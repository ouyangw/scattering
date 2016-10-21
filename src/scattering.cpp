#include "scattering.hpp"
#include "exception.hpp"
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <ostream>

using std::vector;
using utility::Exception;
using std::stringstream;
using std::size_t;
using std::sqrt;
using std::min;
using std::sin;
using std::cos;

namespace scattering_1d
{
Scattering::Scattering(const Conf &conf)
    : m_conf(conf)
    , m_totalE(0)
    , m_dx((conf.x_max - conf.x_min) / (conf.num_xgrid - 1.0))
    , m_k_inc(sqrt(2.0 * conf.mass * conf.kinetic_energy))
    , m_all_eigenvals(conf.num_xgrid * conf.num_states)
    , m_all_eigenvecs(conf.num_xgrid * conf.num_states * conf.num_states)
    , m_energy_refl(0)
    , m_energy_tran(0)
    , m_k_refl(0)
    , m_k_tran(0)
    , m_x(conf.num_xgrid)
    , m_x_fd2(0)
    , m_A()
    , m_B()
    , m_answer()
    , m_aux_vec(0)
{
  for (size_t i(0); i < m_conf.num_xgrid; ++i)
    m_x[i] = m_conf.x_min + i * m_dx;

  // second order finite difference coefficients
  // kinetic operator (-p^2/2m)
  // assuming 5-stencil for now
  const double fd2_factor(1.0 / (2.0 * m_conf.mass * m_dx * m_dx));
  m_x_fd2.resize(3);
  m_x_fd2[0] = 2.5 * fd2_factor;
  m_x_fd2[1] = -4.0 / 3.0 * fd2_factor;
  m_x_fd2[2] = 1.0 / 12.0 * fd2_factor;
}

////////////////////////////////////////////////////////////////////////////////

void Scattering::cal_adiab_states()
{
  typedef Eigen::SelfAdjointEigenSolver<MatrixAdaptor::el_mat_type>
      eigen_solver_type;
  typedef eigen_solver_type::RealVectorType eval_type;
  typedef eigen_solver_type::EigenvectorsType evec_type;
  const size_t H_size(m_conf.num_states * m_conf.num_states);
  MatrixAdaptor::el_mat_type H(m_conf.num_states, m_conf.num_states);
  MatrixAdaptor mat(H);
  eigen_solver_type solver;
  eval_type vals(m_conf.num_states);
  evec_type vecs(m_conf.num_states, m_conf.num_states);

  // calculate adiabatic states for each x
  for (size_t ix(0); ix < m_conf.num_xgrid; ++ix) {
    H.setZero();
    m_conf.eh_builder_ptr->build_H(mat, m_x[ix]);
    solver.compute(H);
    if (solver.info() != Eigen::Success) {
      stringstream error_ss;
      error_ss << "Error: SelfAdjointEigenSolver error: " << solver.info()
               << '\n';
      throw Exception(error_ss.str());
    }
    vals = solver.eigenvalues();
    vecs = solver.eigenvectors();

    // store the adiabitc state eigensystem
    for (size_t is(0); is < m_conf.num_states; ++is) {
      m_all_eigenvals[ix * m_conf.num_states + is] = vals(is);
      for (size_t js(0); js < m_conf.num_states; ++js)
        m_all_eigenvecs[ix * H_size + is * m_conf.num_states + js] =
            vecs(js, is);
    }
  }

  // make wavefunction consistent
  for (size_t ix(1); ix < m_conf.num_xgrid; ++ix) {
    for (size_t is(0); is < m_conf.num_states; ++is) {
      double overlap(0);
      const size_t idx_offset(ix * H_size + is * m_conf.num_states);
      for (size_t js(0); js < m_conf.num_states; ++js)
        overlap += m_all_eigenvecs[idx_offset + js] *
                   m_all_eigenvecs[idx_offset - H_size + js];
      if (overlap < 0)
        for (size_t js(0); js < m_conf.num_states; ++js)
          m_all_eigenvecs[idx_offset + js] = -m_all_eigenvecs[idx_offset + js];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Scattering::setup_channels()
{
  m_totalE = m_conf.kinetic_energy + m_all_eigenvals[m_conf.inc_state];
  for (size_t is(0); is < m_conf.num_states; ++is) {
    if (m_all_eigenvals[is] < m_totalE) {
      m_energy_refl.push_back(m_all_eigenvals[is]);
      m_k_refl.push_back(
          sqrt(2.0 * m_conf.mass * (m_totalE - m_all_eigenvals[is])));
    }
  }
  const size_t idx_offset((m_conf.num_xgrid - 1) * m_conf.num_states);
  for (size_t is(0); is < m_conf.num_states; ++is) {
    if (m_all_eigenvals[idx_offset + is] < m_totalE) {
      m_energy_tran.push_back(m_all_eigenvals[idx_offset + is]);
      m_k_tran.push_back(sqrt(2.0 * m_conf.mass *
                              (m_totalE - m_all_eigenvals[idx_offset + is])));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void Scattering::setup_equation()
{
  // allocate the input for solver
  const size_t main_diagnal_element_count(m_conf.num_states * m_conf.num_xgrid);
  const size_t main_offdiag_element_count(
      m_conf.num_states * m_conf.num_states * 2 * (2 * m_conf.num_xgrid - 3));
  const size_t channel_element_count(5 * (m_k_refl.size() + m_k_tran.size()));
  const size_t A_dim(main_diagnal_element_count + main_offdiag_element_count +
                     channel_element_count);
  const size_t B_dim(m_conf.num_xgrid * m_conf.num_states +
                     m_energy_refl.size() + m_energy_tran.size());
  m_A.resize(B_dim, B_dim);
  m_aux_vec.reserve(A_dim);
  m_B.resize(B_dim);
  m_B.setZero();
  m_answer.resize(B_dim);

  // global indexing
  size_t Aidx(0);

  // position index is slow, electronic state index is fast
  // store in natural format in tmps then convert to CSR3
  
  // setup A
  // diagonal elements
  for (size_t ix(0); ix < m_conf.num_xgrid; ++ix) {
    const size_t idx_offset(ix * m_conf.num_states);
    for (size_t is(0); is < m_conf.num_states; ++is) {
      const size_t idx(idx_offset + is);
      m_aux_vec.push_back(
          triplet_type(idx, idx, m_all_eigenvals[idx] + m_x_fd2[0] - m_totalE));
      ++Aidx;
    }
  }
  // check the number of element inserted
  if (Aidx != main_diagnal_element_count) {
    stringstream error_ss;
    error_ss << "Error setup equation: diagonal element sum is wrong: "
             << Aidx << " expected: " << main_diagnal_element_count
             << '\n';
    throw Exception(error_ss.str());
  }

  // off-diagonal elements
  const size_t H_size(m_conf.num_states * m_conf.num_states);
  for (size_t ix1(0); ix1 < m_conf.num_xgrid - 1; ++ix1) {
    const size_t idx_offset_x1(ix1 * H_size);
    const size_t idx2_wall(min(ix1 + m_x_fd2.size(), m_conf.num_xgrid));
    for (size_t ix2(ix1 + 1); ix2 < idx2_wall; ++ix2) {
      const size_t idx_offset_x2(ix2 * H_size);
      for (size_t is1(0); is1 < m_conf.num_states; ++is1) {
        const size_t idx_offset1(idx_offset_x1 + is1 * m_conf.num_states);
        for (size_t is2(0); is2 < m_conf.num_states; ++is2) {
          const size_t idx_offset2(idx_offset_x2 + is2 * m_conf.num_states);
          double overlap(0);
          for (size_t js(0); js < m_conf.num_states; ++js)
            overlap += m_all_eigenvecs[idx_offset1 + js] *
                       m_all_eigenvecs[idx_offset2 + js];
          const size_t nat_idx1(ix1 * m_conf.num_states + is1);
          const size_t nat_idx2(ix2 * m_conf.num_states + is2);
          const double value(overlap * m_x_fd2[ix2 - ix1]);
          m_aux_vec.push_back(triplet_type(nat_idx1, nat_idx2, value));
          m_aux_vec.push_back(triplet_type(nat_idx2, nat_idx1, value));
          Aidx += 2;
        }
      }
    }
  }
  // check the number of element inserted
  if (Aidx != main_diagnal_element_count + main_offdiag_element_count) {
    stringstream error_ss;
    error_ss << "Error setup equation: off-diagonal element sum is wrong: "
             << Aidx << " expected: "
             << main_diagnal_element_count + main_offdiag_element_count
             << '\n';
    throw Exception(error_ss.str());
  }

  // channels
  size_t channel_offset(m_conf.num_xgrid * m_conf.num_states);
  for (size_t ir(0); ir < m_k_refl.size(); ++ir) {
    m_aux_vec.push_back(triplet_type(
        ir, channel_offset + ir,
        complex_type(m_x_fd2[1] + m_x_fd2[2] * cos(m_k_refl[ir] * m_dx),
                     m_x_fd2[2] * sin(m_k_refl[ir] * m_dx))));
    ++Aidx;

    m_aux_vec.push_back(
        triplet_type(m_conf.num_states + ir, channel_offset + ir, m_x_fd2[2]));
    ++Aidx;

    m_aux_vec.push_back(triplet_type(channel_offset + ir, ir, m_x_fd2[1]));
    ++Aidx;

    m_aux_vec.push_back(
        triplet_type(channel_offset + ir, m_conf.num_states + ir, m_x_fd2[2]));
    ++Aidx;

    m_aux_vec.push_back(triplet_type(
        channel_offset + ir, channel_offset + ir,
        complex_type(m_energy_refl[ir] - m_totalE + m_x_fd2[0] +
                         m_x_fd2[1] * cos(m_k_refl[ir] * m_dx) +
                         m_x_fd2[2] * cos(2.0 * m_k_refl[ir] * m_dx),
                     m_x_fd2[1] * sin(m_k_refl[ir] * m_dx) +
                         m_x_fd2[2] * sin(2.0 * m_k_refl[ir] * m_dx))));
    ++Aidx;
  }
  channel_offset += m_k_refl.size();
  const size_t Nx_offset((m_conf.num_xgrid - 1) * m_conf.num_states);
  const size_t Nx_1_offset((m_conf.num_xgrid - 2) * m_conf.num_states);
  for (size_t it(0); it < m_k_tran.size(); ++it) {
    m_aux_vec.push_back(triplet_type(
        Nx_offset + it, channel_offset + it,
        complex_type(m_x_fd2[1] + m_x_fd2[2] * cos(m_k_tran[it] * m_dx),
                     m_x_fd2[2] * sin(m_k_tran[it] * m_dx))));
    ++Aidx;

    m_aux_vec.push_back(
        triplet_type(Nx_1_offset + it, channel_offset + it, m_x_fd2[2]));
    ++Aidx;

    m_aux_vec.push_back(
        triplet_type(channel_offset + it, Nx_offset + it, m_x_fd2[1]));
    ++Aidx;

    m_aux_vec.push_back(
        triplet_type(channel_offset + it, Nx_1_offset + it, m_x_fd2[2]));
    ++Aidx;

    m_aux_vec.push_back(triplet_type(
        channel_offset + it, channel_offset + it,
        complex_type(m_energy_tran[it] - m_totalE + m_x_fd2[0] +
                         m_x_fd2[1] * cos(m_k_tran[it] * m_dx) +
                         m_x_fd2[2] * cos(2.0 * m_k_tran[it] * m_dx),
                     m_x_fd2[1] * sin(m_k_tran[it] * m_dx) +
                         m_x_fd2[2] * sin(2.0 * m_k_tran[it] * m_dx))));
    ++Aidx;
  }
  // check the number of element inserted
  if (Aidx !=
      main_diagnal_element_count + main_offdiag_element_count +
          channel_element_count) {
    stringstream error_ss;
    error_ss << "Error setup equation: off-diagonal element sum is wrong: "
             << Aidx << " expected: "
             << main_diagnal_element_count + main_offdiag_element_count +
                    channel_element_count
             << '\n';
    throw Exception(error_ss.str());
  }

  m_A.setFromTriplets(m_aux_vec.begin(), m_aux_vec.end());

  // setup B
  m_B(m_conf.inc_state) =
      complex_type(-(m_x_fd2[1] + m_x_fd2[2] * cos(-m_k_inc * m_dx)),
                   -m_x_fd2[2] * sin(-m_k_inc * m_dx));

  m_B(m_conf.num_states + m_conf.inc_state) = -m_x_fd2[2];

  m_B(m_conf.num_xgrid * m_conf.num_states + m_conf.inc_state) =
      complex_type(m_totalE - m_all_eigenvals[m_conf.inc_state] -
                       (m_x_fd2[0] + m_x_fd2[1] * cos(-m_k_inc * m_dx) +
                        m_x_fd2[2] * cos(-2.0 * m_k_inc * m_dx)),
                   -(m_x_fd2[1] * sin(-m_k_inc * m_dx) +
                     m_x_fd2[2] * sin(-2.0 * m_k_inc * m_dx)));
}

////////////////////////////////////////////////////////////////////////////////

void Scattering::solve_equation(vector<Data> &refl, vector<Data> &tran)
{
  // clear cache eigenvectors
  vector<double>(0).swap(m_all_eigenvecs);

  // include all the states including the energetically nonaccessible states
  refl.resize(m_conf.num_states);
  tran.resize(m_conf.num_states);
  const size_t eigenvals_offset((m_conf.num_xgrid - 1) * m_conf.num_states);
  for (size_t is(0); is < m_conf.num_states; ++is) {
    refl[is].energy = m_all_eigenvals[is];
    tran[is].energy = m_all_eigenvals[eigenvals_offset + is];
  }

  // clear cache eigenvalues
  vector<double>(0).swap(m_all_eigenvals);

  // solve the equation
  Eigen::SparseLU<sparse_mat_type> solver;

  solver.analyzePattern(m_A);
  solver.factorize(m_A);
  m_answer = solver.solve(m_B);

  size_t channel_offset(m_conf.num_xgrid * m_conf.num_states);
  for (size_t ir(0); ir < m_k_refl.size(); ++ir) {
    refl[ir].amplitude_real = m_answer(channel_offset + ir).real();
    refl[ir].amplitude_imag = m_answer(channel_offset + ir).imag();
    refl[ir].raw_norm = std::norm(m_answer(channel_offset + ir));
    refl[ir].normalized_norm = refl[ir].raw_norm * m_k_refl[ir] / m_k_inc;
  }
  channel_offset += m_k_refl.size();
  for (size_t it(0); it < m_k_tran.size(); ++it) {
    tran[it].amplitude_real = m_answer(channel_offset + it).real();
    tran[it].amplitude_imag = m_answer(channel_offset + it).imag();
    tran[it].raw_norm = std::norm(m_answer(channel_offset + it));
    tran[it].normalized_norm = tran[it].raw_norm * m_k_tran[it] / m_k_inc;
  }
}

////////////////////////////////////////////////////////////////////////////////

// one-based indexing!
void Scattering::print_full_AB(std::ostream &os) const
{
  typedef Eigen::SparseMatrix<sparse_mat_type::Scalar, Eigen::RowMajor,
                              sparse_mat_type::Index>
      row_mat_type;
  row_mat_type tmp(m_A);
  typedef sparse_mat_type::Index idx_type;
  for (idx_type row(0); row < tmp.outerSize(); ++row) {
    idx_type col(0);
    for (row_mat_type::InnerIterator it(tmp, row); it; ++it) {
      for (; col < it.index(); ++col)
        os << "0 0 ";
      os << it.value().real() << ' ' << it.value().imag() << ' ';
      ++col;
    }
    for (; col < tmp.cols(); ++col)
      os << "0 0 ";
    os << m_B(row).real() << ' ' << m_B(row).imag() << '\n';
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace scattering_1d
