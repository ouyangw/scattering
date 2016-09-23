#include "scattering.hpp"
#include "lapack.hpp"
#include "exception.hpp"
#include <sstream>
#include <cmath>

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
    , m_A(0)
    , m_B(0)
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
  const size_t H_size(m_conf.num_states * m_conf.num_states);
  int info(0), lwork(-1);
  vector<double> work(1), H(H_size),
      vals(m_conf.num_states);
  char charv('V'), charu('U');
  int matrix_dim(static_cast<int>(m_conf.num_states));

  // query optimal work array size for diagonalization
  m_conf.eh_builder_ptr->build_H(H, m_x[0]);
  dsyev_(&charv, &charu, &matrix_dim, &H[0], &matrix_dim, &vals[0], &work[0],
         &lwork, &info);
  if (info) {
    stringstream error_ss;
    error_ss << "Error: dsyev work size query info: " << info << '\n';
    throw Exception(error_ss.str());
  }
  lwork = static_cast<int>(work[0]);
  work.resize(lwork);

  // calculate adiabatic states for each x
  for (size_t ix(0); ix < m_conf.num_xgrid; ++ix) {
    H.assign(H_size, 0);
    m_conf.eh_builder_ptr->build_H(H, m_x[ix]);
    dsyev_(&charv, &charu, &matrix_dim, &H[0], &matrix_dim, &vals[0], &work[0],
           &lwork, &info);
    if (info) {
      stringstream error_ss;
      error_ss << "Error: diagonalization error info: " << info << '\n';
      throw Exception(error_ss.str());
    }
    // store the adiabitc state eigensystem
    for (size_t is(0); is < m_conf.num_states; ++is) {
      m_all_eigenvals[ix * m_conf.num_states + is] = vals[is];
      for (size_t js(0); js < m_conf.num_states; ++js)
        m_all_eigenvecs[ix * H_size + is * m_conf.num_states + js] =
            H[is * m_conf.num_states + js];
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
  const size_t A_dim(m_conf.num_xgrid * m_conf.num_states +
                     m_energy_refl.size() + m_energy_tran.size());
  m_A.resize(A_dim * A_dim);
  m_B.resize(A_dim);

  // position index is slow, electronic state index is fast
  
  // setup A
  // diagonal elements
  for (size_t ix(0); ix < m_conf.num_xgrid; ++ix) {
    const size_t idx_offset(ix * m_conf.num_states);
    for (size_t is(0); is < m_conf.num_states; ++is) {
      const size_t idx(idx_offset + is);
      m_A[idx * A_dim + idx].real =
          m_all_eigenvals[idx] + m_x_fd2[0] - m_totalE;
    }
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
          const size_t A_idx1((ix2 * m_conf.num_states + is2) * A_dim +
                              ix1 * m_conf.num_states + is1);
          const size_t A_idx2((ix1 * m_conf.num_states + is1) * A_dim +
                              ix2 * m_conf.num_states + is2);
          m_A[A_idx1].real = overlap * m_x_fd2[ix2 - ix1];
          m_A[A_idx2].real = overlap * m_x_fd2[ix2 - ix1];
        }
      }
    }
  }
  // channels
  size_t channel_offset(m_conf.num_xgrid * m_conf.num_states);
  for (size_t ir(0); ir < m_k_refl.size(); ++ir) {
    const size_t idx11((channel_offset + ir) * A_dim + ir),
        idx12((channel_offset + ir) * A_dim + m_conf.num_states + ir),
        idx21(ir * A_dim + channel_offset + ir),
        idx22((m_conf.num_states + ir) * A_dim + channel_offset + ir),
        idx00((channel_offset + ir) * A_dim + channel_offset + ir);
    m_A[idx11].real = m_x_fd2[1] + m_x_fd2[2] * cos(m_k_refl[ir] * m_dx);
    m_A[idx11].imag = m_x_fd2[2] * sin(m_k_refl[ir] * m_dx);
    m_A[idx12].real = m_x_fd2[2];
    m_A[idx21].real = m_x_fd2[1];
    m_A[idx22].real = m_x_fd2[2];
    m_A[idx00].real = m_energy_refl[ir] - m_totalE + m_x_fd2[0] +
                      m_x_fd2[1] * cos(m_k_refl[ir] * m_dx) +
                      m_x_fd2[2] * cos(2.0 * m_k_refl[ir] * m_dx);
    m_A[idx00].imag = m_x_fd2[1] * sin(m_k_refl[ir] * m_dx) +
                      m_x_fd2[2] * sin(2.0 * m_k_refl[ir] * m_dx);
  }
  channel_offset += m_k_refl.size();
  const size_t Nx_offset((m_conf.num_xgrid - 1) * m_conf.num_states);
  const size_t Nx_1_offset((m_conf.num_xgrid - 2) * m_conf.num_states);
  for (size_t it(0); it < m_k_tran.size(); ++it) {
    const size_t idx11((channel_offset + it) * A_dim + Nx_offset + it),
        idx12((channel_offset + it) * A_dim + Nx_1_offset + it),
        idx21((Nx_offset + it) * A_dim + channel_offset + it),
        idx22((Nx_1_offset + it) * A_dim + channel_offset + it),
        idx00((channel_offset + it) * A_dim + channel_offset + it);
    m_A[idx11].real = m_x_fd2[1] + m_x_fd2[2] * cos(m_k_tran[it] * m_dx);
    m_A[idx11].imag = m_x_fd2[2] * sin(m_k_tran[it] * m_dx);
    m_A[idx12].real = m_x_fd2[2];
    m_A[idx21].real = m_x_fd2[1];
    m_A[idx22].real = m_x_fd2[2];
    m_A[idx00].real = m_energy_tran[it] - m_totalE + m_x_fd2[0] +
                      m_x_fd2[1] * cos(m_k_tran[it] * m_dx) +
                      m_x_fd2[2] * cos(2.0 * m_k_tran[it] * m_dx);
    m_A[idx00].imag = m_x_fd2[1] * sin(m_k_tran[it] * m_dx) +
                      m_x_fd2[2] * sin(2.0 * m_k_tran[it] * m_dx);
  }

  // setup B
  m_B[m_conf.inc_state].real =
      -(m_x_fd2[1] + m_x_fd2[2] * cos(-m_k_inc * m_dx));
  m_B[m_conf.inc_state].imag = -m_x_fd2[2] * sin(-m_k_inc * m_dx);
  m_B[m_conf.num_states + m_conf.inc_state].real = -m_x_fd2[2];
  m_B[m_conf.num_xgrid * m_conf.num_states + m_conf.inc_state].real =
      m_totalE - m_all_eigenvals[m_conf.inc_state] -
      (m_x_fd2[0] + m_x_fd2[1] * cos(-m_k_inc * m_dx) +
       m_x_fd2[2] * cos(-2.0 * m_k_inc * m_dx));
  m_B[m_conf.num_xgrid * m_conf.num_states + m_conf.inc_state].imag =
      -(m_x_fd2[1] * sin(-m_k_inc * m_dx) +
        m_x_fd2[2] * sin(-2.0 * m_k_inc * m_dx));
}

////////////////////////////////////////////////////////////////////////////////

void Scattering::solve_equation(vector<Data> &refl, vector<Data> &tran)
{
  // include all the states including the energetically nonaccessible states
  refl.resize(m_conf.num_states);
  tran.resize(m_conf.num_states);
  const size_t eigenvals_offset((m_conf.num_xgrid - 1) * m_conf.num_states);
  for (size_t is(0); is < m_conf.num_states; ++is) {
    refl[is].energy = m_all_eigenvals[is];
    tran[is].energy = m_all_eigenvals[eigenvals_offset + is];
  }
  // clear cache for memory
  vector<double>(0).swap(m_all_eigenvals);
  vector<double>(0).swap(m_all_eigenvecs);

  int Neq(static_cast<int>(m_B.size()));
  int Nrhs(1);
  int info(0);
  vector<int> ipiv(Neq);
  zgesv_(&Neq, &Nrhs, &m_A[0], &Neq, &ipiv[0], &m_B[0], &Neq, &info);

  size_t channel_offset(m_conf.num_xgrid * m_conf.num_states);
  for (size_t ir(0); ir < m_k_refl.size(); ++ir) {
    refl[ir].amplitude_real = m_B[channel_offset + ir].real;
    refl[ir].amplitude_imag = m_B[channel_offset + ir].imag;
    refl[ir].raw_norm = m_B[channel_offset + ir].norm();
    refl[ir].normalized_norm = refl[ir].raw_norm * m_k_refl[ir] / m_k_inc;
  }
  channel_offset += m_k_refl.size();
  for (size_t it(0); it < m_k_tran.size(); ++it) {
    tran[it].amplitude_real = m_B[channel_offset + it].real;
    tran[it].amplitude_imag = m_B[channel_offset + it].imag;
    tran[it].raw_norm = m_B[channel_offset + it].norm();
    tran[it].normalized_norm = tran[it].raw_norm * m_k_tran[it] / m_k_inc;
  }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace scattering_1d
