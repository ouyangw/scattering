#include "scattering.hpp"
#include "exception.hpp"
#include <sstream>
#include <cmath>
#include <algorithm>
#include <ostream>
#include <mkl.h>

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
// private classes and functions
namespace
{
using namespace scattering_1d;
struct Coor {
  Complex val;
  int row, col;
  Coor();
};

bool operator<(const Coor &, const Coor &);

struct CSR3Adaptor {
  vector<Complex> &values;
  vector<int> &columns, &row_index;
  CSR3Adaptor(vector<Complex> &val, vector<int> &col, vector<int> &rowidx);
};

void coor_to_csr3(vector<Coor> &, CSR3Adaptor &);

void vec_to_H(vector<double> &, size_t, Scattering::element_vec_type &);
} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

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
    , m_iA(0)
    , m_jA(0)
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
  vector<double> work(1), H(H_size), vals(m_conf.num_states);
  char charv('V'), charu('U');
  int matrix_dim(static_cast<int>(m_conf.num_states));
  element_vec_type mat_elements;

  // query optimal work array size for diagonalization
  m_conf.eh_builder_ptr->build_H(mat_elements, m_x[0]);
  vec_to_H(H, m_conf.num_states, mat_elements);
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
    mat_elements.clear();
    m_conf.eh_builder_ptr->build_H(mat_elements, m_x[ix]);
    vec_to_H(H, m_conf.num_states, mat_elements);
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
  // allocate the input for solver
  const size_t main_diagnal_element_count(m_conf.num_states * m_conf.num_xgrid);
  const size_t main_offdiag_element_count(
      m_conf.num_states * m_conf.num_states * 2 * (2 * m_conf.num_xgrid - 3));
  const size_t channel_element_count(5 * (m_k_refl.size() + m_k_tran.size()));
  const size_t A_dim(main_diagnal_element_count + main_offdiag_element_count +
                     channel_element_count);
  const size_t B_dim(m_conf.num_xgrid * m_conf.num_states +
                     m_energy_refl.size() + m_energy_tran.size());
  m_A.resize(A_dim);
  m_iA.resize(B_dim + 1);
  m_jA.resize(A_dim);
  m_B.resize(B_dim);

  CSR3Adaptor csr3(m_A, m_jA, m_iA);

  // allocate tmp variables
  vector<Coor> coor(A_dim);
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
      coor[Aidx].val.real = m_all_eigenvals[idx] + m_x_fd2[0] - m_totalE;
      coor[Aidx].row = coor[Aidx].col = idx;
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
          coor[Aidx].val.real = overlap * m_x_fd2[ix2 - ix1];
          coor[Aidx].row = nat_idx1;
          coor[Aidx].col = nat_idx2;
          ++Aidx;
          coor[Aidx].val.real = coor[Aidx - 1].val.real;
          coor[Aidx].row = nat_idx2;
          coor[Aidx].col = nat_idx1;
          ++Aidx;
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
    coor[Aidx].val.real = m_x_fd2[1] + m_x_fd2[2] * cos(m_k_refl[ir] * m_dx);
    coor[Aidx].val.imag = m_x_fd2[2] * sin(m_k_refl[ir] * m_dx);
    coor[Aidx].row = ir;
    coor[Aidx].col = channel_offset + ir;
    ++Aidx;

    coor[Aidx].val.real = m_x_fd2[2];
    coor[Aidx].row = m_conf.num_states + ir;
    coor[Aidx].col = channel_offset + ir;
    ++Aidx;

    coor[Aidx].val.real = m_x_fd2[1];
    coor[Aidx].row = channel_offset + ir;
    coor[Aidx].col = ir;
    ++Aidx;

    coor[Aidx].val.real = m_x_fd2[2];
    coor[Aidx].row = channel_offset + ir;
    coor[Aidx].col = m_conf.num_states + ir;
    ++Aidx;

    coor[Aidx].val.real = m_energy_refl[ir] - m_totalE + m_x_fd2[0] +
                          m_x_fd2[1] * cos(m_k_refl[ir] * m_dx) +
                          m_x_fd2[2] * cos(2.0 * m_k_refl[ir] * m_dx);
    coor[Aidx].val.imag = m_x_fd2[1] * sin(m_k_refl[ir] * m_dx) +
                          m_x_fd2[2] * sin(2.0 * m_k_refl[ir] * m_dx);
    coor[Aidx].row = coor[Aidx].col = channel_offset + ir;
    ++Aidx;
  }
  channel_offset += m_k_refl.size();
  const size_t Nx_offset((m_conf.num_xgrid - 1) * m_conf.num_states);
  const size_t Nx_1_offset((m_conf.num_xgrid - 2) * m_conf.num_states);
  for (size_t it(0); it < m_k_tran.size(); ++it) {
    coor[Aidx].val.real = m_x_fd2[1] + m_x_fd2[2] * cos(m_k_tran[it] * m_dx);
    coor[Aidx].val.imag = m_x_fd2[2] * sin(m_k_tran[it] * m_dx);
    coor[Aidx].row = Nx_offset + it;
    coor[Aidx].col = channel_offset + it;
    ++Aidx;

    coor[Aidx].val.real = m_x_fd2[2];
    coor[Aidx].row = Nx_1_offset + it;
    coor[Aidx].col = channel_offset + it;
    ++Aidx;

    coor[Aidx].val.real = m_x_fd2[1];
    coor[Aidx].row = channel_offset + it;
    coor[Aidx].col = Nx_offset + it;
    ++Aidx;

    coor[Aidx].val.real = m_x_fd2[2];
    coor[Aidx].row = channel_offset + it;
    coor[Aidx].col = Nx_1_offset + it;
    ++Aidx;

    coor[Aidx].val.real = m_energy_tran[it] - m_totalE + m_x_fd2[0] +
                          m_x_fd2[1] * cos(m_k_tran[it] * m_dx) +
                          m_x_fd2[2] * cos(2.0 * m_k_tran[it] * m_dx);
    coor[Aidx].val.imag = m_x_fd2[1] * sin(m_k_tran[it] * m_dx) +
                          m_x_fd2[2] * sin(2.0 * m_k_tran[it] * m_dx);
    coor[Aidx].row = coor[Aidx].col = channel_offset + it;
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

  coor_to_csr3(coor, csr3);

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

  // one-based indexing for pardiso input data!
  // pardiso init
  int mtype(3);
  vector<void *> pt(64);
  vector<MKL_INT> iparm(64);
  pardisoinit(&pt[0], &mtype, &iparm[0]);

  // pardiso setup
  // iparm[1] = 3; // for OpenMP parallelization
  // iparm[23] = 1; // for parallel factorization with more than 8 cpu
  iparm[26] = 1; // enable checks for CSR3 correctness

  // pardiso parameters
  MKL_INT maxfct(1), mnum(1), n(m_B.size());
  vector<MKL_INT> perm(n);
  MKL_INT nrhs(1), msglvl(0);
  vector<Complex> x(n);
  MKL_INT error(0);

  // phase: Analysis
  MKL_INT phase(11);
  pardiso(&pt[0], &maxfct, &mnum, &mtype, &phase, &n, &m_A[0], &m_iA[0],
          &m_jA[0], &perm[0], &nrhs, &iparm[0], &msglvl, &m_B[0], &x[0],
          &error);
  if (error) {
    // release memory
    phase = -1;
    MKL_INT error_release_mem(0);
    pardiso(&pt[0], &maxfct, &mnum, &mtype, &phase, &n, &m_A[0], &m_iA[0],
            &m_jA[0], &perm[0], &nrhs, &iparm[0], &msglvl, &m_B[0], &x[0],
            &error_release_mem);
    stringstream error_ss;
    error_ss << "Error Pardiso Phase 1: error code: " << error << '\n';
    if (error_release_mem)
      error_ss << "Error Pardiso Memory Release: error code: "
               << error_release_mem << '\n';
    throw Exception(error_ss.str());
  }

  // phase: Numerical factorization; Solve
  phase = 23;
  pardiso(&pt[0], &maxfct, &mnum, &mtype, &phase, &n, &m_A[0], &m_iA[0],
          &m_jA[0], &perm[0], &nrhs, &iparm[0], &msglvl, &m_B[0], &x[0],
          &error);
  if (error) {
    // release memory
    phase = -1;
    MKL_INT error_release_mem(0);
    pardiso(&pt[0], &maxfct, &mnum, &mtype, &phase, &n, &m_A[0], &m_iA[0],
            &m_jA[0], &perm[0], &nrhs, &iparm[0], &msglvl, &m_B[0], &x[0],
            &error_release_mem);
    stringstream error_ss;
    error_ss << "Error Pardiso Phase 23: error code: " << error << '\n';
    if (error_release_mem)
      error_ss << "Error Pardiso Memory Release: error code: "
               << error_release_mem << '\n';
    throw Exception(error_ss.str());
  }

  size_t channel_offset(m_conf.num_xgrid * m_conf.num_states);
  for (size_t ir(0); ir < m_k_refl.size(); ++ir) {
    refl[ir].amplitude_real = x[channel_offset + ir].real;
    refl[ir].amplitude_imag = x[channel_offset + ir].imag;
    refl[ir].raw_norm = x[channel_offset + ir].norm();
    refl[ir].normalized_norm = refl[ir].raw_norm * m_k_refl[ir] / m_k_inc;
  }
  channel_offset += m_k_refl.size();
  for (size_t it(0); it < m_k_tran.size(); ++it) {
    tran[it].amplitude_real = x[channel_offset + it].real;
    tran[it].amplitude_imag = x[channel_offset + it].imag;
    tran[it].raw_norm = x[channel_offset + it].norm();
    tran[it].normalized_norm = tran[it].raw_norm * m_k_tran[it] / m_k_inc;
  }
}

////////////////////////////////////////////////////////////////////////////////

// one-based indexing!
void Scattering::print_full_AB(std::ostream &os) const
{
  typedef std::vector<int>::const_iterator int_iter_type;
  typedef std::vector<Complex>::const_iterator complex_iter_type;
  complex_iter_type B_iter(m_B.begin());
  Complex zero;
  const int dim(*std::max_element(m_jA.begin(), m_jA.end()));
  size_t idx(0);
  for (int_iter_type rowiter(m_iA.begin() + 1); rowiter != m_iA.end();
       ++rowiter, ++B_iter) {
    int col(0);
    for (; idx < *rowiter - 1; ++idx) {
      for (; col < m_jA[idx] - 1; ++col)
        os << zero << ' ';
      os << m_A[idx] << ' ';
      ++col;
    }
    for (; col < dim; ++col)
      os << zero << ' ';
    os << *B_iter << '\n';
  }
}

////////////////////////////////////////////////////////////////////////////////

void Scattering::print_H(const Conf &conf, element_vec_type &elements,
                         std::ostream &os)
{
  const size_t H_size(conf.num_states * conf.num_states);
  vector<double> H(H_size);
  vec_to_H(H, conf.num_states, elements);
  for (size_t row(0); row < conf.num_states; ++row) {
    for (size_t col(0); col < conf.num_states; ++col)
      os << H[col * conf.num_states + row] << ' ';
    os << '\n';
  }
}

////////////////////////////////////////////////////////////////////////////////
// implementations of the private classes and functions
////////////////////////////////////////////////////////////////////////////////

namespace
{
Coor::Coor()
    : val()
    , row(0)
    , col(0)
{
}

////////////////////////////////////////////////////////////////////////////////

CSR3Adaptor::CSR3Adaptor(vector<Complex> &val, vector<int> &col,
                         vector<int> &rowidx)
    : values(val)
    , columns(col)
    , row_index(rowidx)
{
}

////////////////////////////////////////////////////////////////////////////////

// one-based indexing!
void coor_to_csr3(std::vector<Coor> &coor, CSR3Adaptor &csr3)
{
  typedef std::vector<Complex>::iterator complex_iter_type;
  typedef std::vector<Coor>::iterator coor_iter_type;
  typedef std::vector<int>::iterator int_iter_type;
  std::sort(coor.begin(), coor.end());
  csr3.values.resize(coor.size());
  csr3.columns.resize(coor.size());
  csr3.row_index.resize(coor.back().row + 2);
  int rowidx(1), old_row(-1);
  coor_iter_type coor_it(coor.begin());
  complex_iter_type val_it(csr3.values.begin());
  int_iter_type rowidx_it(csr3.row_index.begin());
  int_iter_type col_it(csr3.columns.begin());
  for (; coor_it != coor.end(); ++coor_it, ++val_it, ++col_it, ++rowidx) {
    *val_it = coor_it->val;
    *col_it = coor_it->col + 1;
    if (coor_it->row != old_row) {
      old_row = coor_it->row;
      *rowidx_it = rowidx;
      ++rowidx_it;
    }
  }
  *rowidx_it = rowidx;
}

////////////////////////////////////////////////////////////////////////////////

bool operator<(const Coor &lhs, const Coor &rhs)
{
  if (lhs.row < rhs.row)
    return true;
  else if (lhs.row == rhs.row) {
    if (lhs.col < rhs.col)
      return true;
    else
      return false;
  } else
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void vec_to_H(vector<double> &H, size_t H_dim,
              Scattering::element_vec_type &elements)
{
  for (vector<double>::iterator it(H.begin()); it != H.end(); ++it)
    *it = 0.0;
  for (Scattering::element_vec_type::iterator it(elements.begin());
       it != elements.end(); ++it)
    H[it->j * H_dim + it->i] = it->value;
}
} // anonymous namespace

} // namespace scattering_1d
