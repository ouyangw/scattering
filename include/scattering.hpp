#ifndef SCATTERING_HPP
#define SCATTERING_HPP
#include "scattering_1d.hpp"
#include <Eigen/SparseCore>
#include <iosfwd>
#include <complex>
namespace scattering_1d
{
class Scattering
{
public:
  Scattering(const Conf &conf);
  void cal_adiab_states();
  void setup_channels();
  void setup_equation();
  void solve_equation(std::vector<Data> &refl, std::vector<Data> &tran);
  void print_full_AB(std::ostream &os) const;

private:
  typedef std::complex<double> complex_type;
  typedef Eigen::SparseMatrix<complex_type> sparse_mat_type;
  typedef Eigen::Triplet<complex_type, sparse_mat_type::Index> triplet_type;
  typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vec_type;
  const Conf &m_conf;
  double m_totalE;
  const double m_dx;
  const double m_k_inc;
  std::vector<double> m_all_eigenvals, m_all_eigenvecs;
  std::vector<double> m_energy_refl, m_energy_tran, m_k_refl, m_k_tran;
  std::vector<double> m_x, m_x_fd2;
  sparse_mat_type m_A;
  vec_type m_B, m_answer;
  std::vector<triplet_type> m_aux_vec;
};
}
#endif
