#ifndef SCATTERING_HPP
#define SCATTERING_HPP
#include "scattering_1d.hpp"
#include "complex.hpp"
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

private:
  const Conf &m_conf;
  double m_totalE;
  const double m_dx;
  const double m_k_inc;
  std::vector<double> m_all_eigenvals, m_all_eigenvecs;
  std::vector<double> m_energy_refl, m_energy_tran, m_k_refl, m_k_tran;
  std::vector<double> m_x;
  std::vector<double> m_x_fd2;
  std::vector<Complex> m_A, m_B;
  std::vector<int> m_iA, m_jA;
};
}
#endif
