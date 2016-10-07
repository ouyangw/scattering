#ifndef TULLY_2_HPP
#define TULLY_2_HPP
#include "scattering_1d.hpp"
#include <cmath>
namespace scattering_1d
{
class Tully_2 : public ElectronicHamiltonianBuilder
{
public:
  Tully_2(double A, double B, double C, double D, double E)
      : m_A(A)
      , m_B(B)
      , m_C(C)
      , m_D(D)
      , m_E(E)
  {
  }

  void build_H(std::vector<MatrixElement> &H_uptri, double x)
  {
    H_uptri.push_back(MatrixElement(0, 0, 0));
    H_uptri.push_back(MatrixElement(1, 1, -m_A * std::exp(-m_B * x * x) + m_E));
    H_uptri.push_back(MatrixElement(0, 1, m_C * std::exp(-m_D * x * x)));
  }

  std::string get_description() const
  {
    return "Tully #2. [J. Chem. Phys. 93, 1061 (1990)]";
  }

private:
  double m_A, m_B, m_C, m_D, m_E;
};
}
#endif
