#ifndef TULLY_1_HPP
#define TULLY_1_HPP
#include "scattering_1d.hpp"
#include <cmath>
namespace scattering_1d
{
class Tully_1 : public ElectronicHamiltonianBuilder
{
public:
  Tully_1(double A, double B, double C, double D)
      : m_A(A)
      , m_B(B)
      , m_C(C)
      , m_D(D)
  {
  }

  void build_H(MatrixAdaptor &H_uptri, double x)
  {
    double H00(0);
    if (x > 0)
      H00 = m_A * (1.0 - std::exp(-m_B * x));
    else
      H00 = -m_A * (1.0 - std::exp(m_B * x));
    H_uptri.assign(0, 0, H00);
    H_uptri.assign(1, 1, -H00);
    H_uptri.assign(0, 1, m_C * std::exp(-m_D * x * x));
  }

  std::string get_description() const
  {
    return "Tully #1.";
  }

private:
  double m_A, m_B, m_C, m_D;
};
}
#endif
