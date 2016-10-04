#ifndef TULLY_3_HPP
#define TULLY_3_HPP
#include "scattering_1d.hpp"
#include <cmath>
namespace scattering_1d
{
class Tully_3 : public ElectronicHamiltonianBuilder
{
public:
  Tully_3(double A, double B, double C)
      : m_A(A)
      , m_B(B)
      , m_C(C)
  {
  }

  void build_H(MatrixAdaptor &H_uptri, double x)
  {
    H_uptri.assign(0, 0, m_A);
    H_uptri.assign(1, 1, -m_A);
    if (x < 0)
      H_uptri.assign(0, 1, m_B * std::exp(m_C * x));
    else
      H_uptri.assign(0, 1, m_B * (2.0 - std::exp(-m_C * x)));
  }

  std::string get_description() const
  {
    return "Tully #3.";
  }

private:
  double m_A, m_B, m_C;
};
}
#endif
