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

  virtual void build_H(std::vector<double> &H, double x)
  {
    if (x > 0)
      H[0] = m_A * (1.0 - std::exp(-m_B * x));
    else
      H[0] = -m_A * (1.0 - std::exp(m_B * x));
    H[1] = H[2] = m_C * std::exp(-m_D * x * x);
    H[3] = -H[0];
  }

  virtual std::string get_description() const
  {
    return "Tully #1.";
  }

private:
  double m_A, m_B, m_C, m_D;
};
}
#endif
