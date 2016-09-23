#ifndef COMPLEX_HPP
#define COMPLEX_HPP
namespace scattering_1d
{
struct Complex {
  double real, imag;
  explicit Complex(double re = 0, double im = 0)
      : real(re)
      , imag(im)
  {
  }
  double norm() const { return real * real + imag * imag; }
};
}
#endif
