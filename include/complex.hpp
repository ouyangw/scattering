///////////////////////////////////////////////////////////////////////////////
//
// Scattering: 1D scattering calculation Library
// Copyright (C) 2016 Joseph Subotnik's group, The University of Pennsylvania
//
// This file is part of Scattering.
//
// Scattering is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// Scattering is distributed in the hope that it will be usefull, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Scattering. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef COMPLEX_HPP
#define COMPLEX_HPP
#include <ostream>
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

////////////////////////////////////////////////////////////////////////////////

inline std::ostream &operator<<(std::ostream &os, const Complex &c)
{
  os << c.real << ' ' << c.imag;
  return os;
}
}
#endif
