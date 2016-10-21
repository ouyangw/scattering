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

  void build_H(std::vector<MatrixElement> &H_uptri, double x)
  {
    double H00(0);
    if (x > 0)
      H00 = m_A * (1.0 - std::exp(-m_B * x));
    else
      H00 = -m_A * (1.0 - std::exp(m_B * x));
    H_uptri.push_back(MatrixElement(0, 0, H00));
    H_uptri.push_back(MatrixElement(1, 1, -H00));
    H_uptri.push_back(MatrixElement(0, 1, m_C * std::exp(-m_D * x * x)));
  }

  std::string get_description() const
  {
    return "Tully #1. [J. Chem. Phys. 93, 1061 (1990)]";
  }

private:
  double m_A, m_B, m_C, m_D;
};
}
#endif
