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

  void build_H(std::vector<MatrixElement> &H_uptri, double x)
  {
    H_uptri.push_back(MatrixElement(0, 0, m_A));
    H_uptri.push_back(MatrixElement(1, 1, -m_A));
    if (x < 0)
      H_uptri.push_back(MatrixElement(0, 1, m_B * std::exp(m_C * x)));
    else
      H_uptri.push_back(MatrixElement(0, 1, m_B * (2.0 - std::exp(-m_C * x))));
  }

  std::string get_description() const
  {
    return "Tully #3. [J. Chem. Phys. 93, 1061 (1990)]";
  }

private:
  double m_A, m_B, m_C;
};
}
#endif
