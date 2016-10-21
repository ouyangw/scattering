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

#include "tully_2.hpp"
#include "scattering_1d.hpp"
#include <iostream>
#include <vector>
#include <iomanip>

using std::setw;

int main()
{
  typedef std::vector<scattering_1d::Data>::iterator iter_type;
  std::cout << scattering_1d::get_description() << '\n';
  scattering_1d::Tully_2 tully2(0.1, 0.28, 0.015, 0.06, 0.05);
  scattering_1d::Conf conf;
  conf.num_states = 2;
  conf.inc_state = 0;
  conf.num_xgrid = 6001;
  conf.x_min = -15;
  conf.x_max = 15;
  conf.kinetic_energy = 21.0 * 21.0 * 0.5 / 2000.0;
  conf.mass = 2000.0;
  conf.eh_builder_ptr = &tully2;
  std::cout << conf.to_string();

  std::vector<scattering_1d::Data> refl, tran;
  scattering_1d::compute(conf, refl, tran);

  double total_prob(0);
  std::cout << "Reflection data:\n";
  for (iter_type it(refl.begin()); it != refl.end(); ++it) {
    std::cout << setw(16)
              << it->energy << ' '
              << setw(16)
              << it->amplitude_real << ' '
              << setw(16)
              << it->amplitude_imag << ' '
              << setw(16)
              << it->raw_norm << ' '
              << setw(16)
              << it->normalized_norm << '\n';
    total_prob += it->normalized_norm;
  }
  std::cout << "\nTransmission data:\n";
  for (iter_type it(tran.begin()); it != tran.end(); ++it) {
    std::cout << setw(16)
              << it->energy << ' '
              << setw(16)
              << it->amplitude_real << ' '
              << setw(16)
              << it->amplitude_imag << ' '
              << setw(16)
              << it->raw_norm << ' '
              << setw(16)
              << it->normalized_norm << '\n';
    total_prob += it->normalized_norm;
  }
  std::cout << "Total propability: " << total_prob << '\n';
  return 0;
}
