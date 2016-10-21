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

#include "../tully_2/tully_2.hpp"
#include "scattering_1d.hpp"
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

using std::setw;

int main()
{
  std::cout << scattering_1d::get_description() << '\n';
  scattering_1d::Tully_2 tully2(0.1, 0.28, 0.015, 0.06, 0.05);
  scattering_1d::Conf conf;
  conf.num_states = 2;
  conf.inc_state = 0;
  conf.num_xgrid = 6001;
  conf.x_min = -15;
  conf.x_max = 15;
  conf.mass = 2000.0;
  conf.eh_builder_ptr = &tully2;

  std::ofstream outfile("tully_2_range.dat");
  if (!outfile.good()) {
    std::cerr << "Error: cannot open file tully_2_range.dat.\n";
    return 1;
  }
  outfile << std::setprecision(8);

  std::vector<scattering_1d::Data> refl, tran;

  for (int i(0); i < 8001; ++i) {
    conf.kinetic_energy = (i * 0.0125 + 10) * (i * 0.0125 + 10) * 0.5 / 2000.0;
    scattering_1d::compute(conf, refl, tran);
    const double total_prob(refl[0].normalized_norm + refl[1].normalized_norm +
                            tran[0].normalized_norm + tran[1].normalized_norm);
    outfile << setw(12) << conf.kinetic_energy << ' '
            << setw(12)
            << std::sqrt(2.0 * conf.mass * conf.kinetic_energy) << ' '
            << setw(12) << refl[0].normalized_norm << ' '
            << setw(12) << refl[1].normalized_norm << ' '
            << setw(12) << tran[0].normalized_norm << ' '
            << setw(12) << tran[1].normalized_norm << ' '
            << setw(12) << total_prob << '\n';
  }

  outfile.close();
  return 0;
}
