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

#include "../tully_1/tully_1.hpp"
#include "scattering_1d.hpp"
#include <iostream>
#include <vector>
#include <fstream>

int main()
{
  std::cout << scattering_1d::get_description() << '\n';
  scattering_1d::Tully_1 tully1(0.01, 1.6, 0.005, 1.0);
  scattering_1d::Conf conf;
  conf.num_states = 2;
  conf.inc_state = 0;
  conf.num_xgrid = 11;
  conf.x_min = -10;
  conf.x_max = 10;
  conf.kinetic_energy = 21.0 * 21.0 * 0.5 / 2000.0;
  conf.mass = 2000.0;
  conf.eh_builder_ptr = &tully1;
  std::cout << conf.to_string();

  std::vector<scattering_1d::Data> refl, tran;
  std::ofstream outfile("tully_1_print.dat");
  if (!outfile.good()) {
    std::cerr << "Warning: cannot open file tully_1_print.dat, "
              << "printing to screen.\n\n";
    scattering_1d::print_full_AB(conf, std::cout);
  } else {
    scattering_1d::print_full_AB(conf, outfile);
    outfile.close();
  }

  return 0;
}
