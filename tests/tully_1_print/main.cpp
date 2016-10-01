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
  std::cout << conf.echo();

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
