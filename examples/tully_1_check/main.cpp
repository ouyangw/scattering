#include "../tully_1/tully_1.hpp"
#include "tully_1_wrong.hpp"
#include "tully_1_wrong2.hpp"
#include "scattering_1d.hpp"
#include <iostream>
#include <exception>

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

  std::cout << "checking good builder:\n";
  try {
    scattering_1d::check_electronic_hamiltonian_builder(conf, 0.5, std::cout);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  std::cout << '\n';

  scattering_1d::Tully_1_Wrong tully1wrong(0.01, 1.6, 0.005, 1.0);
  conf.eh_builder_ptr = &tully1wrong;
  std::cout << conf.echo();

  std::cout << "checking bad builder:\n";
  try {
    scattering_1d::check_electronic_hamiltonian_builder(conf, 0.5, std::cout);
  } catch (const std::exception &e) {
    std::cerr << "Expected exception for building full matrix:\n"
              << e.what() << '\n';
  }

  scattering_1d::Tully_1_Wrong2 tully1wrong2(0.01, 1.6, 0.005, 1.0);
  conf.eh_builder_ptr = &tully1wrong2;
  std::cout << conf.echo();

  std::cout << "checking bad builder:\n";
  try {
    scattering_1d::check_electronic_hamiltonian_builder(conf, 0.5, std::cout);
  } catch (const std::exception &e) {
    std::cerr << "Expected exception of overflow index:\n"
              << e.what() << '\n';
    return 2;
  }


  return 0;
}
