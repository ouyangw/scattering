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
  std::cout << conf.echo();

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
