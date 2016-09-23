#ifndef SCATTERING_1D_HPP
#define SCATTERING_1D_HPP
#include "matrix_adaptor.hpp"
#include <string>
#include <vector>
namespace scattering_1d
{
// Get a brief description of the library (including version info).
std::string get_description();

// Base class for electronic Hamiltonian builder.
// Defined later in the file.
class __attribute__((visibility("default"))) ElectronicHamiltonianBuilder;

// struct of configuration for scattering calculation.
struct Conf {
  std::size_t num_states; // number of electronic states
  std::size_t inc_state; // the index of incoming electronic state (0-based)
  std::size_t num_xgrid; // number of points in position basis
  double x_min, x_max; // the range of position basis
  double kinetic_energy; // incoming nuclear kinetic energy
  double mass; // nuclear mass
  // pointer to hamiltomian builder
  ElectronicHamiltonianBuilder *eh_builder_ptr; 
  Conf();
  // return a string describing the configuration
  // prefix will be added in the beginning of every line
  std::string echo(const std::string &prefix = "#") const;
};

// struct of output data.
struct Data {
  double energy;
  double amplitude_real, amplitude_imag;
  double raw_norm;
  double normalized_norm;
};

// scattering calculation.
// Input:
//   conf - the configuration
// Output:
//   reflection - the reflection data
//   transmission - the transmission data
//
// Every element in the output data vector is a channel
// (reflection or transmission)
void compute(const Conf &conf, std::vector<Data> &reflection,
             std::vector<Data> &transmission);

// Base class for electronic Hamiltonian builder
class ElectronicHamiltonianBuilder
{
public:
  virtual ~ElectronicHamiltonianBuilder();
  // build the upper triangular parts of Hamiltonian
  virtual void build_H(MatrixAdaptor &H_uptri, double x) = 0;
  // short description of the Hamiltonian
  virtual std::string get_description() const = 0;
};

} // namespace scattering_1d
#endif
