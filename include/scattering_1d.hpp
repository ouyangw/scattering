#ifndef SCATTERING_1D_HPP
#define SCATTERING_1D_HPP
#include <string>
#include <vector>
#include <iosfwd>
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
  std::string to_string(const std::string &prefix = "#") const;
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

// print full matrix of A and vector of B in Ax=B
// B is in the last column
void print_full_AB(const Conf &conf, std::ostream &os);

// check if electronic matrix builder follows the rules
//   1. index is not beyond num_sates
//   2. row index is no larger than column index (upper trianglar)
// write the electronic hamiltonian to the stream
// if there is an error, the element vector is written to stream
void check_electronic_hamiltonian_builder(const Conf &conf, double x,
                                          std::ostream &os);

// struct of electronic matrix element
// the (i+1)th row and (j+1)th column of the matrix has a value of "value"
// i and j are zero-based index for row and column
struct MatrixElement {
  std::size_t i, j;
  double value;
  MatrixElement(std::size_t i0, std::size_t j0, double value0)
      : i(i0)
      , j(j0)
      , value(value0)
  {
  }
  MatrixElement()
      : i(0)
      , j(0)
      , value(0)
  {
  }
};

// Base class for electronic Hamiltonian builder
class ElectronicHamiltonianBuilder
{
public:
  virtual ~ElectronicHamiltonianBuilder();
  // build the upper triangular parts of Hamiltonian
  virtual void build_H(std::vector<MatrixElement> &H_uptri, double x) = 0;
  // short description of the Hamiltonian
  virtual std::string get_description() const = 0;
};

} // namespace scattering_1d
#endif
