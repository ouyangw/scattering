#include "scattering_1d.hpp"
#include "scattering.hpp"
#include "exception.hpp"
#include <sstream>

using utility::Exception;

namespace
{
void verify_conf(const scattering_1d::Conf &conf)
{
  if (conf.num_states == 0)
    throw Exception("Configuration Error: num_states is zero.\n");
  if (conf.num_xgrid == 0)
    throw Exception("Configuration Error: num_xgrid is zero.\n");
  if (conf.x_min >= conf.x_max)
    throw Exception("Configuration Error: x_max is not larger than x_min.\n");
  if (conf.kinetic_energy < 0)
    throw Exception("Configuration Error: kinetic energy is negative.\n");
  if (conf.mass <= 0)
    throw Exception("Configuration Error: nuclear mass is negative.\n");
  if (conf.eh_builder_ptr == NULL)
    throw Exception("Configuration Error: eh_builder_ptr is NULL.\n");
}
} // anonymous namespace

namespace scattering_1d
{
__attribute__((visibility("default")))
std::string get_description()
{
  return "Scattering calculation for 1D Hamiltonian (v2.0.0)";
}

__attribute__((visibility("default")))
Conf::Conf()
    : num_states(0)
    , inc_state(0)
    , num_xgrid(0)
    , x_min(0)
    , x_max(0)
    , kinetic_energy(0)
    , mass(0)
    , eh_builder_ptr(NULL)
{
}

////////////////////////////////////////////////////////////////////////////////

__attribute__((visibility("default")))
std::string Conf::echo(const std::string &prefix) const
{
  std::stringstream ss;
  ss << prefix
     << " ********** Scattering Conf Echo **********\n"
     << prefix << ' '
     << "num states: " << num_states << '\n'
     << prefix << ' '
     << "incoming state: " << inc_state << '\n'
     << prefix << ' '
     << "num xgrid: " << num_xgrid << '\n'
     << prefix << ' '
     << "x range: [" << x_min << ", " << x_max << "]\n"
     << prefix << ' '
     << "kinetic energy: " << kinetic_energy << '\n'
     << prefix << ' '
     << "nuclear mass: " << mass << '\n'
     << prefix << ' '
     << "Hamiltonian builder: " << eh_builder_ptr->get_description() << '\n'
     << prefix
     << " **************** End Echo ****************\n";
  return ss.str();
}

////////////////////////////////////////////////////////////////////////////////

__attribute__((visibility("default")))
void compute(const Conf &conf, std::vector<Data> &reflection,
             std::vector<Data> &transmission)
{
  verify_conf(conf);
  Scattering comp(conf);
  comp.cal_adiab_states();
  comp.setup_channels();
  comp.setup_equation();
  comp.solve_equation(reflection, transmission);
}

////////////////////////////////////////////////////////////////////////////////

__attribute__((visibility("default")))
void print_full_AB(const Conf &conf, std::ostream &os)
{
  verify_conf(conf);
  Scattering comp(conf);
  comp.cal_adiab_states();
  comp.setup_channels();
  comp.setup_equation();
  comp.print_full_AB(os);
}

////////////////////////////////////////////////////////////////////////////////

__attribute__((visibility("default")))
ElectronicHamiltonianBuilder::~ElectronicHamiltonianBuilder()
{
}
} // namespace scattering_1d
