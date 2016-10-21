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

#include "scattering_1d.hpp"
#include "scattering.hpp"
#include "exception.hpp"
#include <sstream>
#include <vector>
#include <iomanip>

using utility::Exception;
using std::vector;
using std::stringstream;
using std::string;
using std::ostream;
using std::setw;

namespace
{
void verify_conf(const scattering_1d::Conf &conf)
{
  if (conf.num_states == 0)
    throw Exception("Configuration Error: num_states is zero.\n");
  if (conf.inc_state >= conf.num_states)
    throw Exception(
        "Configuration Error: inc_state is not within num_states.\n");
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
string get_description()
{
  return "Scattering calculation for 1D Hamiltonian (v1.0.0 MKL)";
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
string Conf::to_string(const string &prefix) const
{
  stringstream ss;
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
void compute(const Conf &conf, vector<Data> &reflection,
             vector<Data> &transmission)
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
void print_full_AB(const Conf &conf, ostream &os)
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
void check_electronic_hamiltonian_builder(const Conf &conf, double x,
                                          ostream &os)
{
  typedef vector<MatrixElement> vec_type;
  typedef vec_type::iterator iter_type;
  verify_conf(conf);
  bool is_error(false);
  stringstream error_ss;
  vec_type mat_elements;
  conf.eh_builder_ptr->build_H(mat_elements, x);
  for (iter_type it(mat_elements.begin()); it != mat_elements.end(); ++ it) {
    if (it->i >= conf.num_states) {
      error_ss << "Error (check_electronic_hamiltonian_builder): row index is "
                  "not smaller than num_states.\n";
      is_error = true;
      break;
    } else if (it->i > it->j) {
      error_ss << "Error (check_electronic_hamiltonian_builder): row index is "
                  "larger than column index.\n";
      is_error = true;
      break;
    }
  }
  if (is_error) {
    for (iter_type it(mat_elements.begin()); it != mat_elements.end(); ++it)
      os << setw(8) << it->i << ' '
         << setw(8) << it->j << ' '
         << setw(12) << it->value << '\n';
    throw Exception(error_ss.str());
  }
  Scattering::print_H(conf, mat_elements, os);
}

////////////////////////////////////////////////////////////////////////////////

__attribute__((visibility("default")))
ElectronicHamiltonianBuilder::~ElectronicHamiltonianBuilder()
{
}
} // namespace scattering_1d
