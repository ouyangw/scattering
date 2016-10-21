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

#ifndef SCATTERING_HPP
#define SCATTERING_HPP
#include "scattering_1d.hpp"
#include <Eigen/SparseCore>
#include <iosfwd>
#include <complex>
namespace scattering_1d
{
class Scattering
{
public:
  typedef std::complex<double> complex_type;
  typedef Eigen::SparseMatrix<complex_type> sparse_mat_type;
  typedef Eigen::Triplet<complex_type, sparse_mat_type::Index> triplet_type;
  typedef Eigen::Matrix<complex_type, Eigen::Dynamic, 1> vec_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> el_mat_type;
  typedef std::vector<MatrixElement> element_vec_type;

public:
  Scattering(const Conf &conf);
  void cal_adiab_states();
  void setup_channels();
  void setup_equation();
  void solve_equation(std::vector<Data> &refl, std::vector<Data> &tran);
  void print_full_AB(std::ostream &os) const;

public:
  static void print_H(const Conf &conf, element_vec_type &, std::ostream &);

private:
  const Conf &m_conf;
  double m_totalE;
  const double m_dx;
  const double m_k_inc;
  std::vector<double> m_all_eigenvals, m_all_eigenvecs;
  std::vector<double> m_energy_refl, m_energy_tran, m_k_refl, m_k_tran;
  std::vector<double> m_x, m_x_fd2;
  sparse_mat_type m_A;
  vec_type m_B, m_answer;
  std::vector<triplet_type> m_aux_vec;
};
}
#endif
