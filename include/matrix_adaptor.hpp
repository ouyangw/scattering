#ifndef SCATTERING_1D_MATRIX_ADAPTOR_HPP
#define SCATTERING_1D_MATRIX_ADAPTOR_HPP
#include <Eigen/Core>
namespace scattering_1d
{
class MatrixAdaptor
{
public:
  void assign(std::size_t i, std::size_t j, double value)
  {
    m_mat(i, j) = value;
    m_mat(j, i) = value;
  }

public:
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> el_mat_type;
  MatrixAdaptor(el_mat_type &mat): m_mat(mat) {}

private:
  el_mat_type &m_mat;
};
} // namespace scattering_1d
#endif
