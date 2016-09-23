#ifndef SCATTERING_1D_MATRIX_ADAPTOR_HPP
#define SCATTERING_1D_MATRIX_ADAPTOR_HPP
#include <vector>
namespace scattering_1d
{
class MatrixAdaptor
{
public:
  void assign(std::size_t i, std::size_t j, double value)
  {
    m_vec[j * m_dim + i] = value;
  }

public:
  MatrixAdaptor(std::vector<double> &vec, std::size_t dim): m_vec(vec), m_dim(dim) {}

private:
  std::vector<double> &m_vec;
  const std::size_t m_dim;
};
} // namespace scattering_1d
#endif
