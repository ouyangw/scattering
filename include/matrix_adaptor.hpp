#ifndef SCATTERING_1D_MATRIX_ADAPTOR_HPP
#define SCATTERING_1D_MATRIX_ADAPTOR_HPP
#include <cstddef>
namespace scattering_1d
{
struct MatrixFormat {
  enum format {ROWMAJOR, COLUMNMAJOR};
};

struct MatrixIndexingBase {
  enum base {ZEROBASED = 0ull, ONEBASED = 1ull};
};

namespace
{
template <typename data_type>
class MatrixAdaptorBase
{
public:
  MatrixAdaptorBase(data_type *ptr): ptr_(ptr) {}

protected:
  data_type *ptr_;
};
} // anonymous namespace

template <typename data_type, std::size_t dim, MatrixIndexingBase::base N,
          MatrixFormat::format F>
class MatrixAdaptor;

template <typename data_type, std::size_t dim, MatrixIndexingBase::base N>
class MatrixAdaptor<data_type, dim, N, MatrixFormat::ROWMAJOR>
    : public MatrixAdaptorBase<data_type>
{
public:
  MatrixAdaptor(data_type *ptr): MatrixAdaptorBase<data_type>(ptr) {}
  data_type &operator()(std::size_t i, std::size_t j)
  {
    return MatrixAdaptorBase<data_type>::ptr_[(i - N) * dim + (j - N)];
  }
};

template <typename data_type, std::size_t dim, MatrixIndexingBase::base N>
class MatrixAdaptor<data_type, dim, N, MatrixFormat::COLUMNMAJOR>
    : public MatrixAdaptorBase<data_type>
{
public:
  MatrixAdaptor(data_type *ptr): MatrixAdaptorBase<data_type>(ptr) {}
  data_type &operator()(std::size_t i, std::size_t j)
  {
    return MatrixAdaptorBase<data_type>::ptr_[(j - N) * dim + (i - N)];
  }
};
} // namespace scattering_1d
#endif
