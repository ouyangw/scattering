#ifndef LAPACK_HPP
#define LAPACK_HPP
#include "complex.hpp"
extern "C" {
void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *,
            int *);
void zgesv_(int *, int *, scattering_1d::Complex *, int *, int *,
           scattering_1d::Complex *, int *, int *);
}
#endif
