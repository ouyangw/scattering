#ifndef LAPACK_HPP
#define LAPACK_HPP
extern "C" {
void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *,
            int *);
}
#endif
