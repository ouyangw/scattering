# Examples

This directory contains all the examples for the library.

## Contents

The Tully's model problems are published in his paper
[[J. Chem. Phys. **93**, 1061 (1990)](http://dx.doi.org/10.1063/1.459170)].

Directory | Content | Time<sup>1</sup>
--- | --- | ---
`tully_1` | Tully #1 with one incoming momentum | 0.02s
`tully_1_print` | Tully #1 with one incoming momentum. Print the matrix A and vector B to file `tully_1_print.dat` | 0.007s
`tully_1_range` | Tully #1 with a range of incoming momenta. Reproduce the result in the paper. Results in file `tully_1_range.dat`. | 32s
`tully_1_check` | Tully #1 with correct and wrong implementation | 0.007s
`tully_2` | Tully #2 with one incoming momentum | 0.055s
`tully_2_range` | Tully #2 with a range of incoming momenta. Reproduce the result in the paper. Results in file `tully_2_range.dat`. | 288s
`tully_3` | Tully #3 with one incoming momentum | 0.03s
`tully_3_range` | Tully #3 with a range of incoming momenta. Reproduce the result in the paper. Results in file `tully_3_range.dat`. | 50s

*1. The time is measured on a Mac mini (Mid 2011) with 2.5 GHz Intel Core i5
using `time` commandline utility and the `real` time in the result.*

The data from a range of incoming momenta are written into file with 7 columns:

1. incoming energy
2. incoming momentum
3. reflection coefficient on lower state
4. reflection coefficient on upper state
5. transmission coefficient on lower state
6. transmission coefficient on upper state

The matrix A and vector B have complex numbers. Each complex number is written
as a pair of real value and imaginary value, so each complex number occupies
two columns in the output. The matrix and vector are written such that B is the
last two columns of the total output matrix.

The generated result files are kept in the respective directories, and
the standard output of each program is kept in the file with an extension of
`.stdout`.

## How To Use the Examples

These examples also serve as minimal working codes for users to get started
quickly. A very simple way to write a very simple calculation like the examples
is to make your new directory under `examples/` and copy one of the examples
and modify the electronic Hamiltonian builder and input/output handling. To
build your new code as examples, user need to add an entry at the end of
`CMakeLists.txt` file under `examples/` directory:
```cmake
ADD_EXAMPLE (<your-directory-name>)
```
After these settings, user can run `make examples` in the build directory and
find the new executable of user's own code in `examples/bin/` along with all
examples.
