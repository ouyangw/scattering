# 1D Scattering Calculation Library (v2 Eigen)

This is a implementation of quantum scattering calculation in a one dimension
system. The algorithm can be found in the Appendix of the paper 
[[J. Chem. Phys. **142**, 084110 (2015)](http://dx.doi.org/10.1063/1.4908032)].
Please cite the paper if you use the code.

This version use Eigen 3 ([http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/))
for matrix diagonalization and solving sparse linear systems. The Eigen library
is not included and user should install Eigen before compiling the codes. The
Eigen library is open-source, free and header-only, which mean user only need
to copy the Eigen headers into a directory of user's choice.

### Features:

- Modular design makes the code easy to use and modify.
- Potential energy surface independent and the API is consistent across all
  versions of the codes.
- Library only so users are free to write their own caller function (the
  simplest case is the `main()` function) which deals with input and output in
  the way the most appropriate for the users.
- Several example codes provided to get users started quickly.
- CMake([https://cmake.org/](https://cmake.org/)) configuration files are
  provided for easy building in different environments.
- Different versions using different linear algebra libraries (in other
  branches) (Eigen, [MKL](mkl/), [LAPACK](lapack/)).
- Eigen's sparse linear solver makes the last step calculation very fast.\*
- The Eigen library claims to be faster than any free BLAS and has
  comparable speed to the best BLAS, namely Intel MKL and GOTO. (see Eigen's
  [FAQ](http://eigen.tuxfamily.org/index.php?title=FAQ#How_does_Eigen_compare_to_BLAS.2FLAPACK.3F))
  \*

### Limitations:

- Only works for one dimension system
- The wavefunction on both the left and right boundary is assumed as plane
  wave, which requires that the potential energy surfaces near boundaries are
  flat.
- Calculation is done in adiabatic basis (users must be careful that the
  initial condition and the results are in adiabatic basis)
- The diagonalization is done in dense matrix, which can potentially be very
  slow if the electronic Hamiltonian has many dimensions. A possible way to
  speed up in such case is to parallelize the diagonalizations in different
  configurations.
- The Eigen library is not parallelized for most of its algorithms so it
  doesn't by itself take advantage of parallel hardware (see Eigen's
  [FAQ](http://eigen.tuxfamily.org/index.php?title=FAQ#How_does_Eigen_compare_to_BLAS.2FLAPACK.3F))
  \*

\**These items are specific to this version.*

## Build and Install

*Note: The instruction here assume CMake 3.5 or later is installed. I highly
recommend using CMake because its easiness for users. That being said, users
are free to compile the codes in any way they find appropriate.*

Download the code by cloning the repository:
```bash
git clone https://TODO
```

Optionally, checkout the branch (version) you would like to use:
```bash
git checkout <branch-name>
```

Optionally, you can choose to have a out-of-source build so the generated files
are in another directory, which avoids polluting the source directory:
```bach
cd ..
mkdir build
cd build
```
Note here I make a directory named `build` in the directory one level up to the
source, but users are free to create arbitrary directory in arbitrary location,
users just need to provide the path to the source in the next step.

Configure the build system using `cmake`:
```bash
cmake <path-to-source>
```
Next section has more advanced ways to configure the build system, but they are
optional.

*Note: If user uses the GUI version of CMake, user should fill in the source and
build directory and click "configure", optionally change the variables, then
click "generate".*

*Note: The build process assumes Unix Makefiles. If other build system is used,
please follow the procedure for that system, but the targets are the same.*

Build the library:
```bash
make
```
Optionally, you can build the examples together by issuing:
```bash
make examples
```
The generated executables will be in `examples/bin/` under the build directory.

Optionally, you can install the library to the preset location:
```bash
make install
```
Note that the examples will not be installed but only the public interface
headers and the library.

## Configure the Build System:

This section is the instruction of more advanced configuration of the build
system using CMake. It is optional and safe to skip.

To configure the build system with more control, the following command can be
used:
```bash
cmake [-L] [-G <generator-name>] [-D var=value] <path-to-source>
```
The `-L` option will show all the cache variables which users can change on
their needs. This option only shows the variables but not modifies them. It is
just nicer to see the available options.

The `-G` option sets the build system you would like to use. For instance,
"Unix Makefiles" for Linux and Unix system. The available generators can be
found by running `cmake` with `--help` option:
```bash
cmake --help
```
If you don't provide this option, CMake will figure out the default one for
your system.

The `-Dval=value` is to change the value of the variables. Most variables can
be changed after the first time running `cmake` and the way to use `cmake` is
the same. In another word, user can run `cmake` multiple times to change the
configuration variables.

One very important variable that must be set in the first time run and cannot be
change in later run of `cmake` is the compiler option: `CMAKE_CXX_COMPILER`.
For instance `-D CMAKE_CXX_COMPILER=/usr/bin/clang++` for a Clang compiler.
If you don't set this variable, `cmake` will figure out itself, probably the
default C++ compiler on your system.

If you want to change the generator or the compiler, you just need to delete the
generated cache files (`CMakeCache.txt` and `CMakeFiles/`) to start over. One
advantage of out-of-source build is that you can delete everything in the build
directory without worrying about accidentally deleting any source code.

### Configuration Variables

The variables start with `LIBSCATTERING_` prefix are variables specific to this
library.

#### CMAKE_INSTALL_PREFIX

A notable CMake built-in variable is `CMAKE_INSTALL_PREFIX`. This variable
contains the path the program will be installed to when running `cmake
install`. If you are planning to install the library, you might want to
change this variable because the default value is `/usr/local/` which
requires root access to write.

#### LIBSCATTERING_CHECK_CXXSTANDARD

This boolean variable dictates whether detecting C++11/14 availability. Default
is `ON`.

#### LIBSCATTERING_CXXSTANDARD

This string variable dictates the C++ standard to use. Default is the highest
detected. Can be set as `98`, `11` or `14`.

#### LIBSCATTERING_EIGEN_INCLUDE

This path variable sets the path to the Eigen headers. Default is
`/usr/local/include/`.

#### LIBSCATTERING_WARNING_FLAG

This string variable sets the warning flag when compiling the code. It must
be a semicolumn(;) separated list.

## The Examples

The `examples/` directory contains several examples showing how to use the
library. To read more about the examples and how to use the examples, please
refer to the [README](examples/) in `examples/`.

## Public Interface

All the codes of the library is in the namespace `scattering_1d`.

### Header Files

The public headers users should be using are `scattering_1d.hpp` and
`matrix_adaptor.hpp`. They are the only headers installed when user run `make
install`. The other headers in `include/` are private interface that are
subject to change without notice.

### std::string get_description()

This function returns a string with a brief description of the library. It helps
users to keep track of the version of the library used in the calculation.

### void compute(const Conf &conf, std::vector<Data> &reflection, std::vector<Data> &transmission)

This is the function that does the scattering calculation. 

Variable | Input/Output | Usage
--- | --- | ---
`conf` | Input | The configuration for the calculation
`reflection` | Output | The reflection coefficients. The data are stored with increasing potential energy.
`transmission` | Output | The transmission coefficients. The data are stored with increasing potential energy.

The data structure of `Conf` and `Data` will be given [below](#struct-conf).

The function will first check the validity of `conf`. If there is an error in
`conf`, the function will throw an exception with the explanation of the error.

### void print_full_AB(const Conf &conf, std::ostream &os);

Rather than doing the scattering calculation, this function prints the
resulting matrix and vector in the final equation Ax=B to output stream. This
function serves as a way to debug the code.

Variable | Input/Output | Usage
--- | --- | ---
`conf` | Input | The configuration for the calculation
`os` | Output | The output stream the matrix and vector are written to

The matrix and vector are written such that B is the last column of the total
output matrix.

The function will first check the validity of `conf`. If there is an error in
`conf`, the function will throw an exception with the explanation of the error.

### struct Conf

This is the struct containing all the configurations for the calculation.

```c++
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
  std::string echo(const std::string &prefix = "#") const;
};
```

The member `std::string echo(const std::string &prefix = "#") const` is
provided to output the content of the instance to a `std::string`. Every line
of the output string is prefixed by `prefix` and a space. So users can dump this
string into their output file for a record, and the `prefix` serves as a comment
marker so users can read the data in other program more easily.

### struct Data

This is the struct of the output data from `compute`.

```c++
struct Data {
  double energy;
  double amplitude_real, amplitude_imag;
  double raw_norm;
  double normalized_norm;
};
```

The variables are self-explanatory. The `energy` is the potential energy for
the reflection or transmission channel (state). The only subtlety is the
`raw_norm` and `normalized_norm`. The `raw_norm` is the plain absolute value of
the amplitude from the calculation result. The `normalized_norm` is the absolute
value normalized against the incoming and outgoing momenta. The detail is in the
cited paper. The `normalized_norm` is the coefficient users should be looking
at.

### class ElectronicHamiltonianBuilder

The class `ElectronicHamiltonianBuilder` is the abstract base class users
should inherit to provide the implementation of building the electronic
Hamiltonian. For C++ novice, it is easier to just follow the 
[examples](examples/).

```c++
class ElectronicHamiltonianBuilder
{
public:
  virtual ~ElectronicHamiltonianBuilder();
  // build the upper triangular parts of Hamiltonian
  virtual void build_H(MatrixAdaptor &H_uptri, double x) = 0;
  // short description of the Hamiltonian
  virtual std::string get_description() const = 0;
};
```

#### virtual void build_H(MatrixAdaptor &H_uptri, double x)

User must implement this member function in their class to build the electronic
Hamiltonian given a nuclear position of `x`.

#### virtual std::string get_description() const

User must implement this member function in their class to return a string of
brief description of the electronic Hamiltonian. This string is not involved in
the scattering calculation but only serve as a output in `Conf::echo` so users
can keep track of what Hamiltonian is used in the calculation.

### class MatrixAdaptor

This is a helper adapter to hide the actual container for the matrix, so the
interface can be made consistent across different versions of the code. 

Users are supposed to use the member function to build up the electronic
Hamiltonian in `ElectronicHamiltonianBuilder::build_H`. Users should not try to
make the adapter nor are they responsible for making the adapter.

#### void assign(std::size_t i, std::size_t j, double value)

This member function assign the matrix element at `(i+1)`th row and `(j+1)`th
column with the `value`. The `i` and `j` are zero-based index for row and
column, respectively.