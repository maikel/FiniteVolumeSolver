# Finite Volume Solver

This repository is an example of how to use SAMRAI and AMReX. 
It also serves as a first attempt to talk about more concrete types and algorithms.
Questions which are going to be tackled are 

 * What are the right abstractions which we can define to help us in different projects.
 * How to setup a build system and contiguous integration
 * How to write documentation
 * How does proper unit testing look like

# Physical Project Structure

- `docs/` contains a Doxygen.in file to generate documentation
- `docker/` contains Dockerfiles to create docker images with preinstalled development environments. Read the [Docker Tutorial](https://git.imp.fu-berlin.de/ag-klein/FiniteVolumeSolver/tree/develop/docker) for more information. 
- `examples/` contains exemplary applications using this framework.
- `include/` contains public header files which are consumable via some exportable library.
- `src/` contains private header and source files to the library.
- `tests/` contains unit test files written with [Catch2](https://github.com/catchorg/Catch2)
- `third_party/` contains third party libraries which are included as submodules in this project.
- `.clang-format` Format definitions for this project. Auto-format files with the tool `clang-format`

# Build Instructions

This project has currently the follwing dependencies on third party libraries

- [Eigen3](https://eigen.tuxfamily.org) for basic linear algebra and vectorization.
- [fmtlib](http://fmtlib.net), a library to format string which got partially standardized in C++20.
- [boost](https://www.boost.org) hana for simpler template meta programming and some basic reflection.
- [AMReX](https://amrex-codes.github.io) an optional AMR library.
- [SAMRAI](https://github.com/LLNL/SAMRAI) an optional AMR library.

## How To Build Finite Volume Solver

Given that all dependencies are installed. you have to invoke CMake with all paths configured for your locally installed depedencies.
These paths need to point to directories contiaining the `LibraryConig.cmake` for each `Library` in question.

What follows is an example of how to setup a build of FiniteVolumeSolver for the case that your AMReX installation is contained in a folder `${AMREX_INSTALL_PREFIX}` while Eigen3, boost and fmtlib are installed via a package manager.

```bash
git clone git@git.imp.fu-berlin.de:ag-klein/FiniteVolumeSolver.git FiniteVolumeSolver/
cd FiniteVolumeSolver
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release -DAMReX_DIR="${AMREX_INSTALL_PREFIX}/lib/cmake/AMReX/"
```

Being in your preconfigured build folder, you can see possible configure options and make adjustments via the command `ccmake .` (configure cmake).

As as example, to generate a Xcode project I use the following command line

```bash
cmake -G Xcode ../ -DCMAKE_Fortran_COMPILER=/usr/local/Cellar/gcc/8.2.0/bin/gfortran -DAMReX_DIR=/Users/maikel/amrex/amrex2d-eb-debug/lib/cmake/AMReX/ -DCMAKE_BUILD_TYPE=Debug -DEigen3_DIR=/usr/local/Cellar/eigen/3.3.7/share/eigen3/cmake/ -DCMAKE_EXE_LINKER_FLAGS="-lgfortran -L/usr/local/Cellar/gcc/8.2.0/lib/gcc/8" -DBOOST_ROOT=/usr/local/Cellar/boost/1.69.0/ -DCMAKE_CXX_COMPILER=mpic++
```

## AMReX

To install [AMReX](https://github.com/AMReX-Codes/amrex) use the git develop branch and build it on your local machine. 
AMReX comes in multiple flavors which can be configured via CMake. To use the CutCell methods we need to enable the embedded boundaries extentions to AMReX.

To build AMReX with support for 2D meshes and embedded boundaries and OpenMP parallelisation an invocation can look like

```bash
AMREX_INSTALL_PREFIX="/Your/Path/To/AMReX/Installation/Prefix"
git clone --single-branch -b development https://github.com/AMReX-Codes/amrex AMReX/
cd AMReX/
mkdir build/
cmake ../ -DCMAKE_BUILD_TYPE=Release -DDIM=2 -DCMAKE_INSTALL_PREFIX="${AMREX_INSTALL_PREFIX}" -DENABLE_EB=ON -DENABLE_OMP=ON
make install
```

## SAMRAI

*SAMRAI is currently not maintained*

To build this example you will need an installed [SAMRAI](https://github.com/LLNL/SAMRAI) version with [CMake](https://cmake.org) support.
To install such a SAMRAI version download the newest SAMRAI from git

```bash
git clone --single-branch -b feature/blt https://github.com/LLNL/SAMRAI.git SAMRAI/
cd SAMRAI
git submodule init
git submodule update
```

And build it with commands similar to

```bash
mkdir build && cd build
cmake ../ -DCMAKE_INSTALL_PREFIX="YOUR_LOCAL_INSTALL_PATH" -DCMAKE_BUILD_TYPE="Release" -DHDF5_ROOT="YOUR_HDF5_PATH" -DENABLE_OPENMP=OFF -DCMAKE_CXX_COMPILER=mpic++
make -j4 install
```

Tip: You can use `ccmake .` in the build directory to get an interactive view for configurable variables. Press `t` to toggle advanced options.

*Note: You will need MPI and HDF5.*

Then download this repository and invoke this CMakefile with a proper SAMRAI path.
use following commands in the project source directory.


