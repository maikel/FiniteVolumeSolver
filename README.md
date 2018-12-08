# Reactive Euler Solver

This repository is an example of how to use SAMRAI. It also serves as a first attempt to talk about more concrete types and algorithms.
Questions which are going to be tackled are 

 * What are the right abstractions which we can define to help us in different projects.
 * How to setup a build system and contiguous integration
 * How to write documentation
 * How does proper unit testing look like

# Physical Project Structure

- `docs/` contains a Doxygen.in file to generate documentation
- `examples/` contains exemplary applications using this framework.
- `include/` contains *public* header files which are consumable via some exportable library.
- `src/` contains *private* header and source files to the library.
- `tests/` contains unit test files written with [Catch2](https://github.com/catchorg/Catch2)
- `third_party/` contains third party libraries which are included as submodules in this project.
- `.clang-format` Format definitions for this project. Auto-format files with the tool `clang-format`

# Build Instructions

To build this example you will need an installed [SAMRAI](https://github.com/LLNL/SAMRAI) version with [CMake](https://cmake.org) support.
To install such a SAMRAI version download the newest SAMRAI from git

```
> git clone --single-branch -b feature/blt https://github.com/LLNL/SAMRAI.git SAMRAI/
> cd SAMRAI
> git submodule init
> git submodule update
```

And build it with commands similar to

```
> mkdir build && cd build
> cmake .. -DCMAKE_INSTALL_PREFIX="YOUR_LOCAL_INSTALL_PATH" -DCMAKE_BUILD_TYPE="Release" -DHDF5_ROOT="YOUR_HDF5_PATH" -DENABLE_OPENMP=OFF -DCMAKE_CXX_COMPILER=mpic++
> make -j4 install
```

Tip: You can use `ccmake .` in the build directory to get an interactive view for configurable variables. Press `t` to toggle advanced options.

*Note: You will need MPI and HDF5.*

Then download this repository and invoke this CMakefile with a proper SAMRAI path.
use following commands in the project source directory.

```
> git clone https://git.imp.fu-berlin.de/ag-klein/samrai_example.git Example
> cd Example/
> git submodule init
> git submodule update
> mkdir build
> cd build
> cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH="${SAMRAI_INSTALL_PATH}"
> make
> ./examples/HyperbolicTimeIntegrator
```

Note: The current implementation requires some **C++14** library and language features.
This means you need at least **gcc 6**, **clang 5** or **Xcode 9** and the appropiate standard library.
