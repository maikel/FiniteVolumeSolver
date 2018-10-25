# SAMRAI Example

This repository is an example of how to use SAMRAI. It also serves as a first attempt to talk about more concrete types and algorithms.
Questions which are going to be tackled are 

 * What are the right abstractions which we can define to help us in different projects.
 * How to setup a build system and contiguous integration
 * How to write documentation
 * How does proper unit testing look like

# Install Instructions

To build this example you will need an installed [SAMRAI](https://github.com/LLNL/SAMRAI) version with [CMake](https://cmake.org) support.
To install such a SAMRAI version download the newest SAMRAI from git

```
> git clone --single-branch -b feature/blt https://github.com/LLNL/SAMRAI.git SAMRAI/
> cd SAMRAI
> git submodule init
> git submodule update
```

And build it via something like this

```
> mkdir build && cd build
> cmake .. -DCMAKE_PREFIX_PATH="YOUR_LOCAL_INSTALL_PATH" -DHDF5_ROOT="YOUR_HDF5_PATH" -DENABLE_OPENMP=OFF -DSILO_DIR="YOUR_SILO_PATH" -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc
> make -j4 install
```

*Note: You will need a MPI, HDF5 and SILO.*

Then download this repository and invoke this CMakefile with a proper SAMRAI path.
use following commands in the project source directory.

```
> git clone https://git.imp.fu-berlin.de/ag-klein/samrai_example.git Example
> cd Example/
> git submodule init
> git submodule update
> mkdir build
> cd build
> cmake .. -DCMAKE_PREFIX_PATH="${SAMRAI_INSTALL_PATH}/share/SAMRAI/cmake/"

```