# Install Instructions

To build this example you will need an install SAMRAI version with CMake support.
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

Note that you will need a MPI, HDF5 and SILO.

Then download invoke this CMakefile with a proper SAMRAI path...

```
cmake .. -DCMAKE_PREFIX_PATH="${SAMRAI_INSTALL_PATH}/share/SAMRAI/cmake/"
```