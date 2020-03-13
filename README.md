# Finite Volume Solver

This C++17 project provides a framework to solve hyperbolic equations with finite volume methods.
The framework assumes an external structured grid library like AMReX or SAMRAI to adaptively refine regions of interest.

One of the main design philosophies is to make the numerical methods re-usable and consistent across different grid implementations.
We provide some standard flux methods which can be used out-of-the-box with a large set of hyperbolic equations.
Furthermore lot of helper classes exist to generate distributed grids for data mangement which simplify and modernize the usage of libraries like AMReX or SAMRAI.

At last, this library is also capable of handling embedded boundaries in a dimensionally split setting as defined in [Klein2009].

# Installation

## Conan Configuration

We use the C++ package manager [conan](https://conan.io) to install dependencies for this library. `conan` is a python3 package and installable via the python package manager `pip3`. To install conan open a terminal and use

```bash
> pip3 install --user conan # install in user directory, might need to add $HOME/.local/bin to $PATH
```

or

```bash
> sudo pip3 install conan   # try to install conan system-wide
```

To test whether conan is installed type `> conan` in your command line window.

Expected Output:
```
> conan
Consumer commands
  install    Installs the requirements specified in a recipe (conanfile.py or conanfile.txt).
  config     Manages Conan configuration.
  get        Gets a file or list a directory of a given reference or package.
  info       Gets information about the dependency graph of a recipe.
  search     Searches package recipes and binaries in the local cache or in a remote.
Creator commands
  new        Creates a new package recipe template with a 'conanfile.py' and optionally,
             'test_package' testing files.
  create     Builds a binary package for a recipe (conanfile.py).
  upload     Uploads a recipe and binary packages to a remote.
  export     Copies the recipe (conanfile.py & associated files) to your local cache.
  export-pkg Exports a recipe, then creates a package from local source and build folders.
  test       Tests a package consuming it from a conanfile.py with a test() method.
Package development commands
  source     Calls your local conanfile.py 'source()' method.
  build      Calls your local conanfile.py 'build()' method.
  package    Calls your local conanfile.py 'package()' method.
  editable   Manages editable packages (package that resides in the user workspace, but are
             consumed as if they were in the cache).
  workspace  Manages a workspace (a set of packages consumed from the user workspace that
             belongs to the same project).
Misc commands
  profile    Lists profiles in the '.conan/profiles' folder, or shows profile details.
  remote     Manages the remote list and the package recipes associated to a remote.
  user       Authenticates against a remote with user/pass, caching the auth token.
  imports    Calls your local conanfile.py or conanfile.txt 'imports' method.
  copy       Copies conan recipes and packages to another user/channel.
  remove     Removes packages or binaries matching pattern from local cache or remote.
  alias      Creates and exports an 'alias package recipe'.
  download   Downloads recipe and binaries to the local cache, without using settings.
  inspect    Displays conanfile attributes, like name, version and options. Works locally, in
             local cache and remote.
  help       Shows help for a specific command.
  graph      Generates and manipulates lock files.

Conan commands. Type "conan <command> -h" for help
```

As a first step you have to create a [conan profile](https://docs.conan.io/en/latest/reference/commands/misc/profile.html), which descibes which tool chain you want to use.
To create an auto-detected tool-chain use the command

```
> conan profile new default --detect 
```

In case of GCC make sure to use the C++11 ABI. Unfortunately conan does not choose this by default and we have to adjust the configuration by

```
> conan profile update settings.compiler.libcxx=libstdc++11 default
> conan profile show default
Configuration for profile default:

[settings]
os=Linux
os_build=Linux
arch=x86_64
arch_build=x86_64
compiler=gcc
compiler.version=9
compiler.libcxx=libstdc++11
build_type=Release
[options]
[build_requires]
[env]
```

Note: We need a minimum GCC version of 8 or a minimum LLVM clang version of 5.

We added some custom installation recipes which conan can use to install dependencies like `AMReX` or `HDF5`. These are stored in a conan reposiory and we need to point conan to this repository. This is done via the command line 

```
> conan remote add finite-volume https://api.bintray.com/conan/fub-agklein/finite-volume
```

## MPI Installation

Make sure to have a installed some MPI implementation such that your current conan profile is able to link against it. This can be achieved by adapting the compiler option of your profile. 

## CMake Installation

CMake is a build generation tool and generates build configuration as Makefiles for example. 

Before we start installing all dependencies with conan we need to install a rather new version of `CMake`. `AMReX` requires a minimal `CMake` version of 3.14 but updates that requirement fairly regular. Forunately `CMake` is quite simple to download and binary packages can be found at https://cmake.org/download/.

## Building the Library

First use git and clone the library to a local directory 
 
```
> git clone git@git.imp.fu-berlin.de:ag-klein/FiniteVolumeSolver.git
```

If you want to build the unit tests you need to pull `Catch2` as a git submodule. Enter the source direction and call

```
> cd FiniteVolumeSolver/
./FiniteVolumeSolver> git submodule update --init
```

This will checkout the `develop` branch by default and create folder named with relative path `./FiniteVolumeSolver`. Next, we create a out-of-source build directory where we want to build the library archives and example binaries. 

```
./FiniteVolumeSolver> mkdir build
./FiniteVolumeSolver> cd build
./FiniteVolumeSolver/build> cd build
```

Inside the `build` directory we use conan to install the dependencies with the options which we want to use. `AMReX` for examples has the following configurable build options:

```bash
AMReX:eb = True|False [True] # enable embedded boundaries / cut-cells
AMReX:omp = True|False [True] # enable OpenMP parallelization 
AMReX:dim = 1|2|3 [3] # spatial dimension used by AMReX
```

To install AMReX with embedded boundaries and without OpenMP support (there is no OpenMp support on Apple for example) use within the build directory

```
./FiniteVolumeSolver/build> conan install <Path-to-FiniteVolumeSolver-Source-Dir> -o AMReX:dim=2 -o AMReX:omp=False 
```

In our case

```
./FiniteVolumeSolver/build> conan install ../ -o AMReX:dim=2 -o AMReX:omp=False 
```

This will look into the file `FiniteVolumeSolver/conanfile.txt` and tries to locally install all dependencies which are listed there. After installing these it creates a `conanfile.cmake` in the build directory which will be read by our `CMakeLists.txt` file. This in turn injects all necessary include and library paths which we need to build our application. Now we use `cmake` to configure our specific build, i.e.

```
./FiniteVolumeSolver/build> cmake ../
```

to configure a Debug build, or  

```
./FiniteVolumeSolver/build> cmake ../ -DCMAKE_BUILD_TYPE=Release
```

for a Release build. On Linux and MacOs this will create Makefiles by default which you can invoke by typing

```
./FiniteVolumeSolver/build> make
```

which should build all targets. For other systems a more general call would be

```
./FiniteVolumeSolver/build> cmake --build .
```

If some targets fails to build feel free to raise an issue in GitLab.

# Overview

## Policy Pattern

To provide customization points for the many algorithmic choices in an AMR-enabled simulation we make heavily use of the so called _policy pattern_.

> **_From Wikipedia_**: In computer programming, the strategy pattern (also known as the policy pattern) is a behavioral software design pattern that enables selecting an algorithm at runtime. Instead of implementing a single algorithm directly, code receives run-time instructions as to which in a family of algorithms to use.

For each algorithmic customization point we define a concept, which is simply a set of syntactic requirements.
Any type which satisfy a particular policy concept can be used as a drop-in replacement for existing algorithms in order to adjust the method to your special needs.
The customization points are chosen to be orthogonal and thus enable to freely concentrate on a single aspect of the method. 

<details>
<summary>Click to read an example for an InitialData concept</summary>

Any type which has a member function called `void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom)` satisfies the `InitialData` concept for usage with the AMReX library.
This means in practise that an object of type `T` can be used in code as in the example

```cpp
MyInitialDataPolicy my_intial_data{};
amrex::MultiFab data = /* obtain AMReX MultiFab from somewhere */
amrex::Geometry geom = /* obtain AMReX Geometry from somewhere */
my_intial_data.InitializeData(data, geom); 
```

The notion of concepts will be part of the C++ language as of the C++20 standard. 
In compilers which will support C++20 concepts already this type requirement can be expressed in the code as

```cpp
template <typename I>
concept InitialData = requires (I& initial_data, amrex::MultiFab& data, const amrex::Geometry& geometry) {
  { initial_data.InitializeData(data, geometry) };
};

// This line checks at compile-time if the class MyInitialDataPolicy satisfies the InitialData concept.
static_assert(InitialData<MyInitialDataPolicy>);
```

[Try it on godbolt!](https://godbolt.org/z/pu0Hh4)
</details>

## Hyperbolic Equations



## GriddingAlgorithm, DimensionalSplitSystemSolver, SplitSystemSourceLevelIntegrator

When setting up a simulation we decompose the problem in roughly two parts.

1. Manage a locally refined grid containing simulation data distributed over multiple MPI ranks.
2. Set up a numerical method to solve the problem on a distributed grid.

A `PatchHierarchy` holds simulation data for each refinement level on a `PatchLevel` and can be generated by a `GriddingAlgorithm`.
A gridding algorithm generates hierarchies by using a `TaggingMethod` policy object which masks cells which need additional refinement.
The present gridding algorithms use AMReX or SAMRAI to generate patches and hierarchies to cover such tagged cells.
In addition to the `TaggingMethod` policy the `GriddingAlgorithm` also need boundary and initial conditions, given by `BoundaryCondition` and `InitialData` policies.
The boundary conditions are used in communication routines to fill the ghost cell layer touching the computational domain.

Multiple solvers can share a common gridding algorithm and will hold a member variable of type `std::shared_ptr<GriddingAlgorithm>`.
Given a shared `GriddingAlgorithm` a solver is being able to refine the patch hierarchy or to fill ghost cell boundaries at any given time of its algorithm.
Copying a `GriddingAlgorithm` by value will deeply copy all data on all MPI ranks. This can be usefull to create fall-back scenarios of a simulation.

## Conservative and Complete States

Each equation defines two kind of states: `Conservative` and `Complete` states. 
The conservative states contain variables which see a hyperbolic time update from `fub::amrex::TimeIntegrator`.
The complete state variables are a super set of the conservative ones and may define more auxialiary member variables.
This makes sense for variables which are auxialiary or needed very often but their computation cost is high. 
The complete state is the location to cache those expensive computations. 

For example, the template class `template <int Rank> class PerfectGas` defines the conservative variables

```cpp
struct -implementation-defined-class-name- {
  double density;
  array<double, 2> momentum;
  double energy;
};
```

and the complete state variables are 

```cpp
struct -implementation-defined-class-name- {
  double density;
  array<double, 2> momentum;
  double energy;
  double pressure;
  double speed_of_sound;
};
```

## Implement a new FluxMethod

A class `FM` satisfies the concept `FluxMethod<Equation, StencilSize>` if the following constraints are satisfied.

```cpp
template <typename FM, typename Equation, int StencilSize>
concept FluxMethod = requires (FM& flux_method, 
                               Conservative<Equation> cons, 
                               Complete<Equation> stencil[StencilSize], 
                               double dx, double dt, int dir) {
  { flux_method.ComputeNumericFlux(cons, stencil, dx, dt, dir) };
  { flux_method.ComputeStableDt(stencil, dx, dir) } -> double;
  { static_cast<const FM&>(flux_method).GetStencilSize() } -> int;
};
```

For example the following class `MusclHancockMethod` satisfies this concept.

```cpp
struct MusclHancockMethod {
    using Conservative = ::Conservative<PerfectGas<2>>;
    using Complete = ::Complete<PerfectGas<2>>;

    void ComputeNumericFlux(Conservative& cons, const Complete* stencil, 
                              double dx, double dt, int dir);
    
    double ComputeStableDt(const Complete* stencil, double dx, int dir);

    int GetStencilSize() const noexcept { return 4; }
};
```

The implemented MUSCL-type flux methods for an equation `Equation` satisfy `FluxMethod<Equation, 4>` and are implemented in two steps.
First we reconstruct a state to the left and to the right the face in question.
This gives a reduction of the stencil array from `Complete state[4]` to `Complete reconstruced_states[2]`.
Secondly, step we call a (lower-order) base method `BaseFM` which satisfies `FluxMethod<Equation, 2>`.
