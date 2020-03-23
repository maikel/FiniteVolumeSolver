# Finite Volume Solver

This C++17 project provides a framework to solve hyperbolic equations with finite volume methods.
The framework assumes an external structured grid library like AMReX or SAMRAI to adaptively refine regions of interest.

One of the main design philosophies is to make the numerical methods re-usable and consistent across different grid implementations.
We provide some standard flux methods which can be used out-of-the-box with a large set of hyperbolic equations.
Furthermore, a lot of helper classes exist to generate distributed grids for data management which simplifies and modernizes the usage of libraries like AMReX or SAMRAI.

At last, this library is also capable of handling embedded boundaries in a dimensionally split setting as defined in [Klein2009].

- [Finite Volume Solver](#finite-volume-solver)
  - [Installation](#installation)
    - [Conan Configuration](#conan-configuration)
    - [MPI Installation](#mpi-installation)
    - [CMake Installation](#cmake-installation)
    - [Building the Library](#building-the-library)
  - [Code Overview](#code-overview)
    - [Policy Pattern](#policy-pattern)
    - [Structure of a Simulation](#structure-of-a-simulation)
    - [Data Storage: PatchHierarchy, GriddingAlgorithm and IntegratorContext](#data-storage-patchhierarchy-griddingalgorithm-and-integratorcontext)
    - [Conservative and Complete States](#conservative-and-complete-states)
    - [FluxMethod, TimeIntegrator, CompleteFromCons](#fluxmethod-timeintegrator-completefromcons)
    - [Implement a new FluxMethod (simple case)](#implement-a-new-fluxmethod-simple-case)

## Installation

### Conan Configuration

We use the C++ package manager [conan](https://conan.io) to install dependencies for this library. `conan` is a python3 package and installable via the python package manager `pip3`. To install `conan` open a terminal and use

```bash
> pip3 install --user conan # install in user directory, might need to add $HOME/.local/bin to $PATH
```

or

```bash
> sudo pip3 install conan   # try to install conan system-wide
```

To test whether conan is installed type `> conan` in your command line window.

Expected Output:

```text
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

As a first step, you have to create a [conan profile](https://docs.conan.io/en/latest/reference/commands/misc/profile.html), which describes which toolchain you want to use.
To create an auto-detected tool-chain use the command

```text
> conan profile new default --detect
```

In the case of GCC make sure to use the C++11 ABI. Unfortunately, conan does not choose this by default and we have to adjust the configuration by

```text
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

We added some custom installation recipes which `conan` can use to install dependencies like `AMReX` or `HDF5`. These are stored in a `conan` repository and we need to point `conan` to this repository. This is done via the command line

```text
> conan remote add finite-volume https://api.bintray.com/conan/fub-agklein/finite-volume
```

### MPI Installation

Make sure to have some MPI implementation on your system such that your current `conan` profile can link against it. This can be achieved by adapting the compiler option of your profile.

### CMake Installation

CMake is a build generation tool and generates build configuration as Makefiles for example.

Before we start installing all dependencies with `conan` we need to install a rather new version of `CMake`. `AMReX` requires a minimal `CMake` version of 3.14. Fortunately `CMake` is quite simple to download and binary packages can be found [here](https://cmake.org/download/).

### Building the Library

First use git and clone the library to a local directory

```bash
> git clone git@git.imp.fu-berlin.de:ag-klein/FiniteVolumeSolver.git
```

If you want to build the unit tests you need to pull `Catch2` as a git submodule. Enter the source direction and call

```bash
> cd FiniteVolumeSolver/
./FiniteVolumeSolver> git submodule update --init
```

This will checkout the `develop` branch by default and create a folder named with relative path `./FiniteVolumeSolver`. Next, we create an out-of-source build directory where we want to build the library archives and example binaries.

```bash
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

To install AMReX with embedded boundaries and without OpenMP support (there is no OpenMP support on Apple for example) use within the build directory

```bash
./FiniteVolumeSolver/build> conan install <Path-to-FiniteVolumeSolver-Source-Dir> -o AMReX:dim=2 -o AMReX:omp=False
```

In our case

```bash
./FiniteVolumeSolver/build> conan install ../ -o AMReX:dim=2 -o AMReX:omp=False
```

This will look into the file `FiniteVolumeSolver/conanfile.txt` and tries to locally install all dependencies which are listed there. After installing these it creates a `conanfile.cmake` in the build directory which will be read by our `CMakeLists.txt` file. This, in turn, injects all necessary include and library paths which we need to build our application. Now we use `cmake` to configure our specific build, i.e.

```bash
./FiniteVolumeSolver/build> cmake ../
```

to configure a Debug build, or

```bash
./FiniteVolumeSolver/build> cmake ../ -DCMAKE_BUILD_TYPE=Release
```

for a Release build. On Linux and macOS this will create Makefiles by default which you can invoke by typing

```bash
./FiniteVolumeSolver/build> make
```

which should build all targets. For other systems, a more general call would be

```bash
./FiniteVolumeSolver/build> cmake --build .
```

If some targets fail to build feel free to raise an issue in GitLab.

## Code Overview

### Policy Pattern

To provide customization points for the many algorithmic choices in an AMR-enabled simulation we make heavy use of the so-called _policy pattern_.

> **_From Wikipedia_**: In computer programming, the strategy pattern (also known as the policy pattern) is a behavioral software design pattern that enables selecting an algorithm at runtime. Instead of implementing a single algorithm directly, code receives run-time instructions as to which in a family of algorithms to use.

For each algorithmic customization point, we define a concept, which is simply a set of syntactic requirements.
Any type which satisfies a particular policy concept can be used as a drop-in replacement for existing algorithms to adjust the method for your special needs.
The customization points are chosen to be orthogonal and thus enable them to freely concentrate on a single aspect of the method.

Any type which has a member function called `void InitializeData(amrex::MultiFab& data, const amrex::Geometry& geom)` satisfies the `InitialData` concept for usage with the AMReX library.
This means in practice that an object of type `T` can be used in code as in the example

```cpp
MyInitialDataPolicy my_intial_data{};
amrex::MultiFab data = /* obtain AMReX MultiFab from somewhere */
amrex::Geometry geom = /* obtain AMReX Geometry from somewhere */
my_intial_data.InitializeData(data, geom);
```

The notion of concepts will be part of the C++ language as of the C++20 standard.
In compilers which will support C++20 concepts already, this type requirement can be expressed in the code as

```cpp
template <typename I>
concept InitialData = requires (I& initial_data, amrex::MultiFab& data, const amrex::Geometry& geometry) {
  { initial_data.InitializeData(data, geometry) };
};

// This line checks at compile-time if the class MyInitialDataPolicy satisfies the InitialData concept.
static_assert(InitialData<MyInitialDataPolicy>);
```

[Try it on godbolt!](https://godbolt.org/z/pu0Hh4)

### Structure of a Simulation

When setting up a simulation we decompose the setup into two parts:

1. Manage a locally refined grid containing simulation data distributed over multiple MPI ranks.
2. Set up a numerical method to solve the problem on a distributed grid.

The algorithmic choice in the second part leads in general to additional requirements for the data storage, such as needing a certain amount of ghost cells.

### Data Storage: PatchHierarchy, GriddingAlgorithm and IntegratorContext

A `PatchHierarchy` allocates simulation data for each refinement level on a `PatchLevel` and can be either generated by a `GriddingAlgorithm` or can be initialized from a checkpoint to restart a simulation. Data on a patch hierarchy shall represent a valid spatial state at some time point and does not contain a ghost cell layer.
The allocated data can be either associated to be cell-, node- or face-centered. The required data allocation is described by an object of type `DataDescription`, which is defined as

```cpp
struct DataDescription {
  int n_cell_components{0};
  int n_node_components{0};
  int n_face_components{0};
  int dimension{AMREX_SPACEDIM};
};
```

**Note:** ***a patch hierarchy does not know what the components represent. It does not have an equation object as context and doesn't need it. Its only responsibility is to manage a certain amount of data over multiple levels and distributed over MPI ranks.***

A gridding algorithm adds some algorithmic choices to the patch hierarchy which are needed for the usual distributed methods in the context of adaptive mesh refinement. The present gridding algorithms use the `AmrCore` class from the AMReX library to generate patches and hierarchies to cover such tagged cells. It owns a patch hierarchy and initializes or modifies it using a `TaggingMethod` policy object which masks cells that need additional refinement. In addition to the `TaggingMethod` policy, a `GriddingAlgorithm` also need boundary and initial conditions, given by `BoundaryCondition` and `InitialData` policy objects. The boundary conditions are used in communication routines to fill the ghost cell layer touching the computational domain, whereas the initial conditions are only used once at the beginning of a simulation.

**Example**: The following code example will create an initialized and adaptively refined patch hierarchy.

```cpp
#include "fub/AMReX.hpp"
#include "fub/Solver.hpp"

// This class satisfies the InitialData policy concept
struct CircleData {
  void InitializeData(amrex::MultiFab& data,
                      const amrex::Geometry& geom) const {
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      amrex::FArrayBox& fab = data[mfi];
      fub::amrex::ForEachIndex(box, [&](int i, int j) {
        const double x = geom.CellCenter(i, 0);
        const double y = geom.CellCenter(j, 1);
        const double norm2 = x * x + y * y;
        constexpr double r2 = 0.25 * 0.25;
        amrex::IntVect iv(i, j);
        if (norm2 < r2) {
          fab(iv, 0) = 3.0;
        } else {
          fab(iv, 0) = 1.0;
        }
      });
    });
  }
};

int main() {
  // This is needed to initialize the AMReX library
  fub::amrex::ScopeGuard guard{};

  constexpr int Dim = AMREX_SPACEDIM;
  static_assert(AMREX_SPACEDIM == 2);

  const std::array<int, Dim> n_cells{128, 128};
  const std::array<double, Dim> xlower{-1.0, -1.0};
  const std::array<double, Dim> xupper{+1.0, +1.0};

  fub::Advection2d equation{{1.0, 1.0}};

  fub::amrex::CartesianGridGeometry geometry;
  geometry.cell_dimensions = n_cells;
  geometry.coordinates = amrex::RealBox(xlower, xupper);
  geometry.periodicity = std::array<int, Dim>{1, 1};

  fub::amrex::PatchHierarchyOptions hier_opts{};
  hier_opts.max_number_of_levels = 4;

  // Construct an empty patch hierarchy which allocates only one component
  // for mass.
  fub::amrex::PatchHierarchy hierarchy(equation, geometry, hier_opts);

  // Our tagging method will tag cells where the relative difference to the
  // neighbor cells is below a certain tolerance
  using State = fub::Advection2d::Complete;
  fub::amrex::GradientDetector gradient{equation,
                                        std::pair{&State::mass, 1e-2}};

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, CircleData{}, gradient);

  // This call will allocate and initialize data on the hierarchy
  gridding->InitializeHierarchy(0.0);

  // Write a plot file
  fub::amrex::WritePlotFile("InitialHierarchy/plt00000",
                            gridding->GetPatchHierarchy(), equation);
}
```

Multiple solvers can share a common gridding algorithm and will hold a member variable of type `std::shared_ptr<GriddingAlgorithm>`.
Given a shared `GriddingAlgorithm` a solver is being able to refine the patch hierarchy or to fill ghost cell boundaries at any given time of its algorithm.
Copying a `GriddingAlgorithm` by value will deeply copy all data on all MPI ranks. This can be useful to create fall-back scenarios of a simulation.

An `IntegratorContext` allocates, manages and provides access to additional data that is needed for every conservative AMR scheme, such as cell-centered and face-centered data arrays with ghost cell layers. It also manages face-centered data on the coarse-fine interface between two refinement levels.

In to addition to managing the data, every integrator context defines the following member functions.

```cpp
struct IntegratorContext {
   /// \brief Replaces the underlying gridding algorithm with the specified one.
  void
  ResetHierarchyConfiguration(std::shared_ptr<GriddingAlgorithm> gridding);

  /// \brief Whenever the gridding algorithm changes the data hierarchy this
  /// function will regrid all distributed helper variables managed by the
  /// context.
  ///
  /// \param[in] level  The level number of the coarsest level which changed its
  /// shape. Regrid all levels finer than level.
  void ResetHierarchyConfiguration(int level = 0);

  /// \brief Copy data on the specified level from the underlying PatchHierarchy
  /// to the scratch.
  ///
  /// This function does not fill the ghost cell layer.
  void CopyDataToScratch(int level_num);

  /// \brief Copy data on the specified level from the scratch to the underlying
  /// PatchHierarchy.
  void CopyScratchToData(int level_num);

  /// \brief Applies the boundary condition for the scratch space on level
  /// `level` in direcition `dir`.
  ///
  /// \param level  The refinement level on which the boundary condition shall
  /// be used.
  void ApplyBoundaryCondition(int level, Direction dir);
  void ApplyBoundaryCondition(int level, Direction dir, BoundaryCondition& bc);

  /// @{
  /// \brief Fills the ghost layer of the scratch data and interpolates in the
  /// coarse fine layer.
  void FillGhostLayerTwoLevels(int level, BoundaryCondition& fbc, int coarse,
                               BoundaryCondition& cbc);
  void FillGhostLayerTwoLevels(int level, int coarse);
  /// @}

  /// @{
  /// \brief Fills the ghost layer of the scratch data and does nothing in the
  /// coarse fine layer.
  void FillGhostLayerSingleLevel(int level, BoundaryCondition& bc);
  void FillGhostLayerSingleLevel(int level);
  /// @}

  /// \brief Returns a estimate for a stable time step size which can be taken
  /// for specified level number in direction dir.
  [[nodiscard]] Duration ComputeStableDt(int level, Direction dir);

  /// \brief Fill the flux MultiFab with numeric fluxes based on current states
  /// in scratch.
  void ComputeNumericFluxes(int level, Duration dt, Direction dir);

  /// \brief Apply a conservative time update for each conservative variable on
  /// the specified level number and direction.
  void UpdateConservatively(int level, Duration dt, Direction dir);

  /// \brief Reconstruct complete state variables from conservative ones.
  void CompleteFromCons(int level, Duration dt);

  /// \brief Accumulate fluxes on the coarse fine interfaces for a specified
  /// fine level number.
  void AccumulateCoarseFineFluxes(int level, double time_scale, Direction dir);

  /// \brief Replace the coarse fluxes by accumulated fine fluxes on the coarse
  /// fine interfaces.
  void ApplyFluxCorrection(int fine, int coarse, Duration dt);

  /// \brief Resets all accumulates fluxes to zero.
  void ResetCoarseFineFluxes(int fine, int coarse);

  /// \brief Coarsen scratch data from a fine level number to a coarse level
  /// number.
  void CoarsenConservatively(int fine, int coarse);
};
```

This is the minimum amount of functions needed to build common conservative AMR schemes on top of such an integrator context.

### Conservative and Complete States

Each equation defines two kinds of states: `Conservative` and `Complete` states.
The conservative states contain variables which see a hyperbolic time update from `fub::amrex::TimeIntegrator`.
The complete state variables are a superset of the conservative ones and may define more auxiliary member variables.
This makes sense for variables which are auxiliary or needed very often but their computation cost is high.
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

### FluxMethod, TimeIntegrator, CompleteFromCons

### Implement a new FluxMethod (simple case)

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

For example, the following class `MusclHancockMethod` satisfies this concept.

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
First, we reconstruct a state to the left and to the right the face in question.
This gives a reduction of the stencil array from `Complete state[4]` to `Complete reconstruced_states[2]`.
Secondly, we call a (lower-order) base method `BaseFM` which satisfies `FluxMethod<Equation, 2>`.
