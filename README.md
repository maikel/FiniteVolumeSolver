[![Build on Ubuntu 22.04](https://github.com/maikel/FiniteVolumeSolver/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/maikel/FiniteVolumeSolver/actions/workflows/ubuntu.yml)

# Finite Volume Solver

This C++17 project provides a framework to solve hyperbolic equations with finite volume methods.
The framework assumes an external structured grid library like AMReX or SAMRAI to adaptively refine regions of interest.

One of the main design philosophies is to make the numerical methods re-usable and consistent across different grid implementations.
We provide some standard flux methods which can be used out-of-the-box with a large set of hyperbolic equations.
Furthermore, a lot of helper classes exist to generate distributed grids for data management which simplifies and modernizes the usage of libraries like AMReX or SAMRAI.

At last, this library is also capable of handling embedded boundaries in a dimensionally split setting as defined in [Klein2009].

We provide a [Doxygen documentation](http://page.mi.fu-berlin.de/ghastermann/fvs-agklein-dox/).

This framework includes applications that are used in research of combustion processes within a gas turbine. For example, the application AMReX.EB.ConvergentNozzle computes the shock field of a detonation wave that goes through a convergent nozzle within the plenum of a gas turbine

<img src="docs/img/Figure0318.png" />

Installation {#installation}
============================

You can take a look at the GitHub Actions and docker/Containerfile to see the required steps to build this library. Here, we will briefly explain what needs to be done, to compile and link the library against its dependencies.

Conan Configuration {#conan-config}
-----------------------------------

We use the C++ package manager [conan](https://conan.io) to install dependencies for this library. `conan` is a python3 package and installable via the python package manager `pip3`. To install `conan` open a terminal and use the command

```text
> pip3 install --user conan
```

to make user-wide installation. You might need to add `${HOME}/.local/bin` to your `PATH` environment variable in case of a user-wide installation.
If you have root access on your machine you can also consider doing a system-wide installation via the command

```text
> sudo pip3 install conan
```

As a first step, you have to create a [conan profile](https://docs.conan.io/en/latest/reference/commands/misc/profile.html), which describes which toolchain you want to use. To create an auto-detected tool-chain use the command

```text
> conan profile new default --detect
```

In the case of GCC make sure to use the C++11 ABI. Unfortunately, `conan` does not choose this by default and we have to adjust the configuration by

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

Note: The library needs a minimum GCC version of 7 or a minimum LLVM clang version of 5.

We added some custom installation recipes which `conan` can use to install dependencies like [[AMReX]] or [[HDF5]]. These are stored in the `conan` subdirectory within this repository.

To create conan packages for our third-party dependencies use commands such as

```
> conan create /FiniteVolumeSolver/conan/conan-hdf5 HDF5/1.10@local/stable
> conan create /FiniteVolumeSolver/conan/conan-amrex AMReX/development@local/stable
> conan create /FiniteVolumeSolver/conan/conan-vc Vc/1.4.3@local/stable
> conan create /FiniteVolumeSolver/conan/conan-fmt fmt/9.0.0@local/stable
```

MPI Installation {#mpi-install}
-------------------------------

Make sure to have some MPI implementation on your system such that your current `conan` profile can link against it. This can be achieved by adapting the compiler option of your profile.

CMake Installation {#cmake-install}
-----------------------------------

CMake is a build generation tool and generates build configuration as Makefiles for example.

Before we start installing all dependencies with `conan` we need to install a rather new version of `CMake`. `AMReX` requires a minimal `CMake` version of 3.14. Fortunately `CMake` is quite simple to download and binary packages can be found [here](https://cmake.org/download/).

Building the Library {#build-library}
-------------------------------------

First use git and clone the library to a local directory

```text
> git clone https://github.com/maikel/FiniteVolumeSolver.git
```

If you want to build the unit tests you need to pull `Catch2` as a git submodule. Enter the source direction and call

```text
> cd FiniteVolumeSolver/
./FiniteVolumeSolver> git submodule update --init
```

This will checkout the `develop` branch by default and create a folder named with relative path `./FiniteVolumeSolver`. Next, we create an out-of-source build directory where we want to build the library archives and example binaries.

```text
./FiniteVolumeSolver> mkdir build
./FiniteVolumeSolver> cd build
./FiniteVolumeSolver/build> cd build
```

Inside the `build` directory we use conan to install the dependencies with the options which we want to use. `AMReX` for examples has the following configurable build options:

```text
AMReX:eb = True|False [True] # enable embedded boundaries / cut-cells
AMReX:omp = True|False [True] # enable OpenMP parallelization
AMReX:dim = 1|2|3 [3] # spatial dimension used by AMReX
```

To install [[AMReX]] with embedded boundaries and without OpenMP support (there is no OpenMP support on Apple for example) use within the build directory

```text
./FiniteVolumeSolver/build> conan install <Path-to-FiniteVolumeSolver-Source-Dir> -o AMReX:dim=2 -o AMReX:omp=False
```

In our case

```text
./FiniteVolumeSolver/build> conan install ../ -o AMReX:dim=2 -o AMReX:omp=False
```

This will look into the file `FiniteVolumeSolver/conanfile.txt` and tries to locally install all dependencies which are listed there. After installing these it creates a `conanfile.cmake` in the build directory which will be read by our `CMakeLists.txt` file. This, in turn, injects all necessary include and library paths which we need to build our application. Now we use `cmake` to configure our specific build, i.e.

```text
./FiniteVolumeSolver/build> cmake ../
```

to configure a Debug build, or

```text
./FiniteVolumeSolver/build> cmake ../ -DCMAKE_BUILD_TYPE=Release
```

for a Release build. On Linux and macOS this will create Makefiles by default which you can invoke by typing

```text
./FiniteVolumeSolver/build> make
```

which should build all targets. For other systems, a more general call would be

```text
./FiniteVolumeSolver/build> cmake --build .
```

### Useful build flags

- If you want to link to another python version then the default one you can pass
  ```text
    -DPYBIND11_PYTHON_VERSION=3.8
  ```
  or directly specify the python executable by
  ```text
    -DPYTHON_EXECUTABLE=/path/to/python
  ```


If some targets fail to build feel free to raise an issue in GitLab.

Code Overview {#code-overview}
==============================

We use a layered library design to achieve separation of concerns. 
One of the main goals is to have reusable implementations for algorithms which are independent of the concrete storage backend. 
We use namespaces to indicate the dependency on various structured grid libraries such as AMReX or SAMRAI. 
All components inside of the namespace `fub` are independent of the concrete AMR backend. Classes and functions of the namespaces `fub::amrex` or `fub::samrai` depend on their respective libraries.

This library provides tools to write driver files which form complete applications. Each driver makes its own choices on how to handle program options, input files or its output. 
You can find many example applications in the `examples/` subfolder within the Git project.
Every example composes layer for layer until a solver is built and run. 
An application needs to build up its layers in the following order

  1. [PatchHierarchy](#patch-hierarchy) (data storage)
  2. [GriddingAlgorithm](#gridding-algorithm) (algorithmic choice)
  3. [IntegratorContext](#data-integrator-context) (data storage / algorithmic choice)
  4. one or more composed [LevelIntegrator](#level-integrator)s (algorithmic choice)
  5. [SubcycleSolver](#subcycle-solver) (algorithmic choice)

Note: Each layer in this list is implemented in terms of its predecessor.

Policy Pattern {#policy-pattern}
--------------------------------

To provide customization points for the many algorithmic choices in an AMR-enabled simulation we make heavy use of the so-called *policy pattern*.

> **Wikipedia**: In computer programming, the strategy pattern (also known as the policy pattern) is a behavioral software design pattern that enables selecting an algorithm at runtime. Instead of implementing a single algorithm directly, code receives run-time instructions as to which in a family of algorithms to use.

For each algorithmic customization point, we define a concept, which is simply a set of syntactic requirements.
Any type which satisfies a particular policy concept can be used as a drop-in replacement for existing algorithms to adjust the method for your special needs.
The customization points are chosen to be orthogonal and thus enable them to freely concentrate on a single aspect of the method.

<div class="example">
Any class which has a member function called `UpdateConservatively` with a sufficient function signature satisfies the policy concept for `TimeIntegrator<IntegratorContext>`.
This leads to a parametrization of the algorithm in terms of objects.
One could, in the case of a `TimeIntegrator` for example, have multiple implementations for different parallelization strategies, as in the following code snippet

```cpp
fub::amrex::HypberbolicMethod method{};
if (is_gpu_enabled) {
  method.time_integrator 
      = fub::amrex::EulerForwardTimeIntegrator(fub::execution::gpu);
} else {
  method.time_integrator 
      = fub::amrex::EulerForwardTimeIntegrator(fub::execution::openmp);
}
```

The notion of concepts will be part of the C++ language as of the C++20 standard.
With compilers which will support C++20 concepts already, this type requirement can be expressed in the code as in

```cpp
template <typename TI, typename IntegratorContext>
concept TimeIntegrator = std::copy_assignable<TI> && requires (
    TI& time_integrator, IntegratorContext& context, int level, Duration dt,
    Direction dir) {
  { time_integrator.UpdateConservatively(context, level, dt, dir) };
};
```
</div>


Data Storage: PatchHierarchy {#patch-hierarchy}
-----------------------------------------------

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

Note: A patch hierarchy does not know what the components represent. It does not have an equation object as context and does not need it. Its only responsibility is to manage a certain amount of data over multiple levels and distributed over MPI ranks.

Note: [[AMReX]] represents hierarchies very often as pairs of `(Vector<MultiFab*>, Vector<Geometry>)` where the vector is taken over the number of refinement levels. One can interpret a `PatchHierarchy` as an equivalent representation.

Algorithmic Choice: GriddingAlgorithm {#gridding-algorithm}
-----------------------------------------------------

A gridding algorithm adds some algorithmic choices to the patch hierarchy which are needed for the box generation algorithms in the context of adaptive mesh refinement. 
The present gridding algorithms use the `AmrCore` class from the [[AMReX]] library to generate a distributed hierarchy of index boxes. 
It owns a patch hierarchy and initializes or modifies it.
It uses a `TaggingMethod` policy object to mask cells that need additional refinement. 
In addition to the `TaggingMethod` policy, a `GriddingAlgorithm` also needs boundary and initial conditions, given by `BoundaryCondition` and `InitialData` policy objects. 
The boundary conditions are used in communication routines to fill the ghost cell layer touching the computational domain, whereas the initial conditions are only used once at the beginning of a simulation.

The following code example will create an initialized and adaptively refined patch hierarchy.
```cpp
#include "fub/AMReX.hpp"
#include "fub/Solver.hpp"

struct CircleData {
  void InitializeData(fub::amrex::PatchLevel& patch_level,
                      const fub::amrex::GriddingAlgorithm& grid, int level,
                      fub::Duration /*time*/) const {
    const amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
    amrex::MultiFab& data = patch_level.data;
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

  fub::amrex::PatchHierarchy hierarchy(equation, geometry, hier_opts);

  using State = fub::Advection2d::Complete;
  fub::amrex::GradientDetector gradient{equation,
                                        std::pair{&State::mass, 1e-2}};

  auto gridding = std::make_shared<fub::amrex::GriddingAlgorithm>(
      hierarchy, CircleData{}, gradient);
  gridding->InitializeHierarchy(0.0);

  fub::amrex::WritePlotFile("InitialHierarchy/plt00000",
                            gridding->GetPatchHierarchy(), equation);
}
```
This example produces the following figure and can be visualized either in VisIt, Paraview or python-yt.
<center><img width="60%" height="auto" src="docs/img/initial_hierarchy.png"/></center>

Copying a `GriddingAlgorithm` object will deeply copy its data and its policy objects on each MPI rank. This is useful if one wants to make a backup of the current grid state.
The `GriddingAlgorithm` provides three customization points. 
Their concepts `InitialData<GriddingAlgorithm>`, `TaggingMethod<GriddingAlgorithm>` and `BoundaryCondition<GriddingAlgorithm>` are defined as 

```cpp
template <typename I, typename GriddingAlgorithm>
concept InitialData = std::copy_assignable<I> && 
    requires (const I& initial_data, const GriddingAlgorithm& grid, int level) {
        { initial_data.InitializeData(grid, level) };
    }

template <typename T, typename GriddingAlgorithm>
concept TaggingMethod = std::copy_assignable<T> && 
    requires (const T& tagging_method, ::amrex::TagBoxArray& tags, const GriddingAlgorithm& grid, int level) {
        { tagging_method.TagCellsForRefinement(tags, grid, level) };
    }

template <typename BC, typename GriddingAlgorithm>
concept BoundaryCondition = std::copy_assignable<BC> && 
    requires (const BC& boundary_condition, ::amrex::MultiFab& data, const GriddingAlgorithm& grid, int level) {
        { boundary_condition.FillBoundary(data, grid, level) };
    }
```

You find polymorphic value types `AnyInitialData<GriddingAlgorithm>`, `AnyTaggingMethod<GriddingAlgorithm>` and `AnyBoundaryCondition<GriddingAlgorithm>` in the Finite Volume Solver library to match each one of the concepts. These polymorphic types are used to store any object of a class that satisfies those concepts.

<div class="example">
The class `fub::amrex::GriddingAlgorithm` is internally defined in something like

```cpp
class GriddingAlgorithm : private ::amrex::AmrCore {
private:
    PatchHierarchy hierarchy_;
    AnyInitialData initial_data_;
    AnyTaggingMethod tagging_method_;
    AnyBoundaryCondition boundary_condition_;

public:
    /* etc... */
};
```

</div>

Data Storage: IntegratorContext {#data-integrator-context}
----------------------------------------------------------

An `IntegratorContext` allocates, manages and provides access to additional data that is needed for every conservative AMR scheme, such as cell-centered and face-centered data arrays **with ghost layers**. It also manages face-centered data on the coarse-fine interface between two refinement levels.

<div class="example">
The following code snippet shows a how to access cell and face-based data with ghost cells using the [[AMReX]] toolbox. This implements parallelization for GPU with a fallback to an MPI/OpenMP hybrid. This class satisfies the `TimeIntegrator<fub::amrex::IntegratorContext>` concept.

```cpp
class GpuTimeIntegrator {
  /// \brief This function computes a conservative time update for cells in
  /// direction dir
  ///
  /// This specific implementation uses the AMReX tools to do the job
  void UpdateConservatively(IntegratorContext& context, int level, Duration dt,
                            Direction dir) {
    using namespace amrex;
    MultiFab& states = context.GetScratch(level);
    const MultiFab& fluxes = context.GetFluxes(level, dir);
    const int n_conservative_variables = context.GetNConservativeVariables();
    const double dx = context.GetCellSize(level, dir);
    const double lambda = dt.count() / dx;
  #ifdef _OPENMP
  #pragma omp parallel if (Gpu::notInLaunchRegion())
  #endif
    for (MFIter mfi(cells, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      FArrayBox& qs = cells[mfi];
      const FArrayBox& fs = fluxes[mfi];
      const Box cell_box = mfi.growntilebox();
      Array4<double> q = qs.array();
      Array4<const double> f = fs.array();
      ParallelFor(cell_box, n_conservative_variables,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int var) {
                    std::array<int, 3> iv{i, j, k};
                    std::array<int, 3> iL = iv;
                    std::array<int, 3> iR = Shift(iv, dir, 1);
                    q(iv, var) = q(iv, var) + lambda * (f(iL, var) - f(iR, var));
                  });
    }
  }
};
```
</div>

In addition to managing the additional data, every integrator context defines a comprehensive list of functions to deal with the various data locations. The `IntegratorContext` concept is defined by


```cpp
template <typename I>
concept IntegratorContext = std::copy_assignable<I> && requires (
      I& integrator_context, int level, Duration dt, Direction dir,
      double scale, BoundaryCondition& bc) {
  // Access PatchLevel Data
  { integrator_context.GetGriddingAlgorithm() };
  { integrator_context.GetData(level) };
  { integrator_context.GetScratch(level) };
  { integrator_context.GetFluxes(level, dir) };

  // Coarse-fine interface interaction
  { integrator_context.ApplyFluxCorrection(level, level - 1) };
  { integrator_context.ResetCoarseFineFluxes(level, level - 1) };
  { integrator_context.AccumulateCoarseFineFluxes(level, scale, dir) };

  // Action on scratch
  { integrator_context.CoarsenConservatively(level, level - 1) };
  { integrator_context.CompleteFromCons(level) };

  // Copy data
  { integrator_context.CopyDataToScratch(level) };
  { integrator_context.CopyScratchToData(level) };

  // Fill ghost layers within scratch
  { integrator_context.FillGhostLayerTwoLevels(level, bc, level - 1, bc) };
  { integrator_context.FillGhostLayerTwoLevels(level, level - 1) };
  { integrator_context.FillGhostLayerSingleLevel(level, bc) };
  { integrator_context.FillGhostLayerSingleLevel(level) };

  // Reallocate internal data from level 
  { integrator_context.ResetHierarchyConfiguration(level) };
};
```

This is the minimum amount of functions needed to build common conservative AMR schemes on top of an AMR library like [[AMReX]] or [[SAMRAI]]. The implementation of an `IntegratorContext` will, in general, be very library-specific. The layers after the `IntegratorContext` describe parts of the AMR time integration which are implemented on those concepts only.

Nevertheless, we still refine this concept to better fit the dimensionally split time integration.

```cpp
template <typename I>
concept DimensionalSplitIntegratorContext = IntegratorContext<I> && requires (
      I& integrator_context, int level, Duration dt, Direction dir,
      double scale) {
  { integrator_context.ComputeNumericFluxes(level, dt, dir) };
  { integrator_context.ComputeStableDt(level, dir) };
  { integrator_context.UpdateConservatively(level, dir) };
};
```

Each integrator context which is currently implemented satisfies the `DimensionalSplitIntegratorContext`.
The functions `ComputeStableDt`, `ComputeNumericFluxes`, `UpdateConservatively` and `CompleteFromCons` are customization points for those imlpementations. I. e. one needs to provide objects for the policy concepts `FluxMethod<IntegratorContext>`, `TimeIntegrator<IntegratorContext>` and `CompleteFromConsCalculation<IntegratorContext>` to define the behaviour of those functions.

The list of the currently implemented integrator contexts is

 - `fub::amrex::IntegratorContext`
 - `fub::amrex::CompressibleAdvectionIntegratorContext`
 - `fub::amrex::cutcell::IntegratorContext`
 - `fub::amrex::MultiBlockIntegratorContext`
 - `fub::samrai::IntegratorContext`

The `MultiBlockIntegratorContext` class is special, because it is composed of two vectors `std::vector<fub::amrex::IntegratorContext>` and `std::vector<fub::amrex::cutcell::IntegratorContext>` where each integrator context describes a computational domain. Instead of defining policy objects for the composed context one provides them for each member context.

Algorithmic Choice: IntegratorContext {#algorithm-choice-integrator-context}
----------------------------------------------------------------------------

We now describe the customization points for an integrator context.
Its concepts `FluxMethod<IntegratorContext>`, `TimeIntegrator<IntegratorContext>` and `CompleteFromConsCalculation<IntegratorContext>` are defined as 

```cpp
template <typename FM, typename IntegratorContext>
concept FluxMethod = std::copy_assignable<FM> && requires (
    FM& flux_method, IntegratorContext& context, int level, Duration dt,
    Direction dir) {
  { flux_method.ComputeNumericFlux(context, level, dt, dir) };
  { flux_method.ComputeStableDt(context, level, dir) } -> Duration;
  { FM::GetStencilSize() } -> std::integral_constant<int, StencilSize>;
};

template <typename TI, typename IntegratorContext>
concept TimeIntegrator = std::copy_assignable<TI> && requires (
    TI& time_integrator, IntegratorContext& context, int level, Duration dt,
    Direction dir) {
  { time_integrator.UpdateConservatively(context, level, dt, dir) };
};

template <typename R, typename IntegratorContext>
concept CompleteFromConsCalculation = std::copy_assignable<R> && requires(
    R& reconstruct, IntegratorContext& context, int level) {
  { reconstruct.CompleteFromCons(context, level) };
};
```

The member function `ComputeStableDt` shall return a time step size \( \Delta t \) such that the CFL condition for the time integration is satisfied.
The `ComputeNumericFluxes` shall fill the face-centered flux arrays in the specified direction which will be used by the `TimeIntegrator` to update the conservative state variables.

The `TimeIntegrator` provides the usual conservative scheme for cell-based time updates

![equation](https://latex.codecogs.com/svg.image?u^{n&plus;1}_i&space;=&space;u^n_i&space;&plus;&space;\frac{\Delta&space;t}{\Delta&space;x}&space;\left(&space;F^n_{i&space;-&space;\frac&space;12}&space;-&space;F^n_{i&space;&plus;&space;\frac&space;12}&space;\right).)

The numerical fluxes \( F^n_i \) in this update denote the fluxes on the faces of the specified direction and are computed by the `FluxMethod` which is configured for this integrator context.

An AMR integration scheme may use the `CompleteFromConsCalculation` policy object whenever it changes the conservative variables to establish a valid state across all fields in the patch hierarchy.

We collect the three customization points for the `IntegratorContext` in a class which we call `HyperbolicMethod<IntegratorContext>`.
It is defined as

```cpp
template <typename IntegratorContext>
struct HyperbolicMethod {
  AnyFluxMethod<IntegratorContext> flux_method;
  AnyTimeIntegrator<IntegratorContext> time_integrator;
  AnyCompleteFromConsCalculation<IntegratorContext> reconstruct;
};
```

<div class="example">
To construct an `IntegratorContext` object one needs to provide a `GriddingAlgorithm` and a `HyperbolicMethod`.
The following code snippet shows how a concrete integrator context might be constructed in a driver file.

```cpp
void MyMain(const fub::ProgramOptions& opts) {
  std::chrono::steady_clock::time_point wall_time_reference =
      std::chrono::steady_clock::now();

  fub::amrex::ScopeGuard guard{};

  fub::Advection2d equation{{1.0, 1.0}};
  fub::SeverityLogger log = fub::GetInfoLogger();

  fub::amrex::CartesianGridGeometry geometry = fub::GetOptions(opts, "GridGeometry");
  BOOST_LOG(log) << "GridGeometry:";
  geometry.Print(log);

  fub::amrex::PatchHierarchyOptions hier_opts = fub::GetOptions(opts, "PatchHierarchy");
  BOOST_LOG(log) << "PatchHierarchy:";
  hier_opts.Print(log);
 
  using State = fub::Advection2d::Complete;
  fub::amrex::GradientDetector gradient{equation,
                                        std::pair{&State::mass, 1e-3}};

  std::shared_ptr grid = std::make_shared<fub::amrex::GriddingAlgorithm>(
      fub::amrex::PatchHierarchy(equation, geometry, hier_opts), CircleData{},
      fub::amrex::TagAllOf(gradient, fub::amrex::TagBuffer(2)));
  grid->InitializeHierarchy(0.0);

  fub::amrex::HyperbolicMethod method{
      fub::amrex::FluxMethodAdapter(fub::GodunovMethod{equation}),
      fub::amrex::EulerForwardTimeIntegrator(), fub::amrex::NoReconstruction{}};

  fub::amrex::IntegratorContext context(grid, method);
  
  // Now further use the context to do some things...
}
```
</div>

In addition to computing the mathematical correct formula, an object that satisfies `FluxMethod<IntegratorContext>` needs to repeat the logic of how to parallelize and how to iterate over the multi-dimensional index space.
We recognize that these are two orthogonal concerns which are both hard to get right. 
To ease the burden for authors of a `FluxMethod` we introduced some helper classes which separate these concerns. They will be discussed in more detail in <a href="#flux-method-framework">Chapter 3</a>.

Algorithmic Choice: LevelIntegrator {#level-integrator}
-------------------------------------------------------

Level integrators have the task to advance the solution in time for a single refinement level only.
They will be used by larger subcycle schemes which may, in general, do smaller time steps on finer refinement levels.
We define the concept `LevelIntegrator` such that it provides two functions `AdvanceLevelNonRecursively` and `ComputeStableDt`.


```cpp
template <typename L>
concept LevelIntegrator = std::copy_assignable<L> && requires (const L& level_integrator, int level, Duration dt, std::pair<int, int> subcycle) {
  { level_integrator.AdvanceLevelNonRecursively(level, dt, subcycle) } -> Result<void, TimeStepTooLarge>;
  { level_integrator.ComputeStableDt(level) } -> Duration;
};
```

The name `AdvanceLevelNonRecursively` emphasizes that this function will be called within a recursive subcycle algorithm.
We call level integrators which provide an integrator context `HyperbolicSystemLevelIntegrator`.
An integrator context extents level integrators with functions such as `ApplyFluxCorrection`, `CoarsenConservatively` or `CompleteFromCons`.
These additional functions will be used in a subcycle algorithm to maintain the conservation of the solution across refinement levels.

```cpp
template <typename HS>
concept HyperbolicSystemLevelIntegrator = LevelIntegrator<S> && requires (const L& level_integrator) {
  // Access underlying data storage implementation
  { level_integrator.GetIntegratorContext() } -> IntegratorContext;
};
```

This library includes currently the following list of level integrators

- `DimensionalSplitLevelIntegrator` (satisfies `HyperbolicSystemLevelIntegrator`)
- `SystemSourceLevelIntegrator` (satisfies `HyperboilcSystemLevelIntegrator`)
- `BK19LevelIntegrator` (satisfies `HyperbolicSystemLevelIntegrator`)
- `KineticSourceTerm` (satisfies `LevelIntegrator`)
- `AxialSourceTerm` (satisfies `LevelIntegrator`)
- `IgniteDetonation` (satisfies `LevelIntegrator`)

`SystemSourceLevelIntegrator` and `BK19LevelIntegrator` are so-called composed level integrators.
This means, in order to build one of those level integrators you need to provide at least one `HyperbolicSystemLevelIntegrator`.
The only `HyperbolicSystemLevelIntegrator` currently defined in the library is the `DimensionalSplitLevelIntegrator` which applies dimensional operator splitting ontop of a `DimensionalSplitIntegratorContext`.

A `DimensionalSplitLevelIntegrator` is constructible by providing a `DimensionalSplitIntegratorContext` and `SplittingMethod`, such as `GodunovSplitting` or `StrangSplitting`.
A single dimensionally split time step includes the following steps

```cpp
  // Apply boundary condition for the physical boundary only.
  Context::ApplyBoundaryCondition(this_level, dir);

  // Check stable time step size and if the CFL condition is violated then
  // restart the coarse time step
  const Duration local_dt = Context::ComputeStableDt(this_level, dir);
  MPI_Comm comm = GetMpiCommunicator();
  const Duration level_dt = MinAll(comm, local_dt);
  if (level_dt < split_dt) {
    return TimeStepTooLarge{level_dt, this_level};
  }

  // Compute fluxes in the specified direction
  Context::ComputeNumericFluxes(this_level, split_dt, dir);

  // Use the updated fluxes to update cons variables at the "SCRATCH"
  // context.
  Context::UpdateConservatively(this_level, split_dt, dir);

  // A time update happened on conservative variables.
  // We have to re-establish a valid complete state.
  Context::CompleteFromCons(this_level, split_dt);

  // Accumulate the fluxes on the coarse fine interface
  const double scale = split_dt.count();
  Context::AccumulateCoarseFineFluxes(this_level, scale, dir);
```

Algorithmic Choice: SubcycleSovler {#subcycle-solver}
-----------------------------------------------------

Finally, subcycle solvers define in which order level integrators are applied on the refinement level of a hierarchy.
They define the functions `ComputeStableDt` and `AdvanceLevel`.

```cpp
template <typename S>
concept SubcycleSolver = std::copy_assignable<S> && requires (const S& solver, int level, Duration dt) {
  { solver.ComputeStableDt() } -> Duration;
  { solver.AdvanceLevel(level, dt) } -> Result<void, TimeStepTooLarge>;
};
```

The `AdvanceLevel` function advances all levels finer or equal than the specified level.
Most implementations define this function recursively.

<div class=example>
An examplary implementation of an recursive AdvanceLevel function is given here.

```cpp
template <typename LevelIntegrator>
Result<void, TimeStepTooLarge>
SubcycleFineFirstSolver<LevelIntegrator>::AdvanceLevel(
    int this_level, Duration dt, std::pair<int, int> subcycle) {
  // PreAdvanceLevel might regrid all finer levels.
  // The context must make sure that scratch data is allocated
  PreAdvanceLevel(this_level, dt, subcycle);

  auto scale_dt_on_error = [this](Result<void, TimeStepTooLarge> result) {
    TimeStepTooLarge error = result.error();
    int ratio = GetTotalRefineRatio(error.level);
    error.level = 0;
    error.dt *= ratio;
    return error;
  };

  // If a finer level exists in the hierarchy, we subcycle that finer level
  // multiple times and use the fine fluxes on coarse-fine interfaces
  const int next_level = this_level + 1;
  if (LevelExists(next_level)) {
    ResetCoarseFineFluxes(next_level, this_level);
    const int refine_ratio = GetRatioToCoarserLevel(next_level).max();
    for (int r = 0; r < refine_ratio; ++r) {
      auto result =
          AdvanceLevel(next_level, dt / refine_ratio, {r, refine_ratio});
      if (!result) {
        return scale_dt_on_error(result);
      }
    }
  }

  Result<void, TimeStepTooLarge> result =
      AdvanceLevelNonRecursively(this_level, dt, subcycle);
  if (!result) {
    return scale_dt_on_error(result);
  }
  if (LevelExists(next_level)) {
    ApplyFluxCorrection(next_level, this_level, dt);
  }

  // Coarsen inner regions from next finer level to this level.
  if (LevelExists(next_level)) {
    CoarsenConservatively(next_level, this_level);

    // The conservative update and the coarsening happened on conservative
    // variables. We have to reconstruct the missing variables in the complete
    // state.
    CompleteFromCons(this_level, dt);
  }

  CopyScratchToData(this_level);

  // Apply any further context related work after advancing this level.
  // This function can also indicate if some error occured.
  // For example the context could detect unphysical states and return a
  // TooLargeTimeStep error condition.
  result = PostAdvanceLevel(this_level, dt, subcycle);
  if (!result) {
    return scale_dt_on_error(result);
  }
  return result;
}
```

</div>

The present subcycle solvers depend on a level integrator only.
Future solvers may need additional options that may decsribe varying subcycle ratios between pairs of refinement levels.

Framework for defining FluxMethods {#flux-method-framework}
===========================================================

The `FluxMethod<IntegratorContext>` policy is very generic and allows all kinds of possible logic in user-defined flux methods. 
We provide adapter classes to separate the technical aspects of computing those numerical fluxes from the mathematical ones.
To achieve this separation we begin by introducing the notion of an equation.
An equation gives context to the components of a multi-dimensional array and provides the mathematical dependencies needed.
Using those tools we allow a formulation \( F(q(x, t), dx, dt) \) of flux methods, which depend on a stencil of equation states only.

Equation: Conservative and Complete States {#cons-comp-states}
----------------------------------------------------

The Finite Volume Solver defines so-called equation classes, for example, `Advection2d` or `PerfectGas<Rank>`.
Each equation defines two kinds of states: `Conservative` and `Complete` states.
The conservative states contain variables which see a hyperbolic time update of a time integrator.
The complete state variables are a superset of the conservative ones and may define more auxiliary member variables.
This is reasonable for variables which are needed very often but have a high computation cost.
The complete state is the location to cache those expensive computations.

<div class=example>
For example, the template class `template <int Rank> class PerfectGas` defines the conservative variables

```cpp
template <int Rank> struct PerfectGasConservative {
  double density;
  Array<double, Rank> momentum;
  double energy;
};
```

and the complete state variables are

```cpp
template <int Rank> struct PerfectGasComplete : PerfectGasConservative<Rank> {
  double pressure;
  double speed_of_sound;
};
```
</div>

Given a `MultiFab` object we can iterate over all local `FArrayBox` objects using the `MFIter` iterator as in this simplified example:

```cpp
template <typename Function>
void DoSomethingForEachFArrayBox(amrex::MultiFab& multifab) {
  for (amrex::MFIter mfi(multifab); mfi.isValid(); ++mfi) {
    amrex::FArrayBox& fab = multifab[mfi];
    // do something...
  }
}
```

`amrex::FArrayBox` is a container type in [[AMReX]] that stores a multi-dimensional data array for a single patch.
By using an equation object we can *view* this `amrex::FArrayBox` as a set of multi-dimensional arrays for each variable.

<div class=example>
Lets assume for this example that we are in a context where we are given a `FArrayBox` and we want to access the variables of this patch.
In order to so we need an equation object which interprets the raw data.

```cpp
struct MyClass {
  using Complete = ::fub::Complete<PerfectGas<3>>;
  using Conservative = ::fub::Conservative<PerfectGas<3>>;

  void foo(::amrex::FArrayBox& raw_data) {
    // Transform the raw data object to a View object.
    auto mutable_states = MakeView<Complete>(raw_data, equation);
    auto const_states = MakeView<const Complete>(raw_data, equation);
    auto conservatives = MakeView<Conservative>(raw_data, equation);
    // Iterate through all global grid indices
    ForEachIndex(Box<0>(const_states), [=](int i, int j, int k) {
      // Access the variables
      mutable_states.density(i, j, k) = 1.0;
      conservatives.momentum(i, j, k, 0) = 1.0 + const_states.density(i, j, k);
      conservatives.momentum(i, j, k, 1) = 0.0;
      conservatives.momentum(i, j, k, 2) = 0.0;

      // The objects view the same memory.
      assert(const_states.momentum(i, j, k, 0) == conservatives.momentum(i, j, k, 0));
      assert(const_states.momentum(i, j, k, 1) == conservatives.momentum(i, j, k, 1));
      assert(const_states.momentum(i, j, k, 2) == conservatives.momentum(i, j, k, 2));
    });
  }

  PerfectGas<3> equation;
};
```
</div>

Once we have a `View` object that gives variable-wise access of multi-dimensional data we can load and store state objects, such as `Complete<Equation>` and `Conservative<Equation>`, from and to memory locations inidicated by global space indices `(i, j, k)`.

<div class=example>
We can lift a function that is a defined for state equations, such as `Complete<Equation>`, to a function which is defined on a multi-dimensional array `View<Complete<Equation>>` by applying its logic for each index.

```cpp
struct MyClass {
  using Complete = ::fub::Complete<PerfectGas<3>>;
  using Conservative = ::fub::Conservative<PerfectGas<3>>;

  void Foo(Complete& state) {
    state.density = 1.0;
    state.momentum[0] = 1.0 + state.pressure;
    state.momentum[1] = 0.0;
    state.momentum[2] = 0.0;
  }

  void LiftedFoo(const View<Complete>& view) {
    ForEachIndex(Box<0>(view), [=](int i, int j, int k) {
      Load(state, view, {i, j, k});
      Foo(state);
      Store(view, state, {i, j, k});
    });
  }

  PerfectGas<3> equation;
  Complete state(equation);
};
```

The member function `LiftedFoo` in this example can have multiple equivalent implementations that differ in their efficiency.
Its real implementation is often considered an implementation detail and might change.
</div>

Implement a new FluxMethod, the simple case {#implement-flux-method-simple}
------------------------------------------------------------------------

A class `FM` satisfies the concept `FluxMethod<Equation, StencilSize>` if the following constraints are satisfied.

```cpp
template <typename FM, typename Equation, int StencilSize>
concept FluxMethod = requires (FM& flux_method,
                               Conservative<Equation>& numeric_flux,
                               span<const Complete<Equation>, StencilSize> stencil,
                               double dx, Duration dt, Direction dir) {
  { flux_method.ComputeNumericFlux(numeric_flux, stencil, dx, dt, dir) };
  { flux_method.ComputeStableDt(stencil, dx, dir) } -> Duration;
  { FM::GetStencilSize() } -> std::integral_constant<int, StencilSize>;
};
```

<div class="example">
For example, the following declaration of the class `MusclHancockMethod` satisfies this concept.

```cpp
struct MusclHancockMethod {
    using Conservative = ::Conservative<PerfectGas<2>>;
    using Complete = ::Complete<PerfectGas<2>>;

    void ComputeNumericFlux(Conservative& cons, span<const Complete, 4> stencil,
                              double dx, Duration dt, Direction dir);

    Duration ComputeStableDt(span<const Complete, 4> stencil, double dx, int dir);

    std::integral_constant<int, 2> GetStencilSize() const noexcept { return {}; }
};
```

The implemented MUSCL-type flux methods for an equation `Equation` satisfy `FluxMethod<Equation, 4>` and are implemented in two steps.
First, we reconstruct a state to the left and the right of the face in question.
This gives a reduction of the stencil array from `Complete state[4]` to `Complete reconstruced_states[2]`.
Secondly, we call a (lower-order) base method `BaseFM` which satisfies `FluxMethod<Equation, 2>`.
</div>


Here are some tips for setting up interesting test cases:

  1. If you run long simulation from a remote computer make sure to be able to close the connection without stopping the simulation.
     On Linux one can run docker in a virtual terminal emulator such as `tmux` or `screen` which can be suspended and re-attached without stopping the application.
  
  2. Learn how to use the input files to control output behaviour. If you write your output directly to your host system you can visualize the simulation when it is still running. 
     Beware that this might impact performance on file i/o for windows system if you write data outside of your WSL2 container.
     But if something bad happens you will never lose the output data or checkpoints that are generated up to this point.

  3. Learn how to use checkpoints that restart your simulation from a given time point. 
     This can be either useful to reproduce bugs or to start from special initial conditions that require a stationary flow for example.
     Checkpoints are used by settings the `checkpoint = ''` string in the input file.
     For example, `checkpoint = '/home/amrex/ConvergentNozzle/Checkpoints/000000060'` starts from a checkpoint that is stored in this folder. 
     The number represents the coarse time step number of the simulation.
