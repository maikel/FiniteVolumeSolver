#include "fub/AMReX.hpp"
#include "fub/Solver.hpp"

struct CircleData {
  void InitializeData(amrex::MultiFab& data,
                      const amrex::Geometry& geom) const {
    fub::amrex::ForEachFab(data, [&](const amrex::MFIter& mfi) {
      const ::amrex::Box& box = mfi.tilebox();
      amrex::FArrayBox& fab = data[mfi];
      fub::amrex::ForEachIndex(box, [&](int i, int j) {
            const double x = geom.CellCenter(i, 0);
            const double y = geom.CellCenter(j, 1);
            const double norm2 = x*x + y*y;
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

  fub::amrex::WritePlotFile("InitialHierarchy/plt00000", gridding->GetPatchHierarchy(), equation);
}