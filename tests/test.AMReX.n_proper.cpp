#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_ParmParse.H>

#include <sstream>
#include <stdexcept>

struct AmrCoreTest : amrex::AmrCore {
  AmrCoreTest(const amrex::Box& coarse_fill_box, int n_proper_in,
              const amrex::RealBox& rb, int max_level_in,
              const amrex::Vector<int>& n_cell_in, ::amrex::IntVect ref_ratio0)
  : AmrCore(&rb, max_level_in, n_cell_in, -1, amrex::Vector<amrex::IntVect>(max_level_in + 1, ref_ratio0)),
        coarse_fill_box_{coarse_fill_box}, box_arrays_(max_level_in + 1),
        distributions_(max_level_in + 1) {
    amrex::AmrMesh::verbose = 1;
    n_proper = n_proper_in;
    ref_ratio = ::amrex::Vector<::amrex::IntVect>(max_level_in + 1, ref_ratio0);
  }

  amrex::Box GetRefinedFillBox(int level) {
    amrex::Box fill_box = coarse_fill_box_;
    for (int lvl = 0; lvl < level; ++lvl) {
      fill_box.refine(refRatio(lvl));
    }
    return fill_box;
  }

  // Refine a fixed coarse index range over all refinement levels.
  void ErrorEst(int level, ::amrex::TagBoxArray& tags, double,
                int /* ngrow */) override {
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
      const amrex::Box fill_box = GetRefinedFillBox(level) & mfi.tilebox();
      const amrex::Vector<int> data(fill_box.numPts(), 1);
      tags[mfi].tags(data, fill_box);
    }
  }

  // Reset internal box_arrays and distribution mappings.
  void MakeNewLevelFromScratch(
      int level, double /* time_point */, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override {
    box_arrays_[level] = box_array;
    distributions_[level] = distribution_mapping;
  }

  // Reset internal box_arrays and distribution mappings.
  void MakeNewLevelFromCoarse(
      int level, double /* time_point */, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override {
    box_arrays_[level] = box_array;
    distributions_[level] = distribution_mapping;
  }

  // Reset internal box_arrays and distribution mappings.
  void RemakeLevel(
      int level, double /* time_point */, const ::amrex::BoxArray& box_array,
      const ::amrex::DistributionMapping& distribution_mapping) override {
    box_arrays_[level] = box_array;
    distributions_[level] = distribution_mapping;
  }

  // If this function is called the test failed as we always tag an index range.
  void ClearLevel(int level) override {
    std::ostringstream oss{};
    oss << "Test failed: ref_ratio = " << refRatio(level)
        << ", n_proper = " << n_proper << ", level = " << level << '\n';
    throw std::logic_error{oss.str()};
  }

  amrex::Box coarse_fill_box_;
  std::vector<amrex::BoxArray> box_arrays_;
  std::vector<amrex::DistributionMapping> distributions_;
};

void MyMain() {
  amrex::ParmParse pp("amr");
  pp.add("blocking_factor_x", 8);
  pp.add("blocking_factor_y", 1);
  pp.add("blocking_factor_z", 1);

  const amrex::Box coarse_fill_box{{AMREX_D_DECL(15, 0, 0)}, {AMREX_D_DECL(16, 0, 0)}};
  const amrex::Vector<int> n_cell_in{AMREX_D_DECL(32, 1, 1)};
  const amrex::RealBox rb{{AMREX_D_DECL(-1.0, -1.0, -1.0)},
                          {AMREX_D_DECL(+1.0, +1.0, +1.0)}};
  const int max_level_in = 2;

  const ::amrex::IntVect ref_ratio_0{AMREX_D_DECL(2, 1, 1)};
  const int n_proper_0 = 0;
  AmrCoreTest amr_core_0(coarse_fill_box, n_proper_0, rb, max_level_in,
                         n_cell_in, ref_ratio_0);
  amr_core_0.MakeNewGrids(0.0);
  amr_core_0.regrid(1, 0.0);

  const ::amrex::IntVect ref_ratio_1{AMREX_D_DECL(2, 2, 2)};
  const int n_proper_1 = 1;
  AmrCoreTest amr_core_1(coarse_fill_box, n_proper_1, rb, max_level_in,
                         n_cell_in, ref_ratio_1);
  amr_core_1.MakeNewGrids(0.0);
  amr_core_1.regrid(1, 0.0);

  const ::amrex::IntVect ref_ratio_2{AMREX_D_DECL(2, 1, 1)};
  const int n_proper_2 = 1;
  AmrCoreTest amr_core_2(coarse_fill_box, n_proper_2, rb, max_level_in,
                         n_cell_in, ref_ratio_2);
  amr_core_2.MakeNewGrids(0.0);
  amr_core_2.regrid(1, 0.0);
}

int main(int argc, char** argv) {
  amrex::Initialize(argc, argv);
  try {
    MyMain();
  } catch (...) {
    amrex::Finalize();
    throw;
  }
  amrex::Finalize();
}
