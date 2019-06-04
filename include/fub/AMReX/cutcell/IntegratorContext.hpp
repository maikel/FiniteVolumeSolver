// Copyright (c) 2019 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

namespace fub::amrex::cutcell {

/// \brief This class manages data and the numerical method to perform cut cell simulations with the AMReX library.
class IntegratorContext {
public:
  static constexpr int Rank = AMREX_SPACEDIM;

  IntegratorContext(GriddingAlgorithm gridding, NumericalMethod method);

  

private:
  struct LevelData {
    LevelData() = default;
    LevelData(const LevelData& other) = delete;
    LevelData& operator=(const LevelData& other) = delete;
    LevelData(LevelData&&) noexcept = default;
    LevelData& operator=(LevelData&&) noexcept;
    ~LevelData() noexcept = default;

    /// This eb_factory is shared with the underlying patch hierarchy.
    std::shared_ptr<::amrex::EBFArrayBoxFactory> eb_factory;

    ///////////////////////////////////////////////////////////////////////////
    // [cell-centered]

    /// reference states which are used to compute embedded boundary fluxes
    std::unique_ptr<::amrex::MultiFab> reference_states;

    /// scratch space filled with data in ghost cells
    std::array<::amrex::MultiFab, Rank> scratch;

    /// fluxes for the embedded boundary
    ::amrex::MultiCutFab boundary_fluxes;

    ///////////////////////////////////////////////////////////////////////////
    // [face-centered]

    /// @{
    /// various flux types needed by the numerical scheme
    std::array<::amrex::MultiFab, Rank> fluxes;
    std::array<::amrex::MultiFab, Rank> stabilized_fluxes;
    std::array<::amrex::MultiFab, Rank> shielded_left_fluxes;
    std::array<::amrex::MultiFab, Rank> shielded_right_fluxes;
    std::array<::amrex::MultiFab, Rank> doubly_shielded_fluxes;
    /// @}

    ///////////////////////////////////////////////////////////////////////////
    // [misc]

    /// FluxRegister accumulate fluxes on coarse fine interfaces between
    /// refinement level. These will need to be rebuilt whenever the hierarchy
    /// changes.
    ::amrex::FluxRegister coarse_fine;

    std::array<Duration, Rank> time_point;
    std::array<Duration, Rank> regrid_time_point;
    std::array<std::ptrdiff_t, Rank> cycles;
  };

  int ghost_cell_width_;
  GriddingAlgorithm gridding_;
  std::vector<LevelData> data_;
  NumericalMethod method_;
};


}
