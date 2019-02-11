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

#ifndef FUB_AMREX_PATCH_LEVEL_HPP
#define FUB_AMREX_PATCH_LEVEL_HPP

namespace fub {

struct DataDescription {
  int n_state_components;
  int n_flux_components;
  int num_ghosts_cells;
  ::amrex::BoxArray box_array;
  ::amrex::DistributionMapping distribution_mapping;
};

class PatchLevel {
public:
  PatchLevel(const DataDescrption& description, int level_number,
             double time_point);

  int GetLevelNumber() const noexcept { return level_number_; }

  double GetTimePoint() const noexcept { return time_point_; }

  span<::amrex::MultiFab> old_time() noexcept { return old_time_level_; }
  span<const ::amrex::MultiFab> old_time() const noexcept {
    return old_time_level_;
  }

  span<::amrex::MultiFab> new_time() noexcept { return new_time_level_; }
  span<const ::amrex::MultiFab> new_time() const noexcept {
    return new_time_level_;
  }

  span<::amrex::MultiFab> scratch(int direction) noexcept {
    return scratch_[direction];
  }
  span<const ::amrex::MultiFab> scratch(int direction) const noexcept {
    return scratch_[direction];
  }

  span<::amrex::MultiFab> fluxes(int direction) noexcept {
    return fluxes_[direction];
  }
  span<const ::amrex::MultiFab> fluxes(int direction) const noexcept {
    return fluxes_[direction];
  }

  span<::amrex::FluxRegister> coarse_fine() noexcept;
  span<const ::amrex::FluxRegister> coarse_fine() const noexcept;

private:
  int level_number_;
  double time_point_;
  std::vector<::amrex::MultiFab> states_;
};

} // namespace fub

#endif