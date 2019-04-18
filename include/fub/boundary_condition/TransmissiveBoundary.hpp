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

#ifndef FUB_BOUNDARY_CONDITION_TRANSMISSIVE_BOUNDARY_HPP
#define FUB_BOUNDARY_CONDITION_TRANSMISSIVE_BOUNDARY_HPP

#include "fub/Equation.hpp"
#include "fub/State.hpp"
#include "fub/grid/AMReX/PatchHierarchy.hpp"
#include "fub/grid/AMReX/ViewFArrayBox.hpp"

namespace fub {

template <typename Eq> struct TransmissiveBoundary {
  using Equation = Eq;
  using Complete = fub::Complete<Equation>;
  using Conservative = fub::Conservative<Equation>;

  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  explicit TransmissiveBoundary(const Equation& eq) : equation_{eq} {}

  void operator()(const fub::PatchDataView<double, Rank + 1>& data,
                  fub::amrex::PatchHandle, fub::Location location,
                  int fill_width, fub::Duration) {
    fub::BasicView<Complete> complete =
        fub::amrex::MakeView<fub::BasicView<Complete>>(data, equation_);
    fub::IndexBox<Rank + 1> box = data.Box();
    FUB_ASSERT(location.side == 0 || location.side == 1);
    std::array<std::ptrdiff_t, sRank> lower;
    std::copy_n(box.lower.begin(), Rank, lower.begin());
    std::array<std::ptrdiff_t, sRank> upper;
    std::copy_n(box.upper.begin(), Rank, upper.begin());
    if (location.side == 0) {
      upper[location.direction] = lower[location.direction] + fill_width;
      const fub::IndexBox<Rank> fill_box{lower, upper};
      fub::ForEachIndex(fill_box, [&](auto... is) {
        const std::array<std::ptrdiff_t, sRank> dest_index{is...};
        std::array<std::ptrdiff_t, sRank> source_index = dest_index;
        source_index[location.direction] = upper[location.direction];
        Load(state_, AsConst(complete), source_index);
        Store(complete, state_, dest_index);
      });
    } else {
      lower[location.direction] = upper[location.direction] - fill_width;
      const fub::IndexBox<Rank> fill_box{lower, upper};
      fub::ForEachIndex(fill_box, [&](auto... is) {
        const std::array<std::ptrdiff_t, sRank> dest_index{is...};
        std::array<std::ptrdiff_t, sRank> source_index = dest_index;
        source_index[location.direction] =
            upper[location.direction] - fill_width - 1;
        Load(state_, AsConst(complete), source_index);
        Store(complete, state_, dest_index);
      });
    }
  }

  Equation equation_;
  Complete state_{equation_};
};

} // namespace fub

#endif
