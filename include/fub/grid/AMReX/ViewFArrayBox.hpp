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

#ifndef FUB_AMREX_FARRAYBOX_HPP
#define FUB_AMREX_FARRAYBOX_HPP

#include "fub/Equation.hpp"
#include "fub/PatchDataView.hpp"
#include "fub/State.hpp"

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>

namespace fub {
namespace amrex {
std::array<std::ptrdiff_t, AMREX_SPACEDIM> AsArray(const ::amrex::IntVect& vec);

template <int Rank> IndexBox<Rank> AsIndexBox(const ::amrex::Box& box) {
  const std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(box.smallEnd());
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> upper = AsArray(box.bigEnd());
  std::transform(upper.begin(), upper.end(), upper.begin(),
                 [](std::ptrdiff_t i) { return i + 1; });

  if constexpr (AMREX_SPACEDIM != Rank) {
    std::array<std::ptrdiff_t, Rank> lower1;
    std::array<std::ptrdiff_t, Rank> upper1;
    std::copy_n(lower.begin(), Rank, lower1.begin());
    std::copy_n(upper.begin(), Rank, upper1.begin());
    return IndexBox<Rank>{lower1, upper1};
  } else {
    return IndexBox<Rank>{lower, upper};
  }
}

/// Creates a mdspan which views all components of a mutable Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
mdspan<T, AMREX_SPACEDIM + 1> MakeMdSpan(::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[static_cast<std::size_t>(i)] = length[i];
  }
  extents[AMREX_SPACEDIM] = fab.nComp();
  return mdspan<T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

template <typename T>
PatchDataView<T, AMREX_SPACEDIM + 1>
MakePatchDataView(::amrex::BaseFab<T>& fab) {
  fub::mdspan<T, AMREX_SPACEDIM + 1> mdspan = MakeMdSpan(fab);
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(fab.box().smallEnd());
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> lower_comp;
  std::copy_n(lower.begin(), AMREX_SPACEDIM, lower_comp.begin());
  lower_comp[AMREX_SPACEDIM] = 0;
  return PatchDataView<T, AMREX_SPACEDIM + 1>(mdspan, lower_comp);
}

/// Creates a mdspan which views the specified component of a mutable Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
mdspan<T, AMREX_SPACEDIM> MakeMdSpan(::amrex::BaseFab<T>& fab, int component) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[static_cast<std::size_t>(i)] = length[i];
  }
  return mdspan<T, AMREX_SPACEDIM>{fab.dataPtr(component), extents};
}

template <typename T>
PatchDataView<T, AMREX_SPACEDIM> MakePatchDataView(::amrex::BaseFab<T>& fab,
                                                   int component) {
  fub::mdspan<T, AMREX_SPACEDIM> mdspan = MakeMdSpan(fab, component);
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(fab.box().smallEnd());
  return PatchDataView<T, AMREX_SPACEDIM>(mdspan, lower);
}

/// Creates a mdspan which views all components of a const Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
mdspan<const T, AMREX_SPACEDIM + 1> MakeMdSpan(const ::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[static_cast<std::size_t>(i)] = length[i];
  }
  extents[AMREX_SPACEDIM] = fab.nComp();
  return mdspan<const T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

template <typename T>
PatchDataView<const T, AMREX_SPACEDIM + 1>
MakePatchDataView(const ::amrex::BaseFab<T>& fab) {
  fub::mdspan<const T, AMREX_SPACEDIM + 1> mdspan = MakeMdSpan(fab);
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(fab.box().smallEnd());
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> lower_comp;
  std::copy_n(lower.begin(), AMREX_SPACEDIM, lower_comp.begin());
  lower_comp[AMREX_SPACEDIM] = 0;
  return PatchDataView<const T, AMREX_SPACEDIM + 1>(mdspan, lower_comp);
}

/// Creates a mdspan which views the specified component of a const Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
mdspan<const T, AMREX_SPACEDIM> MakeMdSpan(const ::amrex::BaseFab<T>& fab,
                                           int component) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[static_cast<std::size_t>(i)] = length[i];
  }
  return mdspan<const T, AMREX_SPACEDIM>{fab.dataPtr(component), extents};
}

template <typename T>
PatchDataView<const T, AMREX_SPACEDIM>
MakePatchDataView(const ::amrex::BaseFab<T>& fab, int component) {
  mdspan<const T, AMREX_SPACEDIM> mdspan = MakeMdSpan(fab, component);
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(fab.box().smallEnd());
  return PatchDataView<const T, AMREX_SPACEDIM>(mdspan, lower);
}

template <typename State> struct MakeViewImpl {
  using Equation = typename State::Equation;
  using Depths = typename State::Depths;
  using ValueType = typename State::ValueType;

  static constexpr int Rank = Equation::Rank();
  static constexpr std::size_t sRank = static_cast<std::size_t>(Rank);

  State operator()(const PatchDataView<ValueType, AMREX_SPACEDIM + 1>& fab,
                   const Equation& equation) {
    const auto depths = ::fub::Depths<State, Equation>(equation);
    int counter = 0;
    const std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1>& origin = fab.Origin();
    auto transform = overloaded{
        [&](int n_comps) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> efab =
              AsArray(fab.Extents());
          std::array<std::ptrdiff_t, sRank + 1> e;
          std::copy_n(efab.begin(), Rank, e.begin());
          e[sRank] = n_comps;
          std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
          index[AMREX_SPACEDIM] = counter;
          std::array<std::ptrdiff_t, sRank + 1> this_origin{};
          std::copy_n(origin.begin(), Rank, this_origin.begin());
          this_origin[sRank] = counter;
          counter += n_comps;
          mdspan<ValueType, sRank + 1> mds(&fab.MdSpan()(index), e);
          return PatchDataView<ValueType, Rank + 1>(mds, this_origin);
        },
        [&](const ScalarDepth&) {
          std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> efab =
              AsArray(fab.Extents());
          std::array<std::ptrdiff_t, sRank> e;
          std::copy_n(efab.begin(), Rank, e.begin());
          std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
          index[AMREX_SPACEDIM] = counter;
          std::array<std::ptrdiff_t, sRank> this_origin;
          std::copy_n(origin.begin(), Rank, this_origin.begin());
          counter += 1;
          mdspan<ValueType, sRank> mds(&fab.MdSpan()(index), e);
          return PatchDataView<ValueType, Rank>(mds, this_origin);
        }};
    State pd_views{};
    ForEachVariable<State>(
        [&](auto& pdv, const auto& depth) { pdv = transform(depth); }, pd_views,
        depths);
    return pd_views;
  }
};

template <typename State, typename T, typename Equation>
auto MakeView(const PatchDataView<T, AMREX_SPACEDIM + 1>& fab,
              const Equation& equation) {
  MakeViewImpl<State> make_view;
  return make_view(fab, equation);
}

} // namespace amrex
} // namespace fub

#endif
