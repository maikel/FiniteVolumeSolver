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
#include "fub/ext/ProgramOptions.hpp"

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_RealBox.H>

namespace fub {
template <>
::amrex::IntVect GetOptionOr(const ProgramOptions& map, const std::string& name,
                             const ::amrex::IntVect& value);

template <>
::amrex::Box GetOptionOr(const ProgramOptions& map, const std::string& name,
                         const ::amrex::Box& value);

template <>
::amrex::RealBox GetOptionOr(const ProgramOptions& map, const std::string& name,
                             const ::amrex::RealBox& value);

namespace amrex {

std::array<std::ptrdiff_t, AMREX_SPACEDIM> AsArray(const ::amrex::IntVect& vec);

template <int Rank> IndexBox<Rank> AsIndexBox(const ::amrex::Box& box) {
  const std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(box.smallEnd());
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> upper = AsArray(box.bigEnd());
  std::transform(upper.begin(), upper.end(), upper.begin(),
                 [](std::ptrdiff_t i) { return i + 1; });

  if constexpr (AMREX_SPACEDIM != Rank) {
    constexpr std::size_t sRank = static_cast<std::size_t>(Rank);
    std::array<std::ptrdiff_t, sRank> lower1;
    std::array<std::ptrdiff_t, sRank> upper1;
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
  Index<AMREX_SPACEDIM + 1> extents{};
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
  Index<AMREX_SPACEDIM> lower = AsArray(fab.box().smallEnd());
  Index<AMREX_SPACEDIM + 1> lower_comp{};
  std::copy_n(lower.begin(), AMREX_SPACEDIM, lower_comp.begin());
  return PatchDataView<T, AMREX_SPACEDIM + 1>(mdspan, lower_comp);
}

/// Creates a mdspan which views the specified component of a mutable Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
mdspan<T, AMREX_SPACEDIM> MakeMdSpan(::amrex::BaseFab<T>& fab, int component) {
  Index<AMREX_SPACEDIM> extents{};
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
  Index<AMREX_SPACEDIM> lower = AsArray(fab.box().smallEnd());
  return PatchDataView<T, AMREX_SPACEDIM>(mdspan, lower);
}

/// Creates a mdspan which views all components of a const Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
mdspan<const T, AMREX_SPACEDIM + 1> MakeMdSpan(const ::amrex::BaseFab<T>& fab) {
  Index<AMREX_SPACEDIM + 1> extents{};
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[static_cast<std::size_t>(i)] = length[i];
  }
  extents[static_cast<std::size_t>(AMREX_SPACEDIM)] = fab.nComp();
  return mdspan<const T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

template <typename T>
PatchDataView<const T, AMREX_SPACEDIM + 1>
MakePatchDataView(const ::amrex::BaseFab<T>& fab) {
  fub::mdspan<const T, AMREX_SPACEDIM + 1> mdspan = MakeMdSpan(fab);
  Index<AMREX_SPACEDIM> lower = AsArray(fab.box().smallEnd());
  Index<AMREX_SPACEDIM + 1> lower_comp;
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
  Index<AMREX_SPACEDIM> extents{};
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
  Index<AMREX_SPACEDIM> lower = AsArray(fab.box().smallEnd());
  return PatchDataView<const T, AMREX_SPACEDIM>(mdspan, lower);
}

template <typename T>
PatchDataView<const T, AMREX_SPACEDIM, layout_stride>
MakePatchDataView(const ::amrex::BaseFab<T>& fab, int component,
                  const ::amrex::Box& box) {
  mdspan<const T, AMREX_SPACEDIM> mdspan = MakeMdSpan(fab, component);
  Index<AMREX_SPACEDIM> lower = AsArray(fab.box().smallEnd());
  return PatchDataView<const T, AMREX_SPACEDIM>(mdspan, lower)
      .Subview(AsIndexBox<AMREX_SPACEDIM>(box));
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
    const Index<AMREX_SPACEDIM + 1>& origin = fab.Origin();
    auto transform = overloaded{
        [&](int n_comps) {
          Index<AMREX_SPACEDIM + 1> efab = AsArray(fab.Extents());
          Index<Rank + 1> e{};
          std::copy_n(efab.begin(), Rank, e.begin());
          e[sRank] = n_comps;
          Index<AMREX_SPACEDIM + 1> index{};
          index[static_cast<std::size_t>(AMREX_SPACEDIM)] = counter;
          Index<Rank + 1> this_origin{};
          std::copy_n(origin.begin(), Rank, this_origin.begin());
          // this_origin[sRank] = counter;
          counter += n_comps;
          mdspan<ValueType, sRank + 1> mds(&fab.MdSpan()(index), e);
          return PatchDataView<ValueType, Rank + 1>(mds, this_origin);
        },
        [&](const ScalarDepth&) {
          Index<AMREX_SPACEDIM + 1> efab = AsArray(fab.Extents());
          Index<Rank> e{};
          std::copy_n(efab.begin(), Rank, e.begin());
          Index<AMREX_SPACEDIM + 1> index{};
          index[AMREX_SPACEDIM] = counter;
          Index<Rank> this_origin{};
          std::copy_n(origin.begin(), Rank, this_origin.begin());
          counter += 1;
          mdspan<ValueType, sRank> mds(&fab.MdSpan()(index), e);
          return PatchDataView<ValueType, Rank>(mds, this_origin);
        }};
    State pd_views{};
    ForEachVariable(
        [&](auto& pdv, const auto& depth) { pdv = transform(depth); }, pd_views,
        depths);
    return pd_views;
  }
};

template <typename State, typename T, typename Equation>
auto MakeView(const PatchDataView<T, AMREX_SPACEDIM + 1>& fab,
              const Equation& equation) {
  MakeViewImpl<State> make_view{};
  return make_view(fab, equation);
}

template <typename State, typename Equation>
auto MakeView(::amrex::FArrayBox& fab, const Equation& equation) {
  return MakeView<BasicView<State>>(MakePatchDataView(fab), equation);
}

template <typename State, typename Equation>
auto MakeView(const ::amrex::FArrayBox& fab, const Equation& equation) {
  return MakeView<BasicView<State>>(MakePatchDataView(fab), equation);
}

template <typename State, typename Equation>
auto MakeView(::amrex::FArrayBox& fab, const Equation& eq,
              const IndexBox<Equation::Rank()>& box) {
  return Subview(MakeView<State>(fab, eq), box);
}

template <typename State, typename Equation>
auto MakeView(const ::amrex::FArrayBox& fab, const Equation& eq,
              const IndexBox<Equation::Rank()>& box) {
  return Subview(MakeView<State>(fab, eq), box);
}

template <typename State, typename FAB, typename Equation>
auto MakeView(FAB&& fab, const Equation& eq, const ::amrex::Box& box) {
  return MakeView<State>(std::forward<FAB>(fab), eq,
                         AsIndexBox<Equation::Rank()>(box));
}

std::array<::amrex::Box, 2>
GetCellsAndFacesInStencilRange(const ::amrex::Box& cell_tilebox,
                               const ::amrex::Box& face_validbox,
                               int stencil_width, Direction dir);

} // namespace amrex
} // namespace fub

#endif
