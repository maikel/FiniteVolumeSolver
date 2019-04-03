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

#include "fub/PatchDataView.hpp"
#include "fub/Equation.hpp"
#include "fub/State.hpp"

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>

namespace fub {
namespace amrex {
std::array<std::ptrdiff_t, AMREX_SPACEDIM> AsArray(const ::amrex::IntVect& vec);

IndexBox<AMREX_SPACEDIM> AsIndexBox(const ::amrex::Box& box);

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
  fub::mdspan<const T, AMREX_SPACEDIM> mdspan = MakeMdSpan(fab, component);
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> lower =
      AsArray(fab.box().smallEnd());
  return PatchDataView<const T, AMREX_SPACEDIM>(mdspan, lower);
}

template <typename Equation>
std::integral_constant<StateType, StateType::Complete>
GetStateType(basic_type<Equation>, basic_type<View<Complete<Equation>>>) {
  return {};
}

template <typename Equation>
std::integral_constant<StateType, StateType::Conservative>
GetStateType(basic_type<Equation>, basic_type<View<Conservative<Equation>>>) {
  return {};
}

template <typename State, typename T, typename Equation>
auto MakeView(boost::hana::basic_type<State>,
              const PatchDataView<T, AMREX_SPACEDIM + 1>& fab,
              const Equation& equation) {
  auto shape = equation.Shape(GetStateType(type_c<Equation>, type_c<State>));
  int counter = 0;
  const std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1>& origin = fab.Origin();
  auto transform = overloaded{
      [&](int n_comps) {
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e =
            AsArray(fab.Extents());
        e[AMREX_SPACEDIM] = n_comps;
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
        index[AMREX_SPACEDIM] = counter;
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> this_origin = origin;
        counter += n_comps;
        mdspan<T, AMREX_SPACEDIM + 1> mds(&fab.MdSpan()(index), e);
        return PatchDataView<T, AMREX_SPACEDIM + 1>(mds, this_origin);
      },
      [&](std::integral_constant<int, 1>) {
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e1 =
            AsArray(fab.Extents());
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> e2;
        std::copy_n(e1.begin(), AMREX_SPACEDIM, e2.begin());
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
        index[AMREX_SPACEDIM] = counter;
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> this_origin;
        std::copy_n(origin.begin(), AMREX_SPACEDIM, this_origin.begin());
        counter += 1;
        mdspan<T, AMREX_SPACEDIM> mds(&fab.MdSpan()(index), e2);
        return PatchDataView<T, AMREX_SPACEDIM>(mds, this_origin);
      }};
  return boost::hana::unpack(
      shape, [&](auto... sizes) { return State{transform(sizes)...}; });
}

template <typename State, typename T, typename Equation>
auto MakeView(boost::hana::basic_type<State>,
              const PatchDataView<const T, AMREX_SPACEDIM + 1>& fab,
              const Equation& equation) {
  auto shape = equation.Shape(GetStateType(type_c<Equation>, type_c<State>));
  int counter = 0;
  const std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1>& origin = fab.Origin();
  auto transform = overloaded{
      [&](int n_comps) {
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e =
            AsArray(fab.Extents());
        e[AMREX_SPACEDIM] = n_comps;
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
        index[AMREX_SPACEDIM] = counter;
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> this_origin = origin;
        counter += n_comps;
        mdspan<const T, AMREX_SPACEDIM + 1> mds(&fab.MdSpan()(index), e);
        return PatchDataView<const T, AMREX_SPACEDIM + 1>(mds, this_origin);
      },
      [&](std::integral_constant<int, 1>) {
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e1 =
            AsArray(fab.Extents());
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> e2;
        std::copy_n(e1.begin(), AMREX_SPACEDIM, e2.begin());
        std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
        index[AMREX_SPACEDIM] = counter;
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> this_origin;
        std::copy_n(origin.begin(), AMREX_SPACEDIM, this_origin.begin());
        counter += 1;
        mdspan<const T, AMREX_SPACEDIM> mds(&fab.MdSpan()(index), e2);
        return PatchDataView<const T, AMREX_SPACEDIM>(mds, this_origin);
      }};
  using ConstState = decltype(AsConst(std::declval<State>()));
  return boost::hana::unpack(
      shape, [&](auto... sizes) { return ConstState{transform(sizes)...}; });
}

template <typename State, typename T, typename Equation>
auto MakeView(const PatchDataView<T, AMREX_SPACEDIM + 1>& fab,
              const Equation& equation) {
  return MakeView(type_c<State>, fab, equation);
}

} // namespace amrex
} // namespace fub

#endif
