#ifndef FUB_AMREX_FARRAYBOX_HPP
#define FUB_AMREX_FARRAYBOX_HPP

#include "fub/Equation.hpp"
#include "fub/State.hpp"

#include <AMReX_BaseFab.H>
#include <AMReX_FArrayBox.H>

namespace fub {
namespace amrex {
/// Creates a mdspan which views all components of a mutable Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
dynamic_mdspan<T, AMREX_SPACEDIM + 1> MakeMdSpan(::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  extents[AMREX_SPACEDIM] = fab.nComp();
  return dynamic_mdspan<T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

/// Creates a mdspan which views the specified component of a mutable Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
dynamic_mdspan<T, AMREX_SPACEDIM> MakeMdSpan(::amrex::BaseFab<T>& fab,
                                             int component) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  return dynamic_mdspan<T, AMREX_SPACEDIM>{fab.dataPtr(component), extents};
}

/// Creates a mdspan which views all components of a const Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
dynamic_mdspan<const T, AMREX_SPACEDIM + 1>
MakeMdSpan(const ::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  extents[AMREX_SPACEDIM] = fab.nComp();
  return dynamic_mdspan<const T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

/// Creates a mdspan which views the specified component of a const Fab.
///
/// \param[in] fab  The Fab which owns the data.
template <typename T>
dynamic_mdspan<const T, AMREX_SPACEDIM>
MakeMdSpan(const ::amrex::BaseFab<T>& fab, int component) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  return dynamic_mdspan<const T, AMREX_SPACEDIM>{fab.dataPtr(component),
                                                 extents};
}

template <typename Equation>
constant<StateType::Complete>
GetStateType(basic_type<Equation>, basic_type<View<Complete<Equation>>>) {
  return {};
}

template <typename Equation>
constant<StateType::Cons> GetStateType(basic_type<Equation>,
                                       basic_type<View<Cons<Equation>>>) {
  return {};
}

template <typename State, typename T, typename Equation>
auto MakeView(boost::hana::basic_type<State>,
              dynamic_mdspan<T, AMREX_SPACEDIM + 1> fab,
              const Equation& equation) {
  auto shape = equation.Shape(GetStateType(type_c<Equation>, type_c<State>));
  int counter = 0;
  auto transform =
      overloaded{[&](int n_comps) {
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e =
                       AsArray(fab.extents());
                   e[AMREX_SPACEDIM] = n_comps;
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
                   index[AMREX_SPACEDIM] = counter;
                   counter += n_comps;
                   return dynamic_mdspan<T, AMREX_SPACEDIM + 1>(&fab(index), e);
                 },
                 [&](std::integral_constant<int, 1>) {
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e1 =
                       AsArray(fab.extents());
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM> e2;
                   std::copy_n(e1.begin(), AMREX_SPACEDIM, e2.begin());
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
                   index[AMREX_SPACEDIM] = counter;
                   counter += 1;
                   return dynamic_mdspan<T, AMREX_SPACEDIM>(&fab(index), e2);
                 }};
  return boost::hana::unpack(
      shape, [&](auto... sizes) { return State{transform(sizes)...}; });
}

template <typename State, typename T, typename Equation>
auto MakeView(boost::hana::basic_type<State>,
              dynamic_mdspan<const T, AMREX_SPACEDIM + 1> fab,
              const Equation& equation) {
  auto shape = equation.Shape(GetStateType(type_c<Equation>, type_c<State>));
  int counter = 0;
  auto transform =
      overloaded{[&](int n_comps) {
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e =
                       AsArray(fab.extents());
                   e[AMREX_SPACEDIM] = n_comps;
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
                   index[AMREX_SPACEDIM] = counter;
                   counter += n_comps;
                   return dynamic_mdspan<const T, AMREX_SPACEDIM + 1>(&fab(index), e);
                 },
                 [&](std::integral_constant<int, 1>) {
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> e1 =
                       AsArray(fab.extents());
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM> e2;
                   std::copy_n(e1.begin(), AMREX_SPACEDIM, e2.begin());
                   std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> index{};
                   index[AMREX_SPACEDIM] = counter;
                   counter += 1;
                   return dynamic_mdspan<const T, AMREX_SPACEDIM>(&fab(index), e2);
                 }};
  using ConstState = decltype(AsConst(std::declval<State>()));
  return boost::hana::unpack(
      shape, [&](auto... sizes) { return ConstState{transform(sizes)...}; });
}

} // namespace amrex
} // namespace fub

#endif
