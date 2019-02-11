#ifndef FUB_AMREX_FARRAYBOX_HPP
#define FUB_AMREX_FARRAYBOX_HPP

#include "fub/StateFacade.hpp"

#include <AMReX_BaseFab.h>
#include <AMReX_FArrayBox.h>

namespace fub {
namespace amrex {

template <typename T> auto GetBox(T&& state) {
  constexpr auto as = std::decay_t<T>::accessors();
  return boost::hana::second(boost::hana::at_c<0>(as))(state)->box();
}

struct ScalarTag {};
struct VectorTag {};

template <typename T>
dynamic_mdspan<T, AMREX_SPACEDIM + 1> MakeMdSpan(VectorTag,
                                                 ::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM + 1> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  extents[AMREX_SPACEDIM] = fab.nComps();
  return dynamic_mdspan<T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

template <typename T>
dynamic_mdspan<T, AMREX_SPACEDIM> MakeMdSpan(ScalarTag,
                                             ::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  return dynamic_mdspan<T, AMREX_SPACEDIM>{fab.dataPtr(), extents};
}

template <typename T> auto MakeMdSpan(::amrex::BaseFab<T>& fab) {
  return MakeMdSpan(ScalarTag{}, fab);
}

template <typename T>
dynamic_mdspan<const T, AMREX_SPACEDIM + 1>
MakeMdSpan(VectorTag{}, const ::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  extents[AMREX_SPACEDIM] = fab.nComps();
  return dynamic_mdspan<const T, AMREX_SPACEDIM + 1>{fab.dataPtr(), extents};
}

template <typename T>
dynamic_mdspan<const T, AMREX_SPACEDIM>
MakeMdSpan(ScalarTag{}, const ::amrex::BaseFab<T>& fab) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> extents;
  ::amrex::IntVect length{fab.box().length()};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[i] = length[i];
  }
  return dynamic_mdspan<const T, AMREX_SPACEDIM>{fab.dataPtr(), extents};
}

template <typename T> auto MakeMdSpan(const ::amrex::BaseFab<T>& fab) {
  return MakeMdSpan(ScalarTag{}, fab);
}

template <typename State> auto ViewDataOnFAB(State&& state) {
  using T = std::decay_t<State>;
  constexpr auto sizes =
      GetStaticSizes<typename std::decay_t<T>::Equation::State>();
  constexpr auto accessors = T::accessors();
  constexpr auto zipped = boost::hana::zip(accessors, sizes.as_tuple());
  return boost::hana::unpack(zipped, [&](auto... xs) {
    auto as_mdspan = [&](auto x) {
      auto size = [](auto x) { return boost::hana::at_c<1>(x); };
      auto array = [&](auto x) {
        return boost::hana::second(boost::hana::at_c<0>(x))(state);
      };
      if constexpr (std::decay_t<decltype(size(x))>() == 1) {
        return MakeMdSpan(ScalarTag{}, *array(x));
      } else {
        return MakeMdSpan(VectorTag{}, *array(x));
      }
    };
    return boost::hana::make<typename T::hana_tag>(as_mdspan(xs)...);
  });
}

} // namespace amrex
} // namespace fub

#endif