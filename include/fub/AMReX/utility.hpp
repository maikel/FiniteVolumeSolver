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

#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>

namespace fub {
namespace amrex {

template <typename State, typename Equation>
auto ToMultiFabs(const Equation& eq) {
  return Transform([](auto&&) { return ::amrex::MultiFab{}; }, State());
}

template <typename Equation> auto ToFluxRegisters(const Equation& eq) {
  return Transform(
      [](auto&&) { return std::unique_ptr<::amrex::FluxRegister>{}; },
      typename Equation::Cons());
}

template <typename StateT>
using MultiFabs = decltype(
    ToMultiFabs<StateT>(std::declval<const typename StateT::Equation&>()));

template <typename ConsT>
using CoarseFineFluxes =
    decltype(ToFluxRegisters(std::declval<const typename ConsT::Equation&>()));

template <typename T, std::size_t N>
::amrex::Vector<T> ToAmrexVector(const std::array<T, N>& array) {
  return ::amrex::Vector<T>(array.begin(), array.end());
}

template <typename T> auto Size(const StateFacade<T>& facade) {
  using namespace boost::hana::literals;
  return boost::hana::second(facade.accessors()[0_c])(facade).size();
}

template <typename T> auto Size(const ::amrex::FabArray<T>& arrays) {
  return arrays.size();
}

template <typename T> auto LocalSize(const StateFacade<T>& facade) {
  using namespace boost::hana::literals;
  return boost::hana::second(facade.accessors()[0_c])(facade).local_size();
}

template <typename T> auto LocalSize(const ::amrex::FabArray<T>& arrays) {
  return arrays.local_size();
}


template <typename T> auto LocalIndices(const StateFacade<T>& facade) {
  using namespace boost::hana::literals;
  return boost::hana::second(facade.accessors()[0_c])(facade).IndexArray();
}

template <typename T> auto LocalIndices(const ::amrex::FabArray<T>& arrays) {
  return arrays.IndexArray();
}

template <typename Function, typename MultiFab, typename... MultiFabs>
void ForEachFAB(Function function, MultiFab&& fab, MultiFabs&&... fabs) {
  auto get_fab_ptrs = [&](auto&& state, int i) -> decltype(auto) {
    using T = std::decay_t<decltype(state)>;
    if constexpr (std::is_base_of<::amrex::FabArrayBase, T>::value) {
      return state[i];
    } else {
      return Transform([i](auto&& multi_fab) { return &multi_fab[i]; }, state);
    }
  };
  const ::amrex::Vector<int>& local_indices = LocalIndices(fab);
  for (auto i : local_indices) {
    function(get_fab_ptrs(fab, i), get_fab_ptrs(fabs, i)...);
  }
}

} // namespace amrex
} // namespace fub