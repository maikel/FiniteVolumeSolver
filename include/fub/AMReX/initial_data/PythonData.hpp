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

#ifndef FUB_AMREX_INITIAL_DATA_PYTHON_DATA_HPP
#define FUB_AMREX_INITIAL_DATA_PYTHON_DATA_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"

#include "fub/ext/ProgramOptions.hpp"

namespace fub::amrex {

template <typename Equation> class PythonData {
public:
  PythonData(Equation eq, const ProgramOptions& options);

  void InitializeData(PatchLevel& patch_level, const GriddingAlgorithm& grid,
                      int level, Duration time);

private:
  Equation equation_;
  std::function<void(Complete<Equation>&, std::array<double, 3>, Equation&)>
      initial_data_fn_;
};

template <typename Equation>
PythonData<Equation>::PythonData(Equation eq, const ProgramOptions& options)
    : equation_{std::move(eq)} {
  auto py_function = [&]() -> std::optional<pybind11::function> {
    if (auto iter = options.find("function"); iter != options.end()) {
      pybind11::object obj = iter->second;
      pybind11::function fun(obj);
      return fun;
    }
    return {};
  }();

  initial_data_fn_ = [py_function](Complete<Equation>& state,
                                   std::array<double, 3> x,
                                   Equation& equation) {
    if (py_function) {
      pybind11::function function = *py_function;
      pybind11::object ret = function(x);
      auto [T, p, moles, u] = ret.cast<
          std::tuple<double, double, std::string, std::array<double, 3>>>();
      equation.GetReactor().SetMoleFractions(moles);
      equation.GetReactor().SetTemperature(T);
      equation.GetReactor().SetPressure(p);
      Array<double, Equation::Rank(), 1> velocity;
      velocity.setZero();
      for (int d = 0; d < Equation::Rank(); ++d) {
        velocity[d] = u[d];
      }
      equation.CompleteFromReactor(state, velocity);
    } else {
      equation.GetReactor().SetMoleFractions("O2:21,N2:79");
      equation.GetReactor().SetTemperature(300.);
      equation.GetReactor().SetPressure(101325.0);
      Array<double, Equation::Rank(), 1> velocity;
      velocity.setZero();
      equation.CompleteFromReactor(state, velocity);
    }
  };
}

template <typename Equation>
void PythonData<Equation>::InitializeData(PatchLevel& patch_level,
                                          const GriddingAlgorithm& grid,
                                          int level, Duration /*time*/) {
  ::amrex::MultiFab& mf = patch_level.data;
  const ::amrex::Geometry& geom = grid.GetPatchHierarchy().GetGeometry(level);
  Complete<Equation> state(equation_);
  ForEachFab(execution::seq, mf, [&](const ::amrex::MFIter& mfi) {
    ::amrex::FArrayBox& fab = mf[mfi];
    auto view = MakeView<Complete<Equation>>(fab, equation_, mfi.tilebox());
    ForEachIndex(Box<0>(view), [&](auto... is) {
      std::array<double, 3> x{};
      geom.CellCenter(::amrex::IntVect(int(is)...), x.data());
      initial_data_fn_(state, x, equation_);
      Store(view, state, {is...});
    });
  });
}

} // namespace fub::amrex

#endif
