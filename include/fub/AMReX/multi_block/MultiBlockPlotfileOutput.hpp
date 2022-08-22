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

#ifndef FUB_AMREX_MULTI_BLOCK_PLOT_FILE_OUTPUT_HPP
#define FUB_AMREX_MULTI_BLOCK_PLOT_FILE_OUTPUT_HPP

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"
#include "fub/equations/EulerEquation.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

#include <string>

namespace fub::amrex {

class MultiBlockPlotfileOutput
    : public OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm> {
public:
  MultiBlockPlotfileOutput(const std::map<std::string, pybind11::object>& vm);

  void operator()(const MultiBlockGriddingAlgorithm& grid) override;

private:
  std::string parent_path_;
};

template <typename TubeEquation, typename PlenumEquation>
class MultiBlockPlotfileOutput2
    : public OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
public:
  MultiBlockPlotfileOutput2(const std::map<std::string, pybind11::object>& vm,
                            const TubeEquation& tube_equation,
                            const PlenumEquation& plenum_equation);

  void operator()(const MultiBlockGriddingAlgorithm2& grid) override;

private:
  TubeEquation tube_equation_;
  PlenumEquation plenum_equation_;
  std::string parent_path_;
};

template <typename TubeEquation, typename PlenumEquation>
MultiBlockPlotfileOutput2<TubeEquation, PlenumEquation>::
    MultiBlockPlotfileOutput2(const std::map<std::string, pybind11::object>& vm,
                              const TubeEquation& tube_equation,
                              const PlenumEquation& plenum_equation)
    : OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2>(vm),
      tube_equation_(tube_equation), plenum_equation_(plenum_equation) {
  parent_path_ = GetOptionOr(vm, "directory", std::string("."));
}

template <typename TubeEquation, typename PlenumEquation>
void MultiBlockPlotfileOutput2<TubeEquation, PlenumEquation>::
operator()(const MultiBlockGriddingAlgorithm2& grid) {
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "Plotfile");
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
  for (int i = 0; i < grid.GetPlena().size(); ++i) {
    std::string name =
        fmt::format("{}/Plenum{}/plt{:09}", parent_path_, i, grid.GetCycles());
    BOOST_LOG(log) << fmt::format("Write Plotfile output to '{}'.", name);
    cutcell::WritePlotFile(name, grid.GetPlena()[i]->GetPatchHierarchy(),
                           plenum_equation_);
  }
  for (int i = 0; i < grid.GetTubes().size(); ++i) {
    std::string name =
        fmt::format("{}/Tube{}/plt{:09}", parent_path_, i, grid.GetCycles());
    BOOST_LOG(log) << fmt::format("Write Plotfile output to '{}'.", name);
    WritePlotFile(name, grid.GetTubes()[i]->GetPatchHierarchy(),
                  tube_equation_);
  }
}

} // namespace fub::amrex

#endif