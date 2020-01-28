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

#ifndef FUB_AMREX_CUTCELL_PLOT_FILES_HPP
#define FUB_AMREX_CUTCELL_PLOT_FILES_HPP

#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex::cutcell {

template <typename Equation>
class PlotfileOutput : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
public:
  PlotfileOutput(const Equation& equation, const std::string& path)
      : OutputAtFrequencyOrInterval(), equation_(equation), parent_path_(path) {
  }

  PlotfileOutput(std::vector<std::ptrdiff_t> freqs,
                 std::vector<Duration> intervals, const Equation& equation,
                 const std::string& path)
      : OutputAtFrequencyOrInterval(std::move(freqs), std::move(intervals)),
        equation_(equation), parent_path_(path) {}

  void operator()(const GriddingAlgorithm& grid) override {
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
        log(boost::log::keywords::severity = boost::log::trivial::info);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "Plotfile");
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
    std::string name =
        fmt::format("{}/plt{:09}", parent_path_, grid.GetCycles());
    BOOST_LOG(log) << fmt::format("Write Plotfile output to '{}'.", name);
    WritePlotFile(name, grid.GetPatchHierarchy(), equation_);
  }

private:
  Equation equation_;
  std::string parent_path_;
};

} // namespace fub::amrex::cutcell

#endif