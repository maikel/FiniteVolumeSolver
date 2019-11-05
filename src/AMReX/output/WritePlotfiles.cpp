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

#include "fub/AMReX/output/WritePlotfiles.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex {

  MultiBlockPlotfileOutput::MultiBlockPlotfileOutput(const std::map<std::string, pybind11::object>& vm)
  : OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm>(vm) {
    parent_path_ = GetOptionOr(vm, "directory", std::string("."));
  }

  void
  MultiBlockPlotfileOutput::operator()(const MultiBlockGriddingAlgorithm& grid) {
    Burke2012 burke2012;
    IdealGasMix<AMREX_SPACEDIM> plenum_eq(burke2012);
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>
          log(boost::log::keywords::severity = boost::log::trivial::info);
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "Plotfile");
    BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
    for (int i = 0; i < grid.GetPlena().size(); ++i) {
      std::string name = fmt::format("{}/Plenum{}/plt{:09}", parent_path_, i, grid.GetCycles());
      BOOST_LOG(log) << fmt::format("Write Plotfile output to '{}'.", name);
      cutcell::WritePlotFile(name, grid.GetPlena()[i]->GetPatchHierarchy(), plenum_eq);
    }
    IdealGasMix<1> tube_eq(burke2012);
    for (int i = 0; i < grid.GetTubes().size(); ++i) {
      std::string name = fmt::format("{}/Tube{}/plt{:09}", parent_path_, i, grid.GetCycles());
      BOOST_LOG(log) << fmt::format("Write Plotfile output to '{}'.", name);
      WritePlotFile(name, grid.GetTubes()[i]->GetPatchHierarchy(), tube_eq);
    }
  }


}
