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

#ifndef FUB_AMREX_CUTCELL_LOG_PROBES_HPP
#define FUB_AMREX_CUTCELL_LOG_PROBES_HPP

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

#include <map>
#include <string>
#include <vector>

namespace fub::amrex {

class LogProbesOutput
    : public OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm> {
public:
  LogProbesOutput(const std::map<std::string, pybind11::object>& vm);

  void operator()(const MultiBlockGriddingAlgorithm& grid) override;

private:
  std::string plenum_output_path_;
  std::vector<double> plenum_probes_;

  std::string tube_output_path_;
  std::vector<double> tube_probes_;
};

} // namespace fub::amrex

#endif