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

#include "fub/ext/ProgramOptions.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"
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

}

#endif