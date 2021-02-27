// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
// Copyright (c) 2021 Rupert Klein
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

#ifndef FUB_EQUATIONS_PERFECT_GAS_MIX_PLENA_CONTROL_OUTPUT_HPP
#define FUB_EQUATIONS_PERFECT_GAS_MIX_PLENA_CONTROL_OUTPUT_HPP

#include "fub/AMReX/multi_block/MultiBlockIntegratorContext2.hpp"
#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

namespace fub::perfect_gas_mix::gt {

class ControlOutput
    : public OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2> {
public:
  ControlOutput(const ProgramOptions& options,
                std::shared_ptr<const ControlState> control_state);

  void operator()(const amrex::MultiBlockGriddingAlgorithm2& grid) override;

private:
  std::string file_path_;
  std::shared_ptr<const ControlState> control_state_;
  std::map<std::string, function_ref<double(const ControlState&)>> fields_;
  std::vector<double> data_buffer_;

  void CreateHdf5Database();
  void WriteHdf5Database(span<const double> data, Duration time,
                         std::ptrdiff_t cycle);
};

} // namespace fub::perfect_gas_mix::gt

#endif