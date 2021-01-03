// Copyright (c) 2020 Maikel Nadolski
// Copyright (c) 2020 Stefan Vater
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

#ifndef FUB_AMREX_CUTCELL_DEBUG_OUTPUT_HPP
#define FUB_AMREX_CUTCELL_DEBUG_OUTPUT_HPP

#include "fub/AMReX/output/DebugOutput.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

namespace fub::amrex::cutcell {
class DebugOutput : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
public:
  /// \brief Read program options from opts and enable the storage.
  explicit DebugOutput(const ProgramOptions& opts,
                       const std::shared_ptr<DebugStorage>& storage);

  /// \brief Write out the debug storage on the grid.
  void operator()(const GriddingAlgorithm& grid) override;

private:
  /// \brief This is the base directory where the snapshots will be output to.
  std::string directory_;
};

} // namespace fub::amrex::cutcell

#endif