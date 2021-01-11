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

#include "fub/AMReX/cutcell/output/DebugOutput.hpp"

namespace fub::amrex::cutcell {
DebugOutput::DebugOutput(const ProgramOptions& opts,
                         const std::shared_ptr<DebugStorage>& storage)
    : OutputAtFrequencyOrInterval(opts) {
  OutputAtFrequencyOrInterval::frequencies_ = std::vector<std::ptrdiff_t>{1LL};
  directory_ = GetOptionOr(opts, "directory", directory_);
  storage->Enable();
  SeverityLogger log = GetInfoLogger();
  BOOST_LOG(log) << "DebugOutput configured:";
  BOOST_LOG(log) << fmt::format(" - directory = '{}'", directory_);
  OutputAtFrequencyOrInterval::Print(log);
}

void DebugOutput::operator()(const GriddingAlgorithm& grid) {
  DebugStorage& storage = *grid.GetPatchHierarchy().GetDebugStorage();
  const std::ptrdiff_t cycles = grid.GetPatchHierarchy().GetCycles();
  [[maybe_unused]] const auto int_max =
      static_cast<std::ptrdiff_t>(std::numeric_limits<int>::max());
  FUB_ASSERT(cycles < int_max);
  storage.FlushData(directory_, static_cast<int>(cycles), grid.GetTimePoint());
}

} // namespace fub::amrex::cutcell
