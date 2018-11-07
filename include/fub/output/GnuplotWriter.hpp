// Copyright (c) 2018 Maikel Nadolski
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

#include "fub/core/span.hpp"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <string>
#include <vector>

namespace fub {
enum class VariableType { scalar, vector };

/// This class creates output files which can be read by gnuplot.
class GnuplotWriter {
public:
  /// Constructs the class with a directory which will store the output files.
  explicit GnuplotWriter(std::string directory,
                         span<const int> ids = span<const int>());

  /// Adds a specified patch data to the list of output variables.
  void registerPlotQuantity(std::string name, VariableType type,
                            int patch_data_id);

  /// Writes all quantitites which has been registered with this class to the
  /// output directory.
  void
  writePlotData(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                int cylce, double dt);

private:
  std::unique_ptr<std::ostream> out_;
  std::string directory_;
  std::vector<std::string> names_{};
  std::vector<VariableType> types_{};
  std::vector<int> patch_data_ids_{};
};

} // namespace fub