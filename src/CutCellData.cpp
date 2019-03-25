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

#include "fub/CutCellData.hpp"

namespace fub {

Eigen::Vector2d GetBoundaryNormal(const CutCellData<2>& ccdata,
                                  const std::array<std::ptrdiff_t, 2>& index) {
  Eigen::Vector2d normal;
  normal[0] = -ccdata.boundary_normals(index[0], index[1], 0);
  normal[1] = -ccdata.boundary_normals(index[0], index[1], 1);
  return normal;
}

Eigen::Vector3d GetBoundaryNormal(const CutCellData<3>& ccdata,
                                  const std::array<std::ptrdiff_t, 3>& index) {
  Eigen::Vector3d normal;
  normal[0] = -ccdata.boundary_normals(index[0], index[1], index[2], 0);
  normal[1] = -ccdata.boundary_normals(index[0], index[1], index[2], 1);
  normal[2] = -ccdata.boundary_normals(index[0], index[1], index[2], 2);
  return normal;
}

} // namespace fub