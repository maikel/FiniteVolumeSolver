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

#ifndef FUB_GEOMETRY_EXPAND_TUBE_HPP
#define FUB_GEOMETRY_EXPAND_TUBE_HPP

#include <Eigen/Dense>

namespace fub {

class ExpandTube {
public:
  ExpandTube(const Eigen::Matrix<double, 3, Eigen::Dynamic>& center_points,
             double radius);

  double operator()(double x, double y, double z) const;

  template <typename Vector> double operator()(const Vector& x) const {
    return this->operator()(x[0], x[1], x[2]);
  }

private:
  double radius_;
  Eigen::Matrix<double, 3, Eigen::Dynamic> center_points_;
  Eigen::Matrix<double, 3, Eigen::Dynamic> heights_;
  Eigen::Matrix<double, 3, Eigen::Dynamic> planes_;
};

Eigen::Matrix<double, 3, Eigen::Dynamic>
ReadPointsFromFile(const std::string& filename);

} // namespace fub

#endif
