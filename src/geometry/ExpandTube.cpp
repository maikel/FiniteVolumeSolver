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

#include "fub/core/assert.hpp"
#include "fub/geometry/ExpandTube.hpp"
#include <fstream>
#include <sstream>
#include <vector>

namespace fub {
namespace {
Eigen::Matrix<double, 3, Eigen::Dynamic>
GetHeightsFromPoints_(const Eigen::Matrix<double, 3, Eigen::Dynamic>& points) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> heights(3, points.cols() - 1);
  for (int c = 0; c < points.cols() - 1; ++c) {
    heights.col(c) = (points.col(c + 1) - points.col(c)).normalized();
  }
  return heights;
}

Eigen::Matrix<double, 3, Eigen::Dynamic>
GetPlanes(const Eigen::Matrix<double, 3, Eigen::Dynamic>& center_points,
          const Eigen::Matrix<double, 3, Eigen::Dynamic>& heights) {
  Eigen::Matrix<double, 3, Eigen::Dynamic> planes(3, center_points.cols());
  Eigen::Vector3d norm0{0.0, 1.0, 0.0};// = heights.col(0);
  planes(0, 0) = norm0[0];
  planes(1, 0) = norm0[1];
  planes(2, 0) = norm0.dot(center_points.col(0));
  FUB_ASSERT(heights.cols() == center_points.cols() - 1);
  for (int c = 1; c < heights.cols(); ++c) {
    Eigen::Vector3d norm1 = heights.col(c);
    Eigen::Vector3d avg = (0.5 * (norm0 + norm1)).normalized();
    planes(0, c) = avg[0];
    planes(1, c) = avg[1];
    planes(2, c) = avg.dot(center_points.col(c));
    norm0 = norm1;
  }
  const std::ptrdiff_t last = center_points.cols() - 1;
  norm0 = Eigen::Vector3d{0.0, -1.0, 0.0};
  planes(0, last) = norm0[0];
  planes(1, last) = norm0[1];
  planes(2, last) = norm0.dot(center_points.col(last));
  return planes;
}

double ComputeDistance_(const Eigen::Vector3d& x, const Eigen::Vector3d& origin,
                        const Eigen::Vector3d& h, double radius) {
  const double length = h.dot(x - origin);
  Eigen::Vector3d projected = origin + length * h;
  const double r = (x - projected).norm();
  return r - radius;
}

std::string ReadFile_(const std::string& filename) {
  std::ifstream ifs(filename);
  std::string file_content;
  ifs.seekg(0, std::ios::end);
  std::ptrdiff_t size = ifs.tellg();
  if (size < 0) {
    throw std::runtime_error("Invalid File '" + filename + "'!");
  }
  file_content.reserve(static_cast<std::size_t>(size));
  ifs.seekg(0, std::ios::beg);
  file_content.assign(std::istreambuf_iterator<char>(ifs),
                      std::istreambuf_iterator<char>());
  return file_content;
}
} // namespace

ExpandTube::ExpandTube(
    const Eigen::Matrix<double, 3, Eigen::Dynamic>& center_points,
    double radius)
    : radius_{radius}, center_points_{center_points},
      heights_{GetHeightsFromPoints_(center_points_)}, planes_{GetPlanes(
                                                           center_points_,
                                                           heights_)} {}

double ExpandTube::operator()(double x, double y, double z) const {
  Eigen::Vector2d X(x, y);
  const auto num_cylinder = heights_.cols();
  std::vector<double> distances(static_cast<std::size_t>(num_cylinder));
  for (int i = 0; i < num_cylinder; ++i) {
    Eigen::Vector2d normlo{planes_(0, i), planes_(1, i)};
    Eigen::Vector2d normhi{planes_(0, i + 1), planes_(1, i + 1)};
    distances[static_cast<std::size_t>(i)] = std::max(
        {planes_(2, i) - normlo.dot(X), normhi.dot(X) - planes_(2, i + 1),
          ComputeDistance_({x, y, z}, center_points_.col(i),
                          heights_.col(i), radius_)});
  }
  return -*std::min_element(distances.begin(), distances.end());
}

Eigen::Matrix<double, 3, Eigen::Dynamic>
ReadPointsFromFile(const std::string& filename) {
  const std::string content = ReadFile_(filename);
  const std::ptrdiff_t n_points =
      std::count(content.begin(), content.end(), '\n');
  Eigen::Matrix<double, 3, Eigen::Dynamic> points(3, n_points);
  std::istringstream iss(content);
  std::string line{};
  int col = 0;
  while (std::getline(iss, line) && col < n_points) {
    std::istringstream line_in(line);
    double x;
    if (line_in >> x) {
      points(0, col) = x;
    }
    if (line_in >> x) {
      points(1, col) = x;
    }
    points(2, col) = 0.0;
    if (line_in) {
      ++col;
    }
  }
  points.resize(3, col);
  return points;
}

} // namespace fub
