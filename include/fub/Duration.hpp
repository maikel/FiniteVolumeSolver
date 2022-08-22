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

#ifndef FUB_DURATION_HPP
#define FUB_DURATION_HPP

#include <chrono>
#include <cmath>
#include <ostream>

#include <boost/serialization/access.hpp>

namespace fub {

struct Duration : std::chrono::duration<double> {
  using std::chrono::duration<double>::duration;

  Duration(const std::chrono::duration<double>& other) noexcept : duration(other) {}
};

inline std::ostream& operator<<(std::ostream& out, Duration dur) {
  return out << dur.count() << " [s]";
} 

}

namespace boost::serialization {

template <typename Archive>
void serialize(Archive& ar, ::fub::Duration& t, unsigned int) {
  double count = t.count();
  if (!std::isfinite(count)) {
    count = std::numeric_limits<double>::lowest();
  }
  ar& count;
  t = ::fub::Duration(count);
}

} // namespace boost::serialization

#endif
