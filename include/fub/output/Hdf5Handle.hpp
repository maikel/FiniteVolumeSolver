// Copyright (c) 2020 Maikel Nadolski
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

#include <hdf5.h>

namespace fub {

template <typename Deleter> struct H5Handle {
  H5Handle() = default;
  H5Handle(hid_t id) : id_{id} {}

  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;

  H5Handle(H5Handle&& other) noexcept : id_{std::exchange(other.id_, -1)} {}
  H5Handle& operator=(H5Handle&& other) {
    id_ = std::exchange(other.id_, -1);
    return *this;
  }

  ~H5Handle() {
    if (id_ > 0) {
      Deleter{}(id_);
    }
  }

  operator hid_t() const noexcept { return id_; }

  hid_t id_{};
};

struct H5Fdeleter {
  using pointer = hid_t;
  void operator()(hid_t file) const noexcept { H5Fclose(file); }
};
using H5File = H5Handle<H5Fdeleter>;

struct H5Sdeleter {
  using pointer = hid_t;
  void operator()(hid_t dataspace) const noexcept { H5Sclose(dataspace); }
};
using H5Space = H5Handle<H5Sdeleter>;

struct H5Ddeleter {
  using pointer = hid_t;
  void operator()(hid_t dataset) const noexcept { H5Dclose(dataset); }
};
using H5Dataset = H5Handle<H5Ddeleter>;

struct H5Adeleter {
  using pointer = hid_t;
  void operator()(hid_t attributes) const noexcept { H5Aclose(attributes); }
};
using H5Attribute = H5Handle<H5Adeleter>;

struct H5Tdeleter {
  using pointer = hid_t;
  void operator()(hid_t attributes) const noexcept { H5Tclose(attributes); }
};
using H5Type = H5Handle<H5Tdeleter>;

struct H5Pdeleter {
  using pointer = hid_t;
  void operator()(hid_t properties) const noexcept { H5Pclose(properties); }
};
using H5Properties = H5Handle<H5Pdeleter>;

} // namespace fub