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

#include "fub/ext/ProgramOptions.hpp"
#include "fub/ext/Mpi.hpp"

#include <fmt/format.h>
#include <fstream>
#include <string>

namespace fub {

std::map<std::string, pybind11::object>
ParsePythonScript(const boost::filesystem::path& path, MPI_Comm comm) {
  if (!boost::filesystem::is_regular_file(path)) {
    throw std::runtime_error(
        fmt::format("Path '{}' is not a regular file", path.string()));
  }
  std::string content = ReadAndBroadcastFile(path.string(), comm);
  pybind11::exec(content.c_str());
  std::map<std::string, pybind11::object> options;
  for (const auto& [key, value] : pybind11::globals()) {
    const auto name = key.cast<std::string>();
    if (name.substr(0, 2) != "__") {
      options.emplace(name, pybind11::cast<pybind11::object>(value));
    }
  }
  return options;
}

std::map<std::string, pybind11::object> ToMap(const pybind11::dict& dict) {
  std::map<std::string, pybind11::object> options;
  for (const auto& [k, v] : dict) {
    const auto key = k.cast<std::string>();
    options.emplace(key, pybind11::cast<pybind11::object>(v));
  }
  return options;
}

} // namespace fub