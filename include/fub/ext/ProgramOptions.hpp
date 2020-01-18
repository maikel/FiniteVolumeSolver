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

/// \file This file declares a function to parse program options from a python
/// script file.

#ifndef FUB_EXT_PROGRAM_OPTIONS_HPP
#define FUB_EXT_PROGRAM_OPTIONS_HPP

#include <map>
#include <string>
#include <optional>

#include <boost/filesystem.hpp>
#include <mpi.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace fub {

using ProgramOptions = std::map<std::string, pybind11::object>;

ProgramOptions ParsePythonScript(const boost::filesystem::path& path,
                                 MPI_Comm comm);

ProgramOptions ToMap(const pybind11::dict& dict);

template <typename T>
T GetOptionOr(const ProgramOptions& map, const std::string& name,
              const T& value) {
  auto iter = map.find(name);
  if (iter != map.end()) {
    return iter->second.cast<T>();
  }
  return value;
}

ProgramOptions GetOptions(const ProgramOptions& options, const std::string& name);

std::optional<ProgramOptions> ParseCommandLine(int argc, char** argv);

} // namespace fub

#endif