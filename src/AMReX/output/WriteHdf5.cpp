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

#include "fub/AMReX/output/WriteHdf5.hpp"
#include "src/AMReX/output/WriteHdf5Impl.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/ProgramOptions.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

#include <array>

namespace fub::amrex {

WriteHdf5::WriteHdf5(const ProgramOptions& options,
                     std::vector<std::string> field_names)
    : OutputAtFrequencyOrInterval<GriddingAlgorithm>(options),
      field_names_{std::move(field_names)} {
  path_to_file_ = GetOptionOr(options, "path", std::string("grid.h5"));
  fub::SeverityLogger log = GetInfoLogger();
  BOOST_LOG(log) << "WriteHdf5 configured:";
  BOOST_LOG(log) << fmt::format(" - path: {}", path_to_file_);
  auto it = options.find("box");
  if (it != options.end()) {
    auto box = ToMap(it->second);
    std::array<int, AMREX_SPACEDIM> lo =
        GetOptionOr(box, "lower", std::array<int, AMREX_SPACEDIM>{});
    std::array<int, AMREX_SPACEDIM> hi =
        GetOptionOr(box, "upper", std::array<int, AMREX_SPACEDIM>{});
    output_box_.emplace(::amrex::IntVect(lo), ::amrex::IntVect(hi));
    BOOST_LOG(log) << fmt::format(" - box: [{{{}}}, {{{}}}]",
                                  fmt::join(lo, ", "), fmt::join(hi, ", "));
  } else {
    BOOST_LOG(log) << " - box: everything";
  }
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    if (boost::filesystem::is_regular_file(path_to_file_)) {
      boost::log::sources::severity_logger<boost::log::trivial::severity_level>
          log(boost::log::keywords::severity = boost::log::trivial::info);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "HDF5");
      for (int i = 1; i < std::numeric_limits<int>::max(); ++i) {
        std::string new_name = fmt::format("{}.{}", path_to_file_, i);
        if (!boost::filesystem::exists(new_name)) {
          BOOST_LOG(log) << fmt::format(
              "Old output file '{}' detected. Rename old file to '{}'",
              path_to_file_, new_name);
          boost::filesystem::rename(path_to_file_, new_name);
          break;
        }
      }
    } else if (boost::filesystem::exists(path_to_file_)) {
      boost::log::sources::severity_logger<boost::log::trivial::severity_level>
          log(boost::log::keywords::severity = boost::log::trivial::warning);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "HDF5");
      for (int i = 1; i < std::numeric_limits<int>::max(); ++i) {
        std::string new_name = fmt::format("{}.{}", path_to_file_, i);
        if (!boost::filesystem::exists(new_name)) {
          BOOST_LOG(log) << fmt::format(
              "Path'{}' points to some non-file. Output will be directory to "
              "'{}' instead.",
              path_to_file_, new_name);
          path_to_file_ = new_name;
          break;
        }
      }
    }
  }
}

void WriteHdf5::operator()(const GriddingAlgorithm& grid) {
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "HDF5");
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
  BOOST_LOG(log) << fmt::format("Write Hdf5 output to '{}'.", path_to_file_);
  if (output_box_) {
    WriteHdf5RestrictedToBox(path_to_file_, grid.GetPatchHierarchy(),
                             *output_box_, field_names_);
  } else {
    WriteHdf5UnRestricted(path_to_file_, grid.GetPatchHierarchy(),
                          field_names_);
  }
}

} // namespace fub::amrex
