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

#ifndef FUB_EXT_BOOST_LOG_HPP
#define FUB_EXT_BOOST_LOG_HPP

#include "fub/Duration.hpp"
#include <boost/log/common.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/program_options.hpp>
#include <mpi.h>

namespace fub {

struct LogOptions {
  LogOptions() = default;
  LogOptions(const boost::program_options::variables_map& vm);

  static boost::program_options::options_description GetCommandLineOptions();

  std::string file_template{"{rank:04d}.log"};
  std::vector<int> which_mpi_ranks_do_log{0};
};

void InitializeLogging(MPI_Comm comm, const LogOptions& log = {});


using SeverityLogger =
    boost::log::sources::severity_logger<boost::log::trivial::severity_level>;

inline SeverityLogger GetInfoLogger() {
  return SeverityLogger(boost::log::keywords::severity =
                            boost::log::trivial::info);
}

void Log(std::string message, Duration timepoint,
         boost::log::trivial::severity_level level =
             boost::log::trivial::severity_level::info);

} // namespace fub

#endif // FINITEVOLUMESOLVER_LOG_HPP
