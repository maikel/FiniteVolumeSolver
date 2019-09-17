// Copyright (c) 2018 Maikel Nadolski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "fub/ext/log.hpp"

#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/support/date_time.hpp>

#include "fub/Duration.hpp"

#include <fmt/format.h>
#include <iostream>

namespace fub {

void InitializeLogging(MPI_Comm comm, const LogOptions&) {
  boost::log::core::get()->add_global_attribute("TimeStamp", boost::log::attributes::local_clock());
  int rank = -1;
  MPI_Comm_rank(comm, &rank);
//  if (rank == 0) {
//    // Construct the sink
//    using text_sink = boost::log::sinks::synchronous_sink<boost::log::sinks::text_ostream_backend>;
//    boost::shared_ptr<text_sink> sink = boost::make_shared<text_sink>();
//
//    namespace expr = boost::log::expressions;
//
//    sink->set_formatter(expr::format("%1% -- [%2%]: %3%")
//                        % expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S.%f")
//                        % expr::attr<std::string>("Channel")
//                        % expr::smessage);
//
//    // Add a stream to write log to
//    sink->locked_backend()->add_stream(boost::shared_ptr<std::ostream>(&std::clog, boost::null_deleter()));
//
//    // Register the sink in the logging core
//    boost::log::core::get()->add_sink(sink);
//  }
  namespace expr = boost::log::expressions;
  std::string filename = fmt::format("proccess_{:05}.log", rank);
  boost::log::add_file_log(boost::log::keywords::file_name = filename, boost::log::keywords::format = expr::format("%1% -- [%2%] T = %3%: %4%")
                          % expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S.%f")
                          % expr::attr<std::string>("Channel")
                          % expr::attr<double>("Time")
                          % expr::smessage);
}

}
