// Copyright (c) 2018 Maikel Nadolski
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

#include "fub/ext/Log.hpp"

#include <boost/algorithm/string/replace.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>

#include "fub/Duration.hpp"

#include <fmt/format.h>
#include <fstream>

namespace fub {

namespace {
BOOST_LOG_ATTRIBUTE_KEYWORD(a_timestamp, "TimeStamp",
                            boost::log::attributes::local_clock::value_type)

void FormatLogs_(const boost::log::record_view& rec,
                 boost::log::formatting_ostream& stream) {
  namespace log = boost::log;
  namespace expr = log::expressions;
  std::ostringstream prefix;
  prefix << rec[a_timestamp] << " -- ";
  auto channel = log::extract<std::string>("Channel", rec);
  if (channel) {
    prefix << '[' << channel.get() << "] ";
  }
  auto sev = log::extract<log::trivial::severity_level>("Severity", rec);
  if (sev) {
    prefix << '[' << sev.get() << "] ";
  }
  auto time = log::extract<double>("Time", rec);
  if (time) {
    prefix << fmt::format("[T = {:>12.6g}s] ", time.get());
  }
  auto level = log::extract<int>("Level", rec);
  if (level) {
    prefix << fmt::format("[level = {}] ", level.get());
  }

  const std::string replace_by = std::string("\n") + prefix.str();
  std::string message = rec[expr::smessage].get();
  if (message.size() > 1) {
    boost::replace_all(message, "\n", replace_by);
    stream << prefix.str() << message;
  }
}
} // namespace

void InitializeLogging(MPI_Comm comm, const LogOptions& options) {
  boost::log::core::get()->add_global_attribute(
      "TimeStamp", boost::log::attributes::local_clock());
  int rank = -1;
  MPI_Comm_rank(comm, &rank);
  auto first = options.which_mpi_ranks_do_log.begin();
  auto last = options.which_mpi_ranks_do_log.end();
  auto found = std::find(first, last, rank);
  if (found != last) {
    using text_sink = boost::log::sinks::synchronous_sink<
        boost::log::sinks::text_ostream_backend>;
    boost::shared_ptr<text_sink> file = boost::make_shared<text_sink>();
    namespace expr = boost::log::expressions;
    file->locked_backend()->add_stream(boost::make_shared<std::ofstream>(
        fmt::format(options.file_template, fmt::arg("rank", rank))));
    file->set_formatter(&FormatLogs_);
    boost::log::core::get()->add_sink(file);
    boost::shared_ptr<text_sink> console = boost::make_shared<text_sink>();
    boost::shared_ptr<std::ostream> cout(&std::cout, boost::null_deleter{});
    console->locked_backend()->add_stream(cout);
    console->set_formatter(&FormatLogs_);
    console->set_filter(boost::log::trivial::severity >=
                        boost::log::trivial::info);
    boost::log::core::get()->add_sink(console);
  } else {
    boost::log::core::get()->set_filter([](auto&&) { return false; });
  }
}

} // namespace fub
