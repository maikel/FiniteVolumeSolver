// Copyright (c) 2021 Christian Zenker
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

#include "fub/ext/CopyInputFile.hpp"

namespace fub {

InputFileOptions::InputFileOptions(const fub::ProgramOptions& options) {
  FUB_GET_OPTION_VAR(options, file_template);
  FUB_GET_OPTION_VAR(options, copy_input_file);
}

void CopyInputFile(MPI_Comm comm, const InputFileOptions& options, int argc,
                   char** argv) {
  if (!options.copy_input_file) {
    return;
  }

  int rank = -1;
  MPI_Comm_rank(comm, &rank);

  std::string filename = options.file_template;
  if (rank == 0) {
    // create directory if neccessary
    boost::filesystem::path path(filename);
    boost::filesystem::path dir = boost::filesystem::absolute(
        path.parent_path(), boost::filesystem::current_path());
    if (!boost::filesystem::exists(dir)) {
      boost::filesystem::create_directories(dir);
    }
  }
  namespace po = boost::program_options;
  po::options_description desc{};
  std::string config_path{};
  desc.add_options()("config", po::value<std::string>(&config_path),
                     "Path to the config file which can be parsed.");
  desc.add_options()("args",
                     po::value<std::vector<std::string>>()->multitoken(),
                     "Arguments for the input file");
  desc.add_options()("help", "Print this help message.");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("config")) {
    config_path = vm["config"].as<std::string>();
    boost::filesystem::copy_file(
        config_path, filename,
        boost::filesystem::copy_option::overwrite_if_exists);
  }
}

} // namespace fub