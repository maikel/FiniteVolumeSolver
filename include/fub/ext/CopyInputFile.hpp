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

/// \file
/// \brief This file declares a function to copy the
/// python input script file.

#ifndef FUB_EXT_COPY_INPUTFILE_HPP
#define FUB_EXT_COPY_INPUTFILE_HPP

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "fub/ext/ProgramOptions.hpp"

namespace fub {

struct InputFileOptions {
  InputFileOptions() = default;
  InputFileOptions(const fub::ProgramOptions& options);

  std::string file_template{"input.py"};
  bool copy_input_file{0};
};

void CopyInputFile(MPI_Comm comm, const InputFileOptions& options, int argc,
                   char** argv);

} // namespace fub

#endif // FUB_EXT_COPY_INPUTFILE_HPP