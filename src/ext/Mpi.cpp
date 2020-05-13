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

#include "fub/ext/Mpi.hpp"
#include "fub/core/assert.hpp"

#include <fmt/format.h>
#include <fstream>

namespace fub {

std::string ReadAndBroadcastFile(std::string filepath, MPI_Comm comm) {
  int rank = -1;
  MPI_Comm_rank(comm, &rank);
  std::string buffer{};
  int size = -1;
  if (rank == 0) {
    std::ifstream file(filepath);
    if (!file) {
      throw std::runtime_error(
          fmt::format("Could not open file '{}'.", filepath));
    }
    file.seekg(0, std::ios::end);
    size = static_cast<int>(file.tellg());
  }
  MPI_Bcast(&size, 1, MPI_INT, 0, comm);
  FUB_ASSERT(size >= 0);
  const auto ssize = static_cast<std::size_t>(size);
  buffer.resize(ssize);
  if (rank == 0) {
    std::ifstream file(filepath);
    file.read(buffer.data(), size);
  }
  MPI_Bcast(buffer.data(), size, MPI_CHAR, 0, comm);
  return buffer;
}

Duration MinAll(MPI_Comm comm, Duration local) {
  double local_count = local.count();
  double global_count = 0.0;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_DOUBLE, MPI_MIN, comm);
  return Duration(global_count);
}

} // namespace fub
