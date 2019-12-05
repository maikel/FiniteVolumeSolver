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

#include <p4est.h>

#include "fub/core/function_ref.hpp"

#include <mpi.h>

#include <memory>

namespace fub::p4est {
struct P4estDeleter {
  void operator()(p4est_t* pointer) { p4est_destroy(pointer); }
};

class P4est {
public:
  using RefineFn = function_ref<int(const P4est&, int, const Quadrant&)>;
  using CoarsenFn = function_ref<int(const P4est&, int, span<Quadrant, 4>)>;
  using ReplaceFn = function_ref<void(const P4est&)>;
  using WeightFn = function_ref<void(const P4est&, const Quadrant&)>;

  P4est() = delete;

  P4est(const P4est&);
  P4est& operator=(const P4est&);

  P4est(P4est&&) = default;
  P4est& operator=(P4est&&) = default;

  P4est(MPI_Comm mpi_communicator, Connectivity connectivity);

  P4est(MPI_Comm mpi_communicator, Connectivity connectivity, int min_levels);

  void Refine(RefineFn tag_quadrants, ReplaceFn on_replace);

  void Coarsen(CoarsenFn tag_quadrants, ReplaceFn on_replace);

  void Balance(ReplaceFn on_replace);

  void Partition(WeightFn weight);

  operator p4est_t*() noexcept { return pointer_.get(); }
  operator const p4est_t*() const noexcept { return pointer_.get(); }

private:
  std::unique_ptr<p4est_t, P4estDeleter> pointer_;
};

} // namespace fub::p4est