// Copyright (c) 2021 Maikel Nadolski
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

#ifndef FUB_AMREX_CUTCELL_OUTPUT_WRITE_HDF_IMPL_
#define FUB_AMREX_CUTCELL_OUTPUT_WRITE_HDF_IMPL_

#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/core/span.hpp"

#include <string>

namespace fub::amrex::cutcell {

void WriteHdf5UnRestricted(
    const std::string& name, const PatchHierarchy& hierarchy,
    span<const std::string> fields = span<const std::string>{});

void WriteHdf5RestrictedToBox(
    const std::string& name, const PatchHierarchy& hierarchy,
    const ::amrex::Box& finest_box,
    span<const std::string> fields = span<const std::string>{});

}

#endif