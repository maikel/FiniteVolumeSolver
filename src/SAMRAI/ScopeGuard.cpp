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

#include "fub/SAMRAI/ScopeGuard.hpp"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

namespace fub {

/// This method uses the global communication group, MPI_WORLD_COMM, for calls
/// to SAMRAI routines.
ScopeGuard::ScopeGuard(int argc, char** argv) {
  ::SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
  ::SAMRAI::tbox::SAMRAIManager::initialize();
  ::SAMRAI::tbox::SAMRAIManager::startup();
}

/// Calls SAMRAI routines to do the clean-up.
ScopeGuard::~ScopeGuard() noexcept {
  ::SAMRAI::tbox::SAMRAIManager::shutdown();
  ::SAMRAI::tbox::SAMRAIManager::finalize();
  ::SAMRAI::tbox::SAMRAI_MPI::finalize();
}

} // namespace fub