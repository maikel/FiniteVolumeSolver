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

#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>

namespace fub {
namespace amrex {

template <typename Function, typename MultiFab, typename... MultiFabs>
void ForEachFAB(Function function, MultiFab&& fab, MultiFabs&&... fabs) {
  const ::amrex::Vector<int>& local_indices = fab.IndexArray();
  for (auto i : local_indices) {
    function(fab[i], fabs[i]...);
  }
}

template <typename F>
::amrex::PhysBCFunct<F>
MakePhysBCFunct(const ::amrex::Geometry& geom,
                const ::amrex::Vector<::amrex::BCRec>& bcr, F f) {
  return {geom, bcr, f};
}

} // namespace amrex
} // namespace fub