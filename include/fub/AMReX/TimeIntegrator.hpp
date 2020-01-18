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

#ifndef FUB_AMREX_HYPERBOLIC_SPLIT_TIME_INTEGRATOR_HPP
#define FUB_AMREX_HYPERBOLIC_SPLIT_TIME_INTEGRATOR_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace fub::amrex {

template <typename Tag> struct ForwardIntegrator {
  ForwardIntegrator() = default;
  explicit ForwardIntegrator(Tag) {}

  void UpdateConservatively(::amrex::MultiFab& dest,
                            const ::amrex::MultiFab& src,
                            const ::amrex::MultiFab& fluxes,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir);

  void UpdateConservatively(IntegratorContext& context, int level, Duration dt,
                            Direction dir);
};

ForwardIntegrator() -> ForwardIntegrator<execution::OpenMpSimdTag>;

extern template struct ForwardIntegrator<execution::SequentialTag>;
extern template struct ForwardIntegrator<execution::OpenMpTag>;
extern template struct ForwardIntegrator<execution::SimdTag>;
extern template struct ForwardIntegrator<execution::OpenMpSimdTag>;

inline ForwardIntegrator<execution::OpenMpSimdTag> EulerForwardTimeIntegrator() {
  return {};
}

template <typename Tag>
ForwardIntegrator<Tag> EulerForwardTimeIntegrator(Tag) {
  return {};
}

} // namespace fub::amrex

#endif
