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

#include "fub/AMReX/FluxMethod.hpp"

namespace fub::amrex {

FluxMethod::FluxMethod(const FluxMethod& other)
    : flux_method_{other.flux_method_->Clone()} {}

FluxMethod& FluxMethod::operator=(const FluxMethod& other) {
  flux_method_ = other.flux_method_->Clone();
  return *this;
}

double FluxMethod::ComputeStableDt(const ::amrex::FArrayBox& states,
                                   const ::amrex::Box& box,
                                   const ::amrex::Geometry& geom,
                                   Direction dir) {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  return flux_method_->ComputeStableDt(states, box, geom, dir);
}

void FluxMethod::ComputeNumericFluxes(::amrex::FArrayBox& fluxes,
                                      const ::amrex::Box& box,
                                      const ::amrex::FArrayBox& states,
                                      const ::amrex::Geometry& geom,
                                      Duration dt, Direction dir) {
  if (!flux_method_) {
    throw std::runtime_error("Empty flux method.");
  }
  flux_method_->ComputeNumericFluxes(fluxes, box, states, geom, dt, dir);
}

int FluxMethod::GetStencilWidth() const {
  return flux_method_->GetStencilWidth();
}

} // namespace fub::amrex