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

#ifndef FUB_AMREX_FLUX_METHOD_HPP
#define FUB_AMREX_FLUX_METHOD_HPP

#include "fub/AMReX/PatchHandle.hpp"
#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/State.hpp"

#include <AMReX_Geometry.H>

namespace fub::amrex {

struct FluxMethodBase {
  virtual ~FluxMethodBase() = default;

  virtual std::unique_ptr<FluxMethodBase> Clone() const = 0;

  virtual double ComputeStableDt(const ::amrex::FArrayBox& states,
                                 const ::amrex::Box& box,
                                 const ::amrex::Geometry& geom,
                                 Direction dir) = 0;

  virtual void ComputeNumericFluxes(::amrex::FArrayBox& fluxes,
                                    const ::amrex::Box& box,
                                    const ::amrex::FArrayBox& states,
                                    const ::amrex::Geometry& geom, Duration dt,
                                    Direction dir) = 0;

  virtual int GetStencilWidth() const = 0;
};

template <typename FM> struct FluxMethodWrapper : FluxMethodBase {
  FluxMethodWrapper() = default;

  FluxMethodWrapper(const FM& fm);
  FluxMethodWrapper(FM&& fm) noexcept;

  std::unique_ptr<FluxMethodBase> Clone() const override;

  double ComputeStableDt(const ::amrex::FArrayBox& fab, const ::amrex::Box& box,
                         const ::amrex::Geometry& geom, Direction dir) override;

  void ComputeNumericFluxes(::amrex::FArrayBox& fluxes, const ::amrex::Box& box,
                            const ::amrex::FArrayBox& states,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir) override;

  int GetStencilWidth() const override;

  FM flux_method_{};
};

class FluxMethod {
public:
  FluxMethod() = delete;

  template <typename FM, typename = std::enable_if_t<
                             !std::is_same_v<FluxMethod, std::decay_t<FM>>>>
  FluxMethod(FM&& flux_method)
      : flux_method_{std::make_unique<FluxMethodWrapper<std::decay_t<FM>>>(
            std::forward<FM>(flux_method))} {}

  FluxMethod(const FluxMethod& other);
  FluxMethod& operator=(const FluxMethod&);

  FluxMethod(FluxMethod&&) = default;
  FluxMethod& operator=(FluxMethod&&) = default;

  double ComputeStableDt(const ::amrex::FArrayBox& states,
                         const ::amrex::Box& box, const ::amrex::Geometry& geom,
                         Direction dir);

  void ComputeNumericFluxes(::amrex::FArrayBox& fluxes, const ::amrex::Box& box,
                            const ::amrex::FArrayBox& states,
                            const ::amrex::Geometry& geom, Duration dt,
                            Direction dir);

  int GetStencilWidth() const;

private:
  std::unique_ptr<FluxMethodBase> flux_method_;
};

template <typename FM>
FluxMethodWrapper<FM>::FluxMethodWrapper(const FM& fm) : flux_method_{fm} {}

template <typename FM>
FluxMethodWrapper<FM>::FluxMethodWrapper(FM&& fm) noexcept
    : flux_method_{std::move(fm)} {}

template <typename FM>
std::unique_ptr<FluxMethodBase> FluxMethodWrapper<FM>::Clone() const {
  return std::make_unique<FluxMethodWrapper<FM>>(flux_method_);
}

template <typename FM>
double FluxMethodWrapper<FM>::ComputeStableDt(const ::amrex::FArrayBox& fab,
                                              const ::amrex::Box& box,
                                              const ::amrex::Geometry& geom,
                                              Direction dir) {
  using Eq = std::decay_t<decltype(flux_method_.GetEquation())>;
  static constexpr int Rank = Eq::Rank();
  Eq& equation = flux_method_.GetEquation();
  const double dx = geom.CellSize(int(dir));
  const IndexBox<Rank> cells = AsIndexBox<Rank>(box);
  View<const Complete<Eq>> states =
      MakeView<const Complete<Eq>>(fab, equation, cells);
  return flux_method_.ComputeStableDt(states, dx, dir);
}

template <typename FM>
void FluxMethodWrapper<FM>::ComputeNumericFluxes(
    ::amrex::FArrayBox& fluxes, const ::amrex::Box& box,
    const ::amrex::FArrayBox& states, const ::amrex::Geometry& geom,
    Duration dt, Direction dir) {
  using Eq = std::decay_t<decltype(flux_method_.GetEquation())>;
  static const int Rank = Eq::Rank();
  Eq& equation = flux_method_.GetEquation();
  const double dx = geom.CellSize(int(dir));

  const int gcw = flux_method_.GetStencilWidth() + 1;

  const IndexBox<Rank> cells = Grow(AsIndexBox<Rank>(box), dir, {gcw, gcw});
  View<const Complete<Eq>> state =
      MakeView<const Complete<Eq>>(states, equation, cells);

  const IndexBox<Rank> faces = Grow(
      AsIndexBox<Rank>(::amrex::surroundingNodes(box, int(dir))), dir, {1, 1});
  View<Conservative<Eq>> flux =
      MakeView<Conservative<Eq>>(fluxes, equation, faces);

  flux_method_.ComputeNumericFluxes(flux, state, dt, dx, dir);
}

template <typename FM> int FluxMethodWrapper<FM>::GetStencilWidth() const {
  return flux_method_.GetStencilWidth();
}
} // namespace fub::amrex

#endif
