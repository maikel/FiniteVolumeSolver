// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_AMREX_MULTI_BLOCK_BOUNDARY2_HPP
#define FUB_AMREX_MULTI_BLOCK_BOUNDARY2_HPP

#include "fub/AMReX/IntegratorContext.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"
#include "fub/AMReX/multi_block/MultiBlockBoundary.hpp"
#include "fub/equations/EulerEquation.hpp"

#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/PatchDataView.hpp"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <vector>

namespace fub::amrex {

class MultiBlockGriddingAlgorithm2;

/// \ingroup BoundaryCondition
///
struct MultiBlockBoundaryBase {
  virtual ~MultiBlockBoundaryBase() = default;

  virtual void FillTubeGhostLayer(
      ::amrex::FArrayBox& tube_ghost_data,
      const ::amrex::FArrayBox& plenum_mirror_data,
      const Eigen::Matrix<double, AMREX_SPACEDIM, 1>& normal) = 0;
  virtual void
  FillPlenumGhostLayer(::amrex::FArrayBox& plenum_ghost_data,
                       const ::amrex::FArrayBox& tube_mirror_data) = 0;

  virtual std::unique_ptr<MultiBlockBoundaryBase> Clone() const = 0;
};

template <typename Boundary>
struct MultiBlockBoundaryWrapper : public MultiBlockBoundaryBase {
  MultiBlockBoundaryWrapper(const Boundary& impl) : impl_(impl) {}
  MultiBlockBoundaryWrapper(Boundary&& impl) : impl_(std::move(impl)) {}

  std::unique_ptr<MultiBlockBoundaryBase> Clone() const override final {
    return std::make_unique<MultiBlockBoundaryWrapper>(impl_);
  }

  void FillTubeGhostLayer(
      ::amrex::FArrayBox& tube_ghost_data,
      const ::amrex::FArrayBox& plenum_mirror_data,
      const Eigen::Matrix<double, AMREX_SPACEDIM, 1>& normal) override final {
    impl_.FillTubeGhostLayer(tube_ghost_data, plenum_mirror_data, normal);
  }
  void FillPlenumGhostLayer(
      ::amrex::FArrayBox& plenum_ghost_data,
      const ::amrex::FArrayBox& tube_mirror_data) override final {
    impl_.FillPlenumGhostLayer(plenum_ghost_data, tube_mirror_data);
  }

  Boundary impl_;
};

/// \ingroup BoundaryCondition
///
class AnyMultiBlockBoundary {
public:
  /// Constructs coupled boundary states by pre computing mirror and ghost
  /// states for each of the specified domains.
  ///
  /// This function might grow the specified mirror boxes to an extent which is
  /// required to fulfill the specified ghost cell width requirements.
  template <typename Boundary>
  AnyMultiBlockBoundary(Boundary boundary,
                        const MultiBlockGriddingAlgorithm2& gridding,
                        const BlockConnection& connection, int gcw, int level)
      : impl_(std::make_unique<MultiBlockBoundaryWrapper<Boundary>>(
            std::move(boundary))),
        dir_{connection.direction}, side_{connection.side}, level_{level},
        gcw_{gcw}, connection_{connection} {
    Initialize(gridding, connection, gcw, level);
  }

  AnyMultiBlockBoundary(const AnyMultiBlockBoundary& other);

  void PreAdvanceHierarchy(const MultiBlockGriddingAlgorithm2& grid);

  /// Precompute Boundary states for each domain.
  ///
  /// Subsequent calls to FillBoundary will use these computed boundary states.
  ///
  /// \param[in] plenum  The higher dimensional patch hierarchy with geometry
  /// information. States here will be conservatively averaged and projected
  /// onto a one-dimensional space.
  ///
  /// \param[in] tube The low dimensional tube data.
  void ComputeBoundaryData(const cutcell::PatchHierarchy& plenum,
                           const PatchHierarchy& tube);

  /// Precompute Boundary states for each domain using the scratch.
  ///
  /// Subsequent calls to FillBoundary will use these computed boundary states.
  ///
  /// \param[in] plenum  The higher dimensional patch hierarchy with geometry
  /// information. States here will be conservatively averaged and projected
  /// onto a one-dimesnional space.
  ///
  /// \param[in] tube The low dimensional tube data.
  void ComputeBoundaryData(const cutcell::IntegratorContext& plenum,
                           const IntegratorContext& tube);

  void ComputeBoundaryDataForTube(const cutcell::IntegratorContext& plenum,
                                  const IntegratorContext& tube);

  void ComputeBoundaryDataForPlenum(const cutcell::IntegratorContext& plenum,
                                    const IntegratorContext& tube);

  /// Assuming that mf represents a MultiFab living in the higher dimensional
  /// plenum simulation its ghost layer will be filled with data from the tube
  /// simulation.
  void FillBoundary(::amrex::MultiFab& mf,
                    const cutcell::GriddingAlgorithm& gridding, int level);

  void FillBoundary(::amrex::MultiFab& mf,
                    const cutcell::GriddingAlgorithm& gridding, int level,
                    Direction dir) {
    if (dir == dir_) {
      FillBoundary(mf, gridding, level);
    }
  }

  /// Assuming that mf represents a MultiFab living in the one dimensional tube
  /// simulation its ghost layer will be filled with data from the plenum
  /// simulation.
  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level);

  void FillBoundary(::amrex::MultiFab& mf, const GriddingAlgorithm& gridding,
                    int level, Direction dir) {
    if (dir == dir_) {
      FillBoundary(mf, gridding, level);
    }
  }

  const ::amrex::FArrayBox* GetTubeMirrorData() const noexcept;
  const ::amrex::FArrayBox* GetTubeGhostData() const noexcept;

private:
  void Initialize(const MultiBlockGriddingAlgorithm2& gridding,
                  const BlockConnection& connection, int gcw, int level);

  std::unique_ptr<MultiBlockBoundaryBase> impl_;

  ::amrex::Box plenum_mirror_box_{};
  ::amrex::Box tube_mirror_box_{};

  std::unique_ptr<::amrex::FArrayBox> plenum_mirror_data_{};
  std::unique_ptr<::amrex::FArrayBox> tube_ghost_data_{};

  std::unique_ptr<::amrex::FArrayBox> tube_mirror_data_{};
  std::unique_ptr<::amrex::FArrayBox> plenum_ghost_data_{};

  Direction dir_{};
  int side_{};
  int level_{};
  int gcw_{};
  BlockConnection connection_{};
};

template <typename TubeEquation, typename PlenumEquation>
void ReduceStateDimension(TubeEquation& tube_equation,
                          Complete<TubeEquation>& dest,
                          PlenumEquation& /* plenum_equation */,
                          const Conservative<PlenumEquation>& src) {
  dest.density = src.density;
  for (int i = 0; i < dest.momentum.size(); ++i) {
    dest.momentum[i] = src.momentum[i];
  }
  if constexpr (euler::state_with_species<Complete<TubeEquation>>() &&
                euler::state_with_species<Complete<PlenumEquation>>()) {
    dest.species = src.species;
  }
  if constexpr (euler::state_with_passive_scalars<Complete<TubeEquation>>() &&
                euler::state_with_passive_scalars<Complete<PlenumEquation>>()) {
    dest.passive_scalars = src.passive_scalars;
  }
  dest.energy = src.energy;
  CompleteFromCons(tube_equation, dest, AsCons(dest));
}

template <typename PlenumEquation, typename TubeEquation>
void EmbedState(PlenumEquation& plenum_equation, nodeduce_t<Complete<PlenumEquation>&> dest,
                TubeEquation& /* tube_equation */,
                nodeduce_t<const Conservative<TubeEquation>&> src) {
  dest.density = src.density;
  dest.momentum.setZero();
  for (int i = 0; i < src.momentum.size(); ++i) {
    dest.momentum[i] = src.momentum[i];
  }
  if constexpr (euler::state_with_species<Complete<TubeEquation>>() &&
                euler::state_with_species<Complete<PlenumEquation>>()) {
    dest.species = src.species;
  }
  if constexpr (euler::state_with_passive_scalars<Complete<TubeEquation>>() &&
                euler::state_with_passive_scalars<Complete<PlenumEquation>>()) {
    dest.passive_scalars = src.passive_scalars;
  }
  dest.energy = src.energy;
  CompleteFromCons(plenum_equation, dest, AsCons(dest));
}

template <typename TubeEquation, typename PlenumEquation>
struct MultiBlockBoundary2 {
  static_assert(TubeEquation::Rank() <= PlenumEquation::Rank());

  MultiBlockBoundary2(const TubeEquation& tube_equation,
                      const PlenumEquation& plenum_equation)
      : tube_equation_(tube_equation), plenum_equation_(plenum_equation) {}

  MultiBlockBoundary2(const MultiBlockBoundary2& other)
      : tube_equation_(other.tube_equation_),
        plenum_equation_(other.plenum_equation_) {}

  TubeEquation tube_equation_;
  PlenumEquation plenum_equation_;

  void
  FillTubeGhostLayer(::amrex::FArrayBox& tube_ghost_data,
                     const ::amrex::FArrayBox& plenum_mirror_data,
                     const Eigen::Matrix<double, AMREX_SPACEDIM, 1>& normal) {
    BasicView cons_states = MakeView<const Conservative<PlenumEquation>>(
        plenum_mirror_data, plenum_equation_);
    BasicView complete_states =
        MakeView<Complete<TubeEquation>>(tube_ghost_data, tube_equation_);
    Conservative<PlenumEquation> cons(plenum_equation_);
    Conservative<PlenumEquation> rotated(plenum_equation_);
    Complete<TubeEquation> complete(tube_equation_);
    const std::ptrdiff_t i0 = plenum_mirror_data.box().smallEnd(0);
    const std::ptrdiff_t j0 = tube_ghost_data.box().smallEnd(0);
    ForEachIndex(Box<0>(complete_states), [&](std::ptrdiff_t j) {
      const std::ptrdiff_t k = j - j0;
      const std::ptrdiff_t i = i0 + k;
      if (normal.squaredNorm() == 0.0) {
        Load(cons, cons_states, {i});
      } else {
        Load(cons, cons_states, {0});
      }
      if (cons.density > 0.0) {
        if (normal.squaredNorm() > 0.0) {
          const Eigen::Matrix<double, AMREX_SPACEDIM, 1> unit =
              UnitVector<AMREX_SPACEDIM>(Direction::X);
          Rotate(rotated, cons, MakeRotation(normal, unit), plenum_equation_);
          ReduceStateDimension(tube_equation_, complete, plenum_equation_, rotated);
        } else {
          ReduceStateDimension(tube_equation_, complete, plenum_equation_, cons);
        }
        Store(complete_states, complete, {j});
      }
    });
  }

  void FillPlenumGhostLayer(::amrex::FArrayBox& plenum_ghost_data,
                            const ::amrex::FArrayBox& tube_mirror_data) {
    BasicView cons_states = MakeView<const Conservative<TubeEquation>>(
        tube_mirror_data, tube_equation_);
    BasicView complete_states =
        MakeView<Complete<PlenumEquation>>(plenum_ghost_data, plenum_equation_);
    Conservative<TubeEquation> cons(tube_equation_);
    Complete<PlenumEquation> complete(plenum_equation_);
    const std::ptrdiff_t i0 = Box<0>(cons_states).lower[0];
    const std::ptrdiff_t j0 = Box<0>(complete_states).lower[0];
    ForEachIndex(Box<0>(cons_states), [&](std::ptrdiff_t i) {
      const std::ptrdiff_t k = i - i0;
      const std::ptrdiff_t j = j0 + k;
      Load(cons, cons_states, {i});
      if (cons.density > 0.0) {
        EmbedState(plenum_equation_, complete, tube_equation_, cons);
        Store(complete_states, complete, {j});
      }
    });
  }
};

} // namespace fub::amrex

#endif
