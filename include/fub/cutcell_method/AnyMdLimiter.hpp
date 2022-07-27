#ifndef FUB_CUTCELL_ANY_MD_LIMITER_HPP
#define FUB_CUTCELL_ANY_MD_LIMITER_HPP

#include "fub/ext/Eigen.hpp"
#include "fub/CutCellData.hpp"
#include "fub/PatchDataView.hpp"

#include <array>
#include <type_traits>
#include <memory>

namespace fub
{

template <int Rank> struct AnyLimiterBase;

template <int Rank> class AnyLimiter {
public:
  /// @{
  /// \name Constructors

  /// \brief This constructs a method that does nothing on invocation.
  AnyLimiter() = default;

  /// \brief Stores any object which satisfies the tagging method concept.
  template <typename T,
            typename = std::enable_if_t<!decays_to<T, AnyLimiter>()>>
  AnyLimiter(T&& tag);

  /// \brief Copies the implementation.
  AnyLimiter(const AnyLimiter& other);

  /// \brief Copies the implementation.
  AnyLimiter& operator=(const AnyLimiter& other);

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyLimiter(AnyLimiter&& other) noexcept = default;

  /// \brief Moves the `other` object without allocating and leaves an empty
  /// method.
  AnyLimiter& operator=(AnyLimiter&& other) noexcept = default;
  /// @}

  /// @{
  /// \name Member functions

  /// \brief Mask cells that need further refinement in a regridding procedure.
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
  /// @}

private:
  std::unique_ptr<AnyLimiterBase<Rank>> limiter_;
};

template <int Rank> struct LinearOptimizationLimiter {
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};
extern template struct LinearOptimizationLimiter<2>;

template <int Rank> struct NoMdLimiter {
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};
extern template struct NoMdLimiter<2>;

template <int Rank> struct UpwindMdLimiter {
  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const;
};
extern template struct UpwindMdLimiter<2>;

///////////////////////////////////////////////////////////////////////////////
//                                                              IMPLEMENTATION

template <int Rank> struct AnyLimiterBase {
  virtual ~AnyLimiterBase() = default;

  virtual std::unique_ptr<AnyLimiterBase<Rank>> Clone() const = 0;

  virtual void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const = 0;
};


template <int Rank, typename T>
struct LimiterWrapper_ : public AnyLimiterBase<Rank> {

  LimiterWrapper_(const T& limiter) : limiter_{limiter} {}
  LimiterWrapper_(T&& limiter) : limiter_{std::move(limiter)} {}

  void LimitGradientsAtIndex(
      const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
      StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
      const Index<Rank>& index, const Coordinates<Rank>& dx) const override {
    limiter_.LimitGradientsAtIndex(grad_u, u, geom, index, dx);
  }

  std::unique_ptr<AnyLimiterBase<Rank>> Clone() const override {
    return std::make_unique<LimiterWrapper_<Rank, T>>(limiter_);
  };

  T limiter_;
};

template <int Rank>
AnyLimiter<Rank>::AnyLimiter(const AnyLimiter& other)
    : limiter_(other.limiter_ ? other.limiter_->Clone() : nullptr) {}

template <int Rank>
AnyLimiter<Rank>& AnyLimiter<Rank>::operator=(const AnyLimiter& other) {
  AnyLimiter tmp(other);
  return *this = std::move(tmp);
}

template <int Rank>
template <typename T, typename>
AnyLimiter<Rank>::AnyLimiter(T&& tag)
    : limiter_{std::make_unique<LimiterWrapper_<Rank, remove_cvref_t<T>>>(
          std::forward<T>(tag))} {}

template <int Rank>
void AnyLimiter<Rank>::LimitGradientsAtIndex(
    const std::array<StridedDataView<double, Rank>, Rank>& grad_u,
    StridedDataView<const double, Rank> u, const CutCellData<Rank>& geom,
    const Index<Rank>& index, const Coordinates<Rank>& dx) const {
  if (limiter_) {
    return limiter_->LimitGradientsAtIndex(grad_u, u, geom, index, dx);
  }
}

}

#endif