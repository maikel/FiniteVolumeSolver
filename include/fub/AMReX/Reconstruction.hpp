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

#ifndef FUB_AMREX_RECONSTRUCTION_HPP
#define FUB_AMREX_RECONSTRUCTION_HPP

#include "fub/CompleteFromCons.hpp"
#include "fub/ForEach.hpp"
#include "fub/State.hpp"

#include <AMReX_MultiFab.H>

namespace fub::amrex {

struct ReconstructionBase {
  virtual ~ReconstructionBase() = default;
  virtual void CompleteFromCons(::amrex::MultiFab& dest,
                                const ::amrex::MultiFab& src) = 0;
  virtual std::unique_ptr<ReconstructionBase> Clone() const = 0;
};

template <typename T> struct ReconstructionWrapper : ReconstructionBase {
  ReconstructionWrapper(const T& x) : reconstruction_{x} {}
  ReconstructionWrapper(T&& x) : reconstruction_{x} {}

  void CompleteFromCons(::amrex::MultiFab& dest,
                        const ::amrex::MultiFab& src) override {
    reconstruction_.CompleteFromCons(dest, src);
  }

  std::unique_ptr<ReconstructionBase> Clone() const override {
    return std::make_unique<ReconstructionWrapper<T>>(reconstruction_);
  }

  T reconstruction_;
};

class Reconstruction {
public:
  Reconstruction() = delete;

  Reconstruction(const Reconstruction& other)
      : reconstruction_{other.reconstruction_->Clone()} {}
  Reconstruction& operator=(const Reconstruction& other) {
    reconstruction_ = other.reconstruction_->Clone();
    return *this;
  }

  Reconstruction(Reconstruction&& other) = default;
  Reconstruction& operator=(Reconstruction&& other) = default;

  template <typename R, typename = std::enable_if_t<
                            !std::is_same_v<std::decay_t<R>, Reconstruction>>>
  Reconstruction(R&& rec)
      : reconstruction_{std::make_unique<ReconstructionWrapper<R>>(rec)} {}

  void CompleteFromCons(::amrex::MultiFab& dest, const ::amrex::MultiFab& src) {
    reconstruction_->CompleteFromCons(dest, src);
  }

private:
  std::unique_ptr<ReconstructionBase> reconstruction_;
};

template <typename Equation> class ReconstructEquationStates {
public:
  explicit ReconstructEquationStates(const Equation& eq) : equation_{eq} {}

  void CompleteFromCons(::amrex::MultiFab& dest, const ::amrex::MultiFab& src) {
#if defined(_OPENMP) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
    for (::amrex::MFIter mfi(src, true); mfi.isValid(); ++mfi) {
      View<Complete<Equation>> complete =
          MakeView<Complete<Equation>>(dest[mfi], equation_, mfi.tilebox());
      View<const Conservative<Equation>> conservative =
          MakeView<const Conservative<Equation>>(src[mfi], equation_,
                                                 mfi.tilebox());
      ::fub::CompleteFromCons(equation_, complete, conservative);
    }
  }

private:
  Equation equation_;
};

} // namespace fub::amrex

#endif
