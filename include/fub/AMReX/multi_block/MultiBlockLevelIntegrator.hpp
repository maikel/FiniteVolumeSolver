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

#ifndef FUB_AMREX_MULTI_BLOCK_LEVEL_INTEGRATOR_HPP
#define FUB_AMREX_MULTI_BLOCK_LEVEL_INTEGRATOR_HPP

#include "fub/AMReX/multi_block/MultiBlockIntegratorContext.hpp"

#include <functional>
#include <optional>

namespace fub::amrex {

template<typename T>
class assignable
{
public:
    assignable& operator=(assignable const& other)
    {
        value_.emplace(*other.value_);
        return *this;
    }

    assignable& operator=(assignable&& other) = default;
    assignable(assignable&& other) = default;
    assignable(assignable const& other) = default;
    
    assignable(T const& value) : value_(value) {}
    assignable(T&& value) : value_(std::move(value)) {}
    
    T const& get() const { return value_; }
    T& get() { return value_; }
    
    template<typename... Args>
    decltype(auto) operator()(Args&&... args)
    {
        return (*value_)(std::forward<Args>(args)...);
    }
private:
    std::optional<T> value_;
};

template<typename T>
class assignable<T&>
{
public:
    explicit assignable(T& value) : value_(value) {}
    
    T const& get() const { return value_; }
    T& get() { return value_; }
    
    template<typename... Args>
    decltype(auto) operator()(Args&&... args)
    {
        return value_(std::forward<Args>(args)...);
    }
private:
    std::reference_wrapper<T> value_;
};

template <typename AdapterAdvanceLevel, typename AdapterComputeStableDt,
          typename LevelIntegrator>
class MultiBlockLevelIntegrator {
private:
  assignable<AdapterAdvanceLevel> advance_level_map_;
  assignable<AdapterComputeStableDt> compute_stable_dt_map_;
  LevelIntegrator level_integrator_;

public:
  static constexpr int Rank = LevelIntegrator::Rank;

  MultiBlockLevelIntegrator(AdapterAdvanceLevel adapter1,
                            AdapterComputeStableDt adapter2,
                            LevelIntegrator level_integrator)
      : advance_level_map_(std::move(adapter1)),
        compute_stable_dt_map_(std::move(adapter2)),
        level_integrator_(std::move(level_integrator)) {}

  template <typename IntegratorContext>
  [[nodiscard]] Duration ComputeStableDt(const IntegratorContext& context,
                                         int level) {
    return std::invoke(compute_stable_dt_map_, level_integrator_, context,
                       level);
  }

  /// \brief Integrates the source term for each tube in the specified context
  template <typename IntegratorContext>
  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceLevel(IntegratorContext& context, int level, Duration dt,
               const ::amrex::IntVect& ngrow = ::amrex::IntVect(0)) {
    return std::invoke(advance_level_map_, level_integrator_, context, level,
                       dt, ngrow);
  }
};

} // namespace fub::amrex

#endif
