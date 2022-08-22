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

#ifndef FUB_AMREX_MULTI_BLOCK_SOURCE_TERM_HPP
#define FUB_AMREX_MULTI_BLOCK_SOURCE_TERM_HPP

#include "fub/ext/outcome.hpp"

#include <range/v3/view/zip.hpp>
#include <vector>

namespace fub::amrex {

/// This class manages multiple kinetic source terms which are associated to
/// independend one-dimensional domains
template <typename SourceTerm> class MultiBlockSourceTerm {
public:
  static constexpr int Rank = SourceTerm::Rank;

  explicit MultiBlockSourceTerm(const std::vector<SourceTerm>& src_terms)
      : source_terms_(src_terms) {}

  explicit MultiBlockSourceTerm(std::vector<SourceTerm>&& src_terms)
      : source_terms_(std::move(src_terms)) {}

  template <typename MultiBlockGrid, typename = std::enable_if_t<is_detected<meta::ResetHierarchyConfiguration, SourceTerm&, const std::shared_ptr<MultiBlockGrid>&>::value>>
  void ResetHierarchyConfiguration(
      const std::shared_ptr<MultiBlockGrid>& grid) {
    auto&& tubes = grid->GetTubes();
    for (const auto& [tube, source] : ranges::view::zip(tubes, source_terms_)) {
      source.ResetHierarchyConfiguration(*tube);
    }
  }

  template <typename IntegratorContext>
  [[nodiscard]] Duration ComputeStableDt(const IntegratorContext& ctx, int level) noexcept {
    Duration stable_dt(std::numeric_limits<double>::max());
    auto tubes = ctx.Tubes();
    for (const auto& [tube, source] : ranges::view::zip(tubes, source_terms_)) {
      stable_dt = std::min(stable_dt, source.ComputeStableDt(tube, level));
    }
    return stable_dt;
  }

  /// \brief Integrates the source term for each tube in the specified context
  template <typename IntegratorContext>
  [[nodiscard]] Result<void, TimeStepTooLarge>
  AdvanceLevel(IntegratorContext& context, int level, Duration dt,
               [[maybe_unused]] const ::amrex::IntVect& ngrow = ::amrex::IntVect(0)) {
    auto&& tubes = context.Tubes();
    Result<void, TimeStepTooLarge> result = boost::outcome_v2::success();
    int i = 0;
    FUB_ASSERT(static_cast<std::size_t>(tubes.size()) == source_terms_.size());
    for (auto&& tube : tubes) {
      result = source_terms_[i].AdvanceLevel(tube, level, dt, ngrow);
      if (!result) {
        return result;
      }
      i = i + 1;
    }
    return result;
  }

private:
  std::vector<SourceTerm> source_terms_{};
};

} // namespace fub::amrex

#endif
