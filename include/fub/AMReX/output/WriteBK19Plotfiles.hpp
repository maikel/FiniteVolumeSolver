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

#ifndef FUB_AMREX_BK19_PLOT_FILES_HPP
#define FUB_AMREX_BK19_PLOT_FILES_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

#include <boost/log/common.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/trivial.hpp>

namespace fub::amrex {

template <int Rank, int VelocityRank>
struct WriteBK19Plotfile : public OutputAtFrequencyOrInterval<GriddingAlgorithm> {
  WriteBK19Plotfile(const ProgramOptions& options, const CompressibleAdvection<Rank, VelocityRank>& equation)
      : OutputAtFrequencyOrInterval(options), equation_(equation) {
    parent_path_ = GetOptionOr(options, "directory", parent_path_);
  }

  // void operator()(const CompressibleAdvectionIntegratorContext& context) const;
  void operator()(const GriddingAlgorithm& grid) override;

  private:
  CompressibleAdvection<Rank, VelocityRank> equation_;
  std::string parent_path_;
};

template <int Rank, int VelocityRank>
DebugSnapshot::ComponentNames GetCompleteVariableNames__(const CompressibleAdvection<Rank, VelocityRank>& equation) {
  using Equation = CompressibleAdvection<Rank, VelocityRank>;
  using Traits = StateTraits<Complete<Equation>>;
  constexpr auto names = Traits::names;
  const auto depths = Depths(equation, Type<Complete<Equation>>{});
  const std::size_t n_names =
      std::tuple_size<remove_cvref_t<decltype(names)>>::value;
  DebugSnapshot::ComponentNames varnames;
  varnames.reserve(n_names);
  boost::mp11::tuple_for_each(Zip(names, StateToTuple(depths)), [&](auto xs) {
    const int ncomp = std::get<1>(xs);
    if (ncomp == 1) {
      varnames.push_back(std::get<0>(xs));
    } else {
      for (int i = 0; i < ncomp; ++i) {
        varnames.push_back(fmt::format("{}_{}", std::get<0>(xs), i));
      }
    }
  });
  return varnames;
}

inline void WriteRawField(const std::string& path, const std::string& name,
                   const ::amrex::MultiFab& data, int level) {
  ::amrex::VisMF::SetHeaderVersion(::amrex::VisMF::Header::Version_v1);
  const std::string raw_pltname = fmt::format("{}/raw_fields", path);
  const std::string level_prefix = "Level_";
  const std::string full_prefix =
      ::amrex::MultiFabFileFullPrefix(level, raw_pltname, level_prefix, name);
  ::amrex::VisMF::Write(data, full_prefix);
}

template <int Rank, int VelocityRank>
void WriteBK19Plotfile<Rank, VelocityRank>::operator()(const GriddingAlgorithm& grid) {
  const fub::amrex::PatchHierarchy& hier = grid.GetPatchHierarchy();
  std::string name = fmt::format("{}/plt{:09}", parent_path_, grid.GetCycles());
  const int nlevels = hier.GetNumberOfLevels();
  const double time_point = hier.GetTimePoint().count();
  FUB_ASSERT(nlevels >= 0);
  std::size_t size = static_cast<std::size_t>(nlevels);
  ::amrex::Vector<const ::amrex::MultiFab*> mf(size);
  ::amrex::Vector<::amrex::Geometry> geoms(size);
  ::amrex::Vector<int> level_steps(size);
  ::amrex::Vector<::amrex::IntVect> ref_ratio(size);
  for (std::size_t i = 0; i < size; ++i) {
    mf[i] = &hier.GetPatchLevel(static_cast<int>(i)).data;
    geoms[i] = hier.GetGeometry(static_cast<int>(i));
    level_steps[i] = static_cast<int>(hier.GetCycles(static_cast<int>(i)));
    ref_ratio[i] = hier.GetRatioToCoarserLevel(static_cast<int>(i));
  }

  std::vector<std::string> varnames = GetCompleteVariableNames__(equation_);
  ::amrex::Vector<std::string> vnames(varnames.begin(), varnames.end());

  ::amrex::Vector<std::string> rfs{"raw_fields"};
  ::amrex::WriteMultiLevelPlotfile(name, nlevels, mf, vnames, geoms, time_point,
                                   level_steps, ref_ratio, "HyperCLaw-V1.1",
                                   "Level_", "Cell", rfs);

  // write nodal raw fields
  for (int level = 0; level < int(size); ++level) {
    if (hier.GetPatchLevel(level).nodes) {
      const ::amrex::MultiFab& nodes = *hier.GetPatchLevel(level).nodes;
      WriteRawField(name, "pi", nodes, level);
    }
  }
}

} // namespace fub::amrex

#endif
