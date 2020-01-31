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

#ifndef FUB_AMREX_DEBUG_OUTPUT_HPP
#define FUB_AMREX_DEBUG_OUTPUT_HPP

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

namespace fub::amrex {

class DebugOutput{
public:
  using Hierarchy = std::vector<::amrex::MultiFab>;
  using ComponentNames = std::vector<std::string>;

  DebugOutput();
  explicit DebugOutput(const ProgramOptions& options);

  void SaveData(const ::amrex::MultiFab& mf, const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const ::amrex::MultiFab& mf, const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void SaveData(const std::vector<const ::amrex::MultiFab*>& hierarchy,
                const std::string& name,
                ::amrex::SrcComp component = ::amrex::SrcComp(0));

  void SaveData(const std::vector<const ::amrex::MultiFab*>& hierarchy,
                const ComponentNames& names,
                ::amrex::SrcComp first_component = ::amrex::SrcComp(0));

  void ClearAll();

  std::vector<std::pair<Hierarchy, ComponentNames>> GatherCellCentered() const;

private:
  std::vector<ComponentNames> names_{};
  std::vector<Hierarchy> cells_{};
  std::vector<Hierarchy> faces_x_{};
  std::vector<Hierarchy> faces_y_{};
  std::vector<Hierarchy> faces_z_{};
  std::vector<Hierarchy> nodes_{};
};

} // namespace fub::amrex

#endif