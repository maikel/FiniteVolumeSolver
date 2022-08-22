
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

#ifndef FUB_AMREX_OUTPUT_MULTI_WRITE_HDF5_HPP
#define FUB_AMREX_OUTPUT_MULTI_WRITE_HDF5_HPP

#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm2.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"

#include <optional>
#include <string>

namespace fub::amrex {

class MultiWriteHdf5
    : public OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm> {
public:
  MultiWriteHdf5(const ProgramOptions& options);
  virtual ~MultiWriteHdf5() = default;

  void operator()(const MultiBlockGriddingAlgorithm& grid) override;

private:
  enum class Type { plenum, tube };
  Type type_{Type::tube};
  int grid_id_{0};
  std::string path_to_file_{};
  std::optional<::amrex::Box> output_box_{};

  virtual span<const std::string> PlenumFieldNames() const noexcept {
    return {};
  }
  virtual span<const std::string> TubeFieldNames() const noexcept { return {}; }
};

class MultiWriteHdf52
    : public OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2> {
public:
  MultiWriteHdf52(const ProgramOptions& vm);
  virtual ~MultiWriteHdf52() = default;

  void operator()(const MultiBlockGriddingAlgorithm2& grid) override;

private:
  enum class Type { plenum, tube };
  Type type_{Type::tube};
  int grid_id_{0};
  std::string path_to_file_{};
  std::optional<::amrex::Box> output_box_{};

  virtual span<const std::string> PlenumFieldNames() const noexcept {
    return {};
  }
  virtual span<const std::string> TubeFieldNames() const noexcept { return {}; }
};

class MultiWriteHdf5WithNames : public MultiWriteHdf52, public MultiWriteHdf5 {
public:
  template <typename PEquation, typename TEquation>
  MultiWriteHdf5WithNames(const ProgramOptions& options, const PEquation& peq,
                          const TEquation& teq)
      : MultiWriteHdf52(options), MultiWriteHdf5(options),
        plenum_field_names_{
            VarNames<Complete<PEquation>, std::vector<std::string>>(peq)},
        tube_field_names_{
            VarNames<Complete<TEquation>, std::vector<std::string>>(teq)} {
    plenum_field_names_.emplace_back("vfrac");
  }

private:
  std::vector<std::string> plenum_field_names_;
  std::vector<std::string> tube_field_names_;

  span<const std::string> PlenumFieldNames() const noexcept override final {
    return plenum_field_names_;
  }

  span<const std::string> TubeFieldNames() const noexcept override final {
    return tube_field_names_;
  }
};

} // namespace fub::amrex

#endif