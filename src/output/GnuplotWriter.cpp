// Copyright (c) 2018 Maikel Nadolski
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

#include "fub/output/GnuplotWriter.hpp"

#include "fub/SAMRAI/utility.hpp"
#include "fub/core/span.hpp"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"

#include <boost/filesystem.hpp>

#include <array>
#include <cstring>
#include <fstream>
#include <numeric>
#include <string>

namespace fub {
namespace {
std::size_t strlen_(const std::string& str) { return str.size(); }
template <std::size_t N> std::size_t strlen_(const char (&str)[N]) { return N; }
std::size_t strlen_(const char* str) { return std::strlen(str); }

#ifdef __cpp_fold_expressions
template <typename... Strings> std::string strcat_(Strings&&... strings) {
  std::string result;
  result.reserve((0 + ... + strlen_(strings)));
  (result.append(std::forward<Strings>(strings)), ...);
  return result;
}
#else
void strcat_helper_(std::string& result) {}
template <typename String, typename... Strings>
void strcat_helper_(std::string& result, String&& string,
                    Strings&&... strings) {
  result.append(std::forward<String>(string));
  strcat_helper_(result, std::forward<Strings>(strings)...);
}

template <typename... Strings> std::string strcat_(Strings&&... strings) {
  std::array<std::size_t, sizeof...(Strings)> sizes{strlen_(strings)...};
  std::string result;
  result.reserve(std::accumulate(sizes.begin(), sizes.end(), std::size_t{}));
  strcat_helper_(result, std::forward<Strings>(strings)...);
  return result;
}
#endif

} // namespace

/// Constructs the class with a directory which will store the output files.
GnuplotWriter::GnuplotWriter(std::string directory, span<const int> ids)
    : directory_{std::move(directory)} {
  namespace filesystem = boost::filesystem;
  filesystem::path path = directory_;
  if (filesystem::exists(path) && !filesystem::is_directory(path)) {
    throw std::runtime_error(
        strcat_("The specified path '", directory_,
                "' exists already and is not a directory."));
  }
  filesystem::create_directory(path);
  const int rank = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();
  filesystem::path file =
      path / strcat_("gnuplot_", std::to_string(rank), ".out");
  out_ = std::make_unique<std::ofstream>(file.string());

  for (int id : ids) {
    const SAMRAI::hier::VariableDatabase& db =
        *SAMRAI::hier::VariableDatabase::getDatabase();
    std::shared_ptr<SAMRAI::hier::Variable> variable{};
    db.mapIndexToVariable(id, variable);
    if (variable) {
      const std::string& name = variable->getName();
      registerPlotQuantity(name, VariableType::scalar, id);
    }
  }
}

/// Adds a specified patch data to the list of output variables.
void GnuplotWriter::registerPlotQuantity(std::string name, VariableType type,
                                         int patch_data_id) {
  FUB_ASSERT(names_.size() == types_.size());
  FUB_ASSERT(names_.size() == patch_data_ids_.size());
  names_.emplace_back(std::move(name));
  try {
    types_.emplace_back(type);
  } catch (...) {
    names_.pop_back();

    FUB_ASSERT(names_.size() == types_.size());
    FUB_ASSERT(names_.size() == patch_data_ids_.size());
    throw;
  }
  try {
    patch_data_ids_.emplace_back(patch_data_id);
  } catch (...) {
    types_.pop_back();
    names_.pop_back();

    FUB_ASSERT(names_.size() == types_.size());
    FUB_ASSERT(names_.size() == patch_data_ids_.size());
    throw;
  }
  FUB_ASSERT(names_.size() == types_.size());
  FUB_ASSERT(names_.size() == patch_data_ids_.size());
}

namespace {
void writeCellData_(std::ostream& out, const SAMRAI::hier::Patch& patch,
                    span<const int> ids) {
  const SAMRAI::hier::Box& box = patch.getBox();
  const SAMRAI::geom::CartesianPatchGeometry& geom =
      *GetCartesianPatchGeometry(patch);
  for (const SAMRAI::hier::Index& index : box) {
    Coordinates x = ComputeCellCoordinates(geom, box, index);
    for (double xi : x) {
      out << std::setw(20) << xi;
    }
    for (int id : ids) {
      auto data = std::dynamic_pointer_cast<SAMRAI::pdat::CellData<double>>(
          patch.getPatchData(id));
      if (data) {
        SAMRAI::pdat::CellIndex cell(index);
        const int depth = data->getDepth();
        for (int d = 0; d < depth; ++d) {
          out << std::setw(20) << std::setprecision(12) << (*data)(cell, d);
        }
      }
    }
    out << '\n';
  }
}
} // namespace

/// Writes all quantitites which has been registered with this class to the
/// output directory.
void GnuplotWriter::writePlotData(
    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy, int cylce,
    double time_point) {
  forEachPatch(*hierarchy, [&](const SAMRAI::hier::Patch& patch) {
    auto id_iter = patch_data_ids_.begin();
    (*out_) << std::setw(20) << "# Coordinates ";
    for (const std::string& name : names_) {
      const int pdi = *id_iter++;
      auto data = std::dynamic_pointer_cast<SAMRAI::pdat::CellData<double>>(
          patch.getPatchData(pdi));
      if (data) {
        const int depth = data->getDepth();
        if (depth == 1) {
          (*out_) << std::setw(20) << strcat_(name, " ");
        } else {
          for (int d = 0; d < depth; ++d) {
            (*out_) << std::setw(20)
                    << strcat_(name, "_", std::to_string(d), " ");
          }
        }
      }
    }
    (*out_) << '\n';
    writeCellData_(*out_, patch, patch_data_ids_);
  });
  (*out_) << "\n\n";
}
} // namespace fub