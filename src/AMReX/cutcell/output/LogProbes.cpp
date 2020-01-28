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

#include "fub/AMReX/cutcell/output/LogProbes.hpp"
#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include <boost/filesystem.hpp>
#include <fmt/format.h>
#include <fstream>

namespace fub::amrex {
namespace {
template <typename T>
using ProbesView = basic_mdspan<T, extents<AMREX_SPACEDIM, dynamic_extent>>;

void LogTubeProbes(const std::string& p, ProbesView<const double> probes,
                   const fub::amrex::PatchHierarchy& hierarchy, MPI_Comm comm) {
  int rank = -1;
  MPI_Comm_rank(comm, &rank);

  std::vector<double> buffer = GatherStates(hierarchy, probes, comm);
  if (rank == 0) {
    std::ofstream stream = [&p] {
      boost::filesystem::path path(p);
      if (boost::filesystem::is_regular_file(p)) {
        return std::ofstream(p, std::ofstream::app);
      } else if (boost::filesystem::exists(p)) {
        throw std::runtime_error(fmt::format(
            "The specified path '{}' exists but is not regular file.", p));
      }
      boost::filesystem::path dir = path.parent_path();
      if (!boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
      }
      return std::ofstream(p);
    }();
    fub::mdspan<const double, 2> states(buffer.data(), probes.extent(1),
                                        buffer.size() / probes.extent(1));
    fub::Burke2012 burke2012{};
    std::array<double, 11> molar_mass{};
    burke2012.getMolarMass(molar_mass);
    auto sH2 = fub::Burke2012::sH2;
    auto sO2 = fub::Burke2012::sO2;
    auto sH2O = fub::Burke2012::sH2O;
    for (int i = 0; i < probes.extent(1); ++i) {
      const double rho = states(i, 0);
      const double u = states(i, 1) / rho;
      const double h2 = states(i, 3 + sH2) / molar_mass[size_t(sH2)];
      const double o2 = states(i, 3 + sO2) / molar_mass[size_t(sO2)];
      const double h2o = states(i, 3 + sH2O) / molar_mass[size_t(sH2O)];
      const double T = states(i, 16);
      const double p = states(i, 14);
      const double a = states(i, 15);
      const double x = probes(0, i);
      const double t = hierarchy.GetTimePoint().count();
      stream << fmt::format(
          "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}"
          "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}\n",
          t, x, rho, u, T, p, a, h2, o2, h2o);
    }
  }
}

void LogPlenumProbes(const std::string& p, ProbesView<const double> probes,
                     const fub::amrex::cutcell::PatchHierarchy& hierarchy,
                     MPI_Comm comm) {
  int rank = -1;
  MPI_Comm_rank(comm, &rank);

  std::vector<double> buffer = GatherStates(hierarchy, probes, comm);
  if (rank == 0) {
    std::ofstream stream = [&p] {
      boost::filesystem::path path(p);
      if (boost::filesystem::is_regular_file(p)) {
        return std::ofstream(p, std::ofstream::app);
      } else if (boost::filesystem::exists(p)) {
        throw std::runtime_error(fmt::format(
            "The specified path '{}' exists but is not regular file.", p));
      }
      boost::filesystem::path dir = path.parent_path();
      if (!boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
      }
      return std::ofstream(p);
    }();
    fub::mdspan<const double, 2> states(buffer.data(), probes.extent(1),
                                        buffer.size() / probes.extent(1));
    fub::Burke2012 burke2012{};
    std::array<double, 11> molar_mass{};
    burke2012.getMolarMass(molar_mass);
    auto sH2 = fub::Burke2012::sH2;
    auto sO2 = fub::Burke2012::sO2;
    auto sH2O = fub::Burke2012::sH2O;
    for (int i = 0; i < probes.extent(1); ++i) {
      const double rho = states(i, 0);
      const double u = states(i, 1) / rho;
      const double v = states(i, 2) / rho;
      const double w = states(i, 3) / rho;
      const double h2 = states(i, 5 + sH2) / molar_mass[size_t(sH2)];
      const double o2 = states(i, 5 + sO2) / molar_mass[size_t(sO2)];
      const double h2o = states(i, 5 + sH2O) / molar_mass[size_t(sH2O)];
      const double a = states(i, 17);
      const double T = states(i, 18);
      const double p = states(i, 16);
      const double x = probes(0, i);
      const double y = probes(1, i);
      const double z = probes(2, i);
      const double t = hierarchy.GetTimePoint().count();
      stream << fmt::format(
          "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}"
          "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}"
          "{:< 24.15g}{:< 24.15g}{:< 24.15g}{:< 24.15g}\n",
          t, x, y, z, rho, u, v, w, T, p, a, h2, o2, h2o);
    }
  }
}

} // namespace

LogProbesOutput::LogProbesOutput(
    const std::map<std::string, pybind11::object>& vm)
    : OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm>(vm) {
  plenum_output_path_ =
      GetOptionOr(vm, "plenum_filename", std::string("plenum.dat"));
  tube_output_path_ = GetOptionOr(vm, "tube_filename", std::string("tube.dat"));
  {
    const std::map<std::string, pybind11::object> plenum_vm =
        ToMap(GetOptionOr(vm, "plenum", pybind11::dict()));
    std::vector<double> xs =
        GetOptionOr(plenum_vm, "x_coordinates", std::vector<double>());
    std::vector<double> ys =
        GetOptionOr(plenum_vm, "y_coordinates", std::vector<double>());
    std::vector<double> zs =
        GetOptionOr(plenum_vm, "z_coordinates", std::vector<double>());
    FUB_ASSERT(xs.size() == ys.size());
    FUB_ASSERT(zs.size() == ys.size());
    plenum_probes_.resize(AMREX_SPACEDIM * xs.size());
    ProbesView<double> probes(plenum_probes_.data(), xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
      probes(0, i) = xs[i];
      probes(1, i) = ys[i];
      probes(2, i) = zs[i];
    }
  }
  const std::map<std::string, pybind11::object> tube_vm =
      ToMap(GetOptionOr(vm, "tube", pybind11::dict()));
  std::vector<double> xs =
      GetOptionOr(tube_vm, "x_coordinates", std::vector<double>());
  std::vector<double> ys =
      GetOptionOr(tube_vm, "y_coordinates", std::vector<double>());
  std::vector<double> zs =
      GetOptionOr(tube_vm, "z_coordinates", std::vector<double>());
  FUB_ASSERT(xs.size() == ys.size());
  FUB_ASSERT(zs.size() == ys.size());
  tube_probes_.resize(AMREX_SPACEDIM * xs.size());
  ProbesView<double> probes(tube_probes_.data(), xs.size());
  for (std::size_t i = 0; i < xs.size(); ++i) {
    probes(0, i) = xs[i];
    probes(1, i) = ys[i];
    probes(2, i) = zs[i];
  }
}

void LogProbesOutput::operator()(const MultiBlockGriddingAlgorithm& grid) {
  const std::ptrdiff_t n_tubes = grid.GetTubes().size();
  const std::ptrdiff_t n_tube_probes = static_cast<std::ptrdiff_t>(
      tube_probes_.size() / n_tubes / AMREX_SPACEDIM);
  mdspan<const double, 3> tube_coords(tube_probes_.data(), AMREX_SPACEDIM,
                                      n_tube_probes, n_tubes);
  MPI_Comm comm = ::amrex::ParallelDescriptor::Communicator();
  std::ptrdiff_t i_block = 0;
  for (const auto& tube : grid.GetTubes()) {
    ProbesView<const double> probes(&tube_coords(0, 0, i_block), n_tube_probes);
    LogTubeProbes(tube_output_path_, probes, tube->GetPatchHierarchy(), comm);
    ++i_block;
  }
  i_block = 0;
  const std::ptrdiff_t n_plena = grid.GetPlena().size();
  const std::ptrdiff_t n_plenum_probes = static_cast<std::ptrdiff_t>(
      plenum_probes_.size() / n_plena / AMREX_SPACEDIM);
  mdspan<const double, 3> plenum_coords(plenum_probes_.data(), AMREX_SPACEDIM,
                                        n_plenum_probes, n_plena);
  for (const auto& plenum : grid.GetPlena()) {
    ProbesView<const double> probes(&plenum_coords(0, 0, i_block),
                                    n_plenum_probes);
    LogPlenumProbes(plenum_output_path_, probes, plenum->GetPatchHierarchy(),
                    comm);
    ++i_block;
  }
}

} // namespace fub::amrex
