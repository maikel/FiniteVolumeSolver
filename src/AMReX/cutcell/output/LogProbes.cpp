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

#include "fub/output/Hdf5Handle.hpp"

#include <boost/filesystem.hpp>
#include <fmt/format.h>
#include <fstream>

namespace fub::amrex {
namespace {
template <typename T>
using ProbesView = basic_mdspan<T, extents<AMREX_SPACEDIM, dynamic_extent>>;

template <typename GriddingAlgorithm>
void CreateHdf5Database(const std::string& name, const GriddingAlgorithm& grid,
                        ProbesView<const double> probes,
                        mdspan<const double, 2> states) {
  //SeverityLogger log = GetLogger(boost::log::trivial::info);
  // BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesOutput");
  // BOOST_LOG(log) << "Create HDF5 file  '" << name << "'.";
  boost::filesystem::path path(name);
  boost::filesystem::path dir = boost::filesystem::absolute(
      path.parent_path(), boost::filesystem::current_path());
  if (!boost::filesystem::exists(dir)) {
    boost::filesystem::create_directories(dir);
  }
  H5File file(H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  const hsize_t n_probes = static_cast<hsize_t>(probes.extent(1));
  const hsize_t n_fields = static_cast<hsize_t>(states.extent(1));
  // BOOST_LOG(log) << "Write initial probe data.";
  // create a main data set with probe data for each time step
  {
    constexpr std::size_t size = 3;
    hsize_t dims[size] = {n_fields, n_probes, 1};
    hsize_t maxdims[size] = {dims[0], dims[1], H5S_UNLIMITED};
    H5Properties properties(H5Pcreate(H5P_DATASET_CREATE));
    H5Pset_chunk(properties, size, dims);
    H5Pset_alloc_time(properties, H5D_ALLOC_TIME_EARLY);
    H5Space dataspace(H5Screate_simple(size, dims, maxdims));
    H5Dataset dataset(H5Dcreate(file, "/data", H5T_IEEE_F64LE, dataspace,
                                H5P_DEFAULT, properties, H5P_DEFAULT));
    H5Dwrite(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             states.data());

    // Write attribute data which describes how many probes are in it
    hsize_t scalar_value = 1;
    H5Space spacedim(H5Screate_simple(1, &scalar_value, nullptr));
    H5Attribute nprobes_attr(H5Acreate2(dataset, "n_probes", H5T_STD_I32LE,
                                        spacedim, H5P_DEFAULT, H5P_DEFAULT));
    int nprobes = static_cast<int>(n_probes);
    H5Awrite(nprobes_attr, H5T_STD_I32LE, &nprobes);
  }
  // create a dataset to store all coordinates of the probes
  {
    std::array<hsize_t, 2> dims{AMREX_SPACEDIM, n_probes};
    H5Space dataspace(H5Screate_simple(2, dims.data(), nullptr));
    H5Dataset dataset(H5Dcreate(file, "/coordinates", H5T_IEEE_F64LE, dataspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    H5Dwrite(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             probes.data());
  }
  hsize_t dims[1] = {1};
  hsize_t maxdims[1] = {H5S_UNLIMITED};
  H5Properties properties(H5Pcreate(H5P_DATASET_CREATE));
  H5Pset_chunk(properties, 1, dims);
  H5Space dataspace(H5Screate_simple(1, dims, maxdims));

  H5Dataset times(H5Dcreate(file, "/times", H5T_IEEE_F64LE, dataspace,
                            H5P_DEFAULT, properties, H5P_DEFAULT));
  const double count = grid.GetTimePoint().count();
  H5Dwrite(times, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &count);

  H5Dataset cycles(H5Dcreate(file, "/cycles", H5T_STD_I64LE_g, dataspace,
                             H5P_DEFAULT, properties, H5P_DEFAULT));
  const std::int64_t cycle_count = static_cast<std::int64_t>(grid.GetCycles());
  H5Dwrite(cycles, H5T_STD_I64LE_g, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &cycle_count);
}

template <typename GriddingAlgorithm>
void OpenHdf5Database(const std::string& name, const GriddingAlgorithm& grid,
                      mdspan<const double, 2> states) {
  //SeverityLogger log = GetLogger(boost::log::trivial::info);
  //BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesOutput");
  //BOOST_LOG(log) << "Append probes to HDF5 file  '" << name << "'.";
  H5File file(H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  if (file < 0) {
    return;
  }
  if (H5Dataset dataset(H5Dopen(file, "/data", H5P_DEFAULT)); dataset < 0) {
    return;
  } else {
    auto map_last = [](const std::array<hsize_t, 3>& x, auto f) {
      return std::array<hsize_t, 3>{x[0], x[1], hsize_t(f(x[2]))};
    };
    std::array<hsize_t, 3> dims = {};
    std::array<hsize_t, 3> maxdims = {};
    H5Space dataspace(H5Dget_space(dataset));
    H5Sget_simple_extent_dims(dataspace, dims.data(), maxdims.data());
    FUB_ASSERT(maxdims[2] == H5S_UNLIMITED);
    std::array<hsize_t, 3> new_dims =
        map_last(dims, [](hsize_t dim) { return dim + 1; });
    H5Dset_extent(dataset, new_dims.data());
    H5Space filespace(H5Dget_space(dataset));
    std::array<hsize_t, 3> count = map_last(dims, [](hsize_t) { return 1; });
    std::array<hsize_t, 3> offset{0, 0, dims[2]};
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), nullptr,
                        count.data(), nullptr);
    H5Space memspace(H5Screate_simple(3, count.data(), nullptr));
    H5Dwrite(dataset, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT,
             states.data());
  }
  if (H5Dataset times(H5Dopen(file, "/times", H5P_DEFAULT)); times < 0) {
    return;
  } else {
    hsize_t dims{};
    hsize_t maxdims{};
    H5Space dataspace(H5Dget_space(times));
    H5Sget_simple_extent_dims(dataspace, &dims, &maxdims);
    FUB_ASSERT(maxdims == H5S_UNLIMITED);
    const hsize_t new_dims = dims + 1;
    H5Dset_extent(times, &new_dims);
    H5Space filespace(H5Dget_space(times));
    const hsize_t count = 1;
    const hsize_t offset = dims;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, &count, &count,
                        &count);
    const double tp_count = grid.GetTimePoint().count();
    H5Space memspace(H5Screate_simple(1, &count, nullptr));
    H5Dwrite(times, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT,
             &tp_count);
  }
  if (H5Dataset cycles(H5Dopen(file, "/cycles", H5P_DEFAULT)); cycles < 0) {
    return;
  } else {
    hsize_t dims{};
    hsize_t maxdims{};
    H5Space dataspace(H5Dget_space(cycles));
    H5Sget_simple_extent_dims(dataspace, &dims, &maxdims);
    FUB_ASSERT(maxdims == H5S_UNLIMITED);
    const hsize_t new_dims = dims + 1;
    H5Dset_extent(cycles, &new_dims);
    H5Space filespace(H5Dget_space(cycles));
    const hsize_t count = 1;
    const hsize_t offset = dims;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr, &count,
                        nullptr);
    const int cycle_count = static_cast<int>(grid.GetCycles());
    H5Space memspace(H5Screate_simple(1, &count, nullptr));
    H5Dwrite(cycles, H5T_STD_I32LE, memspace, filespace, H5P_DEFAULT,
             &cycle_count);
  }
}

void RenameFileIfExists(const std::string& filename) {
  SeverityLogger log = GetLogger(boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesOutput");
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    if (boost::filesystem::exists(filename)) {
      int old_num = 1;
      std::string backup_name =
          fmt::format("{}.old.{}", filename, old_num);
      while (boost::filesystem::exists(backup_name)) {
        old_num += 1;
        backup_name = fmt::format("{}.old.{}", filename, old_num);
      }
      BOOST_LOG(log) << fmt::format(
          "File '{}' exists already. Rename to '{}'.", filename,
          backup_name);
      boost::filesystem::rename(filename, backup_name);
    }
  }
}

} // namespace

LogProbesOutput::LogProbesOutput(const ProgramOptions& vm)
    : OutputAtFrequencyOrInterval<MultiBlockGriddingAlgorithm2>(vm) {
  using namespace std::literals;
  {
    const ProgramOptions plenum_vm = GetOptions(vm, "Plenum");
    plenum_output_path_ = GetOptionOr(plenum_vm, "filename", "Plenum.h5"s);
    RenameFileIfExists(plenum_output_path_);
    std::vector<std::array<double, AMREX_SPACEDIM>> xs{};
    xs = GetOptionOr(plenum_vm, "coordinates", xs);
    plenum_probes_.resize(AMREX_SPACEDIM * xs.size());
    ProbesView<double> probes(plenum_probes_.data(), xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i) {
      probes(0, i) = xs[i][0];
      probes(1, i) = xs[i][1];
      probes(2, i) = xs[i][2];
    }
  }
  {
    const ProgramOptions tube_vm = GetOptions(vm, "Tube");
    tube_output_path_ = GetOptionOr(tube_vm, "filename", "Tube"s);
    n_tubes_ = GetOptionOr(tube_vm, "n_tubes", 1);
    RenameFileIfExists(tube_output_path_);
    std::vector<std::array<double, AMREX_SPACEDIM>> xs{};
    xs = GetOptionOr(tube_vm, "coordinates", xs);
    tube_probes_.resize(AMREX_SPACEDIM * xs.size());
    const int n_tube_probes = xs.size() / n_tubes_;
    mdspan<double, 3> tube_coords(tube_probes_.data(), AMREX_SPACEDIM,
                                      n_tube_probes, n_tubes_);
    for (int tube = 0; tube < n_tubes_; ++tube) {
      ProbesView<double> probes(&tube_coords(0, 0, tube), n_tube_probes);
      for (int i = 0; i < n_tube_probes; ++i) {
        probes(0, i) = xs[i][0];
        probes(1, i) = xs[i][1];
        probes(2, i) = xs[i][2];
      }
    }
  }
}

void LogProbesOutput::operator()(const MultiBlockGriddingAlgorithm2& grid) {
  const std::ptrdiff_t n_tube_probes =
      static_cast<std::ptrdiff_t>(tube_probes_.size()) / n_tubes_ /
      AMREX_SPACEDIM;
  mdspan<const double, 3> tube_coords(tube_probes_.data(), AMREX_SPACEDIM,
                                      n_tube_probes, n_tubes_);
  MPI_Comm comm = ::amrex::ParallelDescriptor::Communicator();
  int rank = -1;
  MPI_Comm_rank(comm, &rank);
  std::ptrdiff_t i_block = 0;
  for (const auto& tube : grid.GetTubes()) {
    ProbesView<const double> probes(&tube_coords(0, 0, i_block), n_tube_probes);
    std::vector<double> buffer = GatherStates(tube->GetPatchHierarchy(), probes, comm);
    if (rank == 0 && buffer.size() > 0) {
      fub::mdspan<const double, 2> states(
          buffer.data(), probes.extent(1),
          static_cast<std::ptrdiff_t>(buffer.size()) / probes.extent(1));
      std::string filename = fmt::format("{}_{}.h5", tube_output_path_, i_block);
      if (!boost::filesystem::exists(filename)) {
        CreateHdf5Database(filename, *tube, probes, states);
      } else if (boost::filesystem::is_regular_file(filename)) {
        OpenHdf5Database(filename, *tube, states);
      } else {
        SeverityLogger log = GetLogger(boost::log::trivial::warning);
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesOutput");
        BOOST_LOG(log) << fmt::format(
            "Path '{}' points to some non regular file. No probes output.",
            filename);
      }
    }
    ++i_block;
  }
  i_block = 0;
  const std::ptrdiff_t n_plena = grid.GetPlena().size();
  const std::ptrdiff_t n_plenum_probes =
      static_cast<std::ptrdiff_t>(plenum_probes_.size()) / n_plena /
      AMREX_SPACEDIM;
  mdspan<const double, 3> plenum_coords(plenum_probes_.data(), AMREX_SPACEDIM,
                                        n_plenum_probes, n_plena);
  for (const auto& plenum : grid.GetPlena()) {
    ProbesView<const double> probes(&plenum_coords(0, 0, i_block),
                                    n_plenum_probes);
    std::vector<double> buffer = GatherStates(plenum->GetPatchHierarchy(), probes, comm);
    if (rank == 0 && buffer.size() > 0) {
      fub::mdspan<const double, 2> states(
          buffer.data(), probes.extent(1),
          static_cast<std::ptrdiff_t>(buffer.size()) / probes.extent(1));
      if (!boost::filesystem::exists(plenum_output_path_)) {
        CreateHdf5Database(plenum_output_path_, *plenum, probes, states);
      } else if (boost::filesystem::is_regular_file(plenum_output_path_)) {
        OpenHdf5Database(plenum_output_path_, *plenum, states);
      } else {
        SeverityLogger log = GetLogger(boost::log::trivial::warning);
        BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ProbesOutput");
        BOOST_LOG(log) << fmt::format(
            "Path '{}' points to some non regular file. No probes output.",
            plenum_output_path_);
      }
    }
    ++i_block;
  }
}

} // namespace fub::amrex
