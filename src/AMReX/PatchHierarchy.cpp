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

#include "fub/AMReX/PatchHierarchy.hpp"
#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#include <boost/filesystem.hpp>
#include <hdf5.h>

#ifdef AMREX_USE_EB
#include "fub/AMReX/cutcell/AllRegularIndexSpace.hpp"
#include <AMReX_EB2.H>
#endif

namespace fub {
namespace amrex {

PatchLevel::PatchLevel(const PatchLevel& other)
    : level_number(other.level_number), time_point(other.time_point),
      cycles(other.cycles), box_array(other.box_array.boxList()),
      distribution_mapping(other.distribution_mapping.ProcessorMap()) {
  if (other.data.ok()) {
    data.define(box_array, distribution_mapping, other.data.nComp(),
                other.data.nGrowVect(), ::amrex::MFInfo(),
                other.data.Factory());
    data.copy(other.data);
  }
  if (other.nodes) {
    nodes = std::make_unique<::amrex::MultiFab>(
        box_array.surroundingNodes(), distribution_mapping,
        other.nodes->nComp(), other.nodes->nGrowVect(), ::amrex::MFInfo(),
        other.nodes->Factory());
    nodes->copy(*other.nodes);
  }
  if (other.faces[0]) {
    faces[0] = std::make_unique<::amrex::MultiFab>(
        box_array.convert({AMREX_D_DECL(1, 0, 0)}), distribution_mapping,
        other.faces[0]->nComp(), other.faces[0]->nGrowVect(), ::amrex::MFInfo(),
        other.faces[0]->Factory());
    faces[0]->copy(*other.faces[0]);
    if (other.faces[2]) {
      faces[2] = std::make_unique<::amrex::MultiFab>(
          box_array.convert({AMREX_D_DECL(0, 0, 1)}), distribution_mapping,
          other.faces[2]->nComp(), other.faces[2]->nGrowVect(),
          ::amrex::MFInfo(), other.faces[2]->Factory());
      faces[2]->copy(*other.faces[2]);
    }
    if (other.faces[1]) {
      faces[1] = std::make_unique<::amrex::MultiFab>(
          box_array.convert({AMREX_D_DECL(0, 1, 0)}), distribution_mapping,
          other.faces[1]->nComp(), other.faces[1]->nGrowVect(),
          ::amrex::MFInfo(), other.faces[1]->Factory());
      faces[1]->copy(*other.faces[1]);
    }
  }
}

PatchLevel& PatchLevel::operator=(const PatchLevel& other) {
  PatchLevel tmp(other);
  return *this = std::move(tmp);
}

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm, int n_components)
    : level_number{level}, time_point{tp}, box_array{ba},
      distribution_mapping{dm}, data{box_array, distribution_mapping,
                                     n_components, 0} {}

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm,
                       const DataDescription& desc)
    : level_number{level}, time_point{tp}, box_array{ba}, distribution_mapping{
                                                              dm} {
  if (desc.n_state_components) {
    data.define(box_array, distribution_mapping, desc.n_state_components, 0);
  }
  if (desc.n_node_components) {
    ::amrex::BoxArray on_nodes = box_array;
    on_nodes.surroundingNodes();
    nodes = std::make_unique<::amrex::MultiFab>(on_nodes,
                                                distribution_mapping,
                                                desc.n_node_components, 0);
  }
  if (desc.n_face_components) {
    ::amrex::BoxArray on_faces = box_array;
    switch (desc.dimension) {
    case 3:
      on_faces.convert({AMREX_D_DECL(0, 0, 1)});
      faces[2] = std::make_unique<::amrex::MultiFab>(
          on_faces, distribution_mapping,
          desc.n_face_components, 0);
      [[fallthrough]];
    case 2:
      on_faces.convert({AMREX_D_DECL(0, 1, 0)});
      faces[1] = std::make_unique<::amrex::MultiFab>(
          on_faces, distribution_mapping,
          desc.n_face_components, 0);
      [[fallthrough]];
    default:
      on_faces.convert({AMREX_D_DECL(1, 0, 0)});
      faces[0] = std::make_unique<::amrex::MultiFab>(
          on_faces, distribution_mapping,
          desc.n_face_components, 0);
    }
  }
}

PatchLevel::PatchLevel(int level, Duration tp, const ::amrex::BoxArray& ba,
                       const ::amrex::DistributionMapping& dm, int n_components,
                       const ::amrex::FabFactory<::amrex::FArrayBox>& factory)
    : level_number{level}, time_point{tp}, box_array{ba},
      distribution_mapping{dm}, data{box_array,         distribution_mapping,
                                     n_components,      0,
                                     ::amrex::MFInfo(), factory} {}

PatchHierarchy::PatchHierarchy(DataDescription desc,
                               const CartesianGridGeometry& geometry,
                               const PatchHierarchyOptions& options)
    : description_{std::move(desc)}, grid_geometry_{geometry},
      options_{options}, patch_level_{}, patch_level_geometry_{} {
  patch_level_.reserve(static_cast<std::size_t>(options.max_number_of_levels));
  patch_level_geometry_.resize(
      static_cast<std::size_t>(options.max_number_of_levels));
  ::amrex::IntVect lower{};
  ::amrex::IntVect upper{};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    upper[i] = geometry.cell_dimensions[static_cast<std::size_t>(i)] - 1;
  }
  ::amrex::Box level_box(lower, upper);
  for (::amrex::Geometry& geom : patch_level_geometry_) {
    geom = ::amrex::Geometry(level_box, &grid_geometry_.coordinates, -1,
                             grid_geometry_.periodicity.data());
    level_box.refine(options_.refine_ratio);
  }
#ifdef AMREX_USE_EB
  index_spaces_ =
      cutcell::BuildRegularSpace(patch_level_geometry_[0], options.refine_ratio,
                                 options.max_number_of_levels);
#endif
}

int PatchHierarchy::GetRatioToCoarserLevel(int level, Direction dir) const
    noexcept {
  if (level == 0) {
    return 1;
  }
  return options_.refine_ratio[static_cast<int>(dir)];
}

::amrex::IntVect PatchHierarchy::GetRatioToCoarserLevel(int level) const
    noexcept {
  if (level == 0) {
    return ::amrex::IntVect::TheUnitVector();
  }
  return options_.refine_ratio;
}

const ::amrex::Geometry& PatchHierarchy::GetGeometry(int level) const {
  return patch_level_geometry_[static_cast<std::size_t>(level)];
}

const PatchHierarchyOptions& PatchHierarchy::GetOptions() const noexcept {
  return options_;
}

const CartesianGridGeometry& PatchHierarchy::GetGridGeometry() const noexcept {
  return grid_geometry_;
}

std::ptrdiff_t PatchHierarchy::GetCycles(int level) const {
  return GetPatchLevel(level).cycles;
}

Duration PatchHierarchy::GetTimePoint(int level) const {
  return GetPatchLevel(level).time_point;
}

int PatchHierarchy::GetNumberOfLevels() const noexcept {
  return static_cast<int>(patch_level_.size());
}

int PatchHierarchy::GetMaxNumberOfLevels() const noexcept {
  return GetOptions().max_number_of_levels;
}

PatchLevel& PatchHierarchy::GetPatchLevel(int level) {
  return patch_level_.at(static_cast<std::size_t>(level));
}

span<const ::amrex::EB2::IndexSpace*>
PatchHierarchy::GetIndexSpaces() noexcept {
  return index_spaces_;
}

const PatchLevel& PatchHierarchy::GetPatchLevel(int level) const {
  return patch_level_[static_cast<std::size_t>(level)];
}

void PatchHierarchy::PushBack(const PatchLevel& level) {
  patch_level_.push_back(level);
}

void PatchHierarchy::PushBack(PatchLevel&& level) {
  patch_level_.push_back(std::move(level));
}

void PatchHierarchy::PopBack() { patch_level_.pop_back(); }

const DataDescription& PatchHierarchy::GetDataDescription() const noexcept {
  return description_;
}

void WriteCheckpointFile(const std::string checkpointname,
                         const fub::amrex::PatchHierarchy& hier) {
  const int nlevels = hier.GetNumberOfLevels();
  ::amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);
  // write Header file
  if (::amrex::ParallelDescriptor::IOProcessor()) {
    const std::string header(checkpointname + "/Header");
    ::amrex::VisMF::IO_Buffer io_buffer(::amrex::VisMF::IO_Buffer_Size);
    std::ofstream hout;
    hout.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    hout.open(header.c_str(), std::ofstream::out | std::ofstream::trunc |
                                  std::ofstream::binary);
    if (!hout.good()) {
      ::amrex::FileOpenFailed(header);
    }

    hout.precision(17);

    // write out title line
    hout << "Checkpoint File\n";

    // write out finest_level
    hout << nlevels - 1 << '\n';

    // write out array of istep
    for (int level = 0; level < nlevels; ++level) {
      hout << hier.GetCycles(level) << ' ';
    }
    hout << '\n';

    // write out array of t_new
    for (int level = 0; level < nlevels; ++level) {
      hout << hier.GetTimePoint(level).count() << ' ';
    }
    hout << '\n';

    // write the BoxArray at each level
    for (int level = 0; level < nlevels; ++level) {
      hier.GetPatchLevel(level).data.boxArray().writeOn(hout);
      hout << '\n';
    }
  }

  // write the MultiFab data to, e.g., chk00010/Level_0/
  for (int level = 0; level < nlevels; ++level) {
    ::amrex::VisMF::Write(hier.GetPatchLevel(level).data,
                          ::amrex::MultiFabFileFullPrefix(level, checkpointname,
                                                          "Level_", "data"));
  }
}

PatchHierarchy ReadCheckpointFile(const std::string& checkpointname,
                                  DataDescription desc,
                                  const CartesianGridGeometry& geometry,
                                  const PatchHierarchyOptions& options) {
  std::string File(checkpointname + "/Header");
  ::amrex::VisMF::IO_Buffer io_buffer(
      static_cast<std::size_t>(::amrex::VisMF::GetIOBufferSize()));
  ::amrex::Vector<char> fileCharPtr;
  ::amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);

  PatchHierarchy hierarchy(desc, geometry, options);

  // read in title line
  std::string line;
  std::getline(is, line);

  // read in finest_level
  std::getline(is, line);
  const int finest_level = std::stoi(line);
  std::vector<std::ptrdiff_t> cycles(
      static_cast<std::size_t>(finest_level + 1));
  std::vector<double> time_points(static_cast<std::size_t>(finest_level + 1));

  // read in array of istep
  std::getline(is, line);
  {
    std::istringstream lis(line);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(finest_level); ++i) {
      lis >> cycles[i];
    }
  }

  // read in array of t_new
  std::getline(is, line);
  {
    std::istringstream lis(line);
    for (std::size_t i = 0; i <= static_cast<std::size_t>(finest_level); ++i) {
      lis >> time_points[i];
    }
  }

  for (int lev = 0; lev <= finest_level; ++lev) {
    // read in level 'lev' BoxArray from Header
    ::amrex::BoxArray ba;
    ba.readFrom(is);

    // create a distribution mapping
    ::amrex::DistributionMapping dm{ba, ::amrex::ParallelDescriptor::NProcs()};

    hierarchy.PushBack(
        PatchLevel(lev, Duration(time_points[static_cast<std::size_t>(lev)]),
                   ba, dm, desc.n_state_components));
    hierarchy.GetPatchLevel(lev).cycles = cycles[static_cast<std::size_t>(lev)];
  }

  // read in the MultiFab data
  for (int lev = 0; lev <= finest_level; ++lev) {
    ::amrex::VisMF::Read(
        hierarchy.GetPatchLevel(lev).data,
        ::amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "data"));
  }

  return hierarchy;
}

namespace {
template <typename Deleter> struct H5Handle {
  H5Handle() = default;
  H5Handle(hid_t id) : id_{id} {}

  ~H5Handle() {
    if (id_ > 0) {
      Deleter{}(id_);
    }
  }

  operator hid_t() const noexcept { return id_; }

  hid_t id_{};
};

std::array<hsize_t, AMREX_SPACEDIM + 1>
GetExtents(const ::amrex::FArrayBox& fab) {
  std::array<hsize_t, AMREX_SPACEDIM + 1> extents = {hsize_t(fab.nComp())};
  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
    extents[size_t(i + 1)] = hsize_t(fab.box().length(i));
  }
  return extents;
}

struct H5Fdeleter {
  using pointer = hid_t;
  void operator()(hid_t file) const noexcept { H5Fclose(file); }
};
using H5File = H5Handle<H5Fdeleter>;

struct H5Sdeleter {
  using pointer = hid_t;
  void operator()(hid_t dataspace) const noexcept { H5Sclose(dataspace); }
};
using H5Space = H5Handle<H5Sdeleter>;

struct H5Ddeleter {
  using pointer = hid_t;
  void operator()(hid_t dataset) const noexcept { H5Dclose(dataset); }
};
using H5Dataset = H5Handle<H5Ddeleter>;

struct H5Adeleter {
  using pointer = hid_t;
  void operator()(hid_t attributes) const noexcept { H5Aclose(attributes); }
};
using H5Attribute = H5Handle<H5Adeleter>;

struct H5Pdeleter {
  using pointer = hid_t;
  void operator()(hid_t properties) const noexcept { H5Pclose(properties); }
};
using H5Properties = H5Handle<H5Pdeleter>;

void CreateHdf5Database(const std::string& name, const ::amrex::FArrayBox& fab,
                        const ::amrex::Geometry& geom, Duration tp,
                        std::ptrdiff_t cycle) {
  H5File file(H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  // Create Data Set
  {
    std::array<hsize_t, AMREX_SPACEDIM + 1> extents = GetExtents(fab);
    std::array<hsize_t, AMREX_SPACEDIM + 2> dims = {
        1, extents[0], AMREX_D_DECL(extents[1], extents[2], extents[3])};
    std::array<hsize_t, AMREX_SPACEDIM + 2> maxdims = {
        H5S_UNLIMITED, extents[0],
        AMREX_D_DECL(extents[1], extents[2], extents[3])};

    H5Properties properties(H5Pcreate(H5P_DATASET_CREATE));
    H5Pset_chunk(properties, AMREX_SPACEDIM + 2, dims.data());
    H5Pset_alloc_time(properties, H5D_ALLOC_TIME_EARLY);
    H5Space dataspace(
        H5Screate_simple(AMREX_SPACEDIM + 2, dims.data(), maxdims.data()));
    H5Dataset dataset(H5Dcreate(file, "/data", H5T_IEEE_F64LE, dataspace,
                                H5P_DEFAULT, properties, H5P_DEFAULT));
    H5Dwrite(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             fab.dataPtr());
    hsize_t spacedim_value = AMREX_SPACEDIM;
    H5Space spacedim(H5Screate_simple(1, &spacedim_value, nullptr));

    H5Attribute xlower(H5Acreate2(dataset, "xlower", H5T_IEEE_F64LE, spacedim,
                                  H5P_DEFAULT, H5P_DEFAULT));
    H5Awrite(xlower, H5T_IEEE_F64LE, geom.ProbLo());

    H5Attribute cell_size(H5Acreate2(dataset, "cell_size", H5T_IEEE_F64LE,
                                     spacedim, H5P_DEFAULT, H5P_DEFAULT));
    H5Awrite(cell_size, H5T_IEEE_F64LE, geom.CellSize());
  }
  hsize_t dims[1] = {1};
  hsize_t maxdims[1] = {H5S_UNLIMITED};
  H5Properties properties(H5Pcreate(H5P_DATASET_CREATE));
  H5Pset_chunk(properties, 1, dims);
  H5Space dataspace(H5Screate_simple(1, dims, maxdims));

  H5Dataset times(H5Dcreate(file, "/times", H5T_IEEE_F64LE, dataspace,
                            H5P_DEFAULT, properties, H5P_DEFAULT));
  const double count = tp.count();
  H5Dwrite(times, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &count);
  H5Dataset cycles(H5Dcreate(file, "/cycles", H5T_STD_I64LE_g, dataspace,
                             H5P_DEFAULT, properties, H5P_DEFAULT));
  const std::int64_t cycle_count = static_cast<std::int64_t>(cycle);
  H5Dwrite(cycles, H5T_STD_I64LE_g, H5S_ALL, H5S_ALL, H5P_DEFAULT,
           &cycle_count);
}

void OpenHdf5Database(const std::string& name, const ::amrex::FArrayBox& fab,
                      Duration tp, std::ptrdiff_t cycle) {
  // Open or Create Data Set
  H5File file(H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  if (file < 0) {
    return;
  }
  if (H5Dataset dataset(H5Dopen(file, "/data", H5P_DEFAULT)); dataset < 0) {
    return;
  } else {
    auto map_first = [](const std::array<hsize_t, AMREX_SPACEDIM + 2>& x,
                        auto f) {
      return std::array<hsize_t, AMREX_SPACEDIM + 2>{
          hsize_t(f(x[0])), x[1], AMREX_D_DECL(x[2], x[3], x[4])};
    };
    std::array<hsize_t, AMREX_SPACEDIM + 2> dims = {};
    std::array<hsize_t, AMREX_SPACEDIM + 2> maxdims = {};
    H5Space dataspace(H5Dget_space(dataset));
    H5Sget_simple_extent_dims(dataspace, dims.data(), maxdims.data());
    FUB_ASSERT(maxdims[0] == H5S_UNLIMITED);
    std::array<hsize_t, AMREX_SPACEDIM + 2> new_dims =
        map_first(dims, [](hsize_t dim) { return dim + 1; });
    H5Dset_extent(dataset, new_dims.data());
    H5Space filespace(H5Dget_space(dataset));
    std::array<hsize_t, AMREX_SPACEDIM + 2> count =
        map_first(dims, [](hsize_t) { return 1; });
    std::array<hsize_t, AMREX_SPACEDIM + 2> offset{dims[0]};
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), nullptr,
                        count.data(), nullptr);
    H5Space memspace(
        H5Screate_simple(AMREX_SPACEDIM + 2, count.data(), nullptr));
    H5Dwrite(dataset, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT,
             fab.dataPtr());
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
    const double tp_count = tp.count();
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
    const std::int64_t cycle_count = static_cast<std::int64_t>(cycle);
    H5Space memspace(H5Screate_simple(1, &count, nullptr));
    H5Dwrite(cycles, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT,
             &cycle_count);
  }
}
} // namespace

void WriteToHDF5(const std::string& name, const ::amrex::FArrayBox& fab,
                 const ::amrex::Geometry& geom, Duration time_point,
                 std::ptrdiff_t cycle) noexcept {
  if (!boost::filesystem::exists(name)) {
    CreateHdf5Database(name, fab, geom, time_point, cycle);
  } else if (boost::filesystem::is_regular_file(name)) {
    OpenHdf5Database(name, fab, time_point, cycle);
  }
}

void WriteMatlabData(std::ostream& out, const ::amrex::FArrayBox& fab,
                     const fub::IdealGasMix<1>& eq,
                     const ::amrex::Geometry& geom) {
  auto view = MakeView<const Complete<IdealGasMix<1>>>(fab, eq);
  out << fmt::format("{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}", "X",
                     "Density", "VelocityX", "SpeedOfSound", "Temperature",
                     "Pressure", "Gamma", "HeatCapacityAtConstantPressure");
  const int nspecies = eq.GetReactor().GetNSpecies();
  span<const std::string> names = eq.GetReactor().GetSpeciesNames();
  for (int s = 0; s < nspecies - 1; ++s) {
    out << names[s] << ' ';
  }
  out << names[nspecies - 1] << '\n';
  ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i) {
    double x[3] = {0.0, 0.0, 0.0};
    ::amrex::IntVect iv{AMREX_D_DECL(int(i), 0, 0)};
    geom.CellCenter(iv, x);
    const double density = view.density(i);
    const double velocity_x =
        density > 0.0 ? view.momentum(i, 0) / density : 0.0;
    const double temperature = view.temperature(i);
    const double speed_of_sound = view.speed_of_sound(i);
    const double gamma = view.gamma(i);
    const double c_p = view.c_p(i);
    const double pressure = view.pressure(i);
    out << fmt::format("{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< "
                       "24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}",
                       x[0], density, velocity_x, speed_of_sound, temperature,
                       pressure, gamma, c_p);
    for (int s = 0; s < nspecies; ++s) {
      out << fmt::format("{:< 24.15e}", view.species(i, s));
    }
    out << '\n';
  });
}

#if AMREX_SPACEDIM == 3
void WriteMatlabData(std::ostream& out, const ::amrex::FArrayBox& fab,
                     const fub::IdealGasMix<3>& eq,
                     const ::amrex::Geometry& geom) {
  auto view = MakeView<const Complete<IdealGasMix<3>>>(fab, eq);
  out << fmt::format("{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:24s}{:"
                     "24s}{:24s}{:24s}",
                     "X", "Y", "Z", "Density", "VelocityX", "VelocityY",
                     "VelocityZ", "SpeedOfSound", "Temperature", "Pressure",
                     "Gamma", "HeatCapacityAtConstantPressure");
  const int nspecies = eq.GetReactor().GetNSpecies();
  span<const std::string> names = eq.GetReactor().GetSpeciesNames();
  for (int s = 0; s < nspecies - 1; ++s) {
    out << fmt::format("{:24s}", names[s]);
  }
  out << fmt::format("{:24s}\n", names[nspecies - 1]);
  ForEachIndex(Box<0>(view), [&](std::ptrdiff_t i, std::ptrdiff_t j,
                                 std::ptrdiff_t k) {
    double x[3] = {0.0, 0.0, 0.0};
    ::amrex::IntVect iv{int(i), int(j), int(k)};
    geom.CellCenter(iv, x);
    const double density = view.density(i, j, k);
    const double velocity_x =
        density > 0.0 ? view.momentum(i, j, k, 0) / density : 0.0;
    const double velocity_y =
        density > 0.0 ? view.momentum(i, j, k, 1) / density : 0.0;
    const double velocity_z =
        density > 0.0 ? view.momentum(i, j, k, 2) / density : 0.0;
    const double temperature = view.temperature(i, j, k);
    const double speed_of_sound = view.speed_of_sound(i, j, k);
    const double gamma = view.gamma(i, j, k);
    const double c_p = view.c_p(i, j, k);
    const double pressure = view.pressure(i, j, k);
    out << fmt::format("{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}"
                       "{:< 24.15e}{:< 24.15e}{:< 24.15e}{:< 24.15e}"
                       "{:< 24.15e}{:< 24.15e}{:< 24.15e}",
                       x[0], x[1], x[2], density, velocity_x, velocity_y,
                       velocity_z, speed_of_sound, temperature, pressure, gamma,
                       c_p);
    for (int s = 0; s < nspecies; ++s) {
      out << fmt::format("{:< 24.15e}", view.species(i, j, k, s));
    }
    out << '\n';
  });
}
#endif

std::vector<double> GatherStates(
    const PatchHierarchy& hierarchy,
    basic_mdspan<const double, extents<AMREX_SPACEDIM, dynamic_extent>> xs,
    MPI_Comm comm) {
  const int nlevel = hierarchy.GetNumberOfLevels();
  const int finest_level = nlevel - 1;
  const int ncomp = hierarchy.GetDataDescription().n_state_components;
  std::vector<double> buffer(xs.extent(1) * ncomp * nlevel);
  mdspan<double, 3> states(buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 0; level < nlevel; ++level) {
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(level).data;
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(level);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      ForEachIndex(mfi.tilebox(), [&](auto... is) {
        double lo[AMREX_SPACEDIM]{};
        double hi[AMREX_SPACEDIM]{};
        const ::amrex::IntVect iv{int(is)...};
        level_geom.LoNode(iv, lo);
        level_geom.HiNode(iv, hi);
        for (int k = 0; k < xs.extent(1); ++k) {
          bool is_in_range = true;
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            is_in_range = is_in_range && lo[d] <= xs(d, k) && xs(d, k) < hi[d];
          }
          if (is_in_range) {
            for (int comp = 0; comp < level_data.nComp(); ++comp) {
              states(k, comp, level) = level_data[mfi](iv, comp);
            }
          }
        }
      });
    });
  }
  std::vector<double> global_buffer(buffer.size());
  ::MPI_Allreduce(buffer.data(), global_buffer.data(), global_buffer.size(),
                  MPI_DOUBLE, MPI_SUM, comm);
  states = mdspan<double, 3>(global_buffer.data(), xs.extent(1), ncomp, nlevel);
  for (int level = 1; level < nlevel; ++level) {
    for (int comp = 0; comp < ncomp; ++comp) {
      for (int i = 0; i < xs.extent(1); ++i) {
        if (states(i, comp, level) == 0.0) {
          states(i, comp, level) = states(i, comp, level - 1);
        }
      }
    }
  }
  std::vector<double> result(&states(0, 0, finest_level),
                             &states(0, 0, finest_level) +
                                 xs.extent(1) * ncomp);
  return result;
}

void WriteTubeData(const std::string& name, const PatchHierarchy& hierarchy,
                   const IdealGasMix<1>&, fub::Duration time_point,
                   std::ptrdiff_t cycle_number, MPI_Comm comm) {
  int rank = -1;
  ::MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    boost::filesystem::path path(name);
    boost::filesystem::path dir = path.parent_path();
    boost::filesystem::create_directories(dir);
  }
  const std::size_t n_level =
      static_cast<std::size_t>(hierarchy.GetNumberOfLevels());
  std::vector<::amrex::FArrayBox> fabs{};
  fabs.reserve(n_level);
  for (std::size_t level = 0; level < n_level; ++level) {
    const int ilvl = static_cast<int>(level);
    const ::amrex::Geometry& level_geom = hierarchy.GetGeometry(ilvl);
    ::amrex::Box domain = level_geom.Domain();
    const ::amrex::MultiFab& level_data = hierarchy.GetPatchLevel(ilvl).data;
    ::amrex::FArrayBox local_fab(domain, level_data.nComp());
    local_fab.setVal(0.0);
    ForEachFab(level_data, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::FArrayBox& patch_data = level_data[mfi];
      local_fab.copy(patch_data);
    });
    if (rank == 0) {
      ::amrex::FArrayBox& fab = fabs.emplace_back(domain, level_data.nComp());
      fab.setVal(0.0);
      ::MPI_Reduce(local_fab.dataPtr(), fab.dataPtr(), local_fab.size(),
                   MPI_DOUBLE, MPI_SUM, 0, comm);
      if (level > 0) {
        for (int comp = 1; comp < level_data.nComp(); ++comp) {
          for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
            ::amrex::IntVect fine_i{
                AMREX_D_DECL(i, domain.smallEnd(1), domain.smallEnd(2))};
            ::amrex::IntVect coarse_i = fine_i;
            coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
            if (fab(fine_i, 0) == 0.0) {
              fab(fine_i, comp) = fabs[level - 1](coarse_i, comp);
            }
          }
        }
        for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
          ::amrex::IntVect fine_i{
              AMREX_D_DECL(i, domain.smallEnd(1), domain.smallEnd(2))};
          ::amrex::IntVect coarse_i = fine_i;
          coarse_i.coarsen(hierarchy.GetRatioToCoarserLevel(level));
          if (fab(fine_i, 0) == 0.0) {
            fab(fine_i, 0) = fabs[level - 1](coarse_i, 0);
          }
        }
      }
      if (level == n_level - 1) {
        WriteToHDF5(name, fab, level_geom, time_point, cycle_number);
      }
    } else {
      ::MPI_Reduce(local_fab.dataPtr(), nullptr, local_fab.size(), MPI_DOUBLE,
                   MPI_SUM, 0, comm);
    }
  }
}

} // namespace amrex
} // namespace fub
