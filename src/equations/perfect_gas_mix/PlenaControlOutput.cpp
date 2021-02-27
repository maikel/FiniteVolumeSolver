// Copyright (c) 2021 Maikel Nadolski
// Copyright (c) 2021 Christian Zenker
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

#include "fub/equations/perfect_gas_mix/PlenaControlOutput.hpp"

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/map.hpp>

#include "fub/output/Hdf5Handle.hpp"

namespace fub::perfect_gas_mix::gt {

namespace {
auto GetFieldMap() {
  using namespace std::literals;
  using projection = function_ref<double(const ControlState&)>;
  return std::map<std::string, projection>{
      std::pair<std::string, projection>{
          "current_rpm"s, [](const ControlState& s) { return s.current_rpm; }},
      std::pair<std::string, projection>{
          "power_out"s, [](const ControlState& s) { return s.power_out; }},
      std::pair<std::string, projection>{
          "fuel_consumption"s,
          [](const ControlState& s) { return s.fuel_consumption; }},
      std::pair<std::string, projection>{
          "fuel_consumption_rate"s,
          [](const ControlState& s) { return s.fuel_consumption_rate; }},
      std::pair<std::string, projection>{
          "efficiency"s, [](const ControlState& s) { return s.efficiency; }},
      std::pair<std::string, projection>{
          "compressor_pressure"s,
          [](const ControlState& s) { return s.compressor.pressure; }},
      std::pair<std::string, projection>{
          "compressor_temperature"s,
          [](const ControlState& s) { return s.compressor.temperature; }},
      std::pair<std::string, projection>{
          "compressor_power"s,
          [](const ControlState& s) { return s.compressor.power; }},
      std::pair<std::string, projection>{
          "compressor_mass_flow_in"s,
          [](const ControlState& s) { return s.compressor.mass_flow_in; }},
      std::pair<std::string, projection>{
          "compressor_mass_flow_out"s,
          [](const ControlState& s) { return s.compressor.mass_flow_out; }},
      std::pair<std::string, projection>{"compressor_SEC_Mode"s,
                                         [](const ControlState& s) {
                                           return static_cast<double>(
                                               s.compressor.SEC_Mode);
                                         }},
      std::pair<std::string, projection>{
          "turbine_pressure"s,
          [](const ControlState& s) { return s.turbine.pressure; }},
      std::pair<std::string, projection>{
          "turbine_temperature"s,
          [](const ControlState& s) { return s.turbine.temperature; }},
      std::pair<std::string, projection>{
          "turbine_power"s,
          [](const ControlState& s) { return s.turbine.power; }},
      std::pair<std::string, projection>{
          "turbine_mass_flow_in"s,
          [](const ControlState& s) { return s.turbine.mass_flow_in; }},
      std::pair<std::string, projection>{
          "turbine_mass_flow_out"s,
          [](const ControlState& s) { return s.turbine.mass_flow_out; }}};
}
} // namespace

void ControlOutput::CreateHdf5Database() {
  H5File file_ =
      H5Fcreate(file_path_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  // Create Data Set
  {
    const hsize_t field_size = static_cast<hsize_t>(fields_.size());
    std::array<hsize_t, 2> dims{1, field_size};
    std::array<hsize_t, 2> maxdims{H5S_UNLIMITED, field_size};
    // write each time step immediately
    H5Properties properties(H5Pcreate(H5P_DATASET_CREATE));
    H5Pset_chunk(properties, dims.size(), dims.data());
    H5Pset_alloc_time(properties, H5D_ALLOC_TIME_EARLY);
    // create an empty dataset
    dims[0] = 0;
    H5Space data_dataspace_ =
        H5Screate_simple(dims.size(), dims.data(), maxdims.data());
    H5Dataset data_dataset_ =
        H5Dcreate(file_, "/data", H5T_IEEE_F64LE, data_dataspace_, H5P_DEFAULT,
                  properties, H5P_DEFAULT);
  }
  // Create Field names
  if (!fields_.empty()) {
    std::array<hsize_t, 1> dims = {static_cast<hsize_t>(fields_.size())};
    H5Space fields_dataspace(
        H5Screate_simple(dims.size(), dims.data(), nullptr));
    H5Type datatype(H5Tcopy(H5T_C_S1));
    H5Tset_size(datatype, H5T_VARIABLE);
    H5Dataset fields_dataset(H5Dcreate2(file_, "/fields", datatype,
                                        fields_dataspace, H5P_DEFAULT,
                                        H5P_DEFAULT, H5P_DEFAULT));
    hsize_t count[] = {1};
    H5Space memspace(H5Screate_simple(dims.size(), count, nullptr));
    // Now we write each field name into the fields_dataset
    for (auto&& [i, field] :
         ranges::view::enumerate(ranges::view::keys(fields_))) {
      // H5Space filespace(H5Dget_space(fields_dataset));
      hsize_t offset[] = {i};
      H5Sselect_hyperslab(fields_dataspace, H5S_SELECT_SET, offset, nullptr,
                          count, nullptr);
      const char* c_str = field.c_str();
      H5Dwrite(fields_dataset, datatype, memspace, fields_dataspace,
               H5P_DEFAULT, &c_str);
    }
  }
  // Create times and cycles datasets
  std::array<hsize_t, 1> dims = {1};
  std::array<hsize_t, 1> maxdims = {H5S_UNLIMITED};
  H5Properties properties(H5Pcreate(H5P_DATASET_CREATE));
  H5Pset_chunk(properties, dims.size(), dims.data());
  H5Pset_alloc_time(properties, H5D_ALLOC_TIME_EARLY);
  dims[0] = 0;
  H5Space times_dataspace_ =
      H5Screate_simple(dims.size(), dims.data(), maxdims.data());
  H5Dataset times_dataset_ =
      H5Dcreate(file_, "/times", H5T_IEEE_F64LE, times_dataspace_, H5P_DEFAULT,
                properties, H5P_DEFAULT);
  H5Space cycles_dataspace_ =
      H5Screate_simple(dims.size(), dims.data(), maxdims.data());
  H5Dataset cycles_dataset_ =
      H5Dcreate(file_, "/cycles", H5T_STD_I64LE_g, cycles_dataspace_,
                H5P_DEFAULT, properties, H5P_DEFAULT);
}

void ControlOutput::WriteHdf5Database(span<const double> data, Duration time,
                                      std::ptrdiff_t cycle) {
  // Write to /data dataset
  H5File file_ = H5Fopen(file_path_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_ > 0) {
    std::array<hsize_t, 2> dims{};
    std::array<hsize_t, 2> maxdims{};
    H5Dataset data_dataset_ = H5Dopen(file_, "/data", H5P_DEFAULT);
    if (data_dataset_ > 0) {
      H5Space data_dataspace_ = H5Dget_space(data_dataset_);
      H5Sget_simple_extent_dims(data_dataspace_, dims.data(), maxdims.data());
      FUB_ASSERT(maxdims[0] == H5S_UNLIMITED);
      std::array<hsize_t, 2> new_dims{dims[0] + 1, dims[1]};
      H5Dset_extent(data_dataset_, new_dims.data());
      H5Space filespace(H5Dget_space(data_dataset_));
      std::array<hsize_t, 2> count{1, dims[1]};
      std::array<hsize_t, 2> offset{dims[0], 0};
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), nullptr,
                          count.data(), nullptr);
      H5Space memspace(H5Screate_simple(2, count.data(), nullptr));
      H5Dwrite(data_dataset_, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT,
               data.data());
    }
  }
  // Write to /times dataset
  if (file_ > 0) {
    hsize_t dims{};
    hsize_t maxdims{};
    H5Dataset times_dataset_ = H5Dopen(file_, "/times", H5P_DEFAULT);
    if (times_dataset_ > 0) {
      H5Space times_dataspace_ = H5Dget_space(times_dataset_);
      H5Sget_simple_extent_dims(times_dataspace_, &dims, &maxdims);
      FUB_ASSERT(maxdims == H5S_UNLIMITED);
      const hsize_t new_dims = dims + 1;
      H5Dset_extent(times_dataset_, &new_dims);
      H5Space filespace(H5Dget_space(times_dataset_));
      const hsize_t count = 1;
      const hsize_t offset = dims;
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, &count, &count,
                          &count);
      const double tp_count = time.count();
      H5Space memspace(H5Screate_simple(1, &count, nullptr));
      H5Dwrite(times_dataset_, H5T_IEEE_F64LE, memspace, filespace, H5P_DEFAULT,
               &tp_count);
    }
  }
  // Write to /cycles dataset
  if (file_ > 0) {
    hsize_t dims{};
    hsize_t maxdims{};
    H5Dataset cycles_dataset_ = H5Dopen(file_, "/cycles", H5P_DEFAULT);
    if (cycles_dataset_ > 0) {
      H5Space cycles_dataspace_ = H5Dget_space(cycles_dataset_);
      H5Sget_simple_extent_dims(cycles_dataspace_, &dims, &maxdims);
      FUB_ASSERT(maxdims == H5S_UNLIMITED);
      const hsize_t new_dims = dims + 1;
      H5Dset_extent(cycles_dataset_, &new_dims);
      cycles_dataspace_ = H5Dget_space(cycles_dataset_);
      H5Space filespace(H5Dget_space(cycles_dataset_));
      const hsize_t count = 1;
      const hsize_t offset = dims;
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr, &count,
                          nullptr);
      const std::int64_t cycle_count = static_cast<std::int64_t>(cycle);
      H5Space memspace(H5Screate_simple(1, &count, nullptr));
      H5Dwrite(cycles_dataset_, H5T_STD_I64LE_g, memspace, filespace,
               H5P_DEFAULT, &cycle_count);
    }
  }
}

ControlOutput::ControlOutput(const ProgramOptions& options,
                             std::shared_ptr<const ControlState> control_state)
    : OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2>(options),
      control_state_(std::move(control_state)), fields_{GetFieldMap()},
      data_buffer_(fields_.size()) {
  file_path_ = GetOptionOr(options, "path", std::string("ControlState.h5"));
  fub::SeverityLogger log = GetInfoLogger();
  BOOST_LOG(log) << "ControlOutput configured:";
  BOOST_LOG(log) << fmt::format(" - path: {}", file_path_);
  OutputAtFrequencyOrInterval<amrex::MultiBlockGriddingAlgorithm2>::Print(log);
  int rank = -1;
  MPI_Comm_rank(::amrex::ParallelDescriptor::Communicator(), &rank);
  if (rank == 0) {
    boost::filesystem::path path(file_path_);
    boost::filesystem::path dir = boost::filesystem::absolute(
        path.parent_path(), boost::filesystem::current_path());
    if (!boost::filesystem::exists(dir)) {
      boost::filesystem::create_directories(dir);
    }
    if (boost::filesystem::is_regular_file(path)) {
      for (int i = 1; i < std::numeric_limits<int>::max(); ++i) {
        std::string new_name = fmt::format("{}.{}", file_path_, i);
        if (!boost::filesystem::exists(new_name)) {
          BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ControlOutput");
          BOOST_LOG(log) << fmt::format(
              "Old output file '{}' detected. Rename old file to '{}'",
              file_path_, new_name);
          boost::filesystem::rename(file_path_, new_name);
          break;
        }
      }
    } else if (boost::filesystem::exists(path)) {
      boost::log::sources::severity_logger<boost::log::trivial::severity_level>
          log(boost::log::keywords::severity = boost::log::trivial::warning);
      BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ControlOutput");
      for (int i = 1; i < std::numeric_limits<int>::max(); ++i) {
        std::string new_name = fmt::format("{}.{}", file_path_, i);
        if (!boost::filesystem::exists(new_name)) {
          BOOST_LOG(log) << fmt::format(
              "Path'{}' points to some non-file. Output will be directory to "
              "'{}' instead.",
              file_path_, new_name);
          file_path_ = new_name;
          break;
        }
      }
    }
    CreateHdf5Database();
  }
}

void ControlOutput::operator()(
    const amrex::MultiBlockGriddingAlgorithm2& grid) {
  int rank = -1;
  MPI_Comm_rank(::amrex::ParallelDescriptor::Communicator(), &rank);
  boost::log::sources::severity_logger<boost::log::trivial::severity_level> log(
      boost::log::keywords::severity = boost::log::trivial::info);
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Channel", "ControlOutput");
  BOOST_LOG_SCOPED_LOGGER_TAG(log, "Time", grid.GetTimePoint().count());
  BOOST_LOG(log) << fmt::format("Write Hdf5 output to '{}'.", file_path_);
  if (rank == 0) {
    auto values = ranges::view::values(fields_) |
                  ranges::view::transform(
                      [&](auto&& a) -> double { return a(*control_state_); });
    ranges::copy(values, data_buffer_.begin());
    WriteHdf5Database(data_buffer_, grid.GetTimePoint(), grid.GetCycles());
  }
}

} // namespace fub::perfect_gas_mix::gt
