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

#include "fub/equations/perfect_gas_mix/PlenaControl.hpp"

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/map.hpp>

namespace fub::perfect_gas_mix::gt {

namespace {
/// Notify the control with a new turbine plenum value and update its current
/// power.
///
/// The formula is given in GT_Notes equation (35).
///
/// \return Returns the new power value for the turbine
void UpdateTurbinePlenum(ControlState& state, const PerfectGasConstants& eq,
                         const ControlOptions& options, double mdot,
                         const PlenumState& turbine_plenum_current,
                         double flux_rho_last, Duration dt) noexcept {
  const double pressure_half =
      0.5 * (state.turbine.pressure + turbine_plenum_current.pressure);
  const double T_half =
      0.5 * (state.turbine.temperature + turbine_plenum_current.temperature);
  const double gmoog = eq.gamma_minus_one_over_gamma;
  const double cp = eq.heat_capacity_at_constant_pressure;
  const double ootau = 1.0;
  double dt_over_tau = dt.count() * ootau;
  state.turbine.mass_flow_in =
      (1.0 - dt_over_tau) * state.turbine.mass_flow_in +
      dt_over_tau * flux_rho_last * options.surface_area_tube_outlet;
  state.turbine.mass_flow_out =
      (1.0 - dt_over_tau) * state.turbine.mass_flow_out + dt_over_tau * mdot;
  state.turbine.power = mdot * cp * options.efficiency_turbine * T_half *
                        (1.0 - std::pow(pressure_half, -gmoog));
  state.turbine.pressure = turbine_plenum_current.pressure;
  state.turbine.temperature = turbine_plenum_current.temperature;
}

/// The formula is given in GT_Notes equation (20).
///
/// \return Returns the compressure pressure ratio
double CompressorPressureRatio(const ControlOptions& options,
                               double rpm) noexcept {
  double pratio = options.pratiomean +
                  options.c_0 * options.pratiovar *
                      std::atan(M_PI * (rpm / options.rpmmean - 1.0)) / M_PI;
  return pratio;
}

/// \brief Compute the compressor state by given speed of the compressor.
void UpdateCompressorFromRPM(ControlState& state, const PerfectGasConstants& eq,
                             const ControlOptions& options) {
  double new_pressure = CompressorPressureRatio(options, state.current_rpm);
  state.compressor.pressure = new_pressure;
  const double gmoog = eq.gamma_minus_one_over_gamma;

  // see GT_Notes equation 21
  state.compressor.temperature = 1.0 + (std::pow(new_pressure, gmoog) - 1.0) /
                                           options.efficiency_compressor;
}

/// Notify the control with a new compressor plenum value and update its
/// current power.
///
/// The formula is given in GT_Notes equation (24).
///
/// \return Returns the new power value
void UpdateCompressorPlenum(ControlState& state, const PerfectGasConstants& eq,
                            const ControlOptions& options, double flux_rho,
                            double flux_spec, Duration dt) noexcept {
  double temperature_old = state.compressor.temperature;
  double rho_old =
      state.compressor.pressure / state.compressor.temperature * eq.ooRspec;

  // update the compressor_plenum state with the new speed of the compressor
  // so we get new pressure and temperature
  UpdateCompressorFromRPM(state, eq, options);

  double rho_new =
      state.compressor.pressure / state.compressor.temperature * eq.ooRspec;

  const double mdot = flux_rho * options.surface_area_tube_inlet +
                      options.volume_compressor_plenum / dt.count() *
                          (rho_new - rho_old); // equation (23) GT_Notes

  double T_v_half = 0.5 * (temperature_old + state.compressor.temperature);
  double T_ref = 1.0;
  const double cp = eq.heat_capacity_at_constant_pressure;

  double power_compressor_increase =
      mdot * (T_v_half - T_ref) * cp / options.efficiency_compressor;

  const double ootau = 1.0;
  double dt_over_tau = dt.count() * ootau;

  state.compressor.mass_flow_in =
      (1.0 - dt_over_tau) * state.compressor.mass_flow_in + dt_over_tau * mdot;

  state.compressor.mass_flow_out =
      (1.0 - dt_over_tau) * state.compressor.mass_flow_out +
      dt_over_tau * flux_rho * options.surface_area_tube_inlet;

  state.compressor.power = (1.0 - dt_over_tau) * state.compressor.power +
                           dt_over_tau * power_compressor_increase;
  state.fuel_consumption +=
      flux_spec * options.surface_area_tube_inlet * dt.count();

  state.fuel_consumption_rate =
      (1.0 - dt_over_tau) * state.fuel_consumption_rate +
      dt_over_tau * flux_spec * options.surface_area_tube_inlet;
}

/// \brief Compute the new rotation speed of the compressor.
///
/// see GT_notes equation (42 - 44)
void ChangeRPM(ControlState& state, const ControlOptions& options,
               Duration dt) noexcept {
  double pressure = state.compressor.pressure;
  const double pressure_ratio =
      (options.target_pressure_compressor - pressure) /
      options.target_pressure_compressor;
  const double exp_term = std::exp(-options.mu * pressure_ratio);
  const double comp_rate = 0.015 * (1.0 - exp_term);
  constexpr double pi2 = M_PI * M_PI;
  const double Ieff = 2.0 * pi2 * options.inertial_moment;

  const double power_netto = state.turbine.power - state.compressor.power;
  const double d_Erot_dt = (pressure < options.target_pressure_compressor)
                               ? comp_rate * power_netto
                               : 0.0;
  const double Erot =
      Ieff * state.current_rpm * state.current_rpm + d_Erot_dt * dt.count();
  const double rpm_new =
      std::clamp(std::sqrt(Erot / Ieff), options.rpmmin, options.rpmmax);
  // FUB_ASSERT(options.rpmmin <= rpm_new && options.rpmmax >= rpm_new);
  state.current_rpm = rpm_new;

  // values to track:
  state.power_out = power_netto - d_Erot_dt; // Power output

  // efficiency of the machine
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  state.efficiency = std::max(0.0, state.power_out) /
                     (options.Q * state.fuel_consumption_rate + eps);

  if (!state.compressor.SEC_Mode) {
    // 0.9995 just shortcut because floating point comparison
    double target_pressure = 0.9995 * options.target_pressure_compressor;
    state.compressor.SEC_Mode = pressure >= target_pressure;
    // state.compressor.SEC_Mode = (pressure < target_pressure) ? false : true;
  } else {
    double fall_back_pressure = 0.9 * options.target_pressure_compressor;
    state.compressor.SEC_Mode = pressure > fall_back_pressure;
    // state.compressor.SEC_Mode = (pressure > fall_back_pressure) ? true : false;
  }
}
} // namespace

#ifdef FUB_GT_CONTROL_GET_OPTION
#undef FUB_GT_CONTROL_GET_OPTION
#endif
#define FUB_GT_CONTROL_GET_OPTION(var) var = GetOptionOr(options, #var, var)
ControlOptions::ControlOptions(const ProgramOptions& options) {
  FUB_GT_CONTROL_GET_OPTION(efficiency_turbine);
  FUB_GT_CONTROL_GET_OPTION(length_tube);
  FUB_GT_CONTROL_GET_OPTION(surface_area_tube_inlet);
  FUB_GT_CONTROL_GET_OPTION(surface_area_tube_outlet);
  FUB_GT_CONTROL_GET_OPTION(surface_area_compressor_to_compressor_plenum);
  FUB_GT_CONTROL_GET_OPTION(volume_compressor_plenum);
  FUB_GT_CONTROL_GET_OPTION(surface_area_turbine_plenum_to_turbine);
  FUB_GT_CONTROL_GET_OPTION(volume_turbine_plenum);
  FUB_GT_CONTROL_GET_OPTION(efficiency_compressor);
  FUB_GT_CONTROL_GET_OPTION(rpmmin);
  FUB_GT_CONTROL_GET_OPTION(rpmmax);
  FUB_GT_CONTROL_GET_OPTION(pratiomin);
  FUB_GT_CONTROL_GET_OPTION(pratiomax);
  FUB_GT_CONTROL_GET_OPTION(c_0);
  FUB_GT_CONTROL_GET_OPTION(inertial_moment);
  FUB_GT_CONTROL_GET_OPTION(mu);
  FUB_GT_CONTROL_GET_OPTION(target_pressure_compressor);
  FUB_GT_CONTROL_GET_OPTION(Q);
  initial_turbine_state.pressure = GetOptionOr(
      options, "initial_turbine_pressure", initial_turbine_state.pressure);
  initial_turbine_state.temperature =
      GetOptionOr(options, "initial_turbine_temperature",
                  initial_turbine_state.temperature);

  rpmmean = 0.5 * (rpmmax + rpmmin);
  pratiomean = 0.5 * (pratiomax + pratiomin);
  pratiovar = pratiomax - pratiomin;
}
#undef FUB_GT_CONTROL_GET_OPTION

#ifdef FUB_GT_CONTROL_PRINT
#undef FUB_GT_CONTROL_PRINT
#endif
#define FUB_GT_CONTROL_PRINT(variable)                                         \
  BOOST_LOG(log) << " - " #variable " = " << variable;

void ControlOptions::Print(SeverityLogger& log) {
  FUB_GT_CONTROL_PRINT(efficiency_turbine);
  FUB_GT_CONTROL_PRINT(length_tube);
  FUB_GT_CONTROL_PRINT(surface_area_tube_inlet);
  FUB_GT_CONTROL_PRINT(surface_area_tube_outlet);
  FUB_GT_CONTROL_PRINT(surface_area_compressor_to_compressor_plenum);
  FUB_GT_CONTROL_PRINT(volume_compressor_plenum);
  FUB_GT_CONTROL_PRINT(surface_area_turbine_plenum_to_turbine);
  FUB_GT_CONTROL_PRINT(volume_turbine_plenum);
  FUB_GT_CONTROL_PRINT(efficiency_compressor);
  FUB_GT_CONTROL_PRINT(rpmmin);
  FUB_GT_CONTROL_PRINT(rpmmax);
  FUB_GT_CONTROL_PRINT(pratiomin);
  FUB_GT_CONTROL_PRINT(pratiomax);
  FUB_GT_CONTROL_PRINT(c_0);
  FUB_GT_CONTROL_PRINT(inertial_moment);
  FUB_GT_CONTROL_PRINT(mu);
  FUB_GT_CONTROL_PRINT(target_pressure_compressor);
  FUB_GT_CONTROL_PRINT(Q);
  BOOST_LOG(log) << " - initial_turbine_temperature = "
                 << initial_turbine_state.temperature;
  BOOST_LOG(log) << " - initial_turbine_pressure = "
                 << initial_turbine_state.pressure;
}
#undef FUB_GT_CONTROL_PRINT

Control::Control(const PerfectGasConstants& eq, const ControlOptions& options)
    : equation_(eq), options_{options}, state_{
                                            std::make_shared<ControlState>()} {
  state_->current_rpm = options_.rpmmin;
  state_->fuel_consumption = 0.0;
  state_->power_out = 0.0;
  state_->turbine = options_.initial_turbine_state;
  UpdateCompressorFromRPM(*state_, equation_, options_);
}

void Control::UpdatePlena(double mdot_turbine,
                          const PlenumState& turbine_boundary_state,
                          double flux_rho, double flux_spec,
                          double flux_rho_last, Duration dt) {
  UpdateTurbinePlenum(*state_, equation_, options_, mdot_turbine,
                      turbine_boundary_state, flux_rho_last, dt);
  UpdateCompressorPlenum(*state_, equation_, options_, flux_rho, flux_spec, dt);
  ChangeRPM(*state_, options_, dt);
}

// struct ControlState {
//   double current_rpm{};
//   PlenumState compressor{};
//   PlenumState turbine{};
//   double power_out{};
//   double fuel_consumption{};
// };

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
      std::pair<std::string, projection>{
          "compressor_SEC_Mode"s,
          [](const ControlState& s) {
            return static_cast<double>(s.compressor.SEC_Mode);
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
