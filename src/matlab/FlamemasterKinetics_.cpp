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

#include "fub/core/mdspan.hpp"
#include "fub/ideal_gas/FlameMasterReactor.hpp"
#include "fub/ideal_gas/mechanism/Burke2012.hpp"
#include "fub/ideal_gas/mechanism/Gri30.hpp"
#include "fub/ideal_gas/mechanism/Zhao2008Dme.hpp"
#include "fub/ode_solver/CVodeSolver.hpp"
#include "fub/ode_solver/RadauSolver.hpp"

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmacro-redefined"
#elif __CLANG__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmacro-redefined"
#endif
#include "mex.hpp"
#include "mexAdapter.hpp"
#ifdef __GNUC__
#pragma GCC diagnostic pop
#elif __CLANG__
#pragma clang diagnostic pop
#endif

#include <boost/algorithm/string.hpp>

#define TBB_PREVIEW_WAITING_FOR_WORKERS 1
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include <map>
#include <memory>
#include <numeric>
#include <string>

class MexFunction : public matlab::mex::Function {
public:
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    DispatchMemberFunction_(outputs, inputs);
  }

  ~MexFunction() {
    scheduler_.blocking_terminate();
  }

private:
  enum class MemberFunction : int {
    kSetMechanism,
    kSetPressureIsentropic,
    kSetTPX,
    kRecPUTY,
    kAdvance,
    kGetMolarMasses,
    kGetSpeciesNames,
    kFindSpeciesIndex,
    kGetEnthalpies,
    kSetOdeSolver,
    kUnknown
  };

  using FactoryMap =
      std::map<std::string,
               std::unique_ptr<fub::ideal_gas::FlameMasterMechanism>>;

  matlab::data::optional<MemberFunction>
  GetMemberFunction_(matlab::mex::ArgumentList& inputs) {
    if (inputs.size() == 0) {
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>({factory_.createScalar(
                         "FlameMasterKinetics_: Not enough arugments to call "
                         "this function.")}));
      return {};
    }
    if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE) {
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>({factory_.createScalar(
                         "FlameMasterKinetics_: First argument has to be an "
                         "integral value.")}));
      return {};
    }
    matlab::data::TypedArray<double> index = inputs[0];
    int value = index[0];
    return static_cast<MemberFunction>(value);
  }

  void DispatchMemberFunction_(matlab::mex::ArgumentList outputs,
                               matlab::mex::ArgumentList inputs) {
    matlab::data::optional<MemberFunction> memfn = GetMemberFunction_(inputs);
    if (!memfn) {
      return;
    }
    switch (*memfn) {
    case MemberFunction::kSetMechanism:
      SetMechanism_(std::move(inputs));
      break;
    case MemberFunction::kSetTPX:
      SetTPX_(std::move(outputs), std::move(inputs));
      break;
    case MemberFunction::kSetPressureIsentropic:
      SetPressureIsentropic_(std::move(outputs), std::move(inputs));
      break;
    case MemberFunction::kRecPUTY:
      RecPUTY_(std::move(outputs), std::move(inputs));
      break;
    case MemberFunction::kAdvance:
      Advance_(std::move(outputs), std::move(inputs));
      break;
    case MemberFunction::kGetSpeciesNames:
      GetSpeciesNames_(std::move(outputs));
      break;
    case MemberFunction::kGetMolarMasses:
      GetMolarMasses_(std::move(outputs));
      break;
    case MemberFunction::kFindSpeciesIndex:
      FindSpeciesIndex_(std::move(outputs), std::move(inputs));
      break;
    case MemberFunction::kGetEnthalpies:
      GetEnthalpies_(std::move(outputs), std::move(inputs));
      break;
    case MemberFunction::kSetOdeSolver:
      SetOdeSolver_(std::move(inputs));
      break;
    default:
      const int value = static_cast<int>(*memfn);
      std::string what = "FlameMasterKinetics_: Unknown command with value '" +
                         std::to_string(value) + "' requested.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      break;
    }
  }

  static FactoryMap MechanismFactory_() {
    FactoryMap factory{};
    factory["Zhao2008Dme"] = std::make_unique<fub::ideal_gas::Zhao2008Dme>();
    factory["Burke2012"] = std::make_unique<fub::ideal_gas::Burke2012>();
    factory["Gri30"] = std::make_unique<fub::ideal_gas::Gri30>();
    return factory;
  }

  static std::unique_ptr<fub::ideal_gas::FlameMasterMechanism>
  GetMechanism_(const std::string& name) {
    static FactoryMap factory = MechanismFactory_();
    return factory.at(name)->Clone();
  }

  void SetMechanism_(matlab::mex::ArgumentList inputs) {
    matlab::data::CharArray array = inputs[1];
    std::string name = array.toAscii();
    mechanism_ = GetMechanism_(name);
    std::size_t p = tbb::task_scheduler_init::default_num_threads();
    reactor_ =
        std::make_unique<fub::ideal_gas::FlameMasterReactor>(*mechanism_);
    tbb_reactor_ = std::make_unique<
        tbb::enumerable_thread_specific<fub::ideal_gas::FlameMasterReactor>>(
        *reactor_);
    fractions_buffer_.resize(reactor_->getNSpecies());
    tbb_fractions_buffer_ =
        tbb::enumerable_thread_specific<std::vector<double>>(fractions_buffer_);
    engine_->feval(u"fprintf", 0,
                   std::vector<matlab::data::Array>({factory_.createScalar(
                       "FlamemasterKinetics_: Using the '" + name +
                       "' mechanism now.\n")}));
  }

  std::ptrdiff_t GetNumberOfVariables_() {
    const int n_vars = static_cast<int>(Variable::VARIABLE_COUNT) - 1;
    const int n_i_vars = static_cast<int>(VariableI::VARIABLE_I_COUNT);
    return n_vars + n_i_vars + reactor_->getNSpecies() - 1;
  }

  void GatherMoles_(fub::span<double> moles,
                    const matlab::data::TypedArray<double>& X,
                    std::ptrdiff_t cell) const {
    const int n_species = moles.size();
    for (int s = 0; s < n_species; ++s) {
      moles[s] = X[s][cell];
    }
  }

  enum class Variable : int {
    density,
    momentum,
    energy,
    species,
    VARIABLE_COUNT
  };

  enum class VariableI : int {
    pressure,
    temperature,
    cp,
    cv,
    VARIABLE_I_COUNT
  };

  template <typename T>
  using MdSpan = fub::basic_mdspan<T, fub::DynamicExtents<2>, fub::layout_left>;

  static double& At_(MdSpan<double> states, Variable var, std::ptrdiff_t cell) {
    return states(static_cast<int>(var), cell);
  }

  static double& At_(MdSpan<double> states, VariableI var,
                     std::ptrdiff_t cell) {
    const int ivar = static_cast<int>(var);
    const int size = states.extent(0);
    FUB_ASSERT(0 <= ivar && ivar < size);
    return states(size - ivar - 1, cell);
  }

  static void
  UpdateStateFromKinetics_(const fub::ideal_gas::FlameMasterReactor& reactor,
                           MdSpan<double> states, std::ptrdiff_t cell,
                           double velocity = 0.0) {
    auto at = [&](auto var) -> double& { return At_(states, var, cell); };
    const double rho = reactor.getDensity();
    at(Variable::density) = rho;
    at(Variable::momentum) = rho * velocity;
    at(Variable::energy) =
        rho * (reactor.getInternalEnergy() + 0.5 * velocity * velocity);
    at(VariableI::cp) = reactor.getCp();
    at(VariableI::cv) = reactor.getCv();
    at(VariableI::temperature) = reactor.getTemperature();
    at(VariableI::pressure) = reactor.getPressure();
    const int rhoY1 = static_cast<int>(Variable::species);
    fub::span<const double> Y = reactor.getMassFractions();
    for (int s = 0; s < reactor.getNSpecies() - 1; ++s) {
      at(Variable(rhoY1 + s)) = rho * Y[s];
    }
  }

  void UpdateKineticsFromState_(MdSpan<double> states,
                                std::ptrdiff_t cell) const {
    auto at = [&](auto var) -> double& { return At_(states, var, cell); };
    const int rhoY = static_cast<int>(Variable::species);
    for (int s = 0; s < reactor_->getNSpecies() - 1; ++s) {
      fractions_buffer_[s] = at(Variable(rhoY + s));
    }
    const double rho = at(Variable::density);
    fractions_buffer_[reactor_->getNSpecies() - 1] =
        std::max(0.0, rho - std::accumulate(fractions_buffer_.begin(),
                                            fractions_buffer_.end() - 1, 0.0));
    reactor_->setMassFractions(fractions_buffer_);
    reactor_->setTemperature(at(VariableI::temperature));
    reactor_->setPressure(at(VariableI::pressure));
  }

  void SetPressureIsentropic_(matlab::mex::ArgumentList outputs,
                              matlab::mex::ArgumentList inputs) {
    matlab::data::TypedArray<double> array = std::move(inputs[1]);
    matlab::data::ArrayDimensions dims = array.getDimensions();
    if (dims.size() != 2) {
      throw std::runtime_error("FlamemasterKinetics_::SetPressureIsentropic: "
                               "input dimension is wrong.");
    }
    matlab::data::buffer_ptr_t<double> buffer = array.release();
    MdSpan<double> states(buffer.get(), dims[0], dims[1]);
    matlab::data::TypedArray<double> P = std::move(inputs[2]);
    for (std::size_t cell = 0; cell < dims[1]; ++cell) {
      UpdateKineticsFromState_(states, cell);
      reactor_->setPressureIsentropic(P[cell]);
      UpdateStateFromKinetics_(*reactor_, states, cell);
    }
    outputs[0] = factory_.createArrayFromBuffer(dims, std::move(buffer));
  }

  void GetSpeciesNames_(matlab::mex::ArgumentList outputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    const std::size_t n_species = reactor_->getNSpecies();
    matlab::data::CellArray names =
        factory_.createArray<matlab::data::Array>({n_species, 1});
    for (std::size_t s = 0; s < n_species; ++s) {
      names[s] = factory_.createCharArray(reactor_->getSpeciesName(s));
    }
    outputs[0] = std::move(names);
  }

  void Advance_(matlab::mex::ArgumentList outputs,
                matlab::mex::ArgumentList inputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    matlab::data::TypedArray<double> array = std::move(inputs[1]);
    matlab::data::ArrayDimensions dims = array.getDimensions();
    if (dims.size() != 2) {
      throw std::runtime_error("FlamemasterKinetics_::SetPressureIsentropic: "
                               "input dimension is wrong.");
    }
    matlab::data::buffer_ptr_t<double> buffer = array.release();
    const std::size_t n_cells = dims[1];
    const int n_species = reactor_->getNSpecies();
    std::vector<double> rates_buffer(n_species * n_cells);
    MdSpan<double> states(buffer.get(), dims[0], dims[1]);
    MdSpan<double> rates(rates_buffer.data(), n_species, n_cells);
    matlab::data::TypedArray<double> times = std::move(inputs[2]);
    const double dt = times[0];
    tbb::parallel_for(std::size_t{0}, n_cells, [&](std::size_t cell) {
      std::vector<double>& fractions_buffer = tbb_fractions_buffer_.local();
      fub::ideal_gas::FlameMasterReactor reactor = tbb_reactor_->local();
      auto at = [&](auto var, std::size_t cell) -> double& {
        return At_(states, var, cell);
      };
      const double rhoE = at(Variable::energy, cell);
      const double rhou = at(Variable::momentum, cell);
      const double rho = at(Variable::density, cell);
      const int rhoY = static_cast<int>(Variable::species);
      for (int s = 0; s < reactor.getNSpecies() - 1; ++s) {
        fractions_buffer[s] = at(Variable(rhoY + s), cell);
      }
      fractions_buffer[reactor.getNSpecies() - 1] =
          std::max(0.0, rho - std::accumulate(fractions_buffer.begin(),
                                              fractions_buffer.end() - 1, 0.0));
      reactor.setMassFractions(fractions_buffer);
      reactor.setDensity(at(Variable::density, cell));
      const double e = (rhoE - 0.5 * rhou * rhou / rho) / rho;
      reactor.setInternalEnergy(e);
      reactor.advance(dt);
      UpdateStateFromKinetics_(reactor, states, cell, rhou / rho);
      fub::span<const double> current_rates = reactor.getProductionRates();
      for (int s = 0; s < n_species; ++s) {
        rates(s, cell) = current_rates[s];
      }
    });
    outputs[0] = factory_.createArrayFromBuffer(dims, std::move(buffer));
    outputs[1] =
        factory_.createArray({static_cast<std::size_t>(n_species), n_cells},
                             rates_buffer.begin(), rates_buffer.end());
  }

  void SetTPX_(matlab::mex::ArgumentList outputs,
               matlab::mex::ArgumentList inputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    matlab::data::TypedArray<double> T = std::move(inputs[1]);
    matlab::data::TypedArray<double> P = std::move(inputs[2]);

    const std::size_t num_variables = GetNumberOfVariables_();
    const std::size_t num_cells = T.getNumberOfElements();
    std::vector<double> buffer(num_variables * num_cells);

    MdSpan<double> states(buffer.data(), num_variables, num_cells);

    if (inputs[3].getType() == matlab::data::ArrayType::CHAR) {
      matlab::data::CharArray X = std::move(inputs[3]);
      std::string Xstr = X.toAscii();
      reactor_->setMoleFractions(Xstr);
      for (std::size_t cell = 0; cell < num_cells; ++cell) {
        reactor_->setTemperature(T[cell]);
        reactor_->setPressure(P[cell]);
        reactor_->advance(0);
        UpdateStateFromKinetics_(*reactor_, states, cell);
      }
      outputs[0] = factory_.createArray<double>({num_variables, num_cells},
                                                states.data(),
                                                states.data() + states.size());
    } else {
      matlab::data::TypedArray<double> X = std::move(inputs[3]);
      for (std::size_t cell = 0; cell < num_cells; ++cell) {
        GatherMoles_(fractions_buffer_, X, cell);
        reactor_->setMoleFractions(fractions_buffer_);
        reactor_->setTemperature(T[cell]);
        reactor_->setPressure(P[cell]);
        reactor_->advance(0);
        UpdateStateFromKinetics_(*reactor_, states, cell);
      }
      outputs[0] = factory_.createArray<double>({num_variables, num_cells},
                                                states.data(),
                                                states.data() + states.size());
    }
  }

  void RecPUTY_(matlab::mex::ArgumentList outputs,
                matlab::mex::ArgumentList inputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    matlab::data::TypedArray<double> array = std::move(inputs[1]);
    matlab::data::ArrayDimensions dims = array.getDimensions();
    if (dims.size() != 2) {
      throw std::runtime_error("FlamemasterKinetics_::SetPressureIsentropic: "
                               "input dimension is wrong.");
    }
    matlab::data::buffer_ptr_t<double> buffer = array.release();
    const std::size_t num_variables = dims[0];
    const std::size_t num_cells = dims[1];
    MdSpan<double> states(buffer.get(), num_variables, num_cells);
    std::vector<double> recbuf(GetNumberOfVariables_() * num_cells);
    MdSpan<double> rec(recbuf.data(), GetNumberOfVariables_(), num_cells);
    for (std::size_t cell = 0; cell < num_cells; ++cell) {
      for (int s = 0; s < reactor_->getNSpecies(); ++s) {
        fractions_buffer_[s] = states(3 + s, cell);
      }
      reactor_->setMassFractions(fractions_buffer_);
      reactor_->setTemperature(states(2, cell));
      reactor_->setPressure(states(0, cell));
      reactor_->advance(0);
      const double velocity = states(1, cell);
      UpdateStateFromKinetics_(*reactor_, rec, cell, velocity);
    }
    std::size_t nvars = GetNumberOfVariables_();
    outputs[0] =
        factory_.createArray({nvars, num_cells}, recbuf.begin(), recbuf.end());
  }

  void GetMolarMasses_(matlab::mex::ArgumentList outputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    fub::span<const double> masses = reactor_->getMolarMasses();
    outputs[0] = factory_.createArray<double>(
        {static_cast<std::size_t>(masses.size()), 1}, masses.begin(),
        masses.end());
  }

  void FindSpeciesIndex_(matlab::mex::ArgumentList outputs,
                         matlab::mex::ArgumentList inputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    matlab::data::ArrayDimensions dims = inputs[0].getDimensions();
    fub::span<const std::string> names = reactor_->getSpeciesNames();
    matlab::data::CharArray input = inputs[1];
    std::string needle = input.toAscii();
    boost::trim_right(needle);
    auto pos = std::find(names.begin(), names.end(), needle);
    if (pos == names.end()) {
      outputs[0] = factory_.createScalar<double>(-1);
    } else {
      outputs[0] = factory_.createScalar<double>(pos - names.begin());
    }
  }

  void GetEnthalpies_(matlab::mex::ArgumentList outputs,
                      matlab::mex::ArgumentList inputs) {
    if (!reactor_) {
      std::string what = "FlameMasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    matlab::data::TypedArray<double> Ts = std::move(inputs[1]);
    const double T = Ts[0];
    reactor_->setTemperature(T);
    fub::span<const double> h = reactor_->getEnthalpies();
    std::size_t n_species = h.size();
    outputs[0] =
        factory_.createArray<double>({n_species, 1}, h.begin(), h.end());
  }

  void SetOdeSolver_(matlab::mex::ArgumentList inputs) {
    if (!reactor_) {
      std::string what = "FlamemasterKinetics_: No mechanism loaded yet.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
      return;
    }
    matlab::data::CharArray name = inputs[1];
    std::string ode_solver_name = name.toAscii();
    if (ode_solver_name == "CVode") {
      reactor_->setOdeSolver(
          std::make_unique<fub::CVodeSolver>(reactor_->getNSpecies() + 1));
      tbb_reactor_ = std::make_unique<
          tbb::enumerable_thread_specific<fub::ideal_gas::FlameMasterReactor>>(
          *reactor_);
      engine_->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>({factory_.createScalar(
                         "FlamemasterKinetics_: Using the '" + ode_solver_name +
                         "' ode solver now.\n")}));
    } else if (ode_solver_name == "RADAU") {
      reactor_->setOdeSolver(std::make_unique<fub::RadauSolver>());
      tbb_reactor_ = std::make_unique<
          tbb::enumerable_thread_specific<fub::ideal_gas::FlameMasterReactor>>(
          *reactor_);
      engine_->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>({factory_.createScalar(
                         "FlamemasterKinetics_: Using the '" + ode_solver_name +
                         "' ode solver now.\n")}));
    } else {
      std::string what =
          "FlamemasterKinetics_: Unknown ode solver '" + ode_solver_name + "'.";
      engine_->feval(u"error", 0,
                     std::vector<matlab::data::Array>(
                         {factory_.createScalar(std::move(what))}));
    }
  }

  std::shared_ptr<matlab::engine::MATLABEngine> engine_ = getEngine();
  matlab::data::ArrayFactory factory_;
  mutable std::vector<double> fractions_buffer_;
  mutable tbb::enumerable_thread_specific<std::vector<double>>
      tbb_fractions_buffer_;
  tbb::task_scheduler_init scheduler_{tbb::task_scheduler_init::automatic};
  std::unique_ptr<fub::ideal_gas::FlameMasterMechanism> mechanism_{};
  std::unique_ptr<fub::ideal_gas::FlameMasterReactor> reactor_{};
  std::unique_ptr<
      tbb::enumerable_thread_specific<fub::ideal_gas::FlameMasterReactor>>
      tbb_reactor_{};
};