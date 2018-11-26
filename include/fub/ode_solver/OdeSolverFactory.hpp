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

#ifndef FUB_ODE_SOLVER_ODE_SOLVER_FACTORY_HPP
#define FUB_ODE_SOLVER_ODE_SOLVER_FACTORY_HPP

#include "fub/ode_solver/OdeSolver.hpp"

#include <functional>
#include <map>
#include <string>

namespace fub {

enum class OdeSolverType { RADAU, CVode };

struct OdeSolverOptions {
  OdeSolverType type;
};

class OdeSolverFactory {
public:
  using UniqueSolver = std::unique_ptr<OdeSolver>;
  using SpecificFactory = std::function<UniqueSolver(const OdeSolverOptions&)>;

  OdeSolverFactory(const OdeSolverFactory&) = delete;
  OdeSolverFactory& operator=(const OdeSolverFactory&) = delete;
  OdeSolverFactory(OdeSolverFactory&&) = delete;
  OdeSolverFactory& operator=(OdeSolverFactory&&) = delete;

  ~OdeSolverFactory() = default;

  static OdeSolverFactory& GetInstance();

  void RegisterSpecificFactory(std::string name, SpecificFactory factory);

  UniqueSolver MakeSolver(const std::string& name,
                          const OdeSolverOptions& options = {}) const;

private:
  OdeSolverFactory() = default;
  std::map<std::string, SpecificFactory> factory_{};
};

struct RegisterSpecificFactory {
  RegisterSpecificFactory(std::string name,
                          OdeSolverFactory::SpecificFactory factory);
};

} // namespace fub

#endif