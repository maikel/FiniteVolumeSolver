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

#include "fub/ode_solver/OdeSolverFactory.hpp"

namespace fub {

OdeSolverFactory& OdeSolverFactory::GetInstance() {
  static OdeSolverFactory factory{};
  return factory;
}

void OdeSolverFactory::RegisterSpecificFactory(std::string name,
                                               SpecificFactory factory) {
  factory_.insert({std::move(name), std::move(factory)});
}

OdeSolverFactory::UniqueSolver
OdeSolverFactory::MakeSolver(const std::string& name,
                             const OdeSolverOptions& options) const {
  return factory_.at(name)(options);
}

RegisterSpecificFactory::RegisterSpecificFactory(
    std::string name, OdeSolverFactory::SpecificFactory factory) {
  OdeSolverFactory& instance = OdeSolverFactory::GetInstance();
  instance.RegisterSpecificFactory(std::move(name), std::move(factory));
}

} // namespace fub