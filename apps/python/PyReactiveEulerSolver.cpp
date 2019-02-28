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

#include "fub/grid/SAMRAI/ideal_gas/KineticDriver.hpp"

#include "fub/ideal_gas/mechanism/Burke2012.hpp"

#ifdef HAVE_SYS_TIMES_H
#undef HAVE_SYS_TIMES_H
#endif

#ifdef HAVE_UNISTD_H
#undef HAVE_UNISTD_H
#endif

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

using Burke2012 = fub::ideal_gas::Burke2012;
using FlamemasterReactor = fub::ideal_gas::FlameMasterReactor;
using Equation = fub::ideal_gas::IdealGasKinetics;
using Solver = fub::ideal_gas::KineticDriver;

using PyArray = pybind11::array_t<double, pybind11::array::c_style |
                                              pybind11::array::forcecast>;

PYBIND11_MODULE(PyReactiveEulerSolver, module) {
  module.doc() =
      "A module providing python types to solver the reactive euler equations.";

  pybind11::class_<Burke2012>(module, "Burke2012").def(pybind11::init<>());

  pybind11::class_<FlamemasterReactor>(module, "FlamemasterReactor")
      .def(pybind11::init<>(
          [](const Burke2012& burke) { return FlamemasterReactor(burke); }))
      .def("GetNSpecies", &FlamemasterReactor::getNSpecies)
      .def("GetPressure", &FlamemasterReactor::getPressure)
      .def("SetPressure", &FlamemasterReactor::setPressure)
      .def("GetTemperature", &FlamemasterReactor::getTemperature)
      .def("SetTemperature", &FlamemasterReactor::setTemperature)
      .def("SetDensity", &FlamemasterReactor::setDensity)
      .def("GetDensity", &FlamemasterReactor::getDensity)
      .def("Advance",
           [](FlamemasterReactor& reactor, double dt) { reactor.advance(dt); })
      .def("SetMassFractions", (void (FlamemasterReactor::*)(std::string))(
                                   &FlamemasterReactor::setMassFractions))
      .def("SetMassFractions",
           [](FlamemasterReactor& reactor, PyArray fractions) {
             pybind11::buffer_info buf = fractions.request();
             if (buf.ndim != 1) {
               throw std::runtime_error("Number of dimensions must be one");
             }
             fub::span<const double> fracs(static_cast<double*>(buf.ptr),
                                           buf.size);
             reactor.setMassFractions(fracs);
           })
      .def("GetMassFractions", [](const FlamemasterReactor& reactor) {
        fub::span<const double> fractions = reactor.getMassFractions();
        PyArray result(fractions.size());
        
        return result;
      });
}