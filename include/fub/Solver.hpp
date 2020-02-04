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

#include "fub/CartesianCoordinates.hpp"
#include "fub/CompleteFromCons.hpp"
#include "fub/Direction.hpp"
#include "fub/Duration.hpp"
#include "fub/Equation.hpp"
#include "fub/Execution.hpp"
#include "fub/ForEach.hpp"
#include "fub/PatchDataView.hpp"

#include "fub/HyperbolicPatchIntegrator.hpp"
#include "fub/solver/DimensionalSplitLevelIntegrator.hpp"
#include "fub/solver/NoSubcycleSolver.hpp"
#include "fub/solver/SplitSystemSourceLevelIntegrator.hpp"
#include "fub/solver/SubcycleFineFirstSolver.hpp"

#include "fub/NewtonIteration.hpp"

#include "fub/RunSimulation.hpp"

#include "fub/equations/Advection.hpp"
#include "fub/equations/Burgers.hpp"
#include "fub/equations/IdealGasMix.hpp"
#include "fub/equations/PerfectGas.hpp"
#include "fub/equations/ShallowWater.hpp"
#include "fub/equations/CompressibleAdvection.hpp"

#include "fub/equations/ideal_gas_mix/mechanism/Burke2012.hpp"

#include "fub/flux_method/GodunovMethod.hpp"
#include "fub/flux_method/HllMethod.hpp"
#include "fub/flux_method/MusclHancockMethod.hpp"

#include "fub/ext/Log.hpp"
#include "fub/ext/Mpi.hpp"
#include "fub/ext/ProgramOptions.hpp"
#include "fub/ext/Vc.hpp"

#include "fub/cutcell_method/KbnStabilisation.hpp"

#include "fub/split_method/GodunovSplitting.hpp"
#include "fub/split_method/StrangSplitting.hpp"
#include "fub/split_method/StrangSplittingLumped.hpp"

#include "fub/tagging/GradientDetector.hpp"

#include "fub/geometry/Cone.hpp"
#include "fub/geometry/ExpandTube.hpp"
#include "fub/geometry/Geometry.hpp"
#include "fub/geometry/Halfspace.hpp"
#include "fub/geometry/Invert.hpp"
#include "fub/geometry/Polygon.hpp"
#include "fub/geometry/RotateAxis.hpp"

#include "fub/output/AsOutput.hpp"
#include "fub/output/BasicOutput.hpp"
#include "fub/output/CounterOutput.hpp"
#include "fub/output/MultipleOutputs.hpp"
#include "fub/output/OutputAtFrequencyOrInterval.hpp"
#include "fub/output/OutputFactory.hpp"