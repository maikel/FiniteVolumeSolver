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

#ifndef FUB_AMREX_HPP
#define FUB_AMREX_HPP

#include <AMReX.H>
#include <AMReX_IntVect.H>

#include "fub/AMReX/GriddingAlgorithm.hpp"
#include "fub/AMReX/IntegratorContext.hpp"

#include "fub/AMReX/ScopeGuard.hpp"

#include "fub/equations/CompressibleAdvection.hpp"

#include "fub/AMReX/tagging_method/ConstantRegion.hpp"
#include "fub/AMReX/tagging_method/GradientDetector.hpp"
#include "fub/AMReX/tagging_method/TagAllOf.hpp"
#include "fub/AMReX/tagging_method/TagBuffer.hpp"

#include "fub/AMReX/boundary_condition/BoundarySet.hpp"
#include "fub/AMReX/boundary_condition/IsentropicPressureBoundary.hpp"
#include "fub/AMReX/boundary_condition/PressureValveBoundary.hpp"
#include "fub/AMReX/boundary_condition/ReflectiveBoundary.hpp"
#include "fub/AMReX/boundary_condition/TransmissiveBoundary.hpp"
#include "fub/AMReX/boundary_condition/GenericPressureValveBoundary.hpp"
#include "fub/AMReX/boundary_condition/IsentropicPressureExpansion.hpp"

#include "fub/AMReX/initial_data/ConstantData.hpp"

#include "fub/AMReX/IgniteDetonation.hpp"

#include "fub/AMReX/ForEachFab.hpp"
#include "fub/AMReX/ForEachIndex.hpp"

#include "fub/AMReX/FluxMethodAdapter.hpp"
#include "fub/AMReX/Reconstruction.hpp"
#include "fub/AMReX/TimeIntegrator.hpp"
#include "fub/AMReX/FluxMethodFactory.hpp"

#include "fub/AMReX/Geometry.hpp"
#include "fub/geometry/Union.hpp"

#include "fub/AMReX/AxialSourceTerm.hpp"
#include "fub/equations/ideal_gas_mix/KineticSourceTerm.hpp"

#include "fub/equations/perfect_gas_mix/IgnitionDelayKinetics.hpp"
#include "fub/equations/perfect_gas_mix/ArrheniusKinetics.hpp"

#include "fub/AMReX/output/DebugOutput.hpp"
#include "fub/AMReX/output/WriteHdf5.hpp"
#include "fub/AMReX/output/WritePlotfiles.hpp"

#endif
