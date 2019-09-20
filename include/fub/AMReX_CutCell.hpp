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

#ifndef FUB_AMREX_CUTCELL_HPP
#define FUB_AMREX_CUTCELL_HPP

#include "fub/AMReX/cutcell/BoundaryCondition.hpp"
#include "fub/AMReX/cutcell/FluxMethod.hpp"
#include "fub/AMReX/cutcell/GriddingAlgorithm.hpp"
#include "fub/AMReX/cutcell/IndexSpace.hpp"
#include "fub/AMReX/cutcell/InitialData.hpp"
#include "fub/AMReX/cutcell/IntegratorContext.hpp"
#include "fub/AMReX/cutcell/PatchHierarchy.hpp"
#include "fub/AMReX/cutcell/Reconstruction.hpp"
#include "fub/AMReX/cutcell/Tagging.hpp"

#include "fub/AMReX/cutcell/TimeIntegrator.hpp"

#include "fub/AMReX/cutcell/initial_data/RiemannProblem.hpp"

#include "fub/AMReX/cutcell/tagging/ConstantRegion.hpp"
#include "fub/AMReX/cutcell/tagging/GradientDetector.hpp"
#include "fub/AMReX/cutcell/tagging/TagBuffer.hpp"
#include "fub/AMReX/cutcell/tagging/TagCutCells.hpp"

#include "fub/AMReX/cutcell/boundary_condition/BoundarySet.hpp"
#include "fub/AMReX/cutcell/boundary_condition/IsentropicPressureBoundary.hpp"
#include "fub/AMReX/cutcell/boundary_condition/MassflowBoundary.hpp"
#include "fub/AMReX/cutcell/boundary_condition/TransmissiveBoundary.hpp"

#include "fub/AMReX/multi_block/MultiBlockBoundary.hpp"
#include "fub/AMReX/multi_block/MultiBlockGriddingAlgorithm.hpp"
#include "fub/AMReX/multi_block/MultiBlockIntegratorContext.hpp"
#include "fub/AMReX/multi_block/MultiBlockKineticSourceTerm.hpp"
#include "fub/AMReX/multi_block/MultiBlockIgniteDetonation.hpp"

#endif
