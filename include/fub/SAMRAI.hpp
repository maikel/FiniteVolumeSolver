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

#ifndef FUB_SAMRAI_HPP
#define FUB_SAMRAI_HPP

#include "fub/SAMRAI/ScopeGuard.hpp"
#include "fub/SAMRAI/RegisterVariables.hpp"
#include "fub/SAMRAI/PatchHierarchy.hpp"
#include "fub/SAMRAI/GriddingAlgorithm.hpp"
#include "fub/SAMRAI/IntegratorContext.hpp"
#include "fub/SAMRAI/BoundaryCondition.hpp"
#include "fub/SAMRAI/ViewPatch.hpp"

#include "fub/SAMRAI/Reconstruction.hpp"
#include "fub/SAMRAI/TimeIntegrator.hpp"
#include "fub/SAMRAI/FluxMethod.hpp"

#include "fub/SAMRAI/tagging/GradientDetector.hpp"
#include "fub/SAMRAI/tagging/ConstantBox.hpp"

#endif