// Copyright (c) 2020 Maikel Nadolski
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

#ifndef FUB_EQUATIONS_PERFECT_GAS_MIX_HLLEM_HPP
#define FUB_EQUATIONS_PERFECT_GAS_MIX_HLLEM_HPP

#include "fub/equations/PerfectGasMix.hpp"
#include "fub/equations/perfect_gas/HllemMethod.hpp"

namespace fub::perfect_gas {

extern template struct Hllem<PerfectGasMix<1>>;
extern template struct Hllem<PerfectGasMix<2>>;
extern template struct Hllem<PerfectGasMix<3>>;

} // namespace fub::perfect_gas

namespace fub {

extern template class FluxMethod<perfect_gas::Hllem<PerfectGasMix<1>>>;
extern template class FluxMethod<perfect_gas::Hllem<PerfectGasMix<2>>>;
extern template class FluxMethod<perfect_gas::Hllem<PerfectGasMix<3>>>;

} // namespace fub

#endif