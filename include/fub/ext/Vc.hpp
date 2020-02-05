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

#ifndef FUB_EXT_VC_HPP
#define FUB_EXT_VC_HPP

#include <Vc/vector.h>
#include <immintrin.h>
#include <xmmintrin.h>

namespace fub {
inline void EnableFloatingPointExceptions() {
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
}

template <typename Abi>
Vc::Vector<double, Abi> mask_load(const double* p, Vc::Mask<double, Abi> mask) {
  if constexpr (std::is_same_v<Abi, Vc::VectorAbi::Sse>) {
#ifdef __AVX__
    __m128i native_mask = *reinterpret_cast<__m128i*>(&mask);
    __m128d xmm = _mm_maskload_pd(p, native_mask);
    Vc::Vector<double, Abi> x(xmm);
    return x;
#else
    Vc::Vector<double, Abi> x;
    if (mask[0])
      x[0] = p[0];
    if (mask[1])
      x[1] = p[1];
    return x;
#endif
  } else if constexpr (std::is_same_v<Abi, Vc::VectorAbi::Avx>) {
    __m256i native_mask = *reinterpret_cast<__m256i*>(&mask);
    __m256d ymm = _mm256_maskload_pd(p, native_mask);
    Vc::Vector<double, Abi> x(ymm);
    return x;
  }
  //  else if constexpr (std::is_same_v<Abi, Vc::VectorAbi::Avx512) {
  //    __m512i native_mask = *reinterpret_cast<__m512i*>(&mask);
  //    __m512d zmm = _mm512_maskload_pd(p, native_mask);
  //    Vc::Vector<double, Abi> x(zmm);
  //    return x;
  //  }
}

} // namespace fub

#endif
