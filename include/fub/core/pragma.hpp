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

#ifndef FUB_CORE_PRAGMA_HPP
#define FUB_CORE_PRAGMA_HPP

#define DIAG_STR(s) #s
#define DIAG_JOINSTR(x, y) DIAG_STR(x##y)
#ifdef _MSC_VER
#define DIAG_DO_PRAGMA(x) __pragma(#x)
#define DIAG_PRAGMA(compiler, x) DIAG_DO_PRAGMA(warning(x))
#else
#define DIAG_DO_PRAGMA(x) _Pragma(#x)
#define DIAG_PRAGMA(compiler, x) DIAG_DO_PRAGMA(compiler diagnostic x)
#endif
#if defined(__clang__)
#define DISABLE_WARNING(gcc_unused, clang_option, msvc_unused)                 \
  DIAG_PRAGMA(clang, push)                                                     \
  DIAG_PRAGMA(clang, ignored DIAG_JOINSTR(-W, clang_option))
#define ENABLE_WARNING(gcc_unused, clang_option, msvc_unused)                  \
  DIAG_PRAGMA(clang, pop)
#elif defined(_MSC_VER)
#define DISABLE_WARNING(gcc_unused, clang_unused, msvc_errorcode)              \
  DIAG_PRAGMA(msvc, push) DIAG_DO_PRAGMA(warning(disable :##msvc_errorcode))
#define ENABLE_WARNING(gcc_unused, clang_unused, msvc_errorcode)               \
  DIAG_PRAGMA(msvc, pop)
#elif defined(__GNUC__)
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#define DISABLE_WARNING(gcc_option, clang_unused, msvc_unused)                 \
  DIAG_PRAGMA(GCC, push) DIAG_PRAGMA(GCC, ignored DIAG_JOINSTR(-W, gcc_option))
#define ENABLE_WARNING(gcc_option, clang_unused, msvc_unused)                  \
  DIAG_PRAGMA(GCC, pop)
#else
#define DISABLE_WARNING(gcc_option, clang_unused, msvc_unused)                 \
  DIAG_PRAGMA(GCC, ignored DIAG_JOINSTR(-W, gcc_option))
#define ENABLE_WARNING(gcc_option, clang_option, msvc_unused)                  \
  DIAG_PRAGMA(GCC, warning DIAG_JOINSTR(-W, gcc_option))
#endif
#endif

#endif