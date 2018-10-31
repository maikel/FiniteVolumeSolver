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

#ifndef FUB_ODE_SOLVER_RADAU_HPP
#define FUB_ODE_SOLVER_RADAU_HPP

#include "fub/core/span.hpp"

#include <chrono>
#include <stdexcept>
#include <string>
#include <vector>

namespace fub {
namespace radau_error {
struct Exception : std::runtime_error {
  using runtime_error::runtime_error;
};

struct InputIsNotConsistent : Exception {
  InputIsNotConsistent()
      : Exception("ode_solver::Radau: Input is not consistent.") {}
};

struct LargerNMaxIsNeeded : Exception {
  LargerNMaxIsNeeded()
      : Exception("ode_solver::Radau: Larger NMAX is needed.") {}
};

struct StepSizeBecomesTooSmall : Exception {
  StepSizeBecomesTooSmall()
      : Exception("ode_solver::Radau: Step size becomes too small.") {}
};

struct MatrixIsRepeatedlySingular : Exception {
  MatrixIsRepeatedlySingular()
      : Exception("ode_solver::Radau: Matrix is repeatedly singular.") {}
};
} // namespace radau_error

namespace radau_detail {
#if __cpp_noexcept_function_type >= 201510
using system_type = void (*)(int*, double*, double*, double*, float*,
                             int*) noexcept;

using feedback_type = void (*)(int*, double*, double*, double*, int*, int*,
                               int*, float*, int*, int*) noexcept;
#else
using system_type = void (*)(int*, double*, double*, double*, float*, int*);
using feedback_type = void (*)(int*, double*, double*, double*, int*, int*,
                               int*, float*, int*, int*);
#endif

/// @brief This function calls the fortran implementation of E. Hairer.
///
/// We fetched the fortran code from http://[...].
/// It is a fully implicit scheme which is apdatively of 5th, 9th and 13th
/// order. Phillip has shown in his thesis that this scheme is one of few which
/// shows great stability with our kinetic ODEs.
int call_radau(system_type system, // RHS function
               int* par,           // See RPAR
               double t_0,         // Initial time
               double t_end,       // Final integration time
               span<double> y,     // Initial y array
               double h,           // Initial step size
               double rotl,        // Tolerances
               double atol,        //
               span<double> work,  // Temporary array of size LWORK
               span<int> iwork     // Temporary array of size LIWORK
);

/// @brief This function calls the fortran implementation of E. Hairer.
///
/// We fetched the fortran code from http://[...].
/// It is a fully implicit scheme which is apdatively of 5th, 9th and 13th
/// order. Phillip has shown in his thesis that this scheme is one of few which
/// shows great stability with our kinetic ODEs.
int call_radau_feedback(system_type system,    // RHS function
                        int* par,              // See RPAR
                        double t_0,            // Initial time
                        double t_end,          // Final integration time
                        span<double> y,        // Initial y array
                        double h,              // Initial step size
                        double rotl,           // Tolerances
                        double atol,           //
                        span<double> work,     // Temporary array of size LWORK
                        span<int> iwork,       // Temporary array of size LIWORK
                        feedback_type feedback // Feedback function
);
} // namespace radau_detail

struct Radau {
  template <typename System, index Size, typename Feedback,
            typename Allocator = std::allocator<double>>
  static void integrate_feedback(System system, span<double, Size> y0,
                                 double t0, double dt, Feedback feedback,
                                 Allocator alloc = Allocator()) {
    struct IntegrateData_ {
      System* system;
      Feedback* feedback;
    };

    IntegrateData_ user_data{&system, &feedback};

    auto radau_rhs = [](int* n, double* x_, double* y_, double* dydt_, float*,
                        int* ipar) noexcept {
      System& rhs = *(reinterpret_cast<IntegrateData_*>(ipar)->system);
      FUB_ASSERT(Size == dynamic_extent || Size == *n);
      span<double> dydt(dydt_, *n);
      span<const double> y(y_, *n);
      double t = *x_;
      rhs(dydt, y, t);
    };

    auto radau_feedback =
        [](int* nr, double* xold_, double* x_, double* y_, int* cont, int* lrc,
           int* n, float* rpar, int* ipar, int* irtrn) noexcept {
      Feedback& feedback = *(reinterpret_cast<IntegrateData_*>(ipar)->feedback);
      FUB_ASSERT(Size == dynamic_extent || Size == *n);
      span<const double> y(y_, *n);
      double t = *x_;
      feedback(y, t);
    };

    // Set-up properties for the fortran code.

    constexpr double h0 =
        1e-6; // Initial step size. Its okay to set this large. If set
              // low, change WORK[2] to a smaller value!
    double t = t0 + dt;                   // final fime
    constexpr double reltol = 1e-9;  // relative tolerances
    constexpr double abstol = 1e-12; // absolute tolerances

    static constexpr int work_len = Size * (Size + 7 * Size + 3 * 7 + 3) + 20;
    std::vector<double, Allocator> work(work_len, alloc);

    using IAllocator =
        typename std::allocator_traits<Allocator>::template rebind_alloc<int>;
    static constexpr int iwork_len = (2 + (7 - 1) / 2) * Size + 20;
    std::vector<int, IAllocator> iwork(iwork_len, IAllocator(alloc));

    iwork[0] = 1; // Use Hessenberg matrix form
    // iwork[1] = 10000; // Increase maximum number of steps
    iwork[12] = 7;  // Start off with order 13
    work[2] = 0.01; // Recompute Jacobian less often (default 0.001)

    int idid = radau_detail::call_radau_feedback(
        radau_rhs,                          // RHS function
        reinterpret_cast<int*>(&user_data), // IPAR
        t0,                                 // Initial time
        t,                                  // Final integration time
        y0,                                 // Initial y array and output
        h0,                                 // Initial step size
        reltol,                             // Tolerances
        abstol,                             //
        work,                               // Temporary array of size LWORK
        iwork,                              // Temporary array of size LIWORK
        radau_feedback);

    //  IDID= 1  COMPUTATION SUCCESSFUL,
    //  IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
    //  IDID=-1  INPUT IS NOT CONSISTENT,
    //  IDID=-2  LARGER NMAX IS NEEDED,
    //  IDID=-3  STEP SIZE BECOMES TOO SMALL,
    //  IDID=-4  MATRIX IS REPEATEDLY SINGULAR.

    switch (idid) {
    case 1:
      return;
    case 2:
      return;
    case -1:
      throw radau_error::InputIsNotConsistent{};
    case -2:
      throw radau_error::LargerNMaxIsNeeded{};
    case -3:
      idid = radau_detail::call_radau_feedback(
          radau_rhs,                          // RHS function
          reinterpret_cast<int*>(&user_data), // IPAR
          t0,                                 // Initial time
          t / 2,                              // Final integration time
          y0,                                 // Initial y array and output
          h0,                                 // Initial step size
          reltol,                             // Tolerances
          abstol,                             //
          work,                               // Temporary array of size LWORK
          iwork,                              // Temporary array of size LIWORK
          radau_feedback);
      if (idid > 0) {
        idid = radau_detail::call_radau_feedback(
            radau_rhs,                          // RHS function
            reinterpret_cast<int*>(&user_data), // IPAR
            t0,                                 // Initial time
            t / 2,                              // Final integration time
            y0,                                 // Initial y array and output
            h0,                                 // Initial step size
            reltol,                             // Tolerances
            abstol,                             //
            work,                               // Temporary array of size LWORK
            iwork, // Temporary array of size LIWORK
            radau_feedback);
        if (idid > 0) {
          break;
        }
      }
      throw radau_error::StepSizeBecomesTooSmall{};
    case -4:
      throw radau_error::MatrixIsRepeatedlySingular{};
    default:
      throw std::logic_error{"Radau5: The fortran kernel returned an unknown "
                             "return code (which is: " +
                             std::to_string(idid) + ")."};
    }
  }

  template <typename System, index Size,
            typename Allocator = std::allocator<double>>
  static void integrate(System system, span<double, Size> y0, double t0,
                        double dt, Allocator alloc = Allocator()) {
    auto radau_rhs = [](int* n, double* x_, double* y_, double* dydt_, float*,
                        int* ipar) noexcept {
      System& rhs = *reinterpret_cast<System*>(ipar);
      FUB_ASSERT(Size == dynamic_extent || Size == *n);
      span<double> dydt(dydt_, *n);
      span<const double> y(y_, *n);
      double t = *x_;
      rhs(dydt, y, t);
    };

    // Set-up properties for the fortran code.

    constexpr double h0 =
        1e-6; // Initial step size. Its okay to set this large. If set
              // low, change WORK[2] to a smaller value!
    double t = t0 + dt;                   // final fime
    constexpr double reltol = 1e-9;  // relative tolerances
    constexpr double abstol = 1e-12; // absolute tolerances

    static constexpr int work_len = Size * (Size + 7 * Size + 3 * 7 + 3) + 20;
    std::vector<double, Allocator> work(work_len, alloc);

    using IAllocator =
        typename std::allocator_traits<Allocator>::template rebind_alloc<int>;
    static constexpr int iwork_len = (2 + (7 - 1) / 2) * Size + 20;
    std::vector<int, IAllocator> iwork(iwork_len, IAllocator(alloc));

    iwork[0] = 1; // Use Hessenberg matrix form
    // iwork[1] = 10000; // Increase maximum number of steps
    iwork[12] = 7;  // Start off with order 13
    work[2] = 0.01; // Recompute Jacobian less often (default 0.001)

    int idid = radau_detail::call_radau(
        radau_rhs,                                      // RHS function
        reinterpret_cast<int*>(std::addressof(system)), // IPAR
        t0,                                             // Initial time
        t,      // Final integration time
        y0,     // Initial y array and output
        h0,     // Initial step size
        reltol, // Tolerances
        abstol, //
        work,   // Temporary array of size LWORK
        iwork   // Temporary array of size LIWORK
    );

    //  IDID= 1  COMPUTATION SUCCESSFUL,
    //  IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
    //  IDID=-1  INPUT IS NOT CONSISTENT,
    //  IDID=-2  LARGER NMAX IS NEEDED,
    //  IDID=-3  STEP SIZE BECOMES TOO SMALL,
    //  IDID=-4  MATRIX IS REPEATEDLY SINGULAR.

    switch (idid) {
    case 1:
      return;
    case 2:
      return;
    case -1:
      throw radau_error::InputIsNotConsistent{};
    case -2:
      throw radau_error::LargerNMaxIsNeeded{};
    case -3:
      idid = radau_detail::call_radau(
          radau_rhs,                                      // RHS function
          reinterpret_cast<int*>(std::addressof(system)), // IPAR
          t0,                                             // Initial time
          t / 2,  // Final integration time
          y0,     // Initial y array and output
          h0,     // Initial step size
          reltol, // Tolerances
          abstol, //
          work,   // Temporary array of size LWORK
          iwork   // Temporary array of size LIWORK
      );
      if (idid > 0) {
        idid = radau_detail::call_radau(
            radau_rhs,                                      // RHS function
            reinterpret_cast<int*>(std::addressof(system)), // IPAR
            t0,                                             // Initial time
            t / 2,  // Final integration time
            y0,     // Initial y array and output
            h0,     // Initial step size
            reltol, // Tolerances
            abstol, //
            work,   // Temporary array of size LWORK
            iwork   // Temporary array of size LIWORK
        );
        if (idid > 0) {
          break;
        }
      }
      throw radau_error::StepSizeBecomesTooSmall{};
    case -4:
      throw radau_error::MatrixIsRepeatedlySingular{};
    default:
      throw std::logic_error{"Radau5: The fortran kernel returned an unknown "
                             "return code (which is: " +
                             std::to_string(idid) + ")."};
    }
  }
};

} // namespace fub
#endif // !RADAU_HPP
