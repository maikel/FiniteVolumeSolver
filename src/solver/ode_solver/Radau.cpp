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

#include "src/solver/ode_solver/Radau.hpp"

namespace {
typedef void (*rhs_t)(int *n, double *x, double *y, double *f, float *rpar,
                      int *ipar);
typedef void (*jac_t)(int *n, double *x, double *y, double *dfy, int *ldfy,
                      float *rpar, int *ipar);
typedef void (*mas_t)(int *n, double *am, int *lmas, float *rpar, int *ipar);
typedef void (*sol_t)(int *nr, double *xold, double *x, double *y, int *cont,
                      int *lrc, int *n, float *rpar, int *ipar, int *irtrn);

extern "C" {
void radau_(int *N,       // Dimension
            rhs_t FCN,    // RHS function
            double *X,    // Initial time
            double *Y,    // Initial y array
            double *XEND, // Final integration time
            double *H,    // Initial step size
            double *RTOL, // Tolerances
            double *ATOL, //
            int *ITOL,    // Whether rtol and atol are vector valued
            jac_t JAC,    // Jacobian function
            int *IJAC,    // Whether to use JAC (1) or finite differences
            int *MLJAC,   // Band-width of the jacobian. N for full.
            int *MUJAC,   // Upper band-width (0 is MLJAC == N)
            mas_t MAS,    // Mass matrix function
            int *IMAS,    // Whether to use the mass matrix function
            int *MLMAS,   // Band-width of the mass matrix
            int *MUMAS,   // Upper band-widh of the mass matrix
            sol_t SOLOUT, // Dense output function
            int *IOUT,    // Wether to call the dense output function
            double *WORK, // Temporary array of size LWORK
            int *LWORK,   // N*(LJAC+LMAS+7*N+3*NSMAX+3)+20
            int *IWORK,   // Temporary array of size LIWORK
            int *LIWORK,  // (2+(NSMAX-1)/2)*N+20
            float *RPAR,  // User-supplied RHS arguments
            int *IPAR,    // See RPAR
            int *IDID     // Return value
                          //  IDID= 1  COMPUTATION SUCCESSFUL,
            //  IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
            //  IDID=-1  INPUT IS NOT CONSISTENT,
            //  IDID=-2  LARGER NMAX IS NEEDED,
            //  IDID=-3  STEP SIZE BECOMES TOO SMALL,
            //  IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
);
}
} // namespace

namespace fub {
namespace radau_detail {
int call_radau(system_type system,   // RHS function
               int *ipar,         // IPAR
               double t_0,        // Initial time
               double t_end,      // Final integration time
               span<double> y,    // Initial y array
               double h,          // Initial step size
               double reltol,     // Tolerances
               double abstol,     //
               span<double> work, // Temporary array of size LWORK
               span<int> iwork    // Temporary array of size LIWORK
) {
  int no[] = {0, 0, 0, 0};
  int n = static_cast<int>(y.size());
  int lwork = static_cast<int>(work.size());
  int liwork = static_cast<int>(iwork.size());
  int nc[] = {n, n, n};
  int bw[] = {0, 0};
  int idid = 0;
  radau_(&nc[0],       // Dimension
         system,       // RHS function
         &t_0,         // Initial time
         y.data(),     // Initial y array and output
         &t_end,       // Final integration time
         &h,           // Initial step size
         &reltol,      // Tolerances
         &abstol,      //
         &no[0],       // Whether rtol and atol are vector valued
         nullptr,      // Jacobian function
         &no[1],       // Whether to use JAC (1) or finite differences
         &nc[1],       // Band-width of the jacobian. N for full.
         &bw[0],       // Upper band-width (0 is MLJAC == N)
         nullptr,      // Mass matrix function
         &no[2],       // Whether to use the mass matrix function
         &nc[2],       // Band-width of the mass matrix
         &bw[1],       // Upper band-widh of the mass matrix
         nullptr,      // Dense output function
         &no[3],       // Wether to call the dense output function
         work.data(),  // Temporary array of size LWORK
         &lwork,       // N*(LJAC+LMAS+7*N+3*NSMAX+3)+20
         iwork.data(), // Temporary array of size LIWORK
         &liwork,      // (2+(NSMAX-1)/2)*N+20
         nullptr,      // User-supplied RHS arguments
         ipar,         // See RPAR
         &idid         // Return value
                       // IDID= 1  COMPUTATION SUCCESSFUL,
                       // IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
                       // IDID=-1  INPUT IS NOT CONSISTENT,
                       // IDID=-2  LARGER NMAX IS NEEDED,
                       // IDID=-3  STEP SIZE BECOMES TOO SMALL,
                       // IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
  );
  return idid;
}

int call_radau_feedback(system_type system,   // RHS function
               int *ipar,         // IPAR
               double t_0,        // Initial time
               double t_end,      // Final integration time
               span<double> y,    // Initial y array
               double h,          // Initial step size
               double reltol,     // Tolerances
               double abstol,     //
               span<double> work, // Temporary array of size LWORK
               span<int> iwork,    // Temporary array of size LIWORK
               feedback_type feedback
) {
  int no[] = {0, 0, 0, 0};
  int yes[] = {1};
  int n = static_cast<int>(y.size());
  int lwork = static_cast<int>(work.size());
  int liwork = static_cast<int>(iwork.size());
  int nc[] = {n, n, n};
  int bw[] = {0, 0};
  int idid = 0;
  radau_(&nc[0],       // Dimension
         system,       // RHS function
         &t_0,         // Initial time
         y.data(),     // Initial y array and output
         &t_end,       // Final integration time
         &h,           // Initial step size
         &reltol,      // Tolerances
         &abstol,      //
         &no[0],       // Whether rtol and atol are vector valued
         nullptr,      // Jacobian function
         &no[1],       // Whether to use JAC (1) or finite differences
         &nc[1],       // Band-width of the jacobian. N for full.
         &bw[0],       // Upper band-width (0 is MLJAC == N)
         nullptr,      // Mass matrix function
         &no[2],       // Whether to use the mass matrix function
         &nc[2],       // Band-width of the mass matrix
         &bw[1],       // Upper band-widh of the mass matrix
         feedback,      // Dense output function
         &yes[0],       // Wether to call the dense output function
         work.data(),  // Temporary array of size LWORK
         &lwork,       // N*(LJAC+LMAS+7*N+3*NSMAX+3)+20
         iwork.data(), // Temporary array of size LIWORK
         &liwork,      // (2+(NSMAX-1)/2)*N+20
         nullptr,      // User-supplied RHS arguments
         ipar,         // See RPAR
         &idid         // Return value
                       // IDID= 1  COMPUTATION SUCCESSFUL,
                       // IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
                       // IDID=-1  INPUT IS NOT CONSISTENT,
                       // IDID=-2  LARGER NMAX IS NEEDED,
                       // IDID=-3  STEP SIZE BECOMES TOO SMALL,
                       // IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
  );
  return idid;
}
} // namespace radau_detail
} // namespace fub
