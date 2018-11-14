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

#include "fub/ode_solver/CVodeSolver.hpp"

extern "C" {
#include <cvode/cvode.h>
#include <cvodes/cvodes_direct.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
}

#include <cmath>

namespace fub {

struct CVodeSolver::CVodeInternalData_ {
  void* cvode_mem;
  N_Vector y;
  N_Vector abstol;
  SUNMatrix A;
  SUNLinearSolver solver;
};

void CVodeSolver::ReleaseCVodeMemory::operator()(CVodeInternalData_* data) const
    noexcept {
  if (data != nullptr) {
    SUNLinSolFree(data->solver);
    SUNMatDestroy(data->A);
    N_VDestroy(data->abstol);
    N_VDestroy(data->y);
    CVodeFree(&data->cvode_mem);
  }
  delete data;
}

CVodeSolver::CVodeSolver(int size)
    : data_(new CVodeInternalData_{CVodeCreate(CV_BDF, CV_NEWTON)}) {
  data_->abstol = N_VNew_Serial(size);
  data_->y = N_VNew_Serial(size);
  data_->A = SUNDenseMatrix(size, size);
  data_->solver = SUNDenseLinearSolver(data_->y, data_->A);
  CVodeSetMaxNumSteps(data_->cvode_mem, 10000);
  CVodeSetMaxOrd(data_->cvode_mem, 5);
}

void CVodeSolver::integrate(system_type system, span<double> y_0, double t0,
                            double dt, feedback_type* feedback,
                            jacobian_type* jacobian) const {
  struct IntegrateData_ {
    system_type* system;
    jacobian_type* jacobian;
    const root_function_type* root;
    int root_size;
  };
  IntegrateData_ data{&system, jacobian, &root_function_, root_size_};
  auto cvode_rhs = [](realtype t, N_Vector x, N_Vector xdot, void* user_data) {
    FUB_ASSERT(user_data != nullptr);
    system_type& system =
        *reinterpret_cast<const IntegrateData_*>(user_data)->system;
    span<double> dydt(NV_DATA_S(xdot), N_VGetLength_Serial(xdot));
    span<const double> y(NV_DATA_S(x), N_VGetLength_Serial(x));
    system(dydt, y, t);
    return CV_SUCCESS;
  };
  auto cvode_jacobian = [](realtype t, N_Vector x, N_Vector fx, SUNMatrix Jac,
                           void* user_data, N_Vector tmp1, N_Vector tmp2,
                           N_Vector tmp3) {
    jacobian_type& jacobian =
        *reinterpret_cast<const IntegrateData_*>(user_data)->jacobian;
    span<double> jac(SUNDenseMatrix_Data(Jac), SUNDenseMatrix_LData(Jac));
    span<const double> y(N_VGetArrayPointer(x), N_VGetLength_Serial(x));
    jacobian(jac, y, t);
    return CV_SUCCESS;
  };
  auto cvode_root = [](realtype t, N_Vector x, realtype* r, void* user_data) {
    const root_function_type& root_function =
        *reinterpret_cast<const IntegrateData_*>(user_data)->root;
    const int root_size =
        reinterpret_cast<const IntegrateData_*>(user_data)->root_size;
    span<double> roots(r, root_size);
    span<const double> y(NV_DATA_S(x), N_VGetLength_Serial(x));
    root_function(roots, y, t);
    return CV_SUCCESS;
  };
  void* cvode_mem = data_->cvode_mem;
  std::copy(y_0.begin(), y_0.end(), N_VGetArrayPointer(data_->y));
  double* abstol = N_VGetArrayPointer(data_->abstol);
  std::transform(y_0.begin(), y_0.end(), abstol, [&](double yi) {
    return std::max({relative_tolerance_ * std::abs(yi), 1e-8});
  });

  CVodeInit(cvode_mem, cvode_rhs, t0, data_->y);
  CVodeSetUserData(cvode_mem, &data);
  CVDlsSetLinearSolver(data_->cvode_mem, data_->solver, data_->A);
  CVodeSVtolerances(cvode_mem, relative_tolerance_, data_->abstol);
  if (root_function_) {
    CVodeRootInit(cvode_mem, root_size_, cvode_root);
  }
  if (jacobian) {
    CVDlsSetJacFn(cvode_mem, cvode_jacobian);
  }
  const double tend = t0 + dt;
  double t = t0;
  if (!feedback) {
    CVode(cvode_mem, tend, data_->y, &t, CV_NORMAL);
  } else {
    while (t < tend) {
      CVode(cvode_mem, tend, data_->y, &t, CV_ONE_STEP);
      (*feedback)(y_0, t);
    }
  }
  std::copy_n(N_VGetArrayPointer(data_->y), y_0.size(), y_0.begin());
}

} // namespace fub
