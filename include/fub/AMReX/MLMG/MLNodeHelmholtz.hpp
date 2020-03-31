// Copyright (c) 2019 Stefan Vater
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

#ifndef AMREX_ML_NODEHELMHOLTZ_HPP
#define AMREX_ML_NODEHELMHOLTZ_HPP

#include <AMReX.H>
#include <AMReX_MLNodeLinOp.H>

namespace amrex {

struct MLNodeHelmholtz : public MLNodeLinOp {

  virtual ~MLNodeHelmholtz(){};

  virtual void setSigma(int amrlev, const MultiFab& a_sigma) = 0;

  virtual void setAlpha(int amrlev, const MultiFab& a_alpha) = 0;

  virtual void compDivergence(const Vector<MultiFab*>& rhs,
                              const Vector<MultiFab*>& vel) = 0;

  virtual void compRHS(const Vector<MultiFab*>& rhs,
                       const Vector<MultiFab*>& vel,
                       const Vector<const MultiFab*>& rhnd,
                       const Vector<MultiFab*>& rhcc) = 0;

};

} // namespace amrex

#endif
