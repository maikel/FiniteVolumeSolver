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

#include <AMReX.H>
#include <AMReX_MultiFab.H>

int main(int argc, char** argv)
{
  amrex::Initialize(argc, argv);
  {
    const amrex::Box box(amrex::IntVect(0), amrex::IntVect(32));
    amrex::BoxArray ba(box);
    amrex::DistributionMapping dm(ba);

    amrex::MultiFab mf(ba, dm, 1, 4);

    for (amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi) {
      amrex::Print() << "----------------------------------------\n";
      amrex::Print() << "validbox(): " << mfi.validbox() << '\n';
      amrex::Print() << "growntilebox() " << mfi.growntilebox() << '\n';
      amrex::Print() << "grownnodaltilebox(0): " << mfi.grownnodaltilebox(0) << '\n';
      amrex::Print() << "grownnodaltilebox(1): " << mfi.grownnodaltilebox(1) << '\n';
      amrex::Print() << "grownnodaltilebox(2): " << mfi.grownnodaltilebox(2) << '\n';
    }
  }
  amrex::Finalize();
}