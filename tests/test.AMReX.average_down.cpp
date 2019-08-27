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

#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>

#include <iostream>

int main(int argc, char** argv) {
  amrex::Initialize(argc, argv);
  {
    // Define test parameters

    const amrex::IntVect refine_ratio(2);
    const int max_box_size = 2;
    const int ghost_cell_width = 1;
    const int n_components = 1;
    const amrex::Box box({AMREX_D_DECL(0, 0, 0)}, {AMREX_D_DECL(3, 3, 3)});
    const amrex::RealBox xbox({AMREX_D_DECL(0.0, 0.0, 0.0)},
                              {AMREX_D_DECL(1.0, 1.0, 1.0)});

    // Define geometries

    const amrex::Geometry coarse_geom(box, &xbox);
    amrex::Geometry fine_geom = coarse_geom;
    fine_geom.refine(refine_ratio);

    // Define BoxArrays and DistributionMapping

    amrex::BoxArray coarse_ba(box);
    coarse_ba.maxSize(max_box_size);
    amrex::BoxArray fine_ba(coarse_ba.boxList());
    fine_ba.refine(refine_ratio);

    amrex::DistributionMapping coarse_dm(coarse_ba);
    amrex::DistributionMapping fine_dm(fine_ba);

    // Define MultiFabs and fill with 0s and 1s

    amrex::MultiFab coarse_mf(coarse_ba, coarse_dm, n_components,
                              ghost_cell_width);
    amrex::MultiFab fine_mf(fine_ba, fine_dm, n_components, ghost_cell_width);

    coarse_mf.setVal(0.0);
    fine_mf.setVal(1.0);

    // Average Down!

    amrex::average_down(fine_mf, coarse_mf, fine_geom, coarse_geom, 0,
                        n_components, refine_ratio);

    // Print FabArrays and watch ghost cell data

    for (amrex::MFIter mfi(coarse_mf); mfi.isValid(); ++mfi) {
      amrex::Print() << coarse_mf[mfi];
    }
  }
  amrex::Finalize();
}
