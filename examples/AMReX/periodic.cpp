#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

int main(int argc, char** argv) {
  Initialize(argc, argv);
  {
    Box domain(IntVect{0, 0}, IntVect{3, 3});
    BoxArray ba(domain);
    DistributionMapping dm{ba};
    int ncomp = 1;
    IntVect ngrow{2, 0};
    MultiFab src(ba, dm, ncomp, 0);
    MultiFab dest(ba, dm, ncomp, ngrow);
    // Set Iner Values
    FArrayBox& fab = src[0];
    int counter = 0;
    for (BoxIterator bit(domain); bit.ok(); ++bit) {
      fab(bit()) = counter++;
    }
    fab.setFormat(FABio::FAB_ASCII);
    fab.writeOn(std::cout);
    std::cout << '\n';

    // Fill Boundaries
    double xlo[] = {-1.0, -1.0};
    double xup[] = {+1.0, +1.0};
    RealBox coords(xlo, xup);
    int is_periodic[] = {1, 1};
    Geometry geom(domain, &coords, -1, is_periodic);
    Vector<MultiFab*> smf{&src};
    Vector<Real> stime{0};
    Vector<BCRec> bcr(AMREX_SPACEDIM * 2);
    auto noop = [](auto&&...) {};
    PhysBCFunct<decltype(noop)> no_condition(geom, bcr, noop);
    FillPatchSingleLevel(dest, 0, smf, stime, 0, 0, ncomp, geom, no_condition, 0);

    // Print Result
    FArrayBox& dfab = dest[0];
    dfab.setFormat(FABio::FAB_ASCII);
    dfab.writeOn(std::cout);
    std::cout << '\n';
  }
  Finalize();
}