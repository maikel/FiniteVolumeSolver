#ifndef FUB_AMREX_BOUNDARY_CONDITION_ISENTROPIC_HPP
#define FUB_AMREX_BOUNDARY_CONDITION_ISENTROPIC_HPP

#include "fub/AMReX/BoundaryCondition.hpp"
#include "fub/equations/IdealGasMix.hpp"

namespace fub::amrex {

class IsentropicBoundary {
public:
  IsentropicBoundary(const IdealGasMix<1>& eq, double outer_pressure,
                     Direction dir, int side);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom,
                    Duration dt, const GriddingAlgorithm&);

  void FillBoundary(::amrex::MultiFab& mf, const ::amrex::Geometry& geom);

private:
  IdealGasMix<1> equation_;
  double outer_pressure_;
  Direction dir_;
  int side_;
};

} // namespace fub::amrex

#endif