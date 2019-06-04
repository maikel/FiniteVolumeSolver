#ifndef FUB_AMREX_GEOMETRY_HPP
#define FUB_AMREX_GEOMETRY_HPP

#include <AMReX.H>

#include <array>

namespace fub::amrex {

template <typename Geom> class Geometry : private Geom {
public:
  Geometry(const Geom& base) : Geom(base) {}
  Geometry(Geom&& base) : Geom(std::move(base)) {}

  Geom& Base() noexcept { return *this; }
  const Geom& Base() const noexcept { return *this; }

  double operator()(AMREX_D_DECL(double x, double y, double z)) const {
    return Base().ComputeDistanceTo(AMREX_D_DECL(x, y, z));
  }

  double operator()(const std::array<double, AMREX_SPACEDIM>& coords) const {
    return Base().ComputeDistanceTo(
        AMREX_D_DECL(coords[0], coords[1], coords[2]));
  }
};

} // namespace fub::amrex

#endif