#include "fub/equations/PerfectGas.hpp"
#include "fub/equations/perfect_gas/ExactRiemannSolver.hpp"
#include "fub/equations/perfect_gas/HllemMethod.hpp"
#include "fub/flux_method/MusclHancockMethod2.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

using namespace fub;

template <typename Eq, typename Limiter>
using MusclHancockPrim =
    MusclHancock2<Eq,
                 PrimitiveGradient<Eq, fub::CentralDifferenceGradient<Limiter>>,
                 PrimitiveReconstruction<Eq>, perfect_gas::Hllem<Eq>>;

TEST_CASE("Construct Primitive State") {
  PerfectGas<1> equation{};
  std::array<Primitive<PerfectGas<1>>, 4> w{};
  std::array<Complete<PerfectGas<1>>, 4> q{};

  w.fill(Primitive<PerfectGas<1>>{equation});
  
  REQUIRE(w[0].density == 0.0);
  REQUIRE(w[0].velocity[0] == 0.0);
  REQUIRE(w[0].pressure == 0.0);

  for (int i = 0; i < 2; ++i) {
    w[i].density = 1.0;
    w[i].velocity[0] = 0.0;
    w[i].pressure = 1.0;
  }

  for (int i = 2; i < 4; ++i) {
    w[i].density = 0.125;
    w[i].velocity[0] = 0.0;
    w[i].pressure = 0.1;
  }

  for (int i = 0; i < 4; ++i) {
    CompleteFromPrim(equation, q[i], w[i]);
  }

  const double dx = 0.1;
  using MusclPrim = MusclHancockPrim<PerfectGas<1>, fub::VanLeerLimiter>;
  FluxMethod<MusclPrim> muscl_prim{equation};
  Duration dt(muscl_prim.ComputeStableDt(q, dx, Direction::X));
  Conservative<PerfectGas<1>> f{equation};
  muscl_prim.ComputeNumericFlux(f, q, dt, dx, Direction::X);
}
