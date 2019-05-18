#include "fub/equations/Advection.hpp"
#include "fub/equations/PerfectGas.hpp"
#include "fub/equations/ShallowWater.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

using namespace fub;

TEST_CASE("Construct pointer to beginning") {
  std::vector<double> buffer(3 * 32 * 32);
  BasicView<ShallowWater::Complete, layout_left> view;
  view.height = PatchDataView<double, 2>{
      mdspan<double, 2>(buffer.data(), dynamic_extents<2>{32, 32}), {}};
  view.momentum = PatchDataView<double, 3>{
      mdspan<double, 3>(buffer.data() + 32 * 32, dynamic_extents<3>{32, 32, 2}),
      {}};

  SECTION("Construction") {
    ViewPointer pointer = Begin(view);
    REQUIRE(pointer.height == view.height.Span().begin());
    REQUIRE(pointer.momentum.first == view.momentum.Span().begin());
    REQUIRE(pointer.momentum.second == view.momentum.Stride(2));
  }

  SECTION("Load state from beginning") {
    ViewPointer pointer = Begin(AsConst(view));
    buffer[0] = 42.0;
    ShallowWater::Complete q;
    Load(q, pointer);
    REQUIRE(q.height == 42.0);
  }
}

TEST_CASE("Construct pointer to beginning on PerfectGas3d") {
  const std::size_t size = 8 * 8 * 8;
  std::vector<double> buffer(7 * size);
  BasicView<PerfectGas<3>::Complete, layout_left> view;
  view.density = PatchDataView<double, 3>{
      mdspan<double, 3>(buffer.data(), dynamic_extents<3>{8, 8, 8}), {}};
  view.momentum = PatchDataView<double, 4>{
      mdspan<double, 4>(buffer.data() + size, dynamic_extents<4>{8, 8, 8, 3}),
      {}};
  view.energy = PatchDataView<double, 3>{
      mdspan<double, 3>(buffer.data() + 4 * size, dynamic_extents<3>{8, 8, 8}),
      {}};
  view.pressure = PatchDataView<double, 3>{
      mdspan<double, 3>(buffer.data() + 5 * size, dynamic_extents<3>{8, 8, 8}),
      {}};
  view.speed_of_sound = PatchDataView<double, 3>{
      mdspan<double, 3>(buffer.data() + 6 * size, dynamic_extents<3>{8, 8, 8}),
      {}};

  SECTION("Construction") {
    ViewPointer pointer = Begin(view);
    REQUIRE(pointer.density == view.density.Span().begin());
    REQUIRE(pointer.momentum.first == view.momentum.Span().begin());
    REQUIRE(pointer.momentum.second == view.momentum.Stride(3));
    REQUIRE(pointer.energy == view.energy.Span().begin());
    REQUIRE(pointer.pressure == view.pressure.Span().begin());
  }

  SECTION("Load state from beginning") {
    ViewPointer pointer = Begin(AsConst(view));
    buffer[0] = 42.0;
    PerfectGas<3>::Complete q;
    Load(q, pointer);
    REQUIRE(q.density == 42.0);
  }
}
