#include "fub/core/assert.hpp"
#include "fub/geometry/Polygon.hpp"

int main()
{
  fub::Polygon polygon({0.1, 0.2, 0.15, 0.1}, {0.1, 0.1, 0.2, 0.1});
  const double d1 = polygon.ComputeDistanceTo(0.09, 0.1);
  FUB_ASSERT(d1 < 0.0);
  const double d2 = polygon.ComputeDistanceTo(0.11, 0.105);
  FUB_ASSERT(d2 > 0.0);
}
