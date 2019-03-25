#include <fmt/format.h>
#include <limits>

int main()
{
  float eps = std::numeric_limits<float>::epsilon();
  float dt = 1.323e-6;
  float dx = 1.25e-5;
  float dt_over_eps = dt / eps;
  float F = 1.23456789e6;
  float epsF_over_dx = eps * F / dx;
  float G = dt_over_eps * epsF_over_dx;
  float H = dt * (F / dx);
  fmt::print("{:.12e}\n", G);
  fmt::print("{:.12e}\n", H);
}