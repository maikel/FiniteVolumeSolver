#include "fub/grid/AMReX/ViewFArrayBox.hpp"

namespace fub {
namespace amrex {

std::array<std::ptrdiff_t, AMREX_SPACEDIM> AsArray(const ::amrex::IntVect& vec) {
  std::array<std::ptrdiff_t, AMREX_SPACEDIM> array;
  for (std::size_t i = 0; i < static_cast<std::size_t>(AMREX_SPACEDIM); ++i) {
    array[i] = vec[static_cast<int>(i)];
  }
  return array;
}

}
}