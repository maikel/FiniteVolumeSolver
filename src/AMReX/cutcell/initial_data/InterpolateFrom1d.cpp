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

#include "fub/AMReX/cutcell/initial_data/InterpolateFrom1d.hpp"

#include "fub/AMReX/ViewFArrayBox.hpp"
#include "fub/AMReX/ForEachIndex.hpp"
#include "fub/AMReX/ForEachFab.hpp"

#include "fub/output/Hdf5Handle.hpp"

namespace fub::amrex::cutcell {

InterpolateFrom1d::InterpolateFrom1d(const PerfectGas<AMREX_SPACEDIM>& equation,
                    const std::string& name)
  : equation_(equation), raw_prim_data_{}
{
  H5File file(H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  if (file < 0) {
    throw std::runtime_error("Could not open file " + name + ".");
  }
  if (H5Dataset dataset(H5Dopen(file, "DG_Solution", H5P_DEFAULT)); dataset < 0) {
    throw std::runtime_error("Could not open dataset name://DGSolution.");
  } else {
    std::array<hsize_t, 2> dims = {};
    H5Space dataspace(H5Dget_space(dataset));
    H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr);
    std::size_t size = dims[0] * dims[1];
    raw_prim_data_.resize(size);
    H5Dread(dataset, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw_prim_data_.data());
  }
}

void InterpolateFrom1d::InitializeData(::amrex::MultiFab& data, const ::amrex::Geometry& geom)
{
  auto& factory = static_cast<const ::amrex::EBFArrayBoxFactory&>(data.Factory());
  const ::amrex::MultiFab& alphas = factory.getVolFrac();
  [[maybe_unused]] static constexpr int n_vars = 6;
//  static constexpr int i_x = 0;
  [[maybe_unused]] static constexpr int i_density = 1;
  [[maybe_unused]] static constexpr int i_velocity_x = 2;
  [[maybe_unused]] static constexpr int i_pressure = 5;
//  const double dx = geom.CellSize(0);
  ForEachFab(fub::execution::openmp, data, [&](const ::amrex::MFIter& mfi) {
    mdspan<const double, 2> raw_prims(raw_prim_data_.data(), raw_prim_data_.size() / n_vars, n_vars);
    span<const double> xs = raw_prims.get_span().subspan(0, raw_prims.extent(0));
    ::amrex::Box tilebox = mfi.growntilebox();
    auto states = MakeView<Complete<PerfectGas<AMREX_SPACEDIM>>>(data[mfi], equation_, tilebox);
    double xlo[AMREX_SPACEDIM]{};
    double xhi[AMREX_SPACEDIM]{};
    const ::amrex::FArrayBox& alpha = alphas[mfi];
    ForEachIndex(tilebox, [&](auto... is) {
      ::amrex::IntVect iv{int(is)...};
      if (alpha(iv) > 0.0) {
      geom.LoNode(iv, xlo);
      geom.HiNode(iv, xhi);
      const double* xs_lo = std::lower_bound(xs.begin(), xs.end(), xlo[0]) + 1;
      const double* xs_up = std::lower_bound(xs.begin(), xs.end(), xhi[0]);
      if (xs.end() <= xs_up) {
        const std::ptrdiff_t ilast = raw_prims.extent(0) - 1;
        const double density = raw_prims(ilast, i_density);
        const double velocity = raw_prims(ilast, i_velocity_x);
        const double pressure = raw_prims(ilast, i_pressure);
        FUB_ASSERT(density > 0.0);
        FUB_ASSERT(pressure > 0.0);
        Complete<PerfectGas<AMREX_SPACEDIM>> state = equation_.CompleteFromPrim(density, {AMREX_D_DECL(velocity, 0, 0)}, pressure);
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> index{is...};
        Store(states, state, index);
      } else if (!(xs_lo < xs_up)) {
        const std::ptrdiff_t ifirst = xs_up - xs.begin();
        const double density = raw_prims(ifirst, i_density);
        const double velocity = raw_prims(ifirst, i_velocity_x);
        const double pressure = raw_prims(ifirst, i_pressure);
        FUB_ASSERT(density > 0.0);
        FUB_ASSERT(pressure > 0.0);
        Complete<PerfectGas<AMREX_SPACEDIM>> state = equation_.CompleteFromPrim(density, {AMREX_D_DECL(velocity, 0, 0)}, pressure);
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> index{is...};
        Store(states, state, index);
      } else {
        const auto n = xs_up - xs_lo;
        const auto i0 = xs_lo - xs.begin();

        const double* density_first = &raw_prims(i0, i_density);
        const double* velocity_first = &raw_prims(i0, i_velocity_x);
        const double* pressure_first = &raw_prims(i0, i_pressure);

        const double density = std::accumulate(density_first, density_first + n, 0.0) / n;
        const double velocity = std::accumulate(velocity_first, velocity_first + n, 0.0) / n;
        const double pressure = std::accumulate(pressure_first, pressure_first + n, 0.0) / n;
        FUB_ASSERT(density > 0.0);
        FUB_ASSERT(pressure > 0.0);

        Complete<PerfectGas<AMREX_SPACEDIM>> state = equation_.CompleteFromPrim(density, {AMREX_D_DECL(velocity, 0, 0)}, pressure);
        std::array<std::ptrdiff_t, AMREX_SPACEDIM> index{is...};
        Store(states, state, index);
      }
      }
    });
  });
}

}
