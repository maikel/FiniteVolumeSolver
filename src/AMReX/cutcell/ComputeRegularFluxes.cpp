namespace fub::amrex::cutcell {

template <typename FluxMethod>
void ComputeRegularFluxes(FluxMethod& fm, ::amrex::MultiFab& fluxes,
                          const ::amrex::MultiFab& cells,
                          const IntegratorContext& context, Duration dt,
                          double dx, Direction dir) {
    ForEachFab(execution::openmp, fluxes, [&](const ::amrex::MFIter& mfi) {
      const ::amrex::Box face_tilebox = mfi.growntilebox();
      const ::amrex::Box cell_validbox = scratch[mfi].box();
      auto [cell_box, face_box] =
          GetCellsAndFacesInStencilRange(cell_validbox, face_tilebox, gcw, dir);
      auto flux =
          MakeView<Conservative<Equation>>(fluxes[mfi], *equation, face_box);
      auto states =
          MakeView<const Complete<Equation>>(src[mfi], *equation, cell_box);
      const ::amrex::FabType type = context.GetFabType(level, mfi);
      if (type == ::amrex::FabType::singlevalued) {
        CutCellData<AMREX_SPACEDIM> geom = hierarchy.GetCutCellData(level, mfi);
        fm.ComputeRegularFluxes(flux, states, geom, dt, dx, dir);
      } else {
        fm.ComputeNumericFluxes(flux, states, dt, dx, dir);
      }
    });
}

} // namespace fub::amrex::cutcell