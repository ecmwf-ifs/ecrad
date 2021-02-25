! radiation_tcrad_lw.F90 - Interface to TCRAD solver
!
! (C) Copyright 2021- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!

module radiation_tcrad_lw

  public

contains

  subroutine solver_tcrad_lw(nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

!    use radiation_io, only             : nulout
    use radiation_config, only         : config_type, IPdfShapeGamma
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    use radiation_cloudless_lw, only   : solver_cloudless_lw
    use tcrad_3region_solver, only     : calc_flux

    implicit none

    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth of each layer at each longwave
    ! g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od

    ! Gas and aerosol single-scattering albedo and asymmetry factor,
    ! only if longwave scattering by aerosols is to be represented
    real(jprb), intent(in), &
         &  dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: ssa, g

    ! Cloud and precipitation optical depth of each layer in each
    ! longwave band
    real(jprb), intent(in) :: od_cloud(config%n_bands_lw,nlev,istartcol:iendcol)

    ! Cloud and precipitation single-scattering albedo and asymmetry
    ! factor, only if longwave scattering by clouds is to be
    ! represented
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &                            nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function (emitted flux from a black body) at half levels
    ! and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Regridded cloud properties
    real(jprb), dimension(config%n_g_lw,nlev) :: od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid

    real(jprb), dimension(config%n_g_lw,nlev+1) :: flux_up, flux_dn

    real(jprb) :: hook_handle

    integer :: jcol

    if (lhook) call dr_hook('radiation_tcrad_lw:solver_tcrad_lw',0,hook_handle)

    ! Compute clear-sky fluxes using the cloudless solver
    call solver_cloudless_lw(nlev,istartcol,iendcol, &
         &  config, od, ssa, g, &
         &  planck_hl, emission, albedo, flux)
 
    flux%lw_up(istartcol:iendcol,:) = 0.0_jprb
    flux%lw_dn(istartcol:iendcol,:) = 0.0_jprb

    if (.not. config%do_lw_cloud_scattering) then
      ! FIX we have allocated memory that is not used
      ssa_cloud_regrid = 0.0_jprb
      g_cloud_regrid   = 0.0_jprb
    end if

    do jcol = istartcol,iendcol

      od_cloud_regrid    = od_cloud(config%i_band_from_reordered_g_lw,:,jcol)
      if (config%do_lw_cloud_scattering) then
        ssa_cloud_regrid = ssa_cloud(config%i_band_from_reordered_g_lw,:,jcol)
        g_cloud_regrid   = g_cloud(config%i_band_from_reordered_g_lw,:,jcol)
      end if

      call calc_flux(config%n_g_lw, nlev, emission(:,jcol), albedo(:,jcol), planck_hl(:,:,jcol), &
           &         cloud%fraction(jcol,:), cloud%fractional_std(jcol,:), &
           &         od(:,:,jcol), od_cloud_regrid, ssa_cloud_regrid, g_cloud_regrid, &
           &         cloud%overlap_param(jcol,:), &
           &         flux_up, flux_dn, &
           &         n_stream_per_hem=config%n_stream_per_hem_lw, do_3d=config%do_3d_effects)
      flux%lw_up(jcol,:) = sum(flux_up,1)
      flux%lw_dn(jcol,:) = sum(flux_dn,1)

      if (config%do_save_spectral_flux) then
        call indexed_sum_profile(flux_up, &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_up_band(:,jcol,:))
        call indexed_sum_profile(flux_dn, &
             &           config%i_spec_from_reordered_g_lw, &
             &           flux%lw_dn_band(:,jcol,:))
      end if
     
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_g(:,jcol) = flux_dn(:,nlev+1)

    end do

    if (lhook) call dr_hook('radiation_tcrad_lw:solver_tcrad_lw',1,hook_handle)

  end subroutine solver_tcrad_lw

end module radiation_tcrad_lw
