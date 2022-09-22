! radiation_disort_lw.F90 - Longwave DISORT solver
!
! (C) Copyright 2022- ECMWF.
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
! Modifications

module radiation_disort_lw

  public

contains

  !---------------------------------------------------------------------
  ! Longwave disort solver, in which clouds are assumed to fill
  ! the gridbox horizontally
  subroutine solver_disort_lw(nlev,istartcol,iendcol, &
       &  config, cloud, & 
       &  od, ssa, g, od_cloud, ssa_cloud, g_cloud, planck_hl, &
       &  emission, albedo, &
       &  flux)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook

    use radiation_constants, only      : Pi
    use radiation_config, only         : config_type
    use radiation_cloud, only          : cloud_type
    use radiation_flux, only           : flux_type, indexed_sum_profile
    !use radiation_lw_derivatives, only : calc_lw_derivatives_ica

    implicit none

    integer, parameter :: nstream = 4
    ! Number of angles for radiance calculation (not needed, so set to
    ! 1)
    integer, parameter :: nmu = nstream
    integer, parameter :: nphi = 1
    
    ! Inputs
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type),        intent(in) :: config
    type(cloud_type),         intent(in) :: cloud

    ! Gas and aerosol optical depth, single-scattering albedo and
    ! asymmetry factor at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, nlev, istartcol:iendcol) :: &
         &  od
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering, nlev, istartcol:iendcol) :: &
         &  ssa, g

    ! Cloud and precipitation optical depth, single-scattering albedo and
    ! asymmetry factor in each longwave band
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol)   :: &
         &  od_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering, &
         &  nlev,istartcol:iendcol) :: ssa_cloud, g_cloud

    ! Planck function at each half-level and the surface
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: &
         &  planck_hl
  
    ! Emission (Planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) &
         &  :: emission, albedo

    ! Output
    type(flux_type), intent(inout):: flux

    ! Local variables

    ! Optical depth and single scattering albedo
    real(jprb), dimension(nlev) :: od_in, ssa_in

    ! Coefficients of a Legendre decompositon of the phase function
    real(jprb), dimension(0:nstream,nlev) :: pf_mom

    ! Dummy arguments to DISORT
    !real(jprb) :: dummy(nlev), dummy2(1,1), dummy3(1,1,1)
    real(jprb) :: phi(nphi), mu(nmu)
    real(jprb) :: h_layer(0:nlev)
    real(jprb) :: rhoq(nstream/2, 0:nstream/2, 0:(nstream-1))
    real(jprb) :: rhou(nmu, 0:nstream/2, 0:(nstream-1))
    real(jprb) :: rho_accurate(nmu,nphi)
    real(jprb) :: bemst(nstream/2), emust(nmu)
    real(jprb) :: dfdt(nlev+1), uavg(nlev+1), uu(nmu,nlev+1,nphi)
    real(jprb) :: albmed(nmu), trnmed(nmu), flux_dn_dir(nlev+1), od_out(nlev+1)

    ! Planck functions per sterad (black-body flux divided by pi)
    real(jprb) :: planck_hl_per_sr(nlev+1), planck_surf_per_sr
    
    character*127 :: header
    
    ! Fluxes per g point
    real(jprb), dimension(config%n_g_lw, nlev+1) :: flux_up, flux_dn

    ! Is there any cloud in the profile?
    logical :: is_cloudy_profile

    ! Number of g points
    integer :: ng

    ! Loop indices for level and column
    integer :: jlev, jcol, jg

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_disort_lw:solver_disort_lw',0,hook_handle)

    ng = config%n_g_lw

    header = 'ecRad profile'
    !dummy = 0.5_jprb
    !dummy2 = 0.5_jprb
    !dummy3 = 0.5_jprb
    phi = 0.0_jprb
    mu = 0.0_jprb
    
    ! Loop through columns
    do jcol = istartcol,iendcol

      ! Is there any cloud in the profile?
      is_cloudy_profile = .false.
      do jlev = 1,nlev
        if (cloud%fraction(jcol,jlev) >= config%cloud_fraction_threshold) then
          is_cloudy_profile = .true.
          exit
        end if
      end do

      ! Clear sky calculation
      do jg = 1,ng
        od_in = od(jg,:,jcol)
        pf_mom = 0.0_jprb;
        pf_mom(0,:) = 1.0_jprb
        if (config%do_lw_aerosol_scattering) then
          ssa_in = ssa(jg,:,jcol)
          pf_mom(1,:) = g(jg,:,jcol)*3.0_jprb;
        else
          ssa_in = 0.0_jprb
        end if

        planck_hl_per_sr = planck_hl(jg,:,jcol) / Pi
        planck_surf_per_sr = emission(jg,jcol)/(Pi*(1.0_jprb-albedo(jg,jcol)))
        
        call disort(nlev, nstream, nstream, nmu, nphi, nlev+1, &
             &  .false., .false., 0, .true., spread(.false.,1,5), &
             &  .true., .true., .false., .false., od_in, &
             &  ssa_in, pf_mom, planck_hl_per_sr, -1.0_jprb, -1.0_jprb, &
             &  od_out, 1.0_jprb, 0.0_jprb, mu, phi, 0.0_jprb, &
             &  0.0_jprb, albedo(jg,jcol), planck_surf_per_sr, &
             &  0.0_jprb, 1.0_jprb, 6371.0_jprb, h_layer, &
             &  rhoq, rhou, rho_accurate, bemst, emust, &
             &  0.0_jprb, header, flux_dn_dir, flux_dn(jg,:), flux_up(jg,:), &
             &  dfdt, uavg, uu, albmed, trnmed)
      
      end do

      ! Sum over g-points to compute broadband fluxes
      flux%lw_up_clear(jcol,:) = sum(flux_up,1)
      flux%lw_dn_clear(jcol,:) = sum(flux_dn,1)
      ! Store surface spectral downwelling fluxes
      flux%lw_dn_surf_clear_g(:,jcol) = flux_dn(:,nlev+1)
      
      ! Save the spectral fluxes if required
      if (config%do_save_spectral_flux) then
        call indexed_sum_profile(flux_up, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_up_clear_band(:,jcol,:))
        call indexed_sum_profile(flux_dn, config%i_spec_from_reordered_g_lw, &
             &                   flux%lw_dn_clear_band(:,jcol,:))
      end if

    end do

    if (lhook) call dr_hook('radiation_disort_lw:solver_disort_lw',1,hook_handle)
    
  end subroutine solver_disort_lw

end module radiation_disort_lw
