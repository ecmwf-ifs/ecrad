! tcrad_clear_sky_solver.h - Calculate clear-sky fluxes/radiances in TCRAD
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

module tcrad_clear_sky_solver

contains

  subroutine calc_clear_sky_flux(nspec, nlev, surf_emission, &
       &  surf_albedo, planck_hl, &
       &  od_clear, flux_up, flux_dn, n_angles_per_hem)

    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook
    use tcrad_layer_solutions, only : calc_clear_sky_trans_source, &
         &  gauss_legendre, LW_DIFFUSIVITY, MAX_GAUSS_LEGENDRE_POINTS

    implicit none

    real(jprb), parameter :: LW_MU = 1.0_jprb / LW_DIFFUSIVITY

    ! Inputs

    ! Number of spectral intervals and levels. Note that all
    ! level-dependent variables count down from top-of-atmosphere.
    integer(jpim), intent(in) :: nspec, nlev

    ! Surface upwards emission, in W m-2 (i.e. emissivity multiplied
    ! by Planck function at the surface skin temperature) integrated
    ! across each spectral interval, and albedo in the same intervals
    real(jprb), intent(in), dimension(nspec) :: surf_emission, surf_albedo

    ! Planck function integrated over each spectral interval at each
    ! half-level, in W m-2 (i.e. the flux emitted by a horizontal
    ! black-body surface)
    real(jprb), intent(in), dimension(nspec,nlev+1) :: planck_hl

    ! Layer optical depth of gas and aerosol
    real(jprb), intent(in), dimension(nspec,nlev) :: od_clear

    ! Outputs

    ! Upwelling and downwelling fluxes in each spectral interval at
    ! each half-level (W m-2)
    real(jprb), intent(out), dimension(nspec,nlev+1) :: flux_up, flux_dn

    ! Optional inputs

    ! Number of angles to compute radiances per hemisphere, for
    ! example, 2 results in the delta-2-plus-4 algorithm recommended
    ! by Fu et al. (1997). A value of 0 (the default) indicates to use
    ! the output from the two-stream Tripleclouds flux calculation
    ! directly.
    integer,    intent(in), optional :: n_angles_per_hem

    ! Local variables

    ! Transmittance of each layer
    real(jprb), dimension(nspec,nlev) :: transmittance

    ! Rate of emission up from the top or down through the base of each
    ! layer (W m-2)
    real(jprb), dimension(nspec,nlev) :: source_up, source_dn

    ! Surface upwelling flux (W m-2)
    real(jprb), dimension(nspec) :: flux_up_surf

    ! Upwelling and downwelling fluxes at the top and base of each
    ! layer in each region, in W m-2
    !  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_base, flux_dn_base
    !  real(jprb), dimension(nspec,NREGION,nlev) :: flux_up_top, flux_dn_top

    ! Gauss-Legendre points and weights for sampling cosine of zenith
    ! angle distribution
    real(jprb), dimension(MAX_GAUSS_LEGENDRE_POINTS) :: mu_list, weight_list

    ! Actual weight used accounts for projection into horizontal area
    real(jprb) ::  weight

    ! Local version of an optional argument
    integer(jpim) :: n_angles_per_hem_local

    ! Loop index for stream
    integer(jpim) :: jstream

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_clear_sky_flux',0,hook_handle)

    ! Store local value for number of angles per hemisphere. Note that
    ! values of 0 and 1 have the same effect.
    if (present(n_angles_per_hem)) then
      n_angles_per_hem_local = min(abs(n_angles_per_hem), MAX_GAUSS_LEGENDRE_POINTS)
    else
      n_angles_per_hem_local = 0
    end if

    if (n_angles_per_hem_local < 2) then
      ! We need source up and down
      call calc_clear_sky_trans_source(nspec, nlev, LW_MU, &
           &  planck_hl, od_clear, transmittance, &
           &  source_up=source_up, source_dn=source_dn)
    else
      ! We need only the source down, in order to do a single downward
      ! calculation to obtain the reflected flux at the surface
      call calc_clear_sky_trans_source(nspec, nlev, LW_MU, &
           &  planck_hl, od_clear, transmittance, &
           &  source_dn=source_dn)
    end if

    flux_up = 0.0_jprb
    flux_dn = 0.0_jprb

    call calc_clear_sky_radiance_dn(nspec, nlev, 1.0_jprb, &
         &  transmittance, source_dn, flux_dn)
    flux_up_surf = surf_emission + surf_albedo*flux_dn(:,nlev+1)

    if (n_angles_per_hem_local > 1) then
      ! Fu et al. (1997) method: pass N beams through the
      ! atmosphere using the two-stream solution as the scattering
      ! source function
      if (present(n_angles_per_hem)) then
        ! Negative input values for n_angles_per_hem lead to
        ! alternative quadrature, but n_angles_per_hem_local has been
        ! forced to be positive
        call gauss_legendre(n_angles_per_hem, mu_list, weight_list)
      else
        call gauss_legendre(n_angles_per_hem_local, mu_list, weight_list)
      end if

      flux_dn = 0.0_jprb

      do jstream = 1,n_angles_per_hem_local
        weight = weight_list(jstream)*mu_list(jstream) &
             &  / sum(weight_list(1:n_angles_per_hem_local) &
             &          * mu_list(1:n_angles_per_hem_local))
        ! Radiances are computed in pairs: up and down with same
        ! absolute zenith angle
        call calc_clear_sky_trans_source(nspec,nlev, mu_list(jstream), &
             &  planck_hl, od_clear, transmittance, source_dn=source_dn, &
             &  source_up=source_up)
        call calc_clear_sky_radiance_dn(nspec, nlev, &
             &  weight, &
             &  transmittance, source_dn, flux_dn)
        call calc_clear_sky_radiance_up(nspec, nlev, &
             &  weight, flux_up_surf, &
             &  transmittance, source_up, flux_up)
      end do

    else ! n_angles_per_hem_local == 0 or 1
      call calc_clear_sky_radiance_up(nspec, nlev, &
           &  1.0_jprb, flux_up_surf, &
           &  transmittance, source_up, flux_up)
    end if

    if (lhook) call dr_hook('tcrad:calc_clear_sky_flux',1,hook_handle)

  end subroutine calc_clear_sky_flux


  !---------------------------------------------------------------------
  ! Compute clear-sky upward radiance profile by solving the
  ! Schwarzschild radiative transfer equation assuming the source term
  ! to vary linearly with optical depth in each layer. This routine adds
  ! to any existing radiance profile, useful if used as part of a
  ! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
  ! et al. (1997).
  subroutine calc_clear_sky_radiance_up(nspec, nlev, &
       &  weight, surf_up, &
       &  transmittance, source_up, radiance_up)

    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev

    ! Weight sources by this amount
    real(jprb), intent(in) :: weight

    ! Surface upwelling flux in W m-2
    real(jprb), intent(in),  dimension(nspec) :: surf_up

    ! Transmittance of each layer in the direction of the radiance; this
    ! does not include diffuse transmittance, i.e. rays that may be
    ! scattered as they pass through the layer.  We use an assumed-shape
    ! array because this argument may have been sliced from a
    ! multi-region array in which case the level dimension will not be
    ! contiguous.
    real(jprb), intent(in),  dimension(:,:) :: transmittance

    ! Upward source from the top of the layer in the direction of the
    ! radiance, which may include Planck emission, and scattering. We
    ! use an assumed-shape array because this argument may have been
    ! sliced from a multi-region array in which case the level dimension
    ! will not be contiguous.
    real(jprb), intent(in),  dimension(:,:) :: source_up ! (nspec,nlev)

    ! Output

    ! Upward radiance profile: note that we add to any existing radiance
    ! profile, useful when summing over multiple angles to get a flux
    real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_up

    ! Local variables

    ! Spectral radiance at an interface
    real(jprb), dimension(nspec) :: radiance

    ! Loop index for level
    integer(jpim) :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_clear_sky_radiance_up',0,hook_handle)

    ! Set surface upward radiance
    radiance = weight * surf_up

    ! Save radiance profile averaged over regions, adding to existing
    ! values
    radiance_up(:,nlev+1) = radiance_up(:,nlev+1) + radiance

    do jlev = nlev,1,-1
      ! Solution to Schwarzschild equation
      radiance = transmittance(:,jlev)*radiance + weight * source_up(:,jlev)
      ! Save radiances
      radiance_up(:,jlev) = radiance_up(:,jlev) + radiance
    end do

    if (lhook) call dr_hook('tcrad:calc_radiance_clear_sky_up',1,hook_handle)

  end subroutine calc_clear_sky_radiance_up


  !---------------------------------------------------------------------
  ! Compute clear-sky downward radiance profile by solving the
  ! Schwarzschild radiative transfer equation assuming the source term
  ! to vary linearly with optical depth in each layer. This routine adds
  ! to any existing radiance profile, useful if used as part of a
  ! multi-stream flux calculation, e.g. the delta-2-plus-4 method of Fu
  ! et al. (1997).
  subroutine calc_clear_sky_radiance_dn(nspec, nlev, &
       &  weight, transmittance, source_dn, radiance_dn)

    use parkind1, only           : jpim, jprb
    use yomhook,  only           : lhook, dr_hook

    implicit none

    ! Inputs

    ! Number of spectral intervals and levels
    integer(jpim), intent(in) :: nspec, nlev

    ! Weight sources by this amount
    real(jprb), intent(in) :: weight

    ! Transmittance of each layer in the direction of the radiance; this
    ! does not include diffuse transmittance, i.e. rays that may be
    ! scattered as they pass through the layer.  We use an assumed-shape
    ! array because this argument may have been sliced from a
    ! multi-region array in which case the level dimension will not be
    ! contiguous.
    real(jprb), intent(in),  dimension(:,:) :: transmittance ! (nspec,nlev)

    ! Down source from the base of the layer in the direction of the
    ! radiance, which may include Planck emission, and scattering.  We
    ! use an assumed-shape array because this argument may have been
    ! sliced from a multi-region array in which case the level dimension
    ! will not be contiguous.
    real(jprb), intent(in),  dimension(:,:) :: source_dn ! (nspec,nlev)

    ! Output

    ! Downward radiance profile: note that we add to any existing radiance
    ! profile, useful when summing over multiple angles to get a flux
    real(jprb), intent(inout), dimension(nspec,nlev+1) :: radiance_dn

    ! Local variables

    ! Spectral radiance at an interface
    real(jprb), dimension(nspec) :: radiance

    ! Loop index for level
    integer(jpim) :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_clear_sky_radiance_dn',0,hook_handle)

    ! Start with zero at TOA
    radiance = 0.0_jprb

    do jlev = 1,nlev
      ! Solution to Schwarzschild equation
      radiance = transmittance(:,jlev)*radiance + weight * source_dn(:,jlev)
      ! Save radiances
      radiance_dn(:,jlev+1) = radiance_dn(:,jlev+1) + radiance
    end do

    if (lhook) call dr_hook('tcrad:calc_clear_sky_radiance_dn',1,hook_handle)

  end subroutine calc_clear_sky_radiance_dn

end module tcrad_clear_sky_solver

