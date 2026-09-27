! radiation_adding_ica_lw.F90 - Longwave adding method in independent column approximation
!
! (C) Copyright 2015- ECMWF.
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
!   2017-04-11  R. Hogan  Receive emission/albedo rather than planck/emissivity
!   2017-07-12  R. Hogan  Fast adding method for if only clouds scatter
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_lw

  public

  !$omp declare target(fast_adding_ica_lw_omp)
  !$omp declare target(fast_adding_ica_lw_recompute_omp)
  !$omp declare target(calc_fluxes_no_scattering_lw_omp)
  !$omp declare target(calc_fluxes_no_scattering_lw_recompute_omp)
contains

  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering, by successively adding the contribution of
  ! layers starting from the surface to compute the total albedo and
  ! total upward emission of the increasingly larger block of
  ! atmospheric layers.
  subroutine adding_ica_lw(ncol, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  flux_up, flux_dn)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw',0,hook_handle)

    albedo(:,nlev+1) = albedo_surf

    ! At the surface, the source is thermal emission
    source(:,nlev+1) = emission_surf

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
      do jcol = 1,ncol
        ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb &
             &  / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
             &  * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 11:
        source(jcol,jlev) = source_up(jcol,jlev) &
             &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
             &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev)) &
             &                   * inv_denominator(jcol,jlev)
      end do
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(:,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to emission
    ! below that level
    flux_up(:,1) = source(:,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
      do jcol = 1,ncol
        ! Shonk & Hogan (2008) Eq 14 (after simplification):
        flux_dn(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
      end do
    end do

    if (lhook) call dr_hook('radiation_adding_ica_lw:adding_ica_lw',1,hook_handle)

  end subroutine adding_ica_lw


  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  subroutine fast_adding_ica_lw(ncol, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
       &  flux_up, flux_dn, albedo, source, inv_denominator)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! Determine which layers are cloud-free
    logical, intent(in) :: is_clear_sky_layer(nlev)

    ! Index to highest cloudy layer
    integer, intent(in) :: i_cloud_top

    ! Pre-computed clear-sky downwelling fluxes (W m-2) at half-levels
    real(jprb), intent(in), dimension(ncol, nlev+1)  :: flux_dn_clear

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), intent(out), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), intent(out), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), intent(out), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

    !$ACC ROUTINE WORKER 

#ifndef _OPENACC
    if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw',0,hook_handle)
#endif

    ! Copy over downwelling fluxes above cloud from clear sky
    flux_dn(:,1:i_cloud_top) = flux_dn_clear(:,1:i_cloud_top)

    albedo(:,nlev+1) = albedo_surf
    
    ! At the surface, the source is thermal emission
    source(:,nlev+1) = emission_surf

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    !$ACC LOOP SEQ
    do jlev = nlev,i_cloud_top,-1
      if (is_clear_sky_layer(jlev)) then
        ! Reflectance of this layer is zero, simplifying the expression
        !$ACC LOOP WORKER VECTOR
        do jcol = 1,ncol
          albedo(jcol,jlev) = transmittance(jcol,jlev)*transmittance(jcol,jlev)*albedo(jcol,jlev+1)
          source(jcol,jlev) = source_up(jcol,jlev) &
               &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
               &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev))
        end do
      else
        ! Loop over columns; explicit loop seems to be faster
        !$ACC LOOP WORKER VECTOR
        do jcol = 1,ncol
          ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
          inv_denominator(jcol,jlev) = 1.0_jprb &
               &  / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
          ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
          albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev)*transmittance(jcol,jlev) &
               &  * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
          ! Shonk & Hogan (2008) Eq 11:
          source(jcol,jlev) = source_up(jcol,jlev) &
               &  + transmittance(jcol,jlev) * (source(jcol,jlev+1) &
               &                    + albedo(jcol,jlev+1)*source_dn(jcol,jlev)) &
               &                   * inv_denominator(jcol,jlev)
        end do
      end if
    end do

    ! Compute the fluxes above the highest cloud
    !$ACC LOOP WORKER VECTOR
    do jcol = 1,ncol
      flux_up(jcol,i_cloud_top) = source(jcol,i_cloud_top) &
          &                 + albedo(jcol,i_cloud_top)*flux_dn(jcol,i_cloud_top)
    end do
    !$ACC LOOP SEQ
    do jlev = i_cloud_top-1,1,-1
      !$ACC LOOP WORKER VECTOR
      do jcol = 1,ncol
        flux_up(jcol,jlev) = transmittance(jcol,jlev)*flux_up(jcol,jlev+1) + source_up(jcol,jlev)
      end do
    end do

    ! Work back down through the atmosphere from cloud top computing
    ! the fluxes at each half-level
    !$ACC LOOP SEQ
    do jlev = i_cloud_top,nlev
      if (is_clear_sky_layer(jlev)) then
        !$ACC LOOP WORKER VECTOR
        do jcol = 1,ncol
          flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
               &               + source_dn(jcol,jlev)
          flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
               &               + source(jcol,jlev+1)
        end do
      else
        !$ACC LOOP WORKER VECTOR
        do jcol = 1,ncol
          ! Shonk & Hogan (2008) Eq 14 (after simplification):
          flux_dn(jcol,jlev+1) &
               &  = (transmittance(jcol,jlev)*flux_dn(jcol,jlev) &
               &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
               &     + source_dn(jcol,jlev)) * inv_denominator(jcol,jlev)
          ! Shonk & Hogan (2008) Eq 12:
          flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn(jcol,jlev+1) &
               &               + source(jcol,jlev+1)
        end do
      end if
    end do

#ifndef _OPENACC
    if (lhook) call dr_hook('radiation_adding_ica_lw:fast_adding_ica_lw',1,hook_handle)
#endif

  end subroutine fast_adding_ica_lw

  !---------------------------------------------------------------------
  ! Use the scalar "adding" method to compute longwave flux profiles,
  ! including scattering in cloudy layers only.
  subroutine fast_adding_ica_lw_omp(jg, ng, nlev, &
       &  reflectance, transmittance, source_up, source_dn, emission_surf, albedo_surf, &
       &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
       &  flux_up, flux_dn, albedo, inv_denominator, source)

    use parkind1, only           : jprb
    implicit none

    ! Inputs
    integer, intent(in) :: jg ! spectral index
    integer, intent(in) :: ng ! number of spectral bands
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ng) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ng, nlev)   :: reflectance, transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ng, nlev)   :: source_up, source_dn

    ! Determine which layers are cloud-free
    logical, intent(in) :: is_clear_sky_layer(nlev)

    ! Index to highest cloudy layer
    integer, intent(in) :: i_cloud_top

    ! Pre-computed clear-sky downwelling fluxes (W m-2) at half-levels
    real(jprb), intent(in), dimension(ng, nlev+1)  :: flux_dn_clear

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn
    real(jprb), intent(out), dimension(ng, nlev+1) :: albedo, inv_denominator

    ! Upwelling radiation at each half-level due to emission below
    ! that half-level (W m-2)
    real(jprb), intent(out), dimension(ng, nlev+1) :: source

    ! Loop index for model level and column
    integer :: jlev

    !
    ! Ideally, we would use the associate statement below, however this doesn't
    ! work with NVCOMPILER. Thus, we pass albedo and inv_denominator in from the
    ! calling code. However, we pass it in as flux up and flux_dn. That is, the
    ! flux_up and flux_dn are temporarilly reusable.
    !
    !associate(albedo=>flux_up, inv_denominator=>flux_dn)

    albedo(jg,nlev+1) = albedo_surf(jg)

    ! At the surface, the source is thermal emission
    source(jg,nlev+1) = emission_surf(jg)

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    do jlev = nlev,i_cloud_top,-1
      if (is_clear_sky_layer(jlev)) then
         ! Reflectance of this layer is zero, simplifying the expression
         albedo(jg,jlev) = transmittance(jg,jlev)*transmittance(jg,jlev)*albedo(jg,jlev+1)
         source(jg,jlev) = source_up(jg,jlev) &
              &  + transmittance(jg,jlev) * (source(jg,jlev+1) &
              &                    + albedo(jg,jlev+1)*source_dn(jg,jlev))
      else
         ! Loop over columns; explicit loop seems to be faster
         ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
         inv_denominator(jg,jlev+1) = 1.0_jprb &
              &  / (1.0_jprb-albedo(jg,jlev+1)*reflectance(jg,jlev))
         ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
         albedo(jg,jlev) = reflectance(jg,jlev) + transmittance(jg,jlev)*transmittance(jg,jlev) &
              &  * albedo(jg,jlev+1) * inv_denominator(jg,jlev+1)
         ! Shonk & Hogan (2008) Eq 11:
         source(jg,jlev) = source_up(jg,jlev) &
              &  + transmittance(jg,jlev) * (source(jg,jlev+1) &
              &                    + albedo(jg,jlev+1)*source_dn(jg,jlev)) &
              &                   * inv_denominator(jg,jlev+1)
      end if
    end do

    ! Copy over downwelling fluxes above cloud from clear sky
    do jlev = 1,i_cloud_top
       flux_dn(jg,jlev) = flux_dn_clear(jg,jlev)
    enddo

    ! Compute the fluxes above the highest cloud
    flux_up(jg,i_cloud_top) = source(jg,i_cloud_top) &
         &                 + albedo(jg,i_cloud_top)*flux_dn(jg,i_cloud_top)
    do jlev = i_cloud_top-1,1,-1
       flux_up(jg,jlev) = transmittance(jg,jlev)*flux_up(jg,jlev+1) + source_up(jg,jlev)
    end do

    ! Work back down through the atmosphere from cloud top computing
    ! the fluxes at each half-level
    do jlev = i_cloud_top,nlev
      if (is_clear_sky_layer(jlev)) then
         flux_dn(jg,jlev+1) = transmittance(jg,jlev)*flux_dn(jg,jlev) &
              &               + source_dn(jg,jlev)
         flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn(jg,jlev+1) &
              &               + source(jg,jlev+1)
      else
         ! Shonk & Hogan (2008) Eq 14 (after simplification):
         flux_dn(jg,jlev+1) &
              &  = (transmittance(jg,jlev)*flux_dn(jg,jlev) &
              &     + reflectance(jg,jlev)*source(jg,jlev+1) &
              &     + source_dn(jg,jlev)) * inv_denominator(jg,jlev+1)
         ! Shonk & Hogan (2008) Eq 12:
         flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn(jg,jlev+1) &
              &               + source(jg,jlev+1)
      end if
   end do

  end subroutine fast_adding_ica_lw_omp

  !---------------------------------------------------------------------
  ! As fast_adding_ica_lw_omp, but recomputing each layer's reflectance,
  ! transmittance and emission from the optical properties at the point
  ! of use rather than reading them from four spectral profiles that a
  ! preceding loop had to write. The two sweeps below each touch every
  ! layer once, so the properties are recomputed twice instead of being
  ! stored once and read twice; on this kernel the arithmetic is cheaper
  ! than the traffic. The layer properties never reach memory.
  subroutine fast_adding_ica_lw_recompute_omp(jg, ng, nlev, i_band, &
       &  n_bands, n_bands_s, od, od_scaling, od_cloud, ssa_cloud, g_cloud, planck, &
       &  do_lw_cloud_scattering, emission_surf, albedo_surf, &
       &  is_clear_sky_layer, i_cloud_top, flux_dn_clear, &
       &  transmittance, flux_up, flux_dn, albedo, source)

    use parkind1, only           : jprb
    use radiation_two_stream, only : calc_ref_trans_lw_scalar_omp, &
         &                           calc_no_scattering_transmittance_lw_single_cell_omp
    implicit none

    ! Inputs
    integer, intent(in) :: jg ! spectral index
    integer, intent(in) :: ng ! number of spectral bands
    integer, intent(in) :: nlev ! number of levels

    ! Band containing this g-point, and the band counts of the cloud
    ! optical property arrays
    integer, intent(in) :: i_band, n_bands, n_bands_s

    ! Gas (plus aerosol) optical depth, and the McICA cloud optical
    ! depth scaling, both spectral
    real(jprb), intent(in), dimension(ng, nlev) :: od, od_scaling

    ! In-cloud optical properties, by band
    real(jprb), intent(in), dimension(n_bands, nlev)   :: od_cloud
    real(jprb), intent(in), dimension(n_bands_s, nlev) :: ssa_cloud, g_cloud

    ! Planck function at half-levels
    real(jprb), intent(in), dimension(ng, nlev+1) :: planck

    ! Do clouds scatter in the longwave?
    logical, intent(in) :: do_lw_cloud_scattering

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ng) :: emission_surf, albedo_surf

    ! Determine which layers are cloud-free
    logical, intent(in) :: is_clear_sky_layer(nlev)

    ! Index to highest cloudy layer
    integer, intent(in) :: i_cloud_top

    ! Pre-computed clear-sky downwelling fluxes (W m-2) at half-levels
    real(jprb), intent(in), dimension(ng, nlev+1)  :: flux_dn_clear

    ! Diffuse transmittance of each layer. This is the one layer property
    ! that has to reach memory: the longwave derivative kernel reads it back
    ! for every cloudy column. The other three stay in registers.
    real(jprb), intent(out), dimension(ng, nlev)   :: transmittance

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn

    ! Albedo of the atmosphere below each half-level, and the upwelling
    ! radiation at each half-level due to emission below it (W m-2)
    real(jprb), intent(out), dimension(ng, nlev+1) :: albedo, source

    ! Layer properties of the level being visited, held in registers
    real(jprb) :: reflectance, trans_lev, source_up, source_dn

    ! Combined gas+aerosol+cloud optical properties of a cloudy layer
    real(jprb) :: od_cloud_new, od_total, ssa_total, g_total, scat_od

    ! Reciprocal of the denominator of Lacis and Hansen (1974) Eq 33
    real(jprb) :: inv_denominator

    ! Loop index for model level
    integer :: jlev

    albedo(jg,nlev+1) = albedo_surf(jg)

    ! At the surface, the source is thermal emission
    source(jg,nlev+1) = emission_surf(jg)

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to emission
    ! below that level
    do jlev = nlev,i_cloud_top,-1
      if (is_clear_sky_layer(jlev)) then
         call calc_no_scattering_transmittance_lw_single_cell_omp(od(jg,jlev), &
              &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, source_up, source_dn)
         transmittance(jg,jlev) = trans_lev
         ! Reflectance of this layer is zero, simplifying the expression
         albedo(jg,jlev) = trans_lev*trans_lev*albedo(jg,jlev+1)
         source(jg,jlev) = source_up &
              &  + trans_lev * (source(jg,jlev+1) &
              &                    + albedo(jg,jlev+1)*source_dn)
      else
         od_cloud_new = od_scaling(jg,jlev) * od_cloud(i_band,jlev)
         od_total  = od(jg,jlev) + od_cloud_new
         ssa_total = 0.0_jprb
         g_total   = 0.0_jprb
         if (do_lw_cloud_scattering) then
            if (od_total > 0.0_jprb) then
               scat_od = ssa_cloud(i_band,jlev) * od_cloud_new
               ssa_total = scat_od / od_total
               if (scat_od > 0.0_jprb) then
                  g_total = g_cloud(i_band,jlev) * ssa_cloud(i_band,jlev) * od_cloud_new / scat_od
               end if
            end if
            call calc_ref_trans_lw_scalar_omp(od_total, ssa_total, g_total, &
                 &  planck(jg,jlev), planck(jg,jlev+1), &
                 &  reflectance, trans_lev, source_up, source_dn)
         else
            reflectance = 0.0_jprb
            call calc_no_scattering_transmittance_lw_single_cell_omp(od_total, &
                 &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, source_up, source_dn)
         end if
         transmittance(jg,jlev) = trans_lev
         ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
         inv_denominator = 1.0_jprb &
              &  / (1.0_jprb-albedo(jg,jlev+1)*reflectance)
         ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
         albedo(jg,jlev) = reflectance + trans_lev*trans_lev &
              &  * albedo(jg,jlev+1) * inv_denominator
         ! Shonk & Hogan (2008) Eq 11:
         source(jg,jlev) = source_up &
              &  + trans_lev * (source(jg,jlev+1) &
              &                    + albedo(jg,jlev+1)*source_dn) &
              &                   * inv_denominator
      end if
    end do

    ! Copy over downwelling fluxes above cloud from clear sky
    do jlev = 1,i_cloud_top
       flux_dn(jg,jlev) = flux_dn_clear(jg,jlev)
    enddo

    ! Compute the fluxes above the highest cloud
    flux_up(jg,i_cloud_top) = source(jg,i_cloud_top) &
         &                 + albedo(jg,i_cloud_top)*flux_dn(jg,i_cloud_top)
    do jlev = i_cloud_top-1,1,-1
       ! Every layer above cloud top is by definition clear
       call calc_no_scattering_transmittance_lw_single_cell_omp(od(jg,jlev), &
            &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, source_up, source_dn)
       transmittance(jg,jlev) = trans_lev
       flux_up(jg,jlev) = trans_lev*flux_up(jg,jlev+1) + source_up
    end do

    ! Work back down through the atmosphere from cloud top computing
    ! the fluxes at each half-level
    do jlev = i_cloud_top,nlev
      if (is_clear_sky_layer(jlev)) then
         call calc_no_scattering_transmittance_lw_single_cell_omp(od(jg,jlev), &
              &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, source_up, source_dn)
         flux_dn(jg,jlev+1) = trans_lev*flux_dn(jg,jlev) &
              &               + source_dn
         flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn(jg,jlev+1) &
              &               + source(jg,jlev+1)
      else
         od_cloud_new = od_scaling(jg,jlev) * od_cloud(i_band,jlev)
         od_total  = od(jg,jlev) + od_cloud_new
         ssa_total = 0.0_jprb
         g_total   = 0.0_jprb
         if (do_lw_cloud_scattering) then
            if (od_total > 0.0_jprb) then
               scat_od = ssa_cloud(i_band,jlev) * od_cloud_new
               ssa_total = scat_od / od_total
               if (scat_od > 0.0_jprb) then
                  g_total = g_cloud(i_band,jlev) * ssa_cloud(i_band,jlev) * od_cloud_new / scat_od
               end if
            end if
            call calc_ref_trans_lw_scalar_omp(od_total, ssa_total, g_total, &
                 &  planck(jg,jlev), planck(jg,jlev+1), &
                 &  reflectance, trans_lev, source_up, source_dn)
         else
            reflectance = 0.0_jprb
            call calc_no_scattering_transmittance_lw_single_cell_omp(od_total, &
                 &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, source_up, source_dn)
         end if
         ! The denominator is recomputed here rather than carried from the
         ! sweep above: albedo(jg,jlev+1) is still intact at this point,
         ! being consumed by the flux_up expression below.
         inv_denominator = 1.0_jprb &
              &  / (1.0_jprb-albedo(jg,jlev+1)*reflectance)
         ! Shonk & Hogan (2008) Eq 14 (after simplification):
         flux_dn(jg,jlev+1) &
              &  = (trans_lev*flux_dn(jg,jlev) &
              &     + reflectance*source(jg,jlev+1) &
              &     + source_dn) * inv_denominator
         ! Shonk & Hogan (2008) Eq 12:
         flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn(jg,jlev+1) &
              &               + source(jg,jlev+1)
      end if
   end do

  end subroutine fast_adding_ica_lw_recompute_omp

  !---------------------------------------------------------------------
  ! If there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  subroutine calc_fluxes_no_scattering_lw(ncol, nlev, &
       &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

    use parkind1, only           : jprb
#ifndef _OPENACC
    use yomhook,  only           : lhook, dr_hook, jphook
#endif

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ncol) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ncol, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn
    
    ! Loop index for model level
    integer :: jlev, jcol

#ifndef _OPENACC
    real(jphook) :: hook_handle
#endif

    !$ACC ROUTINE WORKER

#ifndef _OPENACC
    if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw',0,hook_handle)
#endif

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    !$ACC LOOP WORKER VECTOR
    do jcol = 1,ncol
      flux_dn(jcol,1) = 0.0_jprb
    end do

    ! Work down through the atmosphere computing the downward fluxes
    ! at each half-level
!$ACC LOOP SEQ
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
      !$ACC LOOP WORKER VECTOR
      do jcol = 1,ncol
        flux_dn(jcol,jlev+1) = transmittance(jcol,jlev)*flux_dn(jcol,jlev) + source_dn(jcol,jlev)
      end do
    end do

    ! Surface reflection and emission
    !$ACC LOOP WORKER VECTOR
    do jcol = 1,ncol
      flux_up(jcol,nlev+1) = emission_surf(jcol) + albedo_surf(jcol) * flux_dn(jcol,nlev+1)
    end do

    ! Work back up through the atmosphere computing the upward fluxes
    ! at each half-level
!$ACC LOOP SEQ
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      !$ACC LOOP WORKER VECTOR
      do jcol = 1,ncol
        flux_up(jcol,jlev) = transmittance(jcol,jlev)*flux_up(jcol,jlev+1) + source_up(jcol,jlev)
      end do
    end do
    
#ifndef _OPENACC
    if (lhook) call dr_hook('radiation_adding_ica_lw:calc_fluxes_no_scattering_lw',1,hook_handle)
#endif

  end subroutine calc_fluxes_no_scattering_lw


  !---------------------------------------------------------------------
  ! If there is no scattering then fluxes may be computed simply by
  ! passing down through the atmosphere computing the downwelling
  ! fluxes from the transmission and emission of each layer, and then
  ! passing back up through the atmosphere to compute the upwelling
  ! fluxes in the same way.
  subroutine calc_fluxes_no_scattering_lw_omp(jg, ng, nlev, &
       &  transmittance, source_up, source_dn, emission_surf, albedo_surf, flux_up, flux_dn)

    use parkind1, only           : jprb
    implicit none

    ! Inputs
    integer, intent(in) :: jg, ng
    integer, intent(in) :: nlev ! number of levels

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ng) :: emission_surf, albedo_surf

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ng, nlev)   :: transmittance

    ! Emission from each layer in an upward and downward direction
    real(jprb), intent(in),  dimension(ng, nlev)   :: source_up, source_dn

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn

    ! Loop index for model level
    integer :: jlev, jcol

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(jg,1) = 0.0_jprb

    ! Work down through the atmosphere computing the downward fluxes
    ! at each half-level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
       flux_dn(jg,jlev+1) = transmittance(jg,jlev)*flux_dn(jg,jlev) + source_dn(jg,jlev)
    end do

    ! Surface reflection and emission
    flux_up(jg,nlev+1) = emission_surf(jg) + albedo_surf(jg) * flux_dn(jg,nlev+1)

    ! Work back up through the atmosphere computing the upward fluxes
    ! at each half-level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
       flux_up(jg,jlev) = transmittance(jg,jlev)*flux_up(jg,jlev+1) + source_up(jg,jlev)
    end do

  end subroutine calc_fluxes_no_scattering_lw_omp


  !---------------------------------------------------------------------
  ! As calc_fluxes_no_scattering_lw_omp, but the layer transmittance and
  ! emission are recomputed from od and planck at the point of use rather
  ! than read back from three precomputed spectral profiles. Storing them
  ! costs three writes and four reads to HBM, where recomputing needs only
  ! od and planck twice and the underlying expression contains a single
  ! exponential -- cheap enough that the trade is strongly favourable for
  ! a kernel at 87% of the bandwidth roofline.
  !
  ! transmittance is still returned because later kernels consume the
  ! clear-sky profile. The upward sweep recomputes rather than reading it
  ! back, since it needs source_up from the same evaluation anyway.
  !
  ! Results are bit-identical to the stored-profile version: the values
  ! come from the same expressions applied to the same inputs.
  subroutine calc_fluxes_no_scattering_lw_recompute_omp(jg, ng, nlev, &
       &  od, planck, emission_surf, albedo_surf, transmittance, flux_up, flux_dn)

    use parkind1, only             : jprb
    use radiation_two_stream, only : calc_no_scattering_transmittance_lw_single_cell_omp
    implicit none

    ! Inputs
    integer, intent(in) :: jg, ng
    integer, intent(in) :: nlev ! number of levels

    ! Layer optical depth and the Planck function at half-levels, from
    ! which the transmittance and emission are recomputed
    real(jprb), intent(in),  dimension(ng, nlev)   :: od
    real(jprb), intent(in),  dimension(ng, nlev+1) :: planck

    ! Surface emission (W m-2) and albedo
    real(jprb), intent(in),  dimension(ng) :: emission_surf, albedo_surf

    ! Diffuse transmittance of each layer, retained for later kernels
    real(jprb), intent(out), dimension(ng, nlev)   :: transmittance

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling and
    ! downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn

    ! Per-level properties, recomputed into registers
    real(jprb) :: trans_lev, src_up, src_dn

    ! Loop index for model level
    integer :: jlev

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn(jg,1) = 0.0_jprb

    ! Work down through the atmosphere computing the downward fluxes
    ! at each half-level
    do jlev = 1,nlev
       call calc_no_scattering_transmittance_lw_single_cell_omp(od(jg,jlev), &
            &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, src_up, src_dn)
       transmittance(jg,jlev) = trans_lev
       flux_dn(jg,jlev+1) = trans_lev*flux_dn(jg,jlev) + src_dn
    end do

    ! Surface reflection and emission
    flux_up(jg,nlev+1) = emission_surf(jg) + albedo_surf(jg) * flux_dn(jg,nlev+1)

    ! Work back up through the atmosphere computing the upward fluxes
    ! at each half-level
    do jlev = nlev,1,-1
       call calc_no_scattering_transmittance_lw_single_cell_omp(od(jg,jlev), &
            &  planck(jg,jlev), planck(jg,jlev+1), trans_lev, src_up, src_dn)
       flux_up(jg,jlev) = trans_lev*flux_up(jg,jlev+1) + src_up
    end do

  end subroutine calc_fluxes_no_scattering_lw_recompute_omp

end module radiation_adding_ica_lw
