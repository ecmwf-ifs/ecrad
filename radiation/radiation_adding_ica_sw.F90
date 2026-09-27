! radiation_adding_ica_sw.F90 - Shortwave adding method in independent column approximation
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
!   2017-10-23  R. Hogan  Renamed single-character variables

module radiation_adding_ica_sw

  public

  !$omp declare target(adding_ica_sw_omp)
  !$omp declare target(adding_ica_sw_recompute_omp)
  !$omp declare target(adding_ica_sw_recompute_total_omp)
contains

  subroutine adding_ica_sw(ncol, nlev, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct, &
       &  albedo, source, inv_denominator)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook, jphook

    implicit none

    ! Inputs
    integer, intent(in) :: ncol ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(ncol)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ncol)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! Cosine of the solar zenith angle
    real(jprb), intent(in)                           :: cos_sza

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: reflectance, transmittance

    ! Fraction of direct-beam solar radiation entering the top of a
    ! layer that is reflected back up or scattered forward into the
    ! diffuse stream at the base of the layer
    real(jprb), intent(in),  dimension(ncol, nlev)   :: ref_dir, trans_dir_diff

    ! Direct transmittance, i.e. fraction of direct beam that
    ! penetrates a layer without being scattered or absorbed
    real(jprb), intent(in),  dimension(ncol, nlev)   :: trans_dir_dir

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(ncol, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct
    
    ! Albedo of the entire earth/atmosphere system below each half
    ! level
    real(jprb), intent(out), dimension(ncol, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), intent(out), dimension(ncol, nlev+1) :: source

    ! Equal to 1/(1-albedo*reflectance)
    real(jprb), intent(out), dimension(ncol, nlev)   :: inv_denominator

    ! Loop index for model level and column
    integer :: jlev, jcol

    real(jphook) :: hook_handle

#ifndef _OPENACC
    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',0,hook_handle)
#endif

    !$ACC ROUTINE WORKER

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    flux_dn_direct(:,1) = incoming_toa
    !$ACC LOOP SEQ
    do jlev = 1,nlev
      flux_dn_direct(:,jlev+1) = flux_dn_direct(:,jlev)*trans_dir_dir(:,jlev)
    end do

    albedo(:,nlev+1) = albedo_surf_diffuse

    ! At the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(:,nlev+1) = albedo_surf_direct * flux_dn_direct(:,nlev+1) * cos_sza

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
!$ACC LOOP SEQ
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(:,jlev) = 1.0_jprb / (1.0_jprb-albedo(:,jlev+1)*reflectance(:,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.
      !$ACC LOOP WORKER VECTOR
      do jcol = 1,ncol
        ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
        inv_denominator(jcol,jlev) = 1.0_jprb / (1.0_jprb-albedo(jcol,jlev+1)*reflectance(jcol,jlev))
        ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
        albedo(jcol,jlev) = reflectance(jcol,jlev) + transmittance(jcol,jlev) * transmittance(jcol,jlev) &
             &                                     * albedo(jcol,jlev+1) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 11:
        source(jcol,jlev) = ref_dir(jcol,jlev)*flux_dn_direct(jcol,jlev) &
             &  + transmittance(jcol,jlev)*(source(jcol,jlev+1) &
             &        + albedo(jcol,jlev+1)*trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) &
             &  * inv_denominator(jcol,jlev)
      end do
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn_diffuse(:,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(:,1) = source(:,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
!$ACC LOOP SEQ
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
      !$ACC LOOP WORKER VECTOR
      do jcol = 1,ncol
        ! Shonk & Hogan (2008) Eq 14 (after simplification):
        flux_dn_diffuse(jcol,jlev+1) &
             &  = (transmittance(jcol,jlev)*flux_dn_diffuse(jcol,jlev) &
             &     + reflectance(jcol,jlev)*source(jcol,jlev+1) &
             &     + trans_dir_diff(jcol,jlev)*flux_dn_direct(jcol,jlev)) * inv_denominator(jcol,jlev)
        ! Shonk & Hogan (2008) Eq 12:
        flux_up(jcol,jlev+1) = albedo(jcol,jlev+1)*flux_dn_diffuse(jcol,jlev+1) &
             &            + source(jcol,jlev+1)
        flux_dn_direct(jcol,jlev) = flux_dn_direct(jcol,jlev)*cos_sza
      end do
    end do
    flux_dn_direct(:,nlev+1) = flux_dn_direct(:,nlev+1)*cos_sza

#ifndef _OPENACC
    if (lhook) call dr_hook('radiation_adding_ica_sw:adding_ica_sw',1,hook_handle)
#endif

  end subroutine adding_ica_sw

  subroutine adding_ica_sw_omp(jg, ng, nlev, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct, albedo, inv_denominator,&
       &  source)

    use parkind1, only           : jprb
    implicit none

    ! Inputs
    integer, intent(in) :: jg, ng ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(ng)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ng)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! Cosine of the solar zenith angle
    real(jprb), intent(in)                           :: cos_sza

    ! Diffuse reflectance and transmittance of each layer
    real(jprb), intent(in),  dimension(ng, nlev)   :: reflectance, transmittance

    ! Fraction of direct-beam solar radiation entering the top of a
    ! layer that is reflected back up or scattered forward into the
    ! diffuse stream at the base of the layer
    real(jprb), intent(in),  dimension(ng, nlev)   :: ref_dir, trans_dir_diff

    ! Direct transmittance, i.e. fraction of direct beam that
    ! penetrates a layer without being scattered or absorbed
    real(jprb), intent(in),  dimension(ng, nlev)   :: trans_dir_dir

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct

    real(jprb), intent(out), dimension(ng, nlev+1) :: albedo, inv_denominator

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), intent(out), dimension(ng, nlev+1) :: source

    ! Loop index for model level and column
    integer :: jlev

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    flux_dn_direct(jg,1) = incoming_toa(jg)
    do jlev = 1,nlev
      flux_dn_direct(jg,jlev+1) = flux_dn_direct(jg,jlev)*trans_dir_dir(jg,jlev)
    end do

    !
    ! Ideally, we would use the associate statement below, however this doesn't
    ! work with NVCOMPILER. Thus, we pass albedo and inv_denominator in from the
    ! calling code however, we pass it in as flux up and flux_dn_diffuse. That is,
    ! the flux_up and flux_dn_diffuse are temporarilly reusable.
    !
    !associate(albedo=>flux_up, inv_denominator=>flux_dn_diffuse)

    albedo(jg,nlev+1) = albedo_surf_diffuse(jg)

    ! At the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(jg,nlev+1) = albedo_surf_direct(jg) * flux_dn_direct(jg,nlev+1) * cos_sza

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = nlev,1,-1
      ! Next loop over columns. We could do this by indexing the
      ! entire inner dimension as follows, e.g. for the first line:
      !   inv_denominator(jg,jlev) = 1.0_jprb / (1.0_jprb-albedo(jg,jlev+1)*reflectance(jg,jlev))
      ! and similarly for subsequent lines, but this slows down the
      ! routine by a factor of 2!  Rather, we do it with an explicit
      ! loop.

      ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
       inv_denominator(jg,jlev+1) = 1.0_jprb / (1.0_jprb-albedo(jg,jlev+1)*reflectance(jg,jlev))
       ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
       albedo(jg,jlev) = reflectance(jg,jlev) + transmittance(jg,jlev) * transmittance(jg,jlev) &
            &                                     * albedo(jg,jlev+1) * inv_denominator(jg,jlev+1)
       ! Shonk & Hogan (2008) Eq 11:
       source(jg,jlev) = ref_dir(jg,jlev)*flux_dn_direct(jg,jlev) &
            &  + transmittance(jg,jlev)*(source(jg,jlev+1) &
            &        + albedo(jg,jlev+1)*trans_dir_diff(jg,jlev)*flux_dn_direct(jg,jlev)) &
            &  * inv_denominator(jg,jlev+1)
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn_diffuse(jg,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(jg,1) = source(jg,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
! Added for DWD (2020)
!NEC$ outerloop_unroll(8)
    do jlev = 1,nlev
       ! Shonk & Hogan (2008) Eq 14 (after simplification):
       flux_dn_diffuse(jg,jlev+1) &
            &  = (transmittance(jg,jlev)*flux_dn_diffuse(jg,jlev) &
            &     + reflectance(jg,jlev)*source(jg,jlev+1) &
            &     + trans_dir_diff(jg,jlev)*flux_dn_direct(jg,jlev)) * inv_denominator(jg,jlev+1)
       ! Shonk & Hogan (2008) Eq 12:
       flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn_diffuse(jg,jlev+1) &
            &            + source(jg,jlev+1)
       flux_dn_direct(jg,jlev) = flux_dn_direct(jg,jlev)*cos_sza
    end do

    flux_dn_direct(jg,nlev+1) = flux_dn_direct(jg,nlev+1)*cos_sza
    !end associate

  end subroutine adding_ica_sw_omp


  !---------------------------------------------------------------------
  ! As adding_ica_sw_omp, but the layer two-stream properties are
  ! recomputed from od/ssa/asymmetry at the point of use rather than read
  ! from five precomputed spectral profiles. Those profiles are pure
  ! per-thread scratch with a lifetime of one call, yet at IFS
  ! resolutions they are far too large to keep in registers or LDS, so
  ! storing them costs five writes and eight reads to HBM per call. The
  ! two sweeps below recompute instead, which trades roughly twice the
  ! two-stream arithmetic for about a third less memory traffic -- the
  ! right way round for a kernel sitting at 88% of the bandwidth
  ! roofline with an arithmetic intensity near one.
  !
  ! Results are bit-identical to adding_ica_sw_omp: the recomputed values
  ! come from the same expressions applied to the same inputs.
  subroutine adding_ica_sw_recompute_omp(jg, ng, nlev, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  od, ssa, asymmetry, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct, albedo, inv_denominator, &
       &  source)

    use parkind1, only             : jprb
    use radiation_two_stream, only : calc_ref_trans_sw_scalar_omp
    implicit none

    ! Inputs
    integer, intent(in) :: jg, ng ! number of columns (may be spectral intervals)
    integer, intent(in) :: nlev ! number of levels

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(ng)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ng)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! Cosine of the solar zenith angle
    real(jprb), intent(in)                           :: cos_sza

    ! Layer optical depth, single scattering albedo and asymmetry factor,
    ! from which the two-stream properties are recomputed
    real(jprb), intent(in),  dimension(ng, nlev)   :: od, ssa, asymmetry

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct

    real(jprb), intent(out), dimension(ng, nlev+1) :: albedo, inv_denominator

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), intent(out), dimension(ng, nlev+1) :: source

    ! Layer two-stream properties, recomputed per level into registers
    real(jprb) :: reflectance, transmittance
    real(jprb) :: ref_dir, trans_dir_diff, trans_dir_dir

    ! Loop index for model level
    integer :: jlev

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere. Only the
    ! unscattered transmittance is needed here, so it is evaluated on its
    ! own rather than through the full two-stream solution.
    flux_dn_direct(jg,1) = incoming_toa(jg)
    do jlev = 1,nlev
      trans_dir_dir = max(-max(od(jg,jlev) * (1.0_jprb/cos_sza),0.0_jprb),-1000.0_jprb)
      trans_dir_dir = exp(trans_dir_dir)
      flux_dn_direct(jg,jlev+1) = flux_dn_direct(jg,jlev)*trans_dir_dir
    end do

    albedo(jg,nlev+1) = albedo_surf_diffuse(jg)

    ! At the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(jg,nlev+1) = albedo_surf_direct(jg) * flux_dn_direct(jg,nlev+1) * cos_sza

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
    do jlev = nlev,1,-1
      call calc_ref_trans_sw_scalar_omp(cos_sza, od(jg,jlev), ssa(jg,jlev), &
           &  asymmetry(jg,jlev), reflectance, transmittance, ref_dir, &
           &  trans_dir_diff, trans_dir_dir)

      ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
       inv_denominator(jg,jlev+1) = 1.0_jprb / (1.0_jprb-albedo(jg,jlev+1)*reflectance)
       ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
       albedo(jg,jlev) = reflectance + transmittance * transmittance &
            &                                     * albedo(jg,jlev+1) * inv_denominator(jg,jlev+1)
       ! Shonk & Hogan (2008) Eq 11:
       source(jg,jlev) = ref_dir*flux_dn_direct(jg,jlev) &
            &  + transmittance*(source(jg,jlev+1) &
            &        + albedo(jg,jlev+1)*trans_dir_diff*flux_dn_direct(jg,jlev)) &
            &  * inv_denominator(jg,jlev+1)
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn_diffuse(jg,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(jg,1) = source(jg,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
       call calc_ref_trans_sw_scalar_omp(cos_sza, od(jg,jlev), ssa(jg,jlev), &
            &  asymmetry(jg,jlev), reflectance, transmittance, ref_dir, &
            &  trans_dir_diff, trans_dir_dir)

       ! Shonk & Hogan (2008) Eq 14 (after simplification):
       flux_dn_diffuse(jg,jlev+1) &
            &  = (transmittance*flux_dn_diffuse(jg,jlev) &
            &     + reflectance*source(jg,jlev+1) &
            &     + trans_dir_diff*flux_dn_direct(jg,jlev)) * inv_denominator(jg,jlev+1)
       ! Shonk & Hogan (2008) Eq 12:
       flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn_diffuse(jg,jlev+1) &
            &            + source(jg,jlev+1)
       flux_dn_direct(jg,jlev) = flux_dn_direct(jg,jlev)*cos_sza
    end do

    flux_dn_direct(jg,nlev+1) = flux_dn_direct(jg,nlev+1)*cos_sza

  end subroutine adding_ica_sw_recompute_omp

  !---------------------------------------------------------------------
  ! As adding_ica_sw_omp, but combining the gas, aerosol and cloud
  ! optical properties and solving the two-stream equations at the point
  ! of use, so that the five layer properties never reach memory. The
  ! first sweep needs only the unscattered transmittance, which depends
  ! on the optical depth alone, so it is evaluated directly there rather
  ! than through the full two-stream solution.
  !
  ! This is for the case where delta-Eddington scaling is applied to the
  ! gases before this point; adding_ica_sw_omp is still used when the
  ! scaling is applied to the merged cloud-aerosol-gas mixture, since
  ! then the clear layers of a cloudy profile take delta-scaled values
  ! computed by the clear-sky solver.
  subroutine adding_ica_sw_recompute_total_omp(jg, ng, nlev, i_band, &
       &  n_bands, n_bands_s, incoming_toa, &
       &  albedo_surf_diffuse, albedo_surf_direct, cos_sza, &
       &  od, ssa, asymmetry, od_scaling, od_cloud, ssa_cloud, g_cloud, &
       &  frac, cloud_fraction_threshold, &
       &  flux_up, flux_dn_diffuse, flux_dn_direct, albedo, source)

    use parkind1, only             : jprb
    use radiation_two_stream, only : calc_ref_trans_sw_scalar_omp
    implicit none

    ! Inputs
    integer, intent(in) :: jg, ng ! spectral index and number of g-points
    integer, intent(in) :: nlev ! number of levels

    ! Band containing this g-point, and the band counts of the cloud
    ! optical property arrays
    integer, intent(in) :: i_band, n_bands, n_bands_s

    ! Incoming downwelling solar radiation at top-of-atmosphere (W m-2)
    real(jprb), intent(in),  dimension(ng)         :: incoming_toa

    ! Surface albedo to diffuse and direct radiation
    real(jprb), intent(in),  dimension(ng)         :: albedo_surf_diffuse, &
         &                                              albedo_surf_direct

    ! Cosine of the solar zenith angle
    real(jprb), intent(in)                           :: cos_sza

    ! Gas and aerosol optical depth, single scattering albedo and
    ! asymmetry factor, and the McICA cloud optical depth scaling
    real(jprb), intent(in),  dimension(ng, nlev)   :: od, ssa, asymmetry, od_scaling

    ! In-cloud optical properties, by band
    real(jprb), intent(in),  dimension(n_bands, nlev)   :: od_cloud
    real(jprb), intent(in),  dimension(n_bands_s, nlev) :: ssa_cloud, g_cloud

    ! Cloud fraction of each layer, and the threshold above which a
    ! layer is treated as cloudy
    real(jprb), intent(in),  dimension(nlev)       :: frac
    real(jprb), intent(in)                           :: cloud_fraction_threshold

    ! Resulting fluxes (W m-2) at half-levels: diffuse upwelling,
    ! diffuse downwelling and direct downwelling
    real(jprb), intent(out), dimension(ng, nlev+1) :: flux_up, flux_dn_diffuse, &
         &                                              flux_dn_direct

    ! Albedo of the atmosphere below each half-level
    real(jprb), intent(out), dimension(ng, nlev+1) :: albedo

    ! Upwelling radiation at each half-level due to scattering of the
    ! direct beam below that half-level (W m-2)
    real(jprb), intent(out), dimension(ng, nlev+1) :: source

    ! Layer two-stream properties, recomputed per level into registers
    real(jprb) :: reflectance, transmittance
    real(jprb) :: ref_dir, trans_dir_diff, trans_dir_dir

    ! Combined gas+aerosol+cloud optical properties of the layer
    real(jprb) :: od_cloud_new, od_total, ssa_total, g_total, scat_od

    ! Reciprocal of the denominator of Lacis and Hansen (1974) Eq 33
    real(jprb) :: inv_denominator

    ! Loop index for model level
    integer :: jlev

    ! Compute profile of direct (unscattered) solar fluxes at each
    ! half-level by working down through the atmosphere
    flux_dn_direct(jg,1) = incoming_toa(jg)
    do jlev = 1,nlev
      if (frac(jlev) >= cloud_fraction_threshold) then
         od_total = od(jg,jlev) + od_scaling(jg,jlev) * od_cloud(i_band,jlev)
      else
         od_total = od(jg,jlev)
      end if
      trans_dir_dir = max(-max(od_total * (1.0_jprb/cos_sza),0.0_jprb),-1000.0_jprb)
      trans_dir_dir = exp(trans_dir_dir)
      flux_dn_direct(jg,jlev+1) = flux_dn_direct(jg,jlev)*trans_dir_dir
    end do

    albedo(jg,nlev+1) = albedo_surf_diffuse(jg)

    ! At the surface, the direct solar beam is reflected back into the
    ! diffuse stream
    source(jg,nlev+1) = albedo_surf_direct(jg) * flux_dn_direct(jg,nlev+1) * cos_sza

    ! Work back up through the atmosphere and compute the albedo of
    ! the entire earth/atmosphere system below that half-level, and
    ! also the "source", which is the upwelling flux due to direct
    ! radiation that is scattered below that level
    do jlev = nlev,1,-1
      if (frac(jlev) >= cloud_fraction_threshold) then
         od_cloud_new = od_scaling(jg,jlev) * od_cloud(i_band,jlev)
         od_total  = od(jg,jlev) + od_cloud_new
         ssa_total = 0.0_jprb
         g_total   = 0.0_jprb
         ! In single precision we need to protect against the case that
         ! od_total > 0.0 and ssa_total > 0.0 but od_total*ssa_total == 0
         ! due to underflow
         if (od_total > 0.0_jprb) then
            scat_od = ssa(jg,jlev)*od(jg,jlev) + ssa_cloud(i_band,jlev)*od_cloud_new
            ssa_total = scat_od / od_total
            if (scat_od > 0.0_jprb) then
               g_total = (asymmetry(jg,jlev)*ssa(jg,jlev)*od(jg,jlev) &
                    &     + g_cloud(i_band,jlev)*ssa_cloud(i_band,jlev)*od_cloud_new) &
                    &     / scat_od
            end if
         end if
      else
         od_total  = od(jg,jlev)
         ssa_total = ssa(jg,jlev)
         g_total   = asymmetry(jg,jlev)
      end if
      call calc_ref_trans_sw_scalar_omp(cos_sza, od_total, ssa_total, g_total, &
           &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir)

      ! Lacis and Hansen (1974) Eq 33, Shonk & Hogan (2008) Eq 10:
       inv_denominator = 1.0_jprb / (1.0_jprb-albedo(jg,jlev+1)*reflectance)
       ! Shonk & Hogan (2008) Eq 9, Petty (2006) Eq 13.81:
       albedo(jg,jlev) = reflectance + transmittance * transmittance &
            &                                     * albedo(jg,jlev+1) * inv_denominator
       ! Shonk & Hogan (2008) Eq 11:
       source(jg,jlev) = ref_dir*flux_dn_direct(jg,jlev) &
            &  + transmittance*(source(jg,jlev+1) &
            &        + albedo(jg,jlev+1)*trans_dir_diff*flux_dn_direct(jg,jlev)) &
            &  * inv_denominator
    end do

    ! At top-of-atmosphere there is no diffuse downwelling radiation
    flux_dn_diffuse(jg,1) = 0.0_jprb

    ! At top-of-atmosphere, all upwelling radiation is due to
    ! scattering by the direct beam below that level
    flux_up(jg,1) = source(jg,1)

    ! Work back down through the atmosphere computing the fluxes at
    ! each half-level
    do jlev = 1,nlev
      if (frac(jlev) >= cloud_fraction_threshold) then
         od_cloud_new = od_scaling(jg,jlev) * od_cloud(i_band,jlev)
         od_total  = od(jg,jlev) + od_cloud_new
         ssa_total = 0.0_jprb
         g_total   = 0.0_jprb
         if (od_total > 0.0_jprb) then
            scat_od = ssa(jg,jlev)*od(jg,jlev) + ssa_cloud(i_band,jlev)*od_cloud_new
            ssa_total = scat_od / od_total
            if (scat_od > 0.0_jprb) then
               g_total = (asymmetry(jg,jlev)*ssa(jg,jlev)*od(jg,jlev) &
                    &     + g_cloud(i_band,jlev)*ssa_cloud(i_band,jlev)*od_cloud_new) &
                    &     / scat_od
            end if
         end if
      else
         od_total  = od(jg,jlev)
         ssa_total = ssa(jg,jlev)
         g_total   = asymmetry(jg,jlev)
      end if
      call calc_ref_trans_sw_scalar_omp(cos_sza, od_total, ssa_total, g_total, &
           &  reflectance, transmittance, ref_dir, trans_dir_diff, trans_dir_dir)

       ! The denominator is recomputed here rather than carried from the
       ! sweep above: albedo(jg,jlev+1) is still intact at this point,
       ! being consumed by the flux_up expression below.
       inv_denominator = 1.0_jprb / (1.0_jprb-albedo(jg,jlev+1)*reflectance)
       ! Shonk & Hogan (2008) Eq 14 (after simplification):
       flux_dn_diffuse(jg,jlev+1) &
            &  = (transmittance*flux_dn_diffuse(jg,jlev) &
            &     + reflectance*source(jg,jlev+1) &
            &     + trans_dir_diff*flux_dn_direct(jg,jlev)) * inv_denominator
       ! Shonk & Hogan (2008) Eq 12:
       flux_up(jg,jlev+1) = albedo(jg,jlev+1)*flux_dn_diffuse(jg,jlev+1) &
            &            + source(jg,jlev+1)
       flux_dn_direct(jg,jlev) = flux_dn_direct(jg,jlev)*cos_sza
    end do

    flux_dn_direct(jg,nlev+1) = flux_dn_direct(jg,nlev+1)*cos_sza

  end subroutine adding_ica_sw_recompute_total_omp

end module radiation_adding_ica_sw
