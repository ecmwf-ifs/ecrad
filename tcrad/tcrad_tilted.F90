! tcrad_tilted.F90 - Modify overlap parameter to represent 3D effects by tilting the cloud field
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

module tcrad_tilted

contains

  !---------------------------------------------------------------------
  ! Calculate a modified overlap parameter that approximately
  ! represents 3D effects in the longwave by making the cloud overlap
  ! more random
  subroutine calc_tilted_overlap(nlev, mu, cloud_fraction, overlap_param, &
             &  layer_thickness, inv_cloud_separation_scale, &
             &  overlap_param_tilted)

    use parkind1, only : jprb, jpim
    use yomhook,  only : lhook, dr_hook

    implicit none

    real(jprb), parameter :: PI_OVER_FOUR = acos(-1.0_jprb) * 0.25_jprb
    real(jprb), parameter :: ONE_OVER_PI  = 1.0_jprb / acos(-1.0_jprb)

    ! Inputs

    ! Number of levels
    integer(jpim), intent(in) :: nlev
    
    ! Cosine of zenith angle
    real(jprb),    intent(in) :: mu

    ! Cloud fraction and overlap-parameter profiles
    real(jprb),    intent(in) :: cloud_fraction(nlev)
    real(jprb),    intent(in) :: overlap_param(nlev-1)

    ! Thickness of layers (m) and inverse of cloud separation scale (m-1)
    real(jprb),    intent(in) :: layer_thickness(nlev)
    real(jprb),    intent(in) :: inv_cloud_separation_scale(nlev)

    ! Outputs

    ! Overlap parameter profile with tilting effect to represent 3D effects
    real(jprb),    intent(out):: overlap_param_tilted(nlev-1)

    ! Local variables

    ! Overlap parameter between the base/top and the mid-point of the
    ! cloud in a layer due to tilting
    real(jprb) :: overlap_param_half_layer(nlev)

    ! Inverse of cloud scale (rather than cloud separation scale)
    real(jprb) :: inv_cloud_scale(nlev)

    ! Tangent of zenith angle
    real(jprb) :: tan_za

    real(jprb), dimension(nlev) :: new_fraction, effective_fraction, norm_perim
    real(jprb) :: effective_pair_cover, pair_cover_max, pair_cover_ran, overlap_param_target

    integer(jpim) :: iotherlev, ioverlap

    integer(jpim) :: jlev

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('tcrad:calc_overlap_tilted',0,hook_handle)

    tan_za = sqrt(1.0_jprb - mu*mu)/mu

    ! Eq. 3 of Fielding et al. (QJRMS 2020)
    inv_cloud_scale = inv_cloud_separation_scale &
         &  / (max(1.0e-6_jprb, sqrt(cloud_fraction*(1.0_jprb-cloud_fraction))))

#ifdef WRONG

    overlap_param_half_layer = exp(-(0.5_jprb*PI_OVER_FOUR*tan_za) * layer_thickness * inv_cloud_scale)

    overlap_param_tilted = overlap_param * overlap_param_half_layer(1:nlev-1) * overlap_param_half_layer(2:nlev)

#else

    norm_perim = 4.0_jprb * inv_cloud_scale * cloud_fraction*(1.0_jprb-cloud_fraction)

    overlap_param_tilted = -1.0_jprb

    new_fraction = min(1.0_jprb, (ONE_OVER_PI * tan_za) * layer_thickness * norm_perim)
    effective_fraction = cloud_fraction + new_fraction - cloud_fraction*new_fraction

    do jlev = 1,nlev

      if (cloud_fraction(jlev) > 0.0_jprb .and. cloud_fraction(jlev) < 1.0_jprb) then

      if (jlev == 1) then
        iotherlev = 2
        ioverlap  = 1
      else if (jlev == nlev) then
        iotherlev = nlev-1
        ioverlap  = nlev-1
      else if (cloud_fraction(jlev+1) > cloud_fraction(jlev-1)) then
        iotherlev = jlev+1
        ioverlap  = jlev
      else
        iotherlev = jlev-1
        ioverlap  = jlev-1
      end if
      
      if (cloud_fraction(iotherlev) > 0.0_jprb .and. cloud_fraction(iotherlev) < 1.0_jprb) then

      pair_cover_max = max(cloud_fraction(jlev),cloud_fraction(iotherlev))
      pair_cover_ran = cloud_fraction(jlev) + cloud_fraction(iotherlev) &
           &          -cloud_fraction(jlev) * cloud_fraction(iotherlev)

      effective_pair_cover = overlap_param(ioverlap)*max(effective_fraction(jlev),cloud_fraction(iotherlev)) &
           &  + (1.0_jprb-overlap_param(ioverlap)) * ( effective_fraction(jlev) + cloud_fraction(iotherlev) &
           &                                          -effective_fraction(jlev) * cloud_fraction(iotherlev))

      overlap_param_target = max(0.0_jprb, &
           &  (pair_cover_ran - effective_pair_cover) / (pair_cover_ran - pair_cover_max))

      if (overlap_param_tilted(ioverlap) < 0.0_jprb) then
        overlap_param_tilted(ioverlap) = overlap_param_target
      else
        overlap_param_tilted(ioverlap) = overlap_param_tilted(ioverlap) &
             &   * overlap_param_target / overlap_param(ioverlap)
      end if

    end if

    end if

    end do

    where (overlap_param_tilted < 0.0_jprb)
      overlap_param_tilted = overlap_param
    end where

#endif

!    overlap_param_tilted = overlap_param

    if (lhook) call dr_hook('tcrad:calc_overlap_tilted',1,hook_handle)

  end subroutine calc_tilted_overlap

end module tcrad_tilted
