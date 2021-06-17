! radiation_optimal_columns.F90 - Generate optimal sub-columns for a shortwave solver
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

module radiation_optimal_columns

  public

contains

  !---------------------------------------------------------------------
  ! Generate ncol sub-columns 
  subroutine optimal_columns(ng, nsub, nlev, frac_threshold, frac, overlap_param, &
       &  fractional_std, od_in, weight, od_out, total_cloud_cover)

    use parkind1, only           : jprb
    use yomhook,  only           : lhook, dr_hook
    use radiation_cloud_cover, only : cum_cloud_cover_exp_ran

    implicit none

    ! Inputs
    integer, intent(in) :: ng    ! Number of spectral intervals
    integer, intent(in) :: nsub  ! Number of cloudy sub-columns to generate.
    integer, intent(in) :: nlev  ! Number of model levels


    ! Only cloud fractions above this threshold are considered to be
    ! clouds
    real(jprb), intent(in)  :: frac_threshold

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! Fractional standard deviation at each layer
    real(jprb), intent(in)  :: fractional_std(nlev)

    ! Mean cloud optical depth of cloudy part of gridbox in each layer
    ! and each spectral interval
    real(jprb), intent(in) :: od_in(ng, nlev)

    ! Outputs

    ! Weight of each sub-column, or equivalently the fractional area
    ! of the total gridbox occupied by that sub-column. Since this
    ! code does not simulate a clear column, the weights sum to the
    ! total cloud cover, which might not be 1. It is assumed that the
    ! calling routine simulates a clear column, if required.
    real(jprb), intent(out) :: weight(nsub)

    ! Cloud optical depth in each layer, subcolumn and spectral
    ! interval
    real(jprb), intent(out) :: od_out(ng,nlev,nsub)

    ! Total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(out) :: total_cloud_cover

    ! Local variables

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb) :: cum_cloud_cover(nlev)

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), dimension(nlev-1) :: pair_cloud_cover

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_optimal_columns:optimal_columns',0,hook_handle)

    ! Work out the cloud cover
    call cum_cloud_cover_exp_ran(nlev, frac, overlap_param, &
         &   cum_cloud_cover, pair_cloud_cover)
    total_cloud_cover = cum_cloud_cover(nlev)

    if (total_cloud_cover < frac_threshold) then

      ! No cloud
      weight(:) = 0.0_jprb
      od_scaling(:,:) = 0.0_jprb

    else if (nsub == 1) then

      ! Simple case of one cloudy column: dilute the in-cloud mean
      ! optical depth into the total cloud cover fraction
      weight(1) = total_cloud_cover
      od_scaling(1,:) = frac / total_cloud_cover

    else

      ! Main optimal subcolumn algorithm

      acc_gridbox_od = 0.0_jprb
      acc_z = 0.0_jprb
      acc_var_od = 0.0_jprb

      do jlev = 1,nlev
        if (frac(jlev) > frac_threshold) then
          ! Only consider cloudy layers

          ! Accumulate the gridbox-mean optical depth
          acc_gridbox_od = acc_total_od + od_in(:,jlev) * frac(jlev)

          ! Variance of layer optical depth across the whole gridbox
          gridbox_var_od = frac(jlev) * (fractional_std(jlev)**2 + 1.0_jprb - frac(jlev)) &
               &   * od_in(:,jlev)**2
          gridbox_std_od = sqrt(gridbox_var_od)

          acc_z = (acc_z + gridbox_std_od) * overlap_param(jlev)
          acc_var_od = acc_var_od + gridbox_var_od + 2.0_jprb * acc_z * gridbox_std_od
        end if
      end do

      ! Trivial case of one homogeneous column: we need an
      ! optical-depth scaling that dilutes the in-cloud mean optical
      ! depth into a gridbox mean
      weight(1) = 1.0_jprb
      od_scaling(1,:) = frac
      if (any(frac >= frac_threshold)) then
        total_cloud_cover = 1.0_jprb
      else
        total_cloud_cover = 0.0_jprb
      end if
    else
      ! We have more than one column; work out the cloud cover
      call cum_cloud_cover_exp_ran(nlev, frac, overlap_param, &
           &   cum_cloud_cover, pair_cloud_cover)
      total_cloud_cover = cum_cloud_cover(nlev)
      ! The first column is clear sky
      weight(1) = 1.0_jprb - total_cloud_cover
      od_scaling(1,:) = 0.0_jprb

      if (total_cloud_cover > frac_threshold) then
        if (nsub == 2) then
          


        end if
      else
        ! Clear-sky profile: set other columns to zero
        weight(1) = 1.0_jprb
        weight(2:nsub) = 0.0_jprb
        od_scaling(2:nsub,:) = 0.0_jprb
      end if

    end if

    if (lhook) call dr_hook('radiation_optimal_columns:optimal_columns',1,hook_handle)

  end subroutine optimal_columns

end module radiation_optimal_columns
