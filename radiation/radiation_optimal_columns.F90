! radiation_optimal_columns.F90 - Generate optimal sub-columns for a shortwave radiance solver
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
  ! Generate nsub cloudy sub-columns in ng spectral intervals using
  ! the optimal-column technique in which the standard deviation of
  ! column optical depth and then Generalized Gauss-Laguerre
  ! Quadrature is used to select quadrature points assuming the
  ! distribution of column optical depth (in the cloud-covered part of
  ! the gridbox) follows a gamma distribution.
  subroutine optimal_columns(ng, nsub, nlev, frac_threshold, frac, overlap_param, &
       &  fractional_std, od_in, weight, od_out, total_cloud_cover, cloudy_fsd_od)

    use parkind1, only           : jprb, jpim
    use yomhook,  only           : lhook, dr_hook
    use radiation_cloud_cover, only : cum_cloud_cover_exp_ran
    !use radiation_gen_gauss_laguerre, only : calc_gen_gauss_laguerre
    !use radiation_gauss_lognormal, only : calc_gauss_lognormal
    use radiation_gauss_gamma, only : calc_gauss_gamma

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
    ! calling routine simulates a clear column, if required. Since the
    ! optical depth fractional standard deviation may be different in
    ! each spectral interval, the weights may vary with spectral
    ! interval.
    real(jprb), intent(out) :: weight(ng,nsub)

    ! Cloud optical depth in each layer, subcolumn and spectral
    ! interval
    real(jprb), intent(out) :: od_out(ng,nlev,nsub)

    ! Total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(out) :: total_cloud_cover

    ! Fractional standard deviation of the total-column optical depth
    ! in the cloudy part of the gridbox
    real(jprb), intent(out), optional :: cloudy_fsd_od(ng)

    ! Local variables

    ! Scaling factor to be multiplied by the optical depth of each
    ! layer and each spectral interval to spread the optical depth in
    ! the cloudy part of the layer to the cloud cover
    real(jprb) :: od_scaling(nlev)

    ! Nodes corresponding to "weight" above
    real(jprb) :: node(ng,nsub)

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb) :: cum_cloud_cover(nlev)

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb) :: pair_cloud_cover(nlev-1)

    ! Accumulated gridbox-mean optical depth from top-of-atmosphere
    real(jprb) :: acc_gridbox_od(ng)

    ! Accumulated covariance of the optical depth variations in a
    ! layer with the optical depth of all the layers above, normalized
    ! by the standard deviation of the optical depth of the layer (all
    ! gridbox total variables)
    real(jprb) :: acc_z(ng)

    ! Accumulated variance of optical depth
    real(jprb) :: acc_var_od(ng)

    ! Variance and standard-deviation of layer optical depth (gridbox
    ! total)
    real(jprb) :: gridbox_var_od(ng,nlev), gridbox_std_od(ng,nlev)

    ! Correlation coefficient of gridbox-total optical depth between
    ! adjacent layers
    real(jprb) :: gridbox_corr_coeff(nlev-1)
    
    ! Properties of cloudy part of gridbox: mean column optical depth,
    ! standard deviation of optical depth and fractional standard
    ! deviation of optical depth
    real(jprb) :: cloudy_od(ng), cloudy_std_od(ng), cloudy_fsd_od_local(ng)

    ! Loop index for level and subcolumn
    integer(jpim) :: jlev, jsub

    ! Considering two adjacent layers, fracUL is the fraction of the
    ! gridbox covered in the case that the upper layer U is cloudy (1)
    ! or clear (0) and equivalently for the lower layer L
    real(jprb) :: frac00, frac01, frac10, frac11

    ! Correlation coefficient of in-cloud heterogeneities
    real(jprb) :: cloud_corr_coeff

    logical :: is_first_layer

    real(jprb) :: hook_handle

    if (lhook) call dr_hook('radiation_optimal_columns:optimal_columns',0,hook_handle)

    ! Work out the cloud cover
    call cum_cloud_cover_exp_ran(nlev, frac, overlap_param, &
         &   cum_cloud_cover, pair_cloud_cover)
    total_cloud_cover = cum_cloud_cover(nlev)

    if (present(cloudy_fsd_od)) then
      cloudy_fsd_od = 0.0_jprb
    end if

    if (total_cloud_cover < frac_threshold) then

      ! No cloud
      weight     = 0.0_jprb
      od_out     = 0.0_jprb

    else if (nsub == 1) then

      ! Simple case of one cloudy column: dilute the in-cloud mean
      ! optical depth into the total cloud cover fraction
      weight(:,1) = total_cloud_cover
      od_scaling = frac / total_cloud_cover
      do jlev = 1,nlev
        od_out(:,jlev,1) = od_in(:,jlev) * od_scaling(jlev)
      end do

    else

      ! Variance of layer optical depth across the whole gridbox
      gridbox_var_od = 0.0_jprb
      gridbox_std_od = 0.0_jprb
      do jlev = 1,nlev
        if (frac(jlev) > frac_threshold) then
          gridbox_var_od(:,jlev) = frac(jlev) * (fractional_std(jlev)**2 &
               &  + 1.0_jprb - frac(jlev)) * od_in(:,jlev)**2
          gridbox_std_od(:,jlev) = sqrt(gridbox_var_od(:,jlev))
        end if
      end do

      ! Correlation coefficient
      gridbox_corr_coeff = 0.0_jprb
      do jlev = 1,nlev-1
        if (frac(jlev) > frac_threshold .and. frac(jlev+1) > frac_threshold) then
          ! Cloudy in both layers
          frac11 = frac(jlev) + frac(jlev+1) - pair_cloud_cover(jlev)
          ! Clear in upper layer, cloudy in lower layer
          frac01 = pair_cloud_cover(jlev) - frac(jlev)
          ! Cloudy in upper layer, clear in lower layer
          frac10 = pair_cloud_cover(jlev) - frac(jlev+1)
          ! Clear in both layers
          frac00 = 1.0_jprb - pair_cloud_cover(jlev)

          ! Assume the cloud optical-depth heterogeneities decorrelate
          ! twice as rapidly as the cloud boundaries
          cloud_corr_coeff = overlap_param(jlev) * overlap_param(jlev)

          gridbox_corr_coeff(jlev) = od_in(1,jlev)*od_in(1,jlev+1) *&
               &   ( frac00 * frac(jlev)*frac(jlev+1) &
               &    +frac01 * (frac(jlev)*frac(jlev+1) - frac(jlev)) &
               &    +frac10 * (frac(jlev)*frac(jlev+1) - frac(jlev+1)) &
               &    +frac11 * (1.0_jprb + frac(jlev)*frac(jlev+1)-frac(jlev)-frac(jlev+1) &
               &               + fractional_std(jlev)*fractional_std(jlev+1)*cloud_corr_coeff) ) &
               &  / (gridbox_std_od(1,jlev)*gridbox_std_od(1,jlev+1))
        end if
      end do 

      ! Main optimal subcolumn algorithm

      acc_gridbox_od = 0.0_jprb
      acc_z = 0.0_jprb
      acc_var_od = 0.0_jprb

      is_first_layer = .true.

      do jlev = 1,nlev
        if (frac(jlev) > frac_threshold) then
          ! Only consider cloudy layers

          ! Accumulate the gridbox-mean optical depth
          acc_gridbox_od = acc_gridbox_od + od_in(:,jlev) * frac(jlev)

          if (.not. is_first_layer) then
            acc_z = (acc_z + gridbox_std_od(:,jlev-1)) * gridbox_corr_coeff(jlev-1)
          end if
            
          if (is_first_layer) then
            acc_var_od = gridbox_var_od(:,jlev)
            is_first_layer = .false.
          else
            acc_var_od = acc_var_od + gridbox_var_od(:,jlev) &
                 &     + 2.0_jprb * acc_z * gridbox_std_od(:,jlev)
          end if

        end if
      end do

      ! From the final accumulated values, compute the mean column
      ! optical depth of the cloud-covered part of the gridbox...
      cloudy_od = acc_gridbox_od / total_cloud_cover
      ! ...and the standard deviation of optical depth of the
      ! cloud-covered part of the gridbox...
      cloudy_std_od = sqrt(max(0.0_jprb, (acc_var_od - acc_gridbox_od*acc_gridbox_od &
           &        *(1.0_jprb-total_cloud_cover)/total_cloud_cover)/total_cloud_cover))
      ! ...from which we calculate the fractional standard deviation.
      cloudy_fsd_od_local = cloudy_std_od / cloudy_od
      
      if (present(cloudy_fsd_od)) then
        cloudy_fsd_od = cloudy_fsd_od_local
      end if

      ! Generalized Gauss-Laguerre Quadrature then yields the weights
      ! and nodes to optimally integrate over the cloudy part of the
      ! gridbox
      !call calc_gen_gauss_laguerre(ng, nsub, cloudy_fsd_od_local, weight, node)
      !call calc_gauss_lognormal(ng, nsub, cloudy_fsd_od_local, weight, node)
      call calc_gauss_gamma(ng, nsub, cloudy_fsd_od_local, weight, node)

      ! Scale the results into the output "weight" and "od_out"
      ! variables
      od_scaling = frac / total_cloud_cover
      weight = weight * total_cloud_cover

      do jsub = 1,nsub
        do jlev = 1,nlev
          od_out(:,jlev,jsub) = node(:,jsub) * od_in(:,jlev) * od_scaling(jlev)
        end do
      end do
      
    end if

    if (lhook) call dr_hook('radiation_optimal_columns:optimal_columns',1,hook_handle)

  end subroutine optimal_columns

end module radiation_optimal_columns
