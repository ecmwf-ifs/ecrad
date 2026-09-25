! radiation_cloud_generator_acc.F90 - Generate water-content or optical-depth scalings for McICA
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
! Generate clouds for McICA using a method modified from Raisanen et
! al. (2002)
! This is a copy of the original cloud_generator, that is better suited for OpenACC
!
! Modifications
!   2018-02-22  R. Hogan  Call masked version of PDF sampler for speed
!   2020-03-31  R. Hogan  More vectorizable version of Exp-Ran
!   2022-11-07  D. Hupp adaptation for ACC

module radiation_cloud_generator_acc

  implicit none

  !$omp declare target(cloud_generator_omp)

  public

contains

  !---------------------------------------------------------------------
  ! Generate scaling factors for the cloud optical depth to represent
  ! cloud overlap, the overlap of internal cloud inhomogeneities and
  ! the fractional standard deviation of these inhomogeneities, for
  ! use in a Monte Carlo Independent Column Approximation radiation
  ! scheme. All returned profiles contain cloud, and the total cloud
  ! cover is also returned, so the calling function can then do a
  ! weighted average of clear and cloudy skies; this is a way to
  ! reduce the Monte Carlo noise in profiles with low cloud cover.
  subroutine cloud_generator_acc(ng, nlev, &
    &  iseed, frac_threshold, &
    &  frac, overlap_param, decorrelation_scaling, &
    &  fractional_std, &
    &  sample_ncdf, sample_nfsd, sample_fsd1, &
    &  sample_inv_fsd_interval, sample_val, &
    &  od_scaling, total_cloud_cover, &
    &  ibegin, iend, &
    &  cum_cloud_cover, &
    &  pair_cloud_cover)

    use parkind1,                 only : jprb, jpib
    use radiation_random_numbers, only : initialize_acc, uniform_distribution_acc
    use radiation_random_numbers, only : initialize_acc, uniform_distribution_acc, IMinstdA0, IMinstdA, IMinstdM

    implicit none

    ! Inputs
    integer, intent(in)     :: ng    ! number of g points
    integer, intent(in)     :: nlev  ! number of model levels
    integer, intent(in)     :: iseed ! seed for random number generator

    ! Only cloud fractions above this threshold are considered to be
    ! clouds
    real(jprb), intent(in)  :: frac_threshold

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! Overlap parameter for internal inhomogeneities
    real(jprb), intent(in)  :: decorrelation_scaling

    ! Fractional standard deviation at each layer
    real(jprb), intent(in)  :: fractional_std(nlev)

    ! Object for sampling from a lognormal or gamma distribution
    integer, intent(in)  :: sample_ncdf, sample_nfsd
    real(jprb), intent(in)  :: sample_fsd1, sample_inv_fsd_interval
    real(jprb), intent(in), dimension(:,:)  :: sample_val

    ! Outputs

    ! Cloud optical depth scaling factor, with 0 indicating clear sky
    real(jprb), intent(out) :: od_scaling(ng,nlev)

    ! Total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(in) :: total_cloud_cover

    ! Local variables

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb), intent(in) :: cum_cloud_cover(nlev)

    ! First and last cloudy layers
    integer, intent(in) :: ibegin, iend

    ! Scaled random number for finding cloud
    real(jprb) :: trigger

    ! Uniform deviates between 0 and 1
    real(jprb) :: rand_top

    ! Overlap parameter of inhomogeneities
    real(jprb) :: overlap_param_inhom

    real(jprb) :: rand_cloud, rand_inhom1, rand_inhom2

    integer :: itrigger

    ! Loop index for model level and g-point
    integer :: jlev, jg

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), intent(inout), dimension(nlev-1) :: pair_cloud_cover
    real(jprb) :: overhang

    integer          :: jcloud
    ! Number of contiguous cloudy layers for which to compute optical
    ! depth scaling
    integer :: n_layers_to_scale

    ! Is it time to fill the od_scaling variable?
    logical :: do_fill_od_scaling

    ! variables from manual inlining of sample_vec_from_pdf
    ! Index to look-up table
    integer    :: ifsd, icdf
    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! local variable for the acc rng
    integer(kind=jpib) :: istate

    !$ACC ROUTINE WORKER

    if (total_cloud_cover >= frac_threshold) then
      ! Loop over ng columns (this loop should be paralized as soon as acc RNG is available)
      !$ACC LOOP WORKER VECTOR PRIVATE(istate, rand_top, trigger, itrigger, n_layers_to_scale)
      do jg = 1,ng

        !$ACC LOOP SEQ 
        do jlev = 1,nlev
          od_scaling(jg,jlev) = 0.0_jprb
        end do

        ! begin manuel inline of istate = initialize_acc(iseed, jg)
        istate = REAL(ABS(iseed),jprb)
        ! Use a modified (and vectorized) C++ minstd_rand0 algorithm to populate the state
        istate = nint(mod( istate*jg*(1._jprb-0.05_jprb*jg+0.005_jprb*jg**2)*IMinstdA0, IMinstdM))

        ! One warmup of the C++ minstd_rand algorithm
        istate = mod(IMinstdA * istate, IMinstdM)
        ! end manuel inline of istate = initialize_acc(iseed, jg)
        rand_top = uniform_distribution_acc(istate)

        ! Find the cloud top height corresponding to the current
        ! random number, and store in itrigger
        trigger = rand_top * total_cloud_cover
        itrigger = iend
        !$ACC LOOP SEQ
        do jlev = ibegin,iend
          if (trigger <= cum_cloud_cover(jlev)) then
            itrigger = min(jlev, itrigger)
          end if
        end do

        ! manual inline: call generate_column_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! So far our vertically contiguous cloud contains only one layer
        n_layers_to_scale = 1

        ! Locate the clouds below this layer: first generate some more
        ! random numbers
        ! Loop from the layer below the local cloud top down to the
        ! bottom-most cloudy layer
        !$ACC LOOP SEQ PRIVATE(do_fill_od_scaling, rand_cloud, rand_inhom1, overhang)
        do jlev = itrigger+1,iend+1
          do_fill_od_scaling = .false.
          if (jlev <= iend) then
            rand_cloud = uniform_distribution_acc(istate)
            if (n_layers_to_scale > 0) then
              ! There is a cloud above, in which case the probability
              ! of cloud in the layer below is as follows
              if (rand_cloud*frac(jlev-1) &
                  &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
                ! Add another cloudy layer
                n_layers_to_scale = n_layers_to_scale + 1
              else
                ! Reached the end of a contiguous set of cloudy layers and
                ! will compute the optical depth scaling immediately.
                do_fill_od_scaling = .true.
              end if
            else
              overhang = cum_cloud_cover(jlev)-cum_cloud_cover(jlev-1)
              ! There is clear-sky above, in which case the
              ! probability of cloud in the layer below is as follows
              if (rand_cloud*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
                  &  < pair_cloud_cover(jlev-1) - overhang - frac(jlev-1)) then
                ! A new cloud top
                n_layers_to_scale = 1
              end if
            end if
          else
            ! We are at the bottom of the cloudy layers in the model,
            ! so in a moment need to populate the od_scaling array
            do_fill_od_scaling = .true.
          end if ! (jlev <= iend)

          if (do_fill_od_scaling) then

            rand_inhom1 = uniform_distribution_acc(istate)
            ! Loop through the sequence of cloudy layers
            !$ACC LOOP SEQ PRIVATE(rand_inhom2, overlap_param_inhom, wcdf, icdf, &
            !$ACC   wfsd, ifsd)
            do jcloud = max(2,jlev-n_layers_to_scale),jlev-1

              rand_inhom2 = uniform_distribution_acc(istate)

              ! Set overlap parameter of inhomogeneities
              overlap_param_inhom = overlap_param(jcloud-1)
              if ( ibegin<=jcloud-1 .and. jcloud-1< iend .and. &
                & overlap_param(jcloud-1) > 0.0_jprb) then
                overlap_param_inhom = &
                  & overlap_param(jcloud-1)**(1.0_jprb/decorrelation_scaling)
              end if

              ! Use second random number, and inhomogeneity overlap
              ! parameter, to decide whether the first random number
              ! should be repeated (corresponding to maximum overlap)
              ! or not (corresponding to random overlap)
              if ( jcloud > jlev-n_layers_to_scale .and. &
                & rand_inhom2 >= overlap_param_inhom) then
                rand_inhom1 = uniform_distribution_acc(istate)
              end if
              ! manual inline >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! Bilinear interpolation with bounds
              wcdf = rand_inhom1 * (sample_ncdf-1) + 1.0_jprb
              icdf = max(1, min(int(wcdf), sample_ncdf-1))
              wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

              wfsd = (fractional_std(jcloud)-sample_fsd1) * sample_inv_fsd_interval + 1.0_jprb
              ifsd = max(1, min(int(wfsd), sample_nfsd-1))
              wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

              od_scaling(jg,jcloud) =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * sample_val(icdf  ,ifsd)   &
                  & + (1.0_jprb-wcdf)*          wfsd  * sample_val(icdf  ,ifsd+1) &
                  & +           wcdf *(1.0_jprb-wfsd) * sample_val(icdf+1,ifsd)   &
                  & +           wcdf *          wfsd  * sample_val(icdf+1,ifsd+1)
              ! manual inline <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end do

            n_layers_to_scale = 0

          end if ! (do_fill_od_scaling)
        end do ! jlev = itrigger+1,iend+1
      end do ! jg = 1,ng
    end if ! (total_cloud_cover < frac_threshold)

  end subroutine cloud_generator_acc


    !---------------------------------------------------------------------
  ! Generate scaling factors for the cloud optical depth to represent
  ! cloud overlap, the overlap of internal cloud inhomogeneities and
  ! the fractional standard deviation of these inhomogeneities, for
  ! use in a Monte Carlo Independent Column Approximation radiation
  ! scheme. All returned profiles contain cloud, and the total cloud
  ! cover is also returned, so the calling function can then do a
  ! weighted average of clear and cloudy skies; this is a way to
  ! reduce the Monte Carlo noise in profiles with low cloud cover.
  subroutine cloud_generator_omp(jg, ng, nlev, &
    &  iseed, frac_threshold, &
    &  frac, overlap_param, decorrelation_scaling, &
    &  fractional_std, &
    &  sample_ncdf, sample_nfsd, sample_fsd1, &
    &  sample_inv_fsd_interval, sample_val, &
    &  od_scaling, total_cloud_cover, &
    &  ibegin, iend, &
    &  cum_cloud_cover, &
    &  pair_cloud_cover)

    use parkind1,                 only : jprb, jpib
    use radiation_random_numbers, only : IMinstdA0, IMinstdA, IMinstdM, uniform_distribution_omp
    implicit none

    ! Inputs
    integer, intent(in)     :: jg    ! index into the g direction
    integer, intent(in)     :: ng    ! number of g points
    integer, intent(in)     :: nlev  ! number of model levels
    integer, intent(in)     :: iseed ! seed for random number generator

    ! Only cloud fractions above this threshold are considered to be
    ! clouds
    real(jprb), intent(in)  :: frac_threshold

    ! Cloud fraction on full levels
    real(jprb), intent(in)  :: frac(nlev)

    ! Cloud overlap parameter for interfaces between model layers,
    ! where 0 indicates random overlap and 1 indicates maximum-random
    ! overlap
    real(jprb), intent(in)  :: overlap_param(nlev-1)

    ! Overlap parameter for internal inhomogeneities
    real(jprb), intent(in)  :: decorrelation_scaling

    ! Fractional standard deviation at each layer
    real(jprb), intent(in)  :: fractional_std(nlev)

    ! Object for sampling from a lognormal or gamma distribution
    integer, intent(in)  :: sample_ncdf, sample_nfsd
    real(jprb), intent(in)  :: sample_fsd1, sample_inv_fsd_interval
    real(jprb), intent(in), dimension(:,:)  :: sample_val

    ! Outputs

    ! Cloud optical depth scaling factor, with 0 indicating clear sky
    real(jprb), intent(out) :: od_scaling(ng,nlev)

    ! Total cloud cover using cloud fraction and overlap parameter
    real(jprb), intent(in) :: total_cloud_cover

    ! Local variables

    ! Cumulative cloud cover from TOA to the base of each layer
    real(jprb), intent(in) :: cum_cloud_cover(nlev)

    ! First and last cloudy layers
    integer, intent(in) :: ibegin, iend

    ! Scaled random number for finding cloud
    real(jprb) :: trigger

    ! Uniform deviates between 0 and 1
    real(jprb) :: rand_top

    ! Overlap parameter of inhomogeneities
    real(jprb) :: overlap_param_inhom

    real(jprb) :: rand_cloud, rand_inhom1, rand_inhom2

    integer :: itrigger

    ! Loop index for model level and g-point
    integer :: jlev

    ! Cloud cover of a pair of layers, and amount by which cloud at
    ! next level increases total cloud cover as seen from above
    real(jprb), intent(inout), dimension(nlev-1) :: pair_cloud_cover
    real(jprb) :: overhang

    integer          :: jcloud
    ! Number of contiguous cloudy layers for which to compute optical
    ! depth scaling
    integer :: n_layers_to_scale

    ! Is it time to fill the od_scaling variable?
    logical :: do_fill_od_scaling

    ! variables from manual inlining of sample_vec_from_pdf
    ! Index to look-up table
    integer    :: ifsd, icdf
    ! Weights in bilinear interpolation
    real(jprb) :: wfsd, wcdf

    ! local variable for the acc rng
    integer(kind=jpib) :: istate

    if (total_cloud_cover >= frac_threshold) then
      ! Loop over ng columns (this loop should be paralized as soon as acc RNG is available)
       do jlev = 1,nlev
          od_scaling(jg,jlev) = 0.0_jprb
       end do

       ! begin manuel inline of istate = initialize_acc(iseed, jg)
       istate = REAL(ABS(iseed),jprb)
       ! Use a modified (and vectorized) C++ minstd_rand0 algorithm to populate the state
       istate = nint(mod( istate*jg*(1._jprb-0.05_jprb*jg+0.005_jprb*jg**2)*IMinstdA0, IMinstdM))

       ! One warmup of the C++ minstd_rand algorithm
       istate = mod(IMinstdA * istate, IMinstdM)
       ! end manuel inline of istate = initialize_acc(iseed, jg)
       call uniform_distribution_omp(istate, rand_top)

       ! Find the cloud top height corresponding to the current
       ! random number, and store in itrigger
       trigger = rand_top * total_cloud_cover
       itrigger = iend
       !$ACC LOOP SEQ
       do jlev = ibegin,iend
         if (trigger <= cum_cloud_cover(jlev)) then
           itrigger = min(jlev, itrigger)
         end if
       end do
       ! manual inline: call generate_column_exp_ran >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       ! So far our vertically contiguous cloud contains only one layer
       n_layers_to_scale = 1

       ! Locate the clouds below this layer: first generate some more
       ! random numbers
       ! Loop from the layer below the local cloud top down to the
       ! bottom-most cloudy layer
       do jlev = itrigger+1,iend+1
          do_fill_od_scaling = .false.
          if (jlev <= iend) then
             call uniform_distribution_omp(istate, rand_cloud)
             if (n_layers_to_scale > 0) then
                ! There is a cloud above, in which case the probability
                ! of cloud in the layer below is as follows
                if (rand_cloud*frac(jlev-1) &
                     &  < frac(jlev) + frac(jlev-1) - pair_cloud_cover(jlev-1)) then
                   ! Add another cloudy layer
                   n_layers_to_scale = n_layers_to_scale + 1
                else
                   ! Reached the end of a contiguous set of cloudy layers and
                   ! will compute the optical depth scaling immediately.
                   do_fill_od_scaling = .true.
                end if
             else
                overhang = cum_cloud_cover(jlev)-cum_cloud_cover(jlev-1)
                ! There is clear-sky above, in which case the
                ! probability of cloud in the layer below is as follows
                if (rand_cloud*(cum_cloud_cover(jlev-1) - frac(jlev-1)) &
                     &  < pair_cloud_cover(jlev-1) - overhang - frac(jlev-1)) then
                   ! A new cloud top
                   n_layers_to_scale = 1
                end if
             end if
          else
             ! We are at the bottom of the cloudy layers in the model,
             ! so in a moment need to populate the od_scaling array
             do_fill_od_scaling = .true.
          end if ! (jlev <= iend)

          if (do_fill_od_scaling) then
             call uniform_distribution_omp(istate, rand_inhom1)
             ! Loop through the sequence of cloudy layers
             !$ACC LOOP SEQ PRIVATE(rand_inhom2, overlap_param_inhom, wcdf, icdf, &
             !$ACC   wfsd, ifsd)
             do jcloud = max(2,jlev-n_layers_to_scale),jlev-1

                call uniform_distribution_omp(istate, rand_inhom2)

                ! Set overlap parameter of inhomogeneities
                overlap_param_inhom = overlap_param(jcloud-1)
                if ( ibegin<=jcloud-1 .and. jcloud-1< iend .and. &
                     & overlap_param(jcloud-1) > 0.0_jprb) then
                   overlap_param_inhom = &
                        & overlap_param(jcloud-1)**(1.0_jprb/decorrelation_scaling)
                end if

                ! Use second random number, and inhomogeneity overlap
                ! parameter, to decide whether the first random number
                ! should be repeated (corresponding to maximum overlap)
                ! or not (corresponding to random overlap)
                if ( jcloud > jlev-n_layers_to_scale .and. &
                     & rand_inhom2 >= overlap_param_inhom) then
                  call uniform_distribution_omp(istate, rand_inhom1)
                end if
                ! manual inline >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                ! Bilinear interpolation with bounds
                wcdf = rand_inhom1 * (sample_ncdf-1) + 1.0_jprb
                icdf = max(1, min(int(wcdf), sample_ncdf-1))
                wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

                wfsd = (fractional_std(jcloud)-sample_fsd1) * sample_inv_fsd_interval + 1.0_jprb
                ifsd = max(1, min(int(wfsd), sample_nfsd-1))
                wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

                od_scaling(jg,jcloud) =    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * sample_val(icdf  ,ifsd)   &
                     & + (1.0_jprb-wcdf)*          wfsd  * sample_val(icdf  ,ifsd+1) &
                     & +           wcdf *(1.0_jprb-wfsd) * sample_val(icdf+1,ifsd)   &
                     & +           wcdf *          wfsd  * sample_val(icdf+1,ifsd+1)
                ! manual inline <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             end do

             n_layers_to_scale = 0

          end if ! (do_fill_od_scaling)
       end do ! jlev = itrigger+1,iend+1
    end if ! (total_cloud_cover < frac_threshold)

  end subroutine cloud_generator_omp


#ifdef __NVCOMPILER
  !---------------------------------------------------------------------
  ! As cloud_generator_omp, but the column and g-point loops and the
  ! target region live here rather than in the caller, so the arrays are
  ! passed whole and indexed by jcol internally.
  !
  ! This exists only because nvfortran mis-computes the device base address
  ! of a per-column slice passed to cloud_generator_omp from inside a target
  ! region; see the call sites in the McICA solvers.
  subroutine cloud_generator_block_omp(ncol, istartcol, iendcol, ng, nlev, &
       &  iseed, iseed_offset, frac_threshold, frac, overlap_param, &
       &  decorrelation_scaling, &
       &  fractional_std, sample_ncdf, sample_nfsd, sample_fsd1, &
       &  sample_inv_fsd_interval, sample_val, od_scaling, total_cloud_cover, &
       &  ibegin, iend, cum_cloud_cover, pair_cloud_cover, thread_limit)

    use parkind1,                 only : jprb, jpib
    use radiation_random_numbers, only : IMinstdA0, IMinstdA, IMinstdM, &
         &                               uniform_distribution_omp

    implicit none

    integer,    intent(in) :: ncol, istartcol, iendcol, ng, nlev
    integer,    intent(in) :: iseed(ncol)
    ! Added to iseed so the longwave and shortwave draw different samples
    integer,    intent(in) :: iseed_offset
    real(jprb), intent(in) :: frac_threshold, decorrelation_scaling
    real(jprb), intent(in) :: frac(nlev, istartcol:iendcol)
    real(jprb), intent(in) :: overlap_param(nlev-1, istartcol:iendcol)
    real(jprb), intent(in) :: fractional_std(nlev, istartcol:iendcol)
    integer,    intent(in) :: sample_ncdf, sample_nfsd
    real(jprb), intent(in) :: sample_fsd1, sample_inv_fsd_interval
    real(jprb), intent(in) :: sample_val(sample_ncdf, sample_nfsd)
    real(jprb), intent(inout) :: od_scaling(ng, nlev, istartcol:iendcol)
    real(jprb), intent(in) :: total_cloud_cover(ncol)
    integer,    intent(in) :: ibegin(istartcol:iendcol), iend(istartcol:iendcol)
    real(jprb), intent(in) :: cum_cloud_cover(nlev, istartcol:iendcol)
    real(jprb), intent(in) :: pair_cloud_cover(nlev-1, istartcol:iendcol)
    integer,    intent(in) :: thread_limit

    real(jprb) :: trigger, rand_top, rand_cloud, rand_inhom1, rand_inhom2
    real(jprb) :: overlap_param_inhom, overhang, wfsd, wcdf, cover
    integer    :: jcol, jg, jlev, jcloud, itrigger, n_layers_to_scale
    integer    :: ifsd, icdf, icol_begin, icol_end
    logical    :: do_fill_od_scaling
    integer(kind=jpib) :: istate

    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
    !$OMP& PRIVATE(jcol, jg, jlev, trigger, rand_top, rand_cloud, &
    !$OMP&         rand_inhom1, rand_inhom2, overlap_param_inhom, overhang, &
    !$OMP&         wfsd, wcdf, cover, itrigger, jcloud, n_layers_to_scale, &
    !$OMP&         ifsd, icdf, icol_begin, icol_end, do_fill_od_scaling, istate) &
    !$OMP& FIRSTPRIVATE(istartcol, iendcol, ng, nlev, iseed_offset) &
    !$OMP& NUM_TEAMS(((iendcol-istartcol+1)*ng + 127) / 128) &
    !$OMP& THREAD_LIMIT(thread_limit)
    do jcol = istartcol,iendcol
       do jg = 1, ng
          cover      = total_cloud_cover(jcol)
          icol_begin = ibegin(jcol)
          icol_end   = iend(jcol)

          if (cover >= frac_threshold) then
             do jlev = 1, nlev
                od_scaling(jg,jlev,jcol) = 0.0_jprb
             end do

             istate = REAL(ABS(iseed(jcol) + iseed_offset),jprb)
             istate = nint(mod( istate*jg &
                  &  *(1._jprb-0.05_jprb*jg+0.005_jprb*jg**2)*IMinstdA0, IMinstdM))
             istate = mod(IMinstdA * istate, IMinstdM)
             call uniform_distribution_omp(istate, rand_top)

             trigger  = rand_top * cover
             itrigger = icol_end
             do jlev = icol_begin, icol_end
                if (trigger <= cum_cloud_cover(jlev,jcol)) then
                   itrigger = min(jlev, itrigger)
                end if
             end do

             n_layers_to_scale = 1
             do jlev = itrigger+1, icol_end+1
                do_fill_od_scaling = .false.
                if (jlev <= icol_end) then
                   call uniform_distribution_omp(istate, rand_cloud)
                   if (n_layers_to_scale > 0) then
                      if (rand_cloud*frac(jlev-1,jcol) &
                           &  < frac(jlev,jcol) + frac(jlev-1,jcol) &
                           &    - pair_cloud_cover(jlev-1,jcol)) then
                         n_layers_to_scale = n_layers_to_scale + 1
                      else
                         do_fill_od_scaling = .true.
                      end if
                   else
                      overhang = cum_cloud_cover(jlev,jcol) &
                           &   - cum_cloud_cover(jlev-1,jcol)
                      if (rand_cloud*(cum_cloud_cover(jlev-1,jcol) &
                           &          - frac(jlev-1,jcol)) &
                           &  < pair_cloud_cover(jlev-1,jcol) - overhang &
                           &    - frac(jlev-1,jcol)) then
                         n_layers_to_scale = 1
                      end if
                   end if
                else
                   do_fill_od_scaling = .true.
                end if

                if (do_fill_od_scaling) then
                   call uniform_distribution_omp(istate, rand_inhom1)
                   do jcloud = max(2,jlev-n_layers_to_scale), jlev-1
                      call uniform_distribution_omp(istate, rand_inhom2)

                      overlap_param_inhom = overlap_param(jcloud-1,jcol)
                      if ( icol_begin <= jcloud-1 .and. &
                           & jcloud-1 < icol_end .and. &
                           & overlap_param(jcloud-1,jcol) > 0.0_jprb) then
                         overlap_param_inhom = overlap_param(jcloud-1,jcol) &
                              &  **(1.0_jprb/decorrelation_scaling)
                      end if

                      if ( jcloud > jlev-n_layers_to_scale .and. &
                           & rand_inhom2 >= overlap_param_inhom) then
                         call uniform_distribution_omp(istate, rand_inhom1)
                      end if

                      wcdf = rand_inhom1 * (sample_ncdf-1) + 1.0_jprb
                      icdf = max(1, min(int(wcdf), sample_ncdf-1))
                      wcdf = max(0.0_jprb, min(wcdf - icdf, 1.0_jprb))

                      wfsd = (fractional_std(jcloud,jcol)-sample_fsd1) &
                           &   * sample_inv_fsd_interval + 1.0_jprb
                      ifsd = max(1, min(int(wfsd), sample_nfsd-1))
                      wfsd = max(0.0_jprb, min(wfsd - ifsd, 1.0_jprb))

                      od_scaling(jg,jcloud,jcol) = &
                           &    (1.0_jprb-wcdf)*(1.0_jprb-wfsd) * sample_val(icdf  ,ifsd)   &
                           &  + (1.0_jprb-wcdf)*          wfsd  * sample_val(icdf  ,ifsd+1) &
                           &  +           wcdf *(1.0_jprb-wfsd) * sample_val(icdf+1,ifsd)   &
                           &  +           wcdf *          wfsd  * sample_val(icdf+1,ifsd+1)
                   end do
                   n_layers_to_scale = 0
                end if
             end do
          end if
       end do
    end do
    !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  end subroutine cloud_generator_block_omp
#endif

end module radiation_cloud_generator_acc
