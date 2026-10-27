! (C) Copyright 2026- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

program cloud_cover_ranges
  use parkind1, only : jpim, jprb
  use radiation_io, only : nulout, radiation_abort
  use radiation_cloud_cover, only : compute_cloud_cover, cloud_cover, OverlapName, &
       & IOverlapMaximumRandom, IOverlapExponential
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: ieee_exceptions, only : ieee_divide_by_zero, ieee_get_flag, ieee_set_flag
  implicit none

  integer(jpim), parameter :: ncol=7, nlev=6
  real(jprb), parameter :: tolerance=100.0_jprb*epsilon(1.0_jprb), sentinel=-1.0_jprb
  real(jprb) :: frac(ncol,nlev), overlap(ncol,nlev-1)
  real(jprb) :: cover(ncol), expected(ncol), partial(ncol), full(ncol)
  integer(jpim) :: lev_top(ncol), lev_base(ncol), jcol, jlev, jscheme, jbeta, nlev_range
  logical :: is_beta, divide_by_zero

  ! Include full-column, shortened, single-level and empty ranges.
  lev_top  = [1,2,3,4,5,6,4]
  lev_base = [6,5,3,6,5,6,3]
  do jlev=1,nlev
    do jcol=1,ncol
      frac(jcol,jlev) = real(mod(jcol+2*jlev,7),jprb)/7.0_jprb
    enddo
  enddo
  overlap = 0.6_jprb

  do jscheme=IOverlapMaximumRandom,IOverlapExponential
    do jbeta=0,1
      is_beta = jbeta == 1
      write(nulout,'(a,a,a,l1)') 'Testing ', trim(OverlapName(jscheme)), ', beta overlap: ', is_beta
      call ieee_set_flag(ieee_divide_by_zero,.false.)

      ! Scalar calculations use only the selected layers and do not pad them.
      expected = 0.0_jprb
      do jcol=1,ncol
        nlev_range = lev_base(jcol)-lev_top(jcol)+1
        if (nlev_range > 0) then
          expected(jcol) = cloud_cover(nlev_range,jscheme, &
               & frac(jcol,lev_top(jcol):lev_base(jcol)), &
               & overlap(jcol,lev_top(jcol):lev_base(jcol)-1),is_beta_overlap=is_beta)
        endif
      enddo

      cover = sentinel
      call compute_cloud_cover(1,ncol,ncol,nlev,jscheme,cover,frac,overlap, &
           & is_beta_overlap=is_beta,lev_top=lev_top,lev_base=lev_base)
      if (any(.not.ieee_is_finite(cover))) call radiation_abort(text='Non-finite level-range cloud cover')
      if (any(abs(cover-expected) > tolerance)) call radiation_abort(text='Incorrect level-range cloud cover')

      ! Splitting a call must preserve results and leave inactive columns alone.
      partial = sentinel
      call compute_cloud_cover(1,3,ncol,nlev,jscheme,partial,frac,overlap, &
           & is_beta_overlap=is_beta,lev_top=lev_top,lev_base=lev_base)
      if (any(partial(4:ncol) /= sentinel)) call radiation_abort(text='Inactive cloud-cover columns changed')
      call compute_cloud_cover(4,ncol,ncol,nlev,jscheme,partial,frac,overlap, &
           & is_beta_overlap=is_beta,lev_top=lev_top,lev_base=lev_base)
      if (any(.not.ieee_is_finite(partial))) call radiation_abort(text='Non-finite partial-column cloud cover')
      if (any(abs(partial-cover) > tolerance)) call radiation_abort(text='Partial-column cloud cover differs')

      ! Omitting the new arguments retains the original full-column result.
      full = sentinel
      call compute_cloud_cover(1,ncol,ncol,nlev,jscheme,full,frac,overlap,is_beta_overlap=is_beta)
      do jcol=1,ncol
        expected(jcol) = cloud_cover(nlev,jscheme,frac(jcol,:),overlap(jcol,:),is_beta_overlap=is_beta)
      enddo
      if (any(.not.ieee_is_finite(full))) call radiation_abort(text='Non-finite full-column cloud cover')
      if (any(abs(full-expected) > tolerance)) call radiation_abort(text='Legacy cloud-cover result differs')

      ! Zero padding can raise division by zero in beta2alpha even when the
      ! returned range result is finite. Check sticky flags in optimized builds,
      ! where floating-point exceptions are not necessarily trapped. Do not
      ! inspect invalid-operation flags: vectorized code may speculatively
      ! evaluate unused branches without changing the returned results.
      call ieee_get_flag(ieee_divide_by_zero,divide_by_zero)
      if (divide_by_zero) call radiation_abort(text='Division by zero in cloud-cover calculation')
    enddo
  enddo
end program cloud_cover_ranges
