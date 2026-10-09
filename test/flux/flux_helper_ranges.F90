! (C) Copyright 2026- ECMWF.
! Licensed under the Apache Licence Version 2.0; see COPYING.
program flux_helper_ranges
  use parkind1, only : jprb
  use omp_lib, only : omp_is_initial_device
  use radiation_flux, only : sum_surf_band_omp, canopy_sw_weights_omp, canopy_lw_nearest_omp
  implicit none
  integer :: failures = 0
  logical :: on_host

  on_host = .true.
  !$omp target map(from:on_host)
  on_host = omp_is_initial_device()
  !$omp end target
  if (on_host) error stop 'Flux helper regression test requires GPU execution'

  call test_case(1, 8, 1, 4, 'leading block')
  call test_case(1, 8, 5, 8, 'later block')
  call test_case(1, 8, 7, 8, 'partial last block')
  call test_case(5, 8, 5, 8, 'matching non-unit allocation')
  call test_case(5, 12, 9, 12, 'subrange of non-unit allocation')
  if (failures /= 0) error stop 'Incorrect flux helper column addressing'

contains
  subroutine test_case(lo, hi, first, last, label)
    integer, intent(in) :: lo, hi, first, last
    character(*), intent(in) :: label
    integer, parameter :: ng=4, nb=2, nc=2
    integer :: map_g(ng), map_e(nb), j, g, b
    real(jprb), parameter :: sentinel=-9999.0_jprb
    real(jprb) :: direct(ng,lo:hi), diffuse(ng,lo:hi), lw(ng,lo:hi)
    real(jprb) :: band_direct(nb,lo:hi), band_total(nb,lo:hi), bands(nb,lo:hi), direct_bands(nb,lo:hi)
    real(jprb) :: canopy_direct(nc,lo:hi), canopy_diffuse(nc,lo:hi), canopy_lw(nc,lo:hi)
    real(jprb) :: expected_direct(nb,lo:hi), expected_total(nb,lo:hi)
    real(jprb) :: expected_canopy_direct(nc,lo:hi), expected_canopy_diffuse(nc,lo:hi)
    real(jprb) :: expected_lw(nc,lo:hi), weights(nc,nb)

    map_g = [1,1,2,2]
    map_e = [2,1]
    weights = reshape([1.0_jprb,0.0_jprb,0.0_jprb,1.0_jprb], [nc,nb])
    do j=lo,hi
      do g=1,ng
        direct(g,j) = real(100*j+g,jprb)
        diffuse(g,j) = real(10*j+g,jprb)
        lw(g,j) = real(1000*j+g,jprb)
      end do
      do b=1,nb
        bands(b,j) = real(100*j+10*b,jprb)
        direct_bands(b,j) = real(10*j+b,jprb)
      end do
    end do
    band_direct=sentinel
    band_total=sentinel
    canopy_direct=sentinel
    canopy_diffuse=sentinel
    canopy_lw=sentinel
    expected_direct=sentinel
    expected_total=sentinel
    expected_canopy_direct=sentinel
    expected_canopy_diffuse=sentinel
    expected_lw=sentinel
    do j=first,last
      expected_direct(:,j)=0.0_jprb
      expected_total(:,j)=0.0_jprb
      expected_lw(:,j)=0.0_jprb
      do g=1,ng
        b=map_g(g)
        expected_direct(b,j)=expected_direct(b,j)+direct(g,j)
        expected_total(b,j)=expected_total(b,j)+diffuse(g,j)
        expected_lw(map_e(b),j)=expected_lw(map_e(b),j)+lw(g,j)
      end do
      expected_total(:,j)=expected_total(:,j)+expected_direct(:,j)
      expected_canopy_direct(:,j)=matmul(weights,direct_bands(:,j))
      expected_canopy_diffuse(:,j)=matmul(weights,bands(:,j))-expected_canopy_direct(:,j)
    end do

    !$omp target data map(to:direct,diffuse,lw,map_g,map_e,weights,bands,direct_bands) &
    !$omp& map(tofrom:band_direct,band_total,canopy_direct,canopy_diffuse,canopy_lw)
    call sum_surf_band_omp(first,last,ng,nb,lo,hi,map_g,direct,diffuse,band_direct,band_total)
    call canopy_sw_weights_omp(first,last,nc,nb,nc,lo,hi,weights,bands,direct_bands,canopy_diffuse,canopy_direct)
    call canopy_lw_nearest_omp(first,last,ng,nb,nc,lo,hi,map_e,map_g,lw,canopy_lw)
    !$omp end target data

    ! Compare the entire allocation: columns outside the range must retain the sentinel.
    call report('sum_surf_band_omp',label, &
         max(maxval(abs(band_direct-expected_direct)),maxval(abs(band_total-expected_total))))
    call report('canopy_sw_weights_omp',label, &
         max(maxval(abs(canopy_direct-expected_canopy_direct)),maxval(abs(canopy_diffuse-expected_canopy_diffuse))))
    call report('canopy_lw_nearest_omp',label,maxval(abs(canopy_lw-expected_lw)))
  end subroutine test_case

  subroutine report(helper,label,error)
    character(*), intent(in) :: helper,label
    real(jprb), intent(in) :: error
    if (error /= 0.0_jprb) failures=failures+1
    print *, trim(helper), ': ', trim(label), ', max error=', error
  end subroutine report
end program flux_helper_ranges
