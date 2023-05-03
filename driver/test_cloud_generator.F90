! (C) Copyright 2017- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

program test_cloud_generator

  use parkind1,                  only : jprb
  use radiation_cloud_generator, only : cloud_generator
  use radiation_pdf_sampler,     only : pdf_sampler_type
  use radiation_cloud_cover,     only : &
       & IOverlapMaximumRandom, IOverlapExponentialRandom, IOverlapExponential

  implicit none

  integer,    parameter :: ncol = 2000
  integer,    parameter :: nlev = 137
  integer,    parameter :: i_overlap_scheme = IOverlapExponential
  real(jprb), parameter :: scale_height = 8000.0_jprb
  real(jprb), parameter :: cloud_inhom_decorr_scaling = 0.5_jprb
  real(jprb), parameter :: frac_threshold = 1.0e-6_jprb

  real(jprb) :: cloud_fraction(nlev), overlap_param(nlev-1), fractional_std(nlev)
!  real(jprb) :: pressure_hl(nlev+1)

!  real(jprb) :: decorrelation_length

  real(jprb) :: od_scaling(ncol,nlev)
  real(jprb) :: total_cloud_cover

  integer :: iseed 

  integer :: jcol, jlev

  type(pdf_sampler_type) :: pdf_sampler

  iseed = 1
  cloud_fraction = 0.0_jprb
  overlap_param  = 0.9_jprb
  fractional_std = 1.0_jprb ! Value up to 45R1

  ! Insert cloud layers
  cloud_fraction(115:125) = 0.1_jprb !0.5_jprb
  cloud_fraction(20:100) = 0.1_jprb !0.75_jprb

  call pdf_sampler%setup('data/mcica_gamma.nc', iverbose=0)

  call cloud_generator(ncol, nlev, i_overlap_scheme, &
       &  iseed, frac_threshold, &
       &  cloud_fraction, overlap_param, &
       &  cloud_inhom_decorr_scaling, &
       &  fractional_std, pdf_sampler, &
       &  od_scaling, total_cloud_cover)

  do jlev = 1,nlev
    do jcol = 1,ncol
      write(6,'(f5.2,a)',advance='no') od_scaling(jcol,jlev), ' '
    end do
    write(6,*)
  end do


end program test_cloud_generator
