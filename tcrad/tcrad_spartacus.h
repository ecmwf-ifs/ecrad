! tcrad_spartacus.h - Layer solutions using SPARTACUS method for TCRAD -*- f90 -*-
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
! This file is included in the modules specifying the NREGION
! parameter (typically 2 or 3) which makes this routine use a
! "doubleclouds" or "tripleclouds" assumption.
!


!---------------------------------------------------------------------
! Calculate the transmittance of each layer along a path with consine
! of zenith angle "mu", including exchange between regions, as well as
! (optionally) the emission up from the top of the layer and down
! through its base
subroutine calc_radiance_trans_source_3d(nspec, nlev, &
     &  mu, region_fracs, region_edge_area, od, &
     &  transmittance, &
     &  rate_up_top, rate_up_base, rate_dn_top, rate_dn_base, &
     &  source_up, source_dn)

  use parkind1, only           : jprb
  use yomhook,  only           : lhook, dr_hook

  real(jprb), parameter :: PI = acos(-1.0_jprb)

  ! Inputs

  ! Number of spectral intervals and levels
  integer(jpim), intent(in) :: nspec, nlev

  ! Cosine of the zenith angle (positive)
  real(jprb) :: mu
  
  ! Fraction of the gridbox occupied by each region (summing to 1)
  ! at each level
  real(jprb), intent(in), dimension(NREGION,nlev) :: region_fracs

  ! Area of the vertical interface between each pair of regions,
  ! divided by the horizontal area of the domain. For 3 regions there
  ! are two areas: between regions 1 and 2 and between regions 2 and 3
  ! (regions 1 and 3 are assumed not to touch).
  real(jprb) :: region_edge_area(NREGION-1,nlev)

  ! Optical depth in each region and layer
  real(jprb), intent(in), dimension(nspec,NREGION,nlev) :: od
  
  ! Outputs

  ! Layer transmittance at the requested zenith angle
  real(jprb), intent(out), dimension(nspec,NREGION,NREGION,nlev) :: transmittance

  ! Optional inputs

  ! Rate of emission/scattering upwards or downwards at the top or
  ! base of each layer and region, in Watts of power per square metre
  ! of the entire gridbox. These terms are available instead of
  ! source_up and source_dn for the 3D radiance calculation. Note that
  ! this is the rate of emission/scattering along the radiance
  ! direction per unit optical depth in that layer region, but since
  ! optical depth is dimensionless, these rates still have units of W
  ! m-2.
  real(jprb), intent(in), dimension(nspec,NREGION,nlev), optional &
       &  :: rate_up_top, rate_up_base, rate_dn_top, rate_dn_base

  ! Optional outputs

  ! Source term up from the top of the layer or down from its base, in
  ! Watts of power per square metre of the entire gridbox, so the
  ! energy is scaled by the size of each region. Since the user may
  ! only require a radiance up or down, these output arguments are
  ! optional.
  real(jprb), intent(out), dimension(nspec,NREGION,nlev), optional &
       &  :: source_up, source_dn

  ! Local variables

  ! SPARTACUS "gamma" matrix multiplied by the layer depth, describing
  ! the rate of exchange between regions multiplied by the layer depth
  real(jprb), dimension(nspec,NREGION,NREGION) :: gamma_z1

  ! Inverse of gamma_z1 matrix
  real(jprb), dimension(nspec,NREGION,NREGION) :: inv_gamma_z1

  ! Working variable
  real(jprb), dimension(nspec,NREGION) :: d0, d_prime, source_diff, source_far
  real(jprb) :: coeff

  ! Secant and tangent of zenith angle
  real(jprb) :: secant_za, tan_za

  ! Exchange rate - the f_ab type terms 
  real(jprb) :: exchange_rate
 
  ! Loop indices
  integer(jpim) :: jlev, jreg, jspec
    
  real(jprb) :: hook_handle

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_3d',0,hook_handle)

  secant_za = 1.0_jprb / mu
  tan_za    = sqrt(1.0_jprb - mu*mu) / mu

  do jlev = 1,nlev

    if (region_fracs(1,jlev) < 1.0_jprb) then
      ! Cloudy layer: create gamma_z1 matrix

      gamma_z1 = 0.0_jprb
      ! Loss due to attenuation
      do jreg = 1,NREGION
        gamma_z1(:,jreg,jreg) = -od(:,jreg,jlev)*secant_za
      end do
!      print *, 'jlev = ', jlev
!      print *, 'gamma_z1 = ',gamma_z1(1,:,:)
      ! Loss and gain due to exchange between regions, the loop here
      ! being on the interfaces between regions, which is one less
      ! than the number of regions
      do jreg = 1,NREGION-1
!        print *,'region_edge_area=',region_edge_area(jreg,jlev)
        if (region_edge_area(jreg,jlev) > 0.0_jprb) then
          ! Radiation lost from region jreg and gained in region jreg+1
          exchange_rate = region_edge_area(jreg,jlev) * tan_za / (PI * region_fracs(jreg,jlev))
          gamma_z1(:,jreg,jreg) = gamma_z1(:,jreg,jreg) - exchange_rate
          gamma_z1(:,jreg+1,jreg) = gamma_z1(:,jreg+1,jreg) + exchange_rate
          ! Radiation lost from region jreg+1 and gained in region jreg
          exchange_rate = region_edge_area(jreg,jlev) * tan_za / (PI * region_fracs(jreg+1,jlev))
          gamma_z1(:,jreg+1,jreg+1) = gamma_z1(:,jreg+1,jreg+1) - exchange_rate
          gamma_z1(:,jreg,jreg+1) = gamma_z1(:,jreg,jreg+1) + exchange_rate
!          print *,'exchange_rate = ', exchange_rate
        end if
      end do
!      print *, 'gamma_z1=',gamma_z1(1,:,:)
      call expm_tridiagonal(nspec, gamma_z1, transmittance(:,:,:,jlev))
      !transmittance(:,:,:,jlev) = 0.0_jprb
      !transmittance(:,1,1,jlev) = exp(gamma_z1(:,1,1))
      !transmittance(:,2,2,jlev) = exp(gamma_z1(:,2,2))
      !transmittance(:,3,3,jlev) = exp(gamma_z1(:,3,3))
      call inv_tridiagonal (nspec, gamma_z1, inv_gamma_z1)

!      print *,'gamma_z1='
!      print *,'  ',gamma_z1(1,1,:)
!      print *,'  ',gamma_z1(1,2,:)
!      print *,'  ',gamma_z1(1,3,:)

      if (present(source_up)) then
        source_far  = rate_up_base(:,:,jlev) * od(:,:,jlev) * secant_za
        source_diff = (rate_up_base(:,:,jlev)-rate_up_top(:,:,jlev)) * od(:,:,jlev) * secant_za
        d_prime = mat_x_vec(nspec, inv_gamma_z1, source_diff)
        d0 = mat_x_vec(nspec, inv_gamma_z1, d_prime - source_far)
        source_up(:,:,jlev) = d0 + d_prime - mat_x_vec(nspec, transmittance(:,:,:,jlev), d0)
        if (source_up(1,1,jlev) < 0.0_jprb) then
          print *, 'level = ', jlev
          print *, 'source_up = ', source_up(1,:,jlev)
          print *, 'rate_up_base = ', rate_up_base(1,:,jlev)
          print *, 'rate_up_top = ', rate_up_top(1,:,jlev)
          print *, 'od = ', od(1,:,jlev)
          print *, 'secant = ', secant_za
          print *, 'transmittance = ', transmittance(1,:,:,jlev)
          print *, 'inv_gamma_z1 = ', inv_gamma_z1(1,:,:)
          print *, 'gamma_z1 = '
          do jreg = 1,NREGION
            print *, '  ', gamma_z1(1,jreg,:)
          end do
        end if
      end if

      if (present(source_dn)) then
        source_far  = rate_dn_top(:,:,jlev) * od(:,:,jlev) * secant_za
        source_diff = (rate_dn_top(:,:,jlev)-rate_dn_base(:,:,jlev)) * od(:,:,jlev) * secant_za
        d_prime = mat_x_vec(nspec, inv_gamma_z1, source_diff)
        d0 = mat_x_vec(nspec, inv_gamma_z1, d_prime - source_far)
        source_dn(:,:,jlev) = d0 + d_prime - mat_x_vec(nspec, transmittance(:,:,:,jlev), d0)
      end if

    else
      ! Clear layer...
      transmittance(:,:,:,jlev) = 0.0_jprb
      transmittance(:,1,1,jlev) = exp(-od(:,1,jlev)*secant_za)
      if (present(source_dn)) then
        do jspec = 1,nspec
          coeff = (rate_dn_top(jspec,1,jlev) - rate_dn_base(jspec,1,jlev)) &
               &   * mu / od(jspec,1,jlev)
          source_dn(jspec,1,jlev) = coeff + rate_dn_base(jspec,1,jlev) &
               - transmittance(jspec,1,1,jlev) * (coeff + rate_dn_top(jspec,1,jlev))
        end do
        source_dn(:,2:NREGION,jlev) = 0.0_jprb
      end if
      if (present(source_up)) then
        do jspec = 1,nspec
          coeff = (rate_up_base(jspec,1,jlev) - rate_up_top(jspec,1,jlev)) &
               &   * mu / od(jspec,1,jlev)
          source_up(jspec,1,jlev) = coeff + rate_up_top(jspec,1,jlev) &
               - transmittance(jspec,1,1,jlev) * (coeff + rate_up_base(jspec,1,jlev))
        end do
        source_up(:,2:NREGION,jlev) = 0.0_jprb
      end if
    end if

  end do

  if (lhook) call dr_hook('tcrad:calc_radiance_trans_source_3d',1,hook_handle)

end subroutine calc_radiance_trans_source_3d

