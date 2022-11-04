! radiation_gaussian_quadrature.F90 - Return nodes and weights for quadrature rule for angular integration
!
! (C) Copyright 2022- ECMWF.
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

module radiation_gaussian_quadrature

  use parkind1, only : jprb

  implicit none
  
  public

  enum, bind(c)
    enumerator IQuadratureLegendre, IQuadratureLaguerre, IQuadratureJacobi, &
         &     IQuadratureElsasser, IQuadratureChandrasekhar, IQuadratureTransOpt, &
         &     IQuadratureOptimal
  end enum

  character(len=*), parameter :: QuadratureName(0:6) = [ 'Legendre     ', &
       &                                                 'Laguerre     ', &
       &                                                 'Jacobi       ', &
       &                                                 'Elsasser     ', &
       &                                                 'Chandrasekhar', &
       &                                                 'TransOpt     ', &
       &                                                 'Optimal      ' ]
  
contains
  
  !---------------------------------------------------------------------
  ! Return the nodes and weights for choosing angles for the discrete
  ! ordinate method. The nodes returned are values of the cosine of
  ! the zenith angle in the range 0-1
  subroutine gaussian_quadrature(quadrature_type, norder, nodes, weights, jacobi_power)

    use radiation_io, only : nulerr, radiation_abort
    
    integer, intent(in) :: quadrature_type
    integer, intent(in) :: norder
    real(jprb), intent(out) :: nodes(norder)
    real(jprb), intent(out) :: weights(norder)
    integer, intent(in), optional :: jacobi_power

    integer :: jacobi_power_local
    real(jprb) :: tmp(norder)
    real(jprb) :: nodes_chand(2*norder)
    real(jprb) :: weights_chand(2*norder)

    if (norder < 1) then
      write(nulerr,'(a)') '*** Error: number of nodes in Gaussian quadrature must be >= 1'
      call radiation_abort()
    end if

    nodes = 0.0_jprb
    weights = 0.0_jprb

    if (quadrature_type == IQuadratureLegendre) then
      ! This is the standard "double-Gauss" quadrature of Sykes (1951)
      ! in which each hemisphere is divided as if there is no weight
      ! in the range 0<mu<1 (even though fluxes are computed by
      ! integrating radiances weighted by mu), i.e. Gauss-Legendre
      ! quadrature
      call cgqf(norder, 1, 0.0_jprb ,0.0_jprb, 0.0_jprb, 1.0_jprb, nodes, weights)
      
    else if (quadrature_type == IQuadratureLaguerre) then
      ! Change of variables t=-2ln(mu), results in
      ! mu*dmu=constant*exp(-t)*dt and weighting by exp(-t) is
      ! Gauss-Laguerre quadrature)
      call cgqf(norder, 5, 0.0_jprb ,0.0_jprb, 0.0_jprb, 1.0_jprb, nodes, weights)
      nodes = exp(-0.5_jprb * nodes)
      ! Flip the nodes & weights to be in order of increasing mu
      tmp = nodes(norder:1:-1)
      nodes = tmp
      tmp = weights(norder:1:-1)
      weights = tmp
      ! Assuming the downstream application will multiply the weights
      ! by mu, we divide by mu then rescale the weights
      weights = weights * 0.5_jprb / sum(weights)
      weights = weights / nodes
      
    else if (quadrature_type == IQuadratureElsasser) then
      ! Elsasser (1942) used a diffusivity of 1.66 in two-stream
      ! radiative transfer, i.e. mu=1/1.66
      if (norder /= 1) then
        write(nulerr,'(a)') '*** Error: Elsasser quadrature requires norder==1'
        call radiation_abort()
      end if
      nodes(1) = 1.0_jprb / 1.66_jprb
      weights(1) = 0.5_jprb / nodes(1)
      
    else if (quadrature_type == IQuadratureJacobi) then
      ! Change of variables to s=mu^(2/(jacobi_power+1)), results in
      ! mu*dmu=constant*s^jacobi_power*ds, and weighting by
      ! s^jacobi_power is Gauss-Jacobi quadrature
      if (present(jacobi_power)) then
        jacobi_power_local = jacobi_power
      else
        jacobi_power_local = 0
      end if
      call cgqf(norder, 4, 0.0_jprb, real(jacobi_power_local,jprb), &
           &    0.0_jprb, 1.0_jprb, nodes, weights)
      nodes = nodes**(0.5_jprb*real(jacobi_power+1,jprb))
      ! Assuming the downstream application will multiply the weights
      ! by mu, we divide by mu then rescale the weights
      weights = weights * 0.5_jprb / sum(weights)
      weights = weights / nodes
      !weights = weights / sum(weights)

    else if (quadrature_type == IQuadratureChandrasekhar) then
      ! Chandrasekhar used Gauss-Legendre quadrature acros the entire
      ! sphere (-1<mu<1)
      call cgqf(norder*2, 1, 0.0_jprb ,0.0_jprb, -1.0_jprb, 1.0_jprb, nodes_chand, weights_chand)
      nodes = nodes_chand(norder+1:2*norder)
      weights = weights_chand(norder+1:2*norder)
      weights = weights * 0.5_jprb / sum(nodes*weights)

    else if (quadrature_type == IQuadratureTransOpt) then
      ! Optimized values to minimize the error in transmission over a
      ! wide range of paths
      if (norder == 1) then
        nodes(1) = 1.0_jprb / 1.6145_jprb
        weights(1) = 0.5_jprb / nodes(1)
        !nodes(1) = 0.619392_jprb
        !weights(1) = 0.807244_jprb
      else if (norder == 2) then
        nodes(1:2)   = [ 0.2599_jprb, 0.7797_jprb ]
        weights(1:2) = [ 0.42894_jprb, 0.49829_jprb ]
        !nodes(1:2) = [0.225988_jprb, 0.756382_jprb]
        !weights(1:2) = [0.399858_jprb, 0.541575_jprb]
      else if (norder == 3) then
        nodes(1:3) = [0.0906266_jprb, 0.363077_jprb, 0.817142_jprb]
        weights(1:3) = [0.162957_jprb, 0.39496_jprb, 0.418325_jprb]
      else if (norder == 4) then
        nodes(1:4) = [0.041851_jprb, 0.168645_jprb, 0.454197_jprb, 0.852579_jprb]
        weights(1:4) = [0.0742774_jprb, 0.193929_jprb, 0.377014_jprb, 0.343602_jprb]
      else if (norder == 5) then
        nodes(1:5) = [0.0213478_jprb, 0.0845132_jprb, 0.233272_jprb, 0.520621_jprb, 0.876068_jprb]
        weights(1:5) = [0.0374927_jprb, 0.0963651_jprb, 0.21195_jprb, 0.355751_jprb, 0.292674_jprb]
      else if (norder == 6) then
        nodes(1:6) = [0.0116689_jprb, 0.0455395_jprb, 0.124417_jprb, &
             &        0.288027_jprb, 0.571625_jprb, 0.892872_jprb]
        weights(1:6) = [0.0203675_jprb, 0.0512055_jprb, 0.113293_jprb, &
             &          0.221709_jprb, 0.334774_jprb, 0.255481_jprb]
      else if (norder == 7) then
        nodes(1:7) = [0.0067115_jprb, 0.0259702_jprb, 0.0700749_jprb, 0.161642_jprb, &
             &        0.335306_jprb, 0.612455_jprb, 0.905608_jprb]
        weights(1:7) = [0.0116868_jprb, 0.0289031_jprb, 0.0629935_jprb, 0.126164_jprb, &
             &          0.226467_jprb, 0.315125_jprb, 0.226839_jprb]
      else if (norder == 8) then
        nodes(1:8) = [0.00411735_jprb, 0.0158923_jprb, 0.0425107_jprb, 0.096991_jprb, &
             &        0.201017_jprb, 0.383466_jprb, 0.652213_jprb, 0.917566_jprb]
        weights(1:8) = [0.0071849_jprb, 0.0175886_jprb, 0.0377575_jprb, 0.0747951_jprb, &
             &          0.138424_jprb, 0.229277_jprb, 0.294212_jprb, 0.199656_jprb]
      else
        write(nulerr,'(a)') '*** Error: Transmission-optimized quadrature requires norder <= 8'
        call radiation_abort()
      end if
    else if (quadrature_type == IQuadratureOptimal) then
      ! Optimized values to minimize the error in fluxes and heating
      ! rates for the MLS standard atmosphere
      if (norder == 1) then
        nodes(1:1) = [0.606205]
        weights(1:1) = [0.824804]
      else if (norder == 2) then
        nodes(1:2) = [0.189901, 0.738542]
        weights(1:2) = [0.383466, 0.578409]
      else if (norder == 3) then
        nodes(1:3) = [0.0783442, 0.373655, 0.832172]
        weights(1:3) = [0.170496, 0.424737, 0.394074]
      else if (norder == 4) then
        nodes(1:4) = [0.0366585, 0.178189, 0.490449, 0.872313]
        weights(1:4) = [0.0820073, 0.219265, 0.391545, 0.304811]
      else if (norder == 5) then
        nodes(1:5) = [0.0255436, 0.101581, 0.268106, 0.56644, 0.893844]
        weights(1:5) = [0.0542268, 0.109763, 0.232149, 0.349368, 0.254327]
      else if (norder == 6) then
        nodes(1:6) = [0.0251896, 0.0978779, 0.251806, 0.498906, 0.681858, 0.908207]
        weights(1:6) = [0.053093, 0.103407, 0.210043, 0.24434, 0.174694, 0.214303]
      else if (norder == 7) then
        nodes(1:7) = [0.0246248, 0.0912663, 0.224044, 0.428732, 0.617081, 0.762264, 0.932376]
        weights(1:7) = [0.0511669, 0.0917827, 0.176662, 0.212719, 0.160546, 0.151536, 0.155521]
      else if (norder == 8) then
        nodes(1:8) = [0.0249335, 0.094519, 0.234396, 0.42446, 0.544502, 0.67376, 0.814058, 0.935843]
        weights(1:8) = [0.0522084, 0.0969639, 0.184164, 0.149288, 0.13727, 0.125729, 0.118452, 0.135833]
      else
        write(nulerr,'(a)') '*** Error: Optimized quadrature requires norder <= 8'
        call radiation_abort()
      end if
    end if
    
  end subroutine gaussian_quadrature


  subroutine self_test(maxorder)

    use radiation_io, only : nulout

    integer, parameter :: NQuad = 11
    integer, parameter :: QuadList(NQuad) = [ IQuadratureLegendre, IQuadratureLaguerre, &
         &  IQuadratureChandrasekhar, &
         &  IQuadratureJacobi, IQuadratureJacobi, IQuadratureJacobi, IQuadratureJacobi, &
         &  IQuadratureJacobi, IQuadratureJacobi, IQuadratureJacobi, IQuadratureJacobi ]
    integer, parameter :: JacobiPower(NQuad) = [0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7]
    
    integer, intent(in) :: maxorder
    integer :: jorder, jquad
    real(jprb) :: nodes(maxorder), weights(maxorder)

    do jorder = 1,maxorder
      write(nulout,'(a,i0)') 'Number of nodes ', jorder
      do jquad = 1,NQuad
        call gaussian_quadrature(QuadList(jquad), jorder, nodes, weights, &
             &  jacobi_power=JacobiPower(jquad))
        write(nulout,*) nodes(1:jorder), '|', weights(1:jorder), '|', trim(QuadratureName(QuadList(jquad))), JacobiPower(jquad)
      end do
    end do

  end subroutine self_test
  

end module radiation_gaussian_quadrature
