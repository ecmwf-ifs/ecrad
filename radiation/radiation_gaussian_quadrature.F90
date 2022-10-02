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
         &     IQuadratureElsasser, IQuadratureChandrasekhar
  end enum

  character(len=*), parameter :: QuadratureName(0:4) = [ 'Legendre     ', &
       &                                                 'Laguerre     ', &
       &                                                 'Jacobi       ', &
       &                                                 'Elsasser     ', &
       &                                                 'Chandrasekhar' ]
  
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
      weights(1) = 1.0_jprb
      
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
