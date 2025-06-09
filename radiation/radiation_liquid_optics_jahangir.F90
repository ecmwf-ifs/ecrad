! radiation_liquid_optics_jahangir.F90
!
! (C) Copyright 2014- ECMWF.
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
! License: see the COPYING file for details
!

module radiation_liquid_optics_jahangir

  implicit none

  integer, parameter :: NLiqOpticsCoeffsJahangir = 33

contains

  !---------------------------------------------------------------------
  subroutine calc_liq_optics_jahangir(nb, coeff, lwp, re, od, scat_od, g)

    use parkind1, only : jprb
    !use yomhook,  only : lhook, dr_hook

    ! Number of bands
    integer, intent(in)  :: nb
    ! Coefficients read from a data file
    real(jprb), intent(in) :: coeff(:,:)
    ! Liquid water path (kg m-2) and effective radius (m)
    real(jprb), intent(in) :: lwp, re
    ! Total optical depth, scattering optical depth and asymmetry factor
    real(jprb), intent(out) :: od(nb), scat_od(nb), g(nb)
    integer :: nn
    !real(jprb) :: hook_handle

    real(jprb) :: bands(14)
    real(jprb) :: thres(14)
    !bands =(/0.7,1.01,1.26,1.45,1.78,2.045,2.325,2.78,3.45/)
    bands= (/3.45,2.78,2.325,2.045,1.78,1.45,1.26,1.01,0.7,0.533,0.39,0.3,0.23,3.92/)
    thres(1:nb) = re*1E6 /bands(1:nb)
    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',0,hook_handle)
    do nn=1,10
      if (thres(nn)<2) then
        od(nn) = lwp *  (coeff(nn,1) + re*(coeff(nn,2) + re*coeff(nn,3))) &
             &  / ((1.0_jprb + re*(coeff(nn,4) + re*(coeff(nn,5) &
             &  + re*coeff(nn,6))))*1E6)
        g(nn) = (coeff(nn,12) + re*(coeff(nn,13) + re*coeff(nn,14))) &
             &  / (1.0_jprb + re*(coeff(nn,15) + re*coeff(nn,16)))
      else
        od(nn) = lwp *  (coeff(nn,17) + re*(coeff(nn,18) + re*coeff(nn,19))) &
             &  / ((1.0_jprb + re*(coeff(nn,20) + re*(coeff(nn,21) &
             &  + re*coeff(nn,22))))*1E6)
        g(nn) = (coeff(nn,28) + re*(coeff(nn,29) + re*coeff(nn,30))) &
             &  / (1.0_jprb + re*(coeff(nn,31) + re*coeff(nn,32)))
      end if

      if (coeff(nn,33)>0 .and. re*1E6>coeff(nn,33)) then
        scat_od(nn) = od(nn) * (1.0_jprb &
             &  - (coeff(nn,23) + re*(coeff(nn,24) + re*coeff(nn,25))) &
             &  / (1.0_jprb + re*(coeff(nn,26) + re*coeff(nn,27))))
      else if (coeff(nn,33)>0 .and. re*1E6<=coeff(nn,33)) then
        scat_od(nn) = od(nn) * (1.0_jprb &
             &  - (coeff(nn,7) + re*(coeff(nn,8) + re*coeff(nn,9))) &
             &  / (1.0_jprb + re*(coeff(nn,10) + re*coeff(nn,11))))
      else
        scat_od(nn) = od(nn) * (1.0_jprb &
             &  - (coeff(nn,7) + re*(coeff(nn,8) + re*coeff(nn,9))) &
             &  / (1.0_jprb + re*(coeff(nn,10) + re*coeff(nn,11))))
      end if

    end do


    do nn=11,14
      od(nn) = lwp *  (coeff(nn,1) + re*(coeff(nn,2) + re*coeff(nn,3))) &
           &  / ((1.0_jprb + re*(coeff(nn,4) + re*(coeff(nn,5) &
           &  + re*coeff(nn,6))))*1E6)
      g(nn) = (coeff(nn,12) + re*(coeff(nn,13) + re*coeff(nn,14))) &
           &  / (1.0_jprb + re*(coeff(nn,15) + re*coeff(nn,16)))

      if (coeff(nn,33)>0 .and. re*1E6>coeff(nn,33)) then
        scat_od(nn) = od(nn) * (1.0_jprb &
             &  - (coeff(nn,23) + re*(coeff(nn,24) + re*coeff(nn,25))) &
             &  / (1.0_jprb + re*(coeff(nn,26) + re*coeff(nn,27))))
      else if (coeff(nn,33)>0 .and. re*1E6<=coeff(nn,33)) then
        scat_od(nn) = od(nn) * (1.0_jprb &
             &  - (coeff(nn,7) + re*(coeff(nn,8) + re*coeff(nn,9))) &
             &  / (1.0_jprb + re*(coeff(nn,10) + re*coeff(nn,11))))
      else
        scat_od(nn) = od(nn) * (1.0_jprb &
             &  - (coeff(nn,7) + re*(coeff(nn,8) + re*coeff(nn,9))) &
             &  / (1.0_jprb + re*(coeff(nn,10) + re*coeff(nn,11))))
      end if

    end do
   ! do nn=1,14
   !   write(0,*)'YS coeffs Jahan=',nn,re,lwp,coeff(nn,:)
   !   write(20,*)'YS coeffs Jahan=',nn,re,lwp,coeff(nn,:)
   !   write(0,*)'od g scat_od',od(nn),g(nn),scat_od(nn)
   !   write(20,*)'od g scat_od',od(nn),g(nn),scat_od(nn)
   ! enddo

!IF (re< 1.E-6.OR. re>50.E-6 ) THEN
!    write(0,*)'pb RE'
!    write(20,*)'pb RE'
!ENDIF

    !if (lhook) call dr_hook('radiation_liquid_optics_socrates:calc_liq_optics_socrates',1,hook_handle)

  end subroutine calc_liq_optics_jahangir

end module radiation_liquid_optics_jahangir
