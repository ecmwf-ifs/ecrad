! interpolation.F90 - linear interpolation routine

module interpolation

contains
  ! 1D linear interpolation: the original data "data_orig" at
  ! wavenumber points "wn_orig" are interpolated to wavenumbers "wn"
  ! and the result put in "data". Both wn and wn_orig should increase
  ! monotonically.
  subroutine interpolate(wn_orig, data_orig, wn, data, out_of_bounds_val)

    use parkind1, only : jprb
    
    implicit none    

    real(jprb), intent(in),  dimension(:) :: wn_orig, data_orig, wn
    real(jprb), intent(out), dimension(:) :: data
    real(jprb), intent(in),  optional     :: out_of_bounds_val
    
    integer :: jwn, nwn, nwn_orig, iwn_orig, iwn

    nwn = size(wn)
    nwn_orig = size(wn_orig)

    iwn_orig = 1
    iwn = 1

    if (present(out_of_bounds_val)) then
      ! Fix values outside bounds of interpolation
      do while (wn(iwn) <= wn_orig(1))
        data(iwn) = out_of_bounds_val
        iwn = iwn+1
      end do
      
      do while (wn(nwn) >= wn_orig(nwn_orig))
        data(nwn) = out_of_bounds_val
        nwn = nwn-1
      end do
    else
      ! Clamp values outside bounds of interpolation
      do while (wn(iwn) <= wn_orig(1))
        data(iwn) = data_orig(1)
        iwn = iwn+1
      end do
      
      do while (wn(nwn) >= wn_orig(nwn_orig))
        data(nwn) = data_orig(nwn_orig)
        nwn = nwn-1
      end do
    end if

    ! Linear interpolation
    do jwn = iwn,nwn
      do while (wn_orig(iwn_orig+1) < wn(jwn))
        iwn_orig = iwn_orig + 1
      end do
      data(jwn) = (data_orig(iwn_orig)*(wn_orig(iwn_orig+1)-wn(jwn)) &
           &     + data_orig(iwn_orig+1)*(wn(jwn)-wn_orig(iwn_orig))) &
           &     / (wn_orig(iwn_orig+1) - wn_orig(iwn_orig))
    end do

  end subroutine interpolate

end module interpolation
