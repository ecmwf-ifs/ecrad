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

! Program to test the functionality of the
! radiation_aerosol_optics_description module. Usage is for example:
! "test_aerosol_optics_description aerosol_ifs_49R1.nc"
program test_aerosol_optics_description

  use radiation_aerosol_optics_description

  type(aerosol_optics_description_type) :: aer_desc

  character(len=512) :: file_name

  integer :: istatus

  call get_command_argument(1, file_name, status=istatus)
  if (istatus /= 0) then
    stop 'Usage: test_aerosol_optics_description <aerosol_file.nc>'
  end if

  call aer_desc%read(file_name=file_name)

  ! These two should issue a warning that the aerosol type is ambiguous
  print *, 'DD: ', aer_desc%get_index('DD',.false.)
  print *, 'DD (bin=2): ', aer_desc%get_index('DD',.false.,ibin=2)
  ! Indicate preferred aerosol optical model, after which further
  ! calls will prefer one particular model
  print *, 'preferred_optical_model(DD,Fouquart)'
  call aer_desc%preferred_optical_model('DD','Fouquart')
  print *, 'DD (bin=2): ', aer_desc%get_index('DD',.false.,ibin=2)
  print *, 'DD (bin=2,model=Woodward): ', aer_desc%get_index('DD',.false.,ibin=2,optical_model_str="Woodward")
  print *, 'DD (bin=2,model=Woodward2001): ', aer_desc%get_index('DD',.false.,ibin=2,optical_model_str="Woodward2001")
  ! This should fail to find a match, returning zero
  print *, 'DD (model=Nobody): ', aer_desc%get_index('DD',.false.,optical_model_str="Nobody")
  print *, 'SS (bin=3): ', aer_desc%get_index('SS',.true.,ibin=3)

end program test_aerosol_optics_description
