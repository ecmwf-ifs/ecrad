! tcrad_matrix.h - Batched matrix operations for TCRAD solver -*- f90 -*-
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
! Treat A as an m-by-m square matrix and b as n NREGION-element vectors
! (with the n dimension varying fastest), and perform n matrix-vector
! multiplications
pure function singlemat_x_vec(n,A,b)

  use parkind1, only : jpim, jprb

  implicit none

  integer,    intent(in)                             :: n
  real(jprb), intent(in), dimension(NREGION,NREGION) :: A
  real(jprb), intent(in), dimension(n,NREGION)       :: b
  real(jprb),             dimension(n,NREGION)       :: singlemat_x_vec
  
  integer    :: j1, j2
  
  ! Array-wise assignment
  singlemat_x_vec = 0.0_jprb
  
  do j1 = 1,NREGION
    do j2 = 1,NREGION
      singlemat_x_vec(:,j1) = singlemat_x_vec(:,j1) + A(j1,j2)*b(:,j2)
    end do
  end do
  
end function singlemat_x_vec

