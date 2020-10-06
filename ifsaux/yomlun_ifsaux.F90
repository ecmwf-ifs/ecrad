! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

! This is taken from yomlun_ifsaux in the IFS

MODULE YOMLUN_IFSAUX

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

PUBLIC

SAVE

!     ------------------------------------------------------------------

!*    Logical units used by code

!     NULOUT :   output unit
!     NULERR :   unit number for comparison with reference run

INTEGER(KIND=JPIM) :: NULOUT = 6
INTEGER(KIND=JPIM) :: NULERR = 0

!     ------------------------------------------------------------------
END MODULE YOMLUN_IFSAUX
