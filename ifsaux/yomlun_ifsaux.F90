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

#ifdef HAVE_FIAT
USE EC_LUN    ,ONLY : NULOUT, NULERR
#else
USE PARKIND1  ,ONLY : JPIM
#endif

IMPLICIT NONE

SAVE
PRIVATE
PUBLIC :: NULOUT, NULERR

!     ------------------------------------------------------------------

!*    Logical units used by code

!     NULOUT :   output unit
!     NULERR :   unit number for comparison with reference run

#ifndef HAVE_FIAT
INTEGER(KIND=JPIM) :: NULOUT = 6
INTEGER(KIND=JPIM) :: NULERR = 0
#endif

!     ------------------------------------------------------------------
END MODULE YOMLUN_IFSAUX
