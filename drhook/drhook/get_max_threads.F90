! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

FUNCTION GET_MAX_THREADS() RESULT(IMAXT)
USE PARKIND1  ,ONLY : JPIM
USE OML_MOD, ONLY : OML_MAX_THREADS
IMPLICIT NONE
INTEGER(KIND=JPIM) :: IMAXT
IMAXT = 1
!$ imaxt = OML_MAX_THREADS()
END FUNCTION GET_MAX_THREADS
