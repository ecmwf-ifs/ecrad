! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

FUNCTION GET_NUM_THREADS() RESULT(INUMT)
USE PARKIND1  ,ONLY : JPIM
USE OML_MOD, ONLY : OML_NUM_THREADS
IMPLICIT NONE
INTEGER(KIND=JPIM) :: INUMT
INUMT = 1
!$ inumt = OML_NUM_THREADS()
END FUNCTION GET_NUM_THREADS
