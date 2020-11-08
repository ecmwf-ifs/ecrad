! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

FUNCTION GET_PROC_ID() RESULT(PID)

USE PARKIND1  ,ONLY : JPIM
USE MPL_MODULE, ONLY : MPL_RANK
IMPLICIT NONE
INTEGER(KIND=JPIM) :: PID
PID = MPL_RANK

END FUNCTION GET_PROC_ID
