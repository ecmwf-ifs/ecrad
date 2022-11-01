! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOECLD

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOECLD* - CONTROL PARAMETERS FOR DIAGNOSTIC CLOUDS
!     -----------------------------------------------------------------

REAL(KIND=JPRB) :: RDECORR_CF = 2.0_JPRB
REAL(KIND=JPRB) :: RDECORR_CW = 1.0_JPRB

END MODULE YOECLD
