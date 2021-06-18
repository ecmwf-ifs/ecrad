! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMHOOKSTACK

! Used by dr_hook_util to monitor thread stack usage 
! Need "export STACKCHECK=yes"

USE PARKIND1  ,ONLY : JPIM     ,JPRB,      JPIB

IMPLICIT NONE

SAVE

INTEGER(KIND=JPIB), ALLOCATABLE :: ISAVE(:) 
INTEGER(KIND=JPIB), ALLOCATABLE :: IMAXSTACK(:) 
LOGICAL,   ALLOCATABLE :: LL_THREAD_FIRST(:)
CHARACTER(LEN=3)       :: CSTACK

END MODULE YOMHOOKSTACK

