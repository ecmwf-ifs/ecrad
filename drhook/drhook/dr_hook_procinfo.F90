! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE DR_HOOK_PROCINFO(KMYPROC, KNPROC)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE MPL_MODULE, ONLY : MPL_WORLD_RANK, MPL_WORLD_SIZE
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(OUT) :: KMYPROC, KNPROC
KMYPROC = MPL_WORLD_RANK + 1
KNPROC = MPL_WORLD_SIZE
END SUBROUTINE DR_HOOK_PROCINFO
