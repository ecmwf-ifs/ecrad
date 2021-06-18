! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMDYNCORE

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

PUBLIC

SAVE

! Aqua planet?
LOGICAL, PARAMETER :: LAQUA = .false.
! Small-planet factor
REAL(KIND=JPRB)    :: RPLRG = 1.0


END MODULE YOMDYNCORE
