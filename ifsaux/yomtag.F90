! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMTAG

USE PARKIND1  ,ONLY : JPIM

IMPLICIT NONE

PUBLIC

SAVE

! MTAGRAD : tag for communications done in SUECRAD (ECMWF physics)
INTEGER(KIND=JPIM), PARAMETER :: MTAGRAD               =  2800

END MODULE YOMTAG
