! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMABRT
USE PARKIND1, ONLY: JPIM
IMPLICIT NONE
PUBLIC
INTEGER(KIND=JPIM), SAVE :: MAB_CNT = 0 ! Must be used with OMP FLUSH inside the OMP CRITICAL regions
END MODULE YOMABRT

