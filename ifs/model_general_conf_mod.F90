! (C) Copyright 2017- ECMWF.
! (C) Copyright 2017- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE MODEL_GENERAL_CONF_MOD
  ! USE TYPE_GEOMETRY, ONLY : GEOMETRY
  ! USE YOMDIMF      , ONLY : TDIMF
  USE YOM_YGFL     , ONLY : TYPE_GFLD
  USE YOMRIP       , ONLY : TRIP
  ! USE YOMMODERRMOD , ONLY : TMODERR
  ! USE TYPE_ACV     , ONLY : TACVDIM
  IMPLICIT NONE

  TYPE MODEL_GENERAL_CONF_TYPE

    ! TYPE(GEOMETRY), POINTER :: GEOM => NULL()

    ! TYPE(TDIMF)             :: YRDIMF                  !! number of fields
    TYPE(TYPE_GFLD)         :: YGFL                    !! gfl descriptors
    TYPE(TRIP)              :: YRRIP                   !! TEMPORARY TREATMENT OF TIME, SHOULD CHANGE AT CY45
    ! TYPE(TACVDIM)           :: YRDIMACV                !! ACV field
    ! TYPE(TMODERR)           :: YRMODERR                !! Model error config

    CONTAINS

    PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION

  END TYPE MODEL_GENERAL_CONF_TYPE

  !---------------------------------------------------------------------

  CONTAINS

  SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
    IMPLICIT NONE
    CLASS(MODEL_GENERAL_CONF_TYPE), INTENT(IN) :: SELF
    INTEGER                       , INTENT(IN) :: KDEPTH
    INTEGER                       , INTENT(IN) :: KOUTNO

    WRITE(KOUTNO,*) REPEAT(' ',KDEPTH) // 'model%yrml_gconf : '
    !CALL SELF%YRDIMF%PRINT(KDEPTH+2,KOUTNO)
    CALL SELF%YRRIP%PRINT(KDEPTH+2,KOUTNO)
    !CALL SELF%YRMODERR%PRINT(KDEPTH+2,KOUTNO)

  END SUBROUTINE PRINT_CONFIGURATION

END MODULE MODEL_GENERAL_CONF_MOD
