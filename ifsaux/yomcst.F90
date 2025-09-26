! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMCST

USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

PUBLIC

SAVE

! * RPI          : number Pi
REAL(KIND=JPRB), PARAMETER :: RPI = 3.14159265358979323846_JPRB
! * RCLUM        : light velocity
REAL(KIND=JPRB), PARAMETER :: RCLUM = 299792458._JPRB
! * RHPLA        : Planck constant
REAL(KIND=JPRB), PARAMETER :: RHPLA = 6.6260755E-34_JPRB
! * RKBOL        : Bolzmann constant
REAL(KIND=JPRB), PARAMETER :: RKBOL = 1.380658E-23_JPRB
! * RNAVO        : Avogadro number
REAL(KIND=JPRB), PARAMETER :: RNAVO = 6.0221367E+23_JPRB
! * RSIGMA       : Stefan-Bolzman constant
REAL(KIND=JPRB), PARAMETER :: RSIGMA = 5.67037321e-8_JPRB ! W m-2 K-4
! * RG           : gravity constant
REAL(KIND=JPRB), PARAMETER :: RG = 9.80665_JPRB ! m s-2
! * RMD          : dry air molar mass
REAL(KIND=JPRB), PARAMETER :: RMD = 28.9644_JPRB
! * RMV          : vapour water molar mass
REAL(KIND=JPRB), PARAMETER :: RMV = 18.0153_JPRB
! * R            : perfect gas constant
REAL(KIND=JPRB), PARAMETER :: R = RNAVO*RKBOL
! * RD           : R_dry (dry air constant)
REAL(KIND=JPRB), PARAMETER :: RD = 287.058_JPRB! J kg-1 K-1
! * RV           : R_vap (vapour water constant)
REAL(KIND=JPRB), PARAMETER :: RV = 1000._JPRB*R/RMV
! * RMO3         : ozone molar mass
REAL(KIND=JPRB), PARAMETER :: RMO3 = 47.9942_JPRB
! * RTT          : Tt = temperature of water fusion at "pre_n" 
REAL(KIND=JPRB), PARAMETER :: RTT = 273.16_JPRB
! * RLVTT        : RLvTt = vaporisation latent heat at T=Tt
REAL(KIND=JPRB), PARAMETER :: RLVTT = 2.5008E+6_JPRB
! * RLSTT        : RLsTt = sublimation latent heat at T=Tt
REAL(KIND=JPRB), PARAMETER :: RLSTT = 2.8345E+6_JPRB
! * RI0          : solar constant
REAL(KIND=JPRB), PARAMETER :: RI0 = 1366.0_JPRB
! * RETV         : R_vap/R_dry - 1
REAL(KIND=JPRB), PARAMETER :: RETV = RV/RD-1.0_JPRB
! * RMCO2        : CO2 (carbon dioxide) molar mass
REAL(KIND=JPRB), PARAMETER :: RMCO2 = 44.0095_JPRB
! * RMCH4        : CH4 (methane) molar mass
REAL(KIND=JPRB), PARAMETER :: RMCH4 = 16.04_JPRB
! * RMN2O        : N2O molar mass
REAL(KIND=JPRB), PARAMETER :: RMN2O = 44.013_JPRB
! * RMNO2        : NO2 (nitrogen dioxide) molar mass
REAL(KIND=JPRB), PARAMETER :: RMNO2 = 46.01_JPRB
! * RMCFC11      : CFC11 molar mass
REAL(KIND=JPRB), PARAMETER :: RMCFC11 = 137.3686_JPRB
! * RMCFC12      : CFC12 molar mass
REAL(KIND=JPRB), PARAMETER :: RMCFC12 = 120.914_JPRB
! * RMHCFC12     : HCFC22 molar mass
REAL(KIND=JPRB), PARAMETER :: RMHCFC22 = 86.469_JPRB
! * RMCCL4       : CCl4 molar mass
REAL(KIND=JPRB), PARAMETER :: RMCCL4 = 153.823_JPRB

REAL(KIND=JPRB), PARAMETER :: RCPD  = 3.5_JPRB*RD
REAL(KIND=JPRB), PARAMETER :: RLMLT = RLSTT-RLVTT

END MODULE YOMCST
