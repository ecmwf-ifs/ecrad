! ecrad_config.h - Preprocessor definitions to configure compilation ecRad -*- f90 -*-
!
! (C) Copyright 2023- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
!
! This file should be included in Fortran source files that require
! different optimizations for different architectures, the main
! difference being between the long-vector NEC SX and short-vector
! x86-64 architectures. Feel free to maintain a site-specific version
! of it.
  
#if defined (__SX__) || defined (_OPENACC)
#define DWD_TWO_STREAM_OPTIMIZATIONS
#endif
  
#if defined (__SX__)
#define DWD_REDUCTION_OPTIMIZATIONS
#endif
  
#if defined (__SX__)
#define DWD_VECTOR_OPTIMIZATIONS
#endif
