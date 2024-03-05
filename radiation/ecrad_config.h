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
! different optimizations or settings for different architectures and
! platforms.  Feel free to maintain a site-specific version of it.

! The following settings turn on optimizations specific to the
! long-vector NEC SX (the short-vector x86-64 architecture is assumed
! otherwise). 

#if defined (__SX__) || defined (_OPENACC)
#define DWD_TWO_STREAM_OPTIMIZATIONS 1
#endif
  
#if defined (__SX__)
#define DWD_REDUCTION_OPTIMIZATIONS 1
#endif
  
#if defined (__SX__)
#define DWD_VECTOR_OPTIMIZATIONS 1
#endif

! In the IFS, an MPI version of easy_netcdf capability is used so that
! only one MPI task reads the data files and shares with the other
! tasks. The MPI version is not used for writing files.

!#define EASY_NETCDF_READ_MPI 1
