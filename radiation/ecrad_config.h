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


! A requirement of the generator is that the operation mod(A*X,M) is
! performed with no loss of precision, so type used for A and X must
! be able to hold the largest possible value of A*X without
! overflowing, going negative or losing precision. The largest
! possible value is 48271*2147483647 = 103661183124337. This number
! can be held in either a double-precision real number, or an 8-byte
! integer. Either may be used, but on some hardwares it has been
! found that operations on double-precision reals are faster. Select
! which you prefer by defining USE_REAL_RNG_STATE for double
! precision, or undefining it for an 8-byte integer.
#if defined (__SX__)
#define USE_REAL_RNG_STATE 1
#endif

! Define RNG_STATE_TYPE based on USE_REAL_RNG_STATE, where jprd
! refers to a double-precision number regardless of the working
! precision described by jprb, while jpib describes an 8-byte
! integer
#ifdef USE_REAL_RNG_STATE
#define RNG_STATE_TYPE real(kind=jprd)
#else
#define RNG_STATE_TYPE integer(kind=jpib)
#endif

! In the IFS, an MPI version of easy_netcdf capability is used so that
! only one MPI task reads the data files and shares with the other
! tasks. The MPI version is not used for writing files.

!#define EASY_NETCDF_READ_MPI 1
