! tcrad_3region_solver.F90 - Solve for longwave fluxes or radiances with the Tripleclouds assumption
!
! (C) Copyright 2021- ECMWF.
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

module tcrad_3region_solver

  use parkind1, only : jpim

  ! By making the number of regions a parameter NREGION, the compiler
  ! can optimize better.  We also provide a preprocessor parameter as
  ! in one or two cases the operations on 2x2 and 3x3 matrices
  ! (e.g. matrix exponentials) are completely different algorithms.
#define NUM_REGIONS 3
  integer(jpim), parameter :: NREGION = NUM_REGIONS

contains

! The following header files define routines and functions that make
! use of the NREGION parameter

#include "tcrad_region.h"

#include "tcrad_matrix.h"

#include "tcrad_overlap.h"

#include "tcrad_two_stream_flux.h"

#include "tcrad_radiance.h"

#include "tcrad_flux_interface.h"

#include "tcrad_radiance_interface.h"

#undef NUM_REGIONS

end module tcrad_3region_solver
