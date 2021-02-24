! tripleclouds_solver.F90 - Solve for longwave fluxes or radiances with the Tripleclouds assumption
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

module tripleclouds_solver

  use parkind1, only : jpim

  integer(jpim), parameter :: NREGION = 3

contains

#include "matrix_functions.h"

#include "calc_multiregion_flux.h"

end module tripleclouds_solver
