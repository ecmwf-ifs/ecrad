! radsurf_homogeneous_vegetation.f90 - Compute radiative transfer in homogeneous vegetation canopy
!
! Copyright (C) 2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_homogeneous_vegetation

contains
  subroutine calc_boundary_conditions_sw(config, tile_fraction, &
       &  canopy_depth, vegetation_optical_depth, vegetation_albedo, &
       &  ground_albedo_diffuse, ground_albedo_direct, &
       &  ref_dif, tra_dif, ref_dir, tra_dir_dif, tra_dir_dir, &
       &  albedo_diffuse_reg, albedo_direct_reg, &
       &  albedo_diffuse_out, albedo_direct_out, &
       &  ext_air, ssa_air, g_air)

    use parkind, only :jprb

    implicit none

  end subroutine calc_boundary_conditions_sw

  subroutine calc_boundary_conditions_lw

  end subroutine calc_boundary_conditions_lw


end module radsurf_homogeneous_vegetation
