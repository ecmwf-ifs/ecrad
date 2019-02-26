! radsurf_3d_vegetation.f90 - Compute radiative transfer in 3D vegetation canopy
!
! Copyright (C) 2018 ECMWF
!
! Author:  Robin Hogan
! Email:   r.j.hogan@ecmwf.int
! License: see the COPYING file for details
!

module radsurf_3d_vegetation

contains
  subroutine calc_boundary_conditions_sw(config, n_albedo_bands, tile_fraction, &
       &  canopy_depth, vegetation_optical_depth, vegetation_albedo, &
       &  ground_albedo_diffuse, ground_albedo_direct, &
       &  ref_dif, tra_dif, ref_dir, tra_dir_dif, tra_dir_dir, &
       &  albedo_diffuse_reg, albedo_direct_reg, &
       &  albedo_diffuse_out, albedo_direct_out, &
       &  ext_air, ssa_air, g_air)

    use parkind, only : jprb
    use radiation_config, only : config_type

    implicit none

    ! Number of regions
    integer, parameter :: nreg = 2

    type(config_type), intent(in) :: config

    integer,           intent(in) :: n_albedo_bands

    ! Fraction of gridbox occupied by this tile
    real(kind=jprb),   intent(in) :: tile_fraction

    ! Depth of vegetation canopy in metres
    real(kind=jprb),   intent(in) :: canopy_depth

    ! Optical properties of vegetation
    real(kind=jprb),   intent(in) :: vegetation_optical_depth
    real(kind=jprb),   intent(in) :: vegetation_albedo(:)  ! Spectral interval

    ! Optical properties of the ground (function of spectral interval)
    real(kind=jprb),   intent(in) :: ground_albedo_diffuse(:)
    real(kind=jprb),   intent(in) :: ground_albedo_direct(:)

    ! Geometric properties
    real(kind=jprb),   intent(in) :: vegetation_fraction
    real(kind=jprb),   intent(in) :: vegetation_normalized_perimeter ! m-1


    ! Intermediate properties to store
    real(kind=jprb),   intent(in), dimension(n_albedo_bands,nreg,nreg) :: ref_dif, tra_dif, tra_dir_dif, tra_dir_dir

    ! Outputs
    real(kind=jprb),   intent(inout) :: albedo_diffuse_out, albedo_direct_out



  end subroutine calc_boundary_conditions_sw

  subroutine calc_boundary_conditions_lw

  end subroutine calc_boundary_conditions_lw


end module radsurf_3d_vegetation
