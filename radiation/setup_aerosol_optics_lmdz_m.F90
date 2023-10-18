module setup_aerosol_optics_lmdz_m

  implicit none

contains

  subroutine setup_aerosol_optics_lmdz(ao, file_name, iverbose)

    ! Read aerosol optical properties. Note differences with
    ! "radiation_aerosol_optics_data::setup_aerosol_optics":

    ! -- The input NetCDF file is not flat, it contains NetCDF groups.

    ! -- We do not define ao%ssa_mono_phobic, ao%g_mono_phobic,
    ! ao%lidar_ratio_mono_phobic, ao%ssa_mono_philic,
    ! ao%g_mono_philic, ao%lidar_ratio_mono_philic. They are not in
    ! the input NetCDF file and they are not used by ECRad.

    ! -- We do not define ao%description_phobic_str and
    ! ao%description_philic_str. We just leave the initialization
    ! value, which is a blank.

    ! -- We have to cshift the shortwave fields because the the
    ! shortwave bands are in ascending order in the NetCDF file while
    ! they are not in ECRad.

    use radiation_aerosol_optics_data, only: aerosol_optics_type, &
         IAerosolClassUndefined
    use netcdf95, only: nf95_open, nf95_inq_grp_full_ncid, nf95_close, &
         nf95_inq_dimid, nf95_inq_varid, nf95_inquire_dimension, &
         nf95_get_var, nf95_gw_var, nf95_nowrite

    class(aerosol_optics_type), intent(inout):: ao
    ! Could be intent(out) here but must conform to the interface of
    ! radiation_aerosol_optics_data::setup_aerosol_optics.

    character(len=*), intent(in):: file_name
    ! NetCDF file containing the aerosol optics data

    integer, intent(in), optional:: iverbose
    ! Not used but we must conform to the interface of
    ! radiation_aerosol_optics_data::setup_aerosol_optics.

    ! Local:
    integer ncid, grpid, dimid, varid

    !-----------------------------------------------------------------------

    ao%use_hydrophilic = .true.
    ao%use_monochromatic = .true.
    call nf95_open(file_name, nf95_nowrite, ncid)
    call nf95_inq_grp_full_ncid(ncid, "Hydrophilic", grpid)
    call nf95_inq_dimid(grpid, "hur", dimid)
    call nf95_inquire_dimension(grpid, dimid, nclen = ao%nrh)
    allocate(ao%rh_lower(ao%nrh))
    call nf95_inq_varid(grpid, "hur_bounds", varid)
    call nf95_get_var(grpid, varid, ao%rh_lower, count_nc = [1, ao%nrh])

    ! Hydrophilic/LW_bands:
    call nf95_inq_grp_full_ncid(ncid, "Hydrophilic/LW_bands", grpid)
    call nf95_inq_varid(grpid, "asymmetry", varid)
    call nf95_gw_var(grpid, varid, ao%g_lw_philic)
    call nf95_inq_varid(grpid, "single_scat_alb", varid)
    call nf95_gw_var(grpid, varid, ao%ssa_lw_philic)
    call nf95_inq_varid(grpid, "mass_ext", varid)
    call nf95_gw_var(grpid, varid, ao%mass_ext_lw_philic)

    ! Hydrophilic/SW_bands:
    call nf95_inq_grp_full_ncid(ncid, "Hydrophilic/SW_bands", grpid)
    call nf95_inq_varid(grpid, "asymmetry", varid)
    call nf95_gw_var(grpid, varid, ao%g_sw_philic)
    ao%g_sw_philic = cshift(ao%g_sw_philic, 1)
    call nf95_inq_varid(grpid, "single_scat_alb", varid)
    call nf95_gw_var(grpid, varid, ao%ssa_sw_philic)
    ao%ssa_sw_philic = cshift(ao%ssa_sw_philic, 1)
    call nf95_inq_varid(grpid, "mass_ext", varid)
    call nf95_gw_var(grpid, varid, ao%mass_ext_sw_philic)
    ao%mass_ext_sw_philic = cshift(ao%mass_ext_sw_philic, 1)

    ! Hydrophilic/Monochromatic:
    call nf95_inq_grp_full_ncid(ncid, "Hydrophilic/Monochromatic", grpid)
    call nf95_inq_varid(grpid, "mass_ext", varid)
    call nf95_gw_var(grpid, varid, ao%mass_ext_mono_philic)

    ! Hydrophobic/LW_bands:
    call nf95_inq_grp_full_ncid(ncid, "Hydrophobic/LW_bands", grpid)
    call nf95_inq_varid(grpid, "asymmetry", varid)
    call nf95_gw_var(grpid, varid, ao%g_lw_phobic)
    call nf95_inq_varid(grpid, "single_scat_alb", varid)
    call nf95_gw_var(grpid, varid, ao%ssa_lw_phobic)
    call nf95_inq_varid(grpid, "mass_ext", varid)
    call nf95_gw_var(grpid, varid, ao%mass_ext_lw_phobic)

    ! Hydrophobic/SW_bands:
    call nf95_inq_grp_full_ncid(ncid, "Hydrophobic/SW_bands", grpid)
    call nf95_inq_varid(grpid, "asymmetry", varid)
    call nf95_gw_var(grpid, varid, ao%g_sw_phobic)
    ao%g_sw_phobic = cshift(ao%g_sw_phobic, 1)
    call nf95_inq_varid(grpid, "single_scat_alb", varid)
    call nf95_gw_var(grpid, varid, ao%ssa_sw_phobic)
    ao%ssa_sw_phobic = cshift(ao%ssa_sw_phobic, 1)
    call nf95_inq_varid(grpid, "mass_ext", varid)
    call nf95_gw_var(grpid, varid, ao%mass_ext_sw_phobic)
    ao%mass_ext_sw_phobic = cshift(ao%mass_ext_sw_phobic, 1)

    ! Hydrophobic/Monochromatic:
    call nf95_inq_grp_full_ncid(ncid, "Hydrophobic/Monochromatic", grpid)
    call nf95_inq_varid(grpid, "mass_ext", varid)
    call nf95_gw_var(grpid, varid, ao%mass_ext_mono_phobic)

    call nf95_close(ncid)

    ! Get array sizes
    ao%n_bands_lw = size(ao%mass_ext_lw_phobic, 1)
    ao%n_bands_sw = size(ao%mass_ext_sw_phobic, 1)
    ao%n_mono_wl = size(ao%mass_ext_mono_phobic, 1)
    ao%n_type_phobic = size(ao%mass_ext_lw_phobic, 2)
    ao%n_type_philic = size(ao%mass_ext_lw_philic, 3)

  end subroutine setup_aerosol_optics_lmdz

end module setup_aerosol_optics_lmdz_m
