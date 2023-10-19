INPUT = ecrad_meridian

#DRIVER = totalview ../../bin/ecrad -a
DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh

#CONFIG = configCY43R3.nam
#CONFIG = configCY46R1.nam
CONFIG = configCY47R1.nam
CONFIG_ECCKD = configCY49R1.nam

# Typing "make" will run radiation scheme on IFS profiles
all: test

test: test_orig test_ecckd

# Tests of CY47R1 configuration with RRTMG gas optics model and
# various solvers
test_orig: test_default test_noaer test_expran test_expran_vec test_tripleclouds test_spartacus \
	test_spartacus_maxentr

# Tests of ecCKD gas optics model and various solvers
test_ecckd: test_ecckd_mcica test_ecckd_tc test_ecckd_tc_noaer test_ecckd_spartacus

# Default IFS settings: McICA solver with exponential-exponential
# overlap
test_default:
	$(DRIVER) $(CONFIG) $(INPUT).nc $(INPUT)_default_out.nc

# Turn off aerosols
test_noaer:
	$(CHANGENAM) $(CONFIG) config_noaer.nam \
		use_aerosols=false
	$(DRIVER) config_noaer.nam $(INPUT).nc $(INPUT)_noaer_out.nc

# Exponential-random overlap
test_expran:
	$(CHANGENAM) $(CONFIG) config_expran.nam \
		overlap_scheme_name='"Exp-Ran"'
	$(DRIVER) config_expran.nam $(INPUT).nc $(INPUT)_expran_out.nc

# Tripleclouds solver with exponential-random overlap 
test_tripleclouds:
	$(CHANGENAM) $(CONFIG) config_tc.nam \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"' \
		overlap_scheme_name='"Exp-Ran"'
	$(DRIVER) config_tc.nam $(INPUT).nc $(INPUT)_tc_out.nc

# Longwave scattering; since 46R1 this is the default
test_lwscat:
	$(CHANGENAM) $(CONFIG) config_lwscat.nam \
		do_lw_cloud_scattering="true"
	$(DRIVER) config_lwscat.nam $(INPUT).nc $(INPUT)_lwscat_out.nc

# 3D radiative transfer
test_spartacus:
	$(CHANGENAM) $(CONFIG) config_spartacus.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" overlap_scheme_name='"Exp-Ran"' \
		do_sw_delta_scaling_with_gases="false"
	$(DRIVER) config_spartacus.nam $(INPUT).nc $(INPUT)_spartacus_out.nc

# 3D radiative transfer using the older "maximum entrapment"
test_spartacus_maxentr:
	$(CHANGENAM) $(CONFIG) config_spartacus_maxentr.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" overlap_scheme_name='"Exp-Ran"' \
		sw_entrapment_name='"Maximum"' \
		do_sw_delta_scaling_with_gases="false"
	$(DRIVER) config_spartacus_maxentr.nam $(INPUT).nc $(INPUT)_spartacus_maxentr_out.nc

# "Cloudless" solver
test_cloudless:
	$(CHANGENAM) $(CONFIG) config_cloudless.nam \
		use_aerosols=false \
		sw_solver_name='"Cloudless"' lw_solver_name='"Cloudless"'
	$(DRIVER) config_cloudless.nam $(INPUT).nc $(INPUT)_cloudless_out.nc

# Exponential-random overlap with "vectorizable" cloud generator
test_expran_vec:
	$(CHANGENAM) $(CONFIG) config_expran_vec.nam \
		overlap_scheme_name='"Exp-Ran"' use_vectorizable_generator=true
	$(DRIVER) config_expran_vec.nam $(INPUT).nc $(INPUT)_expran_vec_out.nc


### The following targets use the $CONFIG_ECCKD configuration file ###

# ecCKD gas optics scheme (note that default solver is Tripleclouds)
test_ecckd_mcica:
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_mcica.nam \
		sw_solver_name='"McICA"' lw_solver_name='"McICA"'
	$(DRIVER) config_ecckd_mcica.nam $(INPUT).nc $(INPUT)_ecckd_mcica_out.nc

# ecCKD with Tripleclouds solver (default)
test_ecckd_tc:
	$(DRIVER) $(CONFIG_ECCKD) $(INPUT).nc $(INPUT)_ecckd_tc_out.nc

# ecCKD with no aerosols
test_ecckd_noaer:
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_noaer.nam \
		use_aerosols=false
	$(DRIVER) config_ecckd_noaer.nam $(INPUT).nc $(INPUT)_ecckd_noaer_out.nc

# Test the different ways that aerosol optical properties can be
# averaged, outputing the aerosol properties in each gas-optics
# spectral interval, producing the following:
#   aerosol_optics_rrtmg.nc: RRTMG gas optics, aerosol properties from band-wise file
#   aerosol_optics_ecckd.nc: ecCKD gas optics, aerosol properties from high-res file
#   aerosol_optics_rrtmg.nc: RRTMG gas optics, aerosol properties from high-res file
#   aerosol_optics_ecckd.nc: ecCKD gas optics, aerosol properties from band-wise file
test_aerosol_averaging:
	$(CHANGENAM) $(CONFIG) config_rrtmg_saveaer.nam \
		do_save_aerosol_optics=true
	$(DRIVER) config_rrtmg_saveaer.nam $(INPUT).nc $(INPUT)_rrtmg_saveaer_out.nc
	mv aerosol_optics.nc aerosol_optics_rrtmg.nc
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_saveaer.nam \
		do_save_aerosol_optics=true
	$(DRIVER) config_ecckd_saveaer.nam $(INPUT).nc $(INPUT)_ecckd_saveaer_out.nc
	mv aerosol_optics.nc aerosol_optics_ecckd.nc
	$(CHANGENAM) $(CONFIG) config_rrtmg_gen_saveaer.nam \
		do_save_aerosol_optics=true use_general_aerosol_optics=true \
		aerosol_optics_override_file_name="'aerosol_ifs_48R1.nc'" 
	$(DRIVER) config_rrtmg_gen_saveaer.nam $(INPUT).nc $(INPUT)_rrtmg_gen_saveaer_out.nc
	mv aerosol_optics.nc aerosol_optics_rrtmg_gen.nc
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_band_saveaer.nam \
		do_save_aerosol_optics=true \
		aerosol_optics_override_file_name="'aerosol_ifs_rrtm_46R1_with_NI_AM.nc'" 
	$(DRIVER) config_ecckd_band_saveaer.nam $(INPUT).nc $(INPUT)_ecckd_band_saveaer_out.nc
	mv aerosol_optics.nc aerosol_optics_ecckd_band.nc

# ecCKD gas optics with SPARTACUS solver (not currently correct)
test_ecckd_spartacus:
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_spartacus.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" 
	$(DRIVER) config_ecckd_spartacus.nam $(INPUT).nc $(INPUT)_ecckd_spartacus_out.nc

# ecCKD gas optics with no aerosol
test_ecckd_tc_noaer:
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_tc_noaer.nam \
		use_aerosols=false 
	$(DRIVER) config_ecckd_tc_noaer.nam $(INPUT).nc $(INPUT)_ecckd_tc_noaer_out.nc

# Clean data files
clean:
	rm -f *_out.nc tmp*.nc radiative_properties*.nc inputs.nc \
	config_*.nam gmon.out
