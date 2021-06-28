INPUT = ecrad_meridian

#DRIVER = totalview ../../bin/ecrad -a
DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh

#CONFIG = configCY43R3.nam
#CONFIG = configCY46R1.nam
CONFIG = configCY47R1.nam

# Typing "make" will run radiation scheme on IFS profiles
all: test

test: test_default test_noaer test_expran test_expran_vec test_tripleclouds test_spartacus \
	test_spartacus_maxentr test_ecckd test_ecckd_noaer test_ecckd_spartacus

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

# ecCKD gas optics scheme
test_ecckd:
	$(CHANGENAM) $(CONFIG) config_ecckd.nam \
		gas_model_name='"ECCKD"' aerosol_optics_override_file_name='"aerosol_ifs_48R1.nc"' \
		do_nearest_spectral_lw_emiss=false \
		use_general_cloud_optics=true use_general_aerosol_optics=true
	$(DRIVER) config_ecckd.nam $(INPUT).nc $(INPUT)_ecckd_out.nc
#		sw_solver_name='"Homogeneous"' lw_solver_name='"Homogeneous"'

test_ecckd_tc:
	$(CHANGENAM) $(CONFIG) config_ecckd_tc.nam \
		gas_model_name='"ECCKD"' aerosol_optics_override_file_name='"aerosol_ifs_48R1.nc"' \
		 overlap_scheme_name='"Exp-Ran"' \
		do_nearest_spectral_lw_emiss=false \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"'
		use_general_cloud_optics=true use_general_aerosol_optics=true
	$(DRIVER) config_ecckd_tc.nam $(INPUT).nc $(INPUT)_ecckd_tc_out.nc

# ecCKD gas optics with SPARTACUS solver
test_ecckd_spartacus:
	$(CHANGENAM) $(CONFIG) config_ecckd_spartacus.nam \
		gas_model_name='"ECCKD"' aerosol_optics_override_file_name='"aerosol_ifs_48R1.nc"' \
		do_nearest_spectral_lw_emiss=false \
		use_general_cloud_optics=true use_general_aerosol_optics=true \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" overlap_scheme_name='"Exp-Ran"' \
		do_sw_delta_scaling_with_gases="false"
	$(DRIVER) config_ecckd_spartacus.nam $(INPUT).nc $(INPUT)_ecckd_spartacus_out.nc

# ecCKD gas optics with no aerosol
test_ecckd_noaer:
	$(CHANGENAM) $(CONFIG) config_ecckd_noaer.nam \
		gas_model_name='"ECCKD"' use_aerosols=false \
		do_nearest_spectral_lw_emiss=false \
		use_general_cloud_optics=true use_general_aerosol_optics=true
	$(DRIVER) config_ecckd_noaer.nam $(INPUT).nc $(INPUT)_ecckd_noaer_out.nc

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

# Clean data files
clean:
	rm -f *_out.nc tmp*.nc radiative_properties*.nc inputs.nc \
	config_*.nam gmon.out
