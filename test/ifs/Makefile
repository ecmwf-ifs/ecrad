INPUT = ecrad_meridian

#DRIVER = totalview ../../bin/ecrad -a
DRIVER = ../../bin/ecrad
IFS_DRIVER = ../../bin/ecrad_ifs
IFS_DRIVER_BLOCKED = ../../bin/ecrad_ifs_blocked
CHANGENAM = ../common/change_namelist.sh

# Various configuration files corresponding to cycles of ECMWF's
# Integrated Forecasting System (IFS)
#CONFIG = configCY43R3.nam
#CONFIG = configCY46R1.nam
#CONFIG = configCY47R1.nam
#CONFIG = configCY47R3.nam
CONFIG = configCY49R1.nam
CONFIG_ECCKD = configCY49R1_ecckd.nam
CONFIG_MIXED = configCY49R1_mixed.nam

# Typing "make" will run radiation scheme on IFS profiles
all: test

test: test_orig test_ecckd

# Tests of CY47R1 configuration with RRTMG gas optics model and
# various solvers
test_orig: test_default test_noaer test_expexp test_vec test_tripleclouds test_spartacus \
	test_spartacus_maxentr

# Tests of ecCKD gas optics model and various solvers
test_ecckd: test_ecckd_mcica test_ecckd_tc test_ecckd_tc_noaer test_ecckd_spartacus

# Default IFS settings: McICA solver with exponential-exponential
# overlap
test_default:
	$(DRIVER) $(CONFIG) $(INPUT).nc $(INPUT)_default_out.nc

test_ifsdriver:
	$(CHANGENAM) $(CONFIG) config_net.nam do_save_net_fluxes=true do_write_double_precision=true \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"'
	$(IFS_DRIVER) config_net.nam $(INPUT).nc $(INPUT)_ifsdriver_out.nc | tee $(INPUT)_ifsdriver_out.log
	$(DRIVER) config_net.nam $(INPUT).nc $(INPUT)_net_out.nc | tee $(INPUT)_net_out.log

test_ifsdriver_blocked:
	$(CHANGENAM) $(CONFIG) config_net.nam do_save_net_fluxes=true do_write_double_precision=true \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"'
	$(IFS_DRIVER_BLOCKED) config_net.nam $(INPUT).nc $(INPUT)_ifsdriver_blocked_out.nc | tee $(INPUT)_ifsdriver_blocked_out.log
	$(DRIVER) config_net.nam $(INPUT).nc $(INPUT)_net_out.nc | tee $(INPUT)_net_out.log

# Turn off aerosols
test_noaer:
	$(CHANGENAM) $(CONFIG) config_noaer.nam \
		use_aerosols=false
	$(DRIVER) config_noaer.nam $(INPUT).nc $(INPUT)_noaer_out.nc

# Older exponential-exponential overlap
test_expexp:
	$(CHANGENAM) $(CONFIG) config_expexp.nam \
		overlap_scheme_name='"Exp-Exp"'
	$(DRIVER) config_expexp.nam $(INPUT).nc $(INPUT)_expexp_out.nc

# Tripleclouds solver with exponential-random overlap 
test_tripleclouds:
	$(CHANGENAM) $(CONFIG) config_tc.nam \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"'
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
		do_3d_effects="true" \
		do_sw_delta_scaling_with_gases="false"
	$(DRIVER) config_spartacus.nam $(INPUT).nc $(INPUT)_spartacus_out.nc

# 3D radiative transfer using the older "maximum entrapment"
test_spartacus_maxentr:
	$(CHANGENAM) $(CONFIG) config_spartacus_maxentr.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" \
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
test_vec:
	$(CHANGENAM) $(CONFIG) config_vec.nam use_vectorizable_generator=true
	$(DRIVER) config_vec.nam $(INPUT).nc $(INPUT)_vec_out.nc


### The following targets use the $CONFIG_ECCKD configuration file ###

# ecCKD gas optics scheme (note that default solver is Tripleclouds)
test_ecckd_mcica:
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_mcica.nam \
		sw_solver_name='"McICA"' lw_solver_name='"McICA"'
	$(DRIVER) config_ecckd_mcica.nam $(INPUT).nc $(INPUT)_ecckd_mcica_out.nc

# ecCKD with Tripleclouds solver (default)
test_ecckd_tc:
	$(DRIVER) $(CONFIG_ECCKD) $(INPUT).nc $(INPUT)_ecckd_tc_out.nc

test_mixed_gas:
	$(DRIVER) $(CONFIG_MIXED) $(INPUT).nc $(INPUT)_sw_ecckd_lw_ecckd_out.nc
	$(CHANGENAM) $(CONFIG_MIXED) config_mix.nam lw_gas_model_name='"RRTMG-IFS"' do_cloud_aerosol_per_lw_g_point=false
	$(DRIVER) config_mix.nam $(INPUT).nc $(INPUT)_sw_ecckd_lw_rrtmg_out.nc
	$(CHANGENAM) $(CONFIG_MIXED) config_mix.nam sw_gas_model_name='"RRTMG-IFS"' do_cloud_aerosol_per_sw_g_point=false
	$(DRIVER) config_mix.nam $(INPUT).nc $(INPUT)_sw_rrtmg_lw_ecckd_out.nc
	$(CHANGENAM) $(CONFIG_MIXED) config_mix.nam sw_gas_model_name='"RRTMG-IFS"' lw_gas_model_name='"RRTMG-IFS"' do_cloud_aerosol_per_lw_g_point=false do_cloud_aerosol_per_sw_g_point=false
	$(DRIVER) config_mix.nam $(INPUT).nc $(INPUT)_sw_rrtmg_lw_rrtmg_out.nc

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

# Profiling
TAG = default
profile:
	mkdir -p profiling_$(TAG)
	$(CHANGENAM) $(CONFIG) config_profile.nam \
		nrepeat=100 do_save_spectral_flux=false
	DR_HOOK=1 DR_HOOK_OPT=prof DR_HOOK_PROFILE=profiling_$(TAG)/$(INPUT)_default_prof_out.drhook \
		$(DRIVER) config_profile.nam $(INPUT).nc $(INPUT)_default_prof_out.nc \
		| tee $(INPUT)_default_prof_out.log
	if [[ -f gmon.out ]]; then gprof $(DRIVER) gmon.out > profiling_$(TAG)/$(INPUT)_default_prof_out.gprof; fi
	$(CHANGENAM) $(CONFIG_ECCKD) config_profile.nam \
		nrepeat=100 do_save_spectral_flux=false do_save_radiative_properties=false
	DR_HOOK=1 DR_HOOK_OPT=prof DR_HOOK_PROFILE=profiling_$(TAG)/$(INPUT)_ecckd_tc_prof_out.drhook \
		$(DRIVER) config_profile.nam $(INPUT).nc $(INPUT)_ecckd_tc_prof_out.nc \
		| tee $(INPUT)_ecckd_tc_prof_out.log
	if [[ -f gmon.out ]]; then gprof $(DRIVER) gmon.out > $(INPUT)_ecckd_tc_prof_out.gprof; fi

	mv $(INPUT)_*_prof_out.* profiling_$(TAG)

# ecCKD gas optics with spectral diagnostics every 25nm in the visible
# written to a separate file; do this both with the 32- and 96-gpoint
# models (which have 100nm and 25nm resolution in the visible,
# respectively).
DIAGBOUNDS="sw_diag_wavelength_bound=.2e-6,.225e-6,.25e-6,.275e-6,.3e-6,.325e-6,.35e-6,.375e-6,.4e-6,.425e-6,.45e-6,.475e-6,.5e-6,.525e-6,.55e-6,.575e-6,.6e-6,.625e-6,.65e-6,.675e-6,.7e-6"
test_diag:
	$(CHANGENAM) $(CONFIG_ECCKD) config_diag.nam $(DIAGBOUNDS) sw_diag_file_name='"sw_diag_rrtmg.nc"' gas_model_name='"RRTMG-IFS"'
	$(DRIVER) config_diag.nam $(INPUT).nc $(INPUT)_diag_rrtmg_out.nc
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_diag.nam $(DIAGBOUNDS) sw_diag_file_name='"sw_diag_rgb_orig.nc"'
	$(DRIVER) config_ecckd_diag.nam $(INPUT).nc $(INPUT)_ecckd_diag_rgb_orig_out.nc
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_diag.nam $(DIAGBOUNDS) sw_diag_file_name='"sw_diag_rgb.nc"' \
		gas_optics_sw_override_file_name='"ecckd-1.4_sw_climate_rgb-32b_ckd-definition.nc"'
	$(DRIVER) config_ecckd_diag.nam $(INPUT).nc $(INPUT)_ecckd_diag_rgb_out.nc
	$(CHANGENAM) $(CONFIG_ECCKD) config_ecckd_diag.nam $(DIAGBOUNDS) sw_diag_file_name='"sw_diag_vfine.nc"' \
		gas_optics_sw_override_file_name='"ecckd-1.4_sw_climate_vfine-96b_ckd-definition.nc"'
	$(DRIVER) config_ecckd_diag.nam $(INPUT).nc $(INPUT)_ecckd_diag_vfine_out.nc

# Clean data files
clean:
	rm -f *_out.nc tmp*.nc radiative_properties*.nc inputs.nc sw_diag*nc \
	config_*.nam gmon.out
