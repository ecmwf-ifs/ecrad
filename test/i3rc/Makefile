INPUT_PROFILE = i3rc_mls_cumulus

DRIVER = ../../bin/ecrad

CONFIG = configI3RC.nam

# Typing "make" will run radiation scheme in all possible ways
all: test

test: $(INPUT_PROFILE)_sza.nc i3rc_spartacus i3rc_tripleclouds i3rc_mcica


# Create netcdf file containing duplicated profile but with varying
# solar zenith angle
$(INPUT_PROFILE)_sza.nc:  $(INPUT_PROFILE).nc duplicate_profiles.sh
	./duplicate_profiles.sh $(INPUT_PROFILE).nc $(INPUT_PROFILE)_sza.nc

i3rc_spartacus: $(INPUT_PROFILE)_sza.nc
	./change_namelist.sh $(CONFIG) config_3reg_3d.nam \
		n_regions=3 do_3d_effects=true do_3d_lw_multilayer_effects=true
	./change_namelist.sh $(CONFIG) config_3reg_1d.nam \
		n_regions=3 do_3d_effects=false do_3d_lw_multilayer_effects=false
	./change_namelist.sh $(CONFIG) config_3reg_3d_clustering.nam \
		n_regions=3 do_3d_effects=true effective_size_scaling=1.449
	$(DRIVER) config_3reg_1d.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_1d_out.nc
	$(DRIVER) config_3reg_3d.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_out.nc
	$(DRIVER) config_3reg_3d_clustering.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_clustering_out.nc
i3rc_spartacus2: $(INPUT_PROFILE)_sza.nc
	./change_namelist.sh $(CONFIG) config_2reg_3d.nam \
		n_regions=2 do_3d_effects=true do_3d_lw_multilayer_effects=true
	$(DRIVER) config_2reg_3d.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_2reg_3d_out.nc
i3rc_spartacus_extra:  $(INPUT_PROFILE)_sza.nc
	./change_namelist.sh $(CONFIG) config_3reg_3d_computed.nam \
		sw_encroachment_name='"Computed"'
	./change_namelist.sh $(CONFIG) config_3reg_3d_fractal.nam \
		sw_encroachment_name='"Fractal"'
	./change_namelist.sh $(CONFIG) config_3reg_3d_minimum.nam \
		sw_encroachment_name='"Minimum"'
	./change_namelist.sh $(CONFIG) config_3reg_1d_computed.nam \
		do_3d_effects=false sw_encroachment_name='"Computed"'
	./change_namelist.sh $(CONFIG) config_3reg_1d_minimum.nam \
		do_3d_effects=false sw_encroachment_name='"Minimum"'
	./change_namelist.sh $(CONFIG) config_3reg_3d_computedleast.nam \
		sw_encroachment_name='"Computed"' encroachment_scaling=0.0
	$(DRIVER) config_3reg_3d_computed.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_computed_out.nc
	$(DRIVER) config_3reg_3d_fractal.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_fractal_out.nc
	$(DRIVER) config_3reg_3d_minimum.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_minimum_out.nc
	$(DRIVER) config_3reg_1d_computed.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_1d_computed_out.nc
	$(DRIVER) config_3reg_1d_minimum.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_1d_minimum_out.nc
	$(DRIVER) config_3reg_3d_computedleast.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_computedleast_out.nc

i3rc_tripleclouds: $(INPUT_PROFILE)_sza.nc
	./change_namelist.sh $(CONFIG) config_tc.nam n_regions=3 \
		do_3d_effects=false sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"'
	$(DRIVER) config_tc.nam $(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_tc_out.nc

i3rc_mcica: $(INPUT_PROFILE)_sza.nc
	./change_namelist.sh $(CONFIG) config_mcica.nam \
		sw_solver_name='"McICA"' lw_solver_name='"McICA"'
	$(DRIVER) config_mcica.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_mcica_out.nc

i3rc_print_encroachment: $(INPUT_PROFILE)_sza.nc
	./change_namelist.sh $(CONFIG) config_print_encr.nam \
		sw_encroachment_name='"Computed"' \
		gas_model_name='"Monochromatic"' \
		sw_albedo=0.2
	$(DRIVER) config_print_encr.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_mono_encr_out.nc

# Clean data files
clean:
	rm -f $(INPUT_PROFILE)*_out.nc tmp*.nc radiative_properties*.nc \
	$(INPUT_PROFILE)*_sza.nc inputs.nc config_*.nam gmon.out
