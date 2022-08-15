INPUT_PROFILE = i3rc_mls_cumulus

DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh
CONFIG = configI3RC.nam

# Typing "make" will run radiation scheme in all possible ways
all: test

test: $(INPUT_PROFILE)_sza.nc i3rc_spartacus i3rc_tripleclouds i3rc_mcica


# Create netcdf file containing duplicated profile but with varying
# solar zenith angle
$(INPUT_PROFILE)_sza.nc:  $(INPUT_PROFILE).nc duplicate_profiles.sh
	./duplicate_profiles.sh $(INPUT_PROFILE).nc $(INPUT_PROFILE)_sza.nc

i3rc_spartacus: $(INPUT_PROFILE)_sza.nc
	$(CHANGENAM) $(CONFIG) config_3reg_3d.nam \
		n_regions=3 do_3d_effects=true do_3d_lw_multilayer_effects=true
	$(CHANGENAM) $(CONFIG) config_3reg_1d.nam \
		n_regions=3 do_3d_effects=false do_3d_lw_multilayer_effects=false
	$(CHANGENAM) $(CONFIG) config_3reg_3d_clustering.nam \
		n_regions=3 do_3d_effects=true effective_size_scaling=1.449
	$(DRIVER) config_3reg_1d.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_1d_out.nc
	$(DRIVER) config_3reg_3d.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_out.nc
	$(DRIVER) config_3reg_3d_clustering.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_clustering_out.nc
i3rc_spartacus2: $(INPUT_PROFILE)_sza.nc
	$(CHANGENAM) $(CONFIG) config_2reg_3d.nam \
		n_regions=2 do_3d_effects=true do_3d_lw_multilayer_effects=true
	$(DRIVER) config_2reg_3d.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_2reg_3d_out.nc
i3rc_spartacus_extra:  $(INPUT_PROFILE)_sza.nc
	$(CHANGENAM) $(CONFIG) config_3reg_3d_explicit.nam \
		sw_entrapment_name='"Explicit"'
	$(CHANGENAM) $(CONFIG) config_3reg_3d_nonfractal.nam \
		sw_entrapment_name='"Non-fractal"'
	$(CHANGENAM) $(CONFIG) config_3reg_3d_edgeonly.nam \
		sw_entrapment_name='"Edge-only"'
	$(CHANGENAM) $(CONFIG) config_3reg_1d_explicit.nam \
		do_3d_effects=false sw_entrapment_name='"Explicit"'
	$(CHANGENAM) $(CONFIG) config_3reg_1d_edgeonly.nam \
		do_3d_effects=false sw_entrapment_name='"Edge-only"'
	$(CHANGENAM) $(CONFIG) config_3reg_3d_explicit_ohf1.nam \
		sw_entrapment_name='"Explicit"' overhang_factor=1.0
	$(DRIVER) config_3reg_3d_explicit.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_explicit_out.nc
	$(DRIVER) config_3reg_3d_nonfractal.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_nonfractal_out.nc	$(DRIVER) config_3reg_3d_edgeonly.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_edgeonly_out.nc
	$(DRIVER) config_3reg_1d_explicit.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_1d_explicit_out.nc
	$(DRIVER) config_3reg_1d_edgeonly.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_1d_edgeonly_out.nc
	$(DRIVER) config_3reg_3d_explicit_ohf1.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_3reg_3d_explicit_ohf1_out.nc

i3rc_tripleclouds: $(INPUT_PROFILE)_sza.nc
	$(CHANGENAM) $(CONFIG) config_tc.nam n_regions=3 \
		do_3d_effects=false sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"'
	$(DRIVER) config_tc.nam $(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_tc_out.nc

i3rc_mcica: $(INPUT_PROFILE)_sza.nc
	$(CHANGENAM) $(CONFIG) config_mcica.nam \
		sw_solver_name='"McICA"' lw_solver_name='"McICA"'
	$(DRIVER) config_mcica.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_mcica_out.nc

i3rc_print_entrapment: $(INPUT_PROFILE)_sza.nc
	$(CHANGENAM) $(CONFIG) config_print_entr.nam \
		sw_entrapment_name='"Explicit"' \
		gas_model_name='"Monochromatic"' \
		sw_albedo=0.2
	$(DRIVER) config_print_entr.nam \
		$(INPUT_PROFILE)_sza.nc $(INPUT_PROFILE)_mono_entr_out.nc

# Clean data files
clean:
	rm -f $(INPUT_PROFILE)*_out.nc tmp*.nc radiative_properties*.nc \
	$(INPUT_PROFILE)*_sza.nc inputs.nc config_*.nam gmon.out
