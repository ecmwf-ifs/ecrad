PROFILE = tro_cirrus
PROFILE = mls_liquid

DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh
CONFIG = configDISORT.nam
CONFIG_TMP = config_tmp.nam

test: all

all: disort_ref disort_2stream disort_4stream disort_4stream_jac disort_elsasser ecrad

disort_ref:
	$(CHANGENAM) $(CONFIG) $(CONFIG_TMP) n_angles_per_hemisphere_lw=8
	$(DRIVER) $(CONFIG_TMP) $(PROFILE).nc $(PROFILE)_disort_ref_out.nc \
		|tee $(PROFILE)_disort_ref_out.log
	mv radiative_properties.nc radiative_properties_ref.nc

disort_2stream:
	$(DRIVER) $(CONFIG) $(PROFILE).nc $(PROFILE)_disort_2stream_out.nc \
		|tee $(PROFILE)_disort_2stream_out.log

disort_4stream:
	$(CHANGENAM) $(CONFIG) $(CONFIG_TMP) n_angles_per_hemisphere_lw=2
	$(DRIVER) $(CONFIG_TMP) $(PROFILE).nc $(PROFILE)_disort_4stream_out.nc \
		|tee $(PROFILE)_disort_4stream_out.log

disort_4stream_jac:
	$(CHANGENAM) $(CONFIG) $(CONFIG_TMP) n_angles_per_hemisphere_lw=2 lw_quadrature_name='"Jacobi"'
	$(DRIVER) $(CONFIG_TMP) $(PROFILE).nc $(PROFILE)_disort_4stream_jac_out.nc \
		|tee $(PROFILE)_disort_4stream_jac_out.log

disort_elsasser:
	$(CHANGENAM) $(CONFIG) $(CONFIG_TMP) lw_quadrature_name='"Elsasser"'
	$(DRIVER) $(CONFIG_TMP) $(PROFILE).nc $(PROFILE)_disort_elsasser_out.nc \
		|tee $(PROFILE)_disort_elsasser_out.log

ecrad:
	$(CHANGENAM) $(CONFIG) $(CONFIG_TMP) lw_solver_name='"Homogeneous"'
	$(DRIVER) $(CONFIG_TMP) $(PROFILE).nc $(PROFILE)_ecrad_out.nc \
		|tee $(PROFILE)_ecrad_out.log

clean:
	rm -f $(PROFILE)*_out.nc $(PROFILE)*_out.log $(CONFIG_TMP)
