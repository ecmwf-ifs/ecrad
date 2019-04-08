INPUT = ecrad_meridian

#DRIVER = totalview ../../bin/ecrad -a
DRIVER = ../../bin/ecrad

#CONFIG = configCY43R3.nam
#CONFIG = configCY46R1.nam
CONFIG = configCY47R1.nam

# Typing "make" will run radiation scheme on IFS profiles
all: test

test: test_default test_noaer test_expran test_tripleclouds test_spartacus \
	test_spartacus_maxentr

# (a) Default IFS settings for cycle 43R3: McICA solver with
# exponential-exponential overlap
test_default:
	$(DRIVER) $(CONFIG) $(INPUT).nc $(INPUT)_default_out.nc

# (b) Turn off aerosols
test_noaer:
	./change_namelist.sh $(CONFIG) config_noaer.nam \
		use_aerosols=false
	$(DRIVER) config_noaer.nam $(INPUT).nc $(INPUT)_noaer_out.nc

# (c) Exponential-random overlap
test_expran:
	./change_namelist.sh $(CONFIG) config_expran.nam \
		overlap_scheme_name='"Exp-Ran"'
	$(DRIVER) config_expran.nam $(INPUT).nc $(INPUT)_expran_out.nc

# (d) Tripleclouds solver with exponential-random overlap 
test_tripleclouds:
	./change_namelist.sh $(CONFIG) config_tc.nam \
		sw_solver_name='"Tripleclouds"' lw_solver_name='"Tripleclouds"' \
		overlap_scheme_name='"Exp-Ran"'
	$(DRIVER) config_tc.nam $(INPUT).nc $(INPUT)_tc_out.nc

# (e) Longwave scattering; since 46R1 this is the default
test_lwscat:
	./change_namelist.sh $(CONFIG) config_lwscat.nam \
		do_lw_cloud_scattering="true"
	$(DRIVER) config_lwscat.nam $(INPUT).nc $(INPUT)_lwscat_out.nc

# (f) 3D radiative transfer
test_spartacus:
	./change_namelist.sh $(CONFIG) config_spartacus.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" overlap_scheme_name='"Exp-Ran"' \
		do_sw_delta_scaling_with_gases="false"
	$(DRIVER) config_spartacus.nam $(INPUT).nc $(INPUT)_spartacus_out.nc

# (g) 3D radiative transfer using the older "maximum entrapment"
test_spartacus_maxentr:
	./change_namelist.sh $(CONFIG) config_spartacus_maxentr.nam \
		sw_solver_name='"SPARTACUS"' lw_solver_name='"SPARTACUS"' \
		do_3d_effects="true" overlap_scheme_name='"Exp-Ran"' \
		sw_entrapment_name='"Maximum"' \
		do_sw_delta_scaling_with_gases="false"
	$(DRIVER) config_spartacus_maxentr.nam $(INPUT).nc $(INPUT)_spartacus_maxentr_out.nc

# Clean data files
clean:
	rm -f *_out.nc tmp*.nc radiative_properties*.nc inputs.nc \
	config_*.nam gmon.out
