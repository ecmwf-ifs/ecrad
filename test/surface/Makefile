INPUT = mls_london1

#DEBUG = totalview
#DRIVER = $(DEBUG) ../../bin/ecrad -a
DRIVER = ../../bin/ecrad

CONFIG = configCY43R3.nam

DUPPROF = ../i3rc/duplicate_profiles.sh
DUPPROFSURF = ./duplicate_profiles_surf.sh

# Typing "make" will run radiation scheme on IFS profiles
all: test

test: test_default test_sza test_bheight test_bheight_vac

# Single profile
test_default:
	$(DRIVER) $(CONFIG) $(INPUT).nc $(INPUT)_default_out.nc $(INPUT)_default_surf_out.nc

# Multiple profiles with varying solar zenith angle
test_sza: $(INPUT)_sza.nc
	$(DRIVER) $(CONFIG) $(INPUT)_sza.nc $(INPUT)_sza_out.nc $(INPUT)_sza_surf_out.nc

# Multiple profiles with varying building height
test_bheight: $(INPUT)_bheight.nc
	$(DRIVER) $(CONFIG) $(INPUT)_bheight.nc $(INPUT)_bheight_out.nc $(INPUT)_bheight_surf_out.nc

test_bheight_vac: $(INPUT)_bheight.nc
	./change_namelist.sh $(CONFIG) config_vac.nam \
		do_canopy_gases_lw="false" do_canopy_gases_sw="false"
	$(DRIVER) config_vac.nam $(INPUT)_bheight.nc $(INPUT)_bheight_vac_out.nc $(INPUT)_bheight_surf_vac_out.nc

$(INPUT)_sza.nc: $(INPUT).nc $(DUPPROF)
	$(DUPPROF) $(INPUT).nc $(INPUT)_sza.nc

$(INPUT)_bheight.nc: $(INPUT).nc $(DUPPROFSURF)
	$(DUPPROFSURF) $(INPUT).nc $(INPUT)_bheight.nc

# Clean data files
clean:
	rm -f $(INPUT)_*_out.nc tmp*.nc radiative_properties*.nc inputs.nc \
	config_*.nam $(INPUT)_sza.nc $(INPUT)_bheight.nc gmon.out
