INPUT = ckdmip_evaluation1_concentrations_present_reduced.nc
DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh
CONFIG = config.nam
MKRECDIM = ncks -O --mk_rec_dmn column 

# Function to run ecRad shortwave-only on the 50 evaluation profiles
# for a single value of the solar zenith angle
CALCSW = $(CHANGENAM) $(CONFIG) config_sw$(1).nam \
		do_sw=true cos_solar_zenith_angle=$(1); \
	$(DRIVER) config_sw$(1).nam $(INPUT) tmp_out.nc; \
	$(MKRECDIM) tmp_out.nc ckdmip_evaluation1_sw$(1)_out.nc

all: test

test: sw lw

# Shortwave ecRad calculation: run the five solar zenith angles then
# concatenate the NetCDF files
sw:
	$(call CALCSW,0.1)
	$(call CALCSW,0.3)
	$(call CALCSW,0.5)
	$(call CALCSW,0.7)
	$(call CALCSW,0.9)
	ncrcat -O ckdmip_evaluation1_sw0.1_out.nc \
		ckdmip_evaluation1_sw0.3_out.nc \
		ckdmip_evaluation1_sw0.5_out.nc \
		ckdmip_evaluation1_sw0.7_out.nc \
		ckdmip_evaluation1_sw0.9_out.nc \
		ckdmip_evaluation1_sw_out.nc \

# Longwave ecRad calculation
lw:
	$(CHANGENAM) $(CONFIG) config_lw.nam do_lw=true
	$(DRIVER) config_lw.nam $(INPUT) ckdmip_evaluation1_lw_out.nc

# Clean data files
clean:
	rm -f *_out.nc config_*.nam

.PHONY: check_nco
