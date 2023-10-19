# Makefile for running ecRad on clear-sky CKDMIP profiles. The basic
# test uses the single input file INPUT below. You can also run ecRad
# on CKDMIP scenarios (see lower in file). You need the NCO tools
# ncap2, ncrcat and ncpdq.

# --- BASIC TEST ---

# Type "make" to use the default ecckd gas-optics models, or "make
# GASMODEL=rrtmg" for the RRTMG gas-optics model.

INPUT     = ckdmip_evaluation1_concentrations_present_reduced.nc
DRIVER    = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh
GASMODEL  = ecckd
CONFIG    = config-$(GASMODEL).nam

all: test

# By default run the sw and lw tests
test: sw lw

# Longwave ecRad calculation: simply run ecRad with an appropriate
# config file
lw: config-$(GASMODEL)_lw.nam
	$(DRIVER) config-$(GASMODEL)_lw.nam $(INPUT) ckdmip_evaluation1_lw_out.nc

# Shortwave ecRad calculation: run the five solar zenith angles then
# concatenate the NetCDF files. This uses the GNU make feature "call",
# here applied to CALCSWSZA which takes two arguments, the input and
# output files
sw: config_sw_all
	$(call CALCSWSZA,$(INPUT),ckdmip_evaluation1_sw_out.nc)

# Create a longwave config file (basically turn of the shortwave)
config-$(GASMODEL)_lw.nam: $(CONFIG)
	$(CHANGENAM) $(CONFIG) config-$(GASMODEL)_lw.nam do_sw=false

# Create the namelist files for a particular solar zenith angle,
# e.g. "make config-rrtmg_sw_0.5.nam" creates a namelist file for the
# RRTMG gas optics model and a solar zenith angle of 60 degrees.
config-$(GASMODEL)_sw_%.nam: $(CONFIG)
	$(CHANGENAM) $(CONFIG) $@ do_lw=false cos_solar_zenith_angle=$*

# Produce config files for each of the five solar zenith angles
config_sw_all: config-$(GASMODEL)_sw_0.1.nam config-$(GASMODEL)_sw_0.3.nam \
	config-$(GASMODEL)_sw_0.5.nam config-$(GASMODEL)_sw_0.7.nam config-$(GASMODEL)_sw_0.9.nam


# --- HELPER FUNCTIONS ---

# Function to run ecRad shortwave-only on the 50 evaluation profiles
# for a single value of the solar zenith angle (works only with GNU
# make), and then to add a mu0 variable storing the cosine of the
# solar zenith angle. The arguments are (1) the input file name and
# (2) the cosine of the solar zenith angle.
CALCSW = \
	$(DRIVER) config-$(GASMODEL)_sw_$(2).nam $(1) tmp_out.nc; \
	ncap2 -O -s 'defdim("mu0",1,0);mu0[mu0]={'$(2)'f};mu0@long_name="Cosine of solar zenith angle";flux_dn_sw[mu0,column,half_level]=flux_dn_sw;flux_up_sw[mu0,column,half_level]=flux_up_sw;flux_dn_direct_sw[mu0,column,half_level]=flux_dn_direct_sw' tmp_out.nc tmp_$(2)_out.nc

# Function to run ecRad shortwave-only for all five solar zenith
# angles, then to contatenate the results into a single file. The
# arguments are (1) the input file name and (2) the output file name.
CALCSWSZA = \
	$(call CALCSW,$(1),0.1); \
	$(call CALCSW,$(1),0.3); \
	$(call CALCSW,$(1),0.5); \
	$(call CALCSW,$(1),0.7); \
	$(call CALCSW,$(1),0.9); \
	ncrcat -O tmp_0.1_out.nc tmp_0.3_out.nc tmp_0.5_out.nc tmp_0.7_out.nc tmp_0.9_out.nc tmp_out2.nc; \
	ncpdq -O -a column,mu0,half_level tmp_out2.nc $(2)


# --- CKDMIP SCENARIOS ---

# To run the CKDMIP scenarios you first need to copy the concentration
# files from http://aux.ecmwf.int/ecpds/home/ckdmip/concentrations/,
# or ftp from aux.ecmwf.int (username ckdmip, password any non-empty
# string), which you can do with "make download".  Then do:
#
#    make scenarios
#
# As in the basic test you can add the GASMODEL=rrtmg to use the RRTMG
# gas-optics model instead of ecCKD. To use a different configuration
# of ecCKD than the default, copy config-ecckd.nam to something like
# config-custom.nam, then add the gas_optics_sw_override_file_name and
# gas_optics_lw_override_file_name namelist arguments to specify your
# preferred CKD definition file, and call make with GASMODEL=custom.

CONCDIR    = concentrations
DATASET    = evaluation1
FLUXDIR    = fluxes
LBLFLUXDIR = lbl_fluxes

# Web/FTP location of concentration files (no trailing slash)
#WEBSOURCE = http://aux.ecmwf.int/ecpds/home/ckdmip
WEBSOURCE = ftp://ckdmip:blank@aux.ecmwf.int

# List of scenarios to run
SCENARIOS = present co2-180 co2-280 co2-560 co2-1120 co2-2240 ch4-350 ch4-700 ch4-1200 ch4-2600 ch4-3500 n2o-190 n2o-270 n2o-405 n2o-540

# A longwave and shortwave make rule for each scenario
SCENARIO_RULES = $(addprefix lw_scenario_,$(SCENARIOS)) $(addprefix sw_scenario_,$(SCENARIOS))
# A download rule for each scenario
DOWNLOAD_RULES = $(addprefix download_,$(SCENARIOS))

# Main target of this section runs each scenario longwave and shortwave
scenarios: $(SCENARIO_RULES)

# Pattern matches a shortwave scenario, invoked by calling the
# CALCSWSZA function.
sw_scenario_%: config_sw_all
	mkdir -p $(FLUXDIR)
	$(call CALCSWSZA,$(CONCDIR)/ckdmip_$(DATASET)_concentrations_$*.nc,$(FLUXDIR)/ecrad-$(GASMODEL)_$(DATASET)_sw_fluxes_$*.nc)

# Pattern matches a longwave scenario.
lw_scenario_%: config-$(GASMODEL)_lw.nam
	mkdir -p $(FLUXDIR)
	$(DRIVER) config-$(GASMODEL)_lw.nam $(CONCDIR)/ckdmip_$(DATASET)_concentrations_$*.nc $(FLUXDIR)/ecrad-$(GASMODEL)_$(DATASET)_lw_fluxes_$*.nc

# Target to download concentration files from the web to the CONCDIR
# directory
download: $(DOWNLOAD_RULES)

# Download one scenario, explicitly overwriting if it already
# exists. Note that the "h5" files are NetCDF4/HDF5 format so can be
# read with NetCDF tools
download_%:
	mkdir -p $(CONCDIR)
	wget -nv -O $(CONCDIR)/ckdmip_$(DATASET)_concentrations_$*.nc \
		$(WEBSOURCE)/concentrations/ckdmip_$(DATASET)_concentrations_$*.nc
	mkdir -p $(LBLFLUXDIR)
	wget -nv -O $(LBLFLUXDIR)/ckdmip_$(DATASET)_lw_fluxes_$*.nc \
		$(WEBSOURCE)/lw_fluxes/$(DATASET)/ckdmip_$(DATASET)_lw_fluxes_$*.h5
	wget -nv -O $(LBLFLUXDIR)/ckdmip_$(DATASET)_sw_fluxes_$*.nc \
		$(WEBSOURCE)/sw_fluxes/$(DATASET)/ckdmip_$(DATASET)_sw_fluxes_$*.h5


# Clean data files
clean:
	rm -f *_out.nc config-*_sw_*.nam config-*_lw.nam tmp*.nc wget-log

distclean: clean
	rm -rf $(CONCDIR)/ckdmip_$(DATASET)_concentrations_present.nc \
		$(CONCDIR)/ckdmip_$(DATASET)_*-*.nc \
		$(FLUXDIR)/ecrad-*_$(DATASET)_lw_fluxes_*.nc \
		$(FLUXDIR)/ecrad-*_$(DATASET)_sw_fluxes_*.nc \
		$(LBLFLUXDIR)/ckdmip_$(DATASET)_lw_fluxes_*.nc \
		$(LBLFLUXDIR)/ckdmip_$(DATASET)_sw_fluxes_*.nc

.PHONY: config_sw_all scenarios sw lw all test clean distclean
