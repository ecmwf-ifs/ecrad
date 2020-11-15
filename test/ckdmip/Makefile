INPUT = ckdmip_evaluation1_concentrations_present_reduced.nc
DRIVER = ../../bin/ecrad
CHANGENAM = ../common/change_namelist.sh
CONFIG = config.nam

all: test

test: sw lw

sw:
	$(CHANGENAM) config.nam config_sw0.1.nam \
		do_sw=true cos_solar_zenith_angle=0.1
	$(DRIVER) config_sw0.1.nam $(INPUT) ckdmip_evaluation1_sw0.1_out.nc
	$(CHANGENAM) config.nam config_sw0.3.nam \
		do_sw=true cos_solar_zenith_angle=0.3
	$(DRIVER) config_sw0.3.nam $(INPUT) ckdmip_evaluation1_sw0.3_out.nc
	$(CHANGENAM) config.nam config_sw0.5.nam \
		do_sw=true cos_solar_zenith_angle=0.5
	$(DRIVER) config_sw0.5.nam $(INPUT) ckdmip_evaluation1_sw0.5_out.nc
	$(CHANGENAM) config.nam config_sw0.7.nam \
		do_sw=true cos_solar_zenith_angle=0.7
	$(DRIVER) config_sw0.7.nam $(INPUT) ckdmip_evaluation1_sw0.7_out.nc
	$(CHANGENAM) config.nam config_sw0.9.nam \
		do_sw=true cos_solar_zenith_angle=0.9
	$(DRIVER) config_sw0.9.nam $(INPUT) ckdmip_evaluation1_sw0.9_out.nc

lw:
	$(CHANGENAM) config.nam config_lw.nam \
		do_lw=true
	$(DRIVER) config_lw.nam $(INPUT) ckdmip_evaluation1_lw_out.nc

# Clean data files
clean:
	rm -f *_out.nc config_*.nam
