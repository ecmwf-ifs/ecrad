all: compare

compare: control experiment
	python3 ../common/nccmp.py control.nc experiment.nc

clean:
	rm -f control.nc experiment.nc

control:
	ecrad config.nam era5slice.nc control.nc

experiment:
	# Run ecrad
	driver.py config.nam era5slice.nc experiment.nc
	# Add the history attribute for the nccmp.py script
	python3 -c "import netCDF4; \
	            nc = netCDF4.Dataset('experiment.nc', 'a'); \
	            nc.history='pyecrad test'; \
	            nc.close()"
