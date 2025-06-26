all: compare

compare: control experiment
	compare.py control.nc experiment.nc

clean:
	rm -f control.nc experiment.nc

control:
	ecrad config.nam era5slice.nc control.nc

experiment:
	driver.py config.nam era5slice.nc experiment.nc
