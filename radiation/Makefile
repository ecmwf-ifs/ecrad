SOURCES = radiation_aerosol.F90 radiation_config.F90 \
	radiation_flux.F90 	radiation_cloud.F90 \
	radiation_thermodynamics.F90 radiation_lw_derivatives.F90 \
	radiation_gas.F90	radiation_single_level.F90 \
	radiation_cloud_optics_data.F90 \
	radiation_interface.F90 radiation_cloud_optics.F90 \
	radiation_overlap.F90 	radiation_two_stream.F90 \
	radiation_save.F90 	radiation_monochromatic.F90 \
	radiation_constants.F90 radiation_matrix.F90 \
	radiation_spartacus_sw.F90 radiation_cloud_cover.F90 \
	radiation_spartacus_lw.F90 radiation_pdf_sampler.F90 \
	radiation_aerosol_optics.F90 radiation_aerosol_optics_data.F90 \
	radiation_ifs_rrtm.F90 	radiation_adding_ica_sw.F90 \
	radiation_mcica_sw.F90	radiation_cloud_generator.F90 \
	radiation_mcica_lw.F90	radiation_adding_ica_lw.F90 \
	radiation_homogeneous_sw.F90 radiation_homogeneous_lw.F90 \
	radiation_ice_optics_fu.F90 radiation_ice_optics_baran.F90 \
	radiation_ice_optics_baran2017.F90 radiation_ice_optics_yi.F90 \
	radiation_liquid_optics_socrates.F90 radiation_liquid_optics_slingo.F90 \
	radiation_tripleclouds_sw.F90 radiation_tripleclouds_lw.F90 \
	radiation_regions.F90 radiation_ecckd.F90 \
	radiation_cloudless_lw.F90 radiation_cloudless_sw.F90 \
	radiation_check.F90   radiation_ecckd_interface.F90 \
	radiation_gas_constants.F90 radiation_ecckd_gas.F90 \
	radiation_spectral_definition.F90 radiation_general_cloud_optics.F90 \
	radiation_general_cloud_optics_data.F90 radiation_random_numbers.F90 \
	radiation_aerosol_optics_description.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBRAD = ../lib/libradiation.a

all: $(LIBRAD)

$(LIBRAD): $(OBJECTS)
	ar r $(LIBRAD) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBRAD)

radiation_flux.o radiation_aerosol.o: radiation_config.o
radiation_interface.o: radiation_aerosol.o radiation_config.o \
	radiation_flux.o radiation_cloud.o radiation_save.o \
	radiation_single_level.o radiation_gas.o radiation_monochromatic.o \
	radiation_thermodynamics.o radiation_spartacus_sw.o \
	radiation_spartacus_lw.o radiation_cloud_optics.o \
	radiation_aerosol_optics.o radiation_config.o \
	radiation_ifs_rrtm.o radiation_mcica_sw.o radiation_mcica_lw.o \
	radiation_homogeneous_lw.o radiation_homogeneous_sw.o \
	radiation_cloudless_lw.o radiation_cloudless_sw.o \
	radiation_tripleclouds_sw.o radiation_tripleclouds_lw.o \
	radiation_ecckd_interface.o radiation_general_cloud_optics.o
radiation_ifs_rrtm.o radiation_monochromatic.o: radiation_config.o \
	radiation_thermodynamics.o radiation_gas.o radiation_single_level.o
radiation_spartacus_sw.o radiation_spartacus_lw.o radiation_mcica_sw.o \
	radiation_mcica_lw.o radiation_homogeneous_sw.o radiation_homogeneous_lw.o \
	radiation_cloudless_sw.o radiation_cloudless_lw.o: \
	radiation_config.o \
	radiation_thermodynamics.o \
	radiation_single_level.o radiation_cloud.o \
	radiation_two_stream.o
radiation_cloud.o: radiation_thermodynamics.o
radiation_save.o: radiation_config.o \
	radiation_single_level.o radiation_thermodynamics.o \
	radiation_cloud.o
radiation_single_level.o: radiation_config.o
radiation_monochromatic.o radiation_aerosol_optics.o radiation_cloud_optics.o \
	radiation_homogeneous_sw.o radiation_spartacus_lw.o \
	radiation_cloudless_sw.o radiation_cloudless_lw.o \
	radiation_ecckd.o radiation_flux.o \
	radiation_spartacus_sw.o radiation_thermodynamics.o: radiation_constants.o
radiation_cloud_optics.o radiation_config.o: \
	radiation_cloud_optics_data.o
radiation_cloud_optics.o: radiation_ice_optics_baran.o radiation_ice_optics_fu.o \
	radiation_liquid_optics_socrates.o radiation_liquid_optics_slingo.o \
	radiation_ice_optics_baran2017.o radiation_ice_optics_yi.o
radiation_aerosol_optics.o radiation_config.o: \
	radiation_aerosol_optics_data.o
radiation_mcica_sw.o radiation_homogeneous_sw.o: radiation_adding_ica_sw.o radiation_cloud_generator.o
radiation_mcica_lw.o radiation_homogeneous_lw.o: radiation_adding_ica_lw.o radiation_cloud_generator.o \
	radiation_lw_derivatives.o
radiation_cloudless_sw.o: radiation_adding_ica_sw.o
radiation_cloudless_lw.o radiation_tripleclouds_lw.o: radiation_adding_ica_lw.o radiation_lw_derivatives.o
radiation_aerosol_optics.o radiation_cloud_optics.o: \
	radiation_config.o
radiation_cloud_generator.o radiation_config.o:	radiation_pdf_sampler.o \
	radiation_cloud_cover.o radiation_random_numbers.o
radiation_lw_derivatives.o: radiation_matrix.o

radiation_config.o: radiation_cloud_generator.o radiation_ecckd.o radiation_general_cloud_optics_data.o

radiation_aerosol_optics.F90 radiation_cloud_optics.F90 radiation_homogeneous_sw.F90 \
	radiation_mcica_sw.F90 radiation_cloudless_sw.o: radiation_delta_eddington.h
radiation_spartacus_sw.o radiation_spartacus_lw.o \
	radiation_tripleclouds_sw.o radiation_tripleclouds_lw.o: \
	radiation_optical_depth_scaling.h radiation_overlap.o \
	radiation_matrix.o radiation_flux.o radiation_regions.o
radiation_ecckd.o: radiation_ecckd_gas.o radiation_spectral_definition.o
radiation_ecckd_gas.o: radiation_gas_constants.o
radiation_ecckd_interface.o: radiation_config.o radiation_single_level.o radiation_thermodynamics.o \
	radiation_gas.o radiation_gas_constants.o

radiation_gas.o radiation_single_level.o radiation_thermodynamics.o \
	radiation_flux.o radiation_aerosol.o: radiation_check.o

radiation_gas.o: radiation_gas_constants.o

radiation_general_cloud_optics_data.o: radiation_spectral_definition.o

radiation_general_cloud_optics.o: radiation_config.o radiation_cloud.o radiation_thermodynamics.o radiation_constants.o

*.o: ecrad_config.h
