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
	radiation_mcica_lw.F90	radiation_cloud_generator_acc.F90	radiation_adding_ica_lw.F90 \
	radiation_mcica_acc_sw.F90	radiation_mcica_acc_lw.F90 \
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

include Makefile_deps
