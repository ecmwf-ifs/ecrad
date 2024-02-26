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

# Name of the inline library (a directory):
inlib_name := ecrad.il

ifdef GPU

# Additional compiler flags enabling the inline library usage:

# List of source files and functions that should be inlined
inlib_SOURCES:=
inlib_FUNCTIONS:=

inlib_SOURCES+=radiation_random_numbers.F90
inlib_FUNCTIONS+=initialize_acc uniform_distribution_acc

inlib_SOURCES+=radiation_lw_derivatives.F90
inlib_FUNCTIONS+=calc_lw_derivatives_ica2 modify_lw_derivatives_ica

inlib_SOURCES+=radiation_adding_ica_lw.F90
inlib_FUNCTIONS+=adding_ica_lw fast_adding_ica_lw calc_fluxes_no_scattering_lw

inlib_SOURCES+=radiation_adding_ica_sw.F90
inlib_FUNCTIONS+=adding_ica_sw

inlib_SOURCES+=radiation_cloud_generator_acc.F90
inlib_FUNCTIONS+=cloud_generator_acc

inlib_SOURCES+=radiation_two_stream.F90
inlib_FUNCTIONS+=calc_ref_trans_sw calc_ref_trans_lw calc_reflectance_transmittance_sw calc_no_scattering_transmittance_lw calc_two_stream_gammas_sw

, := ,
inlib_FUNCTIONS:=$(subst $(eval) ,$(,),$(foreach f,$(inlib_FUNCTIONS),name:$(f)))
inlib_extract_FCFLAGS:= -Mextract=lib:$(inlib_name),reshape,$(inlib_FUNCTIONS)
inlib_inline_FCFLAGS:= -Minline=lib:$(inlib_name),reshape,$(inlib_FUNCTIONS) -Mneginfo=inline

########################################

# Object files that correspond to the source files in $(inlib_SOURCES):
inlib_OBJS:= $(inlib_SOURCES:.F90=.o)

# List of objects that are generated using the inline library, but not part of the inline library
ifneq ($(strip $(inlib_OBJS)),)
consumer_objs:= $(shell grep $(foreach obj,$(inlib_OBJS),-e ":.*\s$(obj)\b") Makefile_deps | cut -d ':' -f 1 | sort | uniq)
consumer_objs:= $(filter-out $(inlib_OBJS),$(consumer_objs))

# Pattern- and target-specific assignments are propagated to the prerequisites
# and override the global assignments. Therefore, we introduce the following
# match-anything pattern assignments to prevent that (the eval/value combination
# is required to keep the original global values of the variables without
# changing their flavors, i.e. keep them recursively expanded):
${eval %: FCFLAGS= $(value FCFLAGS)}

# Target-specific variables for objects that use the inline library:
$(consumer_objs): FCFLAGS+= $(inlib_inline_FCFLAGS)

# All object that can be built with the inline library depend on it:
$(consumer_objs): $(inlib_name)

radiation_mcica_acc_lw.o: FCFLAGS += -gpu=maxregcount:96
radiation_mcica_acc_sw.o: FCFLAGS += -gpu=maxregcount:96
radiation_two_stream.o: FCFLAGS += -gpu=maxregcount:96
radiation_lw_derivatives.o: FCFLAGS += -gpu=maxregcount:96
endif

# The inline library generation rule. Note that the source files are provided
# to the compiler in the topological order. We also have to delete the partially
# generated library if the compiler fails:
$(inlib_name): $(inlib_OBJS) $(inlib_SOURCES)
	rm -rf $@ && $(FC) $(FCFLAGS) $(inlib_extract_FCFLAGS) $(inlib_SOURCES)
endif

clean:
	rm -f *.o $(LIBRAD)
	rm -rf $(inlib_name)

include Makefile_deps
