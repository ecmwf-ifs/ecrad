SOURCES = ice_effective_radius.F90 liquid_effective_radius.F90 \
	radiation_scheme.F90 radiation_setup.F90 yoecldp.F90 yoerdu.F90 \
	yoethf.F90 satur.F90 type_model.F90 yomrip0.F90 yomrip.F90 \
	model_general_conf_mod.F90 model_physics_ecmwf_mod.F90 \
	yoephy.F90 yoecld.F90 model_physics_radiation_mod.F90 \
	yoeaerop.F90 cloud_overlap_decorr_len.F90 yoe_mcica.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBIFS = ../lib/libifs.a

all: $(LIBIFS)

deps: includes

$(LIBIFS): $(OBJECTS)
	ar r $(LIBIFS) $(OBJECTS)

%.o: %.F90 ../mod/*.mod
	$(FC) $(FCFLAGS) -c $<

includes:
	cd .. && fcm make --config-file=bin/fcm-make-interfaces.cfg interfaces.ns-incl=ifs interfaces.source=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/..

clean:
	rm -f *.o $(LIBIFS)

.PHONY: deps includes clean-deps

liquid_effective_radius.o: yoecldp.o yoerdu.o
radiation_scheme.o: type_model.o yomrip0.o radiation_setup.o
type_model.o: model_general_conf_mod.o model_physics_ecmwf_mod.o model_physics_radiation_mod.o
model_general_conf_mod.o: yomrip.o
model_physics_ecmwf_mod.o: yoephy.o yoecld.o yoecldp.o
model_physics_radiation_mod.o: radiation_setup.o yoe_mcica.o
radiation_setup.o: yoeaerop.o
satur.o: yoethf.o
cloud_overlap_decorr_len.o: yoecld.o
