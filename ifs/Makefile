SOURCES = ice_effective_radius.F90 liquid_effective_radius.F90 \
	radiation_scheme.F90 radiation_setup.F90 yoerdu.F90 \
	yomrip.F90 yoephy.F90 yoecld.F90 yoe_spectral_planck.F90 \
	cloud_overlap_decorr_len.F90 yoerad.F90 yoethf.F90 satur.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBIFS = ../lib/libifs.a

all: $(LIBIFS)

deps: includes

$(LIBIFS): $(OBJECTS)
	ar r $(LIBIFS) $(OBJECTS)

%.o: %.F90 ../lib/libifsaux.a ../lib/libradiation.a
	$(FC) $(FCFLAGS) -c $<

includes:
	fcm make --config-file=../bin/fcm-make-interfaces.cfg interfaces.ns-incl=ifs interfaces.source=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))/..

clean:
	rm -f *.o $(LIBIFS)

.PHONY: deps includes clean-deps

liquid_effective_radius.o: yoerdu.o
radiation_scheme.o: radiation_setup.o
radiation_setup.o: yoephy.o
cloud_overlap_decorr_len.o: yoecld.o
cos_sza.o ice_effective_radius.o liquid_effective_radius.o radiation_scheme.o radiation_setup.o: yoerad.o
yoerad.o: yoe_spectral_planck.o
satur.o: yoethf.o
