SOURCES = ice_effective_radius.F90 liquid_effective_radius.F90 radiation_scheme.F90
MAKE_INCLUDES = ../bin/make_intfbl.1.pl
INCLUDE_DIR = ../include

OBJECTS := $(SOURCES:.F90=.o)
LIBIFS = ../lib/libifs.a

all: $(LIBIFS)

deps: includes

$(LIBIFS): $(OBJECTS)
	ar r $(LIBIFS) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

includes:
	LOC_INTFBDIR=$(INCLUDE_DIR) INTFBDIR=$(INCLUDE_DIR) $(MAKE_INCLUDES) $(SOURCES)

clean:
	rm -f *.o $(LIBIFS)

#yomlun_ifsaux.o, yomhook.o, random_numbers_mix.o: parkind1.o
#yomlun.o : yomlun_ifsaux.o

.PHONY: deps includes
