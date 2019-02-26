SOURCES = radsurf_properties.F90 radsurf_flux.F90 radsurf_intermediate.F90 radsurf_save.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBSURF = ../lib/libradsurf.a

all: $(LIBSURF)

$(LIBSURF): $(OBJECTS)
	ar r $(LIBSURF) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBSURF)

radsurf_flux.o: radsurf_properties.o
radsurf_intermediate.o: radsurf_flux.o radsurf_properties.o 
radsurf_save.o: radsurf_flux.o
