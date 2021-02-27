SOURCES = tcrad_2region_solver.F90 tcrad_3region_solver.F90 \
	tcrad_layer_solutions.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBMR = ../lib/libtcrad.a

all: $(LIBMR)

$(LIBMR): $(OBJECTS)
	ar r $(LIBMR) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBMR)

tcrad_2region_solver.o tcrad_3region_solver.o: tcrad_layer_solutions.o

tcrad_2region_solver.o tcrad_3region_solver.o: tcrad_matrix.h \
	tcrad_two_stream_flux.h \
	tcrad_radiance.h tcrad_overlap.h tcrad_region.h \
	tcrad_flux_interface.h
