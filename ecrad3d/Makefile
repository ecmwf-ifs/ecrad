SOURCES = ecrad3d_config.F90 ecrad3d_geometry.F90 ecrad3d_layer_solutions.F90 ecrad3d_solver_sw.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBECRAD3D = ../lib/libecrad3d.a

all: $(LIBECRAD3D)

$(LIBECRAD3D): $(OBJECTS)
	ar r $(LIBECRAD3D) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(LIBECRAD3D)

ecrad3d_solver_sw.o: ecrad3d_geometry.o ecrad3d_layer_solutions.o ecrad3d_config.o
