#TEST = test_solver
#TEST = test_cloud_generator
TEST = test_spartacus_math

# In dependency order:
SOURCES = ecrad_driver_read_input.F90 \
	ecrad_driver_config.F90 \
	ecrad_driver.F90 \

OBJECTS := $(SOURCES:.F90=.o)
EXECUTABLE = ../bin/ecrad

all: $(EXECUTABLE) $(TEST)

$(EXECUTABLE): $(OBJECTS) ../lib/*.a
	$(FC) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(TEST): $(TEST).F90 ../lib/*.a
	$(FC) $(FCFLAGS) $(TEST).F90 $(LIBS) -o $(TEST)

# Note that the dependence on mod files can mean that rerunning "make"
# recreates the executable
%.o: %.F90 ../mod/*.mod
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(OBJECTS) $(EXECUTABLE) $(TEST)

ecrad_driver.o: ecrad_driver_config.o ecrad_driver_read_input.o
ecrad_driver_read_input.o: ecrad_driver_config.o
