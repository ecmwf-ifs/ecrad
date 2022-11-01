#TEST = test_solver
#TEST = test_cloud_generator
TEST = test_spartacus_math
TEST = test_aerosol_optics_description

# In dependency order:
SOURCES = ecrad_driver_read_input.F90 \
	ecrad_driver_config.F90

OBJECTS := $(SOURCES:.F90=.o)
EXECUTABLE = ../bin/ecrad
IFS_EXECUTABLE = ../bin/ecrad_ifs

all: $(EXECUTABLE) $(TEST)

$(EXECUTABLE): $(OBJECTS) ../lib/*.a ecrad_driver.o
	$(FC) $(FCFLAGS) ecrad_driver.o $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

ecrad_ifs_driver: $(OBJECTS) ../lib/*.a ecrad_ifs_driver.o
	$(FC) $(FCFLAGS) ecrad_ifs_driver.o $(OBJECTS) -lifs $(LIBS) -o $(IFS_EXECUTABLE)

$(TEST): $(TEST).F90 ../lib/*.a
	$(FC) $(FCFLAGS) $(TEST).F90 $(LIBS) -o $(TEST)

# Note that the dependence on mod files can mean that rerunning "make"
# recreates the executable
%.o: %.F90 ../mod/*.mod
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(EXECUTABLE) $(IFS_EXECUTABLE) $(TEST)

ecrad_driver.o: ecrad_driver_config.o ecrad_driver_read_input.o
ecrad_driver_read_input.o: ecrad_driver_config.o
