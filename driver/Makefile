#TEST = test_solver
#TEST = test_cloud_generator
TEST = test_spartacus_math

# In dependency order:
SOURCES = ecrad_driver_read_input.F90 \
	ecrad_driver_config.F90 \
	ecrad_standalone.F90

OBJECTS := $(SOURCES:.F90=.o)
EXECUTABLE = ../bin/ecrad
IFSEXE = ../bin/ecrad_ifs

all: $(EXECUTABLE) $(TEST) $(IFSEXE)

# Note that the dependence on mod files can mean that rerunning "make"
# recreates the executable
%.o: %.F90 ../mod/*.mod
	$(FC) $(FCFLAGS) -c $<

$(EXECUTABLE): ecrad_driver.o $(OBJECTS) ../lib/*.a
	$(FC) $(FCFLAGS) ecrad_driver.o $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(TEST): $(TEST).F90 ../lib/*.a
	$(FC) $(FCFLAGS) $(TEST).F90 $(LIBS) -o $(TEST)

$(IFSEXE): ecrad_ifs_driver.o $(OBJECTS) ../lib/*.a
	$(FC) $(FCFLAGS) ecrad_ifs_driver.o ecrad_ifs.o $(OBJECTS) -lifs $(LIBS) -o $(IFSEXE)


clean:
	rm -f *.o ecrad_driver.o ecrad_ifs_driver.o $(OBJECTS) $(EXECUTABLE) $(TEST) $(IFSEXE)

ecrad_driver.o: ecrad_driver_config.o ecrad_driver_read_input.o ecrad_standalone.o
ecrad_ifs_driver.o: ecrad_ifs.o ecrad_driver_config.o ecrad_driver_read_input.o
ecrad_driver_read_input.o: ecrad_driver_config.o
ecrad_standalone.o: ecrad_driver_config.o ecrad_driver_read_input.o
