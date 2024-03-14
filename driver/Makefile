TEST_PROGRAMS = test_cloud_generator test_spartacus_math \
	test_aerosol_optics_description test_random_number_generator \
	test_fast_expm

# In dependency order:
SOURCES = ecrad_driver_read_input.F90 \
	ecrad_driver_config.F90

OBJECTS := $(SOURCES:.F90=.o)
EXECUTABLE = ../bin/ecrad
IFS_EXECUTABLE = ../bin/ecrad_ifs
IFS_BLOCKED_EXECUTABLE = ../bin/ecrad_ifs_blocked

all: driver ifs_driver test_programs

driver: $(EXECUTABLE)

ifs_driver: $(IFS_EXECUTABLE) $(IFS_BLOCKED_EXECUTABLE)

test_programs: $(TEST_PROGRAMS)

# Link ecrad executable; add "-lifs" if you want to use the "satur"
# routine in ecrad_driver.F90
$(EXECUTABLE): $(OBJECTS) ../lib/*.a ecrad_driver.o
	$(FC) $(FCFLAGS) ecrad_driver.o $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(IFS_EXECUTABLE): $(OBJECTS) ../lib/*.a ecrad_ifs_driver.o
	$(FC) $(FCFLAGS) ecrad_ifs_driver.o $(OBJECTS) -lifs $(LIBS) -o $(IFS_EXECUTABLE)

$(IFS_BLOCKED_EXECUTABLE): $(OBJECTS) ../lib/*.a ecrad_ifs_driver_blocked.o
	$(FC) $(FCFLAGS) ifs_blocking.o ecrad_ifs_driver_blocked.o $(OBJECTS) -lifs $(LIBS) -o $(IFS_BLOCKED_EXECUTABLE)

test_%: test_%.F90 ../lib/*.a
	$(FC) $(FCFLAGS) $< $(LIBS) -o $@

#$(TEST): $(TEST).F90 ../lib/*.a
#	$(FC) $(FCFLAGS) $(TEST).F90 $(LIBS) -o $(TEST)

# Note that the dependence on mod files can mean that rerunning "make"
# recreates the executable
%.o: %.F90 ../lib/*.a
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o $(EXECUTABLE) $(IFS_EXECUTABLE) $(IFS_BLOCKED_EXECUTABLE) \
		$(TEST_PROGRAMS)

ecrad_driver.o: ecrad_driver_config.o ecrad_driver_read_input.o
ecrad_ifs_driver.o: ecrad_driver_config.o ecrad_driver_read_input.o
ecrad_ifs_driver_blocked.o: ecrad_driver_config.o ecrad_driver_read_input.o ifs_blocking.o
ecrad_driver_read_input.o ifs_blocking.o: ecrad_driver_config.o

.PHONY: driver ifs_driver test_programs all
