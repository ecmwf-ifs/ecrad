LIBDRHOOK = ../lib/libdrhook.a

all: dummy

# Dummy version of Dr Hook creates the YOMHOOK module but the DR_HOOK
# routine does nothing
dummy: $(LIBDRHOOK)

$(LIBDRHOOK): yomhook_dummy.o
	ar -r $(LIBDRHOOK) yomhook_dummy.o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	rm -f *.o $(LIBDRHOOK)

.PHONY: dummy
