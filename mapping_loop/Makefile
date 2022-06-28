
$(SOURCE) ?= .

.DEFAULT_GOAL := library.a

LIBRARY_OBJECTS = library.o

library.o: ${SOURCE}/library.f90
	gfortran -c -Ofast -fPIC -o library.o ${SOURCE}/library.f90

library.a: $(LIBRARY_OBJECTS)
	rm -f library.a
	ar cvr $@ $(LIBRARY_OBJECTS)

.PHONY:	clean

clean:
	-rm -f *.o *.a *.mod *.smod *~ \#*\# .*~ .\#*
