
.PHONY: all clean

FORT=gfortran
AR=ar

#F90OPTS= -Wall -Wextra -fbounds-check 
F90OPTS= -Wall -fbounds-check -funroll-loops -finline -ffast-math -O2

LDLIBS= -lsatlas

#all: compla.a
all: libcompla.so 

#compla.o: compla.f90
#	$(FORT) $(F90OPTS) $(LDLIBS) -c compla.f90
#
#compla.a: compla.o
#	$(AR) rc compla.a compla.o

libcompla.so: compla.f90
	$(FORT) $(F90OPTS) $(LDLIBS) -shared compla.f90 -o libcompla.so -fPIC

test: libcompla.so test.f90  
	$(FORT) $(F90OPTS) -c test.f90
	$(FORT) -o test $(F90OPTS) $(LDLIBS) -lcompla test.o
	./test

clean:
	rm -f libcompla.so test
	rm -f *.o *.mod *.out *.a
