
.PHONY: all clean force

FORT=gfortran
AR=ar

#FORTOPTS= -Wall -Wextra -fbounds-check 
FORTOPTS= -Wall -fbounds-check -funroll-loops -finline -ffast-math -O2

LDLIBS= -lsatlas

#all: compla.a
all: libcompla.so 

#compla.o: compla.f90
#	$(FORT) $(FORTOPTS) $(LDLIBS) -c compla.f90
#
#compla.a: compla.o
#	$(AR) rc compla.a compla.o

libcompla.so: compla.f08
	$(FORT) $(FORTOPTS) $(LDLIBS) -shared compla.f08 -o libcompla.so -fPIC

test: libcompla.so test.f08 force
	$(FORT) $(FORTOPTS) -c test.f08
	$(FORT) -o test $(FORTOPTS) $(LDLIBS) -lcompla test.o
	./test

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f libcompla.so test
	rm -f *.o *.mod *.out *.a
