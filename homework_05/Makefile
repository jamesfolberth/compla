
.PHONY: all clean force

FORT=gfortran
AR=ar

FORTOPTS= -Wall -Wextra -fbounds-check -funroll-loops -finline -ffast-math -O2

LDLIBS= -lsatlas -lcompla -L/opt/hdf5/lib -lhdf5_fortran
INC= -I/usr/include/libcompla -I/opt/hdf5/include
LDADDPATH= -Wl,-rpath=/opt/hdf5/lib # give ld the path to libhdf5_fortran.*

all: p02

p02: p02.f08 force
	$(FORT) $(FORTOPTS) $(INC) -c p02.f08
	$(FORT) -o p02 $(FORTOPTS) $(LDLIBS) $(LDADDPATH) p02.o
	./p02

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f *.o *.mod *.out *.a p02
