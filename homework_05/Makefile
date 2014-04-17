
.PHONY: all clean force

FORT=gfortran
AR=ar

FORTOPTS= -Wall -fbounds-check -funroll-loops -finline -ffast-math -O2

LDLIBS= -lsatlas -lcompla
INC= -I/usr/include/libcompla

all: p02

p02: p02.f08 force
	$(FORT) $(FORTOPTS) $(INC) -c p02.f08
	$(FORT) -o p02 $(FORTOPTS) $(LDLIBS) p02.o
	./p02

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f *.o *.mod *.out *.a p02
