
.PHONY: all clean force

FORT=gfortran
AR=ar

FORTOPTS= -Wall -fbounds-check -funroll-loops -finline -ffast-math -O2

LDLIBS= -lsatlas -lcompla
INC= -I/usr/include/libcompla

all: p03

p03: p03.f08 force
	$(FORT) $(FORTOPTS) $(INC) -c p03.f08
	$(FORT) -o p03 $(FORTOPTS) $(LDLIBS) p03.o
	./p03

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f *.o *.mod *.out *.a
