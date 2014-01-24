
.PHONY: all test clean

FORT=gfortran

F90OPTS= -Wall -Wextra -fbounds-check 
#F90OPTS= -Wall -fbounds-check -funroll-loops -finline -ffast-math -O2

LDLIBS= -lsatlas -lcompla
INCLUDE= -I/usr/include/libcompla # link to location of .mod files

all: main

main: main.f90
	$(FORT) $(F90OPTS) $(INCLUDE) -c main.f90
	$(FORT) -o main $(F90OPTS) $(LDLIBS) main.o

test: main
	./main

clean:
	rm -f main
	rm -f *.o *.mod
