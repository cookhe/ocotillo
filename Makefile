
FC = mpif90

#FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -finit-real=NaN -finit-integer=-2147483648 -g -fbacktrace -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero,overflow,denormal 
FFLAGS = -fbounds-check -Wall -Wunused -O3 -fdefault-real-8 -fdefault-double-8

default: compile

compile: grid.f90 disk.f90 auxiliary.f90 flux_feautrier.f90 input.in
	$(FC) $(FFLAGS) -c grid.f90
	$(FC) $(FFLAGS) -c disk.f90
	$(FC) $(FFLAGS) -c auxiliary.f90
	$(FC) $(FFLAGS) grid.o auxiliary.o disk.o flux_feautrier.f90 -o flux_feautrier.x
clean:
	rm -f *.o *.mod
