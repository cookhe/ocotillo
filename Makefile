
FC = mpifort

FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -finit-real=NaN -finit-integer=-2147483648 -g -fbacktrace -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero,overflow,denormal 
# FFLAGS = -fbounds-check -Wall -Wunused -O3 -fdefault-real-8 -fdefault-double-8

default: compile

compile: common.f90 disk.f90 auxiliary.f90 continuous_opacity.f90 gas_state.f90 flux_feautrier.f90 input.in
	$(FC) $(FFLAGS) -c common.f90
	$(FC) $(FFLAGS) -c disk.f90
	$(FC) $(FFLAGS) -c auxiliary.f90
	$(FC) $(FFLAGS) -c continuous_opacity.f90
	$(FC) $(FFLAGS) -c gas_state.f90
	$(FC) $(FFLAGS) common.o auxiliary.o disk.o continuous_opacity.o gas_state.o flux_feautrier.f90 -o flux_feautrier.x
clean:
	rm -f *.o *.mod
