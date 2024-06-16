
#FC = mpifort
FC=mpif90

FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -finit-real=NaN -finit-integer=-2147483648 -g -fbacktrace -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero,overflow,denormal 
# FFLAGS = -fbounds-check -Wall -Wunused -O3 -fdefault-real-8 -fdefault-double-8

SRC = $(SPECRT_HOME)/source

default: compile

compile: $(SRC)/common.f90 $(SRC)/disk.f90 $(SRC)/auxiliary.f90 $(SRC)/continuous_opacity.f90 $(SRC)/gas_state.f90 $(SRC)/grey.f90 $(SRC)/flux_feautrier.f90 input.in
	cd source; \
	$(FC) $(FFLAGS) -c common.f90
	cd source; \
	$(FC) $(FFLAGS) -c disk.f90
	cd source; \
	$(FC) $(FFLAGS) -c auxiliary.f90
	cd source; \
	$(FC) $(FFLAGS) -c continuous_opacity.f90
	cd source; \
	$(FC) $(FFLAGS) -c gas_state.f90
	cd source; \
	$(FC) $(FFLAGS) -c grey.f90
	cd source; \
	$(FC) $(FFLAGS) common.o auxiliary.o disk.o continuous_opacity.o gas_state.o grey.o flux_feautrier.f90 -o flux_feautrier.x
clean:
	rm -f *.o *.mod
