default: compile
compile:
	mpif90 -fbounds-check -Wall -Wunused -O3 -fdefault-real-8 -fdefault-double-8 -c auxiliary.f90
	mpif90 -fbounds-check -Wall -Wunused -O3 -fdefault-real-8 -fdefault-double-8 auxiliary.o flux_feautrier.f90 -o flux_feautrier.x
clean:
	rm -f *.o *.mod
