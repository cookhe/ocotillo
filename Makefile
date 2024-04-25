default: compile
compile:
	gfortran -c -O3 auxiliary.f90                      
	gfortran -O3 auxiliary.o flux_feautrier.f90 -o flux_feautrier.x	
clean:
	rm -f *.o *.mod
