.PHONY: clean

FC = gfortran
FFLAGS = -O0 -Wall -g -fcheck=array-temps,bounds,do,mem,pointer,recursion -ffpe-trap=invalid,zero,overflow,underflow

graph: output.data 
	gnuplot splot_output_data.gp

output.data: runme
	./runme

runme: main.o const.o vector.o spiral.o
	$(FC) $(FFLAGS) -o runme main.o const.o vector.o spiral.o

main.o: main.f90 const.o vector.o spiral.o
	$(FC) $(FFLAGS) -c main.f90

const.o: const.f90
	$(FC) $(FFLAGS) -c const.f90

vector.o: vector.f90 const.o
	$(FC) $(FFLAGS) -c vector.f90

spiral.o: spiral.f90 const.o vector.o
	$(FC) $(FFLAGS) -c spiral.f90

clean:
	rm -f *.o runme *.mod output.data
