CF=gfortran
CCUDA=nvcc
FFLAGS=-lgfortran
CUDAFLAGS=-lstdc++ -L/usr/local/cuda-7.0/lib64/ -lcuda -lcudart 

OBJS=parameters_dg_2d.o legendre.o limiters.o Cuda_functions.o

exe: $(OBJS) benchmark_2d_dg.f90
	$(CF) -o dg2d benchmark_2d_dg.f90 $(OBJS) $(FFLAGS) $(CUDAFLAGS)

parameters_dg_2d.o: parameters_dg_2d.f90
	$(CF) $(FFLAGS) -c $<

legendre.o: legendre.f90
	$(CF) $(FFLAGS) -c $<

limiters.o: limiters.f90
	$(CF) $(FFLAGS) -c $<
	
Cuda_functions.o: Cuda_functions.cu
	$(CCUDA) -c $<
clean:	
	rm -f *.f90~ 
	rm -f *.mod
	rm dg2d
