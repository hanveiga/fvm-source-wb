all: 
	nvcc -c Cuda_functions.cu
	gfortran -c dg.f90
	gfortran -c legendre.f90
	gfortran dg.o legendre.o Cuda_functions.o -o dg -lstdc++ -L/usr/local/cuda-7.0/lib64/ -lcuda -lcudart 


