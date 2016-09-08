#include <stdio.h>
#include "cuda.h"
#define BLOCK 256
__device__  double legendre(double x, int n, int sq){
  double legendre;
  x=min(max(x,-1.0),1.0);
  switch (n) {
  case 0:
    legendre=1.0;
    break;
  case 1:
    legendre=x;
    break;
  case 2:
    legendre=0.5*(3.0*x*x-1.0);
    break;
  case 3:
    legendre=(2.5*x*x*x-1.5*x);
    break;
  case 4:
    legendre=0.125*(35.0*x*x*x*x-30.0*x*x+3.0);
    break;
  case 5:
    legendre=0.125*(63.0*pow(x,5)-70.0*pow(x,3)+15.0*x);
    break;
  case 6:
    legendre=1.0/16.0*(231.0*pow(x,6)-315.0*pow(x,4)+105.0*pow(x,2)-5.0);
    break;
  }
  if(sq==1)
    legendre=sqrt(2.0*double(n)+1.0)*legendre;
  return legendre;
}

__device__  double legendre_prime(double x, int n){
  double legendre_prime;
  x=min(max(x,-1.0),1.0);
  switch (n) {
  case 0:
    legendre_prime=0.0;
    break;
  case 1:
    legendre_prime=1.0;
    break;
  case 2:
    legendre_prime=3*x;
    break;
  case 3:
    legendre_prime=0.5*(15.0*x*x-3.0);
    break;
  case 4:
    legendre_prime=0.125*(140.0*x*x*x-60.0*x);
    break;
  case 5:
    legendre_prime=0.125*(315.0*pow(x,4)-210.0*pow(x,2)+15.0);
    break;
  case 6:
    legendre_prime=1.0/16.0*(1386.0*pow(x,5)-1260.0*pow(x,3)+210.0*x);
    break;
  }
  //legendre_prime=sqrt((double)(2*n+1))*legendre_prime;
  return legendre_prime;
}

__device__  double minmod(double x, double y, double z){
  int s;
  s=copysign(1.0,x);
  if(copysign(1.0,y) == s && copysign(1.0,z) == s)
    return (double)s*min(fabs(x),min(fabs(y),fabs(z)));
  else
     return 0.0;
}

__global__ void initial_dg (double* rho, double* rho_right, double* rho_left, double* flux_vol, double* flux_right, double* flux_left, double* chsi_quad, double* w_quad, int nx, int n, int nquad)
{
  int icell, quad, i;
  double chsi_left=-1,chsi_right=1, velocity=1.0;
  double flux_quad[20];
  double flux;
  icell = blockDim.x * blockIdx.x + threadIdx.x;
  if( icell < nx ){
    // Compute flux at quadrature points
    for(quad = 0; quad < nquad; quad++){
      flux=0.0;
      for(i=0; i<n; i++)
	flux += rho[icell+i*nx]*legendre(chsi_quad[quad],i,1);
      flux_quad[quad] = velocity * flux * w_quad[quad];
    }
    
    // Compute volume integral DG term, compute left and right fluxes
    rho_left[icell]=0.0;
    rho_right[icell]=0.0;
    for(i=0; i<n; i++){
      flux=0.0;
      for(quad = 0; quad < nquad; quad++)
	flux += legendre_prime(chsi_quad[quad],i)*flux_quad[quad];
      flux_vol[icell+i*nx] = flux*sqrt(2.0*double(i)+1.0);
      rho_left[icell] += rho[icell+i*nx]*legendre(chsi_left,i,1);
      rho_right[icell]+= rho[icell+i*nx]*legendre(chsi_right,i,1);
    }
    //Compute flux at left and right points
    flux_left[icell]=velocity*rho_left[icell];
    flux_right[icell]=velocity*rho_right[icell];;
  }
}

__global__ void Lax_Friedrich (double* flux_face, double* rho_right, double* rho_left, double* flux_right, double* flux_left, int nx)
{
  int iface, ileft, iright;
  double velocity=1.0,cmax;
  iface = blockDim.x * blockIdx.x + threadIdx.x;
  if( iface < nx+1 ){
    ileft=iface-1;
    iright=iface;
    //Periodic boundary conditions
    if(iface==0)ileft=nx-1;
    if(iface==nx)iright=0;
    cmax=fabs(velocity);
    flux_face[iface]=0.5*(flux_right[ileft]+flux_left[iright]-cmax*(rho_left[iright]-rho_right[ileft]));
  }
}

__global__ void final_dg (double* flux_face, double* drhodt, double* flux_vol, int nx, int n, double oneoverdx)
{  
  int icell, i;
  double chsi_left=-1,chsi_right=1;
  icell = blockDim.x * blockIdx.x + threadIdx.x;
  if( icell < nx ){
     for(i=0; i<n; i++)
       drhodt[icell+i*nx]=oneoverdx*(flux_vol[icell+i*nx]-(flux_face[icell+1]*legendre(chsi_right,i,0)-flux_face[icell]*legendre(chsi_left,i,0))*sqrt(2.0*double(i)+1.0));
     
  }
}

__global__ void limiter (double* rho, double* rho_lim, int nx, int n, int use_extremum)
{  
  int icell,i,ileft,iright;
  double rhoL,rhoM,rhoR;
  double rho_min;
  double coeff_i,coeff_ip1,coeff,D2rho;
 
  icell = blockDim.x * blockIdx.x + threadIdx.x;
 
  if( icell < nx){  
    ileft=icell-1;
    iright=icell+1;
    if(icell==0)ileft=nx-1;
    if(icell==nx-1)iright=0;
    for(i=n-1; i>0 ; i--){
      coeff_i=sqrt(2.0*double(i-1)+1.0);
      coeff_ip1=sqrt(2.0*double(i)+1.0)*(2.0*double(i)-1);
      rhoL=(rho[icell+(i-1)*nx]-rho[ileft+(i-1)*nx])*coeff_i;
      rhoR=(rho[iright+(i-1)*nx]-rho[icell+(i-1)*nx])*coeff_i;
      rhoM=rho[icell+i*nx]*coeff_ip1;
      rho_min=minmod(rhoL,rhoM,rhoR)/coeff_ip1;
      if( use_extremum && rho_min != rho[icell+i*nx] && min(rhoL*rhoR,rho[ileft+i*nx]*rho[iright+i+nx]) < 0.0){
	rhoL=2.0*(rho[icell+i*nx]-rho[ileft+i*nx]);
	rhoR=2.0*(rho[iright+i*nx]-rho[icell+i*nx]);
	rhoM=0.5*(rho[iright+i*nx]-rho[ileft+i*nx]);
	D2rho=minmod(rhoL,rhoM,rhoR);
	coeff=1.0;
	if(rhoM != 0.0)coeff=D2rho/rhoM;
	rho_lim[icell+i*nx]=rho_min+(rho[icell+i*nx]-rho_min)*coeff;
      }
      else
	rho_lim[icell+i*nx]=rho_min;
      if(rho_lim[icell+i*nx] == rho[icell+i*nx]) break;  
    }
  
  }
}

__global__ void copy(double* rho, double* rho_lim, int nx, int size)
{  
  int id;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if(id < size)
    rho[id + nx] = rho_lim[id];
}
__global__ void sum3 (double* out, double* A, double* B, double* C, double alpha, double beta, double gamma, int size)
{  
  int id;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size )
    out[id] = A[id]*alpha + B[id]*beta + C[id]*gamma;
}

__global__ void plus_equal (double* out, double* A, double* B, double alpha, double beta, int size)
{  
  int id;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size )
    out[id] += A[id]*alpha + B[id]*beta;
}
__global__ void sum2 (double* out, double* A, double* B, double beta, int size)
{  
  int id;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size )
    out[id] = A[id] + B[id]*beta;
}

extern "C" void device_compute_update_ (double** rho_d,  double**  drhodt_d, double** flux_vol_d, double** rho_right_d, double** rho_left_d, double** flux_right_d, double** flux_left_d, double** flux_face_d, int* Nx, int* N, double** chsi_quad_d, double** w_quad_d, int* Nquad){
     
  int nx = *Nx;
  int n = *N;
  int nquad = *Nquad;
  double dx;
  double oneoverdx;
  cudaError_t error = cudaGetLastError();
  dx=1.0/(double)nx;
  oneoverdx=1.0/dx;
  initial_dg<<<(nx+BLOCK-1)/BLOCK,BLOCK>>>(*rho_d, *rho_right_d, *rho_left_d, *flux_vol_d, *flux_right_d, *flux_left_d, *chsi_quad_d, *w_quad_d, nx, n, nquad);
  Lax_Friedrich<<<(nx+BLOCK)/BLOCK,BLOCK>>>(*flux_face_d, *rho_right_d, *rho_left_d, *flux_right_d, *flux_left_d, nx);
  final_dg<<<(nx+BLOCK-1)/BLOCK,BLOCK>>>(*flux_face_d, *drhodt_d, *flux_vol_d, nx, n, oneoverdx);
  if(error != cudaSuccess){
    printf("CUDA error compute update: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}
extern "C" void device_limiter_(double** rho, double** rho_lim , int* Nx, int* N, int* Use_extremum){
  int n = *N;
  int nx = *Nx;
  int use_extremum = *Use_extremum;
  if(n>1){
    cudaMemcpy( *rho_lim, *rho,  n*nx * sizeof(double) ,cudaMemcpyDeviceToDevice);
    limiter<<<(nx+BLOCK-1)/BLOCK,BLOCK>>>(*rho, *rho_lim, nx, n, use_extremum);
    cudaMemcpy( *rho, *rho_lim,  n*nx * sizeof(double) ,cudaMemcpyDeviceToDevice);
  }
}

extern "C" void device_sum3_ (double** out,  double**  A, double** B, double** C, double* alpha, double* beta, double* gamma, int* n){
  double a = *alpha;
  double b = *beta;
  double c = *gamma;
  int size = *n;
  cudaError_t error = cudaGetLastError();
  sum3<<<(size+BLOCK-1)/BLOCK,BLOCK>>>(*out, *A, *B, *C, a, b, c, size);
  if(error != cudaSuccess){
    printf("CUDA error sum3: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}
extern "C" void device_plus_equal_ (double** out,  double**  A, double** B, double* alpha, double* beta, int* n){
  double a = *alpha;
  double b = *beta;
  int size = *n;
  cudaError_t error = cudaGetLastError();
  plus_equal<<<(size+BLOCK-1)/BLOCK,BLOCK>>>(*out, *A, *B, a, b, size);
  if(error != cudaSuccess){
    printf("CUDA error sum3: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}
extern "C" void device_sum2_ (double** out,  double**  A, double** B, double* beta, int* n){
  double b = *beta;
  int size = *n;
  cudaError_t error = cudaGetLastError();
  sum2<<<(size+BLOCK-1)/BLOCK,BLOCK>>>(*out, *A, *B, b, size);
  if(error != cudaSuccess){
    printf("CUDA error sum2: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void malloc_gpu_ (double **array, int* Size) {
  int size = *Size;
  printf("size %d\n",size);
  cudaError_t error = cudaGetLastError();
  cudaMalloc ( array, size * sizeof(double));
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error malloc: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void h2d_ (double *array, double **darray, int* Size) {
  int size = *Size;
  cudaError_t error = cudaGetLastError();
  cudaMemcpy( *darray, array,  size * sizeof(double) ,cudaMemcpyHostToDevice);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error h2d: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void d2h_ (double **darray, double *array, int* Size) {
  int size = *Size;
  cudaError_t error = cudaGetLastError();
  cudaMemcpy( array, *darray,  size * sizeof(double) ,cudaMemcpyDeviceToHost);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error d2h: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void setdevice_ (int *Device) {
  int device = *Device;
  cudaError_t error = cudaGetLastError();
  cudaSetDevice(device);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error setting device: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  printf("Device Number: %d\n", device);
  printf("Device name: %s\n", prop.name);
}

extern "C" void devices_ () {
  int nDevices;
  cudaGetDeviceCount(&nDevices);
  printf("Devices: %d\n",nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
  }
}
