#include <stdio.h>
#include "cuda.h"
#include "global.h"

#define BLOCK 512

__constant__ double sqrt_mod[5];
__constant__ double xquad[5];
__constant__ double yquad[5];
__constant__ double wxquad[5];
__constant__ double wyquad[5];

__device__  double legendre(double x, int n, int sq){
  double legendre;
  x=min(max(x,-1.0),1.0);
  switch (n) {
  case 0:
    legendre=1.;
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
    legendre *= sqrt_mod[n];
  return legendre;
}

__device__  double legendre_prime(double x, int n, int sq){
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
  if(sq==1)
    legendre_prime *= sqrt_mod[n];
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

__device__  double compute_speed(double rho, double vx, double vy, double p, double gamma){
  double cs, speed;
  cs=sqrt(gamma*max(p,1E-10)/max(rho,1E-10));
  speed=sqrt(vx*vx+vy*vy)+cs;
  return speed;
}

__device__ int BC(int index, int size, int bc){
  if (bc == 1){//periodic
    if (index == -1) 
      index = size-1;
    else if (index == size)
      index = 0;
  }
  else if (bc == 2 || bc == 3){//transmissive or reflective
    if (index == -1) 
      index++;
    else if (index == size)
      index--;
  }
  return index;
}

__global__  void gaussian_quad(double* quads, double* weights, int n){
  int i, id;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  double dpi=acos(-1.0);
  double xx,lp;
  if (id < n){
    xx=(1.0-0.125/n/n+0.125/n/n/n)*
      cos(dpi*(4.0*double(id)-1.0)/(4.0*double(n)+2.0));
    for(i=0; i<100; i++)
      xx=xx-legendre(xx,n,1)/legendre_prime(xx,n,1);
    xx = quads[id] = -xx;
    lp = legendre_prime(xx,n,1);
    weights[id] = 2*(2.0*double(n)+1.0)/(1.0-xx*xx)/(lp*lp);
  }
}


__global__ void get_modes_from_nodes(double* nodes, double* modes, int nx, int ny, int mx, int my, int nvar){

  int id, icell, jcell, imod, jmod, var;
  int xq, yq, idq;
  int a = my;
  int b = mx*a;
  int c = ny*b;
  int d = nx*c;
  int size = nvar*d;
  double val=0.0;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  var = id/d;
  jmod = id - var*d;
  icell = jmod/c;
  jmod -= icell*c;
  jcell = jmod/b;
  jmod -= jcell*b;
  imod = jmod/a;
  jmod -= imod*a;
  idq = jcell*b + icell*c + var*d;
  
  if( id < size ){
    for( xq=0; xq < mx; xq++){
       for( yq=0; yq < my; yq++)
	 val += 0.25*nodes[yq+xq*a+idq]*legendre(xquad[xq],imod,1)*legendre(yquad[yq],jmod,1)
	   * wxquad[xq]*wyquad[yq];
    }
    modes[id] = val;
  }
}

__global__ void get_nodes_from_modes(double* modes, double* nodes, int nx, int ny, int mx, int my, int nvar){

  int id, icell, jcell, imod, jmod, var;
  int xq, yq, idq;
  int a = my;
  int b = mx*a;
  int c = ny*b;
  int d = nx*c;
  int size = nvar*d;
  double val=0.0;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  var = id/d;
  jmod = id - var*d;
  icell = jmod/c;
  jmod -= icell*c;
  jcell = jmod/b;
  jmod -= jcell*b;
  imod = jmod/a;
  jmod -= imod*a;
  idq = jcell*b + icell*c + var*d;
 
  if( id < size ){
    for( xq=0; xq < mx; xq++){
       for( yq=0; yq < my; yq++)
	 val += modes[yq+xq*a+idq]
	   *legendre(xquad[imod],xq,1)*legendre(yquad[jmod],yq,1);
    }
    nodes[id] = val;
  }
}

__global__ void compute_primitive(double* u, double* w, double gamma, int size){

  int id; 
  double rho,vx,vy;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size ){
    rho = w[id] = u[id];
    vx  = w[id+size] = u[id+size]/rho;
    vy  = w[id+size*2] = u[id+size*2]/rho;
    w[id+size*3] = (gamma-1.0)*( u[id+size*3] - 0.5*rho*(vx*vx+vy*vy));
  }
}

__global__ void compute_conservative(double* w, double* u, double gamma, int size){

  int id; 
  double rho,vx,vy;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size ){
    rho = u[id] = w[id];
    vx = w[id+size];
    vy = w[id+size*2];
    u[id+size]   = rho*vx;
    u[id+size*2] = rho*vy;
    u[id+size*3] = w[id+size*3]/(gamma-1.0) + 0.5*rho*(vx*vx+vy*vy);
  }
}

__global__ void compute_flux(double* u, double* w, double* flux1, double* flux2, int size){

  int id; 
  double rhovx,rhovy,e,vx,vy,p;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size ){
    vx  = w[id+size];  
    vy  = w[id+size*2];
    p   = w[id+size*3];  
    rhovx = u[id+size];  
    rhovy = u[id+size*2];
    e = u[id+size*3];
    flux1[id] = rhovx;
    flux2[id] = rhovy;
    flux1[id+size] = rhovx*vx+p;
    flux2[id+size] = rhovx*vy;
    flux1[id+size*2] = rhovx*vy;
    flux2[id+size*2] = rhovy*vy+p;
    flux1[id+size*3] = vx*(e+p);
    flux2[id+size*3] = vy*(e+p);
  }
}

__global__ void flux_vol (double* flux_vol1, double* flux_vol2, double* flux_quad1, double* flux_quad2, 
			   int nx, int ny, int mx, int my, int nvar){
  int id, icell, jcell, imod, jmod, var;
  int xq, yq, idq;
  int a = my;
  int b = mx*a;
  int c = ny*b;
  int d = nx*c;
  int size = nvar*d;
  double val1,val2;
  val1=val2=0.0;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  var = id/d;
  jmod = id - var*d;
  icell = jmod/c;
  jmod -= icell*c;
  jcell = jmod/b;
  jmod -= jcell*b;
  imod = jmod/a;
  jmod -= imod*a;
  idq = jcell*b + icell*c + var*d;
  
  if( id < size ){
    for( xq=0; xq < mx; xq++){
      for( yq=0; yq < my; yq++){
	
	val1 += 0.25*flux_quad1[yq +xq*a +idq]* 
	  legendre_prime(xquad[xq],imod,1)*wxquad[xq]*
	  legendre(yquad[yq],jmod,1)*wyquad[yq];
	  
	val2 +=  0.25*flux_quad2[yq +xq*a +idq]* 
	  legendre_prime(yquad[yq],jmod,1)*wyquad[yq]*
	  legendre(xquad[xq],imod,1)*wxquad[xq];
	
	/*val1 += 0.25*(flux_quad1[yq +xq*a +idq]*legendre_prime(xquad[xq],imod,1)*legendre(yquad[yq],jmod,1)+
	  flux_quad2[yq +xq*a +idq]*legendre_prime(yquad[yq],jmod,1)*legendre(xquad[xq],imod,1))*wxquad[xq]*wyquad[yq];*/
      }
    }
    flux_vol1[id] = val1;
    flux_vol2[id] = val2;
  }
}

__global__ void compute_max_speed(double* w, double* maxvalues, double gamma, int size){
  int id, idofmax, jump; 
  double speed,speed_max;
  __shared__ double maximums[BLOCK];
  __shared__ int ids[BLOCK];
  id = threadIdx.x;
  if(id < size){
    speed_max = compute_speed(w[id],w[id+size],w[id+size*2],w[id+size*3],gamma);
    idofmax = id;
    for (id = threadIdx.x+blockDim.x; id < size; id += blockDim.x){ //This is implemented considering only one block in the reduction launch.
      speed = compute_speed(w[id],w[id+size],w[id+size*2],w[id+size*3],gamma);
	if (speed > speed_max){
	  speed_max= speed;
	  idofmax = id;
	}
    }
    maximums[threadIdx.x] = speed_max;
    ids[threadIdx.x] = idofmax;
  }
  __syncthreads();
  for(jump = blockDim.x/2; jump > 0; jump >>= 1){
    if( threadIdx.x < jump ){
      if (maximums[threadIdx.x+jump] >=	maximums[threadIdx.x]){
	maximums[threadIdx.x] = maximums[threadIdx.x+jump];
	ids[threadIdx.x] = ids[threadIdx.x+jump];
      }
    }
    __syncthreads();
  }
  if(threadIdx.x == 0){
    maxvalues[0] = maximums[0]; //Max speed
    maxvalues[1] = w[ids[0]+size];
    maxvalues[2] = w[ids[0]+size*2];
    maxvalues[3] = sqrt(gamma*max(w[ids[0]+size*3],1E-10)/max(w[ids[0]],1E-10));
  }
  
}

__global__ void compute_AT(double* ul, double* ur, double* ub, double* ut,  double* delta_u, 
			   int my, int mx, int ny, int nx, int nvar){
  int id, var, icell, jcell, imod, jmod, yq, xq, idq, xid, yid;
  double shudl[5],shudr[5],shudb[5],shudt[5];
  double chsi_m = -1, chsi_p = 1, du;
  int a = my;
  int b = mx*a;
  int c = ny*b;
  int d = nx*c;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  var = id/nx*ny;
  jcell = id - var*nx*ny;
  icell = jcell/ny;
  jcell -= icell*ny;
  idq = jcell*b + icell*c + var*d;
  yid = jcell*my + icell*my*ny + var*my*ny*nx;
  xid = jcell*mx + icell*mx*ny + var*mx*ny*nx;
  if( id < nx*ny*nvar ){
    for (xq=0;xq<mx;xq++)
      shudb[xq] = shudt[xq] = 0;
    for (yq=0;yq<my;yq++)
      shudl[yq] = shudr[yq] = 0;
    for (imod=0;imod<mx;imod++){ //This is implemented considering only one block in the reduction launch.
      for (jmod=0;jmod<my;jmod++){
	du = delta_u[jmod+imod*a+idq];
	for (yq=0;yq<my;yq++){
	  shudl[yq] += du*legendre(chsi_m,imod,1)*legendre(yquad[yq],jmod,1);
	  shudr[yq] += du*legendre(chsi_p,imod,1)*legendre(yquad[yq],jmod,1);
	}
	for (xq=0;xq<mx;xq++){
	  shudb[xq] += du*legendre(chsi_m,jmod,1)*legendre(xquad[xq],imod,1);
	  shudt[xq] += du*legendre(chsi_p,jmod,1)*legendre(xquad[xq],imod,1);
	}
      }
    }
    for (yq=0;yq<my;yq++){
      ul[yq+yid]=shudl[yq];
      ur[yq+yid]=shudr[yq];
    }
    for (xq=0;xq<mx;xq++){
      ub[xq+xid]=shudb[xq];
      ut[xq+xid]=shudt[xq];
    }
  }
}

__global__ void compute_FG(double* um, double* up, double* wm, double* wp, double* fm, double* fp, double* FG, 
			   double gamma, int my, int mx, int ny, int nx, int nvar, int dim, int bc, int size){
  
  int id, wid, var, icell, jcell, mod, im, ip, pid, mid, face, a, b, fsize;
  double speed_m, speed_p, cmax;
  id = blockDim.x * blockIdx.x + threadIdx.x; 
  if(dim == 1){
    face = id/(my*ny);
    mod = id-face*my*ny;
    jcell = mod/my;
    mod -= jcell*my;
    wid = mod + jcell*my;
    b = nx*ny*my;
    a = ny*my;
    fsize = nx;
  }
  else if(dim == 0){

    icell = id/(mx*(ny+1));
    mod = id-icell*mx*(ny+1);
    face = mod/mx;
    mod -= face*mx;
    wid = mod + icell*mx*ny;
    b = nx*ny*mx;
    a = mx;
    fsize = ny;
  }
  if(id < size){
    im = face-1;
    ip = face;
    im = BC(im,fsize,bc);
    ip = BC(ip,fsize,bc);
    pid = wid+im*a;
    mid = wid+ip*a;
    speed_p = compute_speed(wp[pid],wp[pid+b],wp[pid+b*2],wp[pid+b*3],gamma);
    speed_m = compute_speed(wm[mid],wm[mid+b],wm[mid+b*2],wm[mid+b*3],gamma);
    cmax=max(speed_m,speed_p);
    for(var = 0; var < nvar; var++)
      FG[id+size*var]=0.5*(fp[pid+b*var]+fm[mid+b*var])+0.5*cmax*(up[pid+b*var]-um[mid+b*var]);
  }
}
__global__ void flux_line_integral(double* edge,double* F,double* G,int my,int mx,int ny,int nx,int nvar){
  int id, icell, jcell, imod, jmod, var;
  int xq, yq, idF, idG;
  int a = my;
  int b = mx*a;
  int c = ny*b;
  int d = nx*c;
  int size = nvar*d;
  double val1,val2;
  double chsi_m = -1, chsi_p = 1;
  val1=val2=0.0;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  var = id/d;
  jmod = id - var*d;
  icell = jmod/c;
  jmod -= icell*c;
  jcell = jmod/b;
  jmod -= jcell*b;
  imod = jmod/a;
  jmod -= imod*a;
  idF = jcell*my + var*my*ny*(nx+1);
  idG = icell*mx*(ny+1) + var*mx*nx*(ny+1);
  if( id < size){
    for(yq = 0; yq < my; yq++){
      val1 += 0.5*F[yq+(icell+1)*my*ny+idF]*legendre(chsi_p,imod,1)*legendre(yquad[yq],jmod,1)*wyquad[yq];
      val2 += 0.5*F[yq+icell*my*ny+idF]*legendre(chsi_m,imod,1)*legendre(yquad[yq],jmod,1)*wyquad[yq];
    }
    edge[id] = val1;
    edge[id+size] = val2;
    val1=val2=0.0;
    for(xq = 0; xq < mx; xq++){
      val1 += 0.5*G[xq+(jcell+1)*mx+idG]*legendre(chsi_p,jmod,1)*legendre(xquad[xq],imod,1)*wxquad[xq];
      val2 += 0.5*G[xq+jcell*mx+idG]*legendre(chsi_m,jmod,1)*legendre(xquad[xq],imod,1)*wxquad[xq];
    }
    edge[id+size*2] = val1;
    edge[id+size*3] = val2;
  }
}

__global__ void compute_dudt(double* dudt, double* flux_vol1, double* flux_vol2, double* edge, double oneoverdx, double oneoverdy, int size){
  int id;
  id = blockDim.x * blockIdx.x + threadIdx.x;
  if( id < size ){
    dudt[id] = 2*( oneoverdx*flux_vol1[id] + oneoverdy*flux_vol2[id] - oneoverdx*(edge[id]-edge[id+size]) - oneoverdy*(edge[id+size*2]-edge[id+size*3]));//+source_vol[id] 
    //dudt[id] =flux_vol1[id];
  }
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

extern "C" void device_get_modes_from_nodes_(double** nodes, double** modes){
  int size = nx*ny*mx*my*nvar;
  cudaError_t error = cudaGetLastError();
  get_modes_from_nodes<<<(size+BLOCK-1)/BLOCK,BLOCK>>>(*nodes,*modes,nx,ny,mx,my,nvar);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error get_modes_from_nodes: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void device_get_nodes_from_modes_(double** modes, double** nodes){
  int size = nx*ny*mx*my*nvar;
  cudaError_t error = cudaGetLastError();
  get_nodes_from_modes<<<(size+BLOCK-1)/BLOCK,BLOCK>>>(*modes,*nodes,nx,ny,mx,my,nvar);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error get_nodes_from_modes: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void device_compute_max_speed_ (double* csmax, double* vxmax, double* vymax, double* cmax){
  double max[4];
  cudaError_t error = cudaGetLastError();
  int size = nx*ny*mx*my;
  compute_primitive<<<(size+BLOCK-1)/BLOCK,BLOCK>>>( u, w, gmma, size);
  compute_max_speed<<<1,BLOCK>>>( w, maxvalues, gmma, size);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error compute max speed: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
  cudaMemcpy(max,maxvalues,4*sizeof(double),cudaMemcpyDeviceToHost);
  *cmax  = max[0];
  *vxmax = max[1];
  *vymax = max[2];
  *csmax = max[3];
 }

extern "C" void device_compute_update_(int* Iter, double* DT){
  double dt = *DT;
  int iter = *Iter;
  //printf("case %d\n",iter);
  switch (iter){
  case 0:
    modes = du;
    break;
  case 1:
    modes = w1;
    break;
  case 2:
    modes = w2;
    break;
  case 3:
    modes = w3;
    break;
  case 4:
    modes = w4;
    break;
  }
  cudaError_t error = cudaGetLastError();
  get_nodes_from_modes<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(modes,u_d_q,nx,ny,mx,my,nvar);
  compute_primitive<<<(usize+BLOCK-1)/BLOCK,BLOCK>>>( u_d_q, w, gmma, usize);
  compute_flux<<<(usize+BLOCK-1)/BLOCK,BLOCK>>>(u_d_q, w, flux_q1, flux_q2, usize);
  flux_vol<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(flux_v1,flux_v2,flux_q1,flux_q2,nx,ny,mx,my,nvar);
  compute_AT<<<(nx*ny*nvar+BLOCK-1)/BLOCK,BLOCK>>>(ul,ur,ub,ut,modes,my,mx,ny,nx,nvar);
  compute_primitive<<<(my*ny*nx+BLOCK-1)/BLOCK,BLOCK>>>( ul, wl, gmma, my*ny*nx);
  compute_primitive<<<(my*ny*nx+BLOCK-1)/BLOCK,BLOCK>>>( ur, wr, gmma, my*ny*nx);
  compute_primitive<<<(mx*ny*nx+BLOCK-1)/BLOCK,BLOCK>>>( ub, wb, gmma, mx*ny*nx);
  compute_primitive<<<(mx*ny*nx+BLOCK-1)/BLOCK,BLOCK>>>( ut, wt, gmma, mx*ny*nx);
  compute_flux<<<(nx*ny*my+BLOCK-1)/BLOCK,BLOCK>>>(ul, wl, fl1, fl2, nx*ny*my);
  compute_flux<<<(nx*ny*my+BLOCK-1)/BLOCK,BLOCK>>>(ur, wr, fr1, fr2, nx*ny*my);
  compute_flux<<<(nx*ny*mx+BLOCK-1)/BLOCK,BLOCK>>>(ub, wb, fb1, fb2, nx*ny*mx);
  compute_flux<<<(nx*ny*mx+BLOCK-1)/BLOCK,BLOCK>>>(ut, wt, ft1, ft2, nx*ny*mx);
  compute_FG<<<(my*ny*(nx+1)+BLOCK-1)/BLOCK,BLOCK>>>(ul,ur,wl,wr,fl1,fr1,F,gmma,my,mx,ny,nx,nvar,1,1,my*ny*(nx+1));
  compute_FG<<<(mx*(ny+1)*nx+BLOCK-1)/BLOCK,BLOCK>>>(ub,ut,wb,wt,fb2,ft2,G,gmma,my,mx,ny,nx,nvar,0,1,mx*(ny+1)*nx);
  flux_line_integral<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(edge,F,G,my,mx,ny,nx,nvar);
  compute_dudt<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(dudt,flux_v1,flux_v2,edge,invdx,invdy,nx*ny*mx*my*nvar);
  switch (iter){
  case 0:
    sum2<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(w1,du,dudt,(double)0.391752226571890*dt, tsize);
    break;
  case 1:
    sum3<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(w2,du,w1,dudt,0.444370493651235,0.555629506348765,0.368410593050371*dt, tsize);
    //get_nodes_from_modes<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(w1,u,nx,ny,mx,my,nvar);
    //cudaDeviceSynchronize();
    //cudaMemcpy(dudt,&edge[tsize*2],tsize*sizeof(double),cudaMemcpyDeviceToDevice);
    break;
  case 2:
    sum3<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(w3,du,w2,dudt,0.620101851488403,0.379898148511597,0.251891774271694*dt, tsize);
    break;
  case 3:
    sum3<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(w4,du,w3,dudt,0.178079954393132,0.821920045606868,0.544974750228521*dt, tsize);
    sum3<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(du,w2,w3,dudt,0.517231671970585,0.096059710526147,0.063692468666290*dt, tsize);
    break;
  case 4:
    plus_equal<<<(tsize+BLOCK-1)/BLOCK,BLOCK>>>(du,w4,dudt,0.386708617503269,0.226007483236906*dt, tsize);
    break;
    }
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error compute update: %s\n", cudaGetErrorString(error));
    exit(-1);
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

extern "C" void gpu_allocation_ (int *Nvar, int* Nx, int* Ny, int* Mx, int* My, double *Bl_x, double *Bl_y, double *CFL, double *Eta, double *Gamma) {
  nvar = *Nvar;
  nx = *Nx;
  ny = *Ny;
  mx = *Mx;
  my = *My;
  mods = max(mx,my);
  tsize = mx*my*nx*ny*nvar;
  usize = mx*my*nx*ny;
  boxlen_x = *Bl_x;
  boxlen_y = *Bl_y;
  dx = boxlen_x/double(nx);
  dy = boxlen_y/double(ny);
  invdx = 1/dx;
  invdy = 1/dy;
  //printf("dx %.14g dy %.14g invdx %.14g invdy %.14g\n",dx,dy,invdx,invdy);
  cfl = *CFL;
  gmma = *Gamma;
  eta = *Eta;
  cudaError_t error = cudaGetLastError();
  cudaMalloc ( &u, tsize * sizeof(double));
  cudaMalloc ( &u_eq, tsize * sizeof(double));
  cudaMalloc ( &u_d_q, tsize * sizeof(double));
  cudaMalloc ( &du, tsize * sizeof(double));
  cudaMalloc ( &w, tsize * sizeof(double));
  cudaMalloc ( &w1, tsize * sizeof(double));
  cudaMalloc ( &w2, tsize * sizeof(double));
  cudaMalloc ( &w3, tsize * sizeof(double));
  cudaMalloc ( &w4, tsize * sizeof(double));
  cudaMalloc ( &dudt, tsize * sizeof(double));
  cudaMalloc ( &ul, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &ur, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &ut, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &ub, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &wl, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &wr, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &wt, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &wb, nvar*nx*ny*mx * sizeof(double));
  cudaMalloc ( &flux_q1, tsize * sizeof(double)); 
  cudaMalloc ( &flux_q2, tsize * sizeof(double)); 
  cudaMalloc ( &flux_v1, tsize * sizeof(double)); 
  cudaMalloc ( &flux_v2, tsize * sizeof(double)); 
  cudaMalloc ( &fl1, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &fr1, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &ft1, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &fb1, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &fl2, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &fr2, nvar*nx*ny*my * sizeof(double)); 
  cudaMalloc ( &ft2, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &fb2, nvar*nx*ny*mx * sizeof(double)); 
  cudaMalloc ( &F, nvar*(nx+1)*ny*my * sizeof(double)); 
  cudaMalloc ( &G, nvar*nx*(ny+1)*mx * sizeof(double)); 
  cudaMalloc ( &edge, tsize*4 * sizeof(double));
  cudaMalloc ( &x, nx*ny*mx*my * sizeof(double));
  cudaMalloc ( &y, nx*ny*mx*my * sizeof(double));
  cudaMalloc ( &maxvalues, 4 * sizeof(double));

  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error gpu init: %s\n", cudaGetErrorString(error));
    exit(-1);
  }
}

extern "C" void gpu_set_pointers_ (double** u_d, double** du_d, double** dudt_d, double** w_d, double** u_eq_d,
			       double** x_d, double** y_d ,double* x_quad,
			       double* y_quad, double* w_x_quad,
				   double* w_y_quad, double* Sqrt_mod) {
  cudaError_t error = cudaGetLastError(); 
  *u_d = u;
  *du_d =du;
  *dudt_d =dudt;
  *w_d = w;
  *u_eq_d = u_eq;
  *x_d = x;
  *y_d = y;
  cudaMemcpyToSymbol(sqrt_mod,Sqrt_mod,sizeof(double)*mods);
  cudaMemcpyToSymbol(xquad,x_quad,sizeof(double)*mx);
  cudaMemcpyToSymbol(yquad,y_quad,sizeof(double)*my);
  cudaMemcpyToSymbol(wxquad,w_x_quad,sizeof(double)*mx);
  cudaMemcpyToSymbol(wyquad,w_y_quad,sizeof(double)*my);
  error = cudaGetLastError();
  if(error != cudaSuccess){
    printf("CUDA error gpu init: %s\n", cudaGetErrorString(error));
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
