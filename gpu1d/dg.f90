program dg
  USE ISO_C_BINDING
  implicit none
  integer,parameter::n=4
  integer,parameter::nx=512
  integer,parameter::nquad=n
  integer,parameter::ninit=3
  logical,parameter::use_limiter=.true.
  logical,parameter::use_extremum=.false.
  Character(LEN=3),parameter::integrator='RK4'
  real(kind=8),dimension(1:nx,1:n)::rho,drhodt,w1,w2,w3,w4
  real(kind=8),dimension(1:nquad)::chsi_quad,w_quad
  integer::iter,icell,i,j
  real(kind=8)::legendre,legendre_prime
  real(kind=8)::rho_init,xcell,dx,x_quad
  real(kind=8)::t,dt,tend,velocity
  character(len=20)::filename

  !Device arrays
  
  type (c_ptr) :: rho_d, drhodt_d, w_quad_d, chsi_quad_d
  type (c_ptr) :: w1_d, w2_d, w3_d, w4_d
  type (c_ptr) :: rho_quad_d,flux_quad_d,flux_vol_d
  type (c_ptr) :: rho_left_d,rho_right_d, flux_left_d
  type (c_ptr) :: flux_right_d,flux_face_d

  !Allocation over device
  !call devices();
  call setdevice(0);
  call malloc_gpu(rho_d,nx*n)
  call malloc_gpu(drhodt_d,nx*n)
  call malloc_gpu(w1_d,nx*n)
  call malloc_gpu(w2_d,nx*n)
  call malloc_gpu(w3_d,nx*n)
  call malloc_gpu(w4_d,nx*n)
  call malloc_gpu(flux_vol_d,nx*n)
  call malloc_gpu(rho_left_d,nx)
  call malloc_gpu(flux_left_d,nx)
  call malloc_gpu(rho_right_d,nx)
  call malloc_gpu(flux_right_d,nx)
  call malloc_gpu(flux_face_d,nx+1)
  call malloc_gpu(chsi_quad_d,nquad)
  call malloc_gpu(w_quad_d,nquad)

  call gl_quadrature(chsi_quad,w_quad,nquad)
  call h2d(chsi_quad,chsi_quad_d,nquad)
  call h2d(w_quad,w_quad_d,nquad)
  
  dx=1.0/dble(nx)
  velocity=1.0
  tend=10.0

  ! Compute initial conditions
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     ! Loop over modes coefficients
     do i=1,n
        rho(icell,i)=0.0
        ! Loop over quadrature points
        do j=1,nquad
           ! Quadrature point in physical space
           x_quad=xcell+dx/2.0*chsi_quad(j)
           ! Perform integration using GL quadrature
           rho(icell,i)=rho(icell,i)+0.5* &
                & rho_init(x_quad,ninit)* &
                & legendre(chsi_quad(j),i-1)* &
                & w_quad(j)
        end do
        !write(*,*)rho(icell,i)
     end do
  end do
  
  call h2d(rho,rho_d,nx*n)
  write(*,*)'Initial conditions'
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(*,'(7(1PE10.3,1X))')xcell,(rho(icell,i),i=1,n)
  end do
    
  
  ! Main time loop
  t=0
  dt=0.5*dx/abs(velocity)/(2.0*dble(n)+1.0)
  iter=0
  write(filename,"(A5,I5.5)")"dens_",iter/10
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(10,'(7(1PE12.5,1X))')xcell,(rho(icell,i),i=1,n)
  end do
  close(10)

  do while(t < tend)

     if(integrator=='RK1')then
        call compute_update(rho,drhodt,nx,n,chsi_quad,w_quad,nquad)
        rho=rho+dt*drhodt
     endif

     if(integrator=='RK2')then
        call compute_update(rho,drhodt,nx,n,chsi_quad,w_quad,nquad)
        w1=rho+dt*drhodt
        if(use_limiter)call limiter(w1,nx,n,use_extremum)
        call compute_update(w1,drhodt,nx,n,chsi_quad,w_quad,nquad)
        rho=0.5*rho+0.5*w1+0.5*dt*drhodt
        if(use_limiter)call limiter(rho,nx,n,use_extremum)
     endif

     if(integrator=='RK3')then
        call compute_update(rho,drhodt,nx,n,chsi_quad,w_quad,nquad)
        w1=rho+dt*drhodt
        if(use_limiter)call limiter(w1,nx,n,use_extremum)
        call compute_update(w1,drhodt,nx,n,chsi_quad,w_quad,nquad)
        w2=0.75*rho+0.25*w1+0.25*dt*drhodt
        if(use_limiter)call limiter(w2,nx,n,use_extremum)
        call compute_update(w2,drhodt,nx,n,chsi_quad,w_quad,nquad)
        rho=1.0/3.0*rho+2.0/3.0*w2+2.0/3.0*dt*drhodt
        if(use_limiter)call limiter(rho,nx,n,use_extremum)
     endif

     if(integrator=='RK4')then
        call device_compute_update(rho_d,drhodt_d,flux_vol_d,rho_right_d,rho_left_d,flux_right_d, &
             & flux_left_d,flux_face_d,nx,n,chsi_quad_d,w_quad_d,nquad)
        
        call device_sum2(w1_d, rho_d, drhodt_d,dble(0.391752226571890)*dt,nx*n)
        if(use_limiter)call device_limiter(w1_d,flux_vol_d,nx,n,use_extremum)
        call device_compute_update(w1_d,drhodt_d,flux_vol_d,rho_right_d,rho_left_d,flux_right_d, &
             & flux_left_d,flux_face_d,nx,n,chsi_quad_d,w_quad_d,nquad)
        
        call device_sum3(w2_d,rho_d,w1_d,drhodt_d,dble(0.444370493651235),dble(0.555629506348765),dble(0.368410593050371)*dt,nx*n)
        if(use_limiter)call device_limiter(w2_d,flux_vol_d,nx,n,use_extremum)
        call device_compute_update(w2_d,drhodt_d,flux_vol_d,rho_right_d,rho_left_d,flux_right_d, &
             & flux_left_d,flux_face_d,nx,n,chsi_quad_d,w_quad_d,nquad)
        
        call device_sum3(w3_d,rho_d,w2_d,drhodt_d,dble(0.620101851488403),dble(0.379898148511597),dble(0.251891774271694)*dt,nx*n)
        if(use_limiter)call device_limiter(w3_d,flux_vol_d,nx,n,use_extremum)
        call device_compute_update(w3_d,drhodt_d,flux_vol_d,rho_right_d,rho_left_d,flux_right_d, & 
             & flux_left_d,flux_face_d,nx,n,chsi_quad_d,w_quad_d,nquad)
        
        call device_sum3(w4_d,rho_d,w3_d,drhodt_d,dble(0.178079954393132),dble(0.821920045606868),dble(0.544974750228521)*dt,nx*n)
        if(use_limiter)call device_limiter(w4_d,flux_vol_d,nx,n,use_extremum)
        call device_sum3(rho_d,w2_d,w3_d,drhodt_d,dble(0.517231671970585),dble(0.096059710526147),dble(0.063692468666290)*dt,nx*n)
        call device_compute_update(w4_d,drhodt_d,flux_vol_d,rho_right_d,rho_left_d,flux_right_d, &
             & flux_left_d,flux_face_d,nx,n,chsi_quad_d,w_quad_d,nquad)
        call device_plus_equal(rho_d,w4_d,drhodt_d,dble(0.386708617503269),dble(0.226007483236906)*dt,nx*n)
        if(use_limiter)call device_limiter(rho_d,flux_vol_d,nx,n,use_extremum)
                
     endif
      
     t=t+dt
     iter=iter+1
     write(*,*)'time=',iter,t,dt

!!$     if(mod(iter,10)==0)then
!!$        write(filename,"(A5,I5.5)")"dens_",iter/10
!!$        open(10,file=TRIM(filename)//".dat",form='formatted')
!!$        do icell=1,nx
!!$           xcell=(dble(icell)-0.5)*dx
!!$           write(10,'(7(1PE12.5,1X))')xcell,(rho(icell,i),i=1,n)
!!$        end do
!!$        close(10)
!!$     endif
  enddo
  call d2h(rho_d,rho,nx*n)  
  ! Write final state
  write(filename,"(A5,I5.5)")"dens_",99999
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(10,'(7(1PE12.5,1X))')xcell,(rho(icell,i),i=1,n)
  end do
  close(10)

  write(*,*)'========================================'
  write(*,*)'time=',t,dt
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(*,'(7(1PE12.5,1X))')xcell,(rho(icell,i),i=1,n)
  end do

end program dg


subroutine limiter(rho,nx,n,use_extremum)
  implicit none
  integer::n,nx
  logical::use_extremum
  real(kind=8),dimension(1:nx,1:n)::rho
  real(kind=8),dimension(1:nx,1:n-1)::rho_lim
  real(kind=8)::rhoL,rhoM,rhoR
  real(kind=8)::minmod,maxmod
  real(kind=8)::rho_min,rho_max
  real(kind=8)::coeff_i,coeff_ip1,coeff,D2rho
  integer::icell,i,j,iface,ileft,iright
  if(n==1)return
  ! Compute classical minmod limiter
  rho_lim=rho(1:nx,2:n)
  do icell=1,nx
     ileft=icell-1
     iright=icell+1
     if(icell==1)ileft=nx
     if(icell==nx)iright=1
     ! Loop only over high-order modes
     do i=n-1,1,-1
        ! Renormalise to get proper Legendre polynomials
        ! and corresponding derivatives
        coeff_i=sqrt(2.0*dble(i-1)+1.0)
        coeff_ip1=sqrt(2.0*dble(i)+1.0)*(2.0*dble(i)-1)
        rhoL=(rho(icell,i)-rho(ileft,i))*coeff_i
        rhoR=(rho(iright,i)-rho(icell,i))*coeff_i
        rhoM=rho(icell,i+1)*coeff_ip1
        rho_min=minmod(rhoL,rhoM,rhoR)/coeff_ip1
        ! If smooth extremum, no limiter applied
        if(use_extremum.AND.rho_min.NE.rho(icell,i+1).AND.&
             & min(rhoL*rhoR,rho(ileft,i+1)*rho(iright,i+1))<0.0)then
           rhoL=2.0*(rho(icell,i+1)-rho(ileft,i+1))
           rhoR=2.0*(rho(iright,i+1)-rho(icell,i+1))
           rhoM=0.5*(rho(iright,i+1)-rho(ileft,i+1))
           D2rho=minmod(rhoL,rhoM,rhoR)
           coeff=1.0
           if(rhoM.NE.0.0)coeff=D2rho/rhoM
           rho_lim(icell,i)=rho_min+(rho(icell,i+1)-rho_min)*coeff
        else
           rho_lim(icell,i)=rho_min
        endif
        ! If no limiter applied, exit loop over modes
        if(rho_lim(icell,i).EQ.rho(icell,i+1)) exit
     end do
  end do
  rho(1:nx,2:n)=rho_lim
end subroutine limiter

function minmod(x,y,z)
  implicit none
  real(kind=8)::x,y,z,s,minmod
  s=sign(1d0,x)
  if(sign(1d0,y)==s.AND.sign(1d0,z)==s)then
     minmod=s*min(abs(x),abs(y),abs(z))
  else
     minmod=0.0
  endif
  return
end function minmod

subroutine compute_update(rho,drhodt,nx,n,chsi_quad,w_quad,nquad)
  
  implicit none
  integer::n,nx,nquad
  real(kind=8),dimension(1:nx,1:n)::rho,drhodt
  real(kind=8),dimension(1:nquad)::chsi_quad,w_quad

  real(kind=8),dimension(1:nquad)::rho_quad,flux_quad
  real(kind=8),dimension(1:nx,1:n)::flux_vol
  real(kind=8),dimension(1:nx)::rho_left,rho_right
  real(kind=8),dimension(1:nx)::flux_left,flux_right
  real(kind=8),dimension(1:nx+1)::flux_face
  
  integer::icell,i,j,iface,ileft,iright
  real(kind=8)::legendre,legendre_prime
  real(kind=8)::chsi_left=-1,chsi_right=+1
  real(kind=8)::dx,x_quad
  real(kind=8)::cmax,oneoverdx,velocity
  
  dx=1.0/dble(nx)
  oneoverdx=1.0/dx
  velocity=1.0

  ! Loop over cells
  do icell=1,nx

     !==================================
     ! Compute flux at quadrature points
     !==================================
     ! Loop over quadrature points
     do j=1,nquad
        rho_quad(j)=0.0
        ! Loop over modes
        do i=1,n
           rho_quad(j)=rho_quad(j)+rho(icell,i)*legendre(chsi_quad(j),i-1)
        end do
        ! Compute flux at quadrature points
        flux_quad(j)=velocity*rho_quad(j)
     end do
     
     !================================
     ! Compute volume integral DG term
     !================================
     ! Loop over modes
     do i=1,n
        flux_vol(icell,i)=0.0
        ! Loop over quadrature points
        do j=1,nquad
           flux_vol(icell,i)=flux_vol(icell,i)+ &
                & flux_quad(j)* &
                & legendre_prime(chsi_quad(j),i-1)* &
                & w_quad(j)
        end do
     end do
     
     !==============================
     ! Compute left and right fluxes
     !==============================
     rho_left(icell)=0.0
     rho_right(icell)=0.0
     ! Loop over modes
     do i=1,n
        rho_left(icell)=rho_left(icell)+rho(icell,i)*legendre(chsi_left,i-1)
        rho_right(icell)=rho_right(icell)+rho(icell,i)*legendre(chsi_right,i-1)
     end do
     ! Compute flux at left and right points
     flux_left(icell)=velocity*rho_left(icell)
     flux_right(icell)=velocity*rho_right(icell)
  end do
  ! End loop over cells
  
  !====================================
  ! Compute Lax Friedrich physical flux
  !====================================
  ! Loop over faces
  do iface=1,nx+1
     ileft=iface-1
     iright=iface
     ! Periodic boundary conditions
     if(iface==1)ileft=nx
     if(iface==nx+1)iright=1
     ! Maximum wave speed
     cmax=abs(velocity)
     flux_face(iface)=0.5*(flux_right(ileft)+flux_left(iright)) &
          & -0.5*cmax*(rho_left(iright)-rho_right(ileft))
  end do
  
  !========================
  ! Compute final DG update
  !========================
  ! Loop over cells
  do icell=1,nx
     ! Loop over modes
     do i=1,n
        drhodt(icell,i)=oneoverdx*(flux_vol(icell,i) &
             & -(flux_face(icell+1)*legendre(chsi_right,i-1) &
             & -flux_face(icell)*legendre(chsi_left,i-1)))
     end do
  end do
  
end subroutine compute_update

function rho_init(x,n)
  integer::n
  real(kind=8)::x,rho_init
  real(kind=8)::dpi=acos(-1d0)
  select case (n)
     case(1)
        rho_init=sin(2.0*dpi*x)
     case(2)
        if(abs(x-0.5)<0.25)then
           rho_init=1.
        else
           rho_init=-1.
        endif
     case(3)
        rho_init=exp(-(x-0.25)**2/2.0/0.05**2)
        if(abs(x-0.7)<0.1)then
           rho_init=rho_init+1.
        endif
     end select
  return
end function rho_init

subroutine gl_quadrature(x_quad,w_quad,n)
  integer::n
  real(kind=8),dimension(1:n)::x_quad,w_quad

  integer::i,iter
  real(kind=8)::dpi=acos(-1.0d0),xx
  real(kind=8)::legendre,legendre_prime

  write(*,*)"Computing Gauss-Legendre quadrature points and weights."
  do i=1,n
     xx=(1.0-0.125/n/n+0.125/n/n/n)* & 
          & cos(dpi*(4.0*dble(i)-1.0)/(4.0*dble(n)+2.0))
     do iter=1,50
        xx=xx-legendre(xx,n)/legendre_prime(xx,n)
     end do
     x_quad(i)=xx
     w_quad(i)=2.0*(2.0*dble(n)+1.0)/(1.0-x_quad(i)**2) &
          & /legendre_prime(x_quad(i),n)**2
  end do
  do i=n/2+1,n
     x_quad(i)=-x_quad(n-i+1)
     w_quad(i)=w_quad(n-i+1)
  end do
!!$  do i=1,n
!!$     write(*,*)i,x_quad(i),w_quad(i)
!!$  end do

end subroutine gl_quadrature

