program dg
  use dg_commons
  implicit none
  !==============================================
  ! This program solves the 1D Euler equations
  ! using the Discontinuous Galerkin method
  ! up to 4th order
  ! (using Legendre polynomials up to x^3).
  ! Romain Teyssier, November 1st 2015
  !==============================================

  ! Main variables
  real(kind=8),dimension(1:nvar,1:n,1:nx)::u,dudt,w1,w2,w3,w4, ureal, uinit

  integer::iter=0,icell,i,j,ivar, intnode
  real(kind=8)::legendre,legendre_prime
  real(kind=8)::xcell,dx,x_quad
  real(kind=8)::t,dt,cmax
  real(kind=8),dimension(1:nvar)::uu,ww
  character(len=20)::filename

  !=====================================================
  ! Compute Gauss-Legendre quadrature points and weights
  !=====================================================
  call gl_quadrature(chsi_quad,w_quad,nquad)

  ! Mesh size
  dx=boxlen/dble(nx)

  !===========================
  ! Compute initial conditions
  !===========================
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     ! Loop over modes coefficients
     do i=1,n
        u(1:nvar,i,icell)=0.0
        ! Loop over quadrature points
        do j=1,nquad
           ! Quadrature point in physical space
           x_quad=xcell+dx/2.0*chsi_quad(j)
           call condinit(x_quad,uu)
           ! Perform integration using GL quadrature
           u(1:nvar,i,icell)=u(1:nvar,i,icell)+0.5* &
                & uu(1:nvar)* &
                & legendre(chsi_quad(j),i-1)* &
                & w_quad(j)
        end do
     end do
  end do

  write(filename,"(A5,I5.5)")"modes",iter/10
  open(7,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(7,'(7(1PE12.5,1X))')xcell,(u(ivar,1,icell),ivar=1,nvar)
  end do
  close(7)

  ureal(:,:,:)=0
  uinit(:,:,:)=0
  ! reconstruct u from modes
  do ivar = 1,nvar
    do icell=1,nx
       xcell=(dble(icell)-0.5)*dx
       do i=1,n
         x_quad=xcell+dx/2.0*chsi_quad(i)
         call condinit(x_quad,uu)
         uinit(ivar,i,icell) = uu(ivar)
         do intnode = 1,n
        ! Loop over quadrature points
          ureal(ivar,i,icell) = ureal(ivar, i, icell) +&
          & u(ivar, intnode ,icell)*legendre(x_quad,intnode-1)
         end do
        end do
      end do
    end do

  write(filename,"(A5,I5.5)")"fullo",iter/10
  open(7,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(7,'(7(1PE12.5,1X))')xcell,(ureal(ivar,1,icell),ivar=1,nvar)
  end do
  close(7)

  write(filename,"(A5,I5.5)")"initi",iter/10
  open(7,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     write(7,'(7(1PE12.5,1X))')xcell,(uinit(ivar,1,icell),ivar=1,nvar)
  end do
  close(7)


  !=================================
  ! Write initial conditions to file
  !=================================
  write(filename,"(A5,I5.5)")"hydroman",iter/10
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     !xcell=(dble(icell)-0.5)*dx
     xcell = (dble(icell)-0.5)*dx + dx/2.0*chsi_quad(1)
     call compute_primitive(uinit(1:nvar,1,icell),ww,gamma,nvar)
     write(10,'(7(1PE12.5,1X))')xcell,(ww(ivar)-exp(-xcell),ivar=3,nvar)
  end do
  close(10)

  !===============
  ! Main  loop
  !===============
  t=0
  iter=0
  do while(t < tend)

     ! Compute time step
     call compute_max_speed(u,cmax)
     dt=0.8*dx/cmax/(2.0*dble(n)+1.0)

     if(integrator=='RK1')then
        call compute_update(u,dudt)
        u=u+dt*dudt
     endif

     if(integrator=='RK2')then
        call compute_update(u,dudt)
        w1=u+dt*dudt
        call limiter(w1)
        call compute_update(w1,dudt)
        u=0.5*u+0.5*w1+0.5*dt*dudt
        call limiter(u)
     endif

     if(integrator=='RK3')then
        call compute_update(u,dudt)
        w1=u+dt*dudt
        call limiter(w1)
        call compute_update(w1,dudt)
        w2=0.75*u+0.25*w1+0.25*dt*dudt
        call limiter(w2)
        call compute_update(w2,dudt)
        u=1.0/3.0*u+2.0/3.0*w2+2.0/3.0*dt*dudt
        call limiter(u)
     endif

     if(integrator=='RK4')then
        call compute_update(u,dudt)
        w1=u+0.391752226571890*dt*dudt
        call limiter(w1)
        call compute_update(w1,dudt)
        w2=0.444370493651235*u+0.555629506348765*w1+0.368410593050371*dt*dudt
        call limiter(w2)
        call compute_update(w2,dudt)
        w3=0.620101851488403*u+0.379898148511597*w2+0.251891774271694*dt*dudt
        call limiter(w3)
        call compute_update(w3,dudt)
        w4=0.178079954393132*u+0.821920045606868*w3+0.544974750228521*dt*dudt
        u=0.517231671970585*w2+0.096059710526147*w3+0.063692468666290*dt*dudt
        call limiter(w4)
        call compute_update(w4,dudt)
        u=u+0.386708617503269*w4+0.226007483236906*dt*dudt
        call limiter(u)
     endif

     t=t+dt
     iter=iter+1
     write(*,*)'time=',iter,t,dt

  enddo


  ureal(:,:,:)=0
  uinit(:,:,:)=0
  ! reconstruct u from modes
  do ivar = 1,nvar
    do icell=1,nx
       xcell=(dble(icell)-0.5)*dx
       do i=1,n
         x_quad=xcell+dx/2.0*chsi_quad(i)
         call condinit(x_quad,uu)
         uinit(ivar,i,icell) = uu(ivar)
         do intnode = 1,n
        ! Loop over quadrature points
          ureal(ivar,i,icell) = ureal(ivar, i, icell) +&
          & u(ivar, intnode ,icell)*legendre(x_quad,intnode-1)
         end do
        end do
      end do
    end do


!==========================
! Write final state to file (reconstructed)
!==========================
  write(filename,"(A5,I5.5)")"recon",99999
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     call compute_primitive(ureal(1:nvar,1,icell),ww,gamma,nvar)
     write(10,'(7(1PE12.5,1X))')xcell,(ww(ivar),ivar=1,nvar)
  end do
  close(7)

!  write(filename,"(A5,I5.5)")"recon",99999
!  open(10,file=TRIM(filename)//".dat",form='formatted')
!  do icell=1,nx
!     xcell=(dble(icell)-0.5)*dx
!     call compute_primitive(u(1:nvar,1,icell),ww,gamma,nvar)
!     write(10,'(7(1PE12.5,1X))')xcell,(ww(ivar),ivar=1,nvar)
!  end do
!  close(10)


  !==========================
  ! Write final state to file
  !==========================
  write(filename,"(A5,I5.5)")"hydro",99999
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     !xcell=(dble(icell)-0.5)*dx
     xcell = (dble(icell)-0.5)*dx + dx/2.0*chsi_quad(1)
     call compute_primitive(ureal(1:nvar,1,icell),ww,gamma,nvar)
     write(10,'(7(1PE12.5,1X))')xcell,(ww(ivar)-exp(-xcell),ivar=3,nvar)
  end do
  close(10)

  write(*,*)'========================================'
  write(*,*)'time=',t,dt
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     call compute_primitive(u(1:nvar,1,icell),ww,gamma,nvar)
     write(*,'(7(1PE12.5,1X))')xcell,(ww(ivar),ivar=1,nvar)
  end do

end program dg
!==============================================
subroutine limiter(u)
  use dg_commons
  implicit none
  real(kind=8),dimension(1:nvar,1:n,1:nx)::u
  !================================================================
  ! This routine applies the moment limiter to the current solution
  ! as in Krivodonova, 2007, JCP, 226, 879
  ! using characteristic variables.
  !================================================================
  real(kind=8),dimension(1:nvar,1:n,1:nx)::u_lim
  real(kind=8),dimension(1:nvar,1:n)::w_lim
  real(kind=8),dimension(1:nvar,1:n)::wL,wM,wR
  real(kind=8),dimension(1:nvar)::w
  real(kind=8),dimension(1:nvar)::uL,uM,uR
  real(kind=8),dimension(1:nvar)::u_left,u_right,w_left,w_right
  real(kind=8)::minmod,maxmod
  real(kind=8)::switch_left,switch_right
  real(kind=8)::u_min,u_max,w_min
  real(kind=8)::coeff_i,coeff_ip1,coeff,D2u
  integer::icell,i,j,iface,ileft,iright,ivar
  if(n==1)return
  ! Compute classical minmod limiter
  u_lim=u
  if(use_limiter)then
  do icell=1,nx
     ileft=icell-1
     iright=icell+1
     switch_left=1.0
     switch_right=1.0
     if(bc==1)then
        if(icell==1)ileft=nx
        if(icell==nx)iright=1
     endif
     if(bc==2)then
        if(icell==1)ileft=1
        if(icell==nx)iright=nx
     endif
     if(bc==3)then
        if(icell==1)then
           ileft=1
           switch_left=-1.0
        endif
        if(icell==nx)then
           iright=nx
           switch_right=-1.0
        endif
     endif
     ! Compute primitive variable for all modes
     call compute_primitive(u(1:nvar,1,icell),w,gamma,nvar)
     ! Loop only over high-order modes
     do i=n-1,1,-1
        ! Renormalise to get proper Legendre polynomials
        ! and corresponding derivatives
        coeff_i=sqrt(2.0*dble(i-1)+1.0)*(2.0*dble(i)-1)
        coeff_ip1=sqrt(2.0*dble(i)+1.0)*(2.0*dble(i)-1)
        uL(1:nvar)=(u(1:nvar,i,icell)-u(1:nvar,i,ileft))*coeff_i/coeff_ip1
        uR(1:nvar)=(u(1:nvar,i,iright)-u(1:nvar,i,icell))*coeff_i/coeff_ip1
        uM(1:nvar)=u(1:nvar,i+1,icell)
        uL(2)=switch_left*uL(2)
        uR(2)=switch_right*uR(2)
        call cons_to_char(uL,wL(1:nvar,i+1),w,nvar,gamma)
        call cons_to_char(uR,wR(1:nvar,i+1),w,nvar,gamma)
        call cons_to_char(uM,wM(1:nvar,i+1),w,nvar,gamma)
     end do
     w_lim=wM
     ! Loop over variables
     do ivar=1,nvar
        ! Loop only over high-order modes
        do i=n-1,1,-1
           w_min=minmod(wL(ivar,i+1),wM(ivar,i+1),wR(ivar,i+1))
           w_lim(ivar,i+1)=w_min
           if(ABS(w_min-wM(ivar,i+1)).LT.0.01*ABS(wM(ivar,i+1)))exit
        end do
        ! End loop over modes
     end do
     ! End loop over variables
     ! Compute conservative variables
     ! Loop only over high-order modes
     do i=n-1,1,-1
        call char_to_cons(w_lim(1:nvar,i+1),u_lim(1:nvar,i+1,icell),w,nvar,gamma)
     end do
  end do
  ! End loop over cells
  endif

  ! Check for unphysical values in the limited states
  do icell=1,nx
     ! Compute primitive variable
     call compute_primitive(u_lim(1:nvar,1,icell),w,gamma,nvar)
     u_left(1:nvar)=0.0; u_right(1:nvar)=0.0
     ! Loop over modes
     do i=1,n
        u_left(1:nvar)=u_left(1:nvar)+u_lim(1:nvar,i,icell)*(-1.0)**(i-1)*sqrt(2.0*dble(i)-1.0)
        u_right(1:nvar)=u_right(1:nvar)+u_lim(1:nvar,i,icell)*sqrt(2.0*dble(i)-1.0)
     end do
     call cons_to_prim(u_left,w_left,w,nvar,gamma)
     call cons_to_prim(u_right,w_right,w,nvar,gamma)
     if(w_left(1)<1d-10.OR.w_right(1)<1d-10.OR.w_left(3)<1d-10.OR.w_left(3)<1d-10)then
        u_lim(1:nvar,2:n,icell)=0.0
     endif
  end do

  ! Update variables with limited states
  u = u_lim

end subroutine limiter
!==============================================
subroutine compute_update(u,dudt)
  use dg_commons
  implicit none
  real(kind=8),dimension(1:nvar,1:n,1:nx)::u,dudt
  !===========================================================
  ! This routine computes the DG update for the input state u.
  !===========================================================
  real(kind=8),dimension(1:nvar,1:nquad)::u_quad,flux_quad, source_quad
  real(kind=8),dimension(1:nvar,1:n,1:nx)::flux_vol, source_vol
  real(kind=8),dimension(1:nvar,1:nx)::u_left,u_right
  real(kind=8),dimension(1:nvar)::flux_riemann,u_tmp
  real(kind=8),dimension(1:nvar,1:nx+1)::flux_face

  integer::icell,i,j,iface,ileft,iright,ivar
  real(kind=8)::legendre,legendre_prime
  real(kind=8)::chsi_left=-1,chsi_right=+1
  real(kind=8)::dx,x_quad
  real(kind=8)::cmax,oneoverdx,c_left,c_right

  dx=boxlen/dble(nx)
  oneoverdx=1.0/dx

  call gl_quadrature(chsi_quad,w_quad,nquad)


  ! Loop over cells
  do icell=1,nx

     !==================================
     ! Compute flux at quadrature points
     !==================================
     ! Loop over quadrature points
     do j=1,nquad
        u_quad(1:nvar,j)=0.0
        ! Loop over modes
        do i=1,n
           u_quad(1:nvar,j)=u_quad(1:nvar,j)+u(1:nvar,i,icell)*legendre(chsi_quad(j),i-1)
        end do
        ! Compute flux at quadrature points
        call compute_flux(u_quad(1:nvar,j),flux_quad(1:nvar,j),gamma,nvar)
     end do

     !================================
     ! Compute volume integral DG term
     !================================
     ! Loop over modes
     do i=1,n
        flux_vol(1:nvar,i,icell)=0.0
        ! Loop over quadrature points
        do j=1,nquad
           flux_vol(1:nvar,i,icell)=flux_vol(1:nvar,i,icell)+ &
                & flux_quad(1:nvar,j)* &
                & legendre_prime(chsi_quad(j),i-1)* &
                & w_quad(j)
        end do
     end do

     !==============================
     ! Compute left and right states
     !==============================
     u_left(1:nvar,icell)=0.0
     u_right(1:nvar,icell)=0.0
     ! Loop over modes
     do i=1,n
        u_left(1:nvar,icell)=u_left(1:nvar,icell)+u(1:nvar,i,icell)*legendre(chsi_left,i-1)
        u_right(1:nvar,icell)=u_right(1:nvar,icell)+u(1:nvar,i,icell)*legendre(chsi_right,i-1)
     end do
  end do
  ! End loop over cells

  !==========================================
  ! Compute physical flux from Riemann solver
  !==========================================
  ! Loop over faces
  do iface=1,nx+1
     ileft=iface-1
     iright=iface
     if(bc==1)then
        ! Periodic boundary conditions
        if(iface==1)ileft=nx
        if(iface==nx+1)iright=1
     endif
     if(bc==2.or.bc==3)then
        ! Zero gradient boundary conditions
        if(iface==1)ileft=1
        if(iface==nx+1)iright=nx
     endif
     ! Compute physical flux using Riemann solver
     select case(riemann)
     case(1)
        call riemann_llf(u_right(1:nvar,ileft),u_left(1:nvar,iright)&
             & ,flux_riemann,gamma,nvar)
     case(2)
        call riemann_hllc(u_right(1:nvar,ileft),u_left(1:nvar,iright)&
             & ,flux_riemann,gamma,nvar)
     end select
     flux_face(1:nvar,iface)=flux_riemann
     ! Compute boundary flux for reflexive BC
     if(bc==3.AND.iface==1)then
        u_tmp(1:nvar)=u_left(1:nvar,iright)
        u_tmp(2)=-u_tmp(2)
        select case(riemann)
        case(1)
           call riemann_llf(u_tmp,u_left(1:nvar,iright)&
                & ,flux_riemann,gamma,nvar)
        case(2)
           call riemann_hllc(u_tmp,u_left(1:nvar,iright)&
                & ,flux_riemann,gamma,nvar)
        end select
        flux_face(1:nvar,iface)=flux_riemann
     endif
     if(bc==3.AND.iface==nx+1)then
        u_tmp(1:nvar)=u_right(1:nvar,ileft)
        u_tmp(2)=-u_tmp(2)
        select case(riemann)
        case(1)
           call riemann_llf(u_right(1:nvar,ileft),u_tmp&
           &,flux_riemann,gamma,nvar)
        case(2)
           call riemann_hllc(u_right(1:nvar,ileft),u_tmp&
                &,flux_riemann,gamma,nvar)
        end select
        flux_face(1:nvar,iface)=flux_riemann
     endif
     if(bc==4.AND.iface==1)then
        u_tmp(1:nvar)= exp(0.)+u_right(1:nvar,1)-exp(-((0.5*dx)+dx/2.0*chsi_quad(1)))
        !u_tmp(2)=-u_tmp(2)
        select case(riemann)
        case(1)
           call riemann_llf(u_tmp,u_left(1:nvar,iright)&
                & ,flux_riemann,gamma,nvar)
        case(2)
           call riemann_hllc(u_tmp,u_left(1:nvar,iright)&
                & ,flux_riemann,gamma,nvar)
        end select
        flux_face(1:nvar,iface)=flux_riemann
     endif
     if(bc==4.AND.iface==nx+1)then
        u_tmp(1:nvar)=exp(-boxlen) + u_left(1:nvar,nx) - exp(-(((nx-0.5)*dx)+dx/2.0*chsi_quad(nquad)))
        !u_tmp(2)=-u_tmp(2)
        select case(riemann)
        case(1)
           call riemann_llf(u_right(1:nvar,ileft),u_tmp&
           &,flux_riemann,gamma,nvar)
        case(2)
           call riemann_hllc(u_right(1:nvar,ileft),u_tmp&
                &,flux_riemann,gamma,nvar)
        end select
        flux_face(1:nvar,iface)=flux_riemann
     endif

  end do

  select case (source)
     case(1)
       source_vol(:,:,:) = 0.0
     case(2)
       ! Compute source
       do icell=1,nx
          !==================================
          ! Compute source at quadrature points
          !==================================
          ! Loop over quadrature points
          do j=1,nquad
             u_quad(1:nvar,j)=0.0
             ! Loop over modes
             do i=1,n
                u_quad(1:nvar,j)=u_quad(1:nvar,j)+u(1:nvar,i,icell)*legendre(chsi_quad(j),i-1)
             end do
             ! Compute source at quadrature points
             call compute_source(u_quad(1:nvar,j),source_quad(1:nvar,j),chsi_quad(j),gamma,nvar)
          end do

          !================================
          ! Compute source volume integral
          !================================
          ! Loop over modes
          do i=1,n
             source_vol(1:nvar,i,icell)=0.0
             do j=1,nquad
                source_vol(1:nvar,i,icell)=source_vol(1:nvar,i,icell)+ &
                     & source_quad(1:nvar,j)* &
                     & legendre(chsi_quad(j),i-1)* &
                     & w_quad(j)
             end do
          end do
        end do
  end select

  !========================
  ! Compute final DG update
  !========================
  ! Loop over cells
  do icell=1,nx
     ! Loop over modes
     do i=1,n
        dudt(1:nvar,i,icell)=oneoverdx*(flux_vol(1:nvar,i,icell) &
             & -(flux_face(1:nvar,icell+1)*legendre(chsi_right,i-1) &
             & -flux_face(1:nvar,icell)*legendre(chsi_left,i-1)) &
             & ) + source_vol(1:nvar,i,icell)
     end do
  end do

  dudt(:,:,1) = 0.
  dudt(:,:,nx) = 0.

end subroutine compute_update
!==============================================
subroutine condinit(x,uu)
  use dg_commons
  real(kind=8)::x
  real(kind=8),dimension(1:nvar)::uu
  !==============================================
  ! This routine computes the initial conditions.
  !==============================================
  real(kind=8),dimension(1:nvar)::ww
  real(kind=8)::dpi=acos(-1d0)

  ! Compute primitive variables
  select case (ninit)
     case(1) ! sine wave (tend=1 or 10)
        ww(1)=1.0+0.5*sin(2.0*dpi*x)
        ww(2)=1.0
        ww(3)=1.0
     case(2) ! step function (tend=1 or 10)
        if(abs(x-0.5)<0.25)then
           ww(1)=2.
        else
           ww(1)=1.0
        endif
        ww(2)=1.0
        ww(3)=1.0
     case(3) ! gaussian + square pulses (tend=1 or 10)
        ww(1)=1.+exp(-(x-0.25)**2/2.0/0.05**2)
        if(abs(x-0.7)<0.1)then
           ww(1)=ww(1)+1.
        endif
        ww(2)=1.0
        ww(3)=1.0
     case(4) ! Sod test (tend=0.245)
        if(abs(x-0.25)<0.25)then
           ww(1)=1.0
           ww(2)=0.0
           ww(3)=1.0
        else
           ww(1)=0.125
           ww(2)=0.0
           ww(3)=0.1
        endif
     case(5) ! blast wave test  (tend=0.038)
        if(x<0.1)then
           ww(1)=1.0
           ww(2)=0.0
           ww(3)=1000.0
        else if(x<0.9)then
           ww(1)=1.0
           ww(2)=0.0
           ww(3)=0.01
        else
           ww(1)=1.0
           ww(2)=0.0
           ww(3)=100.
        endif
     case(6) ! shock entropy interaction (tend=2)
        if(x<10.0)then
           ww(1)=3.857143
           ww(2)=-0.920279
           ww(3)=10.333333
        else
           ww(1)=1.0+0.2*sin(5.0*(x-10.0))
           ww(2)=-3.549648
           ww(3)=1.0
        endif
     case(7) ! hydrostatic eq
           ww(1) = exp(-x)
           ww(2) = 0.
           ww(3) = exp(-x)
     case(8) ! hydro + pert
           ww(1) = exp(-x)
           ww(2) = 0.
           ww(3) = exp(-x) + pert*exp(-100*(x-boxlen/2.)**2)
     end select

     ! Convert primitive to conservative
     uu(1)=ww(1)
     uu(2)=ww(1)*ww(2)
     uu(3)=ww(3)/(gamma-1.0)+0.5*ww(1)*ww(2)**2

  return
end subroutine condinit
!==============================================
subroutine compute_max_speed(u,cmax)
  use dg_commons
  implicit none
  real(kind=8),dimension(1:nvar,1:n,1:nx)::u
  real(kind=8)::cmax
  !==============================================
  ! This routine computes the maximum wave speed.
  !==============================================
  integer::icell
  real(kind=8)::speed
  ! Compute max sound speed
  cmax=0.0
  do icell=1,nx
     call compute_speed(u(1:nvar,1,icell),speed,gamma,nvar)
     cmax=MAX(cmax,speed)
  enddo
end subroutine compute_max_speed

!==============================================
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
!==============================================
subroutine compute_speed(u,speed,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8),dimension(1:nvar)::u
  real(kind=8)::speed,gamma
  real(kind=8),dimension(1:nvar)::w
  real(kind=8)::cs
  ! Compute primitive variables
  call compute_primitive(u,w,gamma,nvar)
  ! Compute sound speed
  cs=sqrt(gamma*max(w(3),1d-10)/max(w(1),1d-10))
  speed=abs(w(2))+cs
end subroutine compute_speed
!==============================================
subroutine compute_flux(u,flux,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8),dimension(1:nvar)::u
  real(kind=8),dimension(1:nvar)::flux
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::w
  ! Compute primitive variables
  call compute_primitive(u,w,gamma,nvar)
  ! Compute flux
  flux(1)=w(2)*u(1)
  flux(2)=w(2)*u(2)+w(3)
  flux(3)=w(2)*u(3)+w(3)*w(2)
end subroutine compute_flux
!==============================================
subroutine compute_source(u,source,x_quad,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8),dimension(1:nvar)::u
  real(kind=8),dimension(1:nvar)::source
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::w
  real(kind=8)::x_quad
  ! Compute primitive variables
  call compute_primitive(u,w,gamma,nvar)
  ! Compute source (linear)
  source(1)=0
  source(2)=-w(1)*0.5 !-rho*phi_x
  source(3)=-w(1)*w(2)*0.5 ! -rho*u*phi_x
  ! hacky part: phi = gx, g = 1
end subroutine compute_source
!==============================================
subroutine compute_primitive(u,w,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::u
  real(kind=8),dimension(1:nvar)::w
  ! Compute primitive variables
  w(1)=u(1)
  w(2)=u(2)/w(1)
  w(3)=(gamma-1.0)*(u(3)-0.5*w(1)*w(2)**2)
end subroutine compute_primitive
!==============================================
subroutine cons_to_prim(du,dw,w,nvar,gamma)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::du,dw,w
  dw(1)=du(1)
  dw(2)=(du(2)-w(2)*du(1))/w(1)
  dw(3)=(gamma-1.0)*(0.5*w(2)**2*du(1)-w(2)*du(2)+du(3))
end subroutine cons_to_prim
!==============================================
subroutine prim_to_cons(dw,du,w,nvar,gamma)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::du,dw,w
  du(1)=dw(1)
  du(2)=w(2)*dw(1)+w(1)*dw(2)
  du(3)=0.5*w(2)**2*dw(1)+w(1)*w(2)*dw(2)+dw(3)/(gamma-1.0)
end subroutine prim_to_cons
!==============================================
subroutine cons_to_char(du,dw,w,nvar,gamma)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::du,dw,w,dp
  real(kind=8)::csq,cs
  csq=gamma*max(w(3),1d-10)/max(w(1),1d-10)
  cs=sqrt(csq)
  call cons_to_prim(du,dp,w,nvar,gamma)
  dw(1)=dp(1)-dp(3)/csq
  dw(2)=0.5*(dp(3)/csq+dp(2)*w(1)/cs)
  dw(3)=0.5*(dp(3)/csq-dp(2)*w(1)/cs)
end subroutine cons_to_char
!==============================================
subroutine char_to_cons(dw,du,w,nvar,gamma)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::du,dw,w,dp
  real(kind=8)::csq,cs
  csq=gamma*max(w(3),1d-10)/max(w(1),1d-10)
  cs=sqrt(csq)
  dp(1)=dw(1)+dw(2)+dw(3)
  dp(2)=(dw(2)-dw(3))*cs/w(1)
  dp(3)=(dw(2)+dw(3))*csq
  call prim_to_cons(dp,du,w,nvar,gamma)
end subroutine char_to_cons
!==============================================
subroutine cons_to_cons(du,dw,w,nvar,gamma)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::du,dw,w,dp
  real(kind=8)::csq,cs
  dw(1)=du(1)
  dw(2)=du(2)
  dw(3)=du(3)
end subroutine cons_to_cons
!==============================================
subroutine riemann_llf(uleft,uright,fgdnv,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::uleft,uright
  real(kind=8),dimension(1:nvar)::fgdnv
  real(kind=8)::cleft,cright,cmax
  real(kind=8),dimension(1:nvar)::fleft,fright
  ! Maximum wave speed
  call compute_speed(uleft,cleft,gamma,nvar)
  call compute_speed(uright,cright,gamma,nvar)
  cmax=max(cleft,cright)
  ! Compute flux at left and right points
  call compute_flux(uleft,fleft,gamma,nvar)
  call compute_flux(uright,fright,gamma,nvar)
  ! Compute Godunox flux
  fgdnv=0.5*(fright+fleft)-0.5*cmax*(uright-uleft)
end subroutine riemann_llf
!==============================================
subroutine riemann_hllc(ul,ur,fgdnv,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::ul,ur
  real(kind=8),dimension(1:nvar)::fgdnv

  real(kind=8)::cl,cr,dl,dr,sl,sr
  real(kind=8),dimension(1:nvar)::wl,wr,wstar,wstarl,wstarr
  real(kind=8),dimension(1:nvar)::ustarl,ustarr,wgdnv,ugdnv
  ! Compute primitive variables
  call compute_primitive(ul,wl,gamma,nvar)
  call compute_primitive(ur,wr,gamma,nvar)
  ! Compute sound speed
  cl=sqrt(gamma*max(wl(3),1d-10)/max(wl(1),1d-10))
  cr=sqrt(gamma*max(wr(3),1d-10)/max(wr(1),1d-10))
  ! Compute HLL wave speed
  SL=min(wl(2),wr(2))-max(cl,cr)
  SR=max(wl(2),wr(2))+max(cl,cr)
  ! Compute Lagrangian sound speed
  DL=wl(1)*(wl(2)-SL)
  DR=wr(1)*(SR-wr(2))
  ! Compute acoustic star state
  wstar(2)=(DR*wr(2)+DL*wl(2)+(wl(3)-wr(3)))/(DL+DR)
  wstar(3)=(DR*wl(3)+DL*wr(3)+DL*DR*(wl(2)-wr(2)))/(DL+DR)
  ! Left star region states
  wstarl(1)=wl(1)*(SL-wl(2))/(SL-wstar(2))
  ustarl(3)=((SL-wl(2))*ul(3)-wl(3)*wl(2)+wstar(3)*wstar(2))/(SL-wstar(2))
  ! Left star region states
  wstarr(1)=wr(1)*(SR-wr(2))/(SR-wstar(2))
  ustarr(3)=((SR-wr(2))*ur(3)-wr(3)*wr(2)+wstar(3)*wstar(2))/(SR-wstar(2))
  ! Sample the solution at x/t=0
  if(SL>0.0)then
     wgdnv(1)=wl(1)
     wgdnv(2)=wl(2)
     wgdnv(3)=wl(3)
     ugdnv(3)=ul(3)
  else if(wstar(2)>0.0)then
     wgdnv(1)=wstarl(1)
     wgdnv(2)=wstar(2)
     wgdnv(3)=wstar(3)
     ugdnv(3)=ustarl(3)
  else if(SR>0.0)then
     wgdnv(1)=wstarr(1)
     wgdnv(2)=wstar(2)
     wgdnv(3)=wstar(3)
     ugdnv(3)=ustarr(3)
  else
     wgdnv(1)=wr(1)
     wgdnv(2)=wr(2)
     wgdnv(3)=wr(3)
     ugdnv(3)=ur(3)
  end if
  fgdnv(1)=wgdnv(1)*wgdnv(2)
  fgdnv(2)=wgdnv(1)*wgdnv(2)*wgdnv(2)+wgdnv(3)
  fgdnv(3)=wgdnv(2)*(ugdnv(3)+wgdnv(3))
end subroutine riemann_hllc
!==============================================
