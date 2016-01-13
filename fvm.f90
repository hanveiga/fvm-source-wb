! run
! fvm method code for euler eqs
! initial conditions
program fvm
  use fvm_commons
  implicit none
  !==============================================
  ! This program solves the 1D Euler equations
  ! using FVM
  !==============================================
  ! Main variables
  real(kind=8),dimension(1:nvar,1:nx)::u,dudt,w1,w2,w3,w4, ureal, uinit, prim
  integer::iter=0,icell,i,j,ivar, intnode
  real(kind=8)::xcell,dx,x_quad
  real(kind=8)::t,dt,cmax
  real(kind=8),dimension(1:nvar)::uu,ww
  character(len=20)::filename
  real(kind=8)::dpi=acos(-1d0)
  real(kind=8)::A

  ! Mesh size
  dx=boxlen/dble(nx)

  !===========================
  ! Compute initial conditions
  !===========================
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     ! Loop over modes coefficients
     call condinit(xcell,uu)
     u(1:nvar,icell) = uu(1:nvar)
  end do

  !=================================
  ! Write initial conditions to file
  !=================================
  write(filename,"(A5,I5.5)")"value",iter/10
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     call compute_primitive(u(1:nvar,icell),ww,gamma,nvar)
     write(10,'(7(1PE12.5,1X))')xcell,(ww(ivar),ivar=1,nvar)
  end do
  close(10)


  !============
  ! main loop
  !============

  ! get time step
  ! loop through time
  t=0
  iter=0
  A = 10**(-6)
  do while(t < tend)
     ! Compute time step
     call compute_max_speed(u,cmax)
     dt=0.8*dx/cmax/(2.0*dble(n)+1.0)

     ! insert perturbation
     !call compute_primitive(u(1:nvar,1), prim(1:nvar,1), gamma, nvar) ! add velocity perturbation
     !prim(2,1)=A*sin(4*dpi*t)
     !call compute_conservative(u(1:nvar,1), prim(1:nvar,1), gamma, nvar)

     ! runge kutta 2nd order
     call compute_update(u,dudt)
     w1=u+dt*dudt
     ! call limiter(w1)
     call compute_update(w1,dudt)
     u=0.5*u+0.5*w1+0.5*dt*dudt
     ! call limiter(u)
     t=t+dt
     iter=iter+1
     write(*,*)'time=',iter,t,dt
  end do

  write(filename,"(A5,I5.5)")"value",99999
  open(10,file=TRIM(filename)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     call compute_primitive(u(1:nvar,icell),ww,gamma,nvar)
     write(10,'(7(1PE12.5,1X))')xcell,(ww(ivar),ivar=1,nvar)
  end do
  close(10)


  write(*,*)'========================================'
  write(*,*)'time=',t,dt
  do icell=1,nx
       xcell=(dble(icell)-0.5)*dx
       call compute_primitive(u(1:nvar,icell),ww,gamma,nvar)
       write(*,'(7(1PE12.5,1X))')xcell,(ww(ivar),ivar=1,nvar)
  end do

end program fvm

!------------------------------------------------
subroutine condinit(x,uu)
  use fvm_commons
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
    case(7) !
          ww(1) = (1-(gamma-1)/(dble(gamma))*x)
          ww(2) = 0
          ww(3) = ww(1)**gamma
     end select
     ! Convert primitive to conservative
     uu(1)=ww(1)
     uu(2)=ww(1)*ww(2)
     uu(3)=ww(3)/(gamma-1.0)+0.5*ww(1)*ww(2)**2
  return
end subroutine condinit
!---------------------------------
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
!-------------------------------
subroutine compute_update(u,dudt)
  use fvm_commons
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::u,dudt, source_term
  real(kind=8),dimension(1:nvar,1:nx)::fluxleft,fluxright
  real(kind=8)::dx,x_quad
  real(kind=8)::cmax,oneoverdx,c_left,c_right, x_minus, x_plus
  integer::icell,i,j,iface,ileft,iright,ivar, u_plus, u_minus, u_c

  dx=boxlen/dble(nx)
  oneoverdx=1.0/dx
  ! dudt = flux approximation using LLF flux
  do iface = 1,nx ! there's gonna be overlapping values probably
    u_minus = iface-1
    u_c = iface
    u_plus = iface+1
    if(bc==1) then     ! periodic BC
      if (iface==1) then
        u_minus = nx
        u_c = iface
        u_plus = iface+1
        ! lala
      endif
      if (iface==nx) then
        u_minus = iface-1
        u_c = iface
        u_plus = 1
      endif
    endif
    if(bc==2) then
      ! zero gradient BC
        if (iface==1) then
          u_minus = 1
          call compute_llflux(u(1:nvar,u_minus),u(1:nvar,u_c),fluxleft(1:nvar,iface),gamma,nvar)
          write (*,*) fluxleft(1:nvar,iface)
        endif
        if (iface==nx) then
          u_plus = nx
        endif
    endif
    call compute_llflux(u(1:nvar,u_minus),u(1:nvar,u_c),fluxleft(1:nvar,iface),gamma,nvar)
    call compute_llflux(u(1:nvar,u_c),u(1:nvar,u_plus),fluxright(1:nvar,iface),gamma,nvar)
  end do

  select case (source)
     case(1)
       source_term(:,:) = 0.0
     case(2)
       do icell=1,nx
          x_minus = (dble(icell-1)-0.5)*dx
          x_plus = (dble(icell+1)-0.5)*dx
          if(icell==1) x_minus = (dble(1)-0.5)*dx
          if(icell==nx) x_plus = (dble(nx)-0.5)*dx
          call compute_source(u(1:nvar,icell),source_term(1:nvar,icell),x_minus,x_plus,dx,gamma,nvar)
       end do
  end select


  do icell = 1,nx
    dudt(1:nvar, icell) = -oneoverdx*(fluxright(1:nvar, icell) - fluxleft(1:nvar, icell)) + source_term(1:nvar,icell)
  end do


end subroutine compute_update
!-------------------------------
subroutine compute_source(u,s,x_minus,x_plus,delta,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8)::gamma, x_minus, x_plus,delta
  real(kind=8),dimension(1:nvar)::u,s
  real(kind=8),dimension(1:nvar)::w
  call compute_primitive(u,w,gamma,nvar)
  s(1) = 0
  s(2) = -w(1)*1*(x_plus-x_minus)/(2*delta)
  s(3) = -w(1)*w(2)*1*(x_plus-x_minus)/(2*delta)

end subroutine compute_source
!-------------------------------
subroutine limiter()

end subroutine limiter

!-------------------------------
subroutine compute_llflux(uleft,uright, fgdnv, gamma, nvar)
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
end subroutine compute_llflux
!-----
subroutine compute_flux(u,flux, gamma, nvar)
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
subroutine compute_max_speed(u,cmax)
  use fvm_commons
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::u
  real(kind=8)::cmax
  !==============================================
  ! This routine computes the maximum wave speed.
  !==============================================
  integer::icell
  real(kind=8)::speed
  ! Compute max sound speed
  cmax=0.0
  do icell=1,nx
     call compute_speed(u(1:nvar,icell),speed,gamma,nvar)
     cmax=MAX(cmax,speed)
  enddo
end subroutine compute_max_speed
!==============================================
subroutine compute_conservative(u,prim,gamma,nvar)
  implicit none
  integer::nvar
  real(kind=8)::gamma
  real(kind=8),dimension(1:nvar)::u
  real(kind=8),dimension(1:nvar)::prim
  ! Convert primitive to conservative
  u(1)=prim(1)
  u(2)=prim(1)*prim(2)
  u(3)=prim(3)/(gamma-1.0)+0.5*prim(1)*prim(2)**2
end subroutine
