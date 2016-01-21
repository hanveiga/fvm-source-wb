!----
! global variables: parameters

program main
  use parameters
  real(kind=8),dimension(1:nx)::x
  real(kind=8),dimension(1:nvar,1:nx)::u, w, u_eq

  call get_x(x,nx)

  call get_initial_conditions(x,u,nx)

  call output_file(x,u,'initial')

  call get_equilibrium_solution(x,u_eq,nx)

  call evolve(u, u_eq)

  call output_file(x,u, 'endofsim')

end program main

!---------
subroutine get_x(x,size)
  use parameters
  integer::size
  real(kind=8),dimension(1:size)::x

  !internal vars
  real(kind=8)::dx, xcell
  integer::i

  dx = 1/dble(nx)
  do i=1,nx
    x(i) = (i-0.5)*dx
  end do

end subroutine get_x
!-------
subroutine get_initial_conditions(x,u,size)
  use parameters
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::u
  real(kind=8),dimension(1:size)::x

  ! internal variables
  real(kind=8),dimension(1:nvar,1:size)::w

  select case (ninit)
    case(1)
      w(1,:) = exp(-x)
      w(2,:) = 0
      w(3,:) = exp(-x)
    case(2)
      w(1,:) = exp(-x)
      w(2,:) = 0
      w(3,:) = exp(-x)+0.0001*exp(-100*(x-0.5)**2)
  end select

  call compute_conservative(w,u,nx)

end subroutine get_initial_conditions

subroutine output_file(x_values, function_values, filen)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::function_values
  real(kind=8),dimension(1:nx)::x_values
  real(kind=8),dimension(1:nvar)::w, w_eq
  character(len=2)::filen

  ! internal vars
  integer::icell, ivar
  real(kind=8)::dx, xcell
  real(kind=8),dimension(1:nvar)::primitives
  integer::one
  dx = 1/dble(nx)
  one = 1
!  write(filen,"(A5,I5.5)")
  open(10,file=TRIM(filen)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     call get_equilibrium_solution([xcell],w_eq,one)
     call compute_primitive(function_values(1:nvar,icell),w,one)
     !write(10,'(7(1PE12.5,1X))')xcell,w(1)
     write(10,'(7(1PE12.5,1X))')xcell,(w(ivar),ivar=3,nvar)
  end do
  close(10)

end subroutine output_file
!--------
subroutine compute_primitive(u,w,size)
  use parameters
  implicit none
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::u
  real(kind=8),dimension(1:nvar,1:size)::w
  ! Compute primitive variables
  w(1,:)=u(1,:)
  w(2,:)=u(2,:)/w(1,:)
  w(3,:)=(gamma-1.0)*(u(3,:)-0.5*w(1,:)*w(2,:)**2)
end subroutine compute_primitive
!-------------
subroutine compute_conservative(ww,u,size)
  use parameters
  implicit none
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::u
  real(kind=8),dimension(1:nvar,1:size)::ww
  ! Compute primitive variables
  u(1,:)=ww(1,:)
  u(2,:)=ww(1,:)*ww(2,:)
  u(3,:)=ww(3,:)/(gamma-1.0)+0.5*ww(1,:)*ww(2,:)**2
end subroutine compute_conservative
!-------------
subroutine get_equilibrium_solution(x,w, size)
  ! get equilibrium solution for primitive variables
  ! use parameters
  ! communicate which initial solution, equilibrium type
  ! boundary conditions and domain properties
  use parameters
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::w
  real(kind=8),dimension(1:size)::x
  ! internal variables

  select case (nequilibrium)
    case(1)
      w(1,:) = exp(-x)
      w(2,:) = 0
      w(3,:) = exp(-x)
    case(2)
      w(1,:) = exp(-x)
      w(2,:) = 0
      w(3,:) = exp(-x)
  end select

end subroutine get_equilibrium_solution
!---------

subroutine compute_max_speed(u,cmax)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::u
  real(kind=8)::cmax
  integer::icell
  real(kind=8)::speed
  ! Compute max sound speed
  cmax=0.0
  do icell=1,nx
     call compute_speed(u(1:nvar,icell),speed)
     cmax=MAX(cmax,speed)
  enddo
end subroutine compute_max_speed
!--------------------
subroutine compute_speed(u,speed)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar)::u
  real(kind=8)::speed
  real(kind=8),dimension(1:nvar)::w
  real(kind=8)::cs
  ! Compute primitive variables
  call compute_primitive(u,w,1)
  ! Compute sound speed
  cs=sqrt(gamma*max(w(3),1d-10)/max(w(1),1d-10))
  speed=abs(w(2))+cs
end subroutine compute_speed
!==============================================
subroutine compute_flux(u,flux,size)
  use parameters
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::u
  real(kind=8),dimension(1:nvar,1:size)::flux
  real(kind=8),dimension(1:nvar,1:size)::w
  ! Compute primitive variables
  call compute_primitive(u,w,size)
  ! Compute flux
  flux(1,:)=w(2,:)*u(1,:)
  flux(2,:)=w(2,:)*u(2,:)+w(3,:)
  flux(3,:)=w(2,:)*u(3,:)+w(3,:)*w(2,:)
end subroutine compute_flux
!==============================
subroutine evolve(u, u_eq)
  use parameters
  implicit none

  real(kind=8),dimension(1:nvar,1:nx)::u,u_new,u_eq

  ! internal variables
  real(kind=8)::t,dt
  real(kind=8)::cmax, dx
  real(kind=8),dimension(1:nvar,1:nx)::dudt, w, w1
  integer::iter, n, i

  dx = 1/dble(nx)

  t=0
  iter=0
  do while(t < tend)
  ! do i = 1,2
     ! Compute time step
     call compute_max_speed(u,cmax)
     dt=0.8*dx/cmax/(2.0*dble(1)+1.0)
     ! runge kutta 2nd order
     call compute_update(u,u_eq, dudt)
     w1=u+dt*dudt
     call compute_update(w1,u_eq,dudt)
     u=0.5*u+0.5*w1+0.5*dt*dudt
     t=t+dt
     iter=iter+1
     write(*,*)'time=',iter,t,dt
  end do

  u_new = u
end subroutine evolve
!--------
subroutine compute_update(u, w_eq, dudt)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::u, w_eq, dudt

  ! internal vars
  real(kind=8),dimension(1:nvar,1:nx)::delta_w, f_left, f_right
  real(kind=8),dimension(1:(nx+1))::x_faces
  real(kind=8),dimension(1:(nx))::x
  real(kind=8),dimension(1:nvar,1:(nx+1))::u_eq_faces, flux_riemann, flux_eq, w_eq_faces
  real(kind=8),dimension(1:nvar, 1:nx)::u_left,u_right, u_eq
  real(kind=8),dimension(1:nvar,1:nx)::s,s_eq
  real(kind=8)::oneoverdx, dx, dt, dx_plus, xcell, zero, xface
  real(kind=8),dimension(1:nvar)::flux_r, u_minus, u_plus, f_plus, f_minus, w_plus, w_minus, u_face, w_face
  integer::one, i, nfaces, iface, ileft, iright
  dx = 1/dble(nx)
  oneoverdx = 1/dx
  nfaces = nx + 1
  dx_plus = 1/dble(nfaces)
  zero = 0
  ! compute perturbation
  call compute_conservative(w_eq,u_eq, nx)

  delta_w = u - u_eq

  ! propagate perturbation to faces
  ! create x_faces
  ! create x
  do i = 1,nx
    x(i) = (i-0.5)*dx
  end do
  do i = 1,nfaces
    x_faces(i) = (i-1)*dx
  end do


  call get_equilibrium_solution(x_faces,w_eq_faces,nfaces)
  call compute_conservative(w_eq_faces, u_eq_faces, nfaces)

  u_left(1:nvar,1:nx) = delta_w(1:nvar,1:nx) + u_eq_faces(1:nvar,1:nx)
  u_right(1:nvar,1:nx) = delta_w(1:nvar,1:nx) + u_eq_faces(1:nvar,2:(nx+1))

  ! compute fluxes
  call compute_flux(u_left, f_left, nx)
  call compute_flux(u_right, f_right, nx)
  call compute_flux(u_eq_faces,flux_eq,nfaces)

  one = 1
  ! compute numerical flux
  do iface = 1,nfaces
    ileft = iface - 1
    iright = iface
    if (bc==1) then ! periodic
      if(iface==1) ileft = nx
      if(iface==nx+1) iright = 1
    end if
    if (bc==2) then ! zero gradient
      if(iface==1) ileft = 1
      if(iface==nx+1) iright = nx
    end if
    if (bc==3) then ! ! include reflexive BC
      if(iface==1) ileft = 1
      if(iface==nx+1) iright = nx
    end if

    call compute_llflux(u_right(1:nvar,ileft),u_left(1:nvar,iright), f_right(1:nvar,ileft), &
                    & f_left(1:nvar,iright), flux_r)

    !call compute_llflux(u_right(1:nvar,ileft),u_left(1:nvar,iright), f_right(1:nvar,ileft), &
    !                    & f_left(1:nvar,iright), flux_r)

    !call compute_llflux(u_right(1:nvar,ileft),u_left(1:nvar,iright), f_right(1:nvar,ileft), &
    !                                    & f_left(1:nvar,iright), flux_l)
    ! treat bc
    flux_riemann(1:nvar,iface) = flux_r
    if (bc==3.and.iface==1) then
      xcell = 0 !dx !(0-0.5)*dx
      !xface = dx
      !w_face = (/ exp(-xface), dble(0), exp(-xface)/)
      w_minus(1:nvar) = (/ exp(-xcell), dble(0), exp(-xcell)/)
      w_face = w_minus
      call compute_conservative(w_minus, u_minus, one)
      call compute_conservative(w_face, u_face, one)
      call compute_flux(u_face+delta_w(1:nvar,1), f_minus, one)
      call compute_llflux(u_face+delta_w(1:nvar,1), u_left(1:nvar,iright), f_minus, &
                      & f_left(1:nvar,iright), flux_riemann(1:nvar,iface))
    end if
    if (bc==3.and.iface==nx+1) then
            xcell = nx*dx !(nx+1-0.5)*dx
            w_plus(1:nvar) = (/ exp(-xcell), dble(0), exp(-xcell)/)
            call compute_conservative(w_plus, u_plus, one)
            call compute_flux(u_plus+delta_w(1:nvar,nx), f_plus, one)
            call compute_llflux(u_right(1:nvar,ileft)+delta_w(1:nvar,nx), u_plus, f_right(1:nvar,ileft), &
                              & f_plus, flux_riemann(1:nvar,iface))
    end if

  end do

  ! compute source
  call get_source(u,s,x, nx)
  call get_source(u_eq,s_eq, x, nx)

  dudt(1:nvar,1:nx) = -(flux_riemann(1:nvar,2:(nx+1))-flux_riemann(1:nvar,1:nx))*oneoverdx + s(1:nvar,1:nx) + &
        & (flux_eq(1:nvar,2:(nx+1))-flux_eq(1:nvar,1:nx))*oneoverdx - s_eq(1:nvar,1:nx)


end subroutine compute_update
!-----
subroutine get_source(w,s,x,size)
  use parameters
  implicit none
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::u,s
  real(kind=8),dimension(1:nvar,1:size)::w
  real(kind=8),dimension(1:size)::x

  !internal
  real(kind=8),dimension(1:size)::x_minus, x_plus
  real(kind=8)::delta

  delta = 1/dble(size)

  x_minus = [x(1),x(1:size-1)] ! zero gradient
  x_plus = [x(2:size),x(size)] ! zero gradient

  s(1,:) = 0
  s(2,:) = -w(1,:)*1*(x_plus-x_minus)/(2*delta)
  s(3,:) = -w(1,:)*w(2,:)*1*(x_plus-x_minus)/(2*delta)

  !s(1) = 0
  !  s(2) = -w(1)*1*(x_plus-x_minus)/(2*delta)
  !  s(3) = -w(1)*w(2)*1*(x_plus-x_minus)/(2*delta)

end subroutine get_source
!--------
subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar)::uleft,uright, f_left, f_right
  real(kind=8),dimension(1:nvar)::fgdnv
  real(kind=8)::cleft,cright,cmax
  real(kind=8),dimension(1:nvar)::fleft,fright
  ! Maximum wave speed
  call compute_speed(uleft,cleft)
  call compute_speed(uright,cright)
  cmax=max(cleft,cright)
  ! Compute Godunox flux
  fgdnv=0.5*(f_right+f_left)-0.5*cmax*(uright-uleft)
end subroutine compute_llflux
!-----
