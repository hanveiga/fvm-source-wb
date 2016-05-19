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

  call evolve(u, u_eq, x)

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

  dx = boxlen/dble(nx)
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
      w(1,:) = 1. !exp(-x)
      w(2,:) = 0
      w(3,:) = 1. !exp(-x)
    case(2)
      w(1,:) = exp(-x)
      w(2,:) = 0
      w(3,:) = exp(-x)!+0.01*exp(-100*(x-0.5)**2)
    case(3)
      w(1,:) = (1 - ((gamma - 1)/gamma)*1*x)**(1/(gamma-1))
      w(2,:) = 0
      w(3,:) = (1 - ((gamma - 1)/gamma)*1*x)**(gamma/(gamma-1))
  end select

  call compute_conservative(w,u,nx)

end subroutine get_initial_conditions


subroutine output_file(x_values, function_values, filen)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::function_values
  real(kind=8),dimension(1:nx)::x_values
  real(kind=8),dimension(1:nvar)::w, w_eq
  character(len=*)::filen

  ! internal vars
  integer::icell, ivar
  real(kind=8)::dx, xcell, maximum_error
  real(kind=8),dimension(1:nvar)::primitives
  integer::one
  dx = boxlen/dble(nx)
  one = 1
  maximum_error = 0
!  write(filen,"(A5,I5.5)")
  open(10,status='REPLACE',file="simul/"//TRIM(filen)//".dat",form='formatted')
  do icell=1,nx
     xcell=(dble(icell)-0.5)*dx
     call get_equilibrium_solution([xcell],w_eq,one)
     call compute_primitive(function_values(1:nvar,icell),w,one)
     !write(10,'(7(1PE12.5,1X))')xcell,w(1)
     write(10,'(7(1PE12.5,1X))')xcell,(w(ivar),ivar=3,nvar)
     if (abs(w(3)-w_eq(3))>maximum_error) then
        maximum_error = abs( w(3) - w_eq(3) )
     end if

     end do
  close(10)
     write(*,*) 'maxerror', maximum_error
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
  !class(:)::x
  ! internal variables

  select case (nequilibrium)
    case(1)
      w(1,:) = 1. !exp(-x)
      w(2,:) = 0
      w(3,:) = 1. !exp(-x)
    case(2)
      w(1,:) = exp(-x)
      w(2,:) = 0
      w(3,:) = exp(-x)
    case(3)
      w(1,:) = (1 - ((gamma - 1)/gamma)*1*x)**(1/(gamma-1))
      w(2,:) = 0
      w(3,:) = (1 - ((gamma - 1)/gamma)*1*x)**(gamma/(gamma-1))
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
subroutine evolve(u, u_eq, x)
  use parameters
  implicit none

  real(kind=8),dimension(1:nvar,1:nx)::u,u_new,u_eq
  real(kind=8),dimension(1:nx)::x
  ! internal variables
  real(kind=8)::t,dt
  real(kind=8)::cmax, dx
  real(kind=8),dimension(1:nvar,1:nx)::dudt, w, w1
  integer::iter, n, i, snap_counter
  character(len=6)::filename

  dx = boxlen/dble(nx)

  t=0
  iter=0
  snap_counter = 0
  do while(t < tend)
  !do while( iter < 10)
     ! Compute time step
     call compute_max_speed(u,cmax)
     dt=0.8*dx/cmax/(2.0*dble(1)+1.0)

     if(solver=='FVM')then
       ! runge kutta 2nd order
       call compute_update_fvm(u,u_eq, dudt)
       w1=u+dt*dudt
       call compute_update_fvm(u,u_eq,dudt)
       u=0.5*u+0.5*w1+0.5*dt*dudt
     endif

     if(solver=='EQL')then
       ! runge kutta 2nd order
       call compute_update(u,u_eq, dudt)
       w1=u+dt*dudt
       call compute_update(u,u_eq,dudt)
       u=0.5*u+0.5*w1+0.5*dt*dudt
     endif

     if(solver=='WB1')then
       ! runge kutta 2nd order
       call compute_update_sr(u,u_eq, dudt)
       w1=u+dt*dudt
       call compute_update_sr(w1,u_eq,dudt)
       u=0.5*u+0.5*w1+0.5*dt*dudt
     endif

     t=t+dt
     iter=iter+1
     write(*,*)'time=',iter,t,dt

     if ((make_movie).and.(MODULO(iter,200)==0)) then
       snap_counter = snap_counter + 1
       write(filename,'(a, i3.3)') 'SIM', snap_counter
       call output_file(x,u,filename)
     end if

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
  real(kind=8)::oneoverdx, dx, dt, dx_plus, zero, xface, xcell
  real(kind=8),dimension(1:nvar)::flux_r, u_minus, u_plus, f_plus, f_minus, w_plus, w_minus, u_face, w_face
  integer::one, i, nfaces, iface, ileft, iright
  dx = boxlen/dble(nx)
  oneoverdx = 1/dx
  nfaces = nx + 1
  dx_plus = boxlen/dble(nfaces)
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
      w_minus(1:nvar) = (/1.,0.,1./)!(/ exp(-xcell), dble(0), exp(-xcell)/)
      !w_minus(1:nvar) = (/ (1 - ((gamma - 1)/gamma)*1*xcell)**(1/(gamma-1)), dble(0),&
      !& (1 - ((gamma - 1)/gamma)*1*xcell)**(gamma/(gamma-1))/)
      !call get_equilibrium_solution(xcell,w_minus, one)
      w_face = w_minus
      call compute_conservative(w_minus, u_minus, one)
      call compute_conservative(w_face, u_face, one)
      call compute_flux(u_face+delta_w(1:nvar,1), f_minus, one)
      call compute_llflux(u_face+delta_w(1:nvar,1), u_left(1:nvar,iright), f_minus, &
                      & f_left(1:nvar,iright), flux_riemann(1:nvar,iface))
    end if
    if (bc==3.and.iface==nx+1) then
            xcell = nx*dx !(nx+1-0.5)*dx
            w_plus(1:nvar) = (/1.,0.,1./) !(/ exp(-xcell), dble(0), exp(-xcell)/)
            !w_plus(1:nvar) = (/ (1 - ((gamma - 1)/gamma)*1*xcell)**(1/(gamma-1)), dble(0)&
            !&, (1 - ((gamma - 1)/gamma)*1*xcell)**(gamma/(gamma-1))/)
            !call get_equilibrium_solution(xcell,w_minus, one)
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

  dudt(1:nvar,1) =  dudt(1:nvar,2)
  dudt(1:nvar,nx) =  dudt(1:nvar,nx-1)

end subroutine compute_update
!-----
subroutine get_source(w,s,x,size)
  use parameters
  implicit none
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::u,s
  real(kind=8),dimension(1:nvar,1:size)::w
  real(kind=8),dimension(1:size)::x
  real(kind=8)::dx

  !internal
  real(kind=8),dimension(1:size)::x_minus, x_plus
  real(kind=8)::delta

  delta = 1/dble(size)

  x_minus = [x(1)-delta,x(1:size-1)] ! zero gradient
  x_plus = [x(2:size),x(size)+delta] ! zero gradient

  s(1,:) = 0
  s(2,:) = -w(1,:)*1*(x_plus-x_minus)/(2*delta)
  s(3,:) = -w(1,:)*w(2,:)*1*(x_plus-x_minus)/(2*delta)

  !s(1) = 0
  !  s(2) = -w(1)*1*(x_plus-x_minus)/(2*delta)
  !  s(3) = -w(1)*w(2)*1*(x_plus-x_minus)/(2*delta)

end subroutine get_source
!--------


subroutine get_source_rg(w,w_left,w_right,s,x,size)
  use parameters
  implicit none
  integer::size
  real(kind=8),dimension(1:nvar,1:size)::w,w_left,w_right,s
  real(kind=8),dimension(1:size)::x
  real(kind=8)::dx

  !internal
  real(kind=8),dimension(1:size)::x_minus, x_plus
  real(kind=8)::delta

  delta = 1./dble(size)

  x_minus = [x(1)-delta,x(1:size-1)] ! zero gradient
  x_plus = [x(2:size),x(size)+delta] ! zero gradient

  s(1,:) = 0
  s(2,:) = (w_right(3,:)-w_left(3,:))/delta
  s(3,:) = -w(1,:)*w(2,:)*1*(x_plus-x_minus)/(2*delta)

  !s(1) = 0
  !  s(2) = -w(1)*1*(x_plus-x_minus)/(2*delta)
  !  s(3) = -w(1)*w(2)*1*(x_plus-x_minus)/(2*delta)

end subroutine get_source_rg

!----

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
!--------
subroutine compute_update_fvm(u, w_eq, dudt)
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
  dx = boxlen/dble(nx)
  oneoverdx = 1/dx
  nfaces = nx + 1
  dx_plus = boxlen/dble(nfaces)
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

  u_left = u
  u_right = u
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
    ! treat bc
    flux_riemann(1:nvar,iface) = flux_r
    if (bc==3.and.iface==1) then
      xcell = 0 !-0.5*dx
      !xface = dx
      !w_face = (/ exp(-xface), dble(0), exp(-xface)/)
      w_minus(1:nvar) = (/1.,0.,1./)!(/ exp(-xcell), dble(0), exp(-xcell)/)
      w_face = w_minus
      call compute_conservative(w_minus, u_minus, one)
      call compute_conservative(w_face, u_face, one)
      call compute_flux(u_face+delta_w(1:nvar,1), f_minus, one)
      call compute_llflux(u_face+delta_w(1:nvar,1), u_left(1:nvar,iright), f_minus, &
                      & f_left(1:nvar,iright), flux_riemann(1:nvar,iface))
    end if
    if (bc==3.and.iface==nx+1) then
            xcell = nx*dx !(nx+0.5)*dx !(nx+1-0.5)*dx
            w_plus(1:nvar) = (/1.,0.,1./)!  (/ exp(-xcell), dble(0), exp(-xcell)/)
            call compute_conservative(w_plus, u_plus, one)
            call compute_flux(u_plus+delta_w(1:nvar,nx), f_plus, one)
            call compute_llflux(u_right(1:nvar,ileft)+delta_w(1:nvar,nx), u_plus, f_right(1:nvar,ileft), &
                              & f_plus, flux_riemann(1:nvar,iface))
    end if

  end do

  ! compute source
  call get_source(u,s,x, nx)

  dudt(1:nvar,1:nx) = -(flux_riemann(1:nvar,2:(nx+1))-flux_riemann(1:nvar,1:nx))*oneoverdx + s(1:nvar,1:nx)
  dudt(1:nvar,1) = dudt(1:nvar,2)
  dudt(1:nvar,nx) = dudt(1:nvar,nx-1)

end subroutine compute_update_fvm
!-----------------------------------------


subroutine compute_update_sr(u, w_eq, dudt)
  use parameters
  implicit none
  real(kind=8),dimension(1:nvar,1:nx)::u, w_eq, dudt
  real(kind=8),dimension(1:nx)::h, Kapp
  real(kind=8),dimension(1:nx)::h0_left, h0_right, phi_faces_left, phi_faces_right, phi_center


  ! internal vars
  real(kind=8),dimension(1:nvar,1:nx)::delta_w, f_left, f_right
  real(kind=8),dimension(1:(nx+1))::x_faces
  real(kind=8),dimension(1:(nx))::x
  real(kind=8),dimension(1:nvar,1:(nx+1))::u_eq_faces, flux_riemann, flux_eq, w_eq_faces
  real(kind=8),dimension(1:nvar, 1:nx)::u_left,u_right, u_eq, w_left, w_right, w
  real(kind=8),dimension(1:nvar,1:nx)::s,s_eq
  real(kind=8)::oneoverdx, dx, dt, dx_plus, xcell, zero, xface, w_temp
  real(kind=8),dimension(1:nvar)::flux_r, u_minus, u_plus, f_plus, f_minus, w_plus, w_minus, u_face, w_face
  integer::one, i, nfaces, iface, ileft, iright
  real(kind=8),dimension(1:nvar)::w_n, w_eq_t, w_plus_t
  real(kind=8)::phi, x_c, x_minus, x_plus
  dx = boxlen/dble(nx)
  oneoverdx = boxlen/dx

  nfaces = nx + 1
  dx_plus = 1/dble(nfaces)

  zero = 0

  ! create x_faces
  ! create x
  do i = 1,nx
    x(i) = (i-0.5)*dx
  end do
  do i = 1,nfaces
    x_faces(i) = (i-1)*dx
  end do
  write(*,*) 'dx', dx
  call compute_primitive(u,w,nx)

  do i=1,nx
    phi_center(i) = phi(x(i))
  end do 

  do i =1,nx
    phi_faces_left(i) = phi(x_faces(i))
  end do 

  write(*,*) 'phiface',maxval(phi_faces_left),minval(phi_faces_left)
  do i =1,nx
    phi_faces_right(i) = phi(x_faces(i+1))
  end do 

  write(*,*) 'phiface',maxval(phi_faces_right),minval(phi_faces_right)

  ! build h:
  h = max(w(3,:),1e-5)/max(1e-5,w(1,:))*(1+1./(gamma-1))
  write(*,*) 'h',maxval(h),minval(h)

  ! build h0_faces:
  h0_left = h + phi_center - phi_faces_left ! repreat h(0) and h(nx)
  h0_right = h + phi_center - phi_faces_right

  ! build k
  Kapp(:) = max(1e-5,w(3,:))/max(1e-5,w(1,:))**gamma
  write(*,*) 'Kapp'
  write(*,*) maxval(Kapp),minval(Kapp)
  ! compute u left, u right with equilibrium formula
  w_left = w
  w_right= w
  !do i = 1,nx
  !  x_c = i*dx
  !  x_minus = (i-1)*dx
  !  x_plus = (i+1)*dx
  ! w_left(3,i) = w(3,i) + w(1,i) * (phi(x_c) - phi(x_minus))/2
  !  w_right(3,i) = w(3,i) - w(1,i) * (phi(x_plus) - phi(x_c))/2
  !end do

  w_left(1,:) = ((1./Kapp)*(gamma-1)/gamma * h0_left)**(1/(gamma-1))
  w_left(3,:) = (1./Kapp)**(1/(gamma-1))*((gamma-1)/gamma * h0_left)**(gamma/(gamma-1))

  w_right(1,:) = ((1./Kapp)*(gamma-1)/gamma * h0_right)**(1/(gamma-1))
  w_right(3,:) = (1./Kapp)**(1/(gamma-1))*((gamma-1)/gamma * h0_right)**(gamma/(gamma-1))
  write (*,*) maxval(w_left), minval(w_left)
  write (*,*) maxval(w_right), minval(w_right)


  call compute_conservative(w_left,u_left,nx)
  call compute_conservative(w_right,u_right,nx)

  ! compute fluxes
  call compute_flux(u_left, f_left, nx)
  call compute_flux(u_right, f_right, nx)

  one = 1
  ! compute numerical flux
  do iface = 2,nfaces-1
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
    ! treat bc
    flux_riemann(1:nvar,iface) = flux_r
    ! if (bc==3.and.iface==1) then
       !xcell = 0 !-0.5*dx
    !   x_c = 0 !-0.5*dx
    !   x_minus = -0.5*dx
    !   w_minus(1:nvar) = (/ (1 - ((gamma - 1)/gamma)*1*x_c)**(1/(gamma-1)), dble(0),&
    !   & dble(0) /)
    !   w_temp = (1 - ((gamma - 1)/gamma)*1*x_c)**(gamma/(gamma-1))
    !   w_minus(3) = w_temp + w_minus(1)*(phi(x_c) - phi(x_minus))/2
    !   call compute_primitive(u(1:nvar,1),w_n,one)
    !   w_eq_t(1:nvar) = (/ (1 - ((gamma - 1)/gamma)*1*(0)*dx)**(1/(gamma-1)),&
    !          &dble(0),&
    !          & (1 - ((gamma - 1)/gamma)*1*(0)*dx)**(gamma/(gamma-1))/)
   !    w_minus = w_minus !+ (w_n - w_eq_t)

    !   call compute_conservative(w_minus, u_minus, one)
    !   call compute_flux(u_minus, f_minus, one)
    !   call compute_llflux(u_minus, u_left(1:nvar,iright), f_minus, &
    !                   & f_left(1:nvar,iright), flux_riemann(1:nvar,iface))
     !end if
     !if (bc==3.and.iface==nx+1) then
      !xcell = nx*dx !(nx+0.5)*dx !(nx+1-0.5)*dx
    !  x_plus = (nx+0.5)*dx
    !  x_c = (nx)*dx
    !  w_plus_t(1:nvar) = (/ (1 - ((gamma - 1)/gamma)*1*x_c)**(1/(gamma-1)), dble(0),&
    !  & dble(0)/)
    !  w_temp = (1 - ((gamma - 1)/gamma)*1*x_c)**(gamma/(gamma-1))
    !  w_plus_t(3) = w_temp - w_plus(1) * (phi(x_plus) - phi(x_c))/2
    !  call compute_primitive(u(1:nvar,nx),w_n,one)
    !  w_eq_t(1:nvar) = !(/ (1 - ((gamma - 1)/gamma)*1*(nx)*dx)**(1/(gamma-1)),&
    !                &dble(0),&
    !                & (1 - ((gamma - 1)/gamma)*1*(nx)*dx)**(gamma/(gamma-1))/)
    !  w_plus = w_plus_t !+ (w_n - w_eq_t)
    !  call compute_conservative(w_plus, u_plus, one)
    !  call compute_flux(u_plus, f_plus, one)
    !  call compute_llflux(u_right(1:nvar,ileft), u_plus, f_right(1:nvar,ileft), &
    !                  & f_plus, flux_riemann(1:nvar,iface))
    ! end if

      if (bc==3.and.iface==1) then
      w_minus(1:nvar) = (/1.,0.,1./)
      x_c = 0.5*dx !-0.5*dx
      x_minus = 0.!-0.5*dx
      w_minus(3) = w_minus(3) + w_minus(1)*(phi(x_c) - phi(x_minus))/2
      w_face = w_minus
      call compute_conservative(w_minus, u_minus, one)
      call compute_conservative(w_face, u_face, one)
      call compute_flux(u_face, f_minus, one)
      call compute_llflux(u_face, u_left(1:nvar,iright), f_minus, &
                      & f_left(1:nvar,iright), flux_riemann(1:nvar,iface))
    end if
    if (bc==3.and.iface==nx+1) then
            !xcell = nx*dx !(nx+1-0.5)*dx
            x_plus = (nx)*dx
            x_c = (nx-1)*dx
            w_plus(1:nvar) = (/1.,0.,1./) !(/ exp(-xcell), dble(0), exp(-xcell)/)
            !w_plus(1:nvar) = (/ (1 - ((gamma - 1)/gamma)*1*xcell)**(1/(gamma-1)), dble(0)&
            !&, (1 - ((gamma - 1)/gamma)*1*xcell)**(gamma/(gamma-1))/)
            !call get_equilibrium_solution(xcell,w_minus, one)
            w_plus(3) = w_plus(3) + w_minus(1)*(phi(x_plus) - phi(x_c))/2
            w_face = w_minus
      
            call compute_conservative(w_plus, u_plus, one)
            call compute_flux(u_plus, f_plus, one)
            call compute_llflux(u_right(1:nvar,ileft), u_plus, f_right(1:nvar,ileft), &
                              & f_plus, flux_riemann(1:nvar,iface))
    end if

  end do

  ! compute source
  call get_source_rg(w,w_left,w_right,s,x, nx)
  write(*,*) 's',minval(s), maxval(s)
  dudt(1:nvar,1:nx) = -(flux_riemann(1:nvar,2:(nx+1))-flux_riemann(1:nvar,1:nx))*oneoverdx + s(1:nvar,1:nx)
  !write(*,*) dudt
  dudt(:,1) =   dudt(:,2) !dudt(:,2)!dudt(:,2)
  dudt(:,nx) =  dudt(:,nx-1) ! dudt(:,nx-1)!dudt(:,nx-1)

  write(*,*) 'flux',minval(flux_riemann), maxval(flux_riemann)
  write(*,*) 'dudt',minval(dudt), maxval(dudt)
end subroutine compute_update_sr

function phi(x) result(phi_m)
  ! assuming linear phi
  real(kind=8)::g, x, phi_m
  g = 1
  phi_m = g*x

end function phi

function phi_sr(x) result(phi_m)
  ! assuming linear phi
  real(kind=8)::g, x, phi_m
  g = 1
  phi_m = g*x

end function phi_sr