!----
! global variables: parameters

program main
  USE ISO_C_BINDING
  use parameters_dg_2d
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::x
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::y
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u, w, u_eq
  
  !device pointers
  type (c_ptr) :: u_d, du_d,dudt_d, w_d, u_eq_d 
  type (c_ptr) :: x_d, y_d
   
  call get_coords(x,y,nx,ny, mx, my)

  call get_initial_conditions(x,y,u,nx,ny, mx, my)
  call output_file(x, y ,u,1,'initcond')

  call get_equilibrium_solution(x, y, u_eq, nx, ny, mx, my)

  !device allocation and memory copies 
  call devices()
  call setdevice(0)
  call gpu_allocation(nvar,nx,ny,mx,my,boxlen_x,boxlen_y,cfl,eta,gamma)
  call gpu_set_pointers(u_d,du_d,dudt_d,w_d,u_eq_d,x_d,y_d &
            & ,x_quad,y_quad,w_x_quad,w_y_quad,sqrt_mod)
  call h2d(u,u_d,nvar*nx*ny*mx*my)
  call h2d(w,w_d,nvar*nx*ny*mx*my)
  call h2d(u_eq,u_eq_d,nvar*nx*ny*mx*my)
  call h2d(x,x_d,nx*ny*mx*my)
  call h2d(y,y_d,nx*ny*mx*my)
  
  call evolve(u, x, y, u_eq, u_d, du_d, dudt_d)

  call output_file(x, y ,u,1,'fintwo')

end program main

!---------
subroutine get_coords(x,y,size_x,size_y, order_x, order_y)
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y, order_x, order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x)::x,y

  !internal vars
  real(kind=8)::dx, dy
  integer::i,j, nj, ni

  ! get quadrature
  call gl_quadrature(x_quad,w_x_quad,order_x)
  call gl_quadrature(y_quad,w_y_quad,order_y)
  do i=1,mods
    sqrt_mod(i)=sqrt((2.0*dble(i-1)+1.0))/sqrt(2.)
  end do  
  dx = boxlen_x/dble(size_x)
  dy = boxlen_y/dble(size_y)
  do i=1,size_x
    do j=1,size_y
      do ni =1,order_x
        do nj = 1,order_y
          x(nj,ni,j,i) = (i-0.5)*dx + dx/2.0*x_quad(ni)
          y(nj,ni,j,i) = (j-0.5)*dy + dy/2.0*y_quad(nj)
        end do
      end do
    end do
  end do

end subroutine get_coords
!-------
subroutine get_initial_conditions(x,y,u,size_x,size_y, order_x, order_y)
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y, i ,j, order_x, order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::u,w
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x)::x,y
  integer::inti,intj, jcell, icell
  ! internal variables
  real(kind=8)::rho_0, p_0, g
  real(kind=8)::dpi=acos(-1d0)

  real(kind=8)::x_dash, y_dash, r, rho_d, delta_r, x_center, y_center

  select case (ninit)
  case(1) !Linear advection of pulse
    w(:,:,:,:,1) = exp(-((x(:,:,:,:)-boxlen_x/2.)**2+(y(:,:,:,:)-boxlen_y/2.)**2)*20)
    w(:,:,:,:,2) = 10.0
    w(:,:,:,:,3) = 10.0
    w(:,:,:,:,4) = minval(w(:,:,:,:,1))
  case(2) !hydrostatic equilibrium with a pulse on pressure
    rho_0 = 1.21
    p_0 = 1.
    g = 1.
    w(:,:,:,:,1) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
    w(:,:,:,:,2) = 0
    w(:,:,:,:,3) = 0
    w(:,:,:,:,4) = p_0*exp(-(rho_0*g/p_0)*(x+y)) + eta*&
    &exp(-100*(rho_0*g/p_0)*((x-0.3)**2+(y-0.3)**2))
  case(3) !2d riemann problem case 4
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y
            if (x(intj,inti,j,i)>=0.5 .and. y(intj,inti,j,i)>=0.5) then
              w(intj,inti,j,i,1) = 1.
              w(intj,inti,j,i,2)  = 0.
              w(intj,inti,j,i,3)  = 0.
              w(intj,inti,j,i,4)  = 1.
            else if (x(intj,inti,j,i)<0.5 .and. y(intj,inti,j,i)>=0.5) then
              w(intj,inti,j,i,1) = 0.5197
              w(intj,inti,j,i,2) = -0.7259
              w(intj,inti,j,i,3) = 0.
              w(intj,inti,j,i,4) = 0.4
            else if (x(intj,inti,j,i)<0.5 .and. y(intj,inti,j,i)<0.5) then
              w(intj,inti,j,i,1) = 0.1072
              w(intj,inti,j,i,2) = -0.7259
              w(intj,inti,j,i,3) = -1.4045
              w(intj,inti,j,i,4) = 0.0439
            else if (x(intj,inti,j,i)>=0.5 .and. y(intj,inti,j,i)<0.5) then
              !else
              w(intj,inti,j,i,1) = 0.2579
              w(intj,inti,j,i,2) = 0.0
              w(intj,inti,j,i,3) = -1.4045
              w(intj,inti,j,i,4) = 0.15
            end if
          end do
        end do
      end do
    end do
  case(4) ! another riemann problem
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y
            if (x(intj,inti,j,i)>=0.5 .and. y(intj,inti,j,i)>=0.5) then
              w(intj,inti,j,i,1) = 1.5
              w(intj,inti,j,i,2)  = 0.
              w(intj,inti,j,i,3)  = 0.
              w(intj,inti,j,i,4)  = 1.5
            else if (x(intj,inti,j,i)<0.5 .and. y(intj,inti,j,i)>=0.5) then
              w(intj,inti,j,i,1) = 0.5323
              w(intj,inti,j,i,2) = 1.206
              w(intj,inti,j,i,3) = 0.
              w(intj,inti,j,i,4) = 0.3
            else if (x(intj,inti,j,i)<0.5 .and. y(intj,inti,j,i)<0.5) then
              w(intj,inti,j,i,1) = 0.138
              w(intj,inti,j,i,2) = 1.206
              w(intj,inti,j,i,3) = 1.206
              w(intj,inti,j,i,4) = 0.029
            else if (x(intj,inti,j,i)>=0.5 .and. y(intj,inti,j,i)<0.5) then
              !else
              w(intj,inti,j,i,1) = 0.5323
              w(intj,inti,j,i,2) = 0.0
              w(intj,inti,j,i,3) = 1.206
              w(intj,inti,j,i,4) = 0.3
            end if
          end do
        end do
      end do
    end do
  case(5)
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y
            if (x(intj,inti,j,i) + y(intj,inti,j,i) >=0.5) then
              w(intj,inti,j,i,1) = 1.
              w(intj,inti,j,i,2)  = 0.
              w(intj,inti,j,i,3)  = 0.
              w(intj,inti,j,i,4)  = 1.
            else if (x(intj,inti,j,i) + y(intj,inti,j,i) < 0.5) then
              w(intj,inti,j,i,1) = 0.125
              w(intj,inti,j,i,2) = 0.
              w(intj,inti,j,i,3) = 0.
              w(intj,inti,j,i,4) = 0.4
            end if
          end do
        end do
      end do
    end do
  case(6) ! vortex http://flash.uchicago.edu/~jbgallag/2012/flash4_ug/node34.html
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y

            w(intj,inti,j,i,1) = 1.*(1. &
            &-(gamma-1.)*5/(8*gamma*dpi**2)*exp(1-&
            &((x(intj,inti,j,i)-5)**2+(y(intj,inti,j,i)-5)**2)))**(1/(gamma-1))

            w(intj,inti,j,i,2) = 2+ &
            &5./(2*dpi)*exp(-1-((x(intj,inti,j,i)-5)**2+&
            &(y(intj,inti,j,i)-5)**2)/2.)*(-y(intj,inti,j,i)+5.)

            w(intj,inti,j,i,3)  = 2+ &
            &5./(2*dpi)*exp(-1-((x(intj,inti,j,i)-5)**2+&
            &(y(intj,inti,j,i)-5)**2)/2.)*&
            &(x(intj,inti,j,i)-5.)

            w(intj,inti,j,i,4)  = w(intj,inti,j,i,1)**gamma
          end do
        end do
      end do
    end do
  case(7) ! smooth rotating disk
    ! initialize variables for smooth rotating disk
    ! boxlen has to be 6x6
    p_0 = 10e-5
    rho_0 = 10e-5
    rho_d = 1
    delta_r = 0.1
    x_center = 3.
    y_center = 3.

    w(:,:,:,:,4)  = p_0
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y
            x_dash = x(intj,inti,j,i) - x_center
            y_dash = y(intj,inti,j,i) - y_center
            r = sqrt(x_dash**2 + y_dash**2)

            if (r<0.5-delta_r/2.) then
              w(intj,inti,j,i,1) = rho_0
            else if ((r<0.5+delta_r/2.).and.(r > 0.5-delta_r/2.)) then
              w(intj,inti,j,i,1) = (rho_d-rho_0)/delta_r * (r - (0.5-delta_r/2.)) + rho_0
            else if ((r >= 0.5+delta_r/2.).and.(r <= 2-delta_r/2.))  then
              w(intj,inti,j,i,1) = rho_d
            else if ((r > 2-delta_r/2.).and.(r < 2+delta_r/2.))  then
              w(intj,inti,j,i,1) = (rho_0-rho_d)/delta_r * (r - (2-delta_r/2.)) + rho_d
            else if (r >= 2+delta_r/2.) then
              w(intj,inti,j,i,1) = rho_0
            end if

            if (r<=0.5-delta_r) then
              w(intj,inti,j,i,2) = 0.
              w(intj,inti,j,i,3) = 0.
            else if ((r<=0.5-delta_r).and.(r > 0.5-2*delta_r)) then
              w(intj,inti,j,i,2) = -y_dash/r**(3./2.)/(delta_r) * (r - (0.5-2*delta_r))
              w(intj,inti,j,i,3) = x_dash/r**(3./2.)/(delta_r) * (r - (0.5-2*delta_r))
            else if ((r > 0.5-delta_r).and.(r <= 2+delta_r))  then
              w(intj,inti,j,i,2) = -y_dash/r**(3./2.)
              w(intj,inti,j,i,3) = x_dash/r**(3./2.)
            else if ((r > 2+delta_r).and.(r <= 2+2*delta_r))  then
              w(intj,inti,j,i,2) = y_dash/r**(3./2.)/delta_r * (r - (2+delta_r)) - y_dash/r**(3./2.)
              w(intj,inti,j,i,3) = -x_dash/r**(3./2.)/delta_r * (r - (2+delta_r)) + x_dash/r**(3./2.)
            else if (r > 2+2*delta_r) then
              w(intj,inti,j,i,2) = 0.
              w(intj,inti,j,i,3) = 0.
            end if

          end do
        end do
      end do
    end do

  case(8) ! square advection
    p_0 = 10e-5
    rho_0 = 10e-5
    rho_d = 1
    delta_r = 0.1
    x_center = 0.5
    y_center = 0.5
    w(:,:,:,:,4) = p_0
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y
            x_dash = x(intj,inti,j,i) - x_center
            y_dash = y(intj,inti,j,i) - y_center

            if ((abs(x_dash)<=0.25).and.(abs(y_dash)<=0.25)) then
              w(intj,inti,j,i,1) = rho_d
            else
              w(intj,inti,j,i,1) = rho_0
            end if

            w(intj,inti,j,i,2) = 1.0 !1.0
            w(intj,inti,j,i,3) = 1.0 !1.0

            w(intj,inti,j,i,4)  = p_0
          end do
        end do
      end do
    end do

  case(9) ! 1d discontinuous pulse advection
    p_0 = 10e-5
    rho_0 = 10e-5
    rho_d = 1
    delta_r = 0.1
    x_center = 0.5
    y_center = 0.5
    w(:,:,:,:,4) = p_0
    do i = 1,size_x
      do j = 1,size_y
        do inti = 1,order_x
          do intj = 1,order_y
            x_dash = x(intj,inti,j,i) - x_center
            y_dash = y(intj,inti,j,i) - y_center

            if ((abs(x_dash)<=0.25)) then
              w(intj,inti,j,i,1) = rho_d
            else
              w(intj,inti,j,i,1) = rho_0
            end if

            w(intj,inti,j,i,2) = 0.
            w(intj,inti,j,i,3) = 1.0

          end do
        end do
      end do
    end do
  case(10) !hydrostatic
    rho_0 = 1.21
    p_0 = 1
    g = 1
    w(:,:,:,:,1) = rho_0*exp(-(rho_0*g/p_0)*(x(:,:,:,:)+y(:,:,:,:)))
    w(:,:,:,:,2) = 0.
    w(:,:,:,:,3) = 0.
    w(:,:,:,:,4) = p_0*exp(-(rho_0*g/p_0)*(x+y))
  end select
  call compute_conservative(w,u,nx,ny,mx,my)
  
end subroutine get_initial_conditions

subroutine output_file(x_values, y_values, function_values, var, filen)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::function_values
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::x_values, y_values
  real(kind=8),dimension(1:nvar)::w, w_eq
  character(len=*)::filen

  ! internal vars
  integer::icell, ivar, jcell
  integer::var
  real(kind=8)::dx, xcell, ycell
  real(kind=8),dimension(1:nvar)::primitives
  integer::one
  one = 1
  open(10,status='REPLACE',file="longrotation2_o1/"//TRIM(filen)//".dat",form='formatted')
  do icell=1,nx
    do jcell=1,ny
      xcell= x_values(1,1,jcell,icell)
      ycell= y_values(1,1,jcell,icell)
      call get_equilibrium_solution([xcell],[ycell],w_eq,one,one,one,one)
      call compute_primitive(function_values(1,1,jcell,icell,1:nvar),w,one,one,one,one)
      write(10,'(7(1PE12.5,1X))')xcell,ycell,(w(ivar)-w_eq(ivar),ivar=var,nvar)
    end do
  end do
  close(10)

end subroutine output_file

subroutine get_modes_from_nodes(nodes, u, size_x, size_y, order_x, order_y)
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y,order_x,order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::nodes,u
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x)::x
 
  ! internal variables
  integer::icell, jcell, i, j, xquad, yquad
  real(kind=8)::legendre
  real(kind=8)::xcell,ycell
  real(kind=8)::dx, dy
  u(:,:,:,:,:) = 0.0

  call gl_quadrature(x_quad,w_x_quad,order_x)
  call gl_quadrature(y_quad,w_y_quad,order_y)

  dx = boxlen_x/dble(size_x)
  dy = boxlen_y/dble(size_y)

  do icell=1,nx
    xcell=(dble(icell)-0.5)*dx
    do jcell=1,ny
      ycell=(dble(jcell)-0.5)*dy
      ! Loop over modes coefficients
      do i=1,mx
        do j=1,my
          ! Loop over quadrature points
          do xquad=1,mx ! for more general use n_x_quad...
            do yquad=1,my
              ! Quadrature point in physical space
              ! Perform integration using GL quadrature
              u(j,i,jcell,icell,1:nvar)=u(j,i,jcell,icell,1:nvar)+ &
              & 0.25*nodes(yquad,xquad,jcell,icell,1:nvar)* &
              & legendre(x_quad(xquad),i-1)* &
              & legendre(y_quad(yquad),j-1)* &
              & w_x_quad(xquad) * &
              & w_y_quad(yquad)
            end do
          end do
        end do
      end do
    end do
  end do
end subroutine get_modes_from_nodes

subroutine get_nodes_from_modes(modes,u,size_x, size_y, order_x, order_y)
  ! reconstruct u from modes
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y,order_x,order_y

  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::modes,u
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x)::x
  
  ! internal variables
  integer::icell, jcell, i, j, intnode, jntnode, ivar
  real::xcell,ycell,xquad,yquad
  real(kind=8)::legendre
  real(kind=8)::dx, dy
  u(:,:,:,:,:) = 0.0

  call gl_quadrature(x_quad,w_x_quad,order_x)
  call gl_quadrature(y_quad,w_y_quad,order_y)

  dx = boxlen_x/dble(size_x)
  dy = boxlen_y/dble(size_y)

  u(:,:,:,:,:) = 0.0

  do ivar = 1,nvar
    do icell=1,nx
      do jcell = 1,ny
        xcell=(dble(icell)-0.5)*dx
        ycell=(dble(jcell)-0.5)*dy
        do i=1,mx
          do j=1,my
            do intnode = 1,mx
              do jntnode = 1,my
                ! Loop over quadrature points
                u(j,i,jcell,icell,ivar) = u(j,i,jcell,icell,ivar) +&
                & modes(jntnode,intnode,jcell,icell,ivar)*&
                &legendre(x_quad(i),intnode-1)*&
                &legendre(y_quad(j),jntnode-1)
              end do
            end do
          end do
        end do
      end do
    end do
  end do

end subroutine get_nodes_from_modes

subroutine get_equilibrium_solution(x,y,w, size_x,size_y, order_x, order_y)
  ! get equilibrium solution for primitive variables
  use parameters_dg_2d
  integer::size_x,size_y, order_x, order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::w
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x)::x,y
 
  real(kind=8)::rho_0,p_0,g

  select case (nequilibrium)
  case(1)
    w(:,:,:,:,1) = exp(-(x+y))
    w(:,:,:,:,2) = 0
    w(:,:,:,:,3) = 0
    w(:,:,:,:,4) = exp(-(x+y))
  case(2)
    rho_0 = 1.21
    p_0 = 1
    g = 1
    w(:,:,:,:,1) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
    w(:,:,:,:,2) = 0
    w(:,:,:,:,3) = 0
    w(:,:,:,:,4) = p_0*exp(-(rho_0*g/p_0)*(x+y))
  case(3)
    w(:,:,:,:,:) = 0.0
  end select

end subroutine get_equilibrium_solution

subroutine evolve(u, x,y, u_eq, u_d, du_d, dudt_d)
  ! eat real nodal values, spit out nodal values
  USE ISO_C_BINDING
  use parameters_dg_2d
  implicit none
  type (c_ptr) :: u_d, du_d, dudt_d
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u,u_eq
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::x,y
  ! internal variables
  real(kind=8)::t,dt
  real(kind=8)::cmax, dx, dy, cs_max,v_xmax,v_ymax
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::dudt, w, w1, w2,w3,w4, delta_u,nodes
  integer::iter, n, i, j
  integer::snap_counter
  integer::var
  character(len=10)::filename

  dx = boxlen_x/dble(nx)
  dy = boxlen_y/dble(ny)
  delta_u(:,:,:,:,:) = 0.0
  !call get_modes_from_nodes(u,delta_u, nx, ny, mx, my)
  !call get_nodes_from_modes(delta_u,u, nx, ny, mx, my)
  !do i=1,1000
  call device_get_modes_from_nodes(u_d,du_d)
  !call device_get_nodes_from_modes(du_d,u_d)
  !end do
  !call d2h(du_d,u,14400)
  !open(unit=8,file='dudt.dat',form='UNFORMATTED',access='direct',recl=14400*8)
  !write(8,rec=1)u
  !stop
 
  t=0
  iter=0
  snap_counter = 0
  call apply_limiter(delta_u)
 
  do while(t < tend)
    ! Compute time step
    call device_compute_max_speed(cs_max,v_xmax,v_ymax,cmax)
    !write(*,*) 'cmax=',cmax
    write(*,*) 'csmax, umax,vmax= ', cs_max, v_xmax, v_ymax
    dt = (cfl/dble(2*4+1))/((abs(v_xmax)+(cs_max))/dx + (abs(v_ymax)+(cs_max))/dx)
    
    if(solver=='EQL')then

      call compute_update(delta_u,x,y,u_eq, dudt)
      w1=delta_u+dt*dudt
      call apply_limiter(w1)
      call compute_update(w1,x,y,u_eq,dudt)
      delta_u=0.5*delta_u+0.5*w1+0.5*dt*dudt
      call apply_limiter(delta_u)

    endif

    if(solver=='RK4')then !4th order SSP RGKT
      !call compute_update(delta_u,x,y,u_eq,dudt)
      call device_compute_update(0,dt)
      !call apply_limiter(w1)
      call device_compute_update(1,dt)
      !call apply_limiter(w2)
      call device_compute_update(2,dt)
      !call apply_limiter(w3)
      call device_compute_update(3,dt)
      !call apply_limiter(w4)
      call device_compute_update(4,dt)
      !call apply_limiter(delta_u)
      call d2h(dudt_d,dudt,14400)
      !open(unit=8,file='dudt.dat',form='UNFORMATTED',access='direct',recl=14400*8)
      !write(8,rec=1)dudt
      !write(*,*) "Here"
      !stop
    endif

    if(solver=='DEB')then
      call compute_update(delta_u,x,y,u_eq,dudt)
      w1=delta_u+0.391752226571890*dt*dudt
      write(*,*) 'maxval w1', maxval(w1)
      call apply_limiter(w1)
      write(*,*) 'dt:', dt
      write(*,*) 'dudt', maxval(dudt)
      write(*,*) maxval(w1)
    endif

    t=t+dt
    iter=iter+1
    write(*,*)'time=',iter,t,dt, cmax
    if(abs(cs_max) > 1.5)then
     pause
    endif
    call device_get_nodes_from_modes(du_d,u_d)
    if ((make_movie).and.(MODULO(iter,interval)==0)) then
      var = 1
      snap_counter = snap_counter + 1
      call d2h(u_d,nodes,nvar*nx*ny*mx*my)
      write(filename,'(a, i5.5)') 'SIM', snap_counter
      call output_file(x,y, nodes,var,filename)

    end if


  end do
  call d2h(du_d,delta_u,nvar*nx*ny*mx*my)
  call get_nodes_from_modes(delta_u, nodes, nx, ny, mx, my)
  u = nodes
end subroutine evolve

subroutine get_boundary_conditions(index, dim)
  use parameters_dg_2d

  integer::index
  integer::dim

  ! if direction = 1 then its x direction, etc.

  if (dim == 1) then
    index = index

    if (bc == 1) then !periodic
      if (index == 0) then
        index = nx
      else if (index == nx+1) then
        index = 1
      end if

    else if ((bc == 2).or.(bc==3)) then !transmissive or reflective
      if (index == 0) then
        index = 1
      else if (index == nx+1) then
        index = nx
      end if

    end if

  else if(dim == 2) then

    if (bc == 1) then !periodic
      if (index == 0) then
        index = ny
      else if (index == ny+1) then
        index = 1
      end if

    else if ((bc == 2).or.(bc==3)) then !transmissive
      if (index == 0) then
        index = 1
      else if (index == ny+1) then
        index = ny
      end if

    end if

  end if

end subroutine get_boundary_conditions

subroutine compute_max_speed(u, cs_max, v_xmax, v_ymax, speed_max)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u
  integer::icell, jcell, j,i
  real(kind=8)::speed, cs, v_x, v_y, cs_max, v_xmax, v_ymax, speed_max
  real(kind=8),dimension(1:nvar)::w
  ! Compute max sound speed
  speed_max=0.0
  do icell=1,nx
    do jcell = 1,ny
      do i=1,mx
        do j=1,my
          call compute_speed(u(j,i,jcell,icell,1:nvar),cs,v_x,v_y,speed)
          if (speed >= speed_max) then
            speed_max=MAX(speed_max,speed)
            v_xmax = v_x
            v_ymax = v_y
            cs_max = cs
            !cs_max = sqrt(1.4)
            call compute_primitive(u(j,i,jcell,icell,1:nvar),w,1,1,1,1)
          end if
        end do
      end do
    end do
  end do
end subroutine compute_max_speed


subroutine compute_speed(u,cs,v_x,v_y,speed)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar)::u
  real(kind=8)::speed
  real(kind=8),dimension(1:nvar)::w
  real(kind=8)::cs, v_x, v_y
  ! Compute primitive variables
  call compute_primitive(u,w,1,1,1,1)
  ! Compute sound speed

  cs=sqrt(gamma*max(w(4),1d-10)/max(w(1),1d-10))
  v_x = w(2)
  v_y = w(3)
  speed=sqrt(w(2)**2+w(3)**2)+cs
end subroutine compute_speed

subroutine compute_primitive(u,w,size_x,size_y,order_x,order_y)
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y, order_x, order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::u,w
  !Compute primitive variables
  w(:,:,:,:,1) = u(:,:,:,:,1)
  w(:,:,:,:,2) = u(:,:,:,:,2)/u(:,:,:,:,1)
  w(:,:,:,:,3) = u(:,:,:,:,3)/u(:,:,:,:,1)
  w(:,:,:,:,4) = (gamma-1.0)*( u(:,:,:,:,4) - 0.5*w(:,:,:,:,1)*(w(:,:,:,:,2)**2+w(:,:,:,:,3)**2) )
end subroutine compute_primitive


subroutine compute_conservative(ww,u,size_x,size_y,order_x,order_y)
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y, order_x, order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::u,ww
  ! Compute primitive variables
  u(:,:,:,:,1)= ww(:,:,:,:,1)
  u(:,:,:,:,2)= ww(:,:,:,:,1)*ww(:,:,:,:,2)
  u(:,:,:,:,3)= ww(:,:,:,:,1)*ww(:,:,:,:,3)
  u(:,:,:,:,4)= ww(:,:,:,:,4)/(gamma-1.) + 0.5*( ww(:,:,:,:,1)*(ww(:,:,:,:,2)**2+ww(:,:,:,:,3)**2) )
  !ww(3,:,:)/(gamma-1.0)+0.5*ww(1,:,:)*ww(2,:,:)**2
end subroutine compute_conservative

subroutine compute_flux(u,flux1, flux2, size_x,size_y,order_x,order_y)
  use parameters_dg_2d
  integer::size_x, size_y,order_x,order_y
  real(kind=8),dimension(1:order_y,1:order_x,1:size_y,1:size_x,1:nvar)::u,flux1,flux2,w
  ! Compute primitive variables
  call compute_primitive(u,w,size_x,size_y,order_x,order_y)
  ! Compute flux
  !    write(*,*) 'cons', maxval(u)

  !    write(*,*) 'prim', maxval(w)

  flux2(:,:,:,:,1)=w(:,:,:,:,1)*w(:,:,:,:,3)
  flux2(:,:,:,:,2)=w(:,:,:,:,1)*w(:,:,:,:,2)*w(:,:,:,:,3)
  flux2(:,:,:,:,3)=w(:,:,:,:,3)*u(:,:,:,:,3)+w(:,:,:,:,4)
  flux2(:,:,:,:,4)=w(:,:,:,:,3)*u(:,:,:,:,4)+w(:,:,:,:,3)*w(:,:,:,:,4)


  flux1(:,:,:,:,1)=w(:,:,:,:,1)*w(:,:,:,:,2)
  flux1(:,:,:,:,2)=w(:,:,:,:,2)*u(:,:,:,:,2)+w(:,:,:,:,4)
  flux1(:,:,:,:,3)=w(:,:,:,:,1)*w(:,:,:,:,2)*w(:,:,:,:,3)
  flux1(:,:,:,:,4)=w(:,:,:,:,2)*u(:,:,:,:,4)+w(:,:,:,:,2)*w(:,:,:,:,4)

end subroutine compute_flux

subroutine compute_flux_int(u,flux)
  use parameters_dg_2d
  real(kind=8),dimension(1:nvar)::u
  real(kind=8),dimension(1:nvar,2)::flux
  real(kind=8),dimension(1:nvar)::w
  ! Compute primitive variables
  call compute_primitive(u,w,1,1,1,1)
  ! Compute flux

  flux(1,1)=w(2)*u(1)
  flux(2,1)=w(2)*u(2)+w(4)
  flux(3,1)=w(1)*w(2)*w(3)
  flux(4,1)=w(2)*u(4)+w(2)*w(4)

  flux(1,2)=w(3)*u(1)
  flux(2,2)=w(1)*w(2)*w(3)
  flux(3,2)=w(3)*u(3)+w(4)
  flux(4,2)=w(3)*u(4)+w(3)*w(4)

end subroutine compute_flux_int


subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv, flag)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar)::uleft,uright, f_left, f_right
  real(kind=8),dimension(1:nvar)::fgdnv
  real(kind=8)::cleft,cs,cmax,speed_left,speed_right,v_x,v_y
  real(kind=8),dimension(1:nvar)::fleft,fright
  integer::flag
  ! Maximum wave speed
  call compute_speed(uleft,cs,v_x,v_y,speed_left)
  call compute_speed(uright,cs,v_x,v_y,speed_right)

  !subroutine compute_speed(u,cs,v_x,v_y,speed)
  cmax=max(speed_left,speed_right)

  ! Compute Godunox flux
  fgdnv=0.5*(f_right+f_left)+0.5*cmax*(uleft-uright)
end subroutine compute_llflux
!-----

subroutine compute_update(delta_u,x,y,u_eq,dudt)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::delta_u,dudt,u_eq
  real(kind=8),dimension(1:nvar,1, 1:mx, 1:nx+1, 1:ny)::F
  real(kind=8),dimension(1:nvar,1:mx, 1, 1:nx, 1:ny+1)::G
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::x,y
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::s, source_vol
  !===========================================================
  ! This routine computes the DG update for the input state u.
  !===========================================================
  real(kind=8),dimension(1:nx,1:ny, 1:mx, 1:my, 1:nx+1,1:ny+1)::x_faces, y_faces
  real(kind=8),dimension(1:nvar,1:mx,1:my)::u_quad,source_quad

  real(kind=8),dimension(1:nvar):: u_temp
  real(kind=8),dimension(1:nvar,2):: flux_temp
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u_delta_quad
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::flux_quad1
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::flux_quad2


  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::flux_vol1, flux_vol2
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1,1:my)::u_left,u_right
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1:mx,1)::u_top, u_bottom
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1:mx,1,2)::flux_top, flux_bottom

  real(kind=8),dimension(1:nvar)::flux_riemann,u_tmp
  real(kind=8),dimension(1:nvar,1:nx+1,1:ny+1)::flux_face,flux_face_eq
  !real(kind=8),dimension(1:nx+1,1:ny+1)::x_faces
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my,2)::flux_left,flux_right
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar,4):: edge

  integer::icell,i,j,iface,ileft,iright,ivar, node, jface
  integer::intnode,jntnode,jcell, edge_num, one
  real(kind=8)::legendre,legendre_prime
  real(kind=8)::chsi_left=-1,chsi_right=+1
  real(kind=8)::chsi_bottom=-1,chsi_top=+1
  real(kind=8)::dx,dy
  real(kind=8)::cmax,oneoverdx,c_left,c_right,oneoverdy
  real(kind=8)::x_right,x_left
  real(kind=8),dimension(1:nvar,1:mx):: u_delta_r, u_delta_l, u_delta_t, u_delta_b
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx)::u_delta_left,u_delta_right,u_delta_top,u_delta_bottom
  real(kind=8),dimension(1:nvar,1:(nx+1))::u_face_eq, w_face_eq

  dx=boxlen_x/dble(nx)
  dy=boxlen_y/dble(ny)
  oneoverdx=1./dble(dx)
  oneoverdy=1./dble(dy)

  call gl_quadrature(x_quad,w_x_quad,mx)
  call gl_quadrature(y_quad,w_y_quad,my)

  ! get equilibrium quantities

  !==================================
  ! Compute volume integral for Flux
  !==================================

  u_delta_quad(:,:,:,:,:) = 0.0
  flux_vol1(:,:,:,:,:) = 0.0
  flux_vol2(:,:,:,:,:) = 0.0
  !flux_quad(:,:,:,:,:,:) = 0.0

  call get_nodes_from_modes(delta_u,u_delta_quad,nx,ny,mx,my)
  call compute_flux(u_delta_quad,flux_quad1,flux_quad2,nx,ny,mx,my)


  do icell=1,nx
    do jcell=1,ny
      !==================================
      ! Compute flux at quadrature points
      !==================================
      !write(*,*) u_delta_quad
      !================================
      ! Compute volume integral DG term
      !================================
      ! Loop over modes
      do i = 1,mx
        do j = 1,my
          flux_vol1(j,i,jcell,icell,1:nvar) = 0.0
          do intnode=1,mx
            do jntnode=1,my
              flux_vol1(j,i,jcell,icell,1:nvar)=flux_vol1(j,i,jcell,icell,1:nvar)+ &
              & 0.25*flux_quad1(jntnode,intnode,jcell,icell,1:nvar)* & ! x direction
              & legendre_prime(x_quad(intnode),i-1)* &
              & w_x_quad(intnode)*&
              & legendre(y_quad(jntnode),j-1)* &
              & w_y_quad(jntnode)! * dx*dy/4.
            end do
          end do

          do intnode=1,mx
            do jntnode=1,my
              flux_vol2(j,i,jcell,icell,1:nvar) = flux_vol2(j,i,jcell,icell,1:nvar)&
              & + 0.25*flux_quad2(jntnode,intnode,jcell,icell,1:nvar)* & !y direction
              & legendre_prime(y_quad(jntnode),j-1)* &
              & w_y_quad(jntnode)*&
              & legendre(x_quad(intnode),i-1)* &
              & w_x_quad(intnode)! * dx*dy/4.
            end do
          end do
        end do
      end do
    end do
  end do


  u_delta_l(:,:) = 0.0
  u_delta_r(:,:) = 0.0
  u_delta_t(:,:) = 0.0
  u_delta_b(:,:) = 0.0


  do icell=1,nx
    do jcell=1,ny
      !==============================
      ! Compute left and right states
      ! computing the value AT the node -1 and 1
      !==============================
      !u_delta(1:nvar,icell)=0.0
      ! Loop over modes
      u_delta_l = 0.
      u_delta_r = 0.
      do i=1,mx
        do j=1,my
          do intnode = 1,my
            u_delta_l(1:nvar, intnode) = u_delta_l(1:nvar, intnode) + delta_u(j,i,jcell,icell,1:nvar)*&
            &legendre(chsi_left,i-1)*legendre(y_quad(intnode),j-1)
            u_delta_r(1:nvar, intnode) = u_delta_r(1:nvar, intnode) + delta_u(j,i,jcell,icell,1:nvar)*&
            &legendre(chsi_right,i-1)*legendre(y_quad(intnode),j-1)
          end do
        end do
      end do
      u_delta_left(1:nvar,icell,jcell,:) = u_delta_l(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.
      u_delta_right(1:nvar,icell,jcell,:) = u_delta_r(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.

      ! compute equilibrium on the fly
      u_left(1:nvar,icell,jcell,1,:) = u_delta_left(1:nvar,icell,jcell,:)
      u_right(1:nvar,icell,jcell,1,:) = u_delta_right(1:nvar,icell,jcell,:)
    end do
  end do

  do icell=1,nx
    do jcell=1,ny
      !==============================
      ! Compute left and right states
      ! computing the value AT the node -1 and 1
      !==============================
      !u_delta(1:nvar,icell)=0.0
      ! Loop over modes
      u_delta_b = 0.
      u_delta_t = 0.
      do i=1,mx
        do j=1,my
          do intnode = 1,mx
            !u_delta_b(1:nvar,intnode) = u_delta_b(1:nvar, intnode) +delta_u(j,i,jcell,icell,1:nvar)*&
            !&legendre(chsi_bottom,j-1)*legendre(x_quad(intnode),i-1)
            u_delta_b(1:nvar,intnode) = u_delta_b(1:nvar, intnode) + delta_u(j,i,jcell,icell,1:nvar)*&
            & legendre(chsi_bottom,j-1)*legendre(x_quad(intnode),i-1)

            u_delta_t(1:nvar,intnode) = u_delta_t(1:nvar, intnode) + delta_u(j,i,jcell,icell,1:nvar)*&
            & legendre(chsi_top,j-1)*legendre(x_quad(intnode),i-1)
          end do
        end do
      end do
      !write(*,*) legendre(chsi_left,1)
      !write(*,*) delta_u(1:nvar,icell,jcell,1,1)
      u_delta_bottom(1:nvar,icell,jcell,:) = u_delta_b(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.
      u_delta_top(1:nvar,icell,jcell,:) = u_delta_t(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.

      ! compute equilibrium on the fly
      u_bottom(1:nvar,icell,jcell,:,1) = u_delta_bottom(1:nvar,icell,jcell,:)
      u_top(1:nvar,icell,jcell,:,1) = u_delta_top(1:nvar,icell,jcell,:)
    end do
  end do

  ! fluxes
  F(:,:,:,:,:) = 0.
  G(:,:,:,:,:) = 0.
  one = 1
  do icell = 1,nx
    do jcell = 1,ny
      do i = 1,mx
        !do j = 1,my
        call compute_flux_int(u_left(1:nvar,icell,jcell,1,i),flux_left(1:nvar,icell,jcell,1,i,:))
        call compute_flux_int(u_right(1:nvar,icell,jcell,1,i),flux_right(1:nvar,icell,jcell,1,i,:))
        call compute_flux_int(u_top(1:nvar,icell,jcell,i,1),flux_top(1:nvar,icell,jcell,i,1,:))
        call compute_flux_int(u_bottom(1:nvar,icell,jcell,i,1),flux_bottom(1:nvar,icell,jcell,i,1,:))
        ! compute artificial flux in x direction
      end do
    end do
  end do

  do j = 1, ny
    do iface = 1, nx+1
      ileft = iface-1
      iright = iface

      call get_boundary_conditions(ileft,2)
      call get_boundary_conditions(iright,2)

      ! subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv)
      do intnode = 1, my
        call compute_llflux(u_right(1:nvar,ileft,j,1,intnode),u_left(1:nvar,iright,j,1,intnode),&
        &flux_right(1:nvar,ileft,j,1,intnode,1),&
        &flux_left(1:nvar,iright,j,1,intnode,1),F(1:nvar,1,intnode,iface,j),1)
      end do

    end do
  end do

  do i = 1,nx
    do jface = 1,ny+1
      ileft = jface-1
      iright = jface

      call get_boundary_conditions(ileft,2)
      call get_boundary_conditions(iright,2)

      do intnode = 1, mx
        call compute_llflux(u_top(1:nvar,i,ileft,intnode,1),u_bottom(1:nvar,i,iright,intnode,1),&
        &flux_top(1:nvar,i,ileft,intnode,1,2),&
        &flux_bottom(1:nvar,i,iright,intnode,1,2),G(1:nvar,intnode,1,i,jface),2)
      end do

    end do
  end do

  !========================
  ! Compute flux line integral
  !========================
  ! Loop over cells]
  edge(:,:,:,:,:,:)=0.0
  do icell = 1,nx
    do jcell = 1,ny
      do i = 1,mx
        do j = 1, my
          do intnode = 1,mx
            edge(j,i,jcell,icell,1:nvar, 1) = &
            &edge(j,i,jcell,icell,1:nvar, 1) + &
            &0.5*F(1:nvar,1, intnode, icell + 1, jcell)*legendre(chsi_right,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode)!&
            !&* dx/2.

            edge(j,i,jcell,icell,1:nvar,2) = &
            &edge(j,i,jcell,icell,1:nvar,2) + &
            &0.5*F(1:nvar,1, intnode, icell, jcell)*legendre(chsi_left,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode)!&
            !&* dx/2.
          end do
        end do
      end do

      do i = 1,mx
        do j =1,my
          do intnode = 1,my
            !edge(j,i,jcell,icell,1:nvar,3) = &
            !&edge(j,i,jcell,icell,1:nvar, 3) + &
            !&G(1:nvar, intnode, 1, icell, jcell+1)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)

            edge(j,i,jcell,icell,1:nvar,3) = &
            &edge(j,i,jcell,icell,1:nvar,3) + &
            &0.5*G(1:nvar,intnode, 1, icell, jcell+1)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)!&
            !&* dx/2.

            edge(j,i,jcell,icell,1:nvar,4) = &
            &edge(j,i,jcell,icell,1:nvar,4) + &
            &0.5*G(1:nvar, intnode, 1, icell, jcell)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)!&
            !&* dx/2.
          end do
        end do
      end do
    end do
  end do


  source_vol(:,:,:,:,:) = 0.0
  select case (source)
  case(1)
    !WRITE(*,*) 'No source'
  case(2) ! sine wave (tend=1 or 10)
    call get_source(u_delta_quad,s,x,y)
    ! evaluate integral
    do icell=1,nx
      do jcell = 1,ny
        do i=1,mx
          do j=1,my
            do intnode=1,mx
              do jntnode=1,my
                source_vol(j,i,jcell,icell,1:nvar) = source_vol(j,i,jcell,icell,1:nvar) + &
                & 0.25*s(jntnode,intnode,jcell,icell,1:nvar)* &
                & legendre(x_quad(intnode),i-1)* &
                & w_x_quad(intnode)*&
                & legendre(y_quad(jntnode),j-1)* &
                & w_y_quad(jntnode)
              end do
            end do
          end do
        end do
      end do
    end do
  end select

  !========================
  ! Compute final DG update
  !========================
  ! Loop over cells
  dudt(:,:,:,:,1:nvar) = 0.0
  do icell=1,nx
    do jcell=1,ny
      ! Loop over modes
      do i=1,mx
        do j=1,my
          dudt(j,i,jcell,icell,1:nvar) = &
          & 2*(oneoverdx*flux_vol1(j,i,jcell,icell,1:nvar) &
          & + oneoverdx*flux_vol2(j,i,jcell,icell,1:nvar) &
          &-oneoverdx*(edge(j,i,jcell,icell,1:nvar,1)&
          &-edge(j,i,jcell,icell,1:nvar,2)) &
          &-oneoverdx*(edge(j,i,jcell,icell,1:nvar,3)&
          &-edge(j,i,jcell,icell,1:nvar,4)) &
          & + source_vol(j,i,jcell,icell,1:nvar))
        end do
      end do
    end do
  end do
  !print*,'max dudt:',maxval(dudt)
  !print*,'min dudt:',minval(dudt)
  !dudt(:,:,:,:,4) = 0.0 ! block the update on pressure

end subroutine compute_update

subroutine apply_limiter(u)
  use parameters_dg_2d
  implicit none

  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u

  if(use_limiter) then
    if (limiter_type == '1OR') then
      call compute_limiter(u)
    else if (limiter_type == 'HIO') then
      call high_order_limiter(u)
    else if (limiter_type=='LOW') then
      call limiter_low_order(u)

    else if (limiter_type=='POS') then
      call limiter_positivity(u)

    else if (limiter_type=='ROS') then
      call limiter_rossmanith(u)

    end if

  end if

end subroutine apply_limiter


subroutine get_source(u,s,x,y)
  use parameters_dg_2d
  implicit none

  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u,s,w
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,2)::grad_p
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::x,y

  call compute_primitive(u,w,nx,ny,mx,my)
  call grad_phi(u,x,y,grad_p)

  s(:,:,:,:,1) = 0.
  s(:,:,:,:,2) = w(:,:,:,:,1)*grad_p(:,:,:,:,1)
  s(:,:,:,:,3) = w(:,:,:,:,1)*grad_p(:,:,:,:,2)
  s(:,:,:,:,4) = w(:,:,:,:,1)*(w(:,:,:,:,2)*grad_p(:,:,:,:,1)&
  &+w(:,:,:,:,3)*grad_p(:,:,:,:,2))

end subroutine get_source
!--------

subroutine grad_phi(u,x,y, grad_p)
  use parameters_dg_2d
  implicit none
  ! this could be a function by the way
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,1:nvar)::u
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx,2)::grad_p
  real(kind=8),dimension(1:my,1:mx,1:ny,1:nx)::x,y

  real(kind=8)::delta_r,x_center,y_center
  real(kind=8)::x_dash,y_dash, r, epsilon

  integer::i,j,icell,jcell

  epsilon = 0.25
  delta_r = 0.1
  x_center = 3.
  y_center = 3.

  do icell = 1,nx
    do jcell = 1,ny
      do i = 1,mx
        do j = 1,my
          x_dash = x(j,i,jcell,icell) - x_center
          y_dash = y(j,i,jcell,icell) - y_center
          r = sqrt(x_dash**2 + y_dash**2)

          if (r > 0.5-0.5*delta_r) then
            grad_p(j,i,jcell,icell,1) = -(x_dash)/(r**3)
            grad_p(j,i,jcell,icell,2) = -(y_dash)/(r**3)
          else if (r <= 0.5-0.5*delta_r) then
            grad_p(j,i,jcell,icell,1) = -(x_dash)/(r*(r**2+epsilon**2))
            grad_p(j,i,jcell,icell,2) = -(y_dash)/(r*(r**2+epsilon**2))
          end if

        end do
      end do
    end do
  end do





end subroutine


