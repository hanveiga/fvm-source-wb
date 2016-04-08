  !----
  ! global variables: parameters

  program main
    use parameters_dg_2d
    real(kind=8),dimension(1:nx,1:ny,1:mx,1:my)::x
    real(kind=8),dimension(1:nx,1:ny,1:mx,1:my)::y
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, w, u_eq
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::modes, nodes


    call get_coords(x,y,nx,ny, mx, my)

    call get_initial_conditions(x,y,u,nx,ny, mx, my)

    !call get_modes_from_nodes(u,modes,nx,ny, mx, my)

    !call get_nodes_from_modes(modes,nodes,nx,ny, mx, my)

    call output_file(x, y ,u,'initwo')

    !call output_file(x, y ,nodes,'nod')

    call get_equilibrium_solution(x, y, u_eq, nx, ny, mx, my)

    call evolve(u, u_eq)

    call output_file(x, y ,u,'fintwo')

  end program main

  !---------
  subroutine get_coords(x,y,size_x,size_y, order_x, order_y)
    use parameters_dg_2d
    integer::size_x,size_y, order_x, order_y
    real(kind=8),dimension(1:size_x,1:size_y,1:order_x,1:order_y)::x,y

    !internal vars
    real(kind=8)::dx, dy
    integer::i,j, nj, ni

    ! get quadrature
    call gl_quadrature(x_quad,w_x_quad,order_x)
    call gl_quadrature(y_quad,w_y_quad,order_y)

    dx = boxlen_x/dble(size_x)
    dy = boxlen_y/dble(size_y)
    do i=1,size_x
      do j=1,size_y
        do ni =1,order_x
          do nj = 1,order_y
            x(i, j, ni, nj) = (i-0.5)*dx + dx/2.0*x_quad(ni)
            y(i, j, ni, nj) = (j-0.5)*dy + dy/2.0*y_quad(nj)
          end do
        end do
      end do
    end do

  end subroutine get_coords
  !-------
  subroutine get_initial_conditions(x,y,u,size_x,size_y, order_x, order_y)
    use parameters_dg_2d
    integer::size_x,size_y, i ,j, order_x, order_y 
    real(kind=8),dimension(1:nvar,1:size_x, 1:size_y, 1:order_x, 1:order_y)::u
    real(kind=8),dimension(1:size_x, size_y,1:order_x, 1:order_y)::x
    real(kind=8),dimension(1:size_x,size_y,1:order_x, 1:order_y)::y
    integer::inti,intj
    ! internal variables
    real(kind=8),dimension(1:nvar,1:size_x, 1:size_y, 1:order_x, 1:order_y)::w
    real(kind=8)::rho_0, p_0, g
    real(kind=8)::dpi=acos(-1d0)

    select case (ninit)
      case(1)
        w(1,:,:,:,:) = exp(-(x+y))
        w(2,:,:,:,:) = 0
        w(3,:,:,:,:) = 0
        w(4,:,:,:,:) = exp(-(x+y)) !+ eta*&
                 !&exp(-100*((x-0.3)**2+(y-0.3)**2))
      case(2)
        rho_0 = 1.21
        p_0 = 1.
        g = 1.
        w(1,:,:,:,:) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
        w(2,:,:,:,:) = 0
        w(3,:,:,:,:) = 0
        w(4,:,:,:,:) = p_0*exp(-(rho_0*g/p_0)*(x+y))
      case(3)
        do i = 1,size_x
          do j = 1,size_y
            do inti = 1,order_x
              do intj = 1,order_y
        if (x(i,j,inti,intj)>=0.5 .and. y(i,j,inti,intj)>=0.5) then        
          w(1,i,j,inti,intj) = 1.
          w(2,i,j,inti,intj)  = 0.
          w(3,i,j,inti,intj)  = 0.
          w(4,i,j,inti,intj)  = 1.
        else if (x(i,j,inti,intj)<0.5 .and. y(i,j,inti,intj)>=0.5) then
          w(1,i,j,inti,intj) = 0.5197
          w(2,i,j,inti,intj) = -0.7259
          w(3,i,j,inti,intj) = 0.
          w(4,i,j,inti,intj) = 0.4
        else if (x(i,j,inti,intj)<0.5 .and. y(i,j,inti,intj)<0.5) then
          w(1,i,j,inti,intj) = 1.
          w(2,i,j,inti,intj) = -0.7259
          w(3,i,j,inti,intj) = -0.7259
          w(4,i,j,inti,intj) = 1.
        else if (x(i,j,inti,intj)>=0.5 .and. y(i,j,inti,intj)<0.5) then
        !else
          w(1,i,j,inti,intj) = 0.5197
          w(2,i,j,inti,intj) = 0.0
          w(3,i,j,inti,intj) = -0.7259
          w(4,i,j,inti,intj) = 0.4
        end if
      end do
    end do
    end do
    end do
      case(4)
        do i = 1,size_x
          do j = 1,size_y
            do inti = 1,order_x
              do intj = 1,order_y
        if (x(i,j,inti,intj)>=0.5 .and. y(i,j,inti,intj)>=0.5) then        
          w(1,i,j,inti,intj) = 1.5
          w(2,i,j,inti,intj)  = 0.
          w(3,i,j,inti,intj)  = 0.
          w(4,i,j,inti,intj)  = 1.5
        else if (x(i,j,inti,intj)<0.5 .and. y(i,j,inti,intj)>=0.5) then
          w(1,i,j,inti,intj) = 0.5323
          w(2,i,j,inti,intj) = 1.206
          w(3,i,j,inti,intj) = 0.
          w(4,i,j,inti,intj) = 0.3
        else if (x(i,j,inti,intj)<0.5 .and. y(i,j,inti,intj)<0.5) then
          w(1,i,j,inti,intj) = 0.138
          w(2,i,j,inti,intj) = 1.206
          w(3,i,j,inti,intj) = 1.206
          w(4,i,j,inti,intj) = 0.029
        else if (x(i,j,inti,intj)>=0.5 .and. y(i,j,inti,intj)<0.5) then
        !else
          w(1,i,j,inti,intj) = 0.5323
          w(2,i,j,inti,intj) = 0.0
          w(3,i,j,inti,intj) = 1.206
          w(4,i,j,inti,intj) = 0.3
        end if
      end do
    end do
    end do
    end do
    end select

    call compute_conservative(w,u,size_x,size_y,mx,my)
    !u = w
  end subroutine get_initial_conditions

  subroutine output_file(x_values, y_values, function_values, filen)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::function_values
    real(kind=8),dimension(1:nx,1:ny, 1:mx, 1:my)::x_values, y_values
    real(kind=8),dimension(1:nvar)::w, w_eq
    character(len=2)::filen

    ! internal vars
    integer::icell, ivar, jcell
    real(kind=8)::dx, xcell, ycell
    real(kind=8),dimension(1:nvar)::primitives
    integer::one
    !dx = boxlen/dble(nx)
    one = 1
  !  write(filen,"(A5,I5.5)")
    open(10,file=TRIM(filen)//".dat",form='formatted')
    do icell=1,nx
      do jcell=1,ny
       xcell= x_values(icell,jcell,1,1)
       ycell = y_values(icell,jcell,1,1)
       call get_equilibrium_solution([xcell],[ycell],w_eq,one,one,one,one)
       !w(1:nvar) = function_values(1:nvar,icell,jcell,1,1)
       call compute_primitive(function_values(1:nvar,icell,jcell,1,1),w,one,one,one,one)
       write(10,'(7(1PE12.5,1X))')xcell,ycell,(w(ivar),ivar=1,nvar)
      end do
    end do
    close(10)

  end subroutine output_file

  !------
  
  subroutine get_modes_from_nodes(nodes, u, size_x, size_y, order_x, order_y)
  use parameters_dg_2d
  implicit none
  integer::size_x,size_y,order_x,order_y
  real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::nodes
  real(kind=8),dimension(1:size_x,1:size_y,1:order_x,1:order_y)::x
  real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::u
 

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
        u(1:nvar,icell,jcell,i,j)=0.0
        ! Loop over quadrature points
        do xquad=1,mx ! for more general use n_x_quad...
          do yquad=1,my
           ! Quadrature point in physical space
           !x_val_quad=xcell+dx/2.0*x_quad(xquad)
           !y_val_quad=ycell+dy/2.0*y_quad(yquad)
           ! Perform integration using GL quadrature
           u(1:nvar,icell,jcell,i,j)=u(1:nvar,icell,jcell,i,j)+0.25* &
                & nodes(1:nvar,icell,jcell,xquad,yquad)* &
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

  real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::modes
  real(kind=8),dimension(1:size_x,1:size_y,1:order_x,1:order_y)::x
  real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::u
 

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

  do ivar = 1,nvar
    do icell=1,nx
      do jcell = 1,ny
       xcell=(dble(icell)-0.5)*dx
       ycell=(dble(jcell)-0.5)*dy
       do i=1,mx
         do j=1,my
         xquad=x_quad(i)
         yquad=y_quad(j)

         u(ivar,icell,jcell,i,j) = 0.0
         do intnode = 1,mx
          do jntnode = 1,my
          ! Loop over quadrature points
            u(ivar,icell,jcell,i,j) = u(ivar, icell, jcell, i, j) +&
            & modes(ivar,icell,jcell,intnode,jntnode)*&
            &legendre(xquad,intnode-1)*&
            &legendre(yquad,jntnode-1)
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
    ! use parameters
    ! communicate which initial solution, equilibrium type
    ! boundary conditions and domain properties
    use parameters_dg_2d
    integer::size_x,size_y, order_x, order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x, 1:order_y)::w
    real(kind=8),dimension(1:size_x,1:size_y, 1:order_x, 1:order_y)::x
    real(kind=8),dimension(1:size_x,1:size_y, 1:order_x, 1:order_y)::y
    
    !class(:)::x
    ! internal variables
    real(kind=8)::rho_0,p_0,g

    select case (nequilibrium)
      case(1)
        w(1,:,:,:,:) = exp(-(x+y))
        w(2,:,:,:,:) = 0
        w(3,:,:,:,:) = 0
        w(4,:,:,:,:) = exp(-(x+y))
      case(2)
        rho_0 = 1.21
        p_0 = 1
        g = 1
        w(1,:,:,:,:) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
        w(2,:,:,:,:) = 0
        w(3,:,:,:,:) = 0
        w(4,:,:,:,:) = p_0*exp(-(rho_0*g/p_0)*(x+y))
    end select

  end subroutine get_equilibrium_solution

  subroutine evolve(u, u_eq)
    ! eat real values, spit out real values
    use parameters_dg_2d
    implicit none

    real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx,1:my)::u,u_eq

    ! internal variables
    real(kind=8)::t,dt
    real(kind=8)::cmax, dx, dy
    real(kind=8),dimension(1:nvar,1:nx, 1:ny,1:mx,1:my)::dudt, w, w1, delta_u,nodes
    integer::iter, n, i, j

    dx = boxlen_x/dble(nx)
    dy = boxlen_y/dble(ny)
    delta_u(:,:,:,:,:) = 0.0
    call get_modes_from_nodes(u,delta_u, nx, ny, mx, my)

    t=0
    iter=0
    do while(t < tend)
    !do while(iter<2)
       ! Compute time step
       call compute_max_speed(delta_u,cmax)
       dt=0.5*dx/cmax*cfl
       
      if(solver=='EQL')then
         ! runge kutta 2nd order
        call compute_update(delta_u,u_eq, dudt)
        w1=delta_u+dt*dudt
        call compute_update(w1,u_eq,dudt)
        delta_u=0.5*delta_u+0.5*w1+0.5*dt*dudt
      endif

      !if(solver=='RKi')then

      !  call compute_update(delta_u,u_eq, dudt)
      !  w1=delta_u+0.391752226571890*dt*dudt
      !  call limiter_cons(delta_u)

      !  call compute_update(w1, u_eq, dudt)
      !  w2=0.444370493651235*delta_u+0.555629506348765*w1+0.368410593050371*dt*dudt
      !  call limiter_cons(w2)

      !  call compute_update(w2,u_eq, dudt)
      !  w3=0.620101851488403*delta_u+0.379898148511597*w2+0.251891774271694*dt*dudt
      !  call limiter_cons(w3)

      !  call compute_update(w3,u_eq, dudt)
      !  w4=0.178079954393132*delta_u+0.821920045606868*w3+0.544974750228521*dt*dudt

      !  delta_u=0.517231671970585*w2+0.096059710526147*w3+0.063692468666290*dt*dudt
      !  call limiter_cons(delta_u)

      !  call compute_update(w4,u_eq, dudt)
      !  delta_u=delta_u+0.386708617503269*w4+0.226007483236906*dt*dudt
      !  call limiter_cons(delta_u)
      !endif

       t=t+dt
       iter=iter+1
       write(*,*)'time=',iter,t,dt, cmax

       !call get_nodes_from_modes(delta_u, nodes, nx, ny, mx, my)
       !u = u_eq + nodes
    end do

    call get_nodes_from_modes(delta_u, nodes, nx, ny, mx, my)

    ! update u
    !u = u_eq + nodes
    u = nodes
  end subroutine evolve

  !-----

  subroutine compute_max_speed(u,cmax)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
    real(kind=8)::cmax
    integer::icell, jcell, j,i
    real(kind=8)::speed
    ! Compute max sound speed
    cmax=0.0
    do icell=1,nx
      do jcell = 1,ny
       do i=1,mx
        do j=1,my
         call compute_speed(u(1:nvar,icell,jcell,i,j),speed)
         cmax=MAX(cmax,speed)
        end do
       end do
      end do
    end do 
  end subroutine compute_max_speed


  subroutine compute_speed(u,speed)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar)::u
    real(kind=8)::speed
    real(kind=8),dimension(1:nvar)::w
    real(kind=8)::cs
    ! Compute primitive variables
    call compute_primitive(u,w,1,1,1,1)
    ! Compute sound speed
    cs=sqrt(gamma*max(w(4),1d-10)/max(w(1),1d-10))
    speed=sqrt(w(2)**2+w(3)**2)+cs
  end subroutine compute_speed



  subroutine compute_primitive(u,w,size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    implicit none
    integer::size_x,size_y, order_x, order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x,1:order_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x,1:order_y)::w
    ! Compute primitive variables
    w(1,:,:,:,:) = u(1,:,:,:,:)
    w(2,:,:,:,:) = u(2,:,:,:,:)/w(1,:,:,:,:)
    w(3,:,:,:,:) = u(3,:,:,:,:)/w(1,:,:,:,:)
    w(4,:,:,:,:) = (gamma-1.0)*( u(4,:,:,:,:) - 0.5*w(1,:,:,:,:)*(w(2,:,:,:,:)**2+w(3,:,:,:,:)**2) )
  end subroutine compute_primitive


  subroutine compute_conservative(ww,u,size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    implicit none
    integer::size_x,size_y, order_x, order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x, 1:order_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x, 1:order_y)::ww
    ! Compute primitive variables
    u(1,:,:,:,:)= ww(1,:,:,:,:)
    u(2,:,:,:,:)= ww(1,:,:,:,:)*ww(2,:,:,:,:)
    u(3,:,:,:,:)= ww(1,:,:,:,:)*ww(3,:,:,:,:)
    u(4,:,:,:,:)= ww(4,:,:,:,:)/(gamma-1.) + 0.5*( ww(1,:,:,:,:)*(ww(2,:,:,:,:)**2+ww(3,:,:,:,:)**2) )
    !ww(3,:,:)/(gamma-1.0)+0.5*ww(1,:,:)*ww(2,:,:)**2
  end subroutine compute_conservative

subroutine compute_flux(u,flux, size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    integer::size_x, size_y,order_x,order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y,2)::flux
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::w
    ! Compute primitive variables
    call compute_primitive(u,w,size_x,size_y,order_x,order_y)
    ! Compute flux

    flux(1,:,:,:,:,1)=w(2,:,:,:,:)*u(1,:,:,:,:)
    flux(2,:,:,:,:,1)=w(2,:,:,:,:)*u(2,:,:,:,:)+w(4,:,:,:,:)
    flux(3,:,:,:,:,1)=w(1,:,:,:,:)*w(2,:,:,:,:)*w(3,:,:,:,:)
    flux(4,:,:,:,:,1)=w(2,:,:,:,:)*u(4,:,:,:,:)+w(2,:,:,:,:)*w(4,:,:,:,:)

    flux(1,:,:,:,:,2) = u(1,:,:,:,:)*w(3,:,:,:,:)
    flux(2,:,:,:,:,2) = u(2,:,:,:,:)*w(3,:,:,:,:)
    flux(3,:,:,:,:,2) = u(3,:,:,:,:)*w(3,:,:,:,:)+w(4,:,:,:,:)
    flux(4,:,:,:,:,2) = w(3,:,:,:,:)*u(4,:,:,:,:)+w(3,:,:,:,:)*w(4,:,:,:,:)

  end subroutine compute_flux

  subroutine get_source(u,s,size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    implicit none
    integer::size_x, size_y, order_x,order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::u,s
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::w
    !real(kind=8),dimension(1:size_x,1:size_y)::x
    !real(kind=8),dimension(1:size_x,1:size_y)::y
    real(kind=8)::phi_x,phi_y
    !internal
    phi_x = 0.0
    phi_y = -0.1
    !x_minus = [x(1),x(1:size-1)] ! zero gradient
    !x_plus = [x(2:size),x(size)] ! zero gradient
    call compute_primitive(u,w,size_x,size_y,order_x,order_y)
    s(1,:,:,:,:) = 0.
    s(2,:,:,:,:) = -w(1,:,:,:,:)*phi_x
    s(3,:,:,:,:) = -w(1,:,:,:,:)*phi_y
    s(4,:,:,:,:) = -w(1,:,:,:,:)*(w(2,:,:,:,:)*phi_x+w(3,:,:,:,:)*phi_y)
    !s(1) = 0
    !  s(2) = -w(1)*1*(x_plus-x_minus)/(2*delta)
    !  s(3) = -w(1)*w(2)*1*(x_plus-x_minus)/(2*delta)

  end subroutine get_source
  !--------

  subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv)
    use parameters_dg_2d
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
    fgdnv=0.5*(f_right+f_left)+0.5*cmax*(uleft-uright)
  end subroutine compute_llflux
  !-----

subroutine compute_update(delta_u,u_eq,dudt)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::delta_u,dudt
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::u_eq
  real(kind=8),dimension(1:nvar,1:mx, 1:nx+1, 1:ny+1)::F,G
  !===========================================================
  ! This routine computes the DG update for the input state u.
  !===========================================================
  real(kind=8),dimension(1:nx,1:ny, 1:mx, 1:my, 1:nx+1,1:ny+1)::x_faces, y_faces
  real(kind=8),dimension(1:nvar,1:mx,1:my)::u_quad,source_quad
  real(kind=8),dimension(1:nvar,1:mx,1:my,2)::flux_quad
  real(kind=8),dimension(1:nvar,1:mx,1:my)::u_quad_eq,flux_quad_eq, source_quad_eq
  real(kind=8),dimension(1:nvar,1:mx,1:my)::u_delta_quad
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::flux_vol, source_vol
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::flux_vol_eq, source_vol_eq
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1:mx,1)::u_left,u_right, u_top, u_bottom
  real(kind=8),dimension(1:nvar)::flux_riemann,u_tmp
  real(kind=8),dimension(1:nvar,1:nx+1,1:ny+1)::flux_face,flux_face_eq
  !real(kind=8),dimension(1:nx+1,1:ny+1)::x_faces
  real(kind=8),dimension(1:nvar, 1:nx, 1:ny, 1:mx, 1, 2)::flux_left,flux_right, flux_top, flux_bottom
  real(kind=8),dimension(1:nvar, 1:nx, 1:ny, 1:mx, 1:my, 4):: edge

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

  !real(kind=8),dimension(1:nvar,1:n):: u_left_bc, u_left_bc_nodes, w_0_eq, w_1_eq
  !real(kind=8),dimension(1:nvar,1:n):: u_0_eq, u_1_eq, u_left_bc_nodes_1, u_left_modes_1
  !real(kind=8),dimension(1:n)::x_0_nodes, x_1_nodes
  !real(kind=8),dimension(1:nvar)::u_left_0, uu
  !real(kind=8),dimension(1:nvar)::w_eq_temp,w_eq_temp_0,u_eq_temp,u_eq_temp_0


  dx=boxlen_x/dble(nx)
  dy=boxlen_y/dble(ny)
  oneoverdx=1./dx
  oneoverdy=1./dy

  call gl_quadrature(x_quad,w_x_quad,mx)
  call gl_quadrature(y_quad,w_y_quad,my)

  ! get equilibrium quantities

  !==================================
  ! Compute volume integral for Flux
  !==================================
  do icell=1,nx
    do jcell=1,ny
     !==================================
     ! Compute flux at quadrature points
     !==================================
     ! Loop over quadrature points
     do i=1,mx
      do j=1,my
        u_delta_quad(1:nvar,i,j)=0.0
        ! Loop over modes
        do intnode=1,mx
          do jntnode=1,my
             u_delta_quad(1:nvar,i,j)=u_delta_quad(1:nvar,i,j)+&
             & delta_u(1:nvar,icell,jcell,i,j)* &
             & legendre(x_quad(i),intnode-1)* &
             & legendre(y_quad(j),jntnode-1)
          end do
        end do 
        ! Compute flux at quadrature points
        call compute_flux(u_delta_quad(1:nvar,i,j),flux_quad(1:nvar,i,j,:),1,1,1,1)
        !call compute_flux(u_eq(1:nvar,icell,jcell,i,j),&
        !  &flux_quad_eq(1:nvar,icell,jcell,i,j,:),1,1,1,1)
        !call compute_flux(u_eq(1:nvar,icell,jcell,i,j),flux_quad_eq(1:nvar,j),gamma,nvar)
      end do
     end do

     !================================
     ! Compute volume integral DG term
     !================================
     ! Loop over modes
     do i = 1,mx
      do j = 1,my
        flux_vol(1:nvar,icell,jcell,i,j)=0.0
        flux_vol_eq(1:nvar,icell,jcell,i,j)=0.0
        ! Loop over quadrature points
        do intnode=1,mx
          do jntnode=1,my
           flux_vol(1:nvar,icell,jcell,i,j)=flux_vol(1:nvar,icell,jcell,i,j)+ &
                & flux_quad(1:nvar,i,j,1)* & ! x direction
                & legendre_prime(x_quad(intnode),i-1)* &
                & w_y_quad(intnode)*&
                & legendre(y_quad(jntnode),j-1)* &
                & w_x_quad(jntnode)&
                & + flux_quad(1:nvar,i,j,2)* & !y direction
                & legendre_prime(y_quad(jntnode),j-1)* &
                & w_y_quad(jntnode)*&
                & legendre(x_quad(intnode),i-1)* &
                & w_x_quad(intnode)


           !flux_vol_eq(1:nvar,i,icell)=flux_vol_eq(1:nvar,i,icell)+ &
           !     & flux_quad_eq(1:nvar,j)* &
           !     & legendre_prime(chsi_quad(j),i-1)* &
           !     & w_quad(j)
          end do
        end do
       end do
      end do
    end do
  end do


!==================================
! Compute volume term of source
!==================================
  select case (source)
     case(1)
       source_vol(:,:,:,:,:) = 0.0
       source_vol_eq(:,:,:,:,:) = 0.0
     case(2)
       ! Compute source
       do icell=1,nx
       do jcell=1,ny
          !==================================
          ! Compute source at quadrature points
          !==================================
          ! Loop over quadrature points
          do i=1,mx
            do j=1,my
             u_quad(1:nvar,i,j)=0.0
             ! Loop over modes
             do intnode=1,mx
              do jntnode=1,my
                u_quad(1:nvar,i,j)=u_quad(1:nvar,i,j)+&
                &delta_u(1:nvar,icell,jcell,i,j)*&
                &legendre(x_quad(intnode),i-1)*&
                &legendre(y_quad(jntnode),j-1)
              end do
             end do
             ! Compute source at quadrature points
             call get_source(u_eq(1:nvar,icell,jcell,i,j)+u_quad(1:nvar,i,j),source_quad(1:nvar,i,j),1,1,1,1)
             call get_source(u_eq(1:nvar,icell,jcell,i,j),source_quad_eq(1:nvar,i,j),1,1,1,1)
            end do
          end do

          !================================
          ! Compute source volume integral
          !================================
          ! Loop over modes
          do i=1,mx
            do j=1,my 
             source_vol(1:nvar,icell,jcell,i,j)=0.0
             source_vol_eq(1:nvar,icell,jcell,i,j)=0.0
             do intnode=1,mx
              do jntnode=1,my
                source_vol(1:nvar,icell,jcell,i,j)=source_vol(1:nvar,icell,jcell,i,j)+ &
                     & source_quad(1:nvar,i,j)!* &
                     !& legendre(x_quad(intnode),i-1)*&
                     !& legendre(y_quad(jntnode),j-1)*&
                     !& w_x_quad(intnode)*0.5&
                     !& w_y_quad(jntnode)*0.5

                !source_vol(1:nvar,icell,jcell,i,j)=source_vol(1:nvar,icell,jcell,i,j)+ &
                !     & source_quad(1:nvar,i,j)* &
                !     & legendre(x_quad(intnode),i-1)*&
                !     & legendre(y_quad(jntnode),j-1)*&
                !     & w_x_quad(intnode)*0.5&
                !     & w_y_quad(jntnode)*0.5

              end do
            end do
          end do
        end do
      end do
      end do
  end select

! compute equilibrium function values at faces
    ! do i=1,nx+1
    !  do j = 1,ny+1
    !    do node = 1,mx
    !      do intnode = 1,my
    !        x_faces(node,intnode,i,j) = (i-1)*dx+dx/2.*x_quad(node)
    !      end do
    !    end do
    ! end do
    !end do

    !do i=1,nx+1
    !  do j = 1,ny+1
    !    do node = 1, mx
       !   do intnode = 1,my
      !       y_faces(node,intnode,i,j) = (j-1)*dx+dx/2.*y_quad(node)
     !     end do
    !    end do
   !   end do
   ! end do

    !call get_equilibrium_solution(x_faces, y, w_x_faces, nx+1, ny+1, mx,my)
    !call get_equilibrium_solution(x, y_faces, w_y_faces, nx+1, ny+1, mx,my)
    !call compute_conservative(w_x_faces, u_x_faces, nx+1, ny+1, mx,my)
   ! call compute_conservative(w_y_faces, u_y_faces, nx+1, ny+1, mx,my)

    !subroutine get_equilibrium_solution(x,y,w, size_x,size_y)
    !write(*,*) w_x_faces
    !write(*,*) u_x_faces
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
            u_delta_l(1:nvar, intnode) = u_delta_l(1:nvar, intnode) +delta_u(1:nvar,icell,jcell,i,j)*&
            &legendre(chsi_left,i-1)*legendre(y_quad(intnode),j-1)
            u_delta_r(1:nvar, intnode) = u_delta_r(1:nvar, intnode) +delta_u(1:nvar,icell,jcell,i,j)*&
            &legendre(chsi_right,i-1)*legendre(y_quad(intnode),j-1)
          end do
        end do
       end do
       u_delta_left(1:nvar,icell,jcell,:) = u_delta_l(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.
       u_delta_right(1:nvar,icell,jcell,:) = u_delta_r(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.
       
       ! compute equilibrium on the fly
       u_left(1:nvar,icell,jcell,:,1) = u_delta_left(1:nvar,icell,jcell,:)
       u_right(1:nvar,icell,jcell,:,1) = u_delta_right(1:nvar,icell,jcell,:)
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
            u_delta_b(1:nvar,intnode) = u_delta_b(1:nvar, intnode) +delta_u(1:nvar,icell,jcell,i,j)*&
            &legendre(chsi_bottom,j-1)*legendre(x_quad(intnode),i-1)
            u_delta_t(1:nvar,intnode) = u_delta_t(1:nvar, intnode) +delta_u(1:nvar,icell,jcell,i,j)*&
            & legendre(chsi_top,j-1)*legendre(x_quad(intnode),i-1)
          end do
        end do
       end do
       u_delta_bottom(1:nvar,icell,jcell,:) = u_delta_b(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.
       u_delta_top(1:nvar,icell,jcell,:) = u_delta_t(1:nvar,:) !(u_delta_r+u_delta_l)/2. - (uu_r+uu_l)/2.
       
       ! compute equilibrium on the fly
       u_bottom(1:nvar,icell,jcell,:,1) = u_delta_bottom(1:nvar,icell,jcell,:)
       u_top(1:nvar,icell,jcell,:,1) = u_delta_top(1:nvar,icell,jcell,:)
      end do
    end do

    !u_left(:,:,:) = u_x_faces(1:nvar,1:nx,1:ny) + delta_u(1:nvar,1:nx,1:ny)
    !u_right(:,:,:) = u_x_faces(1:nvar,2:(nx+1),1:ny) + delta_u(1:nvar,1:nx,1:ny)

    !u_top(:,:,:) = u_y_faces(1:nvar,1:nx,2:(ny+1)) + delta_u(1:nvar,1:nx,1:ny)
    !u_bottom(:,:,:) = u_y_faces(1:nvar,1:nx,1:ny) + delta_u(1:nvar,1:nx,1:ny)

    !write(*,*) 'delta_u'
    !write(*,*) delta_u

    ! fluxes
    F(:,:,:,:) = 0.
    G(:,:,:,:) = 0.
    one = 1
    call compute_flux(u_left,flux_left,nx,ny,mx,one)
    call compute_flux(u_right,flux_right,nx,ny,mx,one)
    call compute_flux(u_top,flux_top,nx,ny,mx,one)
    call compute_flux(u_bottom,flux_bottom,nx,ny,mx,one)
    ! compute artificial flux in x direction
    
    do j = 1, ny
    do iface = 1, nx+1
      ileft = iface-1
      iright = iface
      if (iface == 1) then
        ileft = 1
      end if
      if (iface == nx+1) then
        iright = ny
      end if

        ! subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv)
      do intnode = 1, my
      call compute_llflux(u_right(1:nvar,ileft,j,intnode,1),u_left(1:nvar,iright,j,intnode,1),&
                            &flux_right(1:nvar,ileft,j,intnode,1,1),&
                            &flux_left(1:nvar,iright,j,intnode,1,1),F(1:nvar,intnode,iface,j))
      end do

    end do
    end do 


    do i = 1,nx
      do jface = 1,ny+1
        ileft = jface-1
        iright = jface
        if (jface == 1) then
          ileft = 1
        end if
        if (jface == ny+1) then
          iright = nx
        end if

       !call compute_llflux(u_top(1:nvar,i,ileft),u_bottom(1:nvar,i,iright),&
       !                     &flux_top(1:nvar,i,ileft,2),flux_bottom(1:nvar,i,iright,2),G(1:nvar,i,jface))
      
      do intnode = 1, my
      call compute_llflux(u_top(1:nvar,i,ileft,intnode,1),u_bottom(1:nvar,i,iright,intnode,1),&
                            &flux_top(1:nvar,i,ileft,intnode,1,2),&
                            &flux_bottom(1:nvar,i,iright,intnode,1,2),G(1:nvar,intnode,i,jface))
      end do


      end do
    end do
  !write(*,*) 'F'
  !write(*,*) F



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
               edge(1:nvar,icell,jcell,i,j, 1) = &
               &edge(1:nvar,icell,jcell,i,j, 1) + &
               &0.5*F(1:nvar, intnode, icell + 1, jcell)*legendre(chsi_right,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode)
               
               edge(1:nvar,icell,jcell,i,j,2) = &
               &edge(1:nvar,icell,jcell,i,j,2) + &
               &0.5*F(1:nvar, intnode, icell, jcell)*legendre(chsi_left,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode)
             end do
            end do
        end do

          do i = 1,mx
            do j =1,mx
            do intnode = 1,mx
               edge(1:nvar,icell,jcell,i,j,3) = &
               &edge(1:nvar,icell,jcell,i,j, 3) + &
               &0.5*G(1:nvar, intnode, icell, jcell+1)*legendre(chsi_top,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)
   
               edge(1:nvar,icell,jcell,i,j,4) = &
               &edge(1:nvar,icell,jcell,i,j,4) + &
               &0.5*G(1:nvar, intnode, icell, jcell)*legendre(chsi_bottom,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)
            end do
            end do
        end do
    end do
  end do


  !========================
  ! Compute final DG update
  !========================
  ! Loop over cells
  dudt(1:nvar,:,:,:,:) = 0.0
  do icell=1,nx
    do jcell=1,ny
     ! Loop over modes
     do i=1,mx
      do j=1,my
        dudt(1:nvar,icell,jcell,i,j) =  -(edge(1:nvar,icell,jcell,i,j,1)&
          &-edge(1:nvar,icell,jcell,i,j,2))*oneoverdx &
          &-(edge(1:nvar,icell,jcell,i,j,3)-edge(1:nvar,icell,jcell,i,j,4))*oneoverdx

             ! dudt(1:nvar,icell,jcell,i,j) &
             !& -oneoverdx*(F(1:nvar,j,icell+1,jcell)*&
            ! & legendre(chsi_right,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode) &
            ! & -F(1:nvar,j,icell,jcell)*&
            ! & legendre(chsi_left,i-1)*legendre(x_quad(intnode),j-1))*w_x_quad(intnode)  &
            ! & -oneoverdy*(G(1:nvar,i,icell,jcell+1)*&
            ! & legendre(x_quad(intnode),i-1)*legendre(chsi_top,j-1)*w_x_quad(intnode)  &
            ! & -G(1:nvar,i,icell,jcell)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)  &
            ! & *legendre(chsi_bottom,j-1)) 
               !& oneoverdx*oneoverdy*flux_vol(1:nvar,icell,icell,i,j) &

!             & + source_vol(1:nvar,icell,jcell,i,j)
        ! end do
       end do
     end do
    end do
  end do

!  dudt(1:nvar,1,:,:,:) = 0.0
!  dudt(1:nvar,nx,:,:,:) = 0.0
!  dudt(1:nvar,:,1,:,:) = 0.0
!    dudt(1:nvar,:,ny,:,:) = 0.0

end subroutine compute_update


