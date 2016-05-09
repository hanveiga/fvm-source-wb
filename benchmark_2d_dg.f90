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
        w(1,:,:,:,:) = exp(-((x-0.5)**2+(y-0.5)**2)*10.)
        w(2,:,:,:,:) = 1.
        w(3,:,:,:,:) = 1.
        w(4,:,:,:,:) = 1. !exp(-((x-0.5)**2+(y-0.5)**2)*40.) !+ eta*&
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
          w(1,i,j,inti,intj) = 0.1072
          w(2,i,j,inti,intj) = -0.7259
          w(3,i,j,inti,intj) = -1.4045
          w(4,i,j,inti,intj) = 0.0439
        else if (x(i,j,inti,intj)>=0.5 .and. y(i,j,inti,intj)<0.5) then
        !else
          w(1,i,j,inti,intj) = 0.2579
          w(2,i,j,inti,intj) = 0.0
          w(3,i,j,inti,intj) = -1.4045
          w(4,i,j,inti,intj) = 0.15
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
      case(5)
        do i = 1,size_x
          do j = 1,size_y
            do inti = 1,order_x
              do intj = 1,order_y
        if (x(i,j,inti,intj) + y(i,j,inti,intj) >=0.5) then        
          w(1,i,j,inti,intj) = 1.
          w(2,i,j,inti,intj)  = 0.
          w(3,i,j,inti,intj)  = 0.
          w(4,i,j,inti,intj)  = 1.
        else if (x(i,j,inti,intj) + y(i,j,inti,intj) < 0.5) then
          w(1,i,j,inti,intj) = 0.125
          w(2,i,j,inti,intj) = 0.
          w(3,i,j,inti,intj) = 0.
          w(4,i,j,inti,intj) = 0.4
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

                  w(1,i,j,inti,intj) = 1.*(1. &
                  &+(gamma-1.)*5/(8*gamma*dpi**2)*exp(1-&
                  &((x(i,j,inti,intj)-5)**2+(y(i,j,inti,intj)-5)**2)))**(1/(gamma-1))

                  w(2,i,j,inti,intj)  = 1. + &
                  &5./(2*dpi)*exp(-1-((x(i,j,inti,intj)-5)**2+&
                  &(y(i,j,inti,intj)-5)**2)/2.)*(-y(i,j,inti,intj)+5.)

                  w(3,i,j,inti,intj)  = 1. + &
                  &5./(2*dpi)*exp(-1-((x(i,j,inti,intj)-5)**2+&
                  &(y(i,j,inti,intj)-5)**2)/2.)*&
                  &(x(i,j,inti,intj)-5.)
                  
                  w(4,i,j,inti,intj)  = w(1,i,j,inti,intj)**gamma
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
        ! Loop over quadrature points
        do xquad=1,mx ! for more general use n_x_quad...
          do yquad=1,my
           ! Quadrature point in physical space
           ! Perform integration using GL quadrature
           u(1:nvar,icell,jcell,i,j)=u(1:nvar,icell,jcell,i,j)+ &
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
            u(ivar,icell,jcell,i,j) = u(ivar, icell, jcell, i, j) +&
            & modes(ivar,icell,jcell,intnode,jntnode)*&
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
    real(kind=8)::cmax, dx, dy, cs_max,v_xmax,v_ymax
    real(kind=8),dimension(1:nvar,1:nx, 1:ny,1:mx,1:my)::dudt, w, w1, w2,w3,w4, delta_u,nodes
    integer::iter, n, i, j

    dx = boxlen_x/dble(nx)
    dy = boxlen_y/dble(ny)
    delta_u(:,:,:,:,:) = 0.0
    call get_modes_from_nodes(u,delta_u, nx, ny, mx, my)
    call get_nodes_from_modes(delta_u,nodes, nx, ny, mx, my)
    write(*,*) 'modes', delta_u(1,2,2,:,:)
    t=0
    iter=0

    call apply_limiter(delta_u)
    do while(t < tend)
    !do while(iter<1)
       ! Compute time step
       call compute_max_speed(nodes,cs_max,v_xmax,v_ymax,cmax)
       !dt=cfl*sqrt(dx*dy)/cmax/((2.0*dble(mx)+1.0)*(2.0*dble(my)+1.0))
        dt = (cfl/dble(2*mx+1))/((abs(v_xmax)+cs_max)/dx + (abs(v_ymax)+cs_max)/dx)
      if(solver=='EQL')then
        call compute_update(delta_u,u_eq, dudt)
        
        w1=delta_u+dt*dudt

        call compute_limiter(w1)
        call compute_update(w1,u_eq,dudt)
        delta_u=0.5*delta_u+0.5*w1+0.5*dt*dudt
        call compute_limiter(delta_u)

      endif


     if(solver=='RK4')then
        call compute_update(delta_u,u_eq,dudt)
        w1=delta_u+0.391752226571890*dt*dudt
        call apply_limiter(w1)
        call compute_update(w1,u_eq,dudt)
        w2=0.444370493651235*delta_u+0.555629506348765*w1+0.368410593050371*dt*dudt
        call apply_limiter(w2)
        call compute_update(w2,u_eq,dudt)
        w3=0.620101851488403*delta_u+0.379898148511597*w2+0.251891774271694*dt*dudt

        call apply_limiter(w3)
        call compute_update(w3,u_eq,dudt)
        w4=0.178079954393132*delta_u+0.821920045606868*w3+0.544974750228521*dt*dudt
        call apply_limiter(w4)
        delta_u=0.517231671970585*w2+0.096059710526147*w3+0.063692468666290*dt*dudt
        call compute_update(w4,u_eq,dudt)
        delta_u=delta_u+0.386708617503269*w4+0.226007483236906*dt*dudt
        call apply_limiter(delta_u)
     endif

       t=t+dt
       iter=iter+1
       write(*,*)'time=',iter,t,dt, cmax

       call get_nodes_from_modes(delta_u, nodes, nx, ny, mx, my)
    end do

    call get_nodes_from_modes(delta_u, nodes, nx, ny, mx, my)

    u = nodes
  end subroutine evolve

  !-----

  !subroutine get_characteristic_variables(du,dw,w)

  !  call compute_primitive()

  !end subroutine


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
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
    integer::icell, jcell, j,i
    real(kind=8)::speed, cs, v_x, v_y, cs_max, v_xmax, v_ymax, speed_max
    ! Compute max sound speed
    speed_max=0.0
    do icell=1,nx
      do jcell = 1,ny
       do i=1,mx
        do j=1,my
         call compute_speed(u(1:nvar,icell,jcell,i,j),cs,v_x,v_y,speed)
         if (speed >= speed_max) then
           speed_max=MAX(speed_max,speed)
           v_xmax = v_x
           v_ymax = v_y
           cs_max = cs
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



  subroutine compute_speed_2(u,speed, flag)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar)::u
    real(kind=8)::speed
    real(kind=8),dimension(1:nvar)::w
    real(kind=8)::cs
    integer::flag
    ! Compute primitive variables
    call compute_primitive(u,w,1,1,1,1)
    cs=sqrt(gamma*max(1.,1d-10)/max(w(1),1d-10))
    !write(*,*) cs
    !speed=max(abs(w(2)),abs(w(3)))+cs
    if (flag == 1) then
      speed=sqrt(w(2)**2)+cs
    else if (flag == 2) then
      speed=sqrt(w(3)**2)+cs
    end if
    !speed=100.
  end subroutine compute_speed_2


  subroutine compute_primitive(u,w,size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    implicit none
    integer::size_x,size_y, order_x, order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x,1:order_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y, 1:order_x,1:order_y)::w
    ! Compute primitive variables
    w(1,:,:,:,:) = u(1,:,:,:,:)
    w(2,:,:,:,:) = u(2,:,:,:,:)/u(1,:,:,:,:)
    w(3,:,:,:,:) = u(3,:,:,:,:)/u(1,:,:,:,:)
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

subroutine compute_flux(u,flux1, flux2, size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    integer::size_x, size_y,order_x,order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::flux1
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::flux2
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::w
    ! Compute primitive variables
    call compute_primitive(u,w,size_x,size_y,order_x,order_y)
    ! Compute flux

    flux2(1,:,:,:,:)=w(1,:,:,:,:)*w(3,:,:,:,:)
    flux2(2,:,:,:,:)=w(1,:,:,:,:)*w(2,:,:,:,:)*w(3,:,:,:,:)
    flux2(3,:,:,:,:)=w(3,:,:,:,:)*u(3,:,:,:,:)+w(4,:,:,:,:)
    flux2(4,:,:,:,:)=w(3,:,:,:,:)*u(4,:,:,:,:)+w(3,:,:,:,:)*w(4,:,:,:,:)


    flux1(1,:,:,:,:)=w(1,:,:,:,:)*w(2,:,:,:,:)
    flux1(2,:,:,:,:)=w(2,:,:,:,:)*u(2,:,:,:,:)+w(4,:,:,:,:)
    flux1(3,:,:,:,:)=w(1,:,:,:,:)*w(2,:,:,:,:)*w(3,:,:,:,:)
    flux1(4,:,:,:,:)=w(2,:,:,:,:)*u(4,:,:,:,:)+w(2,:,:,:,:)*w(4,:,:,:,:)

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

  subroutine get_source(u,s,size_x,size_y,order_x,order_y)
    use parameters_dg_2d
    implicit none
    integer::size_x, size_y, order_x,order_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::u,s
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,1:order_x,1:order_y)::w
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

  subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv, flag)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar)::uleft,uright, f_left, f_right
    real(kind=8),dimension(1:nvar)::fgdnv
    real(kind=8)::cleft,cright,cmax
    real(kind=8),dimension(1:nvar)::fleft,fright
    integer::flag
    ! Maximum wave speed
    call compute_speed_2(uleft,cleft,flag)
    call compute_speed_2(uright,cright,flag)
    cmax=max(cleft,cright)
    !if (cmax > 1000) then
    !  write(*,*) 'cmax', cmax
    !  write(*,*) uleft, uright
    !end if
    ! Compute Godunox flux
    fgdnv=0.5*(f_right+f_left)+0.5*cmax*(uleft-uright)
  end subroutine compute_llflux
  !-----

subroutine compute_update(delta_u,u_eq,dudt)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::delta_u,dudt
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::u_eq
  real(kind=8),dimension(1:nvar,1, 1:mx, 1:nx+1, 1:ny)::F
  real(kind=8),dimension(1:nvar,1:mx, 1, 1:nx, 1:ny+1)::G
  !===========================================================
  ! This routine computes the DG update for the input state u.
  !===========================================================
  real(kind=8),dimension(1:nx,1:ny, 1:mx, 1:my, 1:nx+1,1:ny+1)::x_faces, y_faces
  real(kind=8),dimension(1:nvar,1:mx,1:my)::u_quad,source_quad
  !real(kind=8),dimension(1:nvar,1:mx,1:my,2)::flux_quad
  real(kind=8),dimension(1:nvar,1:mx,1:my)::u_quad_eq,flux_quad_eq, source_quad_eq
  !real(kind=8),dimension(1:nvar,1:mx,1:my)::u_delta_quad
  real(kind=8),dimension(1:nvar):: u_temp
  real(kind=8),dimension(1:nvar,2):: flux_temp
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u_delta_quad
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::flux_quad1
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::flux_quad2

  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::flux_vol1, flux_vol2
  real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::flux_vol_eq, source_vol_eq
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1,1:my)::u_left,u_right
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1:mx,1)::u_top, u_bottom
  real(kind=8),dimension(1:nvar, 1:nx,1:ny,1:mx,1,2)::flux_top, flux_bottom

  real(kind=8),dimension(1:nvar)::flux_riemann,u_tmp
  real(kind=8),dimension(1:nvar,1:nx+1,1:ny+1)::flux_face,flux_face_eq
  !real(kind=8),dimension(1:nx+1,1:ny+1)::x_faces
  real(kind=8),dimension(1:nvar, 1:nx, 1:ny, 1, 1:my, 2)::flux_left,flux_right
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
        flux_vol1(1:nvar,icell,jcell,i,j) = 0.0
        do intnode=1,mx
          do jntnode=1,my
           flux_vol1(1:nvar,icell,jcell,i,j)=flux_vol1(1:nvar,icell,jcell,i,j)+ &
                & flux_quad1(1:nvar,icell,jcell,intnode,jntnode)* & ! x direction
                & legendre_prime(x_quad(intnode),i-1)* &
                & w_x_quad(intnode)*&
                & legendre(y_quad(jntnode),j-1)* &
                & w_y_quad(jntnode)
              end do 
            end do

        do intnode=1,mx
          do jntnode=1,my
            flux_vol2(1:nvar,icell,jcell,i,j) = flux_vol2(1:nvar,icell,jcell,i,j)&
                & + flux_quad2(1:nvar,icell,jcell,intnode,jntnode)* & !y direction
                & legendre_prime(y_quad(jntnode),j-1)* &
                & w_y_quad(jntnode)*&
                & legendre(x_quad(intnode),i-1)* &
                & w_x_quad(intnode)
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
            u_delta_l(1:nvar, intnode) = u_delta_l(1:nvar, intnode) + delta_u(1:nvar,icell,jcell,i,j)*&
            &legendre(chsi_left,i-1)*legendre(y_quad(intnode),j-1)
            u_delta_r(1:nvar, intnode) = u_delta_r(1:nvar, intnode) + delta_u(1:nvar,icell,jcell,i,j)*&
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
            !u_delta_b(1:nvar,intnode) = u_delta_b(1:nvar, intnode) +delta_u(1:nvar,icell,jcell,i,j)*&
            !&legendre(chsi_bottom,j-1)*legendre(x_quad(intnode),i-1)
            u_delta_b(1:nvar,intnode) = u_delta_b(1:nvar, intnode) + delta_u(1:nvar,icell,jcell,i,j)*&
            & legendre(chsi_bottom,j-1)*legendre(x_quad(intnode),i-1)

            u_delta_t(1:nvar,intnode) = u_delta_t(1:nvar, intnode) + delta_u(1:nvar,icell,jcell,i,j)*&
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

      call get_boundary_conditions(ileft,1)
      call get_boundary_conditions(iright,1)

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
               edge(1:nvar,icell,jcell,i,j, 1) = &
               &edge(1:nvar,icell,jcell,i,j, 1) + &
               &F(1:nvar,1, intnode, icell + 1, jcell)*legendre(chsi_right,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode)
               
               edge(1:nvar,icell,jcell,i,j,2) = &
               &edge(1:nvar,icell,jcell,i,j,2) + &
               &F(1:nvar,1, intnode, icell, jcell)*legendre(chsi_left,i-1)*legendre(x_quad(intnode),j-1)*w_x_quad(intnode)
             end do
            end do
        end do

          do i = 1,mx 
            do j =1,my
            do intnode = 1,my
               !edge(1:nvar,icell,jcell,i,j,3) = &
               !&edge(1:nvar,icell,jcell,i,j, 3) + &
               !&G(1:nvar, intnode, 1, icell, jcell+1)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)

               edge(1:nvar,icell,jcell,i,j,3) = &
               &edge(1:nvar,icell,jcell,i,j,3) + &
               &G(1:nvar,intnode, 1, icell, jcell+1)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)
   
               edge(1:nvar,icell,jcell,i,j,4) = &
               &edge(1:nvar,icell,jcell,i,j,4) + &
               &G(1:nvar, intnode, 1, icell, jcell)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)*w_x_quad(intnode)
            end do
            end do
        end do
    end do
  end do
  !write(*,*) maxval(edge(1:nvar,:,:,:,:,4)),maxval(edge(1:nvar,:,:,:,:,3))
  !write(*,*) 'end edge'


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
        dudt(1:nvar,icell,jcell,i,j) = & !oneoverdx*flux_vol(1:nvar,icell,jcell,i,j) &
          & oneoverdx*flux_vol1(1:nvar,icell,jcell,i,j) &
          & + oneoverdx*flux_vol2(1:nvar,icell,jcell,i,j) &
          &-oneoverdx*(edge(1:nvar,icell,jcell,i,j,1)&
          &-edge(1:nvar,icell,jcell,i,j,2)) &
          &-oneoverdx*(edge(1:nvar,icell,jcell,i,j,3)&
          &-edge(1:nvar,icell,jcell,i,j,4))
       end do
     end do
    end do
  end do
  !write(*,*) flux_vol(1,1,1,:,:)

end subroutine compute_update

subroutine apply_limiter(u)
  use parameters_dg_2d
  implicit none

  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u

  if(use_limiter) then
    if (limiter_type == '1OR') then
      call compute_limiter(u)
    else if (limiter_type == 'HIO') then
      call high_order_limiter(u)
    else if (limiter_type=='LOW') then
      call limiter_low_order(u)
    end if
  end if

end subroutine apply_limiter