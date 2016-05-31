  !----
  ! global variables: parameters

  program main
    use parameters_rd
    real(kind=8),dimension(1:nx,1:ny,1:mx,1:my)::x
    real(kind=8),dimension(1:nx,1:ny,1:mx,1:my)::y
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, w, u_eq

    call get_coords(x,y,nx,ny, mx, my)

    call get_initial_conditions(x,y,u,nx,ny, mx, my)
    call output_file(x, y ,u,1,'initcond')

    call get_equilibrium_solution(x, y, u_eq, nx, ny, mx, my)

    call evolve(u, x, y,u_eq)

    !call output_file(x, y ,u,'fintwo')

  end program main

  !---------
  subroutine get_coords(x,y,size_x,size_y, order_x, order_y)
    use parameters_rd
    implicit none
    integer::size_x,size_y, order_x, order_y
    real(kind=8),dimension(1:size_x,1:size_y,1:order_x,1:order_y)::x,y

    !internal vars
    real(kind=8)::dx, dy
    integer::i,j, nj, ni

    ! get quadrature
    !x_quad = [-1,1]
    !y_quad = [-1,1]

    dx = boxlen_x/dble(size_x-1)
    dy = boxlen_y/dble(size_y-1)
    do i=1,size_x
      do j=1,size_y
        do ni =1,order_x
          do nj = 1,order_y
            x(i, j, ni, nj) = (i-1)*dx !(i-0.5)*dx + dx/2.0*x_quad(ni)
            y(i, j, ni, nj) = (j-1)*dy !(j-0.5)*dy + dy/2.0*y_quad(nj)
          end do
        end do
      end do
    end do

  end subroutine get_coords
  !-------
  subroutine get_initial_conditions(x,y,u,size_x,size_y, order_x, order_y)
    use parameters_rd
    implicit none
    integer::size_x,size_y, i ,j, order_x, order_y 
    real(kind=8),dimension(1:nvar,1:size_x, 1:size_y, 1:order_x, 1:order_y)::u
    real(kind=8),dimension(1:size_x, size_y,1:order_x, 1:order_y)::x
    real(kind=8),dimension(1:size_x,size_y,1:order_x, 1:order_y)::y
    integer::inti,intj, jcell, icell
    ! internal variables
    real(kind=8),dimension(1:nvar,1:size_x, 1:size_y, 1:order_x, 1:order_y)::w
    real(kind=8)::rho_0, p_0, g
    real(kind=8)::dpi=acos(-1d0)

    ! smooth rotating disk
    real(kind=8)::x_dash, y_dash, r, rho_d, delta_r, x_center, y_center

    select case (ninit)
      case(1)
        w(1,:,:,:,:) = 1+exp(-((x(:,:,:,:)-boxlen_x/2.)**2+(y(:,:,:,:)-boxlen_y/2.)**2)*20)
        w(2,:,:,:,:) = 1.
        w(3,:,:,:,:) = 1.
        w(4,:,:,:,:) = 1.!+exp(-((x(:,:,:,:)-boxlen_x/2.)**2+(y(:,:,:,:)-boxlen_x/2.)**2)*20.) !+ eta*&
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
                  &-(gamma-1.)*5/(8*gamma*dpi**2)*exp(1-&
                  &((x(i,j,inti,intj)-5)**2+(y(i,j,inti,intj)-5)**2)))**(1/(gamma-1))

                  w(2,i,j,inti,intj) = 2+ &
                  &5./(2*dpi)*exp(-1-((x(i,j,inti,intj)-5)**2+&
                  &(y(i,j,inti,intj)-5)**2)/2.)*(-y(i,j,inti,intj)+5.)

                  w(3,i,j,inti,intj)  = 2+ &
                  &5./(2*dpi)*exp(-1-((x(i,j,inti,intj)-5)**2+&
                  &(y(i,j,inti,intj)-5)**2)/2.)*&
                  &(x(i,j,inti,intj)-5.)
                  
                  w(4,i,j,inti,intj)  = w(1,i,j,inti,intj)**gamma
      end do
    end do
    end do
    end do
    case(7) ! smooth rotating disk
      ! initialize variables for smooth rotating disk
      p_0 = 10e-5
      rho_0 = 10e-5
      rho_d = 1
      delta_r = 0.1
      x_center = 3.
      y_center = 3.

      w(:,:,:,:,:) = 0.
      w(1,:,:,:,:) = 1.

      w(4,:,:,:,:)  = p_0
      do i = 1,size_x
          do j = 1,size_y
            do inti = 1,order_x
              do intj = 1,order_y
                x_dash = x(i,j,inti,intj) - x_center
                y_dash = y(i,j,inti,intj) - y_center
                r = sqrt(x_dash**2 + y_dash**2)

                if (r<0.5-delta_r/2.) then
                  w(1,i,j,inti,intj) = rho_0
                else if ((r<0.5+delta_r/2.).and.(r >= 0.5-delta_r/2.)) then
                !else if ((r<=1+delta_r/2.).and.(r>1-delta_r/2.)) then
                  w(1,i,j,inti,intj) = (rho_d-rho_0)/delta_r * (r - (0.5-delta_r/2.)) + rho_0
                else if ((r >= 0.5+delta_r/2.).and.(r <= 2-delta_r/2.))  then
                  w(1,i,j,inti,intj) = rho_d
                else if ((r > 2-delta_r/2.).and.(r < 2+delta_r/2.))  then
                  w(1,i,j,inti,intj) = (rho_0-rho_d)/delta_r * (r - (2-delta_r/2.)) + rho_d
                else if (r >= 2+delta_r/2.) then
                  w(1,i,j,inti,intj) = rho_0
                end if


                !if ((r > 0.5 - 2*delta_r).and.(r < 2+2*delta_r)) then
                !  w(2,i,j,inti,intj) = -y_dash/r**(3./2.)
                !  w(3,i,j,inti,intj) = x_dash/r**(3./2.)
                !else
                !  w(2,i,j,inti,intj) = 0.
                !  w(3,i,j,inti,intj) = 0.
                if (r<=0.5-delta_r) then
                  w(2,i,j,inti,intj) = 0.
                  w(3,i,j,inti,intj) = 0.
                else if ((r<=0.5+delta_r).and.(r > 0.5-delta_r)) then
                  w(2,i,j,inti,intj) = -y_dash/r**(3./2.)/delta_r * (r - (0.5-delta_r))
                  w(3,i,j,inti,intj) = x_dash/r**(3./2.)/delta_r * (r - (0.5-delta_r))
                else if ((r > 0.5+delta_r).and.(r <= 2-delta_r))  then
                  w(2,i,j,inti,intj) = -y_dash/r**(3./2.)
                  w(3,i,j,inti,intj) = x_dash/r**(3./2.)
                else if ((r > 2-delta_r).and.(r <= 2+delta_r))  then
                  w(2,i,j,inti,intj) = y_dash/r**(3./2.)/delta_r * (r - (2+delta_r)) !- y_dash/r**(3./2.)
                  w(3,i,j,inti,intj) = -x_dash/r**(3./2.)/delta_r * (r - (2+delta_r)) !+ x_dash/r**(3./2.)
                else if (r > 2+delta_r) then
                  w(2,i,j,inti,intj) = 0.
                  w(3,i,j,inti,intj) = 0.
                end if
                  
            end do
          end do
        end do
      end do

      case(8)

        w(1,:,:,:,:) = 1.21*exp(-1.21*(x+y))
        w(2,:,:,:,:) = 0
        w(3,:,:,:,:) = 0
        w(4,:,:,:,:) = exp(-1.21*(x+y)) + eta*&
                 &exp(-100*1.21*((x-0.3)**2+(y-0.3)**2))

    end select

    call compute_conservative(w,u,nx,ny,mx,my)
    !u = w
  end subroutine get_initial_conditions

  subroutine output_file(x_values, y_values, function_values, var, filen)
    use parameters_rd
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx, 1:my)::function_values
    real(kind=8),dimension(1:nx,1:ny, 1:mx, 1:my)::x_values, y_values
    real(kind=8),dimension(1:nvar)::w, w_eq
    character(len=*)::filen

    ! internal vars
    integer::icell, ivar, jcell
    integer::var
    real(kind=8)::dx, xcell, ycell
    real(kind=8),dimension(1:nvar)::primitives
    integer::one
    !dx = boxlen/dble(nx)
    one = 1
  !  write(filen,"(A5,I5.5)")
    open(10,status='REPLACE',file="longrotation2_o1/"//TRIM(filen)//".dat",form='formatted')
    do icell=1,nx
      do jcell=1,ny
       xcell= x_values(icell,jcell,1,1)
       ycell = y_values(icell,jcell,1,1)
       call get_equilibrium_solution([xcell],[ycell],w_eq,one,one,one,one)
       call compute_primitive(function_values(1:nvar,icell,jcell,1,1),w,one,one,one,one)
       write(10,'(7(1PE12.5,1X))')xcell,ycell,(w(ivar),ivar=var,nvar)
      end do
    end do
    close(10)

  end subroutine output_file

  !------

  subroutine get_equilibrium_solution(x,y,w, size_x,size_y, order_x, order_y)
    ! get equilibrium solution for primitive variables
    ! use parameters
    ! communicate which initial solution, equilibrium type
    ! boundary conditions and domain properties
    use parameters_rd
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

  subroutine evolve(u, x,y, u_eq)
    ! eat real values, spit out real values
    use parameters_rd
    implicit none

    real(kind=8),dimension(1:nvar,1:nx,1:ny, 1:mx,1:my)::u,u_eq, diff, dudt, up, ua
    real(kind=8),dimension(1:nx,1:ny, 1:mx,1:my)::x,y
    ! internal variables
    real(kind=8)::t,dt
    real(kind=8)::cmax, dx, dy, cs_max,v_xmax,v_ymax
    real(kind=8),dimension(1:nvar,1:nx, 1:ny,1:mx,1:my):: w, w1, w2,w3,w4, delta_u,nodes
    integer::iter, n, i, j

    integer::snap_counter
    integer::var
    character(len=10)::filename

    dx = boxlen_x/dble(nx-1)
    dy = boxlen_y/dble(ny-1)

    t=0
    iter=0
    snap_counter = 0
    do while(t < tend)
    !do while(iter<100)
      CALL compute_max_speed(u,cmax)
      dt=MIN(cfl*dy*dx/cmax/(2.0*DBLE(1)+1.0),tend-t)

      IF(solver=='RD1')THEN
            diff=0
            CALL compute_update_rd(u,diff, dt,x,y,dudt)
            !write(*,*) 'dudt',dudt(3,:,:,:,:)
            !write(*,*) 'u',u(3,:,:,:,:)
            up=u-dudt/(dx*dy) ! this has to be changed for sure
            !ua=0.5*(up+u)
            diff=up-u
            CALL compute_update_rd(u,diff, dt,x,y,dudt)
            u=up-dudt/(dx*dy)
      ENDIF

      t=t+dt
      iter=iter+1
      
     if ((make_movie).and.(MODULO(iter,100)==0)) then
          var = 1
          snap_counter = snap_counter + 1
          write(filename,'(a, i5.5)') 'SIM', snap_counter
          call output_file(x,y, u,var,filename)    
      end if

      write(*,*)'time=',iter,t,dt, cmax


    end do

  end subroutine evolve

  !-----

  !subroutine get_characteristic_variables(du,dw,w)

  !  call compute_primitive()

  !end subroutine


  subroutine get_boundary_conditions(index, dim)
    use parameters_rd

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

  subroutine compute_max_speed(u,speed_max)
    use parameters_rd
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
    integer::icell, jcell, j,i
    real(kind=8)::speed, cs, v_x, v_y, cs_max, v_xmax, v_ymax, speed_max
    real(kind=8),dimension(1:nvar)::w
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
           call compute_primitive(u(1:nvar,icell,jcell,i,j),w,1,1,1,1)
           !write(*,*) 'w1',w(1)
           !write(*,*) 'w4',w(4)
         end if
        end do
       end do
      end do
    end do 
  end subroutine compute_max_speed


  subroutine compute_speed(u,cs,v_x,v_y,speed)
    use parameters_rd
    implicit none
    real(kind=8),dimension(1:nvar)::u
    real(kind=8)::speed
    real(kind=8),dimension(1:nvar)::w
    real(kind=8)::cs, v_x, v_y
    ! Compute primitive variables
    call compute_primitive(u,w,1,1,1,1)
    ! Compute sound speed
    !write(*,*) 'w',w
    cs=sqrt(gamma*max(w(4),1d-10)/max(w(1),1d-10))
    v_x = w(2)
    v_y = w(3)
    speed=sqrt(w(2)**2+w(3)**2)+cs
  end subroutine compute_speed



  subroutine compute_speed_2(u,speed, flag)
    use parameters_rd
    implicit none
    real(kind=8),dimension(1:nvar)::u
    real(kind=8)::speed
    real(kind=8),dimension(1:nvar)::w
    real(kind=8)::cs
    integer::flag
    ! Compute primitive variables
    call compute_primitive(u,w,1,1,1,1)
    cs=sqrt(gamma*max(w(4),1d-10)/max(w(1),1d-10))
    !cs=sqrt(gamma*w(4)/w(1))
    
    !write(*,*) 'cs,', cs
    !speed=max(abs(w(2)),abs(w(3)))+cs
    if (flag == 1) then
      speed=sqrt(w(2)**2)+cs
    else if (flag == 2) then
      speed=sqrt(w(3)**2)+cs
    end if
    !speed=100.
  end subroutine compute_speed_2


  subroutine compute_primitive(u,w,size_x,size_y,order_x,order_y)
    use parameters_rd
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
    use parameters_rd
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
    use parameters_rd
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
    use parameters_rd
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
    use parameters_rd
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
    ! Compute Godunox flux
    fgdnv=0.5*(f_right+f_left)+0.5*cmax*(uleft-uright)
  end subroutine compute_llflux
  !-----

subroutine compute_update_rd(u,diff, dt, x,y, dudt)

  USE parameters_rd
  IMPLICIT NONE
  REAL (kind=8), INTENT(in):: dt
  REAL(kind=8),DIMENSION(1:nvar,1:nx,1:ny,1:mx,1:my), INTENT(in)::u, diff
  REAL(kind=8),DIMENSION(1:nx,1:ny,1:mx,1:my), INTENT(in)::x, y
  REAL(kind=8),DIMENSION(1:nvar,1:nx,1:ny,1:mx,1:my), INTENT(out)::  dudt

  ! internal vars
  REAL(kind=8), DIMENSION(1:nvar):: uu11, uu12, uu21,uu22, ubar
  REAL(kind=8), DIMENSION(1:nvar):: ww11, ww12, ww21,ww22
  REAL(kind=8), DIMENSION(1:nvar)::ff11, ff12, ff21, ff22
  REAL(kind=8), DIMENSION(1:nvar)::g11, g12, g21, g22
  REAL(kind=8), DIMENSION(1:nvar):: bar,totalphi

  REAL(kind=8), DIMENSION(1:nvar,4)::phi, ph
  REAL(kind=8), DIMENSION(1:nvar):: ww1, ww2, w_eq_0, w_eq_1
  REAL(kind=8)                   ::speed1, speed2, speed3, speed4, dx,alpha, voldx, a,b,c
  REAL(kind=8),DIMENSION(1:nvar) :: source11, source21, source12, source22
  REAL(kind=8)::dy
  

  INTEGER:: i1, i2, one, nfaces, i, j, i11, i12, i22, i21
  INTEGER, DIMENSION(2,nx-1):: nubox
  INTEGER, DIMENSION(2,ny-1):: nuboy

  !dx = boxlen/DBLE(nx-1)
  dx = boxlen_x/DBLE(nx-1)
  dy = boxlen_y/DBLE(ny-1)

  !voldx = boxlen/DBLE(nx-1)
  DO i=1,nx-1
     nubox(1,i)=i
     nubox(2,i)=i+1
  ENDDO
  DO i=1,ny-1
     nuboy(1,i)=i
     nuboy(2,i)=i+1
  ENDDO


  dudt=0.
  one = 1.
  nfaces = nx+1

  DO i=1,nx-1
    do j = 1, ny-1
       i11=nubox(1,i)
       i12=nubox(2,i)
       i21=nuboy(1,j)
       i22=nuboy(2,j)
       
       uu11=u(:,i11,i21,1,1)
       uu12=u(:,i12,i21,1,1)
       uu21=u(:,i11,i22,1,1)
       uu22=u(:,i12,i22,1,1)
       
       !     dx=(x(i2)-x(i1))
       CALL compute_flux(uu11,ff11,g11,1,1,1,1)
       CALL compute_flux(uu12,ff12,g12,1,1,1,1)
       CALL compute_flux(uu21,ff21,g21,1,1,1,1)
       CALL compute_flux(uu22,ff22,g22,1,1,1,1)
      
       CALL compute_primitive(uu11,ww11,1,1,1,1)
       CALL compute_primitive(uu12,ww12,1,1,1,1)
       CALL compute_primitive(uu22,ww22,1,1,1,1)
       CALL compute_primitive(uu21,ww21,1,1,1,1)

       CALL get_source(ww11,source11,x(i11,i21,1,1),y(i11,i21,1,1))
       CALL get_source(ww12,source12,x(i12,i21,1,1),y(i12,i21,1,1))
       CALL get_source(ww21,source21,x(i11,i22,1,1),y(i11,i22,1,1))
       CALL get_source(ww22,source22,x(i12,i22,1,1),y(i12,i22,1,1))

       !totalphi=ff2-ff1-0.5*( source1+source2)*dx
       totalphi = -0.5*(ff11+ff21) &
                & +0.5*(ff12+ff22) &
                & -0.5*(g22+g21) &
                & +0.5*(g11+g12) - &
                & 0.25*(source11(:)+source21(:)+source12(:)+source22(:))*dx*dy 

       CALL compute_speed(ww11,a,b,c,speed1)
       CALL compute_speed(ww12,a,b,c,speed2)
       CALL compute_speed(ww21,a,b,c,speed3)
       CALL compute_speed(ww22,a,b,c,speed4)

       alpha=maxval([speed1,speed2, speed3,speed4])
       !write(*,*) alpha
       ubar=0.25*(uu11+uu22+uu21+uu12)
       phi(:,1)=dy*dx*diff(:,i11,i21,1,1)*0.25 + dt*(totalphi*0.25 + alpha*(uu11-ubar)*dy)
       phi(:,2)=dy*dx*diff(:,i12,i21,1,1)*0.25+dt*(totalphi*0.25 + alpha*(uu12-ubar)*dy)
       !alpha=max(speed3,speed4)

       !write(*,*) alpha
       phi(:,3)=dy*dx*diff(:,i11,i22,1,1)*0.25+dt*(totalphi*0.25+ alpha*(uu21-ubar)*dx)
       phi(:,4)=dy*dx*diff(:,i12,i22,1,1)*0.25+dt*(totalphi*0.25+ alpha*(uu22-ubar)*dx)

       !     CALL limit(phi1,phi2)
       if (use_limiter) then
         !CALL lim(uu1,uu2,phi,ph)
       else
         ph = phi
       end if
       ! call limit(phi(:,1),phi(:,2))

       dudt(:,i11,i21,1,1)=dudt(:,i11,i21,1,1)+ph(:,1)
       dudt(:,i12,i21,1,1)=dudt(:,i12,i21,1,1)+ph(:,2)
       dudt(:,i11,i22,1,1)=dudt(:,i11,i22,1,1)+ph(:,3)
       dudt(:,i12,i22,1,1)=dudt(:,i12,i22,1,1)+ph(:,4)

    end do 
   ENDDO

!!$  !
  !boundry conditions
  SELECT CASE(bc)
  CASE(2) ! periodic boundary
     dudt(:,1,:,1,1) = dudt(:,nx,:,1,1)
     dudt(:,nx,:,1,1)= dudt(:,1,:,1,1)
     dudt(:,:,ny,1,1) = dudt(:,:,1,1,1)
     dudt(:,:,1,1,1)= dudt(:,:,ny,1,1)

  CASE(3) ! transmissive
     dudt(:,1,:,1,1) = dudt(:,2,:,1,1)
     dudt(:,nx,:,1,1)= dudt(:,nx-1,:,1,1)
     dudt(:,:,ny,1,1) = dudt(:,:,ny-1,1,1)
     dudt(:,:,1,1,1)= dudt(:,:,2,1,1)
  CASE default
     PRINT*, "boundary condition type ", bc, " not taken into account"
     STOP
  END SELECT


END SUBROUTINE compute_update_rd

subroutine get_source(u,s,x,y)
    use parameters_rd
    implicit none

    real(kind=8),dimension(1:nvar)::u,s
    real(kind=8),dimension(1:nvar)::w
    real(kind=8),dimension(2)::grad_p
    real(kind=8)::x,y

    call compute_primitive(u,w,1,1,1,1)
    !call grad_phi(u,x,y,grad_p)
    grad_p(1) = 1.
    grad_p(2) = 1.
    
    s(1) = 0.
    s(2) = -w(1)*grad_p(1)
    s(3) = -w(1)*grad_p(2)
    s(4) = -w(1)*(w(2)*grad_p(1)&
             &+w(3)*grad_p(2))

end subroutine get_source

subroutine grad_phi(u,x,y, grad_p)
    use parameters_rd
    implicit none
    ! this could be a function by the way
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
    real(kind=8),dimension(1:nx,1:ny,1:mx,1:my,2)::grad_p
    real(kind=8),dimension(1:nx,1:ny,1:mx,1:my)::x,y

    real(kind=8)::delta_r,x_center,y_center
    real(kind=8)::x_dash,y_dash, r, epsilon

    integer::i,j,icell,jcell

    grad_p(:,:,:,:,1) = 1.
    grad_p(:,:,:,:,2) = 1.

  end subroutine

  SUBROUTINE lim(u1,u2,res_in,res_out)
  ! limitation
  REAL, PARAMETER:: gamma=1.4
  REAL(kind=8) ,DIMENSION(3,2), INTENT(in):: res_in
  REAL(kind=8) ,DIMENSION(3,2), INTENT(out):: res_out
  REAL(kind=8), DIMENSION(3), INTENT(in):: u1, u2
  REAL(kind=8)  :: phi 
  REAL(kind=8) ,DIMENSION(3):: r0, r1, r2, vp
  REAL(kind=8) :: ub, hb, ro1, ro2, v1, v2, h1, h2 , den, c2, c
  REAL(kind=8) :: p1, p2, p, c1, pp1, pp2,pp
  REAL(kind=8) :: kappa, x1 , x2, blend
  REAL(kind=8), DIMENSION(2):: aa, bb
  REAL(kind=8), DIMENSION(2,3):: limi, l
  INTEGER :: k
  kappa=gamma-1

  ro1=u1(1)
  ro2=u2(1)
  v1=u1(2)/u1(1)
  v2=u2(2)/u2(1)
  p1=kappa*( u1(3)-0.5*ro1*v1*v1)
  p2=kappa*( u2(3)-0.5*ro2*v2*v2)
  h1=(u1(3)+p1)/ro1
  h2=(u2(3)+p2)/ro2
  c1=SQRT(gamma*p1/ro1)
  c2=SQRT(gamma*p2/ro2)
! this test to have first order in space
!!$    IF (iordre.EQ.1) THEN
!!$       beta=res_in
!!$    ELSE
  !
  ! Compute Roe average
  !


  den=SQRT(ro1)/(SQRT(ro1)+SQRT(ro2))
  ub= den*v1+(1-den)*v2
  hb= den*h1+(1-den)*h2
  c2=kappa*( hb- 0.5*ub*ub )
  c=SQRT(c2)
  vp(1)=ub; vp(2)=ub-c; vp(3)=ub+c

  ! Compute the projection of residual on eigenvectors
  !
  aa(:)=kappa/c2*( res_in(3,:)-ub*res_in(2,:)+0.5*ub*ub*res_in(1,:))
  bb(:)=(res_in(2,:)-ub*res_in(1,:))/c
  l(:,1)=res_in(1,:)-aa(:)
  l(:,2)=0.5*( aa(:)-bb(:) )
  l(:,3)=0.5*( aa(:)+bb(:) )

  r0(1)=1.; r0(2)=ub  ; r0(3)=ub*ub*0.5
  r1(1)=1.; r1(2)=ub-c; r1(3)=hb-ub*c
  r2(1)=1.; r2(2)=ub+c; r2(3)=hb+ub*c
  !
  ! loop over the left eigenvectors for limiting
  !
  DO k=1,3
     pp1=l(1,k)
     pp2=l(2,k)
     pp=pp1 + pp2
     ! psi
     IF (ABS(pp)>0.) THEN 
        den=MAX(pp1/pp,0.)+MAX(pp2/pp,0.)
        limi(1,k)=MAX(pp1/pp,0.)/den*pp
        limi(2,k)=MAX(pp2/pp,0.)/den*pp

        limi(1,k)= limi(1,k)- 0.5*SIGN(1.d0,vp(k))*pp
        limi(2,k)= limi(2,k)+ 0.5*SIGN(1.d0,vp(k))*pp
     ELSE
        limi(1,k)=0.
        limi(1,k) = pp1
        limi(2,k) = pp2 
     ENDIF
!!$! blending between  lxF and SUPG
     blend= (abs(pp)/(abs(pp1)+abs(pp2)+1.e-10*dx))
     limi(1,k)= 0.5*(1.-blend)*(1.-sign(dble(1),vp(k)))*pp+ blend*l(1,k)
     limi(2,k)= 0.5*(1.-blend)*(1.+sign(dble(1),vp(k)))*pp+ blend*l(2,k)
!!$
     ! Here blending between  Galerkin and first order scheme
     !     limi(1,k)= 0.5*(1.-blend)*pp+ blend*l(1,k)
     !     limi(2,k)= 0.5*(1.-blend)*pp+ blend*l(2,k)
     ! spectail treatment of the u wave
!!!$if (k.ne.2) then
!!$!!!$          limi(1,k)= 0.5*(1.-blend)*(1.-sign(1.,vp(k)))*pp+ blend*l(1,k)
!!$!!!$          limi(2,k)= 0.5*(1.-blend)*(1.+sign(1.,vp(k)))*pp+ blend*l(2,k)
!!!$endif
        !  limi(1,k)= 0.5*(1.-blend)*pp+ blend*l(1,k)
       !   limi(2,k)= 0.5*(1.-blend)*pp+ blend*l(2,k)
!!$! blend Galerkin and LxF
      !limi(:,k)= (1-blend)*limi(:,k)+ blend*l(:,k)
  ENDDO

  !limi(1,:)=l(1,:)
  !limi(2,:)=l(2,:)
  res_out(:,1)=limi(1,1)*r0+limi(1,2)*r1+limi(1,3)*r2
  res_out(:,2)=limi(2,1)*r0+limi(2,2)*r1+limi(2,3)*r2
!res_out=res_in
  !res_out(:,1)=limi(1,1)+limi(1,2)+limi(1,3)
  !res_out(:,2)=limi(2,1)+limi(2,2)+limi(2,3)
  !    ENDIF
  !print*, "sum(L*R", limi(1,1)*r0+limi(1,2)*r1+limi(1,3)*r2
  !print*, "res_in", res_in(:,1)
  !print*, "sum(L*R)", limi(2,1)*r0+limi(2,2)*r1+limi(2,3)*r2
  !print*, "res_in", res_in(:,2)
  ! special treatment to avoid sonic glinch
  !
!!$    IF (v1+c1>v2+c2.OR.v1-c1<v2-c2) THEN
!!$       res_out=res_in
!!$    ENDIF
END SUBROUTINE lim
