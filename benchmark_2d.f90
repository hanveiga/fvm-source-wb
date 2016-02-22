  !----
  ! global variables: parameters

  program main
    use parameters_2d
    real(kind=8),dimension(1:nx,1:ny)::x
    real(kind=8),dimension(1:nx,1:ny)::y
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::u, w, u_eq

    call get_coords(x,y,nx,ny)

    call get_initial_conditions(x,y,u,nx,ny)

    call output_file(x, y ,u,'initwo')

    call get_equilibrium_solution(x, y, u_eq, nx, ny)

    call evolve(u, u_eq)

    call output_file(x, y ,u,'fintwo')

  end program main

  !---------
  subroutine get_coords(x,y,size_x,size_y)
    use parameters_2d
    integer::size_x,size_y
    real(kind=8),dimension(1:size_x,1:size_y)::x,y

    !internal vars
    real(kind=8)::dx, dy
    integer::i,j

    dx = boxlen_x/dble(size_x)
    dy = boxlen_y/dble(size_y)
    do i=1,size_x
      do j=1,size_y
        x(i,j) = (i-0.5)*dx
        y(i,j) = (j-0.5)*dy
      end do
    end do

  end subroutine get_coords
  !-------
  subroutine get_initial_conditions(x,y,u,size_x,size_y)
    use parameters_2d
    integer::size_x,size_y, i ,j 
    real(kind=8),dimension(1:nvar,1:size_x, 1:size_y)::u
    real(kind=8),dimension(1:size_x, size_y)::x
    real(kind=8),dimension(1:size_x,size_y)::y

    ! internal variables
    real(kind=8),dimension(1:nvar,1:size_x, 1:size_y)::w
    real(kind=8)::rho_0, p_0, g

    select case (ninit)
      case(1)
        w(1,:,:) = exp(-(x+y))
        w(2,:,:) = 0
        w(3,:,:) = 0
        w(4,:,:) = exp(-(x+y)) !+ eta*&
                 !&exp(-100*((x-0.3)**2+(y-0.3)**2))
      case(2)
        rho_0 = 1.21
        p_0 = 1.
        g = 1.
        w(1,:,:) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
        w(2,:,:) = 0
        w(3,:,:) = 0
        w(4,:,:) = p_0*exp(-(rho_0*g/p_0)*(x+y))
      case(3)
        rho_0 = 1.21
        p_0 = 1
        g = 1
        w(1,:,:) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
        w(2,:,:) = 0
        w(3,:,:) = 0
        w(4,:,:) = p_0*exp(-(rho_0*g/p_0)*(x+y))+eta*&
                 &exp(-(rho_0*g/p_0)*((x-0.3)**2+(y-0.3)**2))

      case(4)
        do i = 1,size_x
          do j = 1,size_y

        if (x(i,j)>=0.5 .and. y(i,j)>=0.5) then        
          w(1,i,j) = 1.5
          w(2,i,j)  = 0.
          w(3,i,j)  = 0.
          w(4,i,j)  = 1.5
        else if (x(i,j)<0.5 .and. y(i,j)>=0.5) then
          w(1,i,j) = 0.5323
          w(2,i,j) = 1.206
          w(3,i,j) = 0.
          w(4,i,j) = 0.3
        else if (x(i,j)<0.5 .and. y(i,j)<0.5) then
          w(1,i,j) = 0.138
          w(2,i,j) = 1.206
          w(3,i,j) = 1.206
          w(4,i,j) = 0.029
        else if (x(i,j)>=0.5 .and. y(i,j)<0.5) then
        !else
          w(1,i,j) = 0.5323
          w(2,i,j) = 0.
          w(3,i,j) = 1.206
          w(4,i,j) = 0.3
        end if
      end do 
    end do 
    end select

    call compute_conservative(w,u,size_x,size_y)

  end subroutine get_initial_conditions

  subroutine output_file(x_values, y_values, function_values, filen)
    use parameters_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny)::function_values
    real(kind=8),dimension(1:nx,1:ny)::x_values, y_values
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
       xcell=x_values(icell,jcell)
       ycell = y_values(icell,jcell)
       call get_equilibrium_solution([xcell],[ycell],w_eq,one,one)
       call compute_primitive(function_values(1:nvar,icell,jcell),w,one,one)
       write(10,'(7(1PE12.5,1X))')xcell,ycell,(w,ivar=4,nvar)
      end do
    end do
    close(10)

  end subroutine output_file
  !--------
  subroutine compute_primitive(u,w,size_x,size_y)
    use parameters_2d
    implicit none
    integer::size_x,size_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::w
    ! Compute primitive variables
    !write(*,*) u(1,:,:)
    w(1,:,:) = u(1,:,:)
    w(2,:,:) = u(2,:,:)/w(1,:,:)
    w(3,:,:) = u(3,:,:)/w(1,:,:)
    w(4,:,:) = (gamma-1.0)*( u(4,:,:) - 0.5*w(1,:,:)*(w(2,:,:)**2+w(3,:,:)**2) )
  end subroutine compute_primitive
  !-------------
  subroutine compute_conservative(ww,u,size_x,size_y)
    use parameters_2d
    implicit none
    integer::size_x,size_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::ww
    ! Compute primitive variables
    u(1,:,:)= ww(1,:,:)
    u(2,:,:)= ww(1,:,:)*ww(2,:,:)
    u(3,:,:)= ww(1,:,:)*ww(3,:,:)
    u(4,:,:)= ww(4,:,:)/(gamma-1.) + 0.5*( ww(1,:,:)*(ww(2,:,:)**2+ww(3,:,:)**2) )
    !ww(3,:,:)/(gamma-1.0)+0.5*ww(1,:,:)*ww(2,:,:)**2
  end subroutine compute_conservative
  !-------------
  subroutine get_equilibrium_solution(x,y,w, size_x,size_y)
    ! get equilibrium solution for primitive variables
    ! use parameters
    ! communicate which initial solution, equilibrium type
    ! boundary conditions and domain properties
    use parameters_2d
    integer::size_x,size_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::w
    real(kind=8),dimension(1:size_x,1:size_y)::x
    real(kind=8),dimension(1:size_x,1:size_y)::y
    
    !class(:)::x
    ! internal variables
    real(kind=8)::rho_0,p_0,g

    select case (nequilibrium)
      case(1)
        w(1,:,:) = exp(-(x+y))
        w(2,:,:) = 0
        w(3,:,:) = 0
        w(4,:,:) = exp(-(x+y))
      case(2)
        rho_0 = 1.21
        p_0 = 1
        g = 1
        w(1,:,:) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
        w(2,:,:) = 0
        w(3,:,:) = 0
        w(4,:,:) = p_0*exp(-(rho_0*g/p_0)*(x+y))
      case(3)
        rho_0 = 1.21
        p_0 = 1
        g = 1
        w(1,:,:) = rho_0*exp(-(rho_0*g/p_0)*(x+y))
        w(2,:,:) = 0
        w(3,:,:) = 0
        w(4,:,:) = p_0*exp(-(rho_0*g/p_0)*(x+y))
      case(4)
        w(1,:,:) = 0.
        w(2,:,:) = 0.
        w(3,:,:) = 0.
        w(4,:,:) = 0.
    end select

  end subroutine get_equilibrium_solution
  !---------

  subroutine evolve(u, u_eq)
    use parameters_2d
    implicit none

    real(kind=8),dimension(1:nvar,1:nx,1:ny)::u,u_new,u_eq

    ! internal variables
    real(kind=8)::t,dt
    real(kind=8)::cmax, dx, dy
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::dudt, w, w1
    integer::iter, n, i, j

    dx = boxlen_x/dble(nx)
    dy = boxlen_y/dble(ny)

    t=0
    iter=0
    do while(t < tend)
       ! Compute time step
       call compute_max_speed(u,cmax)
       dt=0.8*dx/cmax*cfl
       
       !if(solver=='EQL')then
         ! runge kutta 2nd order
        call compute_update(u,u_eq, dudt)
        w1=u+dt*dudt
        !call dirichlet(w1)
        call compute_update(w1,u_eq,dudt)
        u=0.5*u+0.5*w1+0.5*dt*dudt
        !fcall dirichlet(u)
       !endif

       t=t+dt
       iter=iter+1
       write(*,*)'time=',iter,t,dt, cmax
    end do

    u_new = u
  end subroutine evolve

  !-----

  subroutine compute_max_speed(u,cmax)
    use parameters_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny)::u
    real(kind=8)::cmax
    integer::icell, jcell
    real(kind=8)::speed
    ! Compute max sound speed
    cmax=0.0
    do icell=1,nx
      do jcell = 1,ny
       call compute_speed(u(1:nvar,icell,jcell),speed)
       cmax=MAX(cmax,speed)
      end do
    end do 
  end subroutine compute_max_speed

  !--------------------

  subroutine compute_speed(u,speed)
    use parameters_2d
    implicit none
    real(kind=8),dimension(1:nvar)::u
    real(kind=8)::speed
    real(kind=8),dimension(1:nvar)::w
    real(kind=8)::cs
    ! Compute primitive variables
    call compute_primitive(u,w,1,1)
    ! Compute sound speed
    cs=sqrt(gamma*max(w(4),1d-10)/max(w(1),1d-10))
    speed=sqrt(w(2)**2+w(3)**2)+cs
  end subroutine compute_speed

  !--------------

  subroutine compute_flux(u,flux, size_x,size_y)
    use parameters_2d
    integer::size_x, size_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::u
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y,2)::flux
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::w
    ! Compute primitive variables
    call compute_primitive(u,w,size_x,size_y)
    ! Compute flux

    !write(*,*) 'u'
    !write(*,*) u(:,:,:)

    !write(*,*) 'primitive'
    !write(*,*) w(:,:,:)

    flux(1,:,:,1)=w(2,:,:)*u(1,:,:)
    flux(2,:,:,1)=w(2,:,:)*u(2,:,:)+w(4,:,:)
    flux(3,:,:,1)=w(1,:,:)*w(2,:,:)*w(3,:,:)
    flux(4,:,:,1)=w(2,:,:)*u(4,:,:)+w(2,:,:)*w(4,:,:)

    flux(1,:,:,2) = u(1,:,:)*w(3,:,:)
    flux(2,:,:,2) = u(2,:,:)*w(3,:,:)
    flux(3,:,:,2) = u(3,:,:)*w(3,:,:)+w(4,:,:)
    flux(4,:,:,2) = w(3,:,:)*u(4,:,:)+w(3,:,:)*w(4,:,:)

  end subroutine compute_flux

  subroutine get_source(w,s,x,y,size_x,size_y)
    use parameters_2d
    implicit none
    integer::size_x, size_y
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::u,s
    real(kind=8),dimension(1:nvar,1:size_x,1:size_y)::w
    real(kind=8),dimension(1:size_x,1:size_y)::x
    real(kind=8),dimension(1:size_x,1:size_y)::y

    !internal

    !x_minus = [x(1),x(1:size-1)] ! zero gradient
    !x_plus = [x(2:size),x(size)] ! zero gradient

    s(1,:,:) = 0.
    s(2,:,:) = -w(1,:,:)*1
    s(3,:,:) = 0.
    s(4,:,:) = -w(1,:,:)*sqrt(w(2,:,:)**2+w(3,:,:)**2)

    !s(1) = 0
    !  s(2) = -w(1)*1*(x_plus-x_minus)/(2*delta)
    !  s(3) = -w(1)*w(2)*1*(x_plus-x_minus)/(2*delta)

  end subroutine get_source
  !--------

  subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv)
    use parameters_2d
    implicit none
    real(kind=8),dimension(1:nvar)::uleft,uright, f_left, f_right
    real(kind=8),dimension(1:nvar)::fgdnv
    real(kind=8)::cleft,cright,cmax
    real(kind=8),dimension(1:nvar)::fleft,fright
    ! Maximum wave speed
    call compute_speed(uleft,cleft)
    call compute_speed(uright,cright)
    cmax=max(cleft,cright)
    !write(*,*) cmax
    ! Compute Godunox flux
    fgdnv=0.5*(f_right+f_left)+0.5*cmax*(uleft-uright)
  end subroutine compute_llflux
  !-----

  subroutine compute_update(u, w_eq, dudt)
    use parameters_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::u, w_eq, dudt
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::u_top, u_bottom 
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::u_left, u_right
    real(kind=8),dimension(1:nvar,1:nx, 1:ny,2)::flux_top, flux_bottom  
    real(kind=8),dimension(1:nvar,1:nx, 1:ny,2)::flux_right, flux_left
    real(kind=8),dimension(1:nvar,1:nx, 1:ny,2)::flux
    real(kind=8),dimension(1:nvar,1:nx,1:(ny+1))::G
    real(kind=8),dimension(1:nvar,1:(nx+1),1:ny)::F
    real(kind=8)::dx,dy, oneoverdy, oneoverdx
    integer::j,ileft,iright,itop,ibottom,i
    integer::iface, jface

    ! internal vars
    dx = boxlen_x/dble(nx)
    dy = boxlen_y/dble(ny)

    oneoverdx = 1/dx
    oneoverdy = 1/dy

    F(:,:,:) = 0.
    G(:,:,:) = 0.

    ! make u_left, u_right, u_top, u_bottom
    u_left = u
    u_right = u
    u_top = u
    u_bottom = u

    call compute_flux(u_left,flux_left,nx,ny)
    call compute_flux(u_right,flux_right,nx,ny)
    call compute_flux(u_top,flux_top,nx,ny)
    call compute_flux(u_bottom,flux_bottom,nx,ny)
    !write(*,*) 'flux left'
    !write(*,*) flux_left(1,:,:,1)

    ! compute artificial flux in x direction
    
    do j = 1, ny
    do iface = 1, nx+1
      ileft = iface-1
      iright = iface
      if (iface == 1) then
        ileft = nx
      end if
      if (iface == nx+1) then
        iright = 1
      end if

        ! subroutine compute_llflux(uleft,uright, f_left,f_right, fgdnv)
      call compute_llflux(u_right(1:nvar,ileft,j),u_left(1:nvar,iright,j),&
                            &flux_right(1:nvar,ileft,j,1),flux_left(1:nvar,iright,j,1),F(1:nvar,iface,j))
    end do
    end do 


    do i = 1,nx
      do jface = 1,ny+1
        ileft = jface-1
        iright = jface
        if (jface == 1) then
          ileft = ny
        end if
        if (jface == ny+1) then
          iright = 1
        end if

       call compute_llflux(u_top(1:nvar,i,ileft),u_bottom(1:nvar,i,iright),&
                            &flux_top(1:nvar,i,ileft,2),flux_bottom(1:nvar,i,iright,2),G(1:nvar,i,jface))
      end do
    end do

    do i = 1 ,nx
      do j = 1 ,ny 
        dudt(1:nvar,i, j) =  -(F(1:nvar,i+1,j)-F(1:nvar,i,j))*oneoverdx &
          &-(G(1:nvar,i,j+1)-G(1:nvar,i,j))*oneoverdy !+ s(1:nvar,1:nx, 1:ny)
    end do
    end do
    write(*,*) dudt(1,2,2)
  !+ s(1:nvar,1:nx) + &
  !        & (flux_eq(1:nvar,2:(nx+1),1:ny)-flux_eq(1:nvar,1:nx,1:ny))*oneoverdx - s_eq(1:nvar,1:nx) &

  end subroutine compute_update

  subroutine compute_update_old(u, w_eq, dudt)
    use parameters_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::u, w_eq, dudt

    ! internal vars
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::delta_w, f_left, f_right
    real(kind=8),dimension(1:nvar,1:nx, 1:ny)::g_top, g_bottom
    real(kind=8),dimension(1:nvar)::flux_r, flux_g
    real(kind=8),dimension(1:(nx+1),1:(ny+1))::x_faces, y_faces
    real(kind=8),dimension(1:nx,1:ny)::x, y
    real(kind=8),dimension(1:nvar,1:(nx+1),1:ny)::u_eq_faces_x, w_eq_faces_x
    real(kind=8),dimension(1:nvar,1:(nx+1),1:ny)::flux_eq_x, flux_riemann
    real(kind=8),dimension(1:nvar,1:nx,1:(ny+1))::u_eq_faces_y, w_eq_faces_y, flux_eq_y,g_flux

    real(kind=8),dimension(1:nvar, 1:nx, 1:ny)::u_left,u_right, u_top, u_bottom, u_eq
    real(kind=8),dimension(1:nvar,1:nx,1:ny)::s,s_eq
    real(kind=8)::oneoverdx, dx, dt, dx_plus, zero, xface, xcell, oneoverdy, dy
    integer::one, i, nfaces, iface, ileft, iright, jface, itop, ibottom, jcell, icell,j
    
    dx = boxlen_x/dble(nx)
    dy = boxlen_y/dble(ny)
    oneoverdx = 1/dx
    oneoverdy = 1/dy
    nfaces = nx + 1

    zero = 0

    ! compute perturbation
    call compute_conservative(w_eq, u_eq, nx, ny)

    delta_w = u - u_eq

    ! propagate perturbation to faces
    ! create x_faces
    ! create x
    do j = 1,ny
      do i = 1,nx
        x(i,j) = (i-0.5)*dx
        y(j,j) = (j-0.5)*dy
      end do
    end do

    do j = 1,ny
      do i = 1,nfaces
      x_faces(i,j) = (i-1)*dx
      end do
    end do
    do j = 1,nfaces
      do i = 1,nx
      y_faces(i,j) = (j-1)*dy
      end do
    end do

    ! x - direction
    call get_equilibrium_solution(x_faces,y,w_eq_faces_x,nfaces,ny)
    call compute_conservative(w_eq_faces_x, u_eq_faces_x, nfaces, ny)

    !u_left(1:nvar,1:nx,1:ny) = delta_w(1:nvar,1:nx,1:ny) + u_eq_faces_x(1:nvar,1:nx, 1:ny)
    !u_right(1:nvar,1:nx,1:ny) = delta_w(1:nvar,1:nx,1:ny) + u_eq_faces_x(1:nvar,2:(nx+1), 1:ny)
    write(*,*) 'u'
    write(*,*) u(1,:,:)
    do j = 1, ny
      do i = 1 , nx
        if (i == 1) then
          u_left(1:nvar,i,1:ny) = u(1:nvar,i,1:ny)
          u_right(1:nvar,i,1:ny) = u(1:nvar,i+1,1:ny)
          cycle
        end if
         ! u_right(1:nvar,i,1:ny) = u(1:nvar,i,1:ny)
        if (i==nx) then
           u_right(1:nvar,i,1:ny) = u(1:nvar,i,1:ny)
           u_left(1:nvar,i,1:ny) = u(1:nvar, i, 1:ny)
           cycle 
        end if
        u_left(1:nvar,i,1:ny) = u(1:nvar,i-1,1:ny)
        u_right(1:nvar,i,1:ny) = u(1:nvar,i,1:ny)
      end do 
    end do
    write(*,*) 'uleft pressure'
    write(*,*) u_left(1,:,:)

    !do j = 1, ny
      do i = 1 , nx
        if (i == 1) then
          u_bottom(1:nvar,1:nx,i) = u(1:nvar,1:nx,i)
          u_top(1:nvar,1:nx,i) = u(1:nvar,1:nx,i+1)
          cycle
        end if
         ! u_right(1:nvar,i,1:ny) = u(1:nvar,i,1:ny)
        if (i==nx) then
           u_bottom(1:nvar,1:nx,i) = u(1:nvar,1:ny,i-1)
           u_top(1:nvar,1:nx,i) = u(1:nvar, 1:ny,i)
           cycle 
        end if
        u_bottom(1:nvar,1:nx,i) = u(1:nvar,1:nx,i-1)
        u_top(1:nvar,1:nx,i) = u(1:nvar,1:nx,i)
      end do 
    !end do
    write(*,*) 'uleft pressure'
    write(*,*) u_left(1,:,:)


    ! y - direction
    call get_equilibrium_solution(x,y_faces,w_eq_faces_y,nx,nfaces)
    call compute_conservative(w_eq_faces_y, u_eq_faces_y, nx, nfaces)
    !write(*,*) 'u eq faces'
    !write(*,*) u_eq_faces_y(1,:,:)
    !u_top(1:nvar,1:nx,1:ny) = u_eq_faces_y(1:nvar,1:nx,2:(ny+1))!delta_w(1:nvar,1:nx,1:ny) + u_eq_faces_y(1:nvar,1:nx, 2:(ny+1))
    !u_bottom(1:nvar,1:nx,1:ny) = u_eq_faces_y(1:nvar,1:nx,1:ny)!delta_w(1:nvar,1:nx,1:ny) + u_eq_faces_y(1:nvar,1:nx, 1:ny)
    ! compute fluxes
    write(*,*) 'ubottom'
    write(*,*) u_bottom(1,:,:)
    !call compute_flux_x(u_left, f_left, nx,ny)
    !call compute_flux_x(u_right, f_right, nx,ny)
    !call compute_flux_x(u_eq_faces_x,flux_eq_x,nfaces,ny)

    !call compute_flux_y(u_top, g_top, nx,ny)
    !call compute_flux_y(u_bottom, g_bottom, nx,ny)
    !call compute_flux_y(u_eq_faces_y,flux_eq_y,nx,nfaces)

    one = 1
    ! compute numerical flux
    do jcell = 1,ny
      do iface = 1,nfaces
        ileft = iface - 1
        iright = iface
        if (bc==2) then ! zero gradient
          if(iface==1) ileft = 1
          if(iface==nx+1) iright = nx
        end if

        !write(*,*) 'u right, u left'
        !write(*,*) u_right(1:nvar,ileft,jcell),u_left(1:nvar,iright,jcell)

        call compute_llflux(u_right(1:nvar,ileft,jcell),u_left(1:nvar,iright,jcell),&
                        & f_right(1:nvar,ileft,jcell), f_left(1:nvar,iright,jcell), flux_r)

        !write(*,*) 'f right,f left'
        !write(*,*) f_right(1:nvar,ileft,jcell), f_left(1:nvar,iright,jcell)
        !write(*,*) 'riemann flux'
        !write(*,*) flux_r

        flux_riemann(1:nvar,iface,jcell) = flux_r
      end do
    end do

    do icell = 1,nx
      do jface = 1,nfaces
        itop = jface
        ibottom = jface - 1
        if (bc==2) then ! zero gradient
          if(jface==1) ibottom = 1
          if(jface==ny+1) itop = nx
        end if

        call compute_llflux(u_top(1:nvar,icell,ibottom),u_bottom(1:nvar,icell,itop),&
                        & g_top(1:nvar,icell, ibottom), g_bottom(1:nvar,icell, itop), flux_g)

        g_flux(1:nvar,icell,jface) = flux_g
      end do
    end do
    write(*,*) 'flux riemann'
    write (*,*) flux_riemann(1,:,:)
    write(*,*) 'flux g riemann'
    write (*,*) g_flux(1,:,:)
    ! compute source
    call get_source(u,s,x,y, nx,ny)
    !call get_source(u_eq,s_eq, x, nx)

    dudt(1:nvar,1:nx, 1:ny) = -(flux_riemann(1:nvar,2:(nx+1),1:ny)-flux_riemann(1:nvar,1:nx,1:ny))*oneoverdx &
          &-(g_flux(1:nvar,1:nx,2:(ny+1))-g_flux(1:nvar,1:nx,1:ny))*oneoverdy !+ s(1:nvar,1:nx, 1:ny)
    write(*,*) 'dudt'
    write(*,*) dudt(1,:,:)
  !+ s(1:nvar,1:nx) + &
  !        & (flux_eq(1:nvar,2:(nx+1),1:ny)-flux_eq(1:nvar,1:nx,1:ny))*oneoverdx - s_eq(1:nvar,1:nx) &

  end subroutine compute_update_old
