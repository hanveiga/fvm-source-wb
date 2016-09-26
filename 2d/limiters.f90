!function minmod(x,y,z)
!implicit none
!real(kind=8)::x,y,z,s,minmod,dlim,slop
!s=sign(1d0,x)
! slop = min(abs(y),abs(z))
! dlim = slop
!  !if((y*z)<=0.)dlim=0.
!  !minmod = s*min(dlim,abs(x))
! if(sign(1d0,y)==s.AND.sign(1d0,z)==s)then
!    minmod=s*min(abs(x),abs(y),abs(z))
! else
!   minmod=0.0
! endif
! return
!end function minmod


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

function generalized_minmod(x,y,z)
  use parameters_dg_2d
  implicit none
  ! input
  real(kind=8)::x,y,z

  ! output
  real(kind=8)::minmod, dx
  real(kind=8)::generalized_minmod

  ! internal
  real(kind=8)::s

  s = sign(1d0,x)
  dx = boxlen_x/dble(nx)

  if (abs(x) < M*dx**2) then
    generalized_minmod = x
  else
    generalized_minmod = minmod(x,y,z)
  endif

  !write(*,*) minmod
  return
end function generalized_minmod


function minmod2d(u,d_l_x,d_l_y,d_r_x,d_r_y)
  implicit none
  ! input
  real(kind=8)::u,d_l_x,d_l_y,d_r_x,d_r_y

  ! output
  real(kind=8)::minmod2d

  ! internal
  real(kind=8)::s

  s = sign(1d0,u)

  if (sign(1d0,d_l_x) == s .AND. sign(1d0,d_l_y) == s &
  & .and. sign(1d0,d_r_x) == s .and. sign(1d0,d_r_y) == s)  then
  minmod2d = s*min(abs(u),abs(d_l_y),abs(d_l_x),abs(d_r_y),abs(d_r_x))
else
  minmod2d = 0.0
end if
!write(*,*) minmod
return
end function minmod2d



subroutine high_order_limiter(u)
  use parameters_dg_2d
  implicit none

  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, u_new
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes
  real(kind=8)::limited,minmod2d, d_l_x, d_l_y, d_r_x, d_r_y
  real(kind=8)::coeff_plus,coeff_minus,coeff_c,central_u
  integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
  integer::intnode, jntnode

  if (mx==1.and.my==1) then
    ! no limiting when we only have 1st order approx
    !use_limiter = .false.
    return
  end if

  call compute_primitive(u,w,nx,ny,mx,my)
  ! TODO: implement characteristic variable representation
  ! using Roe average
  write(*,*) u(1,1,1,:,:)
  do ivar = 1,nvar
    do icell=1,nx
      do jcell = 1,ny
        ileft = icell - 1
        iright = icell + 1
        itop = jcell + 1
        ibottom = jcell - 1

        call get_boundary_conditions(ileft,1)
        call get_boundary_conditions(iright,1)
        call get_boundary_conditions(itop,2)
        call get_boundary_conditions(ibottom,2)

        do intnode = mx-1,1,-1
          do jntnode = intnode,1,-1
            coeff_c = sqrt(2.0*dble(intnode)+1.0)*sqrt(2.0*dble(jntnode)+1.0)/dble(2.0)
            coeff_minus = sqrt(2.0*dble(intnode)+1.0)*sqrt(2.0*dble(jntnode)+1.0)/dble(2.0)
            !coeff_plus = sqrt(2.0*dble(my-jntnode)+1.0)*sqrt(2.0*dble(mx-intnode)+1.0)/dble(2.0)

            central_u = u(ivar,icell,jcell,intnode+1,jntnode+1)

            d_l_x = u(ivar,icell,jcell,mx-intnode,my-intnode) - &
            & u(ivar,ileft,jcell,mx-intnode,my-intnode)

            d_l_y = u(ivar,icell,jcell,mx-intnode,my-intnode) - &
            & u(ivar,icell,ibottom,mx-intnode,my-intnode)

            d_r_x = u(ivar,iright,jcell,mx-intnode,my-intnode) - &
            & u(ivar,icell,jcell,mx-intnode,my-intnode)

            d_r_y = u(ivar,icell,itop,mx-intnode,my-intnode) - &
            & u(ivar,icell,jcell,mx-intnode,my-intnode)
            write(*,*) 'intnodes', intnode, jntnode
            write(*,*) 'central u, coeff', central_u, coeff_minus*d_l_x

            limited = minmod2d(central_u,coeff_minus*d_l_x,&
            &coeff_minus*d_l_y, coeff_minus*d_r_y, coeff_minus*d_r_x)

            !limited2 = minmod2d(central_u,coeff_minus*d_l_x,&
            !  &coeff_minus*d_l_y, coeff_minus*d_r_y, coeff_minus*d_r_x)


            write (*,*) 'limited', limited
            if (intnode == jntnode) then
              if (abs(limited -   central_u) < 1e-10) then
                exit
              end if
            end if

            u_new(ivar,icell,jcell,mx-intnode+1,my-intnode+1) = limited

          end do
          write (*,*) 'limited', limited
          if (intnode == jntnode) then
            if (abs(limited -   central_u) < 1e-10) then
              exit
            end if
          end if
        end do
        write (*,*) 'limited', limited
        if (intnode == jntnode) then
          if (abs(limited -   central_u) < 1e-10) then
            exit
          end if
        end if
      end do
    end do
  end do

  ! guarantee positivity of density and pressure
  call get_nodes_from_modes(u_new,nodes,nx,ny,mx,my)
  call compute_primitive(nodes,w,nx,ny,mx,my)

  do icell=1,nx
    do jcell = 1,ny
      do intnode = 1,mx
        do jntnode = 1,my

          if (w(1,icell,jcell,intnode,jntnode)<1e-10) then
            w(1,icell,jcell,intnode,jntnode) = 0.0
          end if

          if (w(4,icell,jcell,intnode,jntnode)<1e-10) then
            w(4,icell,jcell,intnode,jntnode) = 0.0
          end if

        end do
      end do
    end do
  end do

  call compute_conservative(w,nodes,nx,ny,mx,my)
  call get_modes_from_nodes(nodes, u_new, nx, ny, mx, my)
  ! Update variables with limited states
  u = u_new

end subroutine high_order_limiter


subroutine compute_limiter(u)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes
  real(kind=8)::limited1,limited2, generalized_minmod
  integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
  integer::intnode, jntnode

  ! look at 1st derivatives, u12, u21
  if (nvar == 1) then
    return
  end if
  if (mx==1.and.my==1) then
    ! no limiting when we only have 1st order approx
    !use_limiter = .false.
    return
  end if

  ! mean is given by u11
  if (use_limiter) then
    do ivar = 1,nvar
      do icell=1,nx
        do jcell = 1,ny
          ileft = icell - 1
          iright = icell + 1
          itop = jcell + 1
          ibottom = jcell - 1

          call get_boundary_conditions(ileft,1)
          call get_boundary_conditions(iright,1)
          call get_boundary_conditions(itop,2)
          call get_boundary_conditions(ibottom,2)


          limited1 = generalized_minmod(u(ivar,icell,jcell,1,2), &
          &u(ivar,iright,jcell,1,1)-u(ivar,icell,jcell,1,1),&
          &u(ivar,icell,jcell,1,1)-u(ivar,ileft,jcell,1,1))

          limited2 = generalized_minmod(u(ivar,icell,jcell,2,1),&
          &u(ivar,icell,itop,1,1)-u(ivar,icell,jcell,1,1),&
          &u(ivar,icell,jcell,1,1)-u(ivar,icell,ibottom,1,1))

          if (abs(limited1 - u(ivar,icell,jcell,1,2))<1e-5) then
            u(ivar,icell,jcell,1,2) = limited1
          else
            u(ivar,icell,jcell,1,2:my) = 0.0
          end if

          if (abs(limited2 - u(ivar,icell,jcell,2,1))<1e-5) then
            u(ivar,icell,jcell,2,1) = limited2
          else
            u(ivar,icell,jcell,2:mx,1) = 0.0
          end if

        end do
      end do
    end do
  end if

  ! guarantee positivity of density and pressure
  call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
  call compute_primitive(nodes,w,nx,ny,mx,my)

  do icell=1,nx
    do jcell = 1,ny
      do intnode = 1,mx
        do jntnode = 1,my
          if (w(1,icell,jcell,intnode,jntnode)<1e-10) then
            w(1,icell,jcell,intnode,jntnode) = 1e-5
            !write(*,*) 'lim'
          end if

          if (w(4,icell,jcell,intnode,jntnode)<1e-10) then
            w(4,icell,jcell,intnode,jntnode) = 1e-5
          end if

        end do
      end do
    end do
  end do

  call compute_conservative(w,nodes,nx,ny,mx,my)
  call get_modes_from_nodes(nodes, u, nx, ny, mx, my)
  ! Update variables with limited states

end subroutine compute_limiter


subroutine limiter_low_order(u)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes
  real(kind=8)::limited1,limited2, generalized_minmod, dx
  integer::icell,jcell, ivar, itop,ibottom,ileft,iright
  integer::intnode, jntnode
  ! convert to nodes
  ! look at 1st derivatives, u12, u21
  ! if crappy then limit
  if (nvar == 1) then
    return
  end if

  if (mx==1.and.my==1) then
    ! no limiting when we only have 1st order approx
    !use_limiter = .false.
    return
  end if

  dx = boxlen_x/dble(nx)
  call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
  !write(*,*) nodes(2,4,4,:,:)
  ! mean is given by u11
  if (use_limiter) then
    do ivar = 1,nvar
      do icell=1,nx
        do jcell = 1,ny
          ileft = icell - 1
          iright = icell + 1
          itop = jcell + 1
          ibottom = jcell - 1

          call get_boundary_conditions(ileft,1)
          call get_boundary_conditions(iright,1)
          call get_boundary_conditions(itop,2)
          call get_boundary_conditions(ibottom,2)

          do intnode = 1,mx
            do jntnode = 1,my
              !    write(*,*) 'lim'
              !limited1 = generalized_minmod(nodes(ivar,icell,jcell,jntnode,intnode), &
              !&(nodes(ivar,iright,jcell,jntnode,jntnode)-nodes(ivar,icell,jcell,jntnode,jntnode))/dx,&
              !&(nodes(ivar,icell,jcell,jntnode,jntnode)-nodes(ivar,ileft,jcell,jntnode,jntnode))/dx)

              !limited2 = generalized_minmod(nodes(ivar,icell,jcell,intnode,jntnode),&
              !&u(ivar,icell,itop,1,1)-u(ivar,icell,jcell,1,1),&
              !&u(ivar,icell,jcell,1,1)-u(ivar,icell,ibottom,1,1))

              !if (abs(limited1 - u(ivar,icell,jcell,jntnode,intnode))<1e-5&
              !  &.and.abs(limited2 - u(ivar,icell,jcell,2,1))<1e-5) then
              !  u(ivar,icell,jcell,jntnode,intnode) = limited1
              !  u(ivar,icell,jcell,intnode,jntnode) = limited2
              !else
              u(ivar,icell,jcell,1,2:my) = 0.0
              u(ivar,icell,jcell,2:mx,1) = 0.0
              !end if
            end do
          end do
        end do
      end do
    end do
  end if
  call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
  !write(*,*) nodes(2,4,4,:,:)
  !call get_nodes_from_modes(u,inodes,nx,ny,mx,my)
  !call compute_primitive(nodes,w,nx,ny,mx,my)

  !do icell=1,nx
  !  do jcell = 1,ny
  !    do intnode = 1,mx
  !      do jntnode = 1,my
  !        if (w(1,icell,jcell,intnode,jntnode)<1e-10) then
  !          w(1,icell,jcell,intnode,jntnode) = 0.0
  !        end if

  !        if (w(4,icell,jcell,intnode,jntnode)<1e-10) then
  !          w(4,icell,jcell,intnode,jntnode) = 0.0
  !        end if

  !      end do
  !    end do
  !  end do
  !end do

  !call compute_conservative(w,nodes,nx,ny,mx,my)
  !call get_modes_from_nodes(nodes, u, nx, ny, mx, my)
  ! Update variables with limited states

end subroutine limiter_low_order


subroutine limiter_positivity(u)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, w_lim,u_lim
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w,w_nodes
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes, nodes_cons
  real(kind=8)::u_left,u_right,u_top,u_bottom, u_center, u_deriv
  real(kind=8)::limited1,limited2, generalized_minmod, minmod, legendre, dx
  integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
  integer::intnode, jntnode
  integer::left,right,bottom,top

  real(kind=8)::chsi_left=-1,chsi_right=+1
  dx = boxlen_x/dble(nx)
  ! look at 1st derivatives, u12, u21

  if (mx==1.and.my==1) then
    ! no limiting when we only have 1st order approx
    !use_limiter = .false.
    return
  end if

  ! mean is given by u11

  ! guarantee positivity of density and pressure

  call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
  call compute_primitive(nodes,w_nodes,nx,ny,mx,my)
  call get_modes_from_nodes(w_nodes,w,nx,ny,mx,my)

  !print*,'r',u(1,:,:,:,:)
  !print*,'p',u(4,:,:,:,:)
  !STOP
  u_lim = u
  !w_lim = w
  do ivar = 1,nvar
    !   if ((ivar == 1).or.(ivar==2).or.(ivar==3)) then
    do icell = 1,nx
      do jcell = 1,nx
        !do intnode = 1,1,-1
        left = icell -1
        right = icell+1
        if (icell == 1) then
          left = 1
        else if (icell == nx) then
          right = nx
        end if
        !print*,'before lim', u(1,icell,jcell,:,:)
        !print*, u(2,icell,jcell,:,:)
        !print*, u(3,icell,jcell,:,:)
        !print*, u(4,icell,jcell,:,:)
        !if (jcell == 1) then
        !   bottom = nx
        !else if (jcell == nx) then
        ! top = 1
        !end if

        u_left = 0.5*u(ivar,left,jcell, 1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_right =0.5*u(ivar,right,jcell, 1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_center = 0.5*u(ivar,icell,jcell, 1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_deriv =  u(ivar,icell,jcell, 2, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1+1)+1.0))
        u_lim(ivar,icell,jcell, 2, 1) = minmod(u_deriv,&
        &(u_center-u_left)/dx,(u_right-u_center)/dx)
        !if (ivar==1) then
        !print*,u_deriv, u_lim(ivar,icell,jcell,2,1)
        !  print*,'before lim', u(ivar,icell,jcell,:,:)
        !end if
        if(ABS(u_lim(ivar,icell,jcell, 2, 1)-u_deriv).GT.0.01*ABS(u_deriv)) then
          u_lim(ivar,icell,jcell,2:mx,1) = 0.0
          !u_lim(ivar,icell,jcell,mx,mx) = 0.0
        end if
      end do
    end do
    !end if
  end do
  ! end do

  do ivar = 1,nvar

    !if ((ivar == 1).or.(ivar==2).or.(ivar==3)) then
    do icell = 1,nx
      do jcell=1,nx
        !do intnode = mx-1,1,-1
        left = icell -1
        right = icell+1
        if (icell == 1) then
          left = 1
        else if (icell == nx) then
          right = nx
        end if
        !if (jcell == 1) then
        ! bottom = nx
        !else if (jcell == nx) then
        ! top = 1
        !end if

        u_left = 0.5*u(ivar,jcell,left,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_right =0.5*u(ivar,jcell,right,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_center = 0.5*u(ivar,jcell,icell,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_deriv = u(ivar,jcell,icell, 1, 2) !*sqrt(dble(2.0))/sqrt((2.0*dble(1+1)+1.0))
        u_lim(ivar,jcell, icell, 1, 2) = minmod(u_deriv,&
        &(u_center-u_left)/dx,(u_right-u_center)/dx)
        !   if  (u_center > 1000) then
        !  print*,(u_center-u_left),(u_right-u_center)
        ! print*, u_deriv, u_lim(ivar,icell,jcell,1,2)
        !  print*, 'beforelim',u(ivar,icell,jcell,1,2)
        ! end if !STOP
        if(ABS(u_lim(ivar,jcell,icell, 1,2)-u_deriv).GT.0.01*ABS(u_deriv)) then
          u_lim(ivar,jcell,icell, 1,2:my) = 0.0
          !u_lim(ivar,jcell,icell, mx,mx) = 0.0
        end if
      end do
      !STOP
    end do
    !end if
  end do
  !end do


  call get_nodes_from_modes(u_lim,nodes,nx,ny,mx,my)
  call compute_primitive(nodes, nodes_cons, nx, ny,mx,my )
  nodes = nodes_cons
  !nodes(4,:,:,:,:) = 10e-5
  !nodes(2,:,:,:,:) = 0.
  !nodes(3,:,:,:,:) = 1.


  !call compute_conservative(nodes,nodes_cons,nx,ny,mx,my)
  !call get_modes_from_nodes(nodes_cons,u_lim,nx,ny,mx,my)

  do ivar=1,nvar
    do icell=1,nx
      do jcell=1,ny
        u_left=0.0; u_right=0.0
        !Loop over modes
        !do intnode=1,mx
        do i = 1, mx
          do j = 1,my
            u_left= nodes(ivar,icell,jcell,1,j) !u_left+u_lim(1,icell,jcell,i,j)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)
            u_right= nodes(ivar,icell,jcell,mx,j)  !u_right+u_lim(1,icell,jcell,i,j)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)
            u_top = nodes(ivar,icell,jcell,i,1)
            u_bottom = nodes(ivar,icell,jcell,i,my)
            !if((u_left<1d-10).or.(u_left<1d-10))then
            if((u_left<1d-10.and.ivar==1).or.(u_left<1d-10.and.ivar==4))then
              nodes(ivar,icell,jcell,i,j)=1e-5
              !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            else if((u_right<1d-10.and.ivar==1).or.(u_right<1d-10.and.ivar==4))then
              !                         else if((u_right<1d-10).or.(u_right<1d-10))then
              nodes(ivar,icell,jcell,i,j)=1e-5
              !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            else if((u_top<1d-10.and.ivar==1).or.(u_top<1d-10.and.ivar==4))then
              !                        else if((u_top<1d-10).or.(u_top<1d-10))then
              nodes(ivar,icell,jcell,i,j)=1e-5
              !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            else if((u_bottom<1d-10.and.ivar==1).or.(u_bottom<1d-10.and.ivar==4))then
              !                          else if((u_bottom<1d-10).or.(u_bottom<1d-10))then
              nodes(ivar,icell,jcell,i,j)=1e-5
              !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            end if
          end do
        end do
      end do
    end do
  end do
  !nodes(4,:,:,:,:) = 10e-5
  !nodes(2,:,:,:,:) = 1.0
  !call get_nodes_from_modes(w_lim,nodes,nx,ny,mx,my)
  call compute_conservative(nodes,nodes_cons,nx,ny,mx,my)
  call get_modes_from_nodes(nodes_cons,u_lim,nx,ny,mx,my)
  !call compute_conservative(nodes,nodes_cons,nx,ny,mx,my)
  !call get_modes_from_nodes(nodes_cons, u_lim, nx, ny, mx, my)
  ! Update variables with limited states
  u = u_lim
end subroutine limiter_positivity

subroutine limiter_rossmanith(u)
  ! eats u, conservative modes and returns limited conservative modes
  ! performs limiting by looking at primitives (as described on paper)
  !
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, u_lim, u_lim_nodes
  real(kind=8),dimension(1:nvar,1:nx,1:ny)::w_average, u_average

  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u_nodes, w_nodes

  real(kind=8)::big_theta, little_theta, theta, w_average_stencil
  real(kind=8),dimension(1:nvar)::minimum_neigh_value, minimum_cell_value, min_temp, max_temp,&
  &maximum_neigh_value, maximum_cell_value, ratio1, ratio2
  real(kind=8),dimension(1:nvar)::maximum_neigh_value_1, minimum_neigh_value_1

  integer,dimension(1:4)::indices_i,indices_j
  real(kind=8)::u_left
  real(kind=8)::limited1,limited2, generalized_minmod, minmod, legendre, dx, ratio, phi_lim
  integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
  integer::intnode, jntnode, indx
  integer::left,right,bottom,top, component

  dx = boxlen_x/dble(nx)

  ! transform to nodes
  call get_nodes_from_modes(u, u_nodes, nx, ny, mx, my)
  call compute_primitive(u_nodes,w_nodes,nx,ny,mx,my)

  call get_modes_from_nodes(w_nodes, w, nx, ny, mx,my)


  ! get averages
  do ivar=1,nvar
    do icell = 1,nx
      do jcell = 1,ny
        w_average(ivar,icell,jcell) =sum(w_nodes(ivar,icell,jcell,:,:))/dble(mx+my)
        u_average(ivar,icell,jcell) =sum( u_nodes(ivar,icell,jcell,:,:))/dble(mx+my)
      end do
    end do
  end do
  ! compute maximum node
  do icell = 1,nx
    do jcell = 1,ny

      do ivar = 1,nvar
        maximum_cell_value(ivar) = maxval(w_nodes(ivar,icell,jcell,:,:))
        minimum_cell_value(ivar) = minval(w_nodes(ivar,icell,jcell,:,:))
        ! periodic bcs
        if (icell == 1) then
          indices_i = (/ nx, icell+1, icell, icell /)
        else if (icell == nx) then
          indices_i = (/ icell - 1, 1, icell, icell /)
        else
          indices_i = (/ icell - 1, icell+1, icell, icell /)
        end if

        if (jcell == 1) then
          indices_j = (/ jcell, jcell, ny, jcell+1 /)
        else if (jcell == ny) then
          indices_j = (/ jcell, jcell, jcell-1, 1 /)
        else
          indices_j = (/ jcell, jcell, jcell-1, jcell+1 /)
        end if

        max_temp(ivar) = 0
        min_temp(ivar) = 0
        maximum_neigh_value(ivar) = 0
        minimum_neigh_value(ivar) = 0

        w_average_stencil = 0.0
        do indx = 1,4
          max_temp(ivar) = maxval(w_nodes(ivar,indices_i(indx),indices_j(indx),:,:))
          min_temp(ivar) = minval(w_nodes(ivar,indices_i(indx),indices_j(indx),:,:))
          w_average_stencil = w_average_stencil + 0.5*w(ivar,indices_i(indx),indices_j(indx),1,1)
          if (max_temp(ivar) > maximum_neigh_value_1(ivar)) then
            maximum_neigh_value_1(ivar) = max_temp(ivar)
          end if
          if (min_temp(ivar) < minimum_neigh_value_1(ivar)) then
            minimum_neigh_value_1(ivar) = min_temp(ivar)
          end if
          ! check the corners too
        end do

        maximum_neigh_value(ivar) = max(w_average(ivar,icell,jcell),maximum_neigh_value_1(ivar))
        minimum_neigh_value(ivar) = min(w_average(ivar,icell,jcell),minimum_neigh_value_1(ivar))

        ! Heuristic f(h) = \alpha h**2 doesn't work very well
        ! Paper states that if f(h) = 0 we do the most agressive limiting and it clips maxima
        ratio1(ivar) = min(1e-6,(minimum_neigh_value(ivar)-w_average(ivar,icell,jcell)))/&
        &min(1e-6,dble(minimum_cell_value(ivar)-w_average(ivar,icell,jcell)))
        ratio2(ivar) = max(1e-6,(maximum_neigh_value(ivar)-w_average(ivar,icell,jcell)))/&
        &max(1e-6,dble(maximum_cell_value(ivar)-w_average(ivar,icell,jcell)))

      end do
      little_theta = 1.
      big_theta = 1.

      do component = 1,nvar
        if (little_theta>phi_lim(ratio1(component))) then
          little_theta=phi_lim(ratio1(component))
        end if
        if (big_theta>phi_lim(ratio2(component))) then
          big_theta=phi_lim(ratio2(component))
        end if
      end do

      theta = max(0.,min(dble(1),min(little_theta,big_theta)))

      ! There's no guarantee that theta is between -1,1
      if (theta>0.90) then
        theta = 0.0
      else if (theta < -.5) then
        theta = 0.0
      end if

      do ivar = 1,nvar
        do intnode = 1,mx
          do jntnode = 1,my
            u_lim_nodes(ivar,icell,jcell,intnode,jntnode) = u_average(ivar,icell,jcell) +&
            &theta*(u_nodes(ivar,icell,jcell,intnode,jntnode)-u_average(ivar,icell,jcell))
          end do
        end do
      end do
    end do
  end do

  call compute_primitive(u_lim_nodes, w_nodes,nx,ny,mx,my)

  !Enforce positivity'

  do ivar=1,nvar
    do icell=1,nx
      do jcell=1,ny
        u_left=0.0!; u_right=0.0
        do i = 1, mx
          do j = 1,my
            u_left= w_nodes(ivar,icell,jcell,i,j) !u_left+u_lim(1,icell,jcell,i,j)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)
            !u_right= nodes(1,icell,jcell,i,j)  !u_right+u_lim(1,icell,jcell,i,j)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)
            if((u_left<1d-6.and.ivar==1).or.(u_left<1d-6.and.ivar==4))then
              w_nodes(ivar,icell,jcell,i,j)=10e-5
              !u_lim(1,icell,jcell,2:mx,2:my)=0.0
            end if
          end do
        end do
      end do
    end do
  end do


  call compute_conservative(w_nodes,u_lim_nodes,nx,ny,mx,my)
  call get_modes_from_nodes(u_lim_nodes,u_lim,nx,ny,mx,my)

  u = u_lim
end subroutine limiter_rossmanith

real(kind=8) function phi_lim(y)
  real(kind=8)::y
  phi_lim = min(1.,y/dble(1.1)) ! as per indication of Rossmanith et al paper.
end function phi_lim

subroutine Krivodonova(modes)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u_lim,modes,nodes,chars,char_modes
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w,w_modes
  real(kind=8),dimension(1:nvar,1:mx,1:my)::uL,uR,uT,uB,uM, w_lim, char_mode
  real(kind=8),dimension(1:nvar,1,1,1:mx,1:my)::char_mode_1,w_lim_1
  real(kind=8),dimension(1:nvar)::w_min

  integer::ileft,iright,icell,jcell,i,j,ivar, itop,ibottom
  real(kind=8)::coeff_i,coeff_ip1, minmod

  call get_nodes_from_modes(modes,nodes,nx,ny,mx,my)
  call compute_primitive(nodes,w,nx,ny,mx,my)
  call compute_characteristics(nodes,chars,nx,ny,mx,my)
  call get_modes_from_nodes(chars,char_modes,nx,ny,mx,my)
  ! ? ? ? ? ?
  print*,chars(1,:,:,1,1)
  print*,char_modes(1,:,:,1,1)
  ! print*,'here'
  ! limiting
  u_lim = modes
  do icell=1,nx
    do jcell=1,ny
      ileft=icell-1
      iright=icell+1
      itop = jcell+1
      ibottom = jcell - 1

      call get_boundary_conditions(ileft,1)
      call get_boundary_conditions(iright,1)
      call get_boundary_conditions(itop,2)
      call get_boundary_conditions(ibottom,2)
      print*,'here'
      ! Compute primitive variable for all modes
      ! Loop only over high-order modes
      do i=mx-1,1,-1
        ! we don't care about cross terms for now
        j = i
        ! Renormalise to get proper Legendre polynomials
        ! and corresponding derivatives
        coeff_i=2/sqrt(2.0*dble(i-1)+1.0)*sqrt(2.0*dble(j-1)+1.0)!*(2.0*dble(i)-1)
        coeff_ip1=2/sqrt(2.0*dble(i)+1.0)*sqrt(2.0*dble(j)+1.0)
        uL(1:nvar,i+1,j+1)=(char_modes(1:nvar,icell,jcell,i,j)-char_modes(1:nvar,ileft,jcell,i,j))
        uR(1:nvar,i+1,j+1)=(char_modes(1:nvar,iright,jcell,i,j)-char_modes(1:nvar,icell,jcell,i,j))
        uT(1:nvar,i+1,j+1)=(char_modes(1:nvar,icell,itop,i,j)-char_modes(1:nvar,ileft,jcell,i,j))
        uB(1:nvar,i+1,j+1)=(char_modes(1:nvar,icell,jcell,i,j)-char_modes(1:nvar,ileft,ibottom,i,j))
        uM(1:nvar,i+1,j+1)=char_modes(1:nvar,icell,jcell,i+1,j+1)
        ! print*,  uL(1:nvar,i+1,j+1),  uR(1:nvar,i+1,j+1),  uB(1:nvar,i+1,j+1)
      end do
      print*,'here'
      w_lim=uM
      ! Loop over variables
      do ivar=1,nvar
        ! Loop only over high-order modes
        do i=mx-1,1,-1
          j = i
          w_min(ivar) = minmod(minmod(uL(ivar,i+1,j+1),uM(ivar,i+1,j+1),&
          &uR(ivar,i+1,j+1)),uT(ivar,i+1,j+1),uB(ivar,i+1,j+1))
          w_lim(ivar,i+1,j+1)=w_min(ivar)

          !(ABS(w_min-uM(ivar,i+1,j+1)).LT.0.01*ABS(uM(ivar,i+1,j+1)))exit
        end do
        ! End loop over modes
      end do
      ! End loop over variables
      ! Compute conservative variables
      ! Loop only over high-order modes
      !do i=mx-1,1,-1
      !    j = i
      !    print*,i,j
      print*,'here before'
      w_lim_1(1:nvar,1,1,:,:) =  w_lim(1:nvar,:,:)
      char_mode_1(1:nvar,1,1,:,:) =  w_lim(1:nvar,:,:)
      print*,w_lim_1
      print*,char_mode_1
      call get_nodes_from_modes(w_lim_1(:,1,1,:,:),char_mode_1(:,1,1,:,:),1,1,mx,my)
      print*,'here after nodes'
      do i=mx-1,1,-1
        j = i
        print*,i,j
        if ((j>2).or.(i>2)) then
          print*,'duh'
        endif
        call compute_cons_from_characteristics(char_mode_1(1:nvar,1,1,i+1,j+1),&
        &nodes(1:nvar,icell,jcell,i+1,j+1),&
        &u_lim(1:nvar,icell,jcell,i+1,j+1),1,1,1,1)
      end do
      print*,'here'
    end do
  end do
  ! End loop over cells

  ! positivity

end subroutine Krivodonova

subroutine limiter_cockburn(u)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, w_lim,u_lim
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w,w_nodes,u_temp,chars,chars_modes
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes, nodes_cons, nodes_new
  real(kind=8)::u_left,u_right,u_top,u_bottom, u_center, u_deriv
  real(kind=8)::limited1,limited2, generalized_minmod, minmod, legendre, dx
  integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
  integer::intnode, jntnode
  integer::left,right,bottom,top

  real(kind=8)::chsi_left=-1,chsi_right=+1
  dx = boxlen_x/dble(nx)
  ! look at 1st derivatives, u12, u21

  if (mx==1.and.my==1) then
    ! no limiting when we only have 1st order approx
    !use_limiter = .false.
    return
  end if

  ! mean is given by u11

  ! guarantee positivity of density and pressure

  call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
  call compute_primitive(nodes,w_nodes,nx,ny,mx,my)
  call compute_characteristics(u,chars,nx,ny,mx,my)
  call get_modes_from_nodes(chars,chars_modes,nx,ny,mx,my)
  call get_modes_from_nodes(w_nodes,w,nx,ny,mx,my)

  u_temp = u
  u = chars_modes
  u_lim = chars_modes
  !w_lim = w
  do ivar = 1,nvar
    !   if ((ivar == 1).or.(ivar==2).or.(ivar==3)) then
    do icell = 1,nx
      do jcell = 1,nx
        !do intnode = 1,1,-1
        left = icell -1
        right = icell+1
        if (icell == 1) then
          left = nx
        else if (icell == nx) then
          right = 1
        end if
        u_left = u(ivar,left,jcell, 1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_right = u(ivar,right,jcell, 1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_center = u(ivar,icell,jcell, 1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_deriv =  u(ivar,icell,jcell, 2, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1+1)+1.0))
        u_lim(ivar,icell,jcell, 2, 1) = 2*minmod(u_deriv,&
        &(u_center-u_left),(u_right-u_center))
        if(ABS(u_lim(ivar,icell,jcell, 2, 1)-u_deriv).GT.0.01*ABS(u_deriv)) then
          u_lim(ivar,icell,jcell,2:mx,1) = 0.0
          u_lim(ivar,icell,jcell,mx,mx) = 0.0
        end if
      end do
    end do
    !end if
  end do
  ! end do

  do ivar = 1,nvar
    do icell = 1,nx
      do jcell=1,nx
        left = icell -1
        right = icell+1
        if (icell == 1) then
          left = nx
        else if (icell == nx) then
          right = 1
        end if
        u_left = u(ivar,jcell,left,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_right = u(ivar,jcell,right,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_center = u(ivar,jcell,icell,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
        u_deriv = u(ivar,jcell,icell, 1, 2) !*sqrt(dble(2.0))/sqrt((2.0*dble(1+1)+1.0))
        u_lim(ivar,jcell, icell, 1, 2) = 2*minmod(u_deriv,&
        &(u_center-u_left),(u_right-u_center))
        if(ABS(u_lim(ivar,jcell,icell, 1,2)-u_deriv).GT.0.01*ABS(u_deriv)) then
          u_lim(ivar,jcell,icell, 1,2:my) = 0.0
          u_lim(ivar,icell,jcell,mx,mx) = 0.0
        end if
      end do
    end do
  end do

  call get_nodes_from_modes(u_lim,chars,nx,ny,mx,my)
  call compute_cons_from_characteristics(chars,u_temp,nodes_new,nx,ny,my,mx)
  call compute_primitive(nodes_new, nodes_cons, nx, ny,mx,my )
  nodes = nodes_cons
  nodes_cons(4,:,:,:,:) = 10e-5
  call compute_conservative(nodes_cons,nodes_new,nx,ny,mx,my)
  call get_modes_from_nodes(nodes_new,u_lim,nx,ny,mx,my)

  u = u_lim
  call get_modes_from_nodes(u,u_lim,nx,ny,mx,my)
  call compute_primitive(u_lim, nodes_cons, nx, ny,mx,my )
  print*,'p',nodes_cons(4,:,:,:,:)
  print*,'u',nodes_cons(2,:,:,:,:)
  print*,'v',nodes_cons(3,:,:,:,:)
  PAUSE

  do ivar=1,nvar
    do icell=1,nx
      do jcell=1,ny
        do i = 1, mx
          do j = 1,my
            u_left= nodes(ivar,icell,jcell,1,j)
            u_right= nodes(ivar,icell,jcell,2,j)
            u_top = nodes(ivar,icell,jcell,i,1)
            u_bottom = nodes(ivar,icell,jcell,i,my)
            !if((u_left<1d-10).or.(u_left<1d-10))then
            if((u_left<1d-10.and.ivar==1).or.(u_left<1d-10.and.ivar==4))then
              nodes(ivar,icell,jcell,i,j)=1e-5
              !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            end if
            if((u_right<1d-10.and.ivar==1).or.(u_right<1d-10.and.ivar==4))then
              !                         else if((u_right<1d-10).or.(u_right<1d-10))then
              nodes(ivar,icell,jcell,i,j)=1e-5
            end if
            !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            if((u_top<1d-10.and.ivar==1).or.(u_top<1d-10.and.ivar==4))then
              !                        else if((u_top<1d-10).or.(u_top<1d-10))then
              nodes(ivar,icell,jcell,i,j)=1e-5
            end if
            !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            if((u_bottom<1d-10.and.ivar==1).or.(u_bottom<1d-10.and.ivar==4))then
              !                          else if((u_bottom<1d-10).or.(u_bottom<1d-10))then
              nodes(ivar,icell,jcell,i,j)=1e-5
              !w_lim(ivar,icell,jcell,2:mx,2:my)=0.0
            end if
          end do
        end do
      end do
    end do
  end do

  call compute_conservative(nodes,nodes_cons,nx,ny,mx,my)
  call get_modes_from_nodes(nodes_cons,u_lim,nx,ny,mx,my)

  !u = u_lim
end subroutine limiter_cockburn
