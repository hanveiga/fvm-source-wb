function minmod(x,y,z)
  implicit none
  ! input
  real(kind=8)::x,y,z

  ! output
  real(kind=8)::minmod

  ! internal
  real(kind=8)::s

  s = sign(1d0,x)

  if (sign(1d0,y) == s .AND. sign(1d0,z) == s) then
     minmod = s*min(abs(x),abs(y),abs(z))
  else
     minmod = 0.0
  endif
  !write(*,*) minmod
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
  call get_modes_from_nodes(nodes, u, nx, ny, mx, my)
  ! Update variables with limited states

  end subroutine compute_limiter


subroutine limiter_low_order(u)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes
    real(kind=8)::limited1,limited2, generalized_minmod
    integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
    integer::intnode, jntnode

    ! look at 1st derivatives, u12, u21

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

          do intnode = 2,mx
            jntnode = intnode-1

            limited1 = generalized_minmod(u(ivar,icell,jcell,jntnode,intnode), &
            &u(ivar,iright,jcell,jntnode,jntnode)-u(ivar,icell,jcell,jntnode,jntnode),&
            &u(ivar,icell,jcell,jntnode,jntnode)-u(ivar,ileft,jcell,jntnode,jntnode))

            limited2 = generalized_minmod(u(ivar,icell,jcell,intnode,jntnode),&
            &u(ivar,icell,itop,1,1)-u(ivar,icell,jcell,1,1),&
            &u(ivar,icell,jcell,1,1)-u(ivar,icell,ibottom,1,1))

            if (abs(limited1 - u(ivar,icell,jcell,jntnode,intnode))<1e-5&
              &.and.abs(limited2 - u(ivar,icell,jcell,2,1))<1e-5) then
              u(ivar,icell,jcell,jntnode,intnode) = limited1
              u(ivar,icell,jcell,intnode,jntnode) = limited2
            else
              u(ivar,icell,jcell,1,2:my) = 0.0
              u(ivar,icell,jcell,2:mx,1) = 0.0
            end if
          end do
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
  call get_modes_from_nodes(nodes, u, nx, ny, mx, my)
  ! Update variables with limited states

  end subroutine limiter_low_order