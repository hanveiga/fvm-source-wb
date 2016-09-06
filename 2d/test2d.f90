program dg
  use parameters_dg_2d
  implicit none
  ! Main variables
  real(kind=8)::legendre,legendre_prime,dx,xcell,x, ycell, y, dy
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::modes
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes

  real(kind=8),dimension(2,1:nx,1:ny,1:mx,1:my)::flux_vol
  integer::ivar, icell,jcell,i,j,intnode,iter,jntnode

  u(:,:,:,:,:) = 0.0
  modes(:,:,:,:,:) = 0.0
  nodes(:,:,:,:,:) = 0.0

  call gl_quadrature(x_quad,w_x_quad,mx)
  call gl_quadrature(y_quad,w_y_quad,my)

  dx = 1./nx
  dy = 1./ny

  do ivar = 1,nvar
  do icell=1,nx
    do jcell =1,ny
      do i = 1,mx
        do j = 1,my
          xcell = (icell-0.5)*dx
          ycell = (jcell-0.5)*dy
          x = xcell + dx/dble(2)*x_quad(i)
          y = ycell + dy/dble(2)*y_quad(j)
          u(ivar,icell,jcell,i,j) = exp(-x+y)
        end do
      end do
    end do
  end do
  end do 

  do ivar =1,nvar
  do icell = 1,nx
    do jcell = 1,ny
      do i = 1,mx
        do j = 1,my
          do intnode = 1,mx
            do jntnode = 1,my
              modes(ivar,icell,jcell,i,j) = modes(ivar,icell,jcell,i,j) + &
                                         & 0.25*u(ivar,icell,jcell,intnode,jntnode)*&
                                         & legendre(x_quad(intnode),i-1)* &
                                         & legendre(y_quad(jntnode),j-1)* &
                                         & w_x_quad(intnode)*w_y_quad(jntnode)
            end do
          end do
        end do
      end do
    end do
  end do
  end do 

  do icell = 1,nx
    do jcell = 1,ny
      do i = 1,mx
        do j = 1,my
          do intnode = 1,mx
            do jntnode = 1,my
              nodes(1:nvar,icell,jcell,i,j) = nodes(1:nvar,icell,jcell,i,j) + &
                                         & modes(1:nvar,icell,jcell,intnode,jntnode)*&
                                         & legendre(x_quad(i),intnode-1)* &
                                         & legendre(y_quad(j),jntnode-1)
            end do
          end do
        end do
      end do
    end do
  end do


  do iter = 1,1000
        call get_modes_from_nodes(nodes, modes, nx,ny,mx,my)
        call get_nodes_from_modes(modes,nodes,nx,ny,mx,my)
  end do


  write(*,*) 'max diff', maxval(u - nodes)
  write(*,*) 'min fiff', minval(u - nodes)
!  write(*,*) modes(1,1,1,:,:)
end program dg
