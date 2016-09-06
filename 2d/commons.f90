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
                & 0.25*nodes(1:nvar,icell,jcell,xquad,yquad)* &
                & legendre(x_quad(xquad),i-1)* &
                & legendre(y_quad(yquad),j-1)* &
                & w_x_quad(xquad) * &
                & w_y_quad(yquad)! *dx*dy/4.
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
