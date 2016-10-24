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



subroutine high_order_limiter_2(u)
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

end subroutine high_order_limiter_2


subroutine compute_limiter(u)
   use parameters_dg_2d
   implicit none
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w, w_nodes
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes
   real(kind=8)::limited1,limited2, generalized_minmod, u_left
   real(kind=8)::beta, norm
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
   beta = 1.
   norm = 3.

   call get_nodes_from_modes(u,nodes, nx,ny,mx,my)
   call compute_primitive(nodes, w_nodes, nx, ny, mx, my)
   call get_modes_from_nodes(w_nodes,w,nx,ny,mx,my)

   modes(:,:,:,:,:)=0
   modes(:,:,:,1,1)=w(:,:,:,1,1)
   u = w


   !call compute_positivity(u)


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


         limited1 = generalized_minmod(norm*u(ivar,icell,jcell,2,1), &
         &(u(ivar,iright,jcell,1,1)-u(ivar,icell,jcell,1,1)),&
         &(u(ivar,icell,jcell,1,1)-u(ivar,ileft,jcell,1,1)))/norm

         limited2 = generalized_minmod(norm*u(ivar,icell,jcell,1,2),&
         &(u(ivar,icell,itop,1,1)-u(ivar,icell,jcell,1,1)),&
         &(u(ivar,icell,jcell,1,1)-u(ivar,icell,ibottom,1,1)))/norm

         if((ABS(limited1 - u(ivar,icell,jcell,2,1)).GT.1E-6) .or. ABS(limited2 - u(ivar,icell,jcell,1,2)).GT.1E-6) then
           modes(ivar,icell,jcell,2,1)=limited1
           modes(ivar,icell,jcell,1,2)=limited2
         else
           modes(ivar,icell,jcell,1,2)=u(ivar,icell,jcell,1,2)
           modes(ivar,icell,jcell,2,1)=u(ivar,icell,jcell,2,1)
         end if
      end do
    end do
  end do
 end if

   call get_nodes_from_modes(modes,nodes,nx,ny,mx,my)

   !call compute_primitive(nodes, w_nodes,nx,ny,mx,my)

   !Enforce positivity'

   !do ivar=1,nvar
    ! do icell=1,nx
    !   do jcell=1,ny
    !     u_left=0.0!; u_right=0.0
    !     do i = 1, mx
    !       do j = 1,my
    !         u_left= w_nodes(ivar,icell,jcell,i,j) !u_left+u_lim(1,icell,jcell,i,j)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)
             !u_right= nodes(1,icell,jcell,i,j)  !u_right+u_lim(1,icell,jcell,i,j)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)
    !         if((u_left<1d-6.and.ivar==1).or.(u_left<1d-6.and.ivar==4))then
    !           w_nodes(ivar,icell,jcell,i,j)=10e-5
               !u_lim(1,icell,jcell,2:mx,2:my)=0.0
    !         end if
    !      end do
    !     end do
    !   end do
    ! end do
   !end do
!w_nodes(2,:,:,:,:) = 0.
!w_nodes(3,:,:,:,:) = 1.
!w_nodes(4,:,:,:,:) = 10e-5
call get_nodes_from_modes(modes,w_nodes,nx,ny,mx,my)
call compute_conservative(w_nodes,nodes,nx,ny,mx,my)

! guarantee positivity of density and pressure
!call get_nodes_from_modes(modes,nodes,nx,ny,mx,my)
call get_modes_from_nodes(nodes, u, nx, ny, mx, my)
! Update variables with limited states
!call limiter_pp(u)

end subroutine compute_limiter


real(kind=8) function solve_for_t(u,u_avg)
  use parameters_dg_2d
  implicit none
  real(kind=8)::p,a,b,c, root2, root1,t1,t2,t,D,q
  real(kind=8)::pa,mxa,mya,ea,pj,mxj,myj,ej,xx
  real(kind=8),dimension(1:nvar)::u, u_avg
  integer::iter
  real(kind=8)::eps
  eps=10e-10
  ! compute with mapple <3
  !a=(gamma-1)*((-pa+pj)*(-Ea+Ej)-(1/2)*(-mxa+mxj)^2-(1/2)*(-mya+myj)^2)
  !b=
  !c=
  !print*,u,u_avg
  pa = u_avg(1); mxa = u_avg(2); mya = u_avg(3); ea = u_avg(4)
  pj = u(1); mxj = u(2); myj = u(3); ej = u(4)
  !a = (gamma-1)*((-pa+pj)*(-Ea+Ej)-(1./2.)*(-mxa+mxj)**2-(1./2.)*(-mya+myj)**2)
  !b = (gamma-1)*(pa*(-Ea+Ej)+(-pa+pj)*Ea-mxa*(-mxa+mxj)-mya*(-mya+myj))-eps*(-pa+pj)
  !c = (gamma-1)*(pa*Ea-(1./2.)*mxa**2-(1./2.)*mya**2)-eps*pa
  !print*,a,b,c
  !a = 2.0*(pj-pa)*(ej-ea) - (mxj-mxa)**2 - (myj-mya)**2
  !b = 2.0*(pj-pa)*(ea-eps/(gamma-1)) + 2.0*pa*(ej-ea) - 2.0*(mxa*(mxj-mxa)+mya*(myj-mya))
  !c = 2.0*pa*ea - (mxa**2+mya**2) - 2.0*eps*pa/(gamma-1.0)

  a = 2.0*(pj-pa)*(ej-ea) - (mxj-mxa)**2 - (myj-mya)**2
  b = 2.0*(pj-pa)*(ea-eps/(gamma-1)) + 2.0*pa*(ej-ea) - 2.0*(mxa*(mxj-mxa)+mya*(myj-mya))
  c = 2.0*pa*ea - (mxa**2+mya**2) - 2.0*eps*pa/(gamma-1.0)
  !q = -1./2.*(b+sign(1d0,b)*sqrt(b**2-4*a*c))
  !t1 = q/max(a,eps)!0.5*(-b-D)
  !t2 = c/max(q,eps)!0.5*(-b+D)
  !print*,t1,t2
  b = b/a
  c = c/a
  D = sqrt(abs(b*b-4*c))
  t1 = 0.5*(-b-D)
  t2 = 0.5*(-b+D)
  !print*,t1,t2
  if((t1 > -eps) .and. (t1 < 1.0 + eps)) then
    t = t1
  else if((t2 > -eps).and.(t2 < 1.0 + eps)) then
    t = t2
  else
    print*,'error, setting t to zero'
    t = 0.0
  end if
  !pause
  t = min(1.0,t)
  t = max(0.0,t)
  solve_for_t = t
  return
end function solve_for_t



real(kind=8) function solve_for_t_iter(u,u_avg)
  use parameters_dg_2d
  implicit none
  real(kind=8)::p,a,b,c, root2, root1,t1,t2,t,D,q, x_old,x_new
  real(kind=8)::pa,mxa,mya,ea,pj,mxj,myj,ej,xx
  real(kind=8),dimension(1:nvar)::u, u_avg
  integer::iter
  real(kind=8)::eps
  eps=10e-10
  ! compute with mapple <3
  !a=(gamma-1)*((-pa+pj)*(-Ea+Ej)-(1/2)*(-mxa+mxj)^2-(1/2)*(-mya+myj)^2)
  !b=
  !c=
  !print*,u,u_avg
  pa = u_avg(1); mxa = u_avg(2); mya = u_avg(3); ea = u_avg(4)
  pj = u(1); mxj = u(2); myj = u(3); ej = u(4)
  !a = (gamma-1)*((-pa+pj)*(-Ea+Ej)-(1./2.)*(-mxa+mxj)**2-(1./2.)*(-mya+myj)**2)
  !b = (gamma-1)*(pa*(-Ea+Ej)+(-pa+pj)*Ea-mxa*(-mxa+mxj)-mya*(-mya+myj))-eps*(-pa+pj)
  !c = (gamma-1)*(pa*Ea-(1./2.)*mxa**2-(1./2.)*mya**2)-eps*pa
  !print*,a,b,c
  !a = 2.0*(pj-pa)*(ej-ea) - (mxj-mxa)**2 - (myj-mya)**2
  !b = 2.0*(pj-pa)*(ea-eps/(gamma-1)) + 2.0*pa*(ej-ea) - 2.0*(mxa*(mxj-mxa)+mya*(myj-mya))
  !c = 2.0*pa*ea - (mxa**2+mya**2) - 2.0*eps*pa/(gamma-1.0)

  a = 2.0*(pj-pa)*(ej-ea) - (mxj-mxa)**2 - (myj-mya)**2
  b = 2.0*(pj-pa)*(ea-eps/(gamma-1)) + 2.0*pa*(ej-ea) - 2.0*(mxa*(mxj-mxa)+mya*(myj-mya))
  c = 2.0*pa*ea - (mxa**2+mya**2) - 2.0*eps*pa/(gamma-1.0)
  !q = -1./2.*(b+sign(1d0,b)*sqrt(b**2-4*a*c))
  !t1 = q/max(a,eps)!0.5*(-b-D)
  !t2 = c/max(q,eps)!0.5*(-b+D)
  !print*,t1,t2
  b = b/a
  c = c/a
  D = sqrt(abs(b*b-4*c))
  t1 = 0.5*(-b-D)
  t2 = 0.5*(-b+D)
  !print*,t1,t2
  if((t1 > -eps) .and. (t1 < 1.0 + eps)) then
    t = t1
  else if((t2 > -eps).and.(t2 < 1.0 + eps)) then
    t = t2
  else
    print*,'error, setting t to zero'
    t = 0.0
  end if
  !pause
  t = min(1.0,t)
  t = max(0.0,t)
  !print*,'t:', t

  x_old = t !use root result as initial guess
  x_new = 0.0
  do while (abs(x_new-x_old) > eps)
    x_new = 2*(x_old**2 + c*a)/(2*x_old+b*a)
    x_old = x_new
    iter = iter + 1
    if (iter > 1000) then
      x_new = t
      x_old = t
      print*,'setting t to zero'
    end if
  end do
  !print*,'x_old 1',x_old

  !x_old = max(x_old,c/(a*x_old))
  !print*,'x_old 2',x_old
  solve_for_t_iter = x_old
  return
end function solve_for_t_iter



subroutine compute_set(modes, node_points)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:mx,1:my)::modes, nodes
  real(kind=8),dimension(1:nvar,1:(mx)*gll+(my)*gll)::node_points
  real(kind=8),dimension(1:nvar)::u_left,u_right, u_top, u_bottom
  real(kind=8)::legendre
  integer::ivar,intnode,i,j, jntnode
  integer::chsi_right,chsi_left
  call get_nodes_from_modes(modes,nodes,1,1,mx,my)

  chsi_right = 1
  chsi_left = -1

  call gl_quadrature(x_quad,w_x_quad,mx)
  call gl_quadrature(y_quad,w_y_quad,my)
  call gll_quadrature(x_gll_quad,w_gll_quad,gll)

  do intnode = 1,my
    do jntnode = 1,gll
      u_left = 0.0
      u_right = 0.0
      do i=1,mx
        do j=1,my
             u_left = u_left + modes(1:nvar,i,j)*legendre(x_gll_quad(jntnode),i-1)*legendre(y_quad(intnode),j-1)
             u_right = u_right + modes(1:nvar,i,j)*legendre(x_quad(intnode),i-1)*legendre(x_gll_quad(jntnode),j-1)
        end do
      end do
      node_points(1:nvar,(intnode-1)*gll+jntnode) = u_left(1:nvar)
      node_points(1:nvar,(intnode-1)*gll+jntnode+my*gll) = u_right(1:nvar)
      !print*,(intnode-1)*gll+jntnode, node_points(1:nvar,(intnode-1)*gll+jntnode)
      !print*,(intnode-1)*gll+jntnode+my*gll, node_points(1:nvar,(intnode-1)*gll+jntnode+my*gll)
    end do
  end do
  !print*,node_points
  !pause

end subroutine compute_set


subroutine compute_positivity(u)
   ! take in conservative modes
   ! implementation of positivity preserving limiter Shu JCP 2010
   ! http://ac.els-cdn.com/S0021999110004535/1-s2.0-S0021999110004535-main.pdf?_tid=a2ea7790-8575-11e6-9c4d-00000aab0f6b&acdnat=1475065277_313c5c41e3c79fb7133c1186b8f179d7

   use parameters_dg_2d
   implicit none
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w, pos_w
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes, pos_nodes, pos_nodes_n, pos_modes
   real(kind=8),dimension(1:nvar,1:nx,1:ny)::u_avg
   real(kind=8),dimension(1:mx,1:my)::minimizing_set
   real(kind=8),dimension(1:nvar)::u_left,u_right,w_left,w_right, u_top, u_bottom
   integer::chsi_right,chsi_left, set_to_zero
   real(kind=8),dimension(1:nvar,1:(mx*gll+my*gll))::face_vals,w_face_vals

   integer::i,j,icell,jcell, ivar, intnode,jntnode
   real(kind=8)::t,t_min,theta, p_min, legendre
   real(kind=8)::eps, solve_for_t_iter
   eps=10e-10
   chsi_right = 1
   chsi_left = -1
   set_to_zero = 0

    call gl_quadrature(x_quad,w_x_quad,mx)
    call gl_quadrature(y_quad,w_y_quad,my)

   if (mx==1.and.my==1) then
     ! no limiting when we only have 1st order approx
     !use_limiter = .false.
     return
   end if

   ! 1. Limit density positivity
   call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
   call compute_primitive(nodes,w,nx,ny,mx,my)

   ! compute average of states
   u_avg(1:nvar,:,:) = u(1:nvar,:,:,1,1)
   !print*,'min mean density entrance', minval(u_avg(1,:,:)),  maxval(u_avg(1,:,:))
   !pause
   pos_modes = u
   pos_nodes = nodes
   pos_nodes_n  = pos_nodes
     do icell=1,nx
       do jcell = 1,ny
         call compute_set(u(1:nvar,icell,jcell,:,:),face_vals)
         ! evaluate polynomial at other quadrature nodes too and extend minizing set


         p_min = minval(face_vals(1,:))
         theta = min(abs((u_avg(1,icell,jcell)-eps)/(u_avg(1,icell,jcell)-p_min)),dble(1))

         do i = 1,mx
           do j = 1,my
             !if ((nodes(1,icell,jcell,i,j)<eps).or.(u_avg(1,icell,jcell)<eps)) then
            !   pos_nodes(1,icell,jcell,i,j) = 1e-5
            ! else
               !pos_nodes(1,icell,jcell,i,j) = u_avg(1,icell,jcell) + theta*(nodes(1,icell,jcell,i,j)-u_avg(1,icell,jcell))
               !if ((icell==31).and.(jcell==51)) then
                 !print*,'AAAA'
                ! print*,'i,j',i,j
                 !print*,'u', u(1:nvar,icell,jcell,i,j)
                 !pause
               !end if
               if ((i.ne.1).or.(j.ne.1)) then
                 !exit
                 !print*,i,j
                 pos_modes(1,icell,jcell,i,j) = theta*u(1,icell,jcell,i,j)
               end if
               !u_avg(1,icell,jcell) + theta*(modes(1,icell,jcell,i,j)-u_avg(1,icell,jcell))

               !if ( pos_nodes(1,icell,jcell,i,j)<eps) then
                 !print*, 'pos nodes negative?',pos_nodes(1,icell,jcell,i,j)
               !end if
            ! end if
           end do
         end do
       end do
     end do
   ! 2. limit pressure

   !pos_nodes_n = pos_nodes
   !if (set_to_zero==1) then
   !pos_nodes_n = pos_nodes
   !call get_modes_from_nodes(pos_nodes,u,nx,ny,mx,my)
   t_min = 1.
   !call compute_primitive(pos_nodes,w,nx,ny,mx,my)
   !pos_modes_n = pos_modes
   u = pos_modes
   !print*,'limited den'
   !do ivar = 1,nvar
     do icell=1,nx
      do jcell = 1,ny

        !do i=1,mx
        !  do j=1,my
        !     if (w(4,icell,jcell,i,j) > eps) then
        !       t = 1.
        !     else
               !t = solve_for_t(pos_nodes(4,icell,jcell,i,j),u_avg(4,icell,jcell))
               !print*,t
        !       t=0.0
        !     end if
        !     if (t_min>t) then
        !       t_min = t
        !     end if
        !   end do
         !end do

         call compute_set(u(1:nvar,icell,jcell,:,:),face_vals)
         do i = 1,gll*mx+gll*my
           call compute_primitive(face_vals(1:nvar,i),w_face_vals(1:nvar,i),1,1,1,1)

           !minimizing_set(:,:) = nodes(1,icell,jcell,:,:)
           !print*,face_vals(1:nvar,i)
           !print*, w_face_vals(1:nvar,i)
            !print*, w_face_vals(4,i)
           if (w_face_vals(4,i) > eps) then
                t = 1.
           else
                t = solve_for_t_iter(face_vals(1:nvar,i),u_avg(1:nvar,icell,jcell))
                !print*,t
                !pause
                !t = 0.0
           end if
           if (t_min>=t) then
               t_min = t
           end if
         end do
        !pause


       do i=1,mx
        do j=1,my
          !pos_nodes_n(1:nvar,icell,jcell,i,j) = u_avg(1:nvar,icell,jcell) + &
          !          &t_min*(pos_nodes(1:nvar,icell,jcell,i,j)-u_avg(1:nvar,icell,jcell))
          if ((i.ne.1).or.(j.ne.1)) then
            pos_modes(1:nvar,icell,jcell,i,j) = t_min*(u(1:nvar,icell,jcell,i,j))
          end if
          !if ((icell==31).and.(jcell==51)) then
              !print*,'AAAA'
            !print*,'jcell,icell,avg,theta',icell,jcell,t_min
            !print*,'rho min', p_min
          !  print*,i,j
          !  print*,'modes before', (u(1:nvar,icell,jcell,i,j))
          !  print*,'modes after', pos_modes(1:nvar,icell,jcell,i,j)
          !  pause
          !end if

         end do
       end do


       t_min = 1.

    end do
  end do
  !call get_modes_from_nodes(pos_nodes_n,modes,nx,ny,mx,my)
  call get_nodes_from_modes(pos_modes,pos_nodes_n,nx,ny,mx,my)

  call compute_primitive(pos_nodes_n,pos_w,nx,ny,mx,my)
  !print*,'maxpressure',maxval(w(4,:,:,:,:)),minval(w(4,:,:,:,:))
  !print*,'maxpressure_limited',maxval(pos_w(4,:,:,:,:)),minval(pos_w(4,:,:,:,:))
  !print*,'velocities',maxval(pos_w(2,:,:,:,:)),minval(pos_w(2,:,:,:,:)),&
  !&maxval(pos_w(3,:,:,:,:)),minval(pos_w(3,:,:,:,:))
  !print*,'density',maxval(pos_w(1,:,:,:,:)),minval(pos_w(1,:,:,:,:))
  !print*,'means',minval(modes(1,:,:,1,1)),minval(modes(4,:,:,1,1))
  !pos_w(3,:,:,:,:) = 1.0
  !pos_w(2,:,:,:,:) = 0.0
  !pos_w(4,:,:,:,:) = 10e-5
  !call compute_conservative(pos_w,nodes,nx,ny,mx,my)
  !call get_modes_from_nodes(nodes,modes,nx,ny,mx,my)
   !print*,'minmax',maxval(abs(u-modes)),minval(abs(modes-u))
   u = pos_modes
   pAUSE
end subroutine compute_positivity


subroutine positivity_on_faces(modes)
  use parameters_dg_2d
  implicit none
  real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my):: modes, nodes, w_nodes
  real(kind=8),dimension(1:nvar)::u_left,u_right,w_left,w_right
  integer::chsi_right,chsi_left, set_to_zero
  real(kind=8),dimension(1:nvar,1:nx,1:ny):: u_avg,w_avg

  integer::i,j,icell,jcell, ivar, intnode,jntnode
  real(kind=8)::legendre
  real(kind=8)::eps
  eps=10e-10
  chsi_right = 1
  chsi_left = -1
  set_to_zero = 0

   call gl_quadrature(x_quad,w_x_quad,mx)
   call gl_quadrature(y_quad,w_y_quad,my)

  u_avg(:,:,:) = modes(:,:,:,1,1)
  call compute_primitive(u_avg,w_avg,nx,ny,1,1)
  call get_nodes_from_modes(modes,nodes,nx,ny,mx,my)
  call compute_primitive(nodes,w_nodes,nx,ny,mx,my)
  do icell=1,nx
   do jcell = 1,ny

     ! check value at faces in x direction
     do intnode = 1,my
       u_left = 0.0
       u_right = 0.0
       do i=1,mx
         do j=1,my
              u_left = u_left + modes(1:nvar,icell,jcell,i,j)*legendre(chsi_left,i-1)*legendre(y_quad(intnode),j-1)
              u_right = u_right + modes(1:nvar,icell,jcell,i,j)*legendre(chsi_right,i-1)*legendre(y_quad(intnode),j-1)
         end do
       end do
       call compute_primitive(u_left,w_left,1,1,1,1)
        call compute_primitive(u_right,w_right,1,1,1,1)
        !print*,u(1:nvar,icell,jcell,:,:)
        !print*,'u_left',u_left
        !print*,'u_right',u_right
        !print*,'wleft',w_left
        !print*,'wright', w_right
        !print*,i,j,icell,jcell
        if ((w_left(1)<eps).or.(w_left(4)<eps).or.(w_right(1)<eps).or.(w_right(4)<eps)) then
          set_to_zero=1
        end if
    end do

    ! check value at faces in y direction

    do intnode = 1,mx
      u_left = 0.0
      u_right = 0.0
      do i=1,mx
        do j=1,my
          u_left = u_left + modes(1:nvar,icell,jcell,i,j)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)
          u_right = u_right + modes(1:nvar,icell,jcell,i,j)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)
        end do
      end do
      call compute_primitive(u_left,w_left,1,1,1,1)
      call compute_primitive(u_right,w_right,1,1,1,1)
      !print*,u(1:nvar,icell,jcell,:,:)
      !print*,'u_left',u_left
      !print*,'u_right',u_right
      !print*,'wleft',w_left
      !print*,'wright', w_right
      !print*,i,j,icell,jcell

      if ((w_left(1)<eps).or.(w_left(4)<eps))then
            set_to_zero=1
      end if
   end do

   do intnode=1,mx
     do jntnode=1,my
         if ((w_nodes(1,icell,jcell,intnode,jntnode)<eps).or.&
         &(w_nodes(4,icell,jcell,intnode,jntnode)<eps)) then
          set_to_zero = 1
         end if
       end do
     end do

    if (set_to_zero==1) then
      !print*,u_avg(ivar,icell,jcell)
      modes(1,icell,jcell,1,1)  = u_avg(1,icell,jcell)
      modes(1,icell,jcell,2:mx,2:my) = 0.0
      modes(1,icell,jcell,1,2:my) = 0.0
      modes(1,icell,jcell,2:mx,1) = 0.0

      modes(4,icell,jcell,1,1)  = 1e-5/(gamma-1) + 0.5*u_avg(1,icell,jcell)*(u_avg(2,icell,jcell)**2/u_avg(1,icell,jcell) +&
      u_avg(3,icell,jcell)**2/u_avg(1,icell,jcell))
      modes(4,icell,jcell,2:mx,2:my) = 0.0
      modes(4,icell,jcell,1,2:my) = 0.0
      modes(4,icell,jcell,2:mx,1) = 0.0
    end if
    !print*,'set to  zero', set_to_zero
    set_to_zero = 0

 end do
end do

!call get_nodes_from_modes(modes, nodes, nx,ny,mx,my)
!call compute_primitive(nodes, w_nodes, nx,ny,mx,my)
!w_nodes(4,:,:,:,:) = 0.1
!call compute_conservative(w_nodes, nodes, nx,ny,mx,my)
!call get_modes_from_nodes(nodes, modes,nx,ny,mx,my)


end subroutine positivity_on_faces


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

  !do ivar=1,nvar
  !  do icell=1,nx
  !    do jcell=1,ny
  !      u_left=0.0!; u_right=0.0
  !      do i = 1, mx
  !        do j = 1,my
  !          u_left= w_nodes(ivar,icell,jcell,i,j) !u_left+u_lim(1,icell,jcell,i,j)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)
            !u_right= nodes(1,icell,jcell,i,j)  !u_right+u_lim(1,icell,jcell,i,j)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)
  !          if((u_left<1d-6.and.ivar==1).or.(u_left<1d-6.and.ivar==4))then
  !            w_nodes(ivar,icell,jcell,i,j)=10e-5
              !u_lim(1,icell,jcell,2:mx,2:my)=0.0
  !          end if
  !        end do
  !      end do
  !    end do
  !  end do
  !end do


  call compute_conservative(w_nodes,u_lim_nodes,nx,ny,mx,my)
  call get_modes_from_nodes(u_lim_nodes,u_lim,nx,ny,mx,my)
  call limiter_pp(u)
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

function limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,intnode,jntnode)
   use parameters_dg_2d
   implicit none
   ! input
   real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u
   integer::ivar,icell,jcell,itop,ibottom,ileft,iright,intnode,jntnode
   !output
   real(kind=8)::limiting, dx, minmod
   !internal
   real(kind=8)::d_l_x, d_l_y, d_r_x, d_r_y
   real(kind=8)::coeff_i=1.,coeff_j=1.,coeff_u=1.,central_u=1.
   real(kind=8)::minmod2d
   !Write(*,*)icell,jcell,intnode,jntnode
   !coeff_j = sqrt(2.0*dble(jntnode-2)+1.0)/sqrt(2.)*sqrt(2.0*dble(intnode-1)+1.0)/sqrt(2.)
   !coeff_i = sqrt(2.0*dble(intnode-2)+1.0)/sqrt(2.)*sqrt(2.0*dble(jntnode-1)+1.0)/sqrt(2.)

   coeff_j = (2.0*dble(intnode-1)+1.0)*(2*dble(jntnode-1)-1)
   coeff_i = (2.0*dble(jntnode-1)+1.0)*(2*dble(intnode-1)-1)
   !coeff_j = 0.5 !0.5 !sqrt(2.0*dble(jntnode-1)-1.0)/sqrt(2.0*dble(jntnode-1)+1.0)
   !coeff_i = 0.5 !0.5 !sqrt(2.0*dble(intnode-1)-1.0)/sqrt(2.0*dble(intnode-1)+1.0)

   !coeff_j = 1/(2.*sqrt(4.0*dble(jntnode-1)**2-1.0))
   !coeff_i = 1/(2.*sqrt(4.0*dble(intnode-1)**2-1.0))
   coeff_u = (2.0*dble(intnode-1)+1.0)*(2.0*dble(jntnode-1)+1.0)
   central_u = u(ivar,icell,jcell,intnode,jntnode)

   d_r_y = (u(ivar,icell,itop,intnode,jntnode-1)-u(ivar,icell,jcell,intnode,jntnode-1))*coeff_j
   d_l_y = (u(ivar,icell,jcell,intnode,jntnode-1)-u(ivar,icell,ibottom,intnode,jntnode-1))*coeff_j
   d_r_x = (u(ivar,iright,jcell,intnode-1,jntnode)-u(ivar,icell,jcell,intnode-1,jntnode))*coeff_i
   d_l_x = (u(ivar,icell,jcell,intnode-1,jntnode)-u(ivar,ileft,jcell,intnode-1,jntnode))*coeff_i

   limiting = minmod2d(central_u*coeff_u, d_r_y, d_l_y, d_r_x, d_l_x)/coeff_u
   !Write(*,*) icell,jcell,intnode,jntnode,d_r_y,d_l_y,d_r_x,d_l_x,central_u,limiting

   return
end function limiting

subroutine high_order_limiter(u)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, u_new
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::chars, chars_m, u_temp
    real(kind=8),dimension(1:nvar,1:nx,1:ny)::u_avg
    real(kind=8)::limited,limited2,limiting,eps,generalized_minmod
    real(kind=8)::d_l_x, d_l_y, d_r_x, d_r_y
    real(kind=8)::coeff,coeff_u,minmod, coeff_x, coeff_y
    integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
    integer::intnode, jntnode
    integer::done
    real(kind=8),dimension(1:nvar)::u_left,u_right,w_left,w_right, u_top, u_bottom,legendre
    integer::chsi_right,chsi_left, set_to_zero

    if (mx==1.and.my==1) then
      ! no limiting when we only have 1st order approx
      !use_limiter = .false.
      return
    end if
    eps = 10e-10
    u_new = u

    do ivar = 1,nvar
      do icell=1,nx
        do jcell = 1,ny
          done = 0
          ileft = icell - 1
          iright = icell + 1
          itop = jcell + 1
          ibottom = jcell - 1
          call get_boundary_conditions(ileft,1)
          call get_boundary_conditions(iright,1)
          call get_boundary_conditions(itop,2)
          call get_boundary_conditions(ibottom,2)
          do intnode = mx,2,-1
            !print*,intnode
            limited = limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,intnode,intnode)
             if(limited .ne. u(ivar,icell,jcell,intnode,intnode))then
                 u_new(ivar,icell,jcell,intnode,intnode)=limited
             else
                 exit
             end if
            ! u_new(ivar,icell,jcell,intnode,intnode)= limited
            do jntnode = intnode-1,2,-1
               limited = limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,intnode,jntnode)
               limited2 = limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,jntnode,intnode)
               if((abs(limited - u(ivar,icell,jcell,intnode,jntnode))<eps) &
               &.and. (abs(limited2 - u(ivar,icell,jcell,jntnode,intnode))<eps)) then
                 done = 1
                 exit
               else
                 u_new(ivar,icell,jcell,intnode,jntnode)=limited
                 u_new(ivar,icell,jcell,jntnode,intnode)=limited2
               end if
             end do
            !
             if (done .eq. 1) then
                   exit
             end if

            !coeff_x = (2*dble(jntnode-2)+1)
            coeff_y = (2*dble(intnode-1)+1)  !sqrt(2.0*dble(intnode-1)-1.0)*sqrt(2.0*dble(intnode-1)+1.0)
            coeff_u = (2*dble(intnode-1)+1)

            d_r_y = u(ivar,icell,itop,intnode-1,1) - &
                & u(ivar,icell,jcell,intnode-1,1)

            d_l_y = u(ivar,icell,jcell,intnode-1,1) - &
                & u(ivar,icell,ibottom,intnode-1,1)

            d_r_x = u(ivar,iright,jcell,1,intnode-1) - &
                & u(ivar,icell,jcell,1,intnode-1)

            d_l_x = u(ivar,icell,jcell,1,intnode-1) - &
                & u(ivar,ileft,jcell,1,intnode-1)

            limited = generalized_minmod(u(ivar,icell,jcell,1,intnode)*coeff_u,d_r_y*coeff_y,d_l_y*coeff_y)/coeff_u
            limited2 = generalized_minmod(u(ivar,icell,jcell,intnode,1)*coeff_u,d_r_x*coeff_y,d_l_x*coeff_y)/coeff_u
            !Write(*,*) icell,jcell,intnode,u(ivar,icell,jcell,1,intnode)*coeff_u,d_r_x*coeff_y,d_l_x*coeff_y,limited2
            if((limited == u(ivar,icell,jcell,1,intnode)) .and. (limited2 == u(ivar,icell,jcell,intnode,1))) then
              !print*,'no limiting'
              exit
            else
              u_new(ivar,icell,jcell,1,intnode)=limited
              u_new(ivar,icell,jcell,intnode,1)=limited2
              !print*,'limiting'
            end if
          end do
        end do
     end do
   end do

  u = u_new
  u_avg = u(1:nvar,:,:,1,1)
  print*,'min mean average on HIO',minval(u_avg(1,:,:))
  call compute_positivity(u)
  !call limiter_pp(u)
  !write(*,*) u(1,1,1,:,:)
  !pausel

  !call positivity_on_faces(u)
  !call compute_positivity(u)
  end subroutine high_order_limiter



  subroutine limiter_positivity_2(u)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, w_lim,u_lim
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w,w_nodes, chars,chars_m,u_temp
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

    call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
    call compute_characteristics(u,chars,nx,ny,mx,my)
    call compute_primitive(nodes,w_nodes,nx,ny,mx,my)
    call get_modes_from_nodes(w_nodes,w,nx,ny,mx,my)
    call get_modes_from_nodes(chars,chars_m,nx,ny,mx,my)

    !print*,'r',u(1,:,:,:,:)
    !print*,'p',u(4,:,:,:,:)
    !STOP
    u_temp = u
    u_lim = chars_m
    u = chars_m
    do ivar = 1,nvar
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
          u_deriv = u(ivar,icell,jcell, 2, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1+1)+1.0))
          u_lim(ivar,icell,jcell, 2, 1) = minmod(u_deriv,&
          &(u_center-u_left),(u_right-u_center))

          if(ABS(u_lim(ivar,icell,jcell, 2, 1)-u_deriv).GT.0.01*ABS(u_deriv)) then
            u_lim(ivar,icell,jcell,2:mx,1) = 0.0
            u_lim(ivar,icell,jcell,mx,mx) = 0.0
          end if
        end do
      end do
    end do
    ! end do

    do ivar = 1,nvar
      do icell = 1,nx
        do jcell=1,nx
          !do intnode = mx-1,1,-1
          left = icell -1
          right = icell+1
          if (icell == 1) then
            left = nx
          else if (icell == nx) then
            right = 1
          end if
          !if (jcell == 1) then
          ! bottom = nx
          !else if (jcell == nx) then
          ! top = 1
          !end if

          u_left = u(ivar,jcell,left,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
          u_right = u(ivar,jcell,right,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
          u_center = u(ivar,jcell,icell,1, 1)!*sqrt(dble(2.0))/sqrt((2.0*dble(1)+1.0))
          u_deriv = u(ivar,jcell,icell, 1, 2) !*sqrt(dble(2.0))/sqrt((2.0*dble(1+1)+1.0))
          u_lim(ivar,jcell, icell, 1, 2) =  minmod(u_deriv,&
          &(u_center-u_left),(u_right-u_center))
          if(ABS(u_lim(ivar,jcell,icell, 1,2)-u_deriv).GT.0.01*ABS(u_deriv)) then
            u_lim(ivar,jcell,icell, 1,2:mx) = 0.0
            u_lim(ivar,jcell,icell, mx,mx) = 0.0
          end if
        end do
      end do
    end do

    call get_nodes_from_modes(u_lim,chars,nx,ny,mx,my)
    call compute_cons_from_characteristics(chars,u_temp,nodes_cons,nx,ny,mx,my)
    call compute_primitive(nodes_cons,nodes,nx,ny,mx,my)

    !call get_nodes_from_modes(u_lim,nodes,nx,ny,mx,my)

    do ivar=1,nvar
      do icell=1,nx
        do jcell=1,ny
          u_left=0.0; u_right=0.0
          !Loop over modes
          !do intnode=1,mx
          do i = 1, mx
            do j = 1,my
              u_left= nodes(ivar,icell,jcell,i,j) !u_left+u_lim(1,icell,jcell,i,j)*legendre(chsi_left,j-1)*legendre(x_quad(intnode),i-1)
              !u_right= nodes(1,icell,jcell,i,j)  !u_right+u_lim(1,icell,jcell,i,j)*legendre(chsi_right,j-1)*legendre(x_quad(intnode),i-1)
              if((u_left<1d-10.and.ivar==1).or.(u_left<1d-10.and.ivar==4))then
                nodes(ivar,icell,jcell,i,j)=1e-5
                !u_lim(1,icell,jcell,2:mx,2:my)=0.0
              end if
            end do
          end do
        end do
      end do
    end do

    call compute_conservative(nodes,nodes_cons,nx,ny,mx,my)
    call get_modes_from_nodes(nodes_cons,u_lim,nx,ny,mx,my)

    u = u_lim
  end subroutine limiter_positivity_2


  subroutine limiter_pp(u)
    ! eats u, conservative modes and returns limited conservative modes
    ! performs limiting by looking at primitives (as described on paper)
    use parameters_dg_2d
    implicit none
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, u_lim, u_lim_nodes
    real(kind=8),dimension(1:nvar,1:nx,1:ny)::w_average, u_avg

    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
    real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u_nodes, w_nodes

    real(kind=8)::big_theta, little_theta, theta, w_average_stencil
    real(kind=8),dimension(1:nvar)::maximum_neigh_value_1, minimum_neigh_value_1

    integer,dimension(1:9)::indices_i,indices_j
    real(kind=8)::u_left, eps, max_temp, min_temp, vs
    real(kind=8)::limited1,limited2, generalized_minmod, minmod, legendre, dx, ratio, phi_lim
    integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
    integer::intnode, jntnode, indx
    integer::left,right,bottom,top, component

    dx = boxlen_x/dble(nx)
    eps = 10e-10

    ! transform to nodes
    !call get_nodes_from_modes(u, u_nodes, nx, ny, mx, my)
    !call compute_primitive(u_nodes,w_nodes,nx,ny,mx,my)
    !call get_modes_from_nodes(w_nodes, w, nx, ny, mx, my)
    !u = w
    u_lim = u
    ! get averages
    do ivar=1,nvar
      do icell = 1,nx
        do jcell = 1,ny
          u_avg(ivar,icell,jcell) = u(ivar,icell,jcell,1,1)!sum(w_nodes(ivar,icell,jcell,:,:))/dble(mx*my)
        end do
      end do
    end do
    print*,'value at entrance'
    print*,maxval(u_avg(1,:,:)),minval(u_avg(1,:,:))
    print*,maxval(u_avg(4,:,:)),minval(u_avg(4,:,:))
    print*,maxval(u_avg(2,:,:)),minval(u_avg(2,:,:))
    print*,maxval(u_avg(3,:,:)),minval(u_avg(3,:,:))


    do icell = 1,nx
      do jcell = 1,ny
        do ivar = 1,nvar

          if (icell == 1) then
            indices_i = (/ nx, icell, icell+1, &
                        &  nx, icell, icell+1, &
                        &  nx, icell, icell+1 /)
          else if (icell == nx) then
            indices_i = (/ icell - 1, icell, 1, &
                        &  icell - 1, icell, 1, &
                        &  icell - 1, icell, 1 /)
          else
            indices_i = (/ icell - 1, icell, icell+1, &
                        &  icell - 1, icell, icell+1, &
                        &  icell - 1, icell, icell+1 /)
          end if

          if (jcell == 1) then
            indices_j = (/ ny , ny , ny, &
                         & jcell,      jcell,         jcell,&
                         & jcell - 1 , jcell - 1 , jcell-1/)
          else if (jcell == ny) then
            indices_j = (/ jcell + 1 , jcell + 1 , jcell+1, &
                         & jcell,      jcell,         jcell,&
                         & 1 , 1 , 1/)
          else
            indices_j = (/ jcell + 1 , jcell + 1 , jcell+1, &
                         & jcell,      jcell,         jcell,&
                         & jcell - 1 , jcell - 1 , jcell-1/)
          end if

          max_temp = 0
          min_temp = 0
          do indx = 1,9
            if (max_temp < u_avg(ivar,indices_i(indx),indices_j(indx))-u_avg(ivar,icell,jcell)) then
              max_temp = u_avg(ivar,indices_i(indx),indices_j(indx))-u_avg(ivar,icell,jcell)
              !print*, max_temp
            end if
            if (min_temp > u_avg(ivar,indices_i(indx),indices_j(indx))-u_avg(ivar,icell,jcell)) then
              min_temp = u_avg(ivar,indices_i(indx),indices_j(indx))-u_avg(ivar,icell,jcell)
              !print*, min_temp
            end if
          end do

          if ((abs(u(ivar,icell,jcell,1,2))+abs(u(ivar,icell,jcell,2,1)))<eps) then
            vs = min(abs(max(max_temp,eps)),abs(min(eps,min_temp)))/(eps)
          else
            vs = min(abs(max(max_temp,eps)),abs(min(eps,min_temp)))/(abs(u(ivar,icell,jcell,1,2))+abs(u(ivar,icell,jcell,2,1)))
          end if
          !print*,vs
          !print*,'min',min(vs,1.)
          !vs = 0.5
          vs = min(vs,1.)
          vs = max(vs,0.)
          vs = 0.0
          u_lim(ivar,icell,jcell,1,2) = min(vs,1.)*u(ivar,icell,jcell,1,2)
          u_lim(ivar,icell,jcell,2,1) = min(vs,1.)*u(ivar,icell,jcell,2,1)
          !if (u_lim(1,icell,jcell,1,1)<0.0) then
          !  u_lim(1,icell,jcell,1,1) = -u_lim(1,icell,jcell,1,1)
          !end if
          !if (u(4,icell,jcell,1,1)<eps) then
          !  u_lim(4,icell,jcell,1,1)=10e-5
          !  u_lim(4,icell,jcell,2,1) = 0.0
          !  u_lim(4,icell,jcell,1,2) = 0.0
          !u_lim(ivar,icell,jcell,2,2) = 0.0 !0.0
            !print*,'press neg',vs
          !end if
          !print*,vs
        end do
      end do
    end do
    u_avg = u_lim(1:ivar,:,:,1,1)
    print*,'lol',maxval(abs(u_lim-u)),minval(abs(u_lim-u))
    print*,maxval(u_avg(1,:,:)),minval(u_avg(1,:,:))
    print*,maxval(u_avg(4,:,:)),minval(u_avg(4,:,:))
    print*,maxval(u_avg(2,:,:)),minval(u_avg(2,:,:))
    print*,maxval(u_avg(3,:,:)),minval(u_avg(3,:,:))
    !u = u_lim
    !call get_nodes_from_modes(u_lim,w_nodes,nx,ny,mx,my)
    !call compute_conservative(w_nodes,u_nodes,nx,ny,mx,my)
    !call get_modes_from_nodes(u_nodes,u_lim,nx,ny,mx,my)
    u = u_lim
    !pause
  end subroutine limiter_pp



  subroutine third_order_limiter(u)
      use parameters_dg_2d
      implicit none
      real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::u, u_new
      real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::w
      real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::nodes, modes
      real(kind=8),dimension(1:nvar,1:nx,1:ny,1:mx,1:my)::chars, chars_m, u_temp
      real(kind=8),dimension(1:nvar,1:nx,1:ny)::u_avg
      real(kind=8)::limited,limited2,limiting
      real(kind=8)::d_l_x, d_l_y, d_r_x, d_r_y
      real(kind=8)::coeff,coeff_u,minmod
      integer::i,j,icell,jcell, ivar, itop,ibottom,ileft,iright
      integer::intnode, jntnode
      integer::done

      if (mx==1.and.my==1) then
        ! no limiting when we only have 1st order approx
        !use_limiter = .false.
        return
      end if
      ! call compute_primitive(u,w,nx,ny,mx,my)
      ! TODO: implement characteristic variable representation
      ! using Roe average
      !write(*,*) u(1,1,1,:,:)

      !call get_nodes_from_modes(u,nodes,nx,ny,mx,my)
      !call compute_characteristics(nodes,chars,nx,ny,mx,my)
      !call get_modes_from_nodes(chars,chars_m,nx,ny,mx,my)
      u_new = u
      !u_temp = u
      !u = chars_m
      !u_new = u
      print*,maxval(u)
      !pause
      do ivar = 1,nvar
        do icell=1,nx
          do jcell = 1,ny
            done = 0
            ileft = icell - 1
            iright = icell + 1
            itop = jcell + 1
            ibottom = jcell - 1
            call get_boundary_conditions(ileft,1)
            call get_boundary_conditions(iright,1)
            call get_boundary_conditions(itop,2)
            call get_boundary_conditions(ibottom,2)
            do intnode = mx,2,-1

              limited = limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,intnode,intnode)
              if(limited .ne. u(ivar,icell,jcell,intnode,intnode))then
                  u_new(ivar,icell,jcell,intnode,intnode)=limited
              else
                  exit
              end if

              do jntnode = intnode-1,2,-1
                limited = limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,intnode,jntnode)
                limited2 = limiting(u,ivar,icell,jcell,itop,ibottom,ileft,iright,jntnode,intnode)
                if((limited .ne. u(ivar,icell,jcell,intnode,jntnode)) .or. (limited2 .ne. u(ivar,icell,jcell,jntnode,intnode))) then
                  u_new(ivar,icell,jcell,intnode,jntnode)=limited
                  u_new(ivar,icell,jcell,jntnode,intnode)=limited2
                else
                  done = 1
                  exit
                end if
              end do

              if (done .eq. 1) then
                    exit
              end if
              !coeff = sqrt(2.0*dble(intnode-2)+1.0)/2.
              coeff = 0.5 !sqrt(2.0*dble(intnode-1)-1.0)/sqrt(2.0*dble(intnode-1)+1.0)
              coeff_u = 1. !sqrt(2.0*dble(intnode-1)+1.0)
              !coeff = 1/(2.*sqrt(4*dble(intnode-1)**2-1))
              d_r_y = u(ivar,icell,itop,1,intnode-1) - &
                  & u(ivar,icell,jcell,1,intnode-1)

              d_l_y = u(ivar,icell,jcell,1,intnode-1) - &
                  & u(ivar,icell,ibottom,1,intnode-1)

              d_r_x = u(ivar,iright,jcell,intnode-1,1) - &
                  & u(ivar,icell,jcell,intnode-1,1)

              d_l_x = u(ivar,icell,jcell,intnode-1,1) - &
                  & u(ivar,ileft,jcell,intnode-1,1)

              limited = minmod(coeff_u*u(ivar,icell,jcell,1,intnode),coeff*d_r_y,coeff*d_l_y)/coeff_u
              limited2 = minmod(coeff_u*u(ivar,icell,jcell,intnode,1),coeff*d_r_x,coeff*d_l_x)/coeff_u
              if((limited .ne. u(ivar,icell,jcell,1,intnode)) .or. (limited2 .ne. u(ivar,icell,jcell,intnode,1))) then
                 u_new(ivar,icell,jcell,1,intnode)=limited
                 u_new(ivar,icell,jcell,intnode,1)=limited2
              else
                 exit
              end if
            end do
          end do
       end do
     end do
    ! Update variables with limited states

    !call get_nodes_from_modes(u_new,chars,nx,ny,mx,my)
    !print*,'min mean average on HIO',minval(chars)
    !call compute_cons_from_characteristics(chars,u_temp,w,nx,ny,mx,my)
    !print*,'min mean average on HIO',minval(w)
    !call get_modes_from_nodes(w,u_new,nx,ny,mx,my)

    u = u_new
    u_avg = u(1:nvar,:,:,1,1)
    print*,'min mean average on HIO',minval(u_avg(1,:,:))
    !call compute_positivity(u)
    call limiter_pp(u)
    !write(*,*) u(1,1,1,:,:)
    !pause
  end subroutine third_order_limiter
