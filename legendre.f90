function legendre(x,n)
  integer::n
  real(kind=8)::x
  real(kind=8)::legendre
  x=min(max(x,-1.0),1.0)
  select case(n)
  case(0)
     legendre=1.0
  case(1)
     legendre=x
  case(2)
     legendre=0.5*(3.0*x**2-1.0)
  case(3)
     legendre=0.5*(5.0*x**3-3.0*x)
  case(4)
     legendre=0.125*(35.0*x**4-30.0*x**2+3.0)
  case(5)
     legendre=0.125*(63.0*x**5-70.0*x**3+15.0*x)
  case(6)
     legendre=1.0/16.0*(231.0*x**6-315.0*x**4+105.0*x**2-5.0)
  end select
  legendre=sqrt((2.0*dble(n)+1.0))*legendre
  return
end function legendre

function legendre_prime(x,n)
  integer::n
  real(kind=8)::x
  real(kind=8)::legendre_prime
  x=min(max(x,-1.0),1.0)
  select case(n)
  case(0)
     legendre_prime=0.0
  case(1)
     legendre_prime=1.0
  case(2)
     legendre_prime=3.0*x
  case(3)
     legendre_prime=0.5*(15.0*x**2-3.0)
  case(4)
     legendre_prime=0.125*(140.0*x**3-60.0*x)
  case(5)
     legendre_prime=0.125*(315.0*x**4-210.0*x**2+15.0)
  case(6)
     legendre_prime=1.0/16.0*(1386.0*x**5-1260.0*x**3+210.0*x)
  end select
  legendre_prime=sqrt((2.0*dble(n)+1.0))*legendre_prime
  return
end function legendre_prime

function legendre_prime_prime(x,n)
  integer::n
  real(kind=8)::x
  real(kind=8)::legendre_prime
  x=min(max(x,-1.0),1.0)
  select case(n)
  case(0)
     legendre_prime=0.0
  case(1)
     legendre_prime=0.0
  case(2)
     legendre_prime=3.0
  case(3)
     legendre_prime=0.5*(30.0*x)
  case(4)
     legendre_prime=0.125*(280.0*x**2-60.0)
  case(5)
     legendre_prime=0.125*(4*315.0*x**3-420.0*x**1)
  case(6)
     legendre_prime=1.0/16.0*(1386.0*5*x**4-1260.0*3*x**2+210.0)
  end select
  legendre_prime=sqrt(2.0*dble(n)+1.0)*legendre_prime
  return
end function legendre_prime_prime

subroutine gl_quadrature(x_quad,w_quad,n)
  integer::n
  real(kind=8),dimension(1:n)::x_quad,w_quad

  integer::i,iter
  real(kind=8)::dpi=acos(-1.0d0),xx
  real(kind=8)::legendre,legendre_prime

  !write(*,*)"Computing Gauss-Legendre quadrature points and weights."

  if (n==1) THEN
    x_quad(1) = 0.0
    w_quad(1) = 2.0
    return
  end if

  do i=1,n
     xx=(1.0-0.125/n/n+0.125/n/n/n)* & 
          & cos(dpi*(4.0*dble(i)-1.0)/(4.0*dble(n)+2.0))
     do iter=1,50
        xx=xx-legendre(xx,n)/legendre_prime(xx,n)
     end do
     x_quad(i)=xx
     w_quad(i)=2.0*(2.0*dble(n)+1.0)/(1.0-x_quad(i)**2) &
          & /legendre_prime(x_quad(i),n)**2
  end do
  do i=n/2+1,n
     x_quad(i)=-x_quad(n-i+1)
     w_quad(i)=w_quad(n-i+1)
  end do
!!$  do i=1,n
!!$     write(*,*)i,x_quad(i),w_quad(i)
!!$  end do

end subroutine gl_quadrature


!subroutine gll_quadrature(x_quad,w_quad,n)
!  integer::n
!  real(kind=8),dimension(1:n)::x_quad,w_quad
!  real(kind=8), dimension(1:(n-2))::x_quad_inner, w_quad_inner
!  integer::i,iter
!  real(kind=8)::dpi=acos(-1.0d0),xx
!  real(kind=8)::legendre,legendre_prime

!  x_quad(1) = -1
!  w_quad(1) = 2/(dble(n)*dble(n-1))
!  x_quad(n) = 1
!  w_quad(n) = 2/(dble(n)*dble(n-1))

!  do i = 2,n-1
!    xx = (1-3*(n-2)/(8*(n-1)**3))* &
!        & cos(dpi*(4.0*dble(i)-3)/(4*(n-1)+1))
!    for iter=1,50
!        xx = x - 2*legendre_prime()*legendre_prime_prime()/&
!        & (2*legendre_prime_prime()**2-lendre
!  write(*,*)'quad=',x_quad_inner,w_quad_inner

!end subroutine gll_quadrature

real(kind=8) function lagrange_poly(x_points, y_points, x, n)
  implicit none
  integer::n
  real(kind=8),dimension(1:n)::x_points
  real(kind=8),dimension(1:n)::y_points
  real(kind=8)::x
  real(kind=8)::L, out
  real(kind=8)::lagrange_basis
  integer::i

  L = 0.0
  do i=1,n
    L = L + dble(y_points(i))*lagrange_basis(x_points, x_points(i), x, n)
  end do
  lagrange_poly = L
end 

real(kind=8) function lagrange_basis(points, base_pt, x, n) 
  implicit none
  integer::n
  real(kind=8),dimension(1:n)::points
  real(kind=8)::x, base_pt
  real(kind=8)::l, temp
  integer::i

  l = 1.0
  temp = 1.0
  DO i=1,n
    !write(*,*),points(i),base_pt,x
    IF (ABS(points(i)-base_pt)<0.00001) THEN
      temp = temp
      !write(*,*),'skip number'
    ELSE
      l = temp*(x-dble(points(i)))/(dble(base_pt)-dble(points(i)))
      temp = l
      !write(*,*),'l=', l
    END IF
  END DO
  lagrange_basis = l
  end

real(kind=8) function lagrange_prime(x_points, y_points, x, n)
  implicit none
  real(kind=8),dimension(1:n)::x_points
  real(kind=8),dimension(1:n)::y_points
  real(kind=8)::x
  real(kind=8)::L
  real(kind=8)::lagrange_prime_basis
  integer::n,i

  L = 0.0
  do i=1,n
    L = L + y_points(i)*lagrange_prime_basis(x_points, x_points(i), x, n)
  end do 
  lagrange_prime = L
end

real(kind=8) function lagrange_prime_basis(points, base_pt, x, n)
  implicit none
  real(kind=8),dimension(1:n)::points, sumpoints
  real(kind=8)::x, base_pt, sumpoint, xpoint
  real(kind=8)::partial, summed
  integer::n,i,j

  summed = 0.0
  partial = 1.0

  do i = 1,n
    partial = 1.0
    IF (ABS(points(i)-base_pt) <0.00001) THEN
      partial = partial
      cycle
    ELSE
      partial = dble(1)/(base_pt - points(i))
    END IF
    do j = 1,n
      IF (ABS(points(j)-base_pt) < 0.00001) THEN
      ELSE IF (ABS(points(i)-points(j)) < 0.00001) THEN
      ELSE
        partial = partial*(x-points(j))/(base_pt-points(j))
      END IF
      end do
    summed = summed + partial
    end do

  lagrange_prime_basis = summed
end function 
