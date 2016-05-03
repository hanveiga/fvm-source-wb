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