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


