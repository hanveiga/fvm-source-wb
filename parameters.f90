module parameters
  ! DG solver parameters
  integer,parameter::nx=200
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.true.

  ! Problem set-up
  integer,parameter::ninit=1
  integer,parameter::bc=2
  integer,parameter::nequilibrium=1

  real(kind=8)::tend=5
  real(kind=8)::boxlen=2.0
  real(kind=8)::gamma=1.4

end module parameters
