module parameters
  ! solver parameters
  integer,parameter::nx=1000
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.true.
  character(LEN=3),parameter::solver='EQL' !or EQL to use the equilibrium solution ^^

  ! Problem set-up
  integer,parameter::ninit=3
  integer,parameter::bc=3
  integer,parameter::nequilibrium=3

  real(kind=8)::tend=1
  real(kind=8)::boxlen=1.0
  real(kind=8)::gamma=1.4

end module parameters
