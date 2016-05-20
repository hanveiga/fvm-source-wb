module parameters
  ! solver parameters
  integer,parameter::nx=128   
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.true.
  logical,parameter::make_movie=.false.
  character(LEN=3),parameter::solver='WB1' !or EQL to use the equilibrium solution ^^

  ! Problem set-up
  integer,parameter::ninit=2
  integer,parameter::bc=2
  integer,parameter::nequilibrium=2

  real(kind=8)::tend=0.20
  real(kind=8)::boxlen=1.0
  real(kind=8)::gamma=1.4
  real(kind=8)::eta = 1e-8

  character(LEN=7)::finaloutput='wb11e-8'

end module parameters
