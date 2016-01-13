module fvm_commons

  ! DG solver parameters
  integer,parameter::n=3
  character(LEN=3),parameter::integrator='RK4'
  integer,parameter::nx=200
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.true.

  ! Problem set-up
  integer,parameter::ninit=4
  integer,parameter::bc=2
  integer,parameter::source=2
  real(kind=8)::tend=0.2
  real(kind=8)::boxlen=1
  real(kind=8)::gamma=1.4

  ! Misc commons

end module fvm_commons
