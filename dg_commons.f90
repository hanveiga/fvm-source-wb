module dg_commons

  ! DG solver parameters
  integer,parameter::n=3
  character(LEN=3),parameter::integrator='RKi'
  integer,parameter::nx=128
  integer,parameter::nquad=n
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.false.

  ! Problem set-up
  integer,parameter::ninit=8
  integer,parameter::bc=5
  integer,parameter::source=2
  real(kind=8)::tend=0.2
  real(kind=8)::boxlen=1.0
  real(kind=8)::gamma=1.4
  real(kind=8)::pert=1e-8

  ! Misc commons
  real(kind=8),dimension(1:nquad)::chsi_quad,w_quad
  character(LEN=7)::finaloutput='dg31e-8'

end module dg_commons
