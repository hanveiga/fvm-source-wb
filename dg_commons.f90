module dg_commons

  ! DG solver parameters
  integer,parameter::n=2
  character(LEN=3),parameter::integrator='RKe'
  integer,parameter::nx=128
  integer,parameter::nquad=n
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.true.

  ! Problem set-up
  integer,parameter::ninit=7
  integer,parameter::bc=4
  integer,parameter::source=2
  real(kind=8)::tend=100
  real(kind=8)::boxlen=1.0
  real(kind=8)::gamma=1.4
  real(kind=8)::pert=0.0000001

  ! Misc commons
  real(kind=8),dimension(1:nquad)::chsi_quad,w_quad

end module dg_commons
