module dg_commons

  ! DG solver parameters
  integer,parameter::n=4
  character(LEN=3),parameter::integrator='RK4'
  integer,parameter::nx=300
  integer,parameter::nquad=n
  integer,parameter::nvar=3
  integer,parameter::riemann=1
  logical,parameter::use_limiter=.true.

  ! Problem set-up
  integer,parameter::ninit=8
  integer,parameter::bc=4
  integer,parameter::source=2
  real(kind=8)::tend=0.10
  real(kind=8)::boxlen=10.0
  real(kind=8)::gamma=1.4
  real(kind=8)::pert=0.01

  ! Misc commons
  real(kind=8),dimension(1:nquad)::chsi_quad,w_quad

end module dg_commons
