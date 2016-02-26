module dg_commons

  ! DG solver parameters
  integer,parameter::n=4
  character(LEN=3),parameter::integrator='RK4'
  integer,parameter::nx=300
  integer,parameter::nquad=n
  integer,parameter::nvar=3
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.false.

  ! Problem set-up
  integer,parameter::ninit=7
  integer,parameter::bc=2
  integer,parameter::source=2
  real(kind=8)::tend=0.25
  real(kind=8)::boxlen=1.0
  real(kind=8)::gamma=1.4
  real(kind=8)::pert=0.01

  ! Misc commons
  real(kind=8),dimension(1:nquad)::chsi_quad,w_quad

end module dg_commons
