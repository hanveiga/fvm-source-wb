module parameters_dg_2d
  ! solver parameters
  integer,parameter::nx=25
  integer,parameter::ny=25
  integer,parameter::mx=3
  integer,parameter::my=3
  
  integer,parameter::nvar=4
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.true.
  character(LEN=3),parameter::solver='EQL' !or EQL to use the equilibrium solution ^^

  ! Problem set-up
  integer,parameter::ninit=1
  integer,parameter::bc=3
  integer,parameter::nequilibrium=1
  integer,parameter::source=2

  real(kind=8)::tend=0.1
  real(kind=8)::boxlen_x=1.0
  real(kind=8)::boxlen_y=1.0
  real(kind=8)::gamma=1.4
  real(kind=8)::cfl = 0.5

  real(kind=8)::eta=0.00001
  
  ! misc commons
  real(kind=8),dimension(1:mx)::x_quad, w_x_quad
  real(kind=8),dimension(1:mx)::y_quad, w_y_quad
end module parameters_dg_2d
