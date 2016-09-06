module parameters_dg_2d
  ! solver parameter
  integer,parameter::nx=30
  integer,parameter::ny=30
  integer,parameter::mx=2
  integer,parameter::my=2

  integer,parameter::nvar=4
  integer,parameter::riemann=2
  logical,parameter::use_limiter=.false.
  logical,parameter::make_movie=.true.
  character(LEN=3),parameter::solver='RK4' !or EQL to use the equilibrium solution ^^
  character(LEN=3),parameter::limiter_type='ROS' !or EQL to use the equilibrium solution ^^
  ! Problem set-up
  integer,parameter::ninit=1
  integer,parameter::bc=1
  integer,parameter::nequilibrium=3
  integer,parameter::source=1

  real(kind=8)::tend=1.0
  real(kind=8)::boxlen_x=1.0
  real(kind=8)::boxlen_y=1.0
  real(kind=8)::gamma=1.4
  real(kind=8)::cfl = 0.5

  real(kind=8)::eta=0.1

  ! misc commons
  real(kind=8),dimension(1:mx)::x_quad, w_x_quad
  real(kind=8),dimension(1:mx)::y_quad, w_y_quad
  real(kind=8)::M=0
  integer,parameter::interval = 50
end module parameters_dg_2d
