module mod_global

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_int), parameter :: NONE = 0

  ! Process distributions
  integer(c_int), parameter :: NORMAL = 1
  integer(c_int), parameter :: GAMMA = 2
  integer(c_int), parameter :: POISSON = 3

  ! Bootstrap types
  integer(c_int), parameter :: PARAM = 1
  integer(c_int), parameter :: RESID = 2
  integer(c_int), parameter :: PAIRS = 3

  ! Residual types
  integer(c_int), parameter :: STANDARDISED = 1
  integer(c_int), parameter :: STUDENTISED = 2
  integer(c_int), parameter :: LOGNORMAL = 3

  ! Simulation types
  integer(c_int), parameter :: SINGLE = 1
  integer(c_int), parameter :: CALENDAR = 2
  integer(c_int), parameter :: ORIGIN = 3

  integer(c_int), parameter :: SUCCESS = 0
  integer(c_int), parameter :: FAILURE = 1

  type(c_ptr) :: rng
  integer(c_int) :: i_thread
  integer(c_int), save :: n_threads
  logical(c_bool), save :: first_call = .true.
  !$omp threadprivate(rng, i_thread)
  
end module mod_global
