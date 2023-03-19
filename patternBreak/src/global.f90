module global

  use, intrinsic :: iso_c_binding

  implicit none

  integer(c_int), parameter :: NORMAL = 1
  integer(c_int), parameter :: GAMMA = 2

  integer(c_int), parameter :: CONDITIONAL = 1
  integer(c_int), parameter :: UNCONDITIONAL = 2
  integer(c_int), parameter :: PAIRS = 3

  integer(c_int), parameter :: NORMAL_STANDARDISED = 1
  integer(c_int), parameter :: NORMAL_MODIFIED = 2
  integer(c_int), parameter :: NORMAL_STUDENTISED = 3
  integer(c_int), parameter :: LOGNORMAL = 4
  integer(c_int), parameter :: PARAMETRIC = 5

  integer(c_int), parameter :: PARAM_RESAMPLE = 1
  integer(c_int), parameter :: NON_PARAM_RESAMPLE = 2
  integer(c_int), parameter :: PAIRS_RESAMPLE = 3

  integer(c_int), parameter :: SINGLE = 1
  integer(c_int), parameter :: CALENDAR = 2
  integer(c_int), parameter :: ORIGIN = 3

  integer(c_int), parameter :: SUCCESS = 0
  integer(c_int), parameter :: FAILURE = 1

  logical(c_bool), parameter :: TRUE = .true.
  logical(c_bool), parameter :: FALSE = .false.

  type(c_ptr) :: rng
  integer(c_int) :: i_thread

  !$omp threadprivate(rng, i_thread)
  
end module global
