interface

! Interfaces to C functions wrapping R's RNG functionality, which mitigates
! the need to use seperate libraries.

subroutine GetRNGstate() bind(C, name="GetRNGstate_fentry_")

end subroutine GetRNGstate

subroutine PutRNGstate() bind(C, name="PutRNGstate_fentry_")

end subroutine PutRNGstate

real(c_double) function rnorm(mean, sd) bind(C, name="rnorm_fentry_")

   import c_double

   real(c_double), intent(in), value :: mean, sd

end function rnorm

real(c_double) function rgamma(shape, scale) bind(C, name="rgamma_fentry_")

   import c_double

   real(c_double), intent(in), value :: shape, scale

end function rgamma

real(c_double) function rpois(lambda) bind(C, name="rpois_fentry_")

   import c_double

   real(c_double), intent(in), value :: lambda

end function rpois

end interface