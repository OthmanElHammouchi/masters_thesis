module interface_mod
  use, intrinsic :: iso_c_binding

  implicit none

  interface

    function init_rng(n_threads, seed) result(res) bind(c)
      import
      integer(c_int), value :: seed
      integer(c_int), value :: n_threads
      type(c_ptr) :: res
    end function init_rng

    function rnorm_par(rng, i_thread, mean, sd) result(sample) bind(c)
      import
      type(c_ptr), value :: rng
      integer(c_int), value :: i_thread
      real(c_double), value :: mean
      real(c_double), value :: sd
      real(c_double) :: sample
    end function rnorm_par

    function rgamma_par(rng, i_thread, shape, scale) result(sample) bind(c)
      import
      type(c_ptr), value :: rng
      integer(c_int), value :: i_thread
      real(c_double), value :: shape
      real(c_double), value :: scale
      real(c_double) :: sample
    end function rgamma_par

    function rpois_par(rng, i_thread, lambda) result(sample) bind(c)
      import
      type(c_ptr), value :: rng
      integer(c_int), value :: i_thread
      real(c_double), value :: lambda
      integer(c_int) :: sample
    end function rpois_par

    function runif_par(rng, i_thread) result(sample) bind(c)
      import
      type(c_ptr), value :: rng
      integer(c_int), value :: i_thread
      real(c_double) :: sample
    end function runif_par

    function pgbar_create(total, freq) result(val) bind(c)
      import
      integer(c_int), value :: total
      integer(c_int), value :: freq
      type(c_ptr) :: val
    end function pgbar_create

    subroutine pgbar_incr(progress_bar) bind(c)
      import
      type(c_ptr), value :: progress_bar
    end subroutine pgbar_incr

    subroutine check_user_input() bind(c)
    end subroutine check_user_input

    subroutine rprint_par(str) bind(c)
      import
      character(kind=c_char), intent(in) :: str(*)
    end subroutine rprint_par

    function raise(sig) bind(C, name="raise")
      use iso_c_binding, only: c_int
      integer(c_int) :: raise
      integer(c_int), value :: sig
    end function

  end interface

end module interface_mod
