module interface

   use, intrinsic :: iso_c_binding

   implicit none

   interface
      
      function init_rng(seed) result(res) bind(c)
         import
         integer(c_int), value :: seed
         type(c_ptr) :: res
      end function init_rng

      function lrng_create(rng, thread_total, thread_num) result(res) bind(c)
         import
         type(c_ptr), value :: rng
         integer(c_int), value :: thread_total
         integer(c_int), value :: thread_num
         type(c_ptr) :: res
      end function lrng_create

      function rnorm_par(lrng, mean, sd) result(sample) bind(c)
         import
         type(c_ptr), value :: lrng
         real(c_double), value :: mean
         real(c_double), value :: sd
         real(c_double) :: sample
      end function rnorm_par

      function rgamma_par(lrng, shape, scale) result(sample) bind(c)
         import
         type(c_ptr), value :: lrng
         real(c_double), value :: shape
         real(c_double), value :: scale
         real(c_double) :: sample
      end function rgamma_par

      function rpois_par(lrng, lambda) result(sample) bind(c)
         import
         type(c_ptr), value :: lrng
         real(c_double), value :: lambda
         integer(c_int) :: sample
      end function rpois_par

      function runif_par(lrng) result(sample) bind(c)
         import
         type(c_ptr), value :: lrng
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

   end interface

end module interface
