program test

   use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
   use :: iso_fortran_env, only: dp => real64
   use dispmodule

   implicit none

   interface
      subroutine reserve_sim(triangle, n_boot, n_dev, config, n_config, results)
         use :: iso_fortran_env, only: dp => real64         

         integer, intent(in) :: n_dev, n_boot, n_config
         ! configuration must be specified in format:
         ! outlier point x, y; perturbation factor; excluded point x, y; residuals type; bootstrap type; distribution
         real(dp), intent(in) :: config(n_config, 8), triangle(n_dev, n_dev)
         real(dp), intent(inout) :: results(n_config, 9)

         integer :: i
         real(dp) :: reserve(n_boot)

         character(:), allocatable :: log_path
         integer :: log_unit
      end subroutine reserve_sim
   end interface

   integer, parameter :: n_dev = 7
   integer, parameter :: n_boot = 1e3
   real(dp) :: triangle(n_dev, n_dev), reserve(n_boot)
   real(dp), allocatable :: config(:, :), results(:, :)
   integer :: fileunit

   triangle = transpose(reshape( &
      [3511._dp, 6726._dp, 8992._dp, 10704._dp, 11763._dp, 12350._dp, 12690._dp, &
      4001._dp, 7703._dp, 9981._dp, 11161._dp, 12117._dp, 12746._dp, 0._dp, &
      4355._dp, 8287._dp, 10233._dp, 11755._dp, 12993._dp, 0._dp, 0._dp, &
      4295._dp, 7750._dp, 9773._dp, 11093._dp, 0._dp, 0._dp, 0._dp, &
      4150._dp, 7897._dp, 10217._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      5102._dp, 9650._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      6283._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp], &
      [n_dev, n_dev]))

   config = reshape([1., 4., 0.75, 1., 2., 1., 1., 1.], [1, 8])

   allocate(results(n_boot, 9))

   call reserve_sim(triangle, n_boot, n_dev, config, 1, results)

   call disp(results)

end program test
