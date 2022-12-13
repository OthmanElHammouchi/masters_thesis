program test

   use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
   use :: iso_fortran_env, only: dp => real64, di => int64
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
   integer :: file_unit
   integer :: i, j, k
   integer, parameter:: n_points = (n_dev**2 - n_dev)/2
   integer :: points(n_points, 2)
   integer :: n_config

   triangle = transpose(reshape( &
      [3511._dp, 6726._dp, 8992._dp, 10704._dp, 11763._dp, 12350._dp, 12690._dp, &
      4001._dp, 7703._dp, 9981._dp, 11161._dp, 12117._dp, 12746._dp, 0._dp, &
      4355._dp, 8287._dp, 10233._dp, 11755._dp, 12993._dp, 0._dp, 0._dp, &
      4295._dp, 7750._dp, 9773._dp, 11093._dp, 0._dp, 0._dp, 0._dp, &
      4150._dp, 7897._dp, 10217._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      5102._dp, 9650._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      6283._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp], &
      [n_dev, n_dev]))

   k = 1
   do j = 2, n_dev
      do i = 1, n_dev + 1 - j
         points(k, :) = [i, j]
         k = k + 1
      end do
   end do

   allocate(config(size(points, 1)**2, 8))

   k = 1
   do j = 1, n_points
      do i = 1, n_points
         config(k, :) = [real(points(j, 1), dp), real(points(j, 2), dp), 0.75_dp, &
         real(points(i, 1), dp), real(points(i, 2), dp), 1._dp, 1._dp, 1._dp]
         k = k + 1
      end do
   end do

   n_config = size(config, 1)

   allocate(results(n_config*n_boot, 9))

   call reserve_sim(triangle, n_boot, n_dev, config, n_config, results)

   open(newunit=file_unit, file="test/test.dat")

   call disp(results, unit=file_unit)

end program test
