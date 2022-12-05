program test

   use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
   use :: iso_fortran_env, only: dp => real64

   implicit none

   interface
      subroutine reserve_boot(n_boot, n_dev, triangle, reserve, distribution_in, resids_type_in, bootstrap_type_in, excl_resids)

        use iso_fortran_env, only: dp => real64

        implicit none

         integer :: i, j, i_diag, i_boot, n_rows, n_resids

         integer, intent(in) :: n_dev, n_boot
         real(dp), intent(in) :: triangle(n_dev, n_dev)
         real(dp), intent(inout) :: reserve(n_boot)

         integer, optional, intent(in) :: excl_resids(:, :)
         integer, optional, intent(in) :: distribution_in, resids_type_in, bootstrap_type_in
         integer :: distribution, resids_type, bootstrap_type

         real(dp) :: dev_facs(n_dev - 1), sigmas(n_dev - 1), latest(n_dev), &
            indiv_dev_facs(n_dev - 1, n_dev - 1), resids(n_dev - 1, n_dev - 1), &
            boot_resids(n_dev - 1, n_dev - 1), boot_indiv_dev_facs(n_dev - 1, n_dev - 1), &
            boot_dev_facs(n_dev - 1), boot_sigmas(n_dev - 1), boot_triangle(n_dev, n_dev)

         real(dp), allocatable :: flat_resids(:)
         integer :: excl_resids_shape(2)

         real :: gamma_shape, gamma_rate
      end subroutine reserve_boot
   end interface


   real(dp) :: nan
   integer, parameter :: n_dev = 7
   integer, parameter :: n_boot = 1e3
   real(dp) :: triangle(n_dev, n_dev), reserve(n_boot)
   integer :: i, j
   integer :: fileunit

   nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

   triangle = &
      transpose(reshape([3511._dp, 6726._dp, 8992._dp, 10704._dp, 11763._dp, 12350._dp, 12690._dp, &
      4001._dp, 7703._dp, 9981._dp, 11161._dp, 12117._dp, 12746._dp, nan, &
      4355._dp, 8287._dp, 10233._dp, 11755._dp, 12993._dp, nan, nan, &
      4295._dp, 7750._dp, 9773._dp, 11093._dp, nan, nan, nan, &
      4150._dp, 7897._dp, 10217._dp, nan, nan, nan, nan, &
      5102._dp, 9650._dp, nan, nan, nan, nan, nan, &
      6283._dp, nan, nan, nan, nan, nan, nan], [n_dev, n_dev]))

   call reserve_boot(n_boot, n_dev, triangle, reserve, 2)

   open(newunit=fileunit, file='test/test.dat')
   write(fileunit, '(f10.2)') reserve

end program
