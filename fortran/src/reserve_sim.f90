subroutine reserve_sim(triangle, n_boot, n_dev, config, n_config, results)

   use iso_fortran_env, only: dp => real64
   use dispmodule
   use progressbar, only: progress_bar

   implicit none

   interface
   subroutine reserve_boot(triangle, n_boot, n_dev, reserve, &
      resids_type_in, boot_type_in, dist_in, excl_resids, n_excl_resids, log_unit)   

         use, intrinsic :: iso_fortran_env, only: dp => real64

         integer, intent(in) :: n_dev, n_boot
         real(dp), intent(in) :: triangle(n_dev, n_dev)
         real(dp), intent(inout) :: reserve(n_boot)
         integer, intent(in) :: n_excl_resids
         integer, intent(in) :: excl_resids(2, n_excl_resids)
         integer, intent(in):: dist_in, resids_type_in, boot_type_in
         integer, intent(in) :: log_unit

      end subroutine reserve_boot

      function single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist) result(sim_triangle)
         
         use, intrinsic :: iso_fortran_env, only: dp => real64

         integer, intent(in):: outlier_rowidx, outlier_colidx
         real(dp) :: factor
         real(dp) :: gamma_shape, gamma_rate
         real(dp), intent(in) :: init_col(:), dev_facs(:), sigmas(:)
         real(dp), allocatable :: sim_triangle(:, :)
         integer, intent(in) :: dist
      end function single_outlier

      end interface

      integer, intent(in) :: n_dev, n_boot, n_config
      ! configuration must be specified in format:
      ! outlier point x, y; perturbation factor; excluded point x, y; residuals type; bootstrap type; distribution
      real(dp), intent(in) :: config(n_config, 8), triangle(n_dev, n_dev)
      real(dp), intent(inout) :: results(n_boot*n_config, 9)

      integer :: i, j, k, n_rows
      real(dp) :: indiv_dev_facs(n_dev - 1, n_dev - 1), dev_facs(n_dev - 1), sigmas(n_dev - 1)
      real(dp) :: reserve(n_boot), init_col(n_dev), sim_triangle(n_dev, n_dev)

      character(:), allocatable :: log_path
      integer :: log_unit

      log_path = "/home/othman/repos/masters_thesis/fortran/log/reserve_boot.log"
      
      open(newunit=log_unit, file=log_path, buffered='yes', blocksize=209715200)

      init_col = triangle(:, 1)

      do j = 1, n_dev - 1

         n_rows = n_dev - j

         indiv_dev_facs(1:n_rows, j) = triangle(1:n_rows, j + 1) / triangle(1:n_rows, j)

         dev_facs(j) = sum(triangle(1:n_rows, j + 1)) / sum(triangle(1:n_rows, j))

         if (j < n_dev - 1) then
            sigmas(j) = sqrt(sum(triangle(1:n_rows, j) * (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) ** 2) / (n_rows - 1))
         else
            sigmas(j) = sqrt(min(sigmas(j - 1) ** 2, sigmas(j - 2) ** 2, sigmas(j - 1) ** 4 / sigmas(j - 2) ** 2))
         end if

      end do

      do i = 1, n_config

         call progress_bar(i, n_config)
         flush(6)

         sim_triangle = single_outlier(int(config(i, 1)), int(config(i, 2)), config(i, 3), &
         init_col, dev_facs, sigmas, int(config(i, 8)))

         call reserve_boot(sim_triangle, n_boot, n_dev, reserve, &
         int(config(i, 6)), int(config(i, 7)), int(config(i, 8)), int(config(i, 4:5)), 1, log_unit)

         results(((i - 1)*n_boot + 1):(i*n_boot), 1:8) = transpose(spread(config(i, :), 2, n_boot))
         results(((i - 1)*n_boot + 1):(i*n_boot), 9) = reserve

      end do

      write(6, '(/)')

      close(log_unit)


   end subroutine reserve_sim

