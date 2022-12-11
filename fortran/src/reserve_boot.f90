! subroutine reserve_sim(triangle, n_dev, config, results, dims)
   
!    use iso_fortran_env, only: dp => real64
   
!    implicit none

!    integer, intent(in) :: dims(2), n_dev
!    ! configuration must be specified in the format:
!    ! outlier point x, y; perturbation factor; excluded point x, y; resids type; bootstrap type; distribution
!    real(dp), intent(in) :: config(dims(1), dims(2)), results(dims(1), dims(2)), triangle(n_dev, n_dev)
!    integer :: i
!    real(dp) :: reserve(int(1e3))
!    integer :: excl_resids(2, 1)

!    do i, dims(1)
!       excl_resids(1, 1) = config(i, 4)
!       excl_resids(1, 2) = config(i, 5)
!       reserve_boot(int(1e3), n_dev, triangle, reserve, config(i, 8), config(i, 6), config(i, 7), )
!       results(i, :) = reserve 
!    end do

! end subroutine reserve_sim

subroutine reserve_boot(n_boot, n_dev, triangle, reserve, &
   distribution, resids_type, bootstrap_type, excl_resids, excl_resids_n_cols)

   use random, only: norm => random_normal, gamma => random_gamma
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

   implicit none

   integer :: i, j, i_diag, i_boot, n_rows, n_resids, k

   integer, intent(in) :: n_dev, n_boot
   real(dp), intent(in) :: triangle(n_dev, n_dev)
   real(dp), intent(inout) :: reserve(n_boot)

   ! Because R's foreign function interface does not allow passing strings to Fortran, we encode the different simulation options as follows:
   ! distribution: 1 = normal, 2 = gamma
   ! resids_type: 1 = raw, 2 = scaled, 3 = parametric
   ! bootstrap_type: 1 = conditional, 2 = unconditional
   ! We also cannot pass lists to Fortran, so the excluded residuals are given as an array whose columns contain the excluded points' coordinates
   integer, intent(in) :: excl_resids_n_cols
   integer, intent(in) :: excl_resids(2, excl_resids_n_cols)
   integer :: distribution, resids_type, bootstrap_type

   real(dp) :: dev_facs(n_dev - 1), sigmas(n_dev - 1), latest(n_dev), &
      indiv_dev_facs(n_dev - 1, n_dev - 1), resids(n_dev - 1, n_dev - 1), &
      scale_facs(n_dev - 1, n_dev - 1), resampled_triangle(n_dev, n_dev), &
      resids_boot(n_dev - 1, n_dev - 1), indiv_dev_facs_boot(n_dev - 1, n_dev - 1), &
      dev_facs_boot(n_dev - 1), sigmas_boot(n_dev - 1), triangle_boot(n_dev, n_dev)

   real(dp), allocatable :: flat_resids(:)
   real(dp) :: gamma_shape, gamma_rate

   real(dp) :: flat_resampled_triangle(n_dev**2)

   n_resids = ((n_dev - 1)**2 + (n_dev - 1))/2 - 1

   do j = 1, n_dev - 1

      n_rows = n_dev - j

      indiv_dev_facs(1:n_rows, j) = triangle(1:n_rows, j + 1) / triangle(1:n_rows, j)

      dev_facs(j) = sum(triangle(1:n_rows, j + 1)) / sum(triangle(1:n_rows, j))

      if (j < n_dev - 1) then
         sigmas(j) = sqrt(sum(triangle(1:n_rows, j) * (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) ** 2) / (n_rows - 1))
      else
         sigmas(j) = sqrt(min(sigmas(j - 1) ** 2, sigmas(j - 2) ** 2, sigmas(j - 1) ** 4 / sigmas(j - 2) ** 2))
      end if

      if (resids_type == 1) then
         resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * sqrt(triangle(1:n_rows, j)) / sigmas(j)
      else if (resids_type == 2) then
         if (j < n_dev - 1) then
            scale_facs(1:n_rows, j) = triangle(1:n_rows, j) / sum(triangle(1:n_rows, j))
         else
            scale_facs(1:n_rows, j) = 1
         end if
         resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * &
            sqrt(triangle(1:n_rows, j)) / (sigmas(j) * scale_facs(1:n_rows, j))
      end if

   end do

   ! remove excluded residuals
   if (excl_resids(1, 1) /= 0) then
      do j = 1, excl_resids_n_cols  
         resids(excl_resids(1, j), excl_resids(2, j) - 1) = 0
         n_resids = n_resids - 1
      end do
   end if

   if (resids_type /= 3) flat_resids = pack(resids, .not. abs(resids - 0) < 1e-10)

   ! main loop
   main_loop: do i_boot = 1, n_boot
      ! resample residual
      if (resids_type == 3) then
         flat_resids = [(norm(), i = 1, n_resids)]
      end if
      do j = 1, n_dev - 1
         do i = 1, n_dev - j
            resids_boot(i, j) = flat_resids(1 + int(n_resids * rand()))
         end do
      end do

      ! compute bootstrapped quantities
      if (bootstrap_type == 1) then
         do j = 1, n_dev - 1

            n_rows = n_dev - j

            indiv_dev_facs_boot(1:n_rows, j) = dev_facs(j) + resids_boot(1:n_rows, j) * sigmas(j) / sqrt(triangle(1:n_rows, j))

            dev_facs_boot(j) = sum(triangle(1:n_rows, j) * indiv_dev_facs_boot(1:n_rows, j)) / sum(triangle(1:n_rows, j))

            if (j < n_dev - 1) then
               sigmas_boot(j) = sqrt(sum(triangle(1:n_rows, j) * &
                  (indiv_dev_facs_boot(1:n_rows, j) - dev_facs_boot(j)) ** 2) / (n_rows - 1))
            else
               sigmas_boot(j) = sqrt(min(sigmas_boot(j - 1) ** 2, &
                  sigmas_boot(j - 2) ** 2, &
                  sigmas_boot(j - 1) ** 4 / sigmas_boot(j - 2) ** 2))
            end if
         end do
      else if (bootstrap_type == 2) then
         resampled_triangle(:, 1) = triangle(:, 1)
         do j = 2, n_dev
            n_rows = n_dev + 1 - j
            resampled_triangle(1:n_rows, j) = dev_facs(j - 1) * resampled_triangle(1:n_rows, j - 1) + &
               sigmas(j - 1) * sqrt(resampled_triangle(1:n_rows, j - 1)) * resids_boot(1:n_rows, j - 1)

            indiv_dev_facs_boot(1:n_rows, j - 1) = dev_facs(j - 1) + resids_boot(1:n_rows, j - 1) * &
               sigmas(j - 1) / sqrt(resampled_triangle(1:n_rows, j - 1))

            dev_facs_boot(j - 1) = sum(resampled_triangle(1:n_rows, j - 1) * &
               indiv_dev_facs_boot(1:n_rows, j - 1)) / sum(resampled_triangle(1:n_rows, j - 1))

            if (j < n_dev) then
               sigmas_boot(j - 1) = sqrt(sum(resampled_triangle(1:n_rows, j - 1) * &
                  (indiv_dev_facs_boot(1:n_rows, j - 1) - dev_facs_boot(j - 1)) ** 2) / (n_rows - 1))
            else
               sigmas_boot(j - 1) = sqrt(min(sigmas_boot(j - 2) ** 2, &
                  sigmas_boot(j - 3) ** 2, &
                  sigmas_boot(j - 2) ** 4 / sigmas_boot(j - 3) ** 2))
            end if
         end do
         
         flat_resampled_triangle = pack(resampled_triangle, .true.)
         if (any(flat_resampled_triangle < 0)) then 
            reserve(i_boot) = 0.
            cycle main_loop
         end if

      end if

      triangle_boot = triangle

      if (distribution == 1) then
         do i_diag = 1, n_dev - 1
            do i = i_diag + 1, n_dev
               j = n_dev + i_diag + 1 - i
               triangle_boot(i, j) = dev_facs_boot(j - 1) * triangle_boot(i, j - 1) + &
               norm()*real((sigmas_boot(j - 1) * sqrt(triangle_boot(i, j - 1) )))
            end do
         end do

         do j = 1, n_dev
            latest(j) = triangle_boot(n_dev + 1 - j, j)
         end do

         reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)

      else if (distribution == 2) then
         do i_diag = 1, n_dev - 1
            do i = i_diag + 1, n_dev
               j = n_dev + i_diag + 1 - i
               gamma_shape = (dev_facs_boot(j - 1)**2 * triangle_boot(i, j - 1)) / sigmas_boot(j - 1) **2
               if (gamma_shape <= tiny(1.0)) then
                  cycle main_loop
               end if
               gamma_rate = dev_facs_boot(j - 1) / sigmas_boot(j - 1) ** 2
               triangle_boot(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp) 
               if (triangle_boot(i, j) <= 0) then
                  cycle main_loop
               end if
            end do
         end do

         do j = 1, n_dev
            latest(j) = triangle_boot(n_dev + 1 - j, j)
         end do

         reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)
      end if

   end do main_loop

end subroutine reserve_boot