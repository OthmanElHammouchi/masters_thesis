include 'mkl_vsl.f90'

subroutine reserve_boot(triangle, n_boot, n_dev, reserve, &
   dist_in, resids_type_in, boot_type_in, excl_resids, n_excl_resids, log_unit)

   !  use random, only: norm => random_normal, gamma => random_gamma
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use dispmodule
   use ifport
   use mkl_vsl

   implicit none

   type (vsl_stream_state) :: stream
   integer(kind=4) errcode
   integer brng,method,seed,n
   real(dp) :: r(n_boot)

   integer :: i, j, i_diag, i_boot, n_rows, n_resids, k

   integer, intent(in) :: n_dev, n_boot
   real(dp), intent(in) :: triangle(n_dev, n_dev)
   real(dp), intent(inout) :: reserve(n_boot)


   ! can't pass lists to R, we have to pass an array of coords
   integer, intent(in) :: n_excl_resids
   integer, intent(in) :: excl_resids(2, n_excl_resids)

   integer, intent(in):: dist_in, resids_type_in, boot_type_in
   character(999) :: dist, resids_type, boot_type

   real(dp) :: dev_facs(n_dev - 1), sigmas(n_dev - 1), latest(n_dev), &
      indiv_dev_facs(n_dev - 1, n_dev - 1), resids(n_dev - 1, n_dev - 1), &
      scale_facs(n_dev - 1, n_dev - 1), resampled_triangle(n_dev, n_dev), &
      resids_boot(n_dev - 1, n_dev - 1), indiv_dev_facs_boot(n_dev - 1, n_dev - 1), &
      dev_facs_boot(n_dev - 1), sigmas_boot(n_dev - 1), triangle_boot(n_dev, n_dev)

   real(dp), allocatable :: flat_resids(:)
   real(dp) :: gamma_shape, gamma_rate

   real(dp) :: flat_resampled_triangle(n_dev**2)

   integer, intent(in) :: log_unit


   brng = VSL_BRNG_MT19937
   method = VSL_RNG_METHOD_GAMMA_GNORM
   seed = 777

   errcode = vslnewstream( stream, brng,  seed )

   ! can't pass strings to R, we have to encode options as integers:
   ! distribution: 1 = normal, 2 = gamma
   ! resids_type: 1 = raw, 2 = scaled, 3 = parametric
   ! bootstrap_type: 1 = conditional, 2 = unconditional

   if (resids_type_in == 1) then
      resids_type = "raw"
   else if (resids_type_in == 2) then
      resids_type = "scaled"
   else if (resids_type_in == 3) then
      resids_type = "parametric"
   else
      stop "Unknown resids type"
   end if

   if (boot_type_in == 1) then
      boot_type = "conditional"
   else if (boot_type_in == 2) then
      boot_type = "unconditional"
   else
      stop "Unknown boot type"
   end if

   if (dist_in == 1) then
      dist = "normal"
   else if (dist_in == 2) then
      dist = "gamma"
   else
      stop "Unknown dist"
   end if

   if (resids_type == "parametric" .and. dist == "gamma") then
      stop "Parametric resids not supported for gamma distribution"
   end if

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

      if (resids_type == "raw") then

         resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * sqrt(triangle(1:n_rows, j)) / sigmas(j)

      else if (resids_type == "scaled") then

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
      do j = 1, n_excl_resids
         resids(excl_resids(1, j), excl_resids(2, j) - 1) = 0
         n_resids = n_resids - 1
      end do
   end if

   if (resids_type /= "parametric") flat_resids = pack(resids, .not. abs(resids - 0) < 1e-10)

   ! main loop
   main_loop: do i_boot = 1, n_boot
      ! resample residual
      if (resids_type == "parametric") then
         flat_resids = [(norm(), i = 1, n_resids)]
      end if
      do j = 1, n_dev - 1
         do i = 1, n_dev - j
            resids_boot(i, j) = flat_resids(1 + int(n_resids * rand()))
         end do
      end do

      ! compute bootstrapped quantities
      if (boot_type == "conditional") then

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

      else if (boot_type == "unconditional") then

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

         if (any(sigmas_boot < 0)) then
            write(log_unit, '("Configuration: ", 3a, /)') trim(resids_type), trim(boot_type), trim(dist)
            write(log_unit, '(a, /)') "sigmas_boot: "
            call disp(sigmas_boot, unit=log_unit)
            write(log_unit, '(/, a, /)') "indiv_devfacs_boot: "
            call disp(indiv_dev_facs_boot, unit=log_unit)
            write(log_unit, '(/, a, /)') "devfacs_boot: "
            call disp(dev_facs_boot, unit=log_unit)
            write(log_unit, '(/, a, /)') "resids_boot: "
            call disp(resids_boot, unit=log_unit)
            stop "Fatal error: negative sigma"
         end if

         flat_resampled_triangle = pack(resampled_triangle, .true.)

      end if

      triangle_boot = triangle

      if (dist == "normal") then
         do i_diag = 1, n_dev - 1
            do i = i_diag + 1, n_dev

               j = n_dev + i_diag + 1 - i

               triangle_boot(i, j) = dev_facs_boot(j - 1) * triangle_boot(i, j - 1) + &
                  norm()*real((sigmas_boot(j - 1) * sqrt(triangle_boot(i, j - 1) )))

               if (triangle_boot(i, j) <= 0) then

                  write(log_unit, '(a, 2i2, //)') "Negative draw at ", i, j
                  write(log_unit, '("Configuration: ", 3a, /)') trim(resids_type), trim(boot_type), trim(dist)
                  write(log_unit, '(a, /)') "sigmas_boot: "
                  call disp(sigmas_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "triangle_boot:"
                  call disp(triangle_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "indiv_devfacs_boot:"
                  call disp(indiv_dev_facs_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "devfacs_boot:"
                  call disp(dev_facs_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "resids_boot:"
                  call disp(resids_boot, unit=log_unit)

                  cycle main_loop
               end if

            end do
         end do

         do j = 1, n_dev
            latest(j) = triangle_boot(n_dev + 1 - j, j)
         end do

         reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)

      else if (dist == "gamma") then
         do i_diag = 1, n_dev - 1
            do i = i_diag + 1, n_dev

               j = n_dev + i_diag + 1 - i

               gamma_shape = (dev_facs_boot(j - 1)**2 * triangle_boot(i, j - 1)) / sigmas_boot(j - 1) **2

               if (gamma_shape <= tiny(1.0)) then

                  write(log_unit, '(a, 2i2, //)') "Negative gamma shape parameter at ", i, j
                  write(log_unit, '("Configuration: ", 3a, /)') trim(resids_type), trim(boot_type), trim(dist)
                  write(log_unit, '(a, /)') "sigmas_boot: "
                  call disp(sigmas_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "triangle_boot:"
                  call disp(triangle_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "indiv_devfacs_boot:"
                  call disp(indiv_dev_facs_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "devfacs_boot:"
                  call disp(dev_facs_boot, unit=log_unit)
                  write(log_unit, '(/, a, /)') "resids_boot:"
                  call disp(resids_boot, unit=log_unit)

                  cycle main_loop
               end if

               gamma_rate = dev_facs_boot(j - 1) / sigmas_boot(j - 1) ** 2

               triangle_boot(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

            end do
         end do

         do j = 1, n_dev
            latest(j) = triangle_boot(n_dev + 1 - j, j)
         end do

         reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)
      end if

   end do main_loop

   errcode = vsldeletestream( stream )

end subroutine reserve_boot

function single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist) result(sim_triangle)

   use random, only: norm => random_normal, gamma => random_gamma
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none

   integer, intent(in):: outlier_rowidx, outlier_colidx
   real(dp) :: factor
   real(dp) :: gamma_shape, gamma_rate
   real(dp), allocatable:: sim_triangle(:, :)
   real(dp), intent(in) :: init_col(:), dev_facs(:), sigmas(:)
   integer, intent(in) :: dist

   integer :: n_dev
   integer :: i, j

   n_dev = size(init_col)

   allocate(sim_triangle(n_dev, n_dev))

   sim_triangle(:, 1) = init_col

   if (dist == 1) then
      do j = 2, n_dev
         do i = 1, n_dev + 1 - j
            if (i == outlier_rowidx) cycle
            sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + &
               norm() * sigmas(j - 1) * sqrt(sim_triangle(i, j - 1))
         end do
      end do

      if (outlier_colidx > 2) then
         do j = 2, outlier_colidx - 1
            sim_triangle(outlier_rowidx, j) = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + &
               norm() * sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
         end do
      end if

      sim_triangle(outlier_rowidx, outlier_colidx) = &
         factor * dev_facs(outlier_colidx - 1) * sim_triangle(outlier_rowidx, outlier_colidx - 1) + &
         norm() * sigmas(outlier_colidx - 1) * sqrt(sim_triangle(outlier_rowidx, outlier_colidx - 1))

      if (outlier_colidx < n_dev) then
         do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
            sim_triangle(outlier_rowidx, j) = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + &
               norm() * sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
         end do
      end if

   else if (dist == 2) then

      do j = 2, n_dev
         do i = 1, n_dev + 1 - j
            if (i == outlier_rowidx) cycle

            gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
            gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

            sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

         end do
      end do

      if (outlier_colidx > 2) then
         do j = 2, outlier_colidx - 1
            gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
            gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

            sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)
         end do
      end if

      gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
      gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

      sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

      if (outlier_colidx < n_dev) then
         do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
            gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
            gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

            sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)
         end do
      end if

   end if
end function single_outlier
