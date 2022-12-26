include 'mkl_vsl.f90'

subroutine reserve_sim(triangle, n_boot, n_dev, config, config_spec, results)

   use iso_fortran_env, only: dp => real64
   use dispmodule
   use progressbar, only: progress_bar

#ifdef __INTEL_COMPILER
   use ifport
   use mkl_vsl
#elif defined __GFORTRAN__
   use random, only: norm => random_normal, gamma => random_gamma
#endif

   implicit none

   integer, intent(in) :: n_dev, n_boot
   integer, intent(in) :: config_spec(3)
   real(dp), intent(in) :: config(config_spec(1), config_spec(2)), triangle(n_dev, n_dev)
   real(dp), intent(inout) :: results(n_boot*config_spec(1), config_spec(2) + 1)

   integer :: i, j, k, i_sim, n_rows, n_config, &
      outlier_rowidx, outlier_colidx, outlier_diagidx, excl_diagidx, excl_rowidx
   integer, allocatable :: excl_resids(:, :)
   real(dp) :: indiv_dev_facs(n_dev - 1, n_dev - 1), dev_facs(n_dev - 1), sigmas(n_dev - 1)
   real(dp) :: reserve(n_boot), init_col(n_dev), sim_triangle(n_dev, n_dev)
   real(dp) :: factor

   character(:), allocatable :: resids_type, boot_type, dist, config_type
   character(:), allocatable :: config_type_key(:), resids_type_key(:), boot_type_key(:), dist_key(:)

   character(:), allocatable :: log_path
   integer :: log_unit

   log_path = "./fortran/log/reserve_boot.log"

#ifdef __INTEL_COMPILER
   open(newunit=log_unit, file=log_path, buffered='yes', blocksize=209715200)
#elif defined __GFORTRAN__
   open(newunit=log_unit, file=log_path)
#endif

   config_type_key = ["single", "calendar", "origin"]
   resids_type_key = ["raw", "scaled", "parametric"]
   boot_type_key = ["conditional", "unconditional"]
   dist_key = ["normal", "gamma"]

   config_type = config_type_key(int(config_spec(3)))

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

   n_config = config_spec(1)

   do i_sim = 1, n_config

      call progress_bar(i_sim, n_config)
      flush(6)

      if (config_type == "single") then

         resids_type = resids_type_key(int(config(i_sim, 6)))
         boot_type = boot_type_key(int(config(i_sim, 7)))
         dist = dist_key(int(config(i_sim, 8)))

         if (resids_type == "parametric" .and. dist == "gamma") then
            stop "Fatal error: parametric resids not supported for gamma distribution"
         end if

         excl_resids(1, :) = int(config(i_sim, 4:5))
         outlier_rowidx = int(config(i_sim, 1))
         outlier_colidx = int(config(i_sim, 2))
         factor = config(i_sim, 3)

         sim_triangle = single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist)

         reserve = reserve_boot(sim_triangle, n_boot, resids_type, boot_type, dist, excl_resids, log_unit)

         results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:config_spec(2)) = transpose(spread(config(i_sim, :), 2, n_boot))
         results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), config_spec(2) + 1) = reserve

      else if (config_type == "calendar") then

         resids_type = resids_type_key(int(config(i_sim, 4)))
         boot_type = boot_type_key(int(config(i_sim, 5)))
         dist = dist_key(int(config(i_sim, 6)))

         if (resids_type == "parametric" .and. dist == "gamma") then
            stop "Fatal error: parametric resids not supported for gamma distribution"
         end if

         outlier_diagidx = int(config(i_sim, 1))
         excl_diagidx = int(config(i_sim, 3))
         factor = config(i_sim, 2)

         allocate(excl_resids(n_dev - excl_diagidx, 2), source=0)

         k = 1
         do j = 2, n_dev
            i = n_dev + 2 - excl_diagidx - j
            if (i <= 0) cycle
            excl_resids(k, :) = [i, j]
            k = k + 1
         end do

         sim_triangle = calendar_outlier(outlier_diagidx, factor, triangle, dev_facs, sigmas, dist)

         reserve = reserve_boot(sim_triangle, n_boot, resids_type, boot_type, dist, excl_resids, log_unit)

         deallocate(excl_resids)

         results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:config_spec(2)) = transpose(spread(config(i_sim, :), 2, n_boot))
         results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), config_spec(2) + 1) = reserve

      else if (config_type == "origin") then

         resids_type = resids_type_key(int(config(i_sim, 4)))
         boot_type = boot_type_key(int(config(i_sim, 5)))
         dist = dist_key(int(config(i_sim, 6)))

         if (resids_type == "parametric" .and. dist == "gamma") then
            stop "Fatal error: parametric resids not supported for gamma distribution"
         end if

         outlier_rowidx = int(config(i_sim, 1))
         excl_rowidx = int(config(i_sim, 3))
         factor = config(i_sim, 2)

         allocate(excl_resids(n_dev - excl_rowidx, 2), source=0)

         k = 1
         do j = 2, n_dev + 1 - excl_rowidx
            excl_resids(k, :) = [excl_rowidx, j]
            k = k + 1
         end do

         sim_triangle = origin_outlier(outlier_rowidx, factor, triangle, dev_facs, sigmas, dist)

         reserve = reserve_boot(sim_triangle, n_boot, resids_type, boot_type, dist, excl_resids, log_unit)

         deallocate(excl_resids)

         results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:config_spec(2)) = transpose(spread(config(i_sim, :), 2, n_boot))
         results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), config_spec(2) + 1) = reserve

      end if

   end do

   write(6, '(/)')

   close(log_unit)

contains

   function reserve_boot(triangle, n_boot, resids_type, boot_type, dist, excl_resids, log_unit) result(reserve)

      integer :: i, j, i_diag, i_boot, n_rows, n_resids, k

      integer, intent(in) :: n_boot
      real(dp), intent(in) :: triangle(:, :)
      integer, intent(in) :: log_unit
      integer, intent(in) :: excl_resids(:, :)

      character(:), allocatable :: dist, resids_type, boot_type

      integer :: n_dev

      real(dp), allocatable :: dev_facs(:), sigmas(:), latest(:), &
         indiv_dev_facs(:, :), resids(:, :), &
         scale_facs(:, :), resampled_triangle(:, :, :)
      real(dp), allocatable :: &
         resids_boot(:, :, :), indiv_dev_facs_boot(:, :, :), &
         dev_facs_boot(:, :), sigmas_boot(:, :), triangle_boot(:, :)

      real(dp) :: reserve(n_boot)
      real(dp), allocatable :: flat_resids(:)
      real(dp) :: gamma_shape, gamma_rate

#ifdef __INTEL_COMPILER
      type (vsl_stream_state) :: stream
      integer(kind=4) :: errcode
      integer :: brng, method, sd, n
      real(dp) :: r(1)

      brng = VSL_BRNG_MT19937
      sd = irand()
      errcode = vslnewstream(stream, brng,  sd)
#endif

      n_dev = size(triangle, dim=1)

      allocate(dev_facs(n_dev - 1), sigmas(n_dev - 1), latest(n_dev), &
         indiv_dev_facs(n_dev - 1, n_dev - 1), resids(n_dev - 1, n_dev - 1), &
         scale_facs(n_dev - 1, n_dev - 1), resampled_triangle(n_dev, n_dev, n_dev - 1), &
         resids_boot(n_dev - 1, n_dev - 1, n_dev - 1), indiv_dev_facs_boot(n_dev - 1, n_dev - 1, n_dev - 1), &
         dev_facs_boot(n_dev - 1, n_dev - 1), sigmas_boot(n_dev - 1, n_dev - 1), triangle_boot(n_dev, n_dev), &
         source=0._dp)


      n_resids = ((n_dev - 1)**2 + (n_dev - 1))/2 - 1 !discard residual from upper right corner

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
               scale_facs(1:n_rows, j) = sqrt(1 - triangle(1:n_rows, j) / sum(triangle(1:n_rows, j)))
            else
               scale_facs(1:n_rows, j) = 1
            end if

            resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * &
               sqrt(triangle(1:n_rows, j)) / (sigmas(j) * scale_facs(1:n_rows, j))

         end if

      end do

      ! remove excluded residuals
      do i = 1, size(excl_resids, dim=1)
         resids(excl_resids(i, 1), excl_resids(i, 2) - 1) = 0
         n_resids = n_resids - 1
      end do

      allocate(flat_resids(n_resids))

      if (resids_type /= "parametric") flat_resids = pack(resids, .not. abs(resids - 0) < 1e-10)

      ! main loop
      main_loop: do i_boot = 1, n_boot
         ! resample residuals
         if (resids_type == "parametric") then
#ifdef __INTEL_COMPILER
            method = VSL_RNG_METHOD_GAUSSIAN_ICDF
            do i = 1, n_dev - 1
               do j = 1, n_dev - 1
                  do k = 1, n_dev - 1
                     errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)
                     resids_boot(i, j, k) = r(1)
                  end do
               end do
            end do
#elif defined __GFORTRAN__
            do i = 1, n_dev - 1
               do i = 1, n_dev - 1
                  do k = 1, n_dev - 1
                     resids_boot(i, j, k) = norm()
                  end do
               end do
            end do
#endif
         else
            do i = 1, n_dev - 1
               do j = 1, n_dev - 1
                  do k = 1, n_dev - 1
                     resids_boot(i, j, k) = flat_resids(1 + int(n_resids * rand()))
                  end do
               end do
            end do
         end if

         indiv_dev_facs_boot = 0
         dev_facs_boot = 0
         sigmas_boot = 0

         ! compute bootstrapped quantities
         if (boot_type == "conditional") then
            do k = 1, n_dev - 1
               do j = 1, n_dev - 1

                  n_rows = n_dev - j

                  indiv_dev_facs_boot(1:n_rows, j, k) = dev_facs(j) + resids_boot(1:n_rows, j, k) * sigmas(j) / sqrt(triangle(1:n_rows, j))

                  dev_facs_boot(j, k) = sum(triangle(1:n_rows, j) * indiv_dev_facs_boot(1:n_rows, j, k)) / sum(triangle(1:n_rows, j))

                  if (j < n_dev - 1) then
                     sigmas_boot(j, k) = sqrt(sum(triangle(1:n_rows, j) * &
                        (indiv_dev_facs_boot(1:n_rows, j, k) - dev_facs_boot(j, k)) ** 2) / (n_rows - 1))
                  else
                     sigmas_boot(j, k) = sqrt(min(sigmas_boot(j - 1, k) ** 2, &
                        sigmas_boot(j - 2, k) ** 2, &
                        sigmas_boot(j - 1, k) ** 4 / sigmas_boot(j - 2, k) ** 2))
                  end if
               end do

               if (any(sigmas_boot < 0) .or. any(isnan(sigmas_boot))) then
#ifdef DEBUG
                  call log_bad_sigma(log_unit, resids_type, boot_type, dist, sigmas_boot, indiv_dev_facs_boot, &
                     dev_facs_boot, resids_boot)
#endif
                  stop "Fatal error: bad sigma"
               end if
            end do

         else if (boot_type == "unconditional") then
            do k = 1, n_dev - 1
               resampled_triangle(:, 1, k) = triangle(:, 1)
               do j = 1, n_dev - 1

                  n_rows = n_dev - j

                  resampled_triangle(1:n_rows, j + 1, k) = dev_facs(j) * resampled_triangle(1:n_rows, j, k) + &
                     sigmas(j) * sqrt(resampled_triangle(1:n_rows, j, k)) * resids_boot(1:n_rows, j, k)

                  if (any(resampled_triangle(1:n_rows, j + 1, k) < 0)) cycle main_loop

                  indiv_dev_facs_boot(1:n_rows, j, k) = dev_facs(j) + resids_boot(1:n_rows, j, k) * &
                     sigmas(j) / sqrt(resampled_triangle(1:n_rows, j, k))

                  dev_facs_boot(j, k) = sum(resampled_triangle(1:n_rows, j, k) * &
                     indiv_dev_facs_boot(1:n_rows, j, k)) / sum(resampled_triangle(1:n_rows, j, k))

                  if (j < n_dev - 1) then
                     sigmas_boot(j, k) = sqrt(sum(resampled_triangle(1:n_rows, j, k) * &
                        (indiv_dev_facs_boot(1:n_rows, j, k) - dev_facs_boot(j, k)) ** 2) / (n_rows - 1))
                  else
                     sigmas_boot(j, k) = sqrt(min(sigmas_boot(j - 1, k) ** 2, &
                        sigmas_boot(j - 2, k) ** 2, &
                        sigmas_boot(j - 1, k) ** 4 / sigmas_boot(j - 2, k) ** 2))
                  end if
               end do

               if (any(sigmas_boot < 0) .or. any(isnan(sigmas_boot))) then
#ifdef DEBUG
                  call log_bad_sigma(log_unit, resids_type, boot_type, dist, sigmas_boot, indiv_dev_facs_boot, &
                     dev_facs_boot, resids_boot)
#endif
                  stop "Fatal error: bad sigma"
               end if

            end do


         end if

         triangle_boot = triangle

         if (dist == "normal") then

#ifdef __INTEL_COMPILER
            method = VSL_RNG_METHOD_GAUSSIAN_ICDF
#endif

            do i_diag = 1, n_dev - 1
               do i = i_diag + 1, n_dev

                  j = n_dev + i_diag + 1 - i

                  errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

                  triangle_boot(i, j) = dev_facs_boot(j - 1, i - 1) * triangle_boot(i, j - 1) + &
                     r(1) * real((sigmas_boot(j - 1, i - 1) * sqrt(triangle_boot(i, j - 1))))

                  if (triangle_boot(i, j) <= 0) then
#ifdef DEBUG

#endif
                     cycle main_loop
                  end if

               end do
            end do

            do j = 1, n_dev
               latest(j) = triangle_boot(n_dev + 1 - j, j)
            end do

            reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)

         else if (dist == "gamma") then

#ifdef __INTEL_COMPILER
            method = VSL_RNG_METHOD_GAMMA_GNORM
#endif

            do i_diag = 1, n_dev - 1
               do i = i_diag + 1, n_dev

                  j = n_dev + i_diag + 1 - i

                  gamma_shape = (dev_facs_boot(j - 1, i - 1)**2 * triangle_boot(i, j - 1)) / sigmas_boot(j - 1, i - 1) **2

                  if (gamma_shape <= tiny(1.0) .or. isnan(gamma_shape)) then
#ifdef DEBUG
                     call log_neg_gamma(log_unit, triangle_boot, resids_type, boot_type, dist, sigmas_boot, indiv_dev_facs_boot, &
                        dev_facs_boot, resids_boot)
#endif
                     cycle main_loop
                  end if

                  gamma_rate = dev_facs_boot(j - 1, i - 1) / sigmas_boot(j - 1, i - 1) ** 2

#ifdef __INTEL_COMPILER
                  errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)
                  triangle_boot(i, j) = r(1)
#elif defined __GFORTRAN__
                  triangle_boot(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)
#endif
               end do
            end do

            do j = 1, n_dev
               latest(j) = triangle_boot(n_dev + 1 - j, j)
            end do

            reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)
         end if

      end do main_loop

#ifdef __INTEL_COMPILER
      errcode = vsldeletestream(stream)
#endif

   end function reserve_boot

   subroutine log_bad_sigma(log_unit, resids_type, boot_type, dist, sigmas_boot, indiv_dev_facs_boot, dev_facs_boot, resids_boot)

      integer, intent(in) :: log_unit
      real(dp), intent(in) :: sigmas_boot(:, :), indiv_dev_facs_boot(:, :, :), dev_facs_boot(:, :), resids_boot(:, :, :)
      character(:), allocatable, intent(in) :: resids_type, boot_type, dist

      write(log_unit, '("Bad sigma at ", i2, /)') k
      write(log_unit, '("Configuration: ", 3(a, x), /)') trim(resids_type), trim(boot_type), trim(dist)
      write(log_unit, '(a, /)') "sigmas_boot: "
      call disp(sigmas_boot(:, k), unit=log_unit)
      write(log_unit, '(/, a, /)') "indiv_dev_facs_boot: "
      call disp(indiv_dev_facs_boot(:, :, k), unit=log_unit)
      write(log_unit, '(/, a, /)') "dev_facs_boot: "
      call disp(dev_facs_boot(:, k), unit=log_unit)
      write(log_unit, '(/, a, /)') "resids_boot: "
      call disp(resids_boot(:, :, k), unit=log_unit)

   end subroutine log_bad_sigma

   subroutine log_neg_gamma(log_unit, triangle_boot, resids_type, boot_type, dist, sigmas_boot, indiv_dev_facs_boot, dev_facs_boot, resids_boot)

      integer, intent(in) :: log_unit
      real(dp), intent(in) :: triangle_boot(:, :), sigmas_boot(:, :), indiv_dev_facs_boot(:, :, :), dev_facs_boot(:, :), resids_boot(:, :, :)
      character(:), allocatable, intent(in) :: resids_type, boot_type, dist

      write(log_unit, '(a, 2i2, /)') "Negative gamma shape parameter ", i_sim, j
      write(log_unit, '("Configuration: ", 3(a, x), /)') trim(resids_type), trim(boot_type), trim(dist)
      write(log_unit, '(a, /)') "sigmas_boot: "
      call disp(sigmas_boot(:, i_sim), unit=log_unit)
      write(log_unit, '(/, a, /)') "triangle_boot:"
      call disp(triangle_boot, unit=log_unit)
      write(log_unit, '(/, a, /)') "indiv_dev_facs_boot:"
      call disp(indiv_dev_facs_boot(:, :, i_sim), unit=log_unit)
      write(log_unit, '(/, a, /)') "dev_facs_boot:"
      call disp(dev_facs_boot(:, i_sim), unit=log_unit)
      write(log_unit, '(/, a, /)') "resids_boot:"
      call disp(resids_boot(:, :, i_sim), unit=log_unit)

   end subroutine log_neg_gamma

   function single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist) result(sim_triangle)

      integer, intent(in):: outlier_rowidx, outlier_colidx
      real(dp), intent(in) :: init_col(:), dev_facs(:), sigmas(:)
      character(:), allocatable, intent(in) :: dist

      real(dp) :: factor
      real(dp) :: gamma_shape, gamma_rate
      real(dp), allocatable:: sim_triangle(:, :)

      integer :: n_dev, i, j

#ifdef __INTEL_COMPILER
      type (vsl_stream_state) :: stream
      integer(kind=4) :: errcode
      integer :: brng, method, sd, n
      real(dp) :: r(1)

      brng = VSL_BRNG_MT19937
      sd = irand()
      errcode = vslnewstream(stream, brng,  sd)
#endif

      n_dev = size(init_col)

      allocate(sim_triangle(n_dev, n_dev))
      sim_triangle = 0

      sim_triangle(:, 1) = init_col

      if (dist == "normal") then

#ifdef __INTEL_COMPILER
         method = VSL_RNG_METHOD_GAUSSIAN_ICDF
#endif

         do j = 2, n_dev
            do i = 1, n_dev + 1 - j
               if (i == outlier_rowidx) cycle

#ifdef __INTEL_COMPILER
               errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

               sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + &
                  r(1) * sigmas(j - 1) * sqrt(sim_triangle(i, j - 1))

#elif defined __GFORTRAN__
               sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + &
                  norm() * sigmas(j - 1) * sqrt(sim_triangle(i, j - 1))

#endif
            end do
         end do

         if (outlier_colidx > 2) then

#ifdef __INTEL_COMPILER
            errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

            do j = 2, outlier_colidx - 1
               sim_triangle(outlier_rowidx, j) = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + &
                  r(1) * sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
            end do
#elif defined __GFORTRAN__
            do j = 2, outlier_colidx - 1
               sim_triangle(outlier_rowidx, j) = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + &
                  norm() * sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
            end do
#endif

         end if

#ifdef __INTEL_COMPILER
         errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

         sim_triangle(outlier_rowidx, outlier_colidx) = &
            factor * dev_facs(outlier_colidx - 1) * sim_triangle(outlier_rowidx, outlier_colidx - 1) + &
            r(1) * sigmas(outlier_colidx - 1) * sqrt(sim_triangle(outlier_rowidx, outlier_colidx - 1))

#elif defined __GFORTRAN__
         sim_triangle(outlier_rowidx, outlier_colidx) = &
            factor * dev_facs(outlier_colidx - 1) * sim_triangle(outlier_rowidx, outlier_colidx - 1) + &
            r(1) * sigmas(outlier_colidx - 1) * sqrt(sim_triangle(outlier_rowidx, outlier_colidx - 1))

#endif

         if (outlier_colidx < n_dev) then
            do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx

#ifdef __INTEL_COMPILER
               errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

               sim_triangle(outlier_rowidx, j) = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + &
                  r(1) * sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
#elif defined __GFORTRAN__

               sim_triangle(outlier_rowidx, j) = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + &
                  r(1) * sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))

#endif

            end do
         end if

      else if (dist == "gamma") then

#ifdef __INTEL_COMPILER
         method = VSL_RNG_METHOD_GAMMA_GNORM
#endif

         do j = 2, n_dev
            do i = 1, n_dev + 1 - j
               if (i == outlier_rowidx) cycle

               gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
               gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

#ifdef __INTEL_COMPILER
               errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

               sim_triangle(i, j) = r(1)

#elif __GFORTRAN__
               sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif

            end do
         end do

         if (outlier_colidx > 2) then
            do j = 2, outlier_colidx - 1
               gamma_shape = dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
               gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

#ifdef __INTEL_COMPILER
               errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

               sim_triangle(outlier_rowidx, j) = r(1)

#elif __GFORTRAN__
               sim_triangle(outlier_rowidx, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif

            end do
         end if

         gamma_shape = dev_facs(outlier_colidx - 1)**2 * sim_triangle(outlier_rowidx, outlier_colidx - 1) / sigmas(outlier_colidx - 1)**2
         gamma_rate = dev_facs(outlier_colidx - 1) / sigmas(outlier_colidx - 1)**2

#ifdef __INTEL_COMPILER
         errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

         sim_triangle(outlier_rowidx, outlier_colidx) = r(1)

#elif __GFORTRAN__
         sim_triangle(outlier_rowidx, outlier_colidx) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif

         if (outlier_colidx < n_dev) then
            do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
               gamma_shape = dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
               gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

#ifdef __INTEL_COMPILER
               errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

               sim_triangle(outlier_rowidx, j) = r(1)

#elif __GFORTRAN__
               sim_triangle(outlier_rowidx, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif

            end do
         end if

      end if

#ifdef __INTEL_COMPILER
      errcode = vsldeletestream(stream)
#endif

   end function single_outlier

   function calendar_outlier(outlier_diagidx, factor, triangle, dev_facs, sigmas, dist) result(sim_triangle)

      integer, intent(in) :: outlier_diagidx
      real(dp), intent(in):: factor
      real(dp), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
      character(:), allocatable, intent(in) :: dist

      integer :: i, j, n_dev, n_cols
      real(dp), allocatable :: sim_triangle(:, :)
      real(dp) :: gamma_shape, gamma_rate

#ifdef __INTEL_COMPILER
      type(vsl_stream_state) :: stream
      integer(kind=4) :: errcode
      integer :: brng, method, sd, n
      real(dp) :: r(1)

      brng = VSL_BRNG_MT19937
      sd = irand()
      errcode = vslnewstream(stream, brng,  sd)
#endif

      n_dev = size(triangle, dim=1)

      allocate(sim_triangle(n_dev, n_dev), source=0._dp)
      sim_triangle(:, 1) = triangle(:, 1)

      do i = 1, n_dev
         n_cols = n_dev + 2 - outlier_diagidx - i
         if (n_cols <= 1) then
            do j = 2, n_dev + 1 - i
               sim_triangle(i, j) = triangle(i, j)
            end do
         else
            do j = 2, n_cols - 1
               sim_triangle(i, j) = triangle(i, j)
            end do

            if (dist == "normal") then

#ifdef __INTEL_COMPILER
               method = VSL_RNG_METHOD_GAUSSIAN_ICDF

               errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

               sim_triangle(i, n_cols) = factor * dev_facs(n_cols - 1) * sim_triangle(i, n_cols - 1) + sigmas(n_cols - 1) * sqrt(sim_triangle(i, n_cols - 1)) * r(1)

#elif defined __GFORTRAN__

               sim_triangle(i, n_cols) = factor * dev_facs(n_cols - 1) * sim_triangle(i, n_cols - 1) + sigmas(n_cols - 1) * sqrt(sim_triangle(i, n_cols - 1)) * norm()

#endif
               do j = n_cols + 1, n_dev + 1 - i
#ifdef __INTEL_COMPILER
                  errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

                  sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(i, j - 1)) * r(1)

#elif defined __GFORTRAN__

                  sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(i, j - 1)) * norm()

#endif
               end do
            else if (dist == "gamma") then

               gamma_shape = factor * dev_facs(n_cols - 1)**2 * sim_triangle(i, n_cols - 1) / sigmas(n_cols - 1)**2
               gamma_rate = factor * dev_facs(n_cols - 1) / sigmas(n_cols - 1)**2

#ifdef __INTEL_COMPILER

               method = VSL_RNG_METHOD_GAMMA_GNORM
               errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

               sim_triangle(i, n_cols) = r(1)

#elif __GFORTRAN__

               sim_triangle(i, n_cols) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif
               do j = n_cols + 1, n_dev + 1 - i

                  gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
                  gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

#ifdef __INTEL_COMPILER

                  errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

                  sim_triangle(i, j) = r(1)

#elif __GFORTRAN__

                  sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif
               end do
            end if
         end if
      end do

#ifdef __INTEL_COMPILER
      errcode = vsldeletestream(stream)
#endif

   end function calendar_outlier

   function origin_outlier(outlier_rowidx, factor, triangle, dev_facs, sigmas, dist) result(sim_triangle)

      integer, intent(in):: outlier_rowidx
      real(dp), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
      character(:), allocatable, intent(in) :: dist
      real(dp), intent(in) :: factor

      real(dp) :: gamma_shape, gamma_rate
      real(dp), allocatable:: sim_triangle(:, :)

      integer :: n_dev, i, j

#ifdef __INTEL_COMPILER
      type(vsl_stream_state) :: stream
      integer(kind=4) :: errcode
      integer :: brng, method, sd, n
      real(dp) :: r(1)

      brng = VSL_BRNG_MT19937
      sd = irand()
      errcode = vslnewstream(stream, brng,  sd)
#endif

      n_dev = size(triangle, dim=1)
      sim_triangle = triangle

      do j = 2, n_dev + 1 - outlier_rowidx

         if (dist == "normal") then

#ifdef __INTEL_COMPILER
            method = VSL_RNG_METHOD_GAUSSIAN_ICDF

            errcode = vdrnggaussian(method, stream, 1, r, 0._dp, 1._dp)

            sim_triangle(outlier_rowidx, j) = factor * dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1)) * r(1)

#elif defined __GFORTRAN__

            sim_triangle(outlier_rowidx, j) = factor * dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1)) * norm()

#endif

         else if (dist == "gamma") then

            gamma_shape = factor * dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
            gamma_rate = factor * dev_facs(j - 1) / sigmas(j - 1)**2

#ifdef __INTEL_COMPILER

            method = VSL_RNG_METHOD_GAMMA_GNORM
            errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._dp, 1/gamma_rate)

            sim_triangle(outlier_rowidx, j) = r(1)

#elif __GFORTRAN__

            sim_triangle(outlier_rowidx, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), dp)

#endif
         end if
      end do

#ifdef __INTEL_COMPILER
      errcode = vsldeletestream(stream)
#endif
   end function origin_outlier

end subroutine reserve_sim


