module mack

   use, intrinsic :: iso_c_binding
   use constants
   use helpers
   use rng_fort_interface
   use print_fort_interface

   implicit none

contains

   subroutine mack_sim(n_dev, triangle, n_config, m_config, config, type, n_boot, results) bind(C, name='mack_sim_')

      integer(c_int), intent(in), value :: n_dev, n_boot, n_config, m_config, type
      real(c_double), intent(in) :: config(n_config, m_config)
      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(out) :: results(n_boot * n_config, m_config + 1)

      integer(c_int) :: i, j, k, i_sim, counter, n_rows
      integer(c_int) :: outlier_rowidx
      integer(c_int) :: outlier_colidx
      integer(c_int) :: outlier_diagidx
      integer(c_int) :: excl_diagidx
      integer(c_int) :: excl_rowidx

      integer(c_int), allocatable :: excl_resids(:, :)
      real(c_double) :: indiv_dev_facs(n_dev - 1, n_dev - 1)
      real(c_double) :: dev_facs(n_dev - 1)
      real(c_double) :: sigmas(n_dev - 1)
      real(c_double) :: reserve(n_boot) 
      real(c_double) :: init_col(n_dev)
      real(c_double) :: triangle_sim(n_dev, n_dev)
      real(c_double) :: factor

      integer(c_int) :: resids_type, boot_type, dist

      real(c_double) :: pct_progress
      real(c_double) :: progress
      character(len=999, kind=c_char) :: progress_str

      init_col = triangle(:, 1)

      call mack_fit(triangle, dev_facs, sigmas)

      counter = 0

      !$OMP PARALLEL DO SHARED(counter) PRIVATE(excl_resids)
      do i_sim = 1, n_config

         if (type == SINGLE) then

            resids_type = int(config(i_sim, 6))
            boot_type = int(config(i_sim, 7))
            dist = int(config(i_sim, 8))

            allocate(excl_resids(1, 2))
            excl_resids(1, :) = int(config(i_sim, 4:5))

            outlier_rowidx = int(config(i_sim, 1))
            outlier_colidx = int(config(i_sim, 2))
            factor = config(i_sim, 3)

            triangle_sim = single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist)

            call mack_boot(n_dev, triangle_sim, resids_type, boot_type, dist, n_boot, reserve, excl_resids)

            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve

            deallocate(excl_resids)

         else if (type == CALENDAR) then

            resids_type = int(config(i_sim, 4))
            boot_type = int(config(i_sim, 5))
            dist = int(config(i_sim, 6))

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

            triangle_sim = calendar_outlier(outlier_diagidx, factor, triangle, dev_facs, sigmas, dist)

            call mack_boot(n_dev, triangle_sim, resids_type, boot_type, dist, n_boot, reserve, excl_resids)

            deallocate(excl_resids)

            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve

         else if (type == ORIGIN) then

            resids_type = int(config(i_sim, 4))
            boot_type = int(config(i_sim, 5))
            dist = int(config(i_sim, 6))

            outlier_rowidx = int(config(i_sim, 1))
            excl_rowidx = int(config(i_sim, 3))
            factor = config(i_sim, 2)

            allocate(excl_resids(n_dev - excl_rowidx, 2), source=0)

            k = 1
            do j = 2, n_dev + 1 - excl_rowidx
               excl_resids(k, :) = [excl_rowidx, j]
               k = k + 1
            end do

            triangle_sim = origin_outlier(outlier_rowidx, factor, triangle, dev_facs, sigmas, dist)

            call mack_boot(n_dev, triangle_sim, resids_type, boot_type, dist, n_boot, reserve, excl_resids)

            deallocate(excl_resids)

            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve

         end if

         progress = real(counter, kind=c_double) / real(n_config, kind=c_double)

         if (abs(mod(progress, 0.005_c_double)) <= 1e-4_c_double) then
            pct_progress = 100*progress
            write(progress_str, "('Progress: ', f5.2)") pct_progress
            call Rprintf(achar(13) // c_null_char)
            call Rprintf(repeat(" ", 118) // c_null_char)
            call Rprintf(achar(13) // c_null_char)
            call Rprintf(trim(progress_str) // c_null_char )
            call Rprintf(achar(13) // c_null_char)
            call R_FlushConsole()
         end if

         counter = counter + 1

      end do
      !$OMP END PARALLEL DO

   end subroutine mack_sim

   ! Subroutine implementing bootstrap of Mack's model for claims reserving.
   subroutine mack_boot(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve, excl_resids)

      integer(c_int), intent(in) :: n_boot, n_dev, dist, resids_type, boot_type
      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      integer(c_int), intent(in), optional :: excl_resids(:, :) ! Long format list of points to exclude.

      integer(c_int) :: i, j, k, i_diag, i_boot, n_rows, n_resids, n_excl_resids ! Bookkeeping variables.

      real(c_double) :: dev_facs(n_dev - 1)
      real(c_double) :: sigmas(n_dev - 1)
      real(c_double) :: latest(n_dev)
      real(c_double) :: indiv_dev_facs(n_dev - 1, n_dev - 1)
      real(c_double) :: resids(n_dev - 1, n_dev - 1)
      real(c_double) :: scale_facs(n_dev - 1, n_dev - 1)
      real(c_double) :: resampled_triangle(n_dev, n_dev, n_dev - 1)

      real(c_double) :: resids_boot(n_dev - 1, n_dev - 1, n_dev - 1)
      real(c_double) :: indiv_dev_facs_boot(n_dev - 1, n_dev - 1, n_dev - 1)
      real(c_double) :: dev_facs_boot(n_dev - 1, n_dev - 1)
      real(c_double) :: sigmas_boot(n_dev - 1, n_dev - 1)
      real(c_double) :: triangle_boot(n_dev, n_dev)

      logical(c_bool) :: triangle_mask(n_dev, n_dev)
      logical(c_bool) :: resids_mask(n_dev - 1, n_dev - 1)

      real(c_double) :: reserve(n_boot)
      real(c_double), allocatable :: flat_resids(:)

      real(c_double) :: mean, sd ! Normal distribution parameters.
      real(c_double) :: shape, scale ! Gamma distribution parameters.

      Call GetRNGstate()

      triangle_mask = .true.
      resids_mask = .true.

      n_resids = ((n_dev - 1)**2 + (n_dev - 1))/2 - 1    ! Discard residual from upper right corner.

      if (present(excl_resids)) then

         n_excl_resids = size(excl_resids, dim=1)
         n_resids = n_resids - n_excl_resids

         do i = 1, n_excl_resids
            triangle_mask(excl_resids(i, 1), excl_resids(i, 2)) = .false.
            resids_mask(excl_resids(i, 1), excl_resids(i, 2) - 1) = .false.
         end do

      end if

      if (resids_type == PARAMETRIC) then
         call mack_fit(triangle, dev_facs, sigmas, triangle_mask=triangle_mask)
      else
         call mack_fit(triangle, dev_facs, sigmas, resids=resids, resids_type=resids_type)
         flat_resids = pack(resids, resids_mask)
      end if


      main_loop: do i_boot = 1, n_boot

         if (resids_type /= PARAMETRIC) then

            do i = 1, n_dev - 1
               do j = 1, n_dev - 1
                  do k = 1, n_dev - 1
                     resids_boot(i, j, k) = flat_resids(1 + int(n_resids * rand()))
                  end do
               end do
            end do

         end if

         ! Parameter error.
         if (boot_type == CONDITIONAL) then
            if (resids_type == PARAMETRIC) then
               do k = 1, n_dev - 1
                  resampled_triangle(:, 1, k) = triangle(:, 1)
               end do

               do k = 1, n_dev - 1
                  do j = 2, n_dev
                     do i = 1, n_dev + 1 - j

                        if (dist == NORMAL) then

                           mean = dev_facs(j - 1) * triangle(i, j - 1)
                           sd = sigmas(j - 1) * sqrt(triangle(i, j - 1))

                           resampled_triangle(i, j, k) = rnorm(mean, sd)

                        else if (dist == GAMMA) then

                           shape = (dev_facs(j - 1)**2 * triangle(i, j - 1)) / sigmas(j - 1) **2
                           scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)

                           resampled_triangle(i, j, k) = rgamma(shape, scale)

                        end if

                     end do
                  end do

                  call mack_fit(resampled_triangle(:, :, k), dev_facs_boot(:, k), sigmas_boot(:, k))

               end do

            else

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

                  if (any(sigmas_boot < 0) .or. any(isnan(sigmas_boot))) cycle main_loop
               end do

            end if

         else if (boot_type == UNCONDITIONAL) then
            if (resids_type == PARAMETRIC) then
               do k = 1, n_dev - 1
                  resampled_triangle(:, 1, k) = triangle(:, 1)
               end do

               do k = 1, n_dev - 1
                  do j = 2, n_dev
                     do i = 1, n_dev + 1 - j

                        if (dist == NORMAL) then

                           mean = dev_facs(j - 1) * resampled_triangle(i, j - 1, k)
                           sd = sigmas(j - 1) * sqrt(resampled_triangle(i, j - 1, k))

                           resampled_triangle(i, j, k) = rnorm(mean, sd)

                        else if (dist == GAMMA) then

                           shape = (dev_facs(j - 1)**2 * resampled_triangle(i, j - 1, k)) / sigmas(j - 1) **2
                           scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)

                           resampled_triangle(i, j, k) = rgamma(shape, scale)

                        end if

                     end do
                  end do

                  call mack_fit(resampled_triangle(:, :, k), dev_facs_boot(:, k), sigmas_boot(:, k))

               end do

            else

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

                  if (any(sigmas_boot < 0) .or. any(isnan(sigmas_boot))) cycle main_loop

               end do

            end if

         end if

         ! Process error.
         triangle_boot = triangle

         if (dist == NORMAL) then

            do i_diag = 1, n_dev - 1
               do i = i_diag + 1, n_dev

                  j = n_dev + i_diag + 1 - i

                  mean = dev_facs_boot(j - 1, i - 1) * triangle_boot(i, j - 1)
                  sd = sigmas_boot(j - 1, i - 1) * sqrt(triangle_boot(i, j - 1))

                  triangle_boot(i, j) = rnorm(mean, sd)

                  if (triangle_boot(i, j) <= 0) then
                     cycle main_loop
                  end if

               end do
            end do

            do j = 1, n_dev
               latest(j) = triangle_boot(n_dev + 1 - j, j)
            end do

            reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)

         else if (dist == GAMMA) then

            do i_diag = 1, n_dev - 1
               do i = i_diag + 1, n_dev

                  j = n_dev + i_diag + 1 - i

                  shape = (dev_facs_boot(j - 1, i - 1)**2 * triangle_boot(i, j - 1)) / sigmas_boot(j - 1, i - 1) **2
                  scale = sigmas_boot(j - 1, i - 1) ** 2 / dev_facs_boot(j - 1, i - 1)

                  if (shape <= tiny(1.0) .or. isnan(shape)) cycle main_loop

                  triangle_boot(i, j) = rgamma(shape, scale)

               end do
            end do

            do j = 1, n_dev
               latest(j) = triangle_boot(n_dev + 1 - j, j)
            end do

            reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)
         end if

      end do main_loop

      Call PutRNGstate()

   end subroutine mack_boot

   subroutine mack_fit(triangle, dev_facs, sigmas, resids, resids_type, triangle_mask)

      real(c_double), intent(in) :: triangle(:, :)
      real(c_double), intent(out) :: dev_facs(:)
      real(c_double), intent(out) :: sigmas(:)
      real(c_double), optional, intent(out) :: resids(:, :)
      integer(c_int), optional, intent(in) :: resids_type
      logical(c_bool), optional, intent(in) :: triangle_mask(:, :)

      integer(c_int) :: i, j, n_rows, n_dev
      logical(c_bool), allocatable :: triangle_mask_(:, :)
      real(c_double), allocatable :: indiv_dev_facs(:, :)
      real(c_double), allocatable :: scale_facs(:, :)


      n_dev = size(triangle, 1)

      allocate(indiv_dev_facs(n_dev, n_dev), source=0._c_double)
      allocate(scale_facs(n_dev, n_dev), source=0._c_double)

      if (present(triangle_mask)) then
         triangle_mask_ = triangle_mask
      else
         allocate(triangle_mask_(n_dev, n_dev))
         triangle_mask_ = .true.
      end if

      do j = 1, n_dev - 1

         n_rows = n_dev - j

         indiv_dev_facs(1:n_rows, j) = triangle(1:n_rows, j + 1) / triangle(1:n_rows, j)

         dev_facs(j) = sum(pack(triangle(1:n_rows, j + 1), triangle_mask_(1:n_rows, j + 1))) / sum(pack(triangle(1:n_rows, j), triangle_mask_(1:n_rows, j)))

         if (j < n_dev - 1) then
            sigmas(j) = sqrt(sum(pack(triangle(1:n_rows, j), triangle_mask_(1:n_rows, j)) * (pack(indiv_dev_facs(1:n_rows, j), triangle_mask_(1:n_rows, j)) - dev_facs(j)) ** 2) / (count(triangle_mask_(1:n_rows, j)) - 1))
         else
            sigmas(j) = sqrt(min(sigmas(j - 1) ** 2, sigmas(j - 2) ** 2, sigmas(j - 1) ** 4 / sigmas(j - 2) ** 2))
         end if

         if (present(resids) .and. present(resids_type)) then

            if (resids_type == RAW) then

               resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * sqrt(triangle(1:n_rows, j)) / sigmas(j)

            else if (resids_type == SCALED) then

               if (j < n_dev - 1) then
                  scale_facs(1:n_rows, j) = sqrt(1 - triangle(1:n_rows, j) / sum(triangle(1:n_rows, j)))
               else
                  scale_facs(1:n_rows, j) = 1
               end if

               resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * &
                  sqrt(triangle(1:n_rows, j)) / (sigmas(j) * scale_facs(1:n_rows, j))

            end if

         end if

      end do

   end subroutine mack_fit


   ! Entry point for C wrapper, to omit excl_resids argument.
   subroutine mack_boot_centry(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve) bind(C, name='mack_boot_')

      integer(c_int), intent(in), value :: n_boot, n_dev, dist, resids_type, boot_type
      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: reserve(n_boot)

      call mack_boot(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve)

   end subroutine mack_boot_centry

end module mack
