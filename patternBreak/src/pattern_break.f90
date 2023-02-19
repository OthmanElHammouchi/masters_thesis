module pattern_break

   use, intrinsic :: iso_c_binding

   implicit none

   integer, parameter :: NORMAL = 1
   integer, parameter :: GAMMA = 2
   integer, parameter :: CONDITIONAL = 1
   integer, parameter :: UNCONDITIONAL = 2
   integer, parameter :: RAW = 1
   integer, parameter :: SCALED = 2
   integer, parameter :: PARAMETRIC = 3
   integer, parameter :: SINGLE = 1
   integer, parameter :: CALENDAR = 2
   integer, parameter :: ORIGIN = 3

   include "c_interface.h"

contains

   subroutine mack_sim(triangle, n_boot, n_dev, config, results)

      integer(c_int), intent(in), value :: n_dev, n_boot
      real(c_double), intent(in) :: config(config_spec(1), config_spec(2)), triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: results(n_boot*config_spec(1), config_spec(2) + 1)

      integer(c_int) :: i, j, k, i_sim, n_rows, n_config
      integer(c_int) :: outlier_rowidx, outlier_colidx, outlier_diagidx, excl_diagidx, excl_rowidx

      integer(c_int), allocatable :: excl_resids(:, :)
      real(c_double) :: indiv_dev_facs(n_dev - 1, n_dev - 1), dev_facs(n_dev - 1), sigmas(n_dev - 1)
      real(c_double) :: reserve(n_boot), init_col(n_dev), sim_triangle(n_dev, n_dev)
      real(c_double) :: factor

      real(c_double) :: pct_completed

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

         pct_completed = 100 * (i_sim / n_config)
         call dblepr1(pct_completed)

         if (config_type == SINGLE) then

            resids_type = int(config(i_sim, 6))
            boot_type = int(config(i_sim, 7))
            dist = int(config(i_sim, 8))

            if (resids_type == "parametric" .and. dist == "gamma") then
               call rexit("Parametric residuals not supported for gamma distribution.")
            end if

            allocate(excl_resids(1, 2))

            excl_resids(1, :) = int(config(i_sim, 4:5))
            outlier_rowidx = int(config(i_sim, 1))
            outlier_colidx = int(config(i_sim, 2))
            factor = config(i_sim, 3)

            sim_triangle = single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist)

            reserve = reserve_boot(sim_triangle, n_dev, n_boot, resids_type, boot_type, dist, excl_resids, log_unit)

            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:config_spec(2)) = transpose(spread(config(i_sim, :), 2, n_boot))
            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), config_spec(2) + 1) = reserve

         else if (config_type == CALENDAR) then

            resids_type = int(config(i_sim, 4))
            boot_type = int(config(i_sim, 5))
            dist = int(config(i_sim, 6))

            if (resids_type == "parametric" .and. dist == "gamma") then
               call rexit("Parametric residuals not supported for gamma distribution.")
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

            reserve = reserve_boot(sim_triangle, n_dev, n_boot, resids_type, boot_type, dist, excl_resids, log_unit)

            deallocate(excl_resids)

            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:config_spec(2)) = transpose(spread(config(i_sim, :), 2, n_boot))
            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), config_spec(2) + 1) = reserve

         else if (config_type == ORIGIN) then

            resids_type = int(config(i_sim, 4))
            boot_type = int(config(i_sim, 5))
            dist = int(config(i_sim, 6))

            if (resids_type == "parametric" .and. dist == "gamma") then
               call rexit("Parametric residuals not supported for gamma distribution.")
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

            reserve = reserve_boot(sim_triangle, n_dev, n_boot, resids_type, boot_type, dist, excl_resids, log_unit)

            deallocate(excl_resids)

            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:config_spec(2)) = transpose(spread(config(i_sim, :), 2, n_boot))
            results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), config_spec(2) + 1) = reserve

         end if

      end do

   end subroutine mack_sim

   ! Subroutine implementing bootstrap of Mack's model for claims reserving.
   subroutine mack_boot(triangle, n_dev, n_boot, resids_type, boot_type, dist, reserve, excl_resids)

      integer(c_int), intent(in) :: n_boot, n_dev, dist, resids_type, boot_type
      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      integer(c_int), intent(in), optional :: excl_resids(:, :)

      integer(c_int) :: i, j, k, i_diag, i_boot, n_rows, n_resids, n_excl_resids ! bookkeeping variables

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

      if (present(excl_resids)) then

         n_excl_resids = size(excl_resids, dim=1)

         do i = 1, n_excl_resids
            triangle_mask(excl_resids(i, 1), excl_resids(i, 2)) = .false.
            resids_mask(excl_resids(i, 1), excl_resids(i, 2) - 1) = .false.
         end do

      end if

      do j = 1, n_dev - 1

         n_rows = n_dev - j

         indiv_dev_facs(1:n_rows, j) = triangle(1:n_rows, j + 1) / triangle(1:n_rows, j)

         if (resids_type == PARAMETRIC) then

            dev_facs(j) = sum(pack(triangle(1:n_rows, j + 1), triangle_mask(1:n_rows, j + 1))) / sum(pack(triangle(1:n_rows, j), triangle_mask(1:n_rows, j)))

            if (j < n_dev - 1) then
               sigmas(j) = sqrt(sum(pack(triangle(1:n_rows, j), triangle_mask(1:n_rows, j)) * (pack(indiv_dev_facs(1:n_rows, j), triangle_mask(1:n_rows, j)) - dev_facs(j)) ** 2) / (count(triangle_mask(1:n_rows, j)) - 1))
            else
               sigmas(j) = sqrt(min(sigmas(j - 1) ** 2, sigmas(j - 2) ** 2, sigmas(j - 1) ** 4 / sigmas(j - 2) ** 2))
            end if

            resids(1:n_rows, j) = (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) * sqrt(triangle(1:n_rows, j)) / sigmas(j)

         else

            dev_facs(j) = sum(triangle(1:n_rows, j + 1)) / sum(triangle(1:n_rows, j))

            if (j < n_dev - 1) then
               sigmas(j) = sqrt(sum(triangle(1:n_rows, j) * (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) ** 2) / (n_rows - 1))
            else
               sigmas(j) = sqrt(min(sigmas(j - 1) ** 2, sigmas(j - 2) ** 2, sigmas(j - 1) ** 4 / sigmas(j - 2) ** 2))
            end if

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

      n_resids = ((n_dev - 1)**2 + (n_dev - 1))/2 - 1    ! Discard residual from upper right corner.
      if (present(excl_resids)) n_resids = n_resids - n_excl_resids

      allocate(flat_resids(n_resids))

      if (resids_type /= PARAMETRIC) flat_resids = pack(resids, resids_mask)

      main_loop: do i_boot = 1, n_boot
         if (resids_type == PARAMETRIC) then

            do i = 1, n_dev - 1
               do j = 1, n_dev - 1
                  do k = 1, n_dev - 1
                     resids_boot(i, j, k) = rnorm(0._c_double, 1._c_double)
                  end do
               end do
            end do

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

         if (boot_type == CONDITIONAL) then
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
                  call rexit("Bad sigma.")
                  cycle main_loop
               end if
            end do

         else if (boot_type == UNCONDITIONAL) then
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
                  call rexit("Bad sigma.")
               end if

            end do

         end if

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

                  if (shape <= tiny(1.0) .or. isnan(shape)) then
#ifdef DEBUG
                     call log_neg_gamma(log_unit, triangle_boot, resids_type, boot_type, dist, sigmas_boot, indiv_dev_facs_boot, &
                        dev_facs_boot, resids_boot)
#endif
                     cycle main_loop
                  end if

                  scale = sigmas_boot(j - 1, i - 1) ** 2 / dev_facs_boot(j - 1, i - 1)
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

   ! Entry point for the C wrapper. Orders the arguments properly and omits excl_resids argument.
   subroutine mack_boot_centry(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve) bind(C, name='mack_boot_')

      integer(c_int), intent(in), value :: n_boot, n_dev, dist, resids_type, boot_type
      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: reserve(n_boot)

      call mack_boot(triangle, n_dev, n_boot, resids_type, boot_type, dist, reserve)

   end subroutine mack_boot_centry


   subroutine glm_boot(triangle, n_dev, n_boot, reserve, excl_resids)

      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: reserve(n_boot)

      integer(c_int), intent(in) :: n_dev, n_boot

      real(c_double), intent(in), optional :: excl_resids(:, :)

      integer(c_int) :: n_pts, n_covs, n_pts_boot
      real(c_double), allocatable :: y(:), X(:, :), b(:), b_boot(:), X_lower(:, :), y_fit_boot(:)
      real(c_double), allocatable :: y_fit(:), triangle_fit(:, :), triangle_boot(:, :), y_boot(:)
      real(c_double), allocatable :: resids(:, :), resids_boot(:, :), flat_resids(:)
      logical(c_bool), allocatable :: resids_mask(:, :)

      integer(c_int) :: i, j, k, i_boot, info

      n_pts = (n_dev**2 + n_dev) / 2
      n_covs = 2*n_dev - 1

      n_pts_boot = n_dev ** 2 - n_pts

      allocate(y(n_pts), X(n_pts, n_covs), X_lower(n_pts_boot, n_covs), b(n_covs), source=0._c_double)

      allocate(resids_boot(n_pts, n_pts), b_boot(n_covs), triangle_boot(n_pts, n_pts), y_fit_boot(n_pts_boot), source = 0._c_double)

      allocate(y_boot(n_pts), y_fit(n_pts), triangle_fit(n_pts, n_pts), source=0._c_double)

      allocate(resids(n_dev, n_dev), flat_resids(n_pts), source=0._c_double)

      allocate(resids_mask(n_dev, n_dev))

      ! Set up feature matrix and response vector.
      X(:, 1) = 1._c_double

      k = 1
      do i = 1, n_dev
         do j = 1, n_dev + 1 - i
            if (i /= 1) X(k, i) = 1
            if (j /= 1) X(k, n_dev + j - 1) = 1
            y(k) = triangle(i, j)
            k = k + 1
         end do
      end do

      y(1) = triangle(1, 1)

      call poisson_fit(X, y, n_pts, n_covs, b)

      ! if (any(isnan(b))) stop "b became NaN"

      y_fit = exp(matmul(X, b))

      k = 1
      do i = 1, n_dev
         do j = 1, n_dev + 1 - i
            triangle_fit(i, j) = y_fit(k)
            k = k + 1
         end do
      end do

      do i = 1, n_dev
         do j = 1, n_dev + 1 - i
            resids(i, j) = (triangle(i, j) - triangle_fit(i, j)) / sqrt(triangle_fit(i, j))
         end do
      end do

      resids_mask = .true.

      flat_resids = pack(resids, resids_mask)

      X_lower(:, 1) = 1._c_double

      k = 1
      do i = 2, n_dev
         do j = n_dev + 1 - i + 1, n_dev
            X_lower(k, i) = 1
            if (j /= 1) X_lower(k, n_dev + j - 1) = 1
            k = k + 1
         end do
      end do

      call GetRNGstate()

      do i_boot = 1, n_boot

         do i = 1, n_dev
            do j = 1, n_dev + 1 - i
               resids_boot(i, j) = flat_resids(1 + int(n_pts * rand()))
               triangle_boot(i, j) = resids_boot(i, j) * sqrt(triangle_fit(i, j)) + triangle_fit(i, j)
            end do
         end do

         k = 1
         do i = 1, n_dev
            do j = 1, n_dev + 1 - i
               y_boot(k) = triangle_boot(i, j)
               k = k + 1
            end do
         end do

         call poisson_fit(X, y_boot, n_pts, n_covs, b_boot)

         y_fit_boot = exp(matmul(X_lower, b_boot))

         do i = 1, n_pts_boot
            y_fit_boot(i) = rpois(y_fit_boot(i))
         end do   
         
         reserve(i_boot) = sum(y_fit_boot)

      end do

      call PutRNGstate()

   end subroutine glm_boot

   subroutine poisson_fit(X, y, N, p, b)

      integer(c_int), intent(in) :: N, p
      real(c_double), intent(in) :: X(N, p), y(N)
      real(c_double), intent(inout) :: b(p)

      real(c_double) :: diff, eps
      real(c_double) :: IPIV(p, p), rhs(p), lhs(p, p)
      real(c_double) :: W(N, N), z(N), eta(N)

      integer(c_int) :: i, info

      ! Initialise coefficients.
      b = spread(0._c_double, 1, p)
      b(1) = log(sum(y) / N)

      ! Apply IRLS to fit Poisson GLM.
      diff = 1E6
      eps = 1E-6
      do while (diff > eps)

         eta = matmul(X, b)

         W = 0

         do i = 1, N
            W(i, i) = exp(eta(i))
         end do

         z = eta + exp(-eta)*y - 1

         lhs = matmul(matmul(transpose(X), W), X)
         rhs = matmul(matmul(transpose(X), W), z)

         call dgesv(p, 1, lhs, p, IPIV, rhs, p, info)

         diff = norm2(b - rhs)

         b = rhs

      end do

   end subroutine poisson_fit

   subroutine glm_boot_centry(n_dev, triangle, n_boot, reserve) bind(C, name='glm_boot_')

      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: reserve(n_boot)

      integer(c_int), intent(in), value :: n_dev, n_boot

      call glm_boot(triangle, n_dev, n_boot, reserve)

   end subroutine glm_boot_centry


end module pattern_break
