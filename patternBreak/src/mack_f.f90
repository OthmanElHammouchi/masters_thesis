module mack

  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  use global
  use helpers
  use interface

  implicit none

contains

  subroutine mack_sim_f(n_dev, triangle, n_config, m_config, config, type, n_boot, results) bind(c)
    integer(c_int), intent(in), value :: n_dev, n_boot, n_config, m_config, type
    real(c_double), intent(in) :: config(n_config, m_config)
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: results(n_boot * n_config, m_config + 1)

    integer(c_int) :: i, j, k, i_sim, counter, n_rows, inc, i_thread, n_threads
    integer(c_int) :: outlier_rowidx, outlier_colidx, outlier_diagidx
    integer(c_int) :: excl_diagidx, excl_rowidx, n_failures

    integer(c_int), allocatable :: excl_resids(:, :)
    real(c_double), allocatable :: reserve(:)
    real(c_double) :: indiv_dev_facs(n_dev - 1, n_dev - 1), dev_facs(n_dev - 1), sigmas(n_dev - 1)
    real(c_double) :: init_col(n_dev), triangle_sim(n_dev, n_dev), factor
    integer(c_int) :: resids_type, boot_type, dist

    integer(c_int) :: status

    type(c_ptr) :: pgbar

    init_col = triangle(:, 1)

    call fit(triangle, dev_facs, sigmas)

    n_threads = init_omp()
    rng = init_rng(n_threads, 42)
    pgbar = pgbar_create(n_config, 1)

    results = 0

    !$omp parallel num_threads(n_threads) &
    !$omp& private(reserve, i_thread, i_sim, resids_type, boot_type, dist, outlier_rowidx, outlier_colidx, &
    !$omp& excl_resids, triangle_sim) &
    !$omp& firstprivate(m_config, n_config, n_boot, n_dev, init_col, sigmas, dev_facs, rng) shared(config, results, pgbar)
    allocate(reserve(n_boot))
    reserve = 0

    i_thread = omp_get_thread_num()

    if (type == SINGLE) then
      allocate(excl_resids(1, 2))
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_config

        call pgbar_incr(pgbar)

        resids_type = int(config(i_sim, 6))
        boot_type = int(config(i_sim, 7))
        dist = int(config(i_sim, 8))

        outlier_rowidx = int(config(i_sim, 1))
        outlier_colidx = int(config(i_sim, 2))
        factor = config(i_sim, 3)

        triangle_sim = single_outlier_mack(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist, rng, i_thread)

        excl_resids(1, :) = int(config(i_sim, 4:5))
        call mack_boot(n_dev, triangle_sim, resids_type, boot_type, dist, n_boot, reserve, excl_resids, status)

        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve
      end do
      !$omp end do

    else if (type == CALENDAR) then
      allocate(excl_resids(n_dev - excl_diagidx, 2))
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_config

        call pgbar_incr(pgbar)

        resids_type = int(config(i_sim, 4))
        boot_type = int(config(i_sim, 5))
        dist = int(config(i_sim, 6))

        outlier_diagidx = int(config(i_sim, 1))
        excl_diagidx = int(config(i_sim, 3))
        factor = config(i_sim, 2)

        excl_resids = 0
        k = 1
        do j = 2, n_dev
          i = n_dev + 2 - excl_diagidx - j
          if (i <= 0) cycle
          excl_resids(k, :) = [i, j]
          k = k + 1
        end do

        triangle_sim = calendar_outlier_mack(outlier_diagidx, factor, triangle, dev_facs, sigmas, dist, rng, i_thread)

        call mack_boot(n_dev, triangle_sim, resids_type, boot_type, dist, n_boot, reserve, excl_resids, status)

        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve
      end do
      !$omp end do

    else if (type == ORIGIN) then
      allocate(excl_resids(n_dev - excl_rowidx, 2), source=0)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_config

        call pgbar_incr(pgbar)

        resids_type = int(config(i_sim, 4))
        boot_type = int(config(i_sim, 5))
        dist = int(config(i_sim, 6))

        outlier_rowidx = int(config(i_sim, 1))
        excl_rowidx = int(config(i_sim, 3))
        factor = config(i_sim, 2)

        k = 1
        do j = 2, n_dev + 1 - excl_rowidx
          excl_resids(k, :) = [excl_rowidx, j]
          k = k + 1
        end do

        triangle_sim = origin_outlier_mack(outlier_rowidx, factor, triangle, dev_facs, sigmas, dist, rng, i_thread)

        call mack_boot(n_dev, triangle_sim, resids_type, boot_type, dist, n_boot, reserve, excl_resids, status)

        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve
      end do
      !$omp end do
    end if
    deallocate(excl_resids)
    deallocate(reserve)
    !$omp end parallel
  end subroutine mack_sim_f

! Subroutine implementing bootstrap of Mack's model for claims reserving.
  subroutine mack_boot(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve, excl_resids, status)

    integer(c_int), intent(in) :: n_boot, n_dev, dist, resids_type, boot_type
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: reserve(n_boot)
    integer(c_int), intent(in), optional :: excl_resids(:, :) ! Long format list of points to exclude.

    integer(c_int) :: i, j, k, i_diag, i_boot, n_rows, n_resids, n_excl_resids, i_thread

    real(c_double) :: indiv_dev_facs(n_dev - 1, n_dev - 1), dev_facs(n_dev - 1)
    real(c_double) :: sigmas(n_dev - 1), resids(n_dev - 1, n_dev - 1)
    logical(c_bool) :: triangle_mask(n_dev, n_dev), resids_mask(n_dev - 1, n_dev - 1)

    real(c_double) :: indiv_dev_facs_boot(n_dev - 1, n_dev - 1), dev_facs_boot(n_dev - 1)
    real(c_double) :: sigmas_boot(n_dev - 1), resids_boot(n_dev - 1, n_dev - 1)
    real(c_double) ::  triangle_boot(n_dev, n_dev), resampled_triangle(n_dev, n_dev)
    real(c_double) :: scale_facs(n_dev, n_dev), sigma_jack(n_dev - 1, n_dev - 1)

    real(c_double) :: latest(n_dev)

    real(c_double) :: mean, sd, shape, scale

    integer(c_int) :: status, max_stuck, stuck_counter

    status = SUCCESS
    max_stuck = 50

    i_thread = omp_get_thread_num()

    n_resids = ((n_dev - 1)**2 + (n_dev - 1))/2 - 1  ! Discard residual from upper right corner.

    triangle_mask = .true.
    resids_mask = .true.

    do j = 1, n_dev
      do i = n_dev + 2 - j, n_dev
        triangle_mask(i, j) = .false.
      end do
    end do

    do j = 1, n_dev - 1
      do i = n_dev + 1 - j, n_dev - 1
        resids_mask(i, j) = .false.
      end do
    end do

    resids_mask(1, n_dev - 1) = .false.

    if (present(excl_resids)) then
      n_excl_resids = size(excl_resids, dim=1)
      n_resids = n_resids - n_excl_resids

      do i = 1, n_excl_resids
        triangle_mask(excl_resids(i, 1), excl_resids(i, 2)) = .false.
        resids_mask(excl_resids(i, 1), excl_resids(i, 2) - 1) = .false.
      end do
    end if

    if (resids_type == PARAMETRIC) then
      call fit(triangle, dev_facs, sigmas, triangle_mask=triangle_mask)
    else if (resids_type == NORMAL_STUDENTISED) then
      call fit(triangle, dev_facs, sigmas, resids=resids, resids_type=resids_type, scale_facs=scale_facs, sigma_jack=sigma_jack)
    else if (resids_type == LOGNORMAL) then
      call fit(triangle, dev_facs, sigmas, resids=resids, resids_type=LOGNORMAL)
    else
      call fit(triangle, dev_facs, sigmas, resids=resids, resids_type=resids_type, scale_facs=scale_facs)
    end if

    stuck_counter = 0
    i_boot = 1
    main_loop: do while (i_boot <= n_boot)

      ! Parameter error.
      if (resids_type == PARAMETRIC) then
        resampled_triangle = resample(triangle, boot_type, dev_facs, sigmas, dist=dist)
        if (any(resampled_triangle < 0)) cycle main_loop
        call fit(resampled_triangle, dev_facs_boot, sigmas_boot)

      else if (resids_type == NORMAL_STUDENTISED) then
        resids_boot = sample(resids, resids_mask, rng)
        resampled_triangle = resample(triangle, boot_type, dev_facs, sigmas, &
          resids=resids_boot, resids_type=resids_type, scale_facs=scale_facs, sigma_jack=sigma_jack)
        ! if (stuck_counter > max_stuck) then
        !   call fit(triangle, dev_facs, sigmas, resids=resids, resids_type=LOGNORMAL)
        !   resampled_triangle = resample(triangle, boot_type, dev_facs, sigmas, &
        !   resids=resids, resids_type=LOGNORMAL)
        !   stuck_counter = 0
        ! end if
        if (any(resampled_triangle < 0)) then
          stuck_counter = stuck_counter + 1
          cycle main_loop
        end if
        call fit(resampled_triangle, dev_facs_boot, sigmas_boot)

      else if (resids_type == LOGNORMAL) then
        resids_boot = sample(resids, resids_mask, rng)
        resampled_triangle = resample(triangle, boot_type, dev_facs, sigmas, &
          resids=resids_boot, resids_type=LOGNORMAL)
        if (any(resampled_triangle < 0)) cycle main_loop
        call fit(resampled_triangle, dev_facs_boot, sigmas_boot)

      else
        resids_boot = sample(resids, resids_mask, rng)
        resampled_triangle = resample(triangle, boot_type, dev_facs, sigmas, &
          resids=resids_boot, resids_type=resids_type, scale_facs=scale_facs)
        if (any(resampled_triangle < 0)) cycle main_loop
        call fit(resampled_triangle, dev_facs_boot, sigmas_boot)
      end if

      ! Process error.
      triangle_boot = triangle

      if (dist == NORMAL) then
        do i_diag = 1, n_dev - 1
          do i = i_diag + 1, n_dev
            j = n_dev + i_diag + 1 - i
            mean = dev_facs_boot(j - 1) * triangle_boot(i, j - 1)
            sd = sigmas_boot(j - 1) * sqrt(triangle_boot(i, j - 1))
            triangle_boot(i, j) = rnorm_par(rng, i_thread, mean, sd)
            if (triangle_boot(i, j) <= 0) cycle main_loop
          end do
        end do

        do j = 1, n_dev
          latest(j) = triangle_boot(n_dev + 1 - j, j)
        end do

        reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)
        i_boot = i_boot + 1

      else if (dist == GAMMA) then
        do i_diag = 1, n_dev - 1
          do i = i_diag + 1, n_dev
            j = n_dev + i_diag + 1 - i
            shape = (dev_facs_boot(j - 1)**2 * triangle_boot(i, j - 1)) / sigmas_boot(j - 1) **2
            scale = sigmas_boot(j - 1) ** 2 / dev_facs_boot(j - 1)
            triangle_boot(i, j) = rgamma_par(rng, i_thread, shape, scale)
          end do
        end do

        do j = 1, n_dev
          latest(j) = triangle_boot(n_dev + 1 - j, j)
        end do

        reserve(i_boot) = sum(triangle_boot(:, n_dev)) - sum(latest)
        i_boot = i_boot + 1

      end if
    end do main_loop
  end subroutine mack_boot

  subroutine fit(triangle, dev_facs, sigmas, resids, resids_type, scale_facs, sigma_jack, triangle_mask)

    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(out) :: dev_facs(:), sigmas(:)
    real(c_double), optional, intent(out) :: resids(:, :), scale_facs(:, :)
    integer(c_int), optional, intent(in) :: resids_type
    logical(c_bool), optional, intent(in) :: triangle_mask(:, :)
    real(c_double), optional, intent(out) :: sigma_jack(:, :)

    integer(c_int) :: i, j, n_rows, n_dev, n_pts_col, n_resids
    real(c_double), allocatable :: indiv_dev_facs(:, :)
    logical(c_bool), allocatable :: col_mask(:)
    real(c_double) :: resids_mean, dev_fac_jack
    real(c_double), allocatable :: shift(:), log_normal_sigmas(:), log_normal_means(:)

    n_dev = size(triangle, 1)
    allocate(indiv_dev_facs(n_dev - 1, n_dev - 1), source=0._c_double)

    if (present(triangle_mask)) then
      do j = 1, n_dev - 1
        n_rows = n_dev - j
        do i = 1, n_rows
          indiv_dev_facs(i, j) = triangle(i, j + 1) / triangle(i, j)
        end do
        col_mask = triangle_mask(1:n_rows, j + 1)
        n_pts_col = count(col_mask)
        dev_facs(j) = sum(triangle(1:n_rows, j + 1), mask=col_mask) / sum(triangle(1:n_rows, j), mask=col_mask)
        if (n_pts_col >= 2) then
          sigmas(j) = sqrt(sum(triangle(1:n_rows, j) * (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) ** 2, mask=col_mask) / n_pts_col)
        else
          sigmas(j) = extrapolate_sigma(sigmas, j)
        end if
      end do
    else
      do j = 1, n_dev - 1
        n_rows = n_dev - j
        do i = 1, n_rows
          indiv_dev_facs(i, j) = triangle(i, j + 1) / triangle(i, j)
        end do
        dev_facs(j) = sum(triangle(1:n_rows, j + 1)) / sum(triangle(1:n_rows, j))
        if (j < n_dev - 1) then
          sigmas(j) = sqrt(sum(triangle(1:n_rows, j) * (indiv_dev_facs(1:n_rows, j) - dev_facs(j)) ** 2) / (n_rows - 1))
        else
          sigmas(j) = extrapolate_sigma(sigmas, j)
        end if
      end do
    end if

    if (present(resids) .and. present(resids_type)) then
      if (present(scale_facs)) then
        scale_facs = 0
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            scale_facs(i, j) = sqrt(1 - triangle(i, j) / sum(triangle(1:n_rows, j)))
          end do
        end do
        scale_facs(1, n_dev) = 1
      end if

      n_resids = (n_dev ** 2 - n_dev) / 2 - 1
      resids = 0
      select case (resids_type)
       case (NORMAL_STANDARDISED)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            resids(1:n_rows, j) = (indiv_dev_facs(i, j) - dev_facs(j)) * sqrt(triangle(i, j)) / (sigmas(j) * scale_facs(i, j))
          end do
        end do
       case (NORMAL_MODIFIED)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            resids(i, j) = (indiv_dev_facs(i, j) - dev_facs(j)) * sqrt(triangle(i, j)) / scale_facs(i, j)
          end do
        end do
       case (NORMAL_STUDENTISED)
        sigma_jack = 0
        allocate(col_mask(n_dev - 1))
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            col_mask = .true.
            col_mask(i) = .false.
            n_pts_col = n_rows - 1
            dev_fac_jack = sum(triangle(1:n_rows, j + 1), mask=col_mask(1:n_rows)) / &
              sum(triangle(1:n_rows, j), mask=col_mask(1:n_rows))
            if (n_pts_col >= 2) then
              sigma_jack(i, j) = sqrt(sum(triangle(1:n_rows, j) * &
                (indiv_dev_facs(1:n_rows, j) - dev_fac_jack) ** 2, mask=col_mask(1:n_rows)) / (n_pts_col - 1))
            else
              sigma_jack(i, j) = extrapolate_sigma(sigma_jack(i, :), j)
            end if
            resids(i, j) = (triangle(i, j + 1) - dev_facs(j) * triangle(i, j)) / &
              (sigma_jack(i, j) * scale_facs(i, j) * sqrt(triangle(i, j)))
          end do
        end do
       case (LOGNORMAL)
        allocate(shift(n_dev - 1))
        allocate(log_normal_means(n_dev - 1))
        allocate(log_normal_sigmas(n_dev - 1))
        do j = 1, n_dev
          n_rows = n_dev - j
          do i = 1, n_rows
            shift(i) = dev_facs(j) * sqrt(triangle(i, j)) / sigmas(j)
            log_normal_sigmas(i) = sqrt(log(1 + 1 / shift(i) ** 2))
            log_normal_means(i) = log(shift(i)) - log_normal_sigmas(i) ** 2 / 2
            resids(i, j) = (triangle(i, j + 1) - dev_facs(j) * triangle(i, j)) / (sigmas(j) * sqrt(triangle(i, j)))
            resids(i, j) = (log(resids(i, j) + shift(i)) - log_normal_means(i)) / log_normal_sigmas(i)
          end do
        end do
      end select

      resids(1, n_dev - 1) = 0

      if (resids_type /= NORMAL_STUDENTISED) then
        resids_mean = 0
        do i = 1, n_dev - 1
          do j = 1, n_dev - i
            resids_mean = resids_mean + resids(i, j)
          end do
        end do
        resids_mean = resids_mean / n_resids
        do i = 1, n_dev - 1
          do j = 1, n_dev - i
            resids(i, j) = resids(i, j) - resids_mean
          end do
        end do
      end if
    end if
    deallocate(indiv_dev_facs)
  end subroutine fit

  ! Extrapolate sigma to linearly to col.
  pure real(c_double) function extrapolate_sigma(sigmas, col)
    real(c_double), intent(in) :: sigmas(:)
    integer(c_int), intent(in) :: col

    extrapolate_sigma = sqrt(min(sigmas(col - 1) ** 2, sigmas(col - 2) ** 2, sigmas(col - 1) ** 4 / sigmas(col - 2) ** 2))
  end function extrapolate_sigma

  ! Resample a cumulative claims triangle conditionally or unconditionally,
  ! given either a set of residuals (nonparametric resample) or a distribution
  ! (parametric resample).
  function resample(triangle, boot_type, dev_facs, sigmas, dist, resids, resids_type, scale_facs, sigma_jack)
    real(c_double), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
    integer(c_int), intent(in) :: boot_type
    integer(c_int), optional, intent(in) :: dist, resids_type
    real(c_double), optional, intent(in) :: resids(:, :), scale_facs(:, :)
    real(c_double), optional, intent(in) :: sigma_jack(:, :)

    real(c_double), allocatable :: log_normal_means(:), log_normal_sigmas(:)
    real(c_double), allocatable :: shift(:)
    real(c_double), allocatable :: resample(:, :)
    integer(c_int) :: i, j, n_rows, n_dev, i_thread
    real(c_double) :: mean, sd, shape, scale

    i_thread = omp_get_thread_num()

    n_dev = size(triangle, 1)
    allocate(resample(n_dev, n_dev), source=0._c_double)

    if


    if (present(dist)) then
      if(boot_type == CONDITIONAL) then
        resample(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle(i, j - 1))
              resample(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              resample(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do
      else if (boot_type == UNCONDITIONAL) then
        resample(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (dist == NORMAL) then
              mean = dev_facs(j - 1) * resample(i, j - 1)
              sd = sigmas(j - 1) * sqrt(resample(i, j - 1))
              resample(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * resample(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              resample(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do
      end if
    else if (present(resids) .and. present(resids_type)) then
      if (resids_type == LOGNORMAL) then
        allocate(log_normal_means(n_dev - 1))
        allocate(log_normal_sigmas(n_dev - 1))
        allocate(shift(n_dev - 1))
      end if
      if (boot_type == CONDITIONAL) then
        resample(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == NORMAL_STANDARDISED) then
              resample(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(triangle(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_MODIFIED) then
              resample(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                scale_facs(i, j) * sqrt(triangle(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_STUDENTISED) then
              resample(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(triangle(i, j)) * resids(i, j)
            else if (resids_type == LOGNORMAL) then
              shift(i) = dev_facs(j) * sqrt(triangle(i, j)) / sigmas(j)
              log_normal_sigmas(i) = sqrt(log(1 + 1 / shift(i) ** 2))
              log_normal_means(i) = log(shift(i)) - log_normal_sigmas(i) ** 2 / 2
              resample(i, j + 1) = exp(resids(i, j) * log_normal_sigmas(i) + log_normal_means(i)) - shift(i)
              resample(i, j + 1) = dev_facs(j) * triangle(i, j) + sigmas(j) * sqrt(triangle(i, j)) * resample(i, j + 1)
            end if
          end do
        end do
      else if (boot_type == UNCONDITIONAL) then
        resample(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == NORMAL_STANDARDISED) then
              resample(i, j + 1) = dev_facs(j) * resample(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(resample(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_MODIFIED) then
              resample(i, j + 1) = dev_facs(j) * resample(i, j) + &
                scale_facs(i, j) * sqrt(resample(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_STUDENTISED) then
              resample(i, j + 1) = dev_facs(j) * resample(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(resample(i, j)) * resids(i, j)
            else if (resids_type == LOGNORMAL) then
              shift(i) = dev_facs(j) * sqrt(resample(i, j)) / sigmas(j)
              log_normal_sigmas(i) = sqrt(log(1 + 1 / shift(i) ** 2))
              log_normal_means(i) = log(shift(i)) - log_normal_sigmas(i) ** 2 / 2
              resample(i, j + 1) = exp(resids(i, j) * log_normal_sigmas(i) + log_normal_means(i)) - shift(i)
              resample(i, j + 1) = dev_facs(j) * resample(i, j) + sigmas(j) * sqrt(resample(i, j)) * resample(i, j + 1)
            end if
          end do
        end do
      end if
    end if
  end function resample

  ! Entry point for C wrapper, to omit excl_resids argument.
  subroutine mack_boot_f(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve) bind(c)
    integer(c_int), intent(in), value :: n_boot, n_dev, dist, resids_type, boot_type
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(inout) :: reserve(n_boot)

    integer(c_int) :: status

    rng = init_rng(1, 42)

    call mack_boot(n_dev, triangle, resids_type, boot_type, dist, n_boot, reserve, status=status)
  end subroutine mack_boot_f

end module mack
