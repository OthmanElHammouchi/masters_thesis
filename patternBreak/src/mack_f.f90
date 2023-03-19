module mack

  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use omp_lib
  use global
  use helpers
  use interface

  implicit none

  integer(c_int) :: boot_type, resids_type, dist
  real(c_double), allocatable :: scale_facs(:, :), sigma_jack(:, :)
  real(c_double), allocatable :: dev_facs_original(:), sigmas_original(:)
  real(c_double), allocatable :: dev_facs(:), sigmas(:)
  real(c_double), allocatable :: resids(:, :)

  !$omp threadprivate(boot_type, resids_type, dist, scale_facs, sigma_jack, &
  !$omp& dev_facs, sigmas, resids)

contains

  subroutine mack_sim_f(n_dev, triangle, n_config, m_config, config, type, n_boot, results) bind(c)
    integer(c_int), intent(in), value :: n_dev, n_boot, n_config, m_config, type
    real(c_double), intent(in) :: config(n_config, m_config)
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: results(n_boot * n_config, m_config + 1)

    integer(c_int) :: i, j, k, i_sim, n_threads, excl_diagidx, excl_rowidx
    integer(c_int) :: outlier_rowidx, outlier_colidx, outlier_diagidx
    real(c_double) :: factor
    integer(c_int), allocatable :: excl_resids(:, :)
    real(c_double), allocatable :: reserve(:)
    real(c_double) :: triangle_sim(n_dev, n_dev)
    type(c_ptr) :: pgbar

    integer(c_int) :: n_duplicates

    allocate(dev_facs_original(n_dev - 1), source=0._c_double)
    allocate(sigmas_original(n_dev - 1), source=0._c_double)

    call fit(triangle, dev_facs_original, sigmas_original)

    n_threads = init_omp()
    rng = init_rng(n_threads, 42)
    pgbar = pgbar_create(n_config, 1)

    results = 0

    !$omp parallel num_threads(n_threads) copyin(rng) &
    !$omp& private(reserve, i_sim, outlier_rowidx, outlier_colidx, excl_resids, triangle_sim, factor) &
    !$omp& firstprivate(m_config, n_config, n_boot, n_dev) &
    !$omp& shared(config, results, pgbar)

    allocate(scale_facs(n_dev, n_dev), source=0._c_double)
    allocate(sigma_jack(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(dev_facs(n_dev - 1), source=0._c_double)
    allocate(sigmas(n_dev - 1), source=0._c_double)
    allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(reserve(n_boot), source=0._c_double)

    i_thread = omp_get_thread_num()

    if (type == SINGLE) then
      allocate(excl_resids(1, 2))
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_config
        call pgbar_incr(pgbar)

        resids_type = int(config(i_sim, 6))
        boot_type = int(config(i_sim, 7))
        dist = int(config(i_sim, 8))
        factor = config(i_sim, 3)
        outlier_rowidx = int(config(i_sim, 1))
        outlier_colidx = int(config(i_sim, 2))
        excl_resids(1, :) = int(config(i_sim, 4:5))

        triangle_sim = single_outlier(triangle, outlier_rowidx, outlier_colidx, factor)
        call mack_boot(n_dev, triangle_sim, n_boot, reserve, excl_resids=excl_resids)

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

        triangle_sim = calendar_outlier(triangle, outlier_diagidx, factor)
        call mack_boot(n_dev, triangle_sim, n_boot, reserve, excl_resids=excl_resids)

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

        triangle_sim = origin_outlier(triangle, outlier_rowidx, factor)

        call mack_boot(n_dev, triangle_sim, n_boot, reserve, excl_resids=excl_resids)

        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_config) = transpose(spread(config(i_sim, :), 2, n_boot))
        results(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_config + 1) = reserve
      end do
      !$omp end do
    end if
    deallocate(excl_resids)
    deallocate(reserve)
    deallocate(scale_facs)
    deallocate(sigma_jack)
    deallocate(dev_facs)
    deallocate(sigmas)
    deallocate(resids)
    !$omp end parallel

    deallocate(dev_facs_original)
    deallocate(sigmas_original)

  end subroutine mack_sim_f

! Subroutine implementing bootstrap of Mack's model for claims reserving.
  subroutine mack_boot(n_dev, triangle, n_boot, reserve, excl_resids)

    integer(c_int), intent(in) :: n_boot, n_dev
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: reserve(n_boot)
    integer(c_int), intent(in), optional :: excl_resids(:, :)

    integer(c_int) :: i, j, i_diag, i_boot, n_resids, n_excl_resids
    logical(c_bool) :: triangle_mask(n_dev, n_dev), resids_mask(n_dev - 1, n_dev - 1)

    real(c_double) :: dev_facs_boot(n_dev - 1), sigmas_boot(n_dev - 1)
    real(c_double) :: resids_boot(n_dev - 1, n_dev - 1)
    real(c_double) :: triangle_boot(n_dev, n_dev), resampled_triangle(n_dev, n_dev)

    real(c_double) :: latest(n_dev)
    real(c_double) :: mean, sd, shape, scale

    integer(c_int) :: max_stuck, stuck_counter

    integer(c_int) :: status

    max_stuck = 50

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
    else
      call fit(triangle, dev_facs, sigmas, ret_resids=TRUE)
    end if

    stuck_counter = 0
    i_boot = 1
    main_loop: do while (i_boot <= n_boot)

      ! Parameter error.
      if (resids_type == PARAMETRIC) then
        call resample(triangle, PARAM_RESAMPLE, dev_facs_boot, sigmas_boot, status)
        if (status == FAILURE) cycle main_loop
      else if (boot_type == PAIRS) then
        call resample(triangle, PAIRS_RESAMPLE, dev_facs_boot, sigmas_boot, status)
        if (status == FAILURE) cycle main_loop
      else
        resids_boot = sample(resids, resids_mask)
        call resample(triangle, NON_PARAM_RESAMPLE, dev_facs_boot, sigmas_boot, status)
        if (status == FAILURE) cycle main_loop
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

  subroutine fit(triangle, dev_facs, sigmas, ret_resids, triangle_mask)

    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(out) :: dev_facs(:), sigmas(:)
    logical(c_bool), optional, intent(in) :: ret_resids
    logical(c_bool), optional, intent(in) :: triangle_mask(:, :)

    integer(c_int) :: i, j, n_rows, n_dev, n_pts_col, n_resids
    logical(c_bool) :: ret_resids_
    real(c_double), allocatable :: indiv_dev_facs(:, :)
    logical(c_bool), allocatable :: col_mask(:)
    real(c_double) :: resids_mean, dev_fac_jack
    real(c_double), allocatable :: shift(:), log_normal_sigmas(:), log_normal_means(:)

    ret_resids_ = .false.

    if (present(ret_resids)) then
      ret_resids_ = ret_resids
    end if

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

    if (ret_resids_) then
      if (resids_type /= LOGNORMAL) then
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
  subroutine resample(triangle, type, dev_facs_boot, sigmas_boot, status)
    real(c_double), intent(in) :: triangle(:, :)
    integer(c_int), intent(in) :: type
    real(c_double), intent(out) :: dev_facs_boot(:), sigmas_boot(:)
    integer(c_int), intent(out) :: status

    real(c_double), allocatable :: log_normal_means(:), log_normal_sigmas(:)
    real(c_double), allocatable :: shift(:)
    real(c_double), allocatable :: triangle_new(:, :)
    integer(c_int) :: i, j, n_rows, n_dev
    real(c_double) :: mean, sd, shape, scale

    logical(c_bool), allocatable :: pairs_mask(:)
    real(c_double), allocatable :: pairs(:, :)
    integer(c_int), allocatable :: pairs_indices(:), resampled_pairs_indices(:)

    n_dev = size(triangle, 1)
    allocate(triangle_new(n_dev, n_dev), source=0._c_double)

    select case (type)
     case (PARAM_RESAMPLE)
      if(boot_type == CONDITIONAL) then
        triangle_new(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle(i, j - 1))
              triangle_new(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              triangle_new(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do
      else if (boot_type == UNCONDITIONAL) then
        triangle_new(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle_new(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle_new(i, j - 1))
              triangle_new(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle_new(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              triangle_new(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do
      end if
      if (any(triangle_new < 0)) then
        status = FAILURE
        return
      end if
      call fit(triangle, dev_facs_boot, sigmas_boot)

     case (NON_PARAM_RESAMPLE)
      if (resids_type == LOGNORMAL) then
        allocate(log_normal_means(n_dev - 1))
        allocate(log_normal_sigmas(n_dev - 1))
        allocate(shift(n_dev - 1))
      end if
      if (boot_type == CONDITIONAL) then
        triangle_new(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == NORMAL_STANDARDISED) then
              triangle_new(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(triangle(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_MODIFIED) then
              triangle_new(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                scale_facs(i, j) * sqrt(triangle(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_STUDENTISED) then
              triangle_new(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(triangle(i, j)) * resids(i, j)
            else if (resids_type == LOGNORMAL) then
              shift(i) = dev_facs(j) * sqrt(triangle(i, j)) / sigmas(j)
              log_normal_sigmas(i) = sqrt(log(1 + 1 / shift(i) ** 2))
              log_normal_means(i) = log(shift(i)) - log_normal_sigmas(i) ** 2 / 2
              triangle_new(i, j + 1) = exp(resids(i, j) * log_normal_sigmas(i) + log_normal_means(i)) - shift(i)
              triangle_new(i, j + 1) = dev_facs(j) * triangle(i, j) + sigmas(j) * sqrt(triangle(i, j)) * triangle_new(i, j + 1)
            end if
          end do
        end do
      else if (boot_type == UNCONDITIONAL) then
        triangle_new(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == NORMAL_STANDARDISED) then
              triangle_new(i, j + 1) = dev_facs(j) * triangle_new(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(triangle_new(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_MODIFIED) then
              triangle_new(i, j + 1) = dev_facs(j) * triangle_new(i, j) + &
                scale_facs(i, j) * sqrt(triangle_new(i, j)) * resids(i, j)
            else if (resids_type == NORMAL_STUDENTISED) then
              triangle_new(i, j + 1) = dev_facs(j) * triangle_new(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(triangle_new(i, j)) * resids(i, j)
            else if (resids_type == LOGNORMAL) then
              shift(i) = dev_facs(j) * sqrt(triangle_new(i, j)) / sigmas(j)
              log_normal_sigmas(i) = sqrt(log(1 + 1 / shift(i) ** 2))
              log_normal_means(i) = log(shift(i)) - log_normal_sigmas(i) ** 2 / 2
              triangle_new(i, j + 1) = exp(resids(i, j) * log_normal_sigmas(i) + log_normal_means(i)) - shift(i)
              triangle_new(i, j + 1) = dev_facs(j) * triangle_new(i, j) + &
                sigmas(j) * sqrt(triangle_new(i, j)) * triangle_new(i, j + 1)
            end if
          end do
        end do
      end if
      if (any(triangle_new < 0)) then
        status = FAILURE
        return
      end if
      call fit(triangle, dev_facs_boot, sigmas_boot)

     case (PAIRS_RESAMPLE)
      allocate(pairs_mask(n_dev - 1), source=.true._c_bool)
      allocate(pairs(n_dev - 1, 2), source=0._c_double)
      allocate(pairs_indices(n_dev - 1))
      allocate(resampled_pairs_indices(n_dev - 1))

      do i = 1, n_dev - 1
        pairs_indices(i) = i
      end do
      do j = 1, n_dev - 1
        n_rows = n_dev - j
        do i = 1, n_rows
          resampled_pairs_indices(i) = pairs_indices(1 + int(n_rows * runif_par(rng, i_thread)))
        end do
        pairs(1:n_rows, :) = triangle(resampled_pairs_indices(1:n_rows), j:(j + 1))
        dev_facs_boot(j) = sum(pairs(1:n_rows, 2)) / sum(pairs(1:n_rows, 1))
        if (j < n_dev - 1) then
          sigmas_boot(j) = sqrt(sum(triangle(1:n_rows, j) * (triangle(1:n_rows, j + 1) / triangle(1:n_rows, j) - &
            dev_facs_boot(j)) ** 2) / (n_rows - 1))
        else
          sigmas_boot(j) = extrapolate_sigma(sigmas_boot, j)
        end if
      end do
      if (any(dev_facs_boot < 1)) then
        status = FAILURE
      end if 
    end select
  end subroutine resample

  ! Entry point for C wrapper, to omit excl_resids argument.
  subroutine mack_boot_f(n_dev, triangle, resids_type_, boot_type_, dist_, n_boot, reserve) bind(c)
    integer(c_int), intent(in), value :: n_boot, n_dev, dist_, resids_type_, boot_type_
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(inout) :: reserve(n_boot)

    boot_type = boot_type_
    resids_type = resids_type_
    dist = dist_
    rng = init_rng(1, 42)

    allocate(scale_facs(n_dev, n_dev), source=0._c_double)
    allocate(sigma_jack(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(dev_facs(n_dev - 1), source=0._c_double)
    allocate(sigmas(n_dev - 1), source=0._c_double)
    allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)

    call mack_boot(n_dev, triangle, n_boot, reserve)
  end subroutine mack_boot_f

  function single_outlier(triangle, outlier_rowidx, outlier_colidx, factor)
    integer(c_int), intent(in):: outlier_rowidx, outlier_colidx
    real(c_double), intent(in) :: triangle(:, :)

    real(c_double) :: factor
    real(c_double), allocatable:: single_outlier(:, :)
    real(c_double) :: shape, scale, mean, sd
    integer(c_int) :: n_dev, i, j

    n_dev = size(triangle, 1)
    allocate(single_outlier(n_dev, n_dev), source=0._c_double)
    single_outlier(:, 1) = triangle(:, 1)

    if (dist == NORMAL) then
      do j = 2, n_dev
        do i = 1, n_dev + 1 - j
          if (i == outlier_rowidx) cycle
          mean = dev_facs_original(j - 1) * single_outlier(i, j - 1)
          sd = sigmas_original(j - 1) * sqrt(single_outlier(i, j - 1))
          single_outlier(i, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end do

      if (outlier_colidx > 2) then
        do j = 2, outlier_colidx - 1
          mean = dev_facs_original(j - 1) * single_outlier(outlier_rowidx, j - 1)
          sd = sigmas_original(j - 1) * sqrt(single_outlier(outlier_rowidx, j - 1))
          single_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end if

      mean = factor * dev_facs_original(outlier_colidx - 1) * single_outlier(outlier_rowidx, outlier_colidx - 1)
      sd = sigmas_original(outlier_colidx - 1) * sqrt(single_outlier(outlier_rowidx, outlier_colidx - 1))
      single_outlier(outlier_rowidx, outlier_colidx) = rnorm_par(rng, i_thread, mean, sd)

      if (outlier_colidx < n_dev) then
        do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
          mean = dev_facs_original(j - 1) * single_outlier(outlier_rowidx, j - 1)
          sd = sigmas_original(j - 1) * sqrt(single_outlier(outlier_rowidx, j - 1))
          single_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end if

    else if (dist == GAMMA) then
      do j = 2, n_dev
        do i = 1, n_dev + 1 - j
          if (i == outlier_rowidx) cycle
          shape = dev_facs_original(j - 1)**2 * single_outlier(i, j - 1) / sigmas_original(j - 1)**2
          scale = sigmas_original(j - 1)**2 / dev_facs_original(j - 1)
          single_outlier(i, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end do

      if (outlier_colidx > 2) then
        do j = 2, outlier_colidx - 1
          shape = dev_facs_original(j - 1)**2 * single_outlier(outlier_rowidx, j - 1) / sigmas_original(j - 1)**2
          scale = sigmas_original(j - 1)**2 / dev_facs_original(j - 1)
          single_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end if

      shape = dev_facs_original(outlier_colidx - 1)**2 * single_outlier(outlier_rowidx, outlier_colidx - 1) / &
        sigmas_original(outlier_colidx - 1)**2
      scale = sigmas_original(outlier_colidx - 1)**2 / dev_facs_original(outlier_colidx - 1)
      single_outlier(outlier_rowidx, outlier_colidx) = rgamma_par(rng, i_thread, shape, scale)

      if (outlier_colidx < n_dev) then
        do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
          shape = dev_facs_original(j - 1)**2 * single_outlier(outlier_rowidx, j - 1) / sigmas_original(j - 1)**2
          scale = sigmas_original(j - 1)**2 / dev_facs_original(j - 1)
          single_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end if
    end if
  end function single_outlier

  function calendar_outlier(triangle, outlier_diagidx, factor)
    integer(c_int), intent(in) :: outlier_diagidx
    real(c_double), intent(in):: factor
    real(c_double), intent(in) :: triangle(:, :)

    integer(c_int) :: i, j, n_dev, n_cols
    real(c_double), allocatable :: calendar_outlier(:, :)
    real(c_double) :: shape, scale, mean, sd

    n_dev = size(triangle, dim=1)

    allocate(calendar_outlier(n_dev, n_dev), source=0._c_double)
    calendar_outlier(:, 1) = triangle(:, 1)

    do i = 1, n_dev
      n_cols = n_dev + 2 - outlier_diagidx - i
      if (n_cols <= 1) then
        do j = 2, n_dev + 1 - i
          calendar_outlier(i, j) = triangle(i, j)
        end do
      else
        do j = 2, n_cols - 1
          calendar_outlier(i, j) = triangle(i, j)
        end do

        if (dist == NORMAL) then
          mean = factor * dev_facs_original(n_cols - 1) * calendar_outlier(i, n_cols - 1)
          sd = sigmas_original(n_cols - 1) * sqrt(calendar_outlier(i, n_cols - 1))
          calendar_outlier(i, n_cols) = rnorm_par(rng, i_thread, mean, sd)
          do j = n_cols + 1, n_dev + 1 - i
            mean = dev_facs_original(j - 1) * calendar_outlier(i, j - 1)
            sd = sigmas_original(j - 1) * sqrt(calendar_outlier(i, j - 1))
            calendar_outlier(i, j) = rnorm_par(rng, i_thread, mean, sd)
          end do

        else if (dist == GAMMA) then
          shape = factor * dev_facs_original(n_cols - 1)**2 * calendar_outlier(i, n_cols - 1) / sigmas_original(n_cols - 1)**2
          scale = sigmas_original(n_cols - 1)**2 / factor * dev_facs_original(n_cols - 1)
          calendar_outlier(i, n_cols) = rgamma_par(rng, i_thread, shape, scale)
          do j = n_cols + 1, n_dev + 1 - i
            shape = dev_facs_original(j - 1)**2 * calendar_outlier(i, j - 1) / sigmas_original(j - 1)**2
            scale = sigmas_original(j - 1)**2 / dev_facs_original(j - 1)
            calendar_outlier(i, j) = rgamma_par(rng, i_thread, shape, scale)
          end do
        end if
      end if
    end do
  end function calendar_outlier

  function origin_outlier(triangle, outlier_rowidx, factor)
    integer(c_int), intent(in):: outlier_rowidx
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in) :: factor

    real(c_double) :: shape, scale, mean, sd
    real(c_double), allocatable:: origin_outlier(:, :)
    integer(c_int) :: n_dev, j

    n_dev = size(triangle, dim=1)
    origin_outlier = triangle

    do j = 2, n_dev + 1 - outlier_rowidx
      if (dist == NORMAL) then
        mean = factor * dev_facs_original(j - 1) * origin_outlier(outlier_rowidx, j - 1)
        sd = sigmas_original(j - 1) * sqrt(origin_outlier(outlier_rowidx, j - 1))
        origin_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)

      else if (dist == GAMMA) then
        shape = factor * dev_facs_original(j - 1)**2 * origin_outlier(outlier_rowidx, j - 1) / sigmas_original(j - 1)**2
        scale = sigmas_original(j - 1)**2 / factor * dev_facs_original(j - 1)
        origin_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
      end if
    end do
  end function origin_outlier

end module mack
