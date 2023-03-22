module mack

  use, intrinsic :: iso_c_binding
  use, intrinsic :: omp_lib
  use global_mod
  use helpers_mod
  use interface_mod
  use triangle_mod
  use sim_mod

  implicit none

contains

  subroutine simulate(sim_type, n_boot, n_dev, triangle_in, factors, boot_types, &
    proc_dists, conds, resids_types, n_res, m_res, res, show_progress) &
    bind(c, name="mack_sim_")
    integer(c_int), intent(in), value :: n_dev, n_boot, sim_type
    integer(c_int), intent(in), value :: n_res, m_res
    logical(c_bool), intent(in), value :: show_progress
    real(c_double), intent(in) :: triangle_in(n_dev, n_dev)
    real(c_double), intent(out) :: res(n_res, m_res)
    integer(c_int) :: boot_types(:), resids_types(:)
    integer(c_int) :: proc_dists(:)
    logical(c_bool) :: conds(:)
    real(c_double) :: factors(:)

    integer(c_int) :: i, j, k, l, i_sim, n_rows, n_threads
    integer(c_int) :: n_outliers, n_excl, n_sim
    integer(c_int) :: excl_diagidx, excl_rowidx, excl_colidx
    integer(c_int) :: outlier_rowidx, outlier_colidx, outlier_diagidx
    real(c_double) :: factor
    real(c_double) :: dev_facs_original(n_dev - 1), sigmas_original(n_dev - 1)
    real(c_double) :: triangle_sim(n_dev, n_dev)
    integer(c_int), allocatable :: excl(:, :), outliers(:, :)
    real(c_double), allocatable :: reserve(:)
    type(c_ptr) :: pgbar

    type(triangle) :: triangle_
    type(single_sim) :: single_sim_
    type(calendar_sim) :: calendar_sim_
    type(origin_sim) :: origin_sim_
    class(mack_sim) :: sim_

    triangle_ = triangle(triangle_in)

    n_threads = init_omp()
    rng = init_rng(n_threads, 42)

    pgbar = pgbar_create(n_sim, 1)

    !$omp parallel num_threads(n_threads) &
    !$omp& copyin(rng) &
    !$omp& private(reserve, i_sim, triangle_sim, triangle_, sim_) &
    !$omp& firstprivate(sim_, n_boot, n_dev, m_res, n_sim) &
    !$omp& shared(res, pgbar)

    allocate(reserve(n_boot), source=0._c_double)

    i_thread = omp_get_thread_num()

    n_res = n_boot * sim_%n_sim
    m_res = sim_%m_sim + 1

    reserve = 0
    select case (sim_type)
     case (SINGLE)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        single_sim_ = single_sim(n_dev, factors, boot_types, proc_dists, conds, resids_types)

        triangle_%mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            triangle_%mask(i, j) = .false.
          end do
        end do
        triangle_%mask(excl_rowidx, excl_colidx) = .false.
        if (boot_type == PARAMETRIC) triangle_%mask(1, n_dev) = .true.

        triangle_sim = single_outlier(triangle_, single_sim_)
        call boot(triangle_sim, single_sim_, n_boot, reserve)
        call single_sim_%update()

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_res) = spread(single_sim_%table(i_sim, 1:(m_res - 1)), 1, n_boot)
        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), m_res) = reserve
      end do
      !$omp end do

     case (CALENDAR)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        calendar_sim_ = calendar_sim(n_dev, factors, boot_types, proc_dists, conds, resids_types)

        triangle_%mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            triangle_%mask(i, j) = .false.
          end do
        end do
        do j = 2, n_dev
          i = n_dev + 2 - excl_diagidx - j
          if (i <= 0) cycle
          triangle_%mask(i, j) = .false.
        end do
        if (boot_type == PARAMETRIC) triangle_%mask(1, n_dev) = .true.

        triangle_sim = calendar_outlier(triangle_, calendar_sim_)
        call boot(triangle_sim, calendar_sim_, n_boot, reserve)
        call calendar_sim_%update()

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_res) = spread(calendar_sim_%table(i_sim, 1:(m_res - 1)), 1, n_boot)
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_res) = reserve
      end do
      !$omp end do

     case (ORIGIN)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        origin_sim_ = origin_sim(n_dev, factors, boot_types, proc_dists, conds, resids_types)

        triangle_%mask = .true.
        do j = 2, n_dev + 1 - excl_rowidx
          triangle_%mask(excl_rowidx, j) = .false.
        end do
        if (boot_type == PARAMETRIC) triangle_%mask(1, n_dev) = .true.

        triangle_sim = origin_outlier(triangle_, sim_)
        call boot(triangle_sim, sim_, n_boot, reserve)
        call origin_sim_%update()

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_res) = spread(origin_sim_%table(i_sim, 1:(m_res - 1)), 1, n_boot)
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_res) = reserve
      end do
      !$omp end do
    end select

    ! Free shared containers
    deallocate(reserve)
    !$omp end parallel
  end subroutine simulate

  subroutine boot(triangle_in, sim_in, n_boot, reserve)
    integer(c_int), intent(in) :: n_boot
    type(triangle), intent(in) :: triangle_in
    type(sim), intent(in) :: sim_in
    real(c_double), intent(out) :: reserve(n_boot)

    integer(c_int) :: i, j, i_diag, i_boot, n_dev, n_resids, n_excl_resids
    real(c_double), allocatable :: dev_facs_boot(:), sigmas_boot(:)
    logical(c_bool), allocatable :: resids_mask(:, :)
    real(c_double), allocatable :: triangle_boot(:, :)
    real(c_double), allocatable :: latest(:)
    real(c_double) :: mean, sd, shape, scale
    integer(c_int) :: status

    integer(c_int) :: boot_type, proc_dist, resids_type
    logical(c_bool) :: cond

    select type (sim_in)
     type is (single_sim)
      boot_type = sim_in%boot_type
      proc_dist = sim_in%proc_dist
      cond = sim_in%cond
      resids_type = sim_in%resids_type
     type is (calendar_sim)
      boot_type = sim_in%boot_type
      proc_dist = sim_in%proc_dist
      cond = sim_in%cond
      resids_type = sim_in%resids_type
     type is (origin_sim)
      boot_type = sim_in%boot_type
      proc_dist = sim_in%proc_dist
      cond = sim_in%cond
      resids_type = sim_in%resids_type
    end select

    n_dev = triangle_in%n_dev

    allocate(dev_facs_boot(n_dev - 1), source=0._c_double)
    allocate(sigmas_boot(n_dev - 1), source=0._c_double)
    allocate(resids_mask(n_dev - 1, n_dev - 1), source=.true._c_bool)
    allocate(latest(n_dev - 1), source=0._c_double)

    if (boot_type == PARAMETRIC) then
      call triangle_in%fit(use_mask=.true., compute_resids=.false.)
    else if (boot_type == RESID) then
      call triangle_in%fit(use_mask=.false., compute_resids=.true.)
    end if

    i_boot = 1
    main_loop: do while (i_boot <= n_boot)
      ! Simulate new parameters
      if (boot_type == RESID) then
        resids_mask = triangle_in%mask(1:(n_dev - 1), 2:n_dev)
        resampled_resids = sample(triangle_in%resids, resids_mask)
      end if

      call resample(triangle_in, sim_in, dev_facs_boot, sigmas_boot, status)
      if (status == FAILURE) cycle main_loop

      ! Simulate process error
      triangle_boot = triangle_

      if (proc_dist == NORMAL) then
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

      else if (proc_dist == GAMMA) then
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
  end subroutine boot

  ! Return simulated development factors and dispersion parameters.
  subroutine resample(triangle_in, sim_in, dev_facs_boot, sigmas_boot, status)
    type(triangle), intent(in) :: triangle_in
    type(sim), intent(in) :: sim_in
    real(c_double), intent(out) :: dev_facs_boot(:), sigmas_boot(:)
    integer(c_int), intent(out) :: status

    real(c_double) :: resampled_triangle(:, :)
    type(triangle) :: resampled_triangle_
    integer(c_int), allocatable :: i, j, n_rows, n_dev
    real(c_double), allocatable :: mean, sd, shape, scale
    real(c_double), allocatable :: triangle_(:, :), dev_facs(:), sigmas(:)
    real(c_double), allocatable :: scale_facs(:, :), sigmas_jack(:, :)
    real(c_double), allocatable :: resids(:, :), log_normal_shift(:)
    real(c_double) :: log_normal_means(:), log_normal_sigmas(:)
    real(c_double), allocatable :: resampled_pairs(:, :)
    integer(c_int), allocatable :: pair_indices(:), resampled_pair_indices(:)

    integer(c_int) :: boot_type, proc_dist, resids_type
    logical(c_bool) :: cond

    select type (sim_in)
     type is (single_sim)
      boot_type = sim_in%boot_type
      proc_dist = sim_in%proc_dist
      cond = sim_in%cond
      resids_type = sim_in%resids_type

      triangle_ = triangle_in%data
      dev_facs = triangle_in%dev_facs
      sigmas = triangle_in%sigmas
      sigmas_jack = triangle_in%sigmas_jack
      scale_facs = triangle_in%scale_facs
      log_normal_shift = triangle_in%log_normal_shift
      log_normal_means = triangle_in%log_normal_means
      log_normal_sigmas = triangle_in%log_normal_sigmas
      resis = triangle_in%resids
     type is (calendar_sim)
      boot_type = sim_in%boot_type
      proc_dist = sim_in%proc_dist
      cond = sim_in%cond
      resids_type = sim_in%resids_type

      triangle_ = triangle_in%data
      dev_facs = triangle_in%dev_facs
      sigmas = triangle_in%sigmas
      sigmas_jack = triangle_in%sigmas_jack
      scale_facs = triangle_in%scale_facs
      log_normal_shift = triangle_in%log_normal_shift
      log_normal_means = triangle_in%log_normal_means
      log_normal_sigmas = triangle_in%log_normal_sigmas
      resis = triangle_in%resids
     type is (origin_sim)
      boot_type = sim_in%boot_type
      proc_dist = sim_in%proc_dist
      cond = sim_in%cond
      resids_type = sim_in%resids_type

      triangle_ = triangle_in%data
      dev_facs = triangle_in%dev_facs
      sigmas = triangle_in%sigmas
      sigmas_jack = triangle_in%sigmas_jack
      scale_facs = triangle_in%scale_facs
      log_normal_shift = triangle_in%log_normal_shift
      log_normal_means = triangle_in%log_normal_means
      log_normal_sigmas = triangle_in%log_normal_sigmas
      resis = triangle_in%resids
    end select

    n_dev = size(triangle_, 1)
    allocate(resampled_triangle(n_dev, n_dev), source=0._c_double)

    status = SUCCESS
    select case (boot_type)
     case (PARAMETRIC)
      if(cond) then
        resampled_triangle(:, 1) = triangle_(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (proc_dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle_(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle_(i, j - 1))
              resampled_triangle(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (proc_dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle_(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              resampled_triangle(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do

      else
        resampled_triangle(:, 1) = triangle_(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (proc_dist == NORMAL) then
              mean = dev_facs(j - 1) * resampled_triangle(i, j - 1)
              sd = sigmas(j - 1) * sqrt(resampled_triangle(i, j - 1))
              resampled_triangle(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (proc_dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * resampled_triangle(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              resampled_triangle(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do
      end if

      do j = 1, n_dev
        n_rows = n_dev + 1 - j
        do i = 1, n_rows
          if (resampled_triangle(i, j) <= 0) then
            status = FAILURE
            return
          end if
        end do
      end do

      resampled_triangle_ = triangle(n_dev, resampled_triangle)
      dev_facs_boot = resampled_triangle_%dev_facs
      sigmas_boot = resampled_triangle_%sigmas

     case (RESID)
      if (cond) then
        resampled_triangle(:, 1) = triangle_(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == STANDARDISED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle_(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(triangle_(i, j)) * resampled_resids(i, j)
            else if (resids_type == MODIFIED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle_(i, j) + &
                scale_facs(i, j) * sqrt(triangle_(i, j)) * resampled_resids(i, j)
            else if (resids_type == STUDENTISED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle_(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(triangle_(i, j)) * resampled_resids(i, j)
            else if (resids_type == LOGNORMAL) then
              resampled_triangle(i, j + 1) = exp(resampled_resids(i, j) * log_normal_sigmas(i) + &
                log_normal_means(i)) - log_normal_shift(i)
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle_(i, j) + &
                sigmas(j) * sqrt(triangle_(i, j)) * resampled_triangle(i, j + 1)
            end if
          end do
        end do

      else
        resampled_triangle(:, 1) = triangle_(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == STANDARDISED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * resampled_triangle(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(resampled_triangle(i, j)) * resampled_resids(i, j)
            else if (resids_type == MODIFIED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * resampled_triangle(i, j) + &
                scale_facs(i, j) * sqrt(resampled_triangle(i, j)) * resampled_resids(i, j)
            else if (resids_type == STUDENTISED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * resampled_triangle(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(resampled_triangle(i, j)) * resampled_resids(i, j)
            else if (resids_type == LOGNORMAL) then
              resampled_triangle(i, j + 1) = exp(resampled_resids(i, j) * log_normal_sigmas(i) + &
                log_normal_means(i)) - log_normal_shift(i)
              resampled_triangle(i, j + 1) = dev_facs(j) * resampled_triangle(i, j) + &
                sigmas(j) * sqrt(resampled_triangle(i, j)) * resampled_triangle(i, j + 1)
            end if
          end do
        end do
      end if

      do j = 1, n_dev
        n_rows = n_dev + 1 - j
        do i = 1, n_rows
          if (resampled_triangle(i, j) <= 0) then
            status = FAILURE
            return
          end if
        end do
      end do

      resampled_triangle_ = triangle(n_dev, resampled_triangle)
      dev_facs_boot = resampled_triangle_%dev_facs
      sigmas_boot = resampled_triangle_%sigmas

     case (PAIRS)
      allocate(resampled_pairs(n_dev - 1, 2), source=0._c_double)
      allocate(pair_indices(n_dev - 1))
      allocate(resampled_pair_indices(n_dev - 1))

      do i = 1, n_dev - 1
        pair_indices(i) = i
      end do

      do j = 1, n_dev - 1
        n_rows = n_dev - j
        do i = 1, n_rows
          resampled_pair_indices(i) = pair_indices(1 + int(n_rows * runif_par(rng, i_thread)))
        end do
        resampled_pairs(1:n_rows, :) = triangle_(resampled_pair_indices(1:n_rows), j:(j + 1))

        dev_facs_boot(j) = sum(resampled_pairs(1:n_rows, 2)) / sum(resampled_pairs(1:n_rows, 1))
        if (j < n_dev - 1) then
          sigmas_boot(j) = sqrt(sum(triangle_(1:n_rows, j) * (triangle_(1:n_rows, j + 1) / triangle_(1:n_rows, j) - &
            dev_facs_boot(j)) ** 2) / (n_rows - 1))
        else
          sigmas_boot(j) = extrapolate_sigma(sigmas_boot, j)
        end if
      end do
    end select
  end subroutine resample

  ! subroutine mack_boot_(n_dev, triangle, boot_type_, proc_dist_, cond_, &
  !   resids_type_, n_boot, reserve) bind(c)
  !   integer(c_int), intent(in), value :: n_boot, n_dev, proc_dist_, resids_type_, boot_type_
  !   logical(c_bool), intent(in), value :: cond_
  !   real(c_double), intent(in) :: triangle(n_dev, n_dev)
  !   real(c_double), intent(out) :: reserve(n_boot)

  !   allocate(dev_facs(n_dev - 1), source=0._c_double)
  !   allocate(sigmas(n_dev - 1), source=0._c_double)
  !   allocate(scale_facs(n_dev, n_dev), source=0._c_double)
  !   allocate(sigma_jack(n_dev - 1, n_dev - 1), source=0._c_double)
  !   allocate(log_normal_shift(n_dev - 1), source=0._c_double)
  !   allocate(log_normal_means(n_dev - 1), source=0._c_double)
  !   allocate(log_normal_sigmas(n_dev - 1), source=0._c_double)
  !   allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)

  !   boot_type = boot_type_
  !   proc_dist = proc_dist_
  !   cond = cond_
  !   resids_type = resids_type_

  !   allocate(mask(n_dev, n_dev))
  !   mask = .true.

  !   rng = init_rng(1, 42)
  !   i_thread = 0

  !   call boot(n_dev, triangle, n_boot, reserve, mask)
  ! end subroutine mack_boot_

  pure type(triangle) function single_outlier(triangle_in, sim_in)
    type(triangle), intent(in) :: triangle_in
    type(sim), intent(in) :: sim_in

    real(c_double), allocatable:: single_outlier_(:, :)
    real(c_double) :: shape, scale, mean, sd
    integer(c_int) :: n_dev, i, j

    real(c_double) :: factor
    real(c_double), allocatable :: dev_facs(:), sigmas(:)
    real(c_double), allocatable :: triangle_(:, :)
    integer(c_int) :: proc_dist
    integer(c_int) :: outlier_rowidx, outlier_colidx

    outlier_rowidx = sim_in%outlier_rowidx
    outlier_colidx = sim_in%outlier_colidx
    proc_dist = sim_in%proc_dist
    factor = sim_in%factor
    triangle_ = triangle_in%data
    dev_facs = triangle_in%dev_facs
    sigmas = triangle_in%sigmas

    if (factor == 1) then
      single_outlier_ = triangle_
      return
    end if

    n_dev = size(triangle_, 1)
    allocate(single_outlier_(n_dev, n_dev), source=0._c_double)
    single_outlier_(:, 1) = triangle_(:, 1)

    if (proc_dist == NORMAL) then
      do j = 2, n_dev
        do i = 1, n_dev + 1 - j
          if (i == outlier_rowidx) cycle
          mean = dev_facs(j - 1) * single_outlier_(i, j - 1)
          sd = sigmas(j - 1) * sqrt(single_outlier_(i, j - 1))
          single_outlier_(i, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end do

      if (outlier_colidx > 2) then
        do j = 2, outlier_colidx - 1
          mean = dev_facs(j - 1) * single_outlier_(outlier_rowidx, j - 1)
          sd = sigmas(j - 1) * sqrt(single_outlier_(outlier_rowidx, j - 1))
          single_outlier_(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end if

      mean = factor * dev_facs(outlier_colidx - 1) * single_outlier_(outlier_rowidx, outlier_colidx - 1)
      sd = sigmas(outlier_colidx - 1) * sqrt(single_outlier_(outlier_rowidx, outlier_colidx - 1))
      single_outlier_(outlier_rowidx, outlier_colidx) = rnorm_par(rng, i_thread, mean, sd)

      if (outlier_colidx < n_dev) then
        do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
          mean = dev_facs(j - 1) * single_outlier_(outlier_rowidx, j - 1)
          sd = sigmas(j - 1) * sqrt(single_outlier_(outlier_rowidx, j - 1))
          single_outlier_(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end if

    else if (proc_dist == GAMMA) then
      do j = 2, n_dev
        do i = 1, n_dev + 1 - j
          if (i == outlier_rowidx) cycle
          shape = dev_facs(j - 1)**2 * single_outlier_(i, j - 1) / sigmas(j - 1)**2
          scale = sigmas(j - 1)**2 / dev_facs(j - 1)
          single_outlier_(i, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end do

      if (outlier_colidx > 2) then
        do j = 2, outlier_colidx - 1
          shape = dev_facs(j - 1)**2 * single_outlier_(outlier_rowidx, j - 1) / sigmas(j - 1)**2
          scale = sigmas(j - 1)**2 / dev_facs(j - 1)
          single_outlier_(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end if

      shape = dev_facs(outlier_colidx - 1)**2 * single_outlier_(outlier_rowidx, outlier_colidx - 1) / &
        sigmas(outlier_colidx - 1)**2
      scale = sigmas(outlier_colidx - 1)**2 / dev_facs(outlier_colidx - 1)
      single_outlier_(outlier_rowidx, outlier_colidx) = rgamma_par(rng, i_thread, shape, scale)

      if (outlier_colidx < n_dev) then
        do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
          shape = dev_facs(j - 1)**2 * single_outlier_(outlier_rowidx, j - 1) / sigmas(j - 1)**2
          scale = sigmas(j - 1)**2 / dev_facs(j - 1)
          single_outlier_(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end if
    end if

    single_outlier = triangle(single_outlier_)
  end function single_outlier

  pure type(triangle) function calendar_outlier(triangle_in, sim_in)
    type(triangle), intent(in) :: triangle_in
    type(sim), intent(in) :: sim_in

    integer(c_int) :: i, j, n_dev, n_cols
    real(c_double), allocatable :: calendar_outlier_(:, :)
    real(c_double) :: shape, scale, mean, sd

    real(c_double) :: factor
    real(c_double), allocatable :: dev_facs(:), sigmas(:)
    real(c_double), allocatable :: triangle_(:, :)
    integer(c_int) :: proc_dist
    integer(c_int) :: outlier_diagidx

    outlier_diagidx = sim_in%outlier_diagidx
    proc_dist = sim_in%proc_dist
    factor = sim_in%factor
    triangle_ = triangle_in%data
    dev_facs = triangle_in%dev_facs
    sigmas = triangle_in%sigmas

    if (factor == 1) then
      calendar_outlier_ = triangle
      return
    end if

    n_dev = size(triangle_, dim=1)

    allocate(calendar_outlier_(n_dev, n_dev), source=0._c_double)
    calendar_outlier_(:, 1) = triangle_(:, 1)

    do i = 1, n_dev
      n_cols = n_dev + 2 - outlier_diagidx - i
      if (n_cols <= 1) then
        do j = 2, n_dev + 1 - i
          calendar_outlier_(i, j) = triangle_(i, j)
        end do
      else
        do j = 2, n_cols - 1
          calendar_outlier_(i, j) = triangle_(i, j)
        end do

        if (proc_dist == NORMAL) then
          mean = factor * dev_facs(n_cols - 1) * calendar_outlier_(i, n_cols - 1)
          sd = sigmas(n_cols - 1) * sqrt(calendar_outlier_(i, n_cols - 1))
          calendar_outlier_(i, n_cols) = rnorm_par(rng, i_thread, mean, sd)
          do j = n_cols + 1, n_dev + 1 - i
            mean = dev_facs(j - 1) * calendar_outlier_(i, j - 1)
            sd = sigmas(j - 1) * sqrt(calendar_outlier_(i, j - 1))
            calendar_outlier_(i, j) = rnorm_par(rng, i_thread, mean, sd)
          end do

        else if (proc_dist == GAMMA) then
          shape = factor * dev_facs(n_cols - 1)**2 * calendar_outlier_(i, n_cols - 1) / sigmas(n_cols - 1)**2
          scale = sigmas(n_cols - 1)**2 / factor * dev_facs(n_cols - 1)
          calendar_outlier_(i, n_cols) = rgamma_par(rng, i_thread, shape, scale)
          do j = n_cols + 1, n_dev + 1 - i
            shape = dev_facs(j - 1)**2 * calendar_outlier_(i, j - 1) / sigmas(j - 1)**2
            scale = sigmas(j - 1)**2 / dev_facs(j - 1)
            calendar_outlier_(i, j) = rgamma_par(rng, i_thread, shape, scale)
          end do
        end if
      end if
    end do

    calendar_outlier = triangle(calendar_outlier_)
  end function calendar_outlier

  pure type(triangle) function origin_outlier(triangle_in, sim_in)
    type(triangle), intent(in) :: triangle_in
    type(sim), intent(in) :: sim_in

    real(c_double) :: shape, scale, mean, sd
    real(c_double), allocatable:: origin_outlier_(:, :)
    integer(c_int) :: n_dev, j

    real(c_double) :: factor
    real(c_double), allocatable :: dev_facs(:), sigmas(:)
    real(c_double), allocatable :: triangle_(:, :)
    integer(c_int) :: proc_dist
    integer(c_int) :: outlier_rowidx

    outlier_rowidx = sim_in%outlier_rowidx
    proc_dist = sim_in%proc_dist
    factor = sim_in%factor
    triangle_ = triangle_in%data
    dev_facs = triangle_in%dev_facs
    sigmas = triangle_in%sigmas

    if (factor == 1) then
      origin_outlier_ = triangle_
      return
    end if

    n_dev = size(triangle_, dim=1)
    origin_outlier_ = triangle_

    do j = 2, n_dev + 1 - outlier_rowidx
      if (proc_dist == NORMAL) then
        mean = factor * dev_facs(j - 1) * origin_outlier_(outlier_rowidx, j - 1)
        sd = sigmas(j - 1) * sqrt(origin_outlier_(outlier_rowidx, j - 1))
        origin_outlier_(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)

      else if (proc_dist == GAMMA) then
        shape = factor * dev_facs(j - 1)**2 * origin_outlier_(outlier_rowidx, j - 1) / sigmas(j - 1)**2
        scale = sigmas(j - 1)**2 / factor * dev_facs(j - 1)
        origin_outlier_(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
      end if
    end do

    origin_outlier = triangle(n_dev, origin_outlier_)
  end function origin_outlier

end module mack
