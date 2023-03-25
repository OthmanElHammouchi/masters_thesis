module mod_mack

  use, intrinsic :: iso_c_binding
  use, intrinsic :: omp_lib
  use mod_global
  use mod_helpers
  use mod_interface

  implicit none

  ! Global simulation variables
  integer(c_int) :: boot_type, proc_dist, resids_type
  logical(c_bool) :: cond
  !$omp threadprivate(boot_type, resids_type, proc_dist, cond)

  ! Shared containers
  real(c_double), allocatable :: scale_facs(:, :), sigma_jack(:, :)
  real(c_double), allocatable :: log_normal_means(:), log_normal_sigmas(:)
  real(c_double), allocatable :: log_normal_shift(:)
  real(c_double), allocatable :: dev_facs(:), sigmas(:)
  real(c_double), allocatable :: resids(:, :), resampled_resids(:, :)
  logical(c_bool), allocatable :: mask(:, :)
  real(c_double), allocatable :: sim_table(:, :)
  !$omp threadprivate(scale_facs, sigma_jack, log_normal_means, log_normal_sigmas, &
  !$omp& log_normal_shift, dev_facs, sigmas, resids, resampled_resids, mask, sim_table)

contains

  subroutine sim(n_dev, triangle, sim_type, n_boot, n_res, m_res, res, show_progress) bind(c, name="mack_sim_")
    integer(c_int), intent(in), value :: n_dev, n_boot, sim_type
    integer(c_int), intent(in), value :: n_res, m_res
    logical(c_bool), intent(in), value :: show_progress
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: res(n_res, m_res)

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

    logical :: test
    integer(c_int) :: status
    logical(c_bool) :: resids_mask(n_dev - 1, n_dev - 1)

    call fit(triangle, dev_facs_original, sigmas_original, use_mask=.false., compute_resids=.false.)

    n_threads = init_omp()
    rng = init_rng(n_threads, 42)
    n_sim = size(sim_table, 1)
    pgbar = pgbar_create(n_sim, 1)

    !$omp parallel num_threads(n_threads) &
    !$omp& copyin(rng, boot_type, proc_dist, cond, resids_type, sim_table) &
    !$omp& private(reserve, i_sim, outlier_rowidx, outlier_colidx, test, status, resids_mask, &
    !$omp& outlier_diagidx, excl_rowidx, excl_diagidx, factor, triangle_sim) &
    !$omp& firstprivate(n_boot, n_dev, dev_facs_original, sigmas_original, &
    !$omp& m_res, n_sim) &
    !$omp& shared(res, pgbar)

    ! Allocate shared containers
    allocate(dev_facs(n_dev - 1), source=0._c_double)
    allocate(sigmas(n_dev - 1), source=0._c_double)
    allocate(scale_facs(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(sigma_jack(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(log_normal_shift(n_dev - 1), source=0._c_double)
    allocate(log_normal_means(n_dev - 1), source=0._c_double)
    allocate(log_normal_sigmas(n_dev - 1), source=0._c_double)
    allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(resampled_resids(n_dev - 1, n_dev - 1), source=0._c_double)

    allocate(mask(n_dev, n_dev))
    allocate(reserve(n_boot), source=0._c_double)

    i_thread = omp_get_thread_num()

    reserve = 0
    select case (sim_type)
     case (SINGLE)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        ! Set simulation variables
        boot_type = int(sim_table(i_sim, 6))
        proc_dist = int(sim_table(i_sim, 7))
        cond = int(sim_table(i_sim, 8))
        resids_type = int(sim_table(i_sim, 9))

        factor = sim_table(i_sim, 3)
        outlier_rowidx = int(sim_table(i_sim, 1))
        outlier_colidx = int(sim_table(i_sim, 2))

        excl_rowidx = sim_table(i_sim, 4)
        excl_colidx = sim_table(i_sim, 5)
        mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            mask(i, j) = .false.
          end do
        end do
        mask(excl_rowidx, excl_colidx) = .false.
        if (boot_type == PARAMETRIC) mask(1, n_dev) = .true.

        triangle_sim = single_outlier(triangle, dev_facs_original, sigmas_original, outlier_rowidx, outlier_colidx, factor)
        call boot(n_dev, triangle_sim, n_boot, reserve, mask)

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_res) = spread(sim_table(i_sim, 1:(m_res - 1)), 1, n_boot)
        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), m_res) = reserve
      end do
      !$omp end do

     case (CALENDAR)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        ! Set simulation variables
        boot_type = int(sim_table(i_sim, 4))
        proc_dist = int(sim_table(i_sim, 5))
        cond = int(sim_table(i_sim, 6))
        resids_type = int(sim_table(i_sim, 7))

        factor = sim_table(i_sim, 2)
        outlier_diagidx = int(sim_table(i_sim, 1))
        excl_diagidx = int(sim_table(i_sim, 3))

        mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            mask(i, j) = .false.
          end do
        end do
        do j = 2, n_dev
          i = n_dev + 2 - excl_diagidx - j
          if (i <= 0) cycle
          mask(i, j) = .false.
        end do
        if (boot_type == PARAMETRIC) mask(1, n_dev) = .true.

        triangle_sim = calendar_outlier(triangle, dev_facs_original, sigmas_original, outlier_diagidx, factor)

        ! test = outlier_diagidx == excl_diagidx .and. outlier_diagidx == 1 .and. factor == 1.5 .and. boot_type == 3

        ! if (test) then
        !   print *, raise(2)
        !   call fit(triangle_sim, dev_facs, sigmas, use_mask = .false., compute_resids = .true.)
        !   call print_vector(dev_facs)
        !   call print_vector(sigmas)
        !   resids_mask = mask(1:(n_dev - 1), 2:n_dev)
        !   resampled_resids = sample(resids, resids_mask)
        !   call resample(triangle_sim, dev_facs, sigmas, status)
        !   call print_vector(dev_facs)
        !   call print_vector(sigmas)
        ! end if

        call boot(n_dev, triangle_sim, n_boot, reserve, mask)

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_res) = spread(sim_table(i_sim, 1:(m_res - 1)), 1, n_boot)
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_res) = reserve
      end do
      !$omp end do

     case (ORIGIN)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        ! Set simulation variables
        boot_type = int(sim_table(i_sim, 4))
        proc_dist = int(sim_table(i_sim, 5))
        cond = int(sim_table(i_sim, 6))
        resids_type = int(sim_table(i_sim, 7))


        factor = sim_table(i_sim, 2)
        outlier_rowidx = int(sim_table(i_sim, 1))
        excl_rowidx = int(sim_table(i_sim, 3))

        mask = .true.
        do j = 2, n_dev + 1 - excl_rowidx
          mask(excl_rowidx, j) = .false.
        end do
        if (boot_type == PARAMETRIC) mask(1, n_dev) = .true.

        triangle_sim = origin_outlier(triangle, dev_facs_original, sigmas_original, outlier_rowidx, factor)
        call boot(n_dev, triangle_sim, n_boot, reserve, mask)

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_res) = spread(sim_table(i_sim, 1:(m_res - 1)), 1, n_boot)
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_res) = reserve
      end do
      !$omp end do
    end select

    ! Free shared containers
    deallocate(log_normal_shift)
    deallocate(log_normal_means)
    deallocate(log_normal_sigmas)
    deallocate(scale_facs)
    deallocate(sigma_jack)
    deallocate(dev_facs)
    deallocate(sigmas)
    deallocate(resids)
    deallocate(resampled_resids)
    deallocate(mask)
    deallocate(reserve)
    !$omp end parallel
    deallocate(sim_table)
  end subroutine sim

  subroutine build_sim_table_(sim_type, n_dev, n_boot, n_factors, factors, &
    n_boot_types, boot_types, n_proc_dists, proc_dists, n_conds, conds, &
    n_resids_types, resids_types, n_res, m_res) bind(c)
    integer(c_int), intent(in), value :: sim_type
    integer(c_int), intent(in), value :: n_dev, n_boot, n_factors, n_boot_types
    integer(c_int), intent(in), value :: n_proc_dists, n_conds, n_resids_types
    integer(c_int), intent(out) :: n_res, m_res
    real(c_double), intent(in) :: factors(n_factors)
    integer(c_int), intent(in) :: boot_types(n_boot_types), proc_dists(n_proc_dists)
    integer(c_int), intent(in) :: resids_types(n_resids_types)
    logical(c_bool), intent(in) :: conds(n_conds)

    integer(c_int) :: n_sim, m_sim, n_sim_, m_sim_, n_outliers, n_excl
    integer(c_int) :: i, j, k, i1, i2, i3, i4, i5, i6, i7
    integer(c_int) :: n_rows
    integer(c_int), allocatable :: outliers(:, :), excl(:, :)
    real(c_double), allocatable :: sim_table_(:, :)

    select case (sim_type)
     case (SINGLE)
      n_outliers = (n_dev ** 2 - n_dev) / 2 - 1
      n_excl = n_outliers
      n_sim_ = n_outliers * n_excl * n_factors * n_boot_types * n_proc_dists * n_conds * n_resids_types
      m_sim_ = 9

      allocate(outliers(n_outliers, 2), source=0)
      allocate(excl(n_excl, 2), source=0)
      allocate(sim_table_(n_sim_, m_sim_), source=0._c_double)

      k = 1
      do j = 2, n_dev - 1
        n_rows = n_dev + 1 - j
        do i = 1, n_rows
          outliers(k, :) = [i, j]
          k = k + 1
        end do
      end do
      excl = outliers

      k = 1
      do i1 = 1, n_outliers
        do i2 = 1, n_factors
          do i3 = 1, n_excl
            do i4 = 1, n_boot_types
              proc_loop_single: do i5 = 1, n_proc_dists
                cond_loop_single: do i6 = 1, n_conds
                  do i7 = 1, n_resids_types
                    sim_table_(k, 1:2) = outliers(i1, :)
                    sim_table_(k, 3) = factors(i2)
                    sim_table_(k, 4:5) = excl(i3, :)
                    sim_table_(k, 6) = boot_types(i4)
                    sim_table_(k, 7) = proc_dists(i5)

                    if (sim_table_(k, 6) == PAIRS) then
                      sim_table_(k, 8) = NONE
                      sim_table_(k, 9) = NONE
                      k = k + 1
                      cycle proc_loop_single

                    else if (sim_table_(k, 6) == PARAMETRIC) then
                      sim_table_(k, 8) = merge(1._c_double, 0._c_double, conds(i6))
                      sim_table_(k, 9) = NONE
                      k = k + 1
                      cycle cond_loop_single
                    else
                      sim_table_(k, 8) = merge(1._c_double, 0._c_double, conds(i6))
                      sim_table_(k, 9) = resids_types(i7)
                      k = k + 1
                    end if
                  end do
                end do cond_loop_single
              end do proc_loop_single 
            end do 
          end do 
        end do 
      end do

     case (CALENDAR)
      n_outliers = n_dev - 1
      n_excl = n_outliers
      n_sim_ = n_outliers * n_excl * n_factors * n_boot_types * n_proc_dists * n_conds * n_resids_types
      m_sim_ = 7

      allocate(sim_table_(n_sim_, m_sim_), source=0._c_double)

      k = 1
      do i1 = 1, n_outliers
        do i2 = 1, n_factors
          do i3 = 1, n_excl
            do i4 = 1, n_boot_types
              proc_loop_cal: do i5 = 1, n_proc_dists
                cond_loop_cal: do i6 = 1, n_conds
                  do i7 = 1, n_resids_types
                    sim_table_(k, 1) = i1
                    sim_table_(k, 2) = factors(i2)
                    sim_table_(k, 3) = i3
                    sim_table_(k, 4) = boot_types(i4)
                    sim_table_(k, 5) = proc_dists(i5)

                    if (sim_table_(k, 4) == PAIRS) then
                      sim_table_(k, 6) = NONE
                      sim_table_(k, 7) = NONE
                      k = k + 1
                      cycle proc_loop_cal

                    else if (sim_table_(k, 4) == PARAMETRIC) then
                      sim_table_(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                      sim_table_(k, 7) = NONE
                      k = k + 1
                      cycle cond_loop_cal

                    else
                      sim_table_(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                      sim_table_(k, 7) = resids_types(i7)
                      k = k + 1
                    end if
                  end do
                end do cond_loop_cal
              end do proc_loop_cal
            end do
          end do 
        end do
      end do

     case (ORIGIN)
      n_outliers = n_dev
      n_excl = n_outliers
      n_sim_ = n_outliers * n_excl * n_factors * n_boot_types * n_proc_dists * n_conds * n_resids_types
      m_sim_ = 7

      allocate(sim_table_(n_sim_, m_sim_), source=0._c_double)

      k = 1
      do i1 = 1, n_outliers
        do i2 = 1, n_factors
          do i3 = 1, n_excl
            do i4 = 1, n_boot_types
              proc_loop_origin: do i5 = 1, n_proc_dists
                cond_loop_origin: do i6 = 1, n_conds
                  do i7 = 1, n_resids_types
                    sim_table_(k, 1) = i1
                    sim_table_(k, 2) = factors(i2)
                    sim_table_(k, 3) = i3
                    sim_table_(k, 4) = boot_types(i4)
                    sim_table_(k, 5) = proc_dists(i5)

                    if (sim_table_(k, 4) == PAIRS) then
                      sim_table_(k, 6) = NONE
                      sim_table_(k, 7) = NONE
                      k = k + 1
                      cycle proc_loop_origin

                    else if (sim_table_(k, 4) == PARAMETRIC) then
                      sim_table_(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                      sim_table_(k, 7) = NONE
                      k = k + 1
                      cycle cond_loop_origin

                    else
                      sim_table_(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                      sim_table_(k, 7) = resids_types(i7)
                      k = k + 1
                    end if
                  end do
                end do cond_loop_origin
              end do proc_loop_origin
            end do
          end do
        end do
      end do
    end select

    n_sim = k - 1
    m_sim = m_sim_

    allocate(sim_table(n_sim, m_sim))
    sim_table = sim_table_(1:n_sim, :)
    deallocate(sim_table_)

    n_res = n_boot * n_sim
    m_res = m_sim + 1
  end subroutine build_sim_table_

  subroutine boot(n_dev, triangle, n_boot, reserve, mask)
    integer(c_int), intent(in) :: n_boot, n_dev
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    logical(c_bool), intent(in) :: mask(:, :)

    real(c_double) :: reserve(n_boot)

    integer(c_int) :: i, j, i_diag, i_boot, n_resids, n_excl_resids
    real(c_double) :: dev_facs_boot(n_dev - 1), sigmas_boot(n_dev - 1)
    logical(c_bool) :: resids_mask(n_dev - 1, n_dev - 1)
    real(c_double) :: triangle_boot(n_dev, n_dev)
    real(c_double) :: latest(n_dev)
    real(c_double) :: mean, sd, shape, scale
    integer(c_int) :: status

    if (boot_type == PARAMETRIC) then
      call fit(triangle, dev_facs, sigmas, use_mask=.true., compute_resids=.false.)
    else if (boot_type == RESID) then
      call fit(triangle, dev_facs, sigmas, use_mask=.false., compute_resids=.true.)
    end if

    i_boot = 1
    main_loop: do while (i_boot <= n_boot)
      ! Simulate new parameters
      if (boot_type == RESID) then
        resids_mask = mask(1:(n_dev - 1), 2:n_dev)
        resampled_resids = sample(resids, resids_mask)
      end if

      call resample(triangle, dev_facs_boot, sigmas_boot, status)
      if (status == FAILURE) cycle main_loop

      ! Simulate process error
      triangle_boot = triangle

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

  ! Return fitted development factors and dispersion parameters, changing
  ! the value of the shared variables in the process.
  subroutine fit(triangle, dev_facs, sigmas, use_mask, compute_resids)
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(out) :: dev_facs(:), sigmas(:)
    logical, intent(in) :: compute_resids, use_mask

    integer(c_int) :: i, j, n_rows, n_dev, n_pts_col, n_resids
    real(c_double), allocatable :: indiv_dev_facs(:, :)
    logical(c_bool), allocatable :: col_mask(:)
    real(c_double) :: resids_mean, dev_fac_jack

    n_dev = size(triangle, 1)
    allocate(indiv_dev_facs(n_dev - 1, n_dev - 1), source=0._c_double)

    if (use_mask) then
      do j = 1, n_dev - 1
        n_rows = n_dev - j
        do i = 1, n_rows
          indiv_dev_facs(i, j) = triangle(i, j + 1) / triangle(i, j)
        end do
        col_mask = mask(1:n_rows, j + 1)
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

    if (compute_resids) then
      if (resids_type /= LOGNORMAL) then
        scale_facs = 0
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            scale_facs(i, j) = sqrt(1 - triangle(i, j) / sum(triangle(1:n_rows, j)))
          end do
        end do
        scale_facs(1, n_dev - 1) = 1
      end if

      n_resids = (n_dev ** 2 - n_dev) / 2 - 1
      resids = 0
      select case (resids_type)
       case (STANDARDISED)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            resids(i, j) = (indiv_dev_facs(i, j) - dev_facs(j)) * sqrt(triangle(i, j)) / (sigmas(j) * scale_facs(i, j))
          end do
        end do

       case (MODIFIED)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            resids(i, j) = (indiv_dev_facs(i, j) - dev_facs(j)) * sqrt(triangle(i, j)) / scale_facs(i, j)
          end do
        end do

       case (STUDENTISED)
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
        do j = 1, n_dev
          n_rows = n_dev - j
          do i = 1, n_rows
            log_normal_shift(i) = dev_facs(j) * sqrt(triangle(i, j)) / sigmas(j)
            log_normal_sigmas(i) = sqrt(log(1 + 1 / log_normal_shift(i) ** 2))
            log_normal_means(i) = log(log_normal_shift(i)) - log_normal_sigmas(i) ** 2 / 2
            resids(i, j) = (triangle(i, j + 1) - dev_facs(j) * triangle(i, j)) / (sigmas(j) * sqrt(triangle(i, j)))
            resids(i, j) = (log(resids(i, j) + log_normal_shift(i)) - log_normal_means(i)) / log_normal_sigmas(i)
          end do
        end do
      end select

      resids(1, n_dev - 1) = 0

      if (resids_type /= STUDENTISED) then
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

  ! Return simulated development factors and dispersion parameters.
  subroutine resample(triangle, dev_facs_boot, sigmas_boot, status)
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(out) :: dev_facs_boot(:), sigmas_boot(:)
    integer(c_int), intent(out) :: status

    real(c_double), allocatable :: resampled_triangle(:, :)
    integer(c_int) :: i, j, n_rows, n_dev
    real(c_double) :: mean, sd, shape, scale

    real(c_double), allocatable :: resampled_pairs(:, :)
    integer(c_int), allocatable :: pair_indices(:), resampled_pair_indices(:)

    n_dev = size(triangle, 1)
    allocate(resampled_triangle(n_dev, n_dev), source=0._c_double)

    status = SUCCESS
    select case (boot_type)
     case (PARAMETRIC)
      if(cond) then
        resampled_triangle(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (proc_dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle(i, j - 1))
              resampled_triangle(i, j) = rnorm_par(rng, i_thread, mean, sd)
            else if (proc_dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              resampled_triangle(i, j) = rgamma_par(rng, i_thread, shape, scale)
            end if
          end do
        end do

      else
        resampled_triangle(:, 1) = triangle(:, 1)
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

      call fit(resampled_triangle, dev_facs_boot, sigmas_boot, use_mask=.false., compute_resids=.false.)

     case (RESID)
      if (cond) then
        resampled_triangle(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == STANDARDISED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigmas(j) * scale_facs(i, j) * sqrt(triangle(i, j)) * resampled_resids(i, j)
            else if (resids_type == MODIFIED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                scale_facs(i, j) * sqrt(triangle(i, j)) * resampled_resids(i, j)
            else if (resids_type == STUDENTISED) then
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigma_jack(i, j) * scale_facs(i, j) * sqrt(triangle(i, j)) * resampled_resids(i, j)
            else if (resids_type == LOGNORMAL) then
              resampled_triangle(i, j + 1) = exp(resampled_resids(i, j) * log_normal_sigmas(i) + &
               log_normal_means(i)) - log_normal_shift(i)
              resampled_triangle(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigmas(j) * sqrt(triangle(i, j)) * resampled_triangle(i, j + 1)
            end if
          end do
        end do

      else
        resampled_triangle(:, 1) = triangle(:, 1)
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

      call fit(resampled_triangle, dev_facs_boot, sigmas_boot, use_mask=.false., compute_resids=.false.)

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
        resampled_pairs(1:n_rows, :) = triangle(resampled_pair_indices(1:n_rows), j:(j + 1))

        dev_facs_boot(j) = sum(resampled_pairs(1:n_rows, 2)) / sum(resampled_pairs(1:n_rows, 1))
        if (j < n_dev - 1) then
          sigmas_boot(j) = sqrt(sum(triangle(1:n_rows, j) * (triangle(1:n_rows, j + 1) / triangle(1:n_rows, j) - &
            dev_facs_boot(j)) ** 2) / (n_rows - 1))
        else
          sigmas_boot(j) = extrapolate_sigma(sigmas_boot, j)
        end if
      end do
    end select
  end subroutine resample

  pure real(c_double) function extrapolate_sigma(sigmas, col)
    real(c_double), intent(in) :: sigmas(:)
    integer(c_int), intent(in) :: col

    extrapolate_sigma = sqrt(min(sigmas(col - 1) ** 2, sigmas(col - 2) ** 2, sigmas(col - 1) ** 4 / sigmas(col - 2) ** 2))
  end function extrapolate_sigma

  subroutine mack_boot_(n_dev, triangle, boot_type_, proc_dist_, cond_, &
    resids_type_, n_boot, reserve) bind(c)
    integer(c_int), intent(in), value :: n_boot, n_dev, proc_dist_, resids_type_, boot_type_
    logical(c_bool), intent(in), value :: cond_
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: reserve(n_boot)

    allocate(dev_facs(n_dev - 1), source=0._c_double)
    allocate(sigmas(n_dev - 1), source=0._c_double)
    allocate(scale_facs(n_dev, n_dev), source=0._c_double)
    allocate(sigma_jack(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(log_normal_shift(n_dev - 1), source=0._c_double)
    allocate(log_normal_means(n_dev - 1), source=0._c_double)
    allocate(log_normal_sigmas(n_dev - 1), source=0._c_double)
    allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)

    boot_type = boot_type_
    proc_dist = proc_dist_
    cond = cond_
    resids_type = resids_type_

    allocate(mask(n_dev, n_dev))
    mask = .true.

    rng = init_rng(1, 42)
    i_thread = 0

    call boot(n_dev, triangle, n_boot, reserve, mask)
  end subroutine mack_boot_

  function single_outlier(triangle, dev_facs, sigmas, outlier_rowidx, outlier_colidx, factor)
    integer(c_int), intent(in):: outlier_rowidx, outlier_colidx
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in):: dev_facs(:), sigmas(:)

    real(c_double) :: factor
    real(c_double), allocatable:: single_outlier(:, :)
    real(c_double) :: shape, scale, mean, sd
    integer(c_int) :: n_dev, i, j

    if (factor == 1) then
      single_outlier = triangle
      return
    end if

    n_dev = size(triangle, 1)
    allocate(single_outlier(n_dev, n_dev), source=0._c_double)
    single_outlier(:, 1) = triangle(:, 1)

    if (proc_dist == NORMAL) then
      do j = 2, n_dev
        do i = 1, n_dev + 1 - j
          if (i == outlier_rowidx) cycle
          mean = dev_facs(j - 1) * single_outlier(i, j - 1)
          sd = sigmas(j - 1) * sqrt(single_outlier(i, j - 1))
          single_outlier(i, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end do

      if (outlier_colidx > 2) then
        do j = 2, outlier_colidx - 1
          mean = dev_facs(j - 1) * single_outlier(outlier_rowidx, j - 1)
          sd = sigmas(j - 1) * sqrt(single_outlier(outlier_rowidx, j - 1))
          single_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end if

      mean = factor * dev_facs(outlier_colidx - 1) * single_outlier(outlier_rowidx, outlier_colidx - 1)
      sd = sigmas(outlier_colidx - 1) * sqrt(single_outlier(outlier_rowidx, outlier_colidx - 1))
      single_outlier(outlier_rowidx, outlier_colidx) = rnorm_par(rng, i_thread, mean, sd)

      if (outlier_colidx < n_dev) then
        do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
          mean = dev_facs(j - 1) * single_outlier(outlier_rowidx, j - 1)
          sd = sigmas(j - 1) * sqrt(single_outlier(outlier_rowidx, j - 1))
          single_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
        end do
      end if

    else if (proc_dist == GAMMA) then
      do j = 2, n_dev
        do i = 1, n_dev + 1 - j
          if (i == outlier_rowidx) cycle
          shape = dev_facs(j - 1)**2 * single_outlier(i, j - 1) / sigmas(j - 1)**2
          scale = sigmas(j - 1)**2 / dev_facs(j - 1)
          single_outlier(i, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end do

      if (outlier_colidx > 2) then
        do j = 2, outlier_colidx - 1
          shape = dev_facs(j - 1)**2 * single_outlier(outlier_rowidx, j - 1) / sigmas(j - 1)**2
          scale = sigmas(j - 1)**2 / dev_facs(j - 1)
          single_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end if

      shape = dev_facs(outlier_colidx - 1)**2 * single_outlier(outlier_rowidx, outlier_colidx - 1) / &
        sigmas(outlier_colidx - 1)**2
      scale = sigmas(outlier_colidx - 1)**2 / dev_facs(outlier_colidx - 1)
      single_outlier(outlier_rowidx, outlier_colidx) = rgamma_par(rng, i_thread, shape, scale)

      if (outlier_colidx < n_dev) then
        do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
          shape = dev_facs(j - 1)**2 * single_outlier(outlier_rowidx, j - 1) / sigmas(j - 1)**2
          scale = sigmas(j - 1)**2 / dev_facs(j - 1)
          single_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
        end do
      end if
    end if
  end function single_outlier

  function calendar_outlier(triangle, dev_facs, sigmas, outlier_diagidx, factor)
    integer(c_int), intent(in) :: outlier_diagidx
    real(c_double), intent(in):: factor
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in):: dev_facs(:), sigmas(:)

    integer(c_int) :: i, j, n_dev, n_cols
    real(c_double), allocatable :: calendar_outlier(:, :)
    real(c_double) :: shape, scale, mean, sd

    if (factor == 1) then
      calendar_outlier = triangle
      return
    end if

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

        if (proc_dist == NORMAL) then
          mean = factor * dev_facs(n_cols - 1) * calendar_outlier(i, n_cols - 1)
          sd = sigmas(n_cols - 1) * sqrt(calendar_outlier(i, n_cols - 1))
          calendar_outlier(i, n_cols) = rnorm_par(rng, i_thread, mean, sd)
          do j = n_cols + 1, n_dev + 1 - i
            mean = dev_facs(j - 1) * calendar_outlier(i, j - 1)
            sd = sigmas(j - 1) * sqrt(calendar_outlier(i, j - 1))
            calendar_outlier(i, j) = rnorm_par(rng, i_thread, mean, sd)
          end do

        else if (proc_dist == GAMMA) then
          shape = factor * dev_facs(n_cols - 1)**2 * calendar_outlier(i, n_cols - 1) / sigmas(n_cols - 1)**2
          scale = sigmas(n_cols - 1)**2 / factor * dev_facs(n_cols - 1)
          calendar_outlier(i, n_cols) = rgamma_par(rng, i_thread, shape, scale)
          do j = n_cols + 1, n_dev + 1 - i
            shape = dev_facs(j - 1)**2 * calendar_outlier(i, j - 1) / sigmas(j - 1)**2
            scale = sigmas(j - 1)**2 / dev_facs(j - 1)
            calendar_outlier(i, j) = rgamma_par(rng, i_thread, shape, scale)
          end do
        end if
      end if
    end do
  end function calendar_outlier

  function origin_outlier(triangle, dev_facs, sigmas, outlier_rowidx, factor)
    integer(c_int), intent(in):: outlier_rowidx
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in) :: factor
    real(c_double), intent(in):: dev_facs(:), sigmas(:)

    real(c_double) :: shape, scale, mean, sd
    real(c_double), allocatable:: origin_outlier(:, :)
    integer(c_int) :: n_dev, j

    if (factor == 1) then
      origin_outlier = triangle
      return
    end if

    n_dev = size(triangle, dim=1)
    origin_outlier = triangle

    do j = 2, n_dev + 1 - outlier_rowidx
      if (proc_dist == NORMAL) then
        mean = factor * dev_facs(j - 1) * origin_outlier(outlier_rowidx, j - 1)
        sd = sigmas(j - 1) * sqrt(origin_outlier(outlier_rowidx, j - 1))
        origin_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)

      else if (proc_dist == GAMMA) then
        shape = factor * dev_facs(j - 1)**2 * origin_outlier(outlier_rowidx, j - 1) / sigmas(j - 1)**2
        scale = sigmas(j - 1)**2 / factor * dev_facs(j - 1)
        origin_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
      end if
    end do
  end function origin_outlier

end module mod_mack
