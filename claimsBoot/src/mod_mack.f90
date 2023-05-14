module mod_mack

  use, intrinsic :: iso_c_binding
  use, intrinsic :: omp_lib
  use mod_global
  use mod_helpers
  use mod_interface

  implicit none

  ! Simulation configuration
  integer(c_int) :: boot_type, dist, resids_type
  logical(c_bool) :: cond
  !$omp threadprivate(boot_type, resids_type, dist, cond)

  ! Shared variables
  real(c_double), allocatable :: scale_facs(:, :), sigmas_jack(:, :)
  real(c_double), allocatable :: ln_means(:, :), ln_sigmas(:, :)
  real(c_double), allocatable :: ln_shifts(:, :)
  real(c_double), allocatable :: dev_facs(:), sigmas(:)
  real(c_double), allocatable :: resids(:, :), resampled_resids(:, :)
  logical(c_bool), allocatable :: obs_mask(:, :), resids_mask(:, :)
  !$omp threadprivate(scale_facs, sigmas_jack, ln_means, ln_sigmas, &
  !$omp& ln_shifts, dev_facs, sigmas, resids, resampled_resids, obs_mask, resids_mask)

contains

  subroutine mack_sim(n_dev, triangle, sim_type, boot_type_, sim_dist, n_boot, &
    n_conf, m_conf, conf, res, show_progress, seed) bind(c, name="mack_sim")
    integer(c_int), intent(in), value :: n_dev, n_boot
    integer(c_int), intent(in), value :: n_conf, m_conf
    integer(c_int), intent(in), value :: sim_type, boot_type_, sim_dist
    logical(c_bool), intent(in), value :: show_progress
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(in) :: conf(n_conf, m_conf)
    real(c_double), intent(out) :: res(n_boot * n_conf, m_conf + 1)
    integer(c_int), intent(in), value :: seed

    integer(c_int) :: i, j, i_sim
    integer(c_int) :: n_excl, n_sim
    integer(c_int) :: excl_diagidx, excl_rowidx, excl_colidx
    integer(c_int) :: outlier_rowidx, outlier_colidx, outlier_diagidx
    real(c_double) :: mean_factor, sd_factor
    real(c_double) :: dev_facs_original(n_dev - 1), sigmas_original(n_dev - 1)
    real(c_double) :: triangle_sim(n_dev, n_dev)
    real(c_double), allocatable :: reserve(:)
    type(c_ptr) :: pgbar

    integer(c_int) :: status
    logical(c_bool) :: print_triangles

    call fit(triangle, dev_facs_original, sigmas_original, use_mask=.false., compute_resids=.false.)

    if (first_call) then
      n_threads = init_omp()
      rng = init_rng(n_threads, seed)
      first_call = .false._c_bool
    end if

    n_sim = n_conf
    pgbar = pgbar_create(n_sim, 1)

    !$omp parallel num_threads(n_threads) &
    !$omp& copyin(rng, boot_type, dist, cond, resids_type) &
    !$omp& private(reserve, i_sim, outlier_rowidx, outlier_colidx, status, &
    !$omp& outlier_diagidx, excl_rowidx, excl_diagidx, mean_factor, sd_factor, triangle_sim) &
    !$omp& firstprivate(n_boot, n_dev, dev_facs_original, sigmas_original, &
    !$omp& n_conf, m_conf, n_sim, sim_dist, sim_type) &
    !$omp& shared(conf, pgbar)

    ! Allocate shared containers
    allocate(dev_facs(n_dev - 1), source=0._c_double)
    allocate(sigmas(n_dev - 1), source=0._c_double)
    allocate(scale_facs(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(sigmas_jack(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(ln_shifts(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(ln_means(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(ln_sigmas(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(resampled_resids(n_dev - 1, n_dev - 1), source=0._c_double)

    allocate(obs_mask(n_dev, n_dev))
    allocate(resids_mask(n_dev - 1, n_dev - 1))
    allocate(reserve(n_boot), source=0._c_double)

    i_thread = omp_get_thread_num()
    boot_type = boot_type_
    print_triangles = .false.

    reserve = 0
    select case (sim_type)
     case (SINGLE)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        outlier_rowidx = int(conf(i_sim, 1))
        outlier_colidx = int(conf(i_sim, 2))
        mean_factor = conf(i_sim, 3)
        sd_factor = conf(i_sim, 4)
        excl_rowidx = int(conf(i_sim, 5))
        excl_colidx = int(conf(i_sim, 6))

        if (boot_type /= PAIRS) then
          if (boot_type == RESID) then
            resids_type = int(conf(i_sim, 7))
          else if (boot_type == PARAM) then
            dist = int(conf(i_sim, 7))
          end if
          cond = int(conf(i_sim, 8))
        end if

        obs_mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            obs_mask(i, j) = .false.
          end do
        end do
        obs_mask(excl_rowidx, excl_colidx) = .false.
        if (boot_type == PARAM) obs_mask(1, n_dev) = .true.

        resids_mask = obs_mask(1:(n_dev - 1), 2:n_dev)
        if (resids_type /= LOGNORMAL) then
          resids_mask(:, n_dev - 1) = .false.
          resids_mask(:, n_dev - 2) = .false.
        end if

        if (boot_type == RESID .and. outlier_colidx == 2 .and. outlier_rowidx == 1 &
        .and. excl_colidx == outlier_colidx .and. excl_rowidx == outlier_rowidx) print_triangles = .true.

        triangle_sim = single_outlier(triangle, dev_facs_original, sigmas_original, &
          outlier_rowidx, outlier_colidx, mean_factor, sd_factor, sim_dist)
        call boot(n_dev, triangle_sim, n_boot, reserve, obs_mask, print_triangles)

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_conf) = spread(conf(i_sim, 1:m_conf), 1, n_boot)
        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), m_conf + 1) = reserve
      end do
      !$omp end do

     case (CALENDAR)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)
        outlier_diagidx = int(conf(i_sim, 1))
        mean_factor = conf(i_sim, 2)
        sd_factor = conf(i_sim, 3)
        excl_diagidx = int(conf(i_sim, 4))

        if (boot_type /= PAIRS) then
          if (boot_type == RESID) then
            resids_type = int(conf(i_sim, 5))
          else if (boot_type == PARAM) then
            dist = int(conf(i_sim, 5))
          end if
          cond = int(conf(i_sim, 6))
        end if

        obs_mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            obs_mask(i, j) = .false.
          end do
        end do
        do j = 2, n_dev
          i = n_dev + 2 - excl_diagidx - j
          if (i <= 0) cycle
          obs_mask(i, j) = .false.
        end do
        if (boot_type == PARAM) obs_mask(1, n_dev) = .true.

        resids_mask = obs_mask(1:(n_dev - 1), 2:n_dev)
        resids_mask(:, n_dev - 1) = .false.
        resids_mask(:, n_dev - 2) = .false.

        triangle_sim = calendar_outlier(triangle, dev_facs_original, sigmas_original, outlier_diagidx, &
          mean_factor, sd_factor, sim_dist)

        call boot(n_dev, triangle_sim, n_boot, reserve, obs_mask, print_triangles)

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_conf) = spread(conf(i_sim, 1:m_conf), 1, n_boot)
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_conf + 1) = reserve
      end do
      !$omp end do

     case (ORIGIN)
      !$omp do schedule(dynamic, 25)
      do i_sim = 1, n_sim
        if (show_progress) call pgbar_incr(pgbar)

        outlier_rowidx = int(conf(i_sim, 1))
        mean_factor = conf(i_sim, 2)
        sd_factor = conf(i_sim, 3)
        excl_rowidx = int(conf(i_sim, 4))

        if (boot_type /= PAIRS) then
          if (boot_type == RESID) then
            resids_type = int(conf(i_sim, 5))
          else if (boot_type == PARAM) then
            dist = int(conf(i_sim, 5))
          end if
          cond = int(conf(i_sim, 6))
        end if

        obs_mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            obs_mask(i, j) = .false.
          end do
        end do
        do j = 2, n_dev + 1 - excl_rowidx
          obs_mask(excl_rowidx, j) = .false.
        end do
        if (boot_type == PARAM) obs_mask(1, n_dev) = .true. ! avoid 0/0

        resids_mask = obs_mask(1:(n_dev - 1), 2:n_dev)
        resids_mask(:, n_dev - 1) = .false.
        resids_mask(:, n_dev - 2) = .false.

        triangle_sim = origin_outlier(triangle, dev_facs_original, sigmas_original, outlier_rowidx, &
          mean_factor, sd_factor, sim_dist)
        call boot(n_dev, triangle_sim, n_boot, reserve, obs_mask, print_triangles)

        res(((i_sim - 1) * n_boot + 1):(i_sim * n_boot), 1:m_conf) = spread(conf(i_sim, 1:m_conf), 1, n_boot)
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_conf + 1) = reserve
      end do
      !$omp end do
    end select

    ! Free shared containers
    deallocate(ln_shifts)
    deallocate(ln_means)
    deallocate(ln_sigmas)
    deallocate(scale_facs)
    deallocate(sigmas_jack)
    deallocate(dev_facs)
    deallocate(sigmas)
    deallocate(resids)
    deallocate(resampled_resids)
    deallocate(obs_mask)
    deallocate(resids_mask)
    deallocate(reserve)
    !$omp end parallel
  end subroutine mack_sim

  subroutine boot(n_dev, triangle, n_boot, reserve, mask, print)
    integer(c_int), intent(in) :: n_boot, n_dev
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    logical(c_bool), intent(in) :: mask(:, :)
    real(c_double), intent(inout) :: reserve(n_boot)
    logical(c_bool), intent(in) :: print

    real(c_double) :: res

    integer(c_int) :: i_boot, n_resids, n_excl_resids, n_fails
    real(c_double) :: dev_facs_boot(n_dev - 1), sigmas_boot(n_dev - 1)
    real(c_double) :: ln_means_boot(n_dev - 1, n_dev - 1), ln_sigmas_boot(n_dev - 1, n_dev - 1)
    real(c_double) :: ln_shifts_boot(n_dev - 1, n_dev - 1)
    integer(c_int) :: status

    character(len=999, kind=c_char) :: info

    n_resids = count(mask) - n_dev

    if (boot_type == PARAM) then
      call fit(triangle, dev_facs, sigmas, use_mask=.true., compute_resids=.false.)
    else if (boot_type == RESID) then
      call fit(triangle, dev_facs, sigmas, use_mask=.false., compute_resids=.true.)
    end if

    i_boot = 1
    main_loop: do while (i_boot <= n_boot)
      ! Simulate new parameters
      if (boot_type == RESID) call sample_resids("boot", ln_shifts_boot, ln_means_boot, ln_sigmas_boot)
      call boot_params(triangle, dev_facs_boot, sigmas_boot, ln_shifts_boot, ln_means_boot, ln_sigmas_boot, status)

      ! Simulate process error
      if (boot_type == RESID) call sample_resids("sim", ln_shifts_boot, ln_means_boot, ln_sigmas_boot)
      res = sim_res(triangle, dev_facs_boot, sigmas_boot, ln_shifts_boot, ln_means_boot, ln_sigmas_boot, status)

      reserve(i_boot) = res
      i_boot = i_boot + 1
    end do main_loop

    write(info, *) "Residual type: ", resids_type, c_new_line, "Triangle: "
    if (print) then
      call print_triangle(triangle, trim(info))
      call print_triangle(resids, "Resids: ")
    end if
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
        col_mask = obs_mask(1:n_rows, j + 1)
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

      resids = 0
      select case (resids_type)
       case (STANDARDISED)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            resids(i, j) = (indiv_dev_facs(i, j) - dev_facs(j)) * sqrt(triangle(i, j)) / (sigmas(j) * scale_facs(i, j))
          end do
        end do

       case (STUDENTISED)
        sigmas_jack = 0
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
              sigmas_jack(i, j) = sqrt(sum(triangle(1:n_rows, j) * &
                (indiv_dev_facs(1:n_rows, j) - dev_fac_jack) ** 2, mask=col_mask(1:n_rows)) / (n_pts_col - 1))
            else
              sigmas_jack(i, j) = extrapolate_sigma(sigmas_jack(i, :), j)
            end if
            resids(i, j) = (triangle(i, j + 1) - dev_facs(j) * triangle(i, j)) / &
              (sigmas_jack(i, j) * scale_facs(i, j) * sqrt(triangle(i, j)))
          end do
        end do

       case (LOGNORMAL)
        do j = 1, n_dev
          n_rows = n_dev - j
          do i = 1, n_rows
            ln_shifts(i, j) = dev_facs(j) * sqrt(triangle(i, j)) / sigmas(j)
            ln_sigmas(i, j) = sqrt(log(1 + 1 / ln_shifts(i, j) ** 2))
            ln_means(i, j) = log(ln_shifts(i, j)) - ln_sigmas(i, j) ** 2 / 2
            resids(i, j) = (triangle(i, j + 1) - dev_facs(j) * triangle(i, j)) / (sigmas(j) * sqrt(triangle(i, j)))
            resids(i, j) = (log(resids(i, j) + ln_shifts(i, j)) - ln_means(i, j)) / ln_sigmas(i, j)
          end do
        end do
      end select

      if (resids_type /= LOGNORMAL) then
        resids(1, n_dev - 1) = 0
        resids(:, n_dev - 2) = 0
      end if
      
      n_resids = (n_dev ** 2 - n_dev) / 2 - 3

      !   if (resids_type /= STUDENTISED) then
      !     resids_mean = 0
      !     do i = 1, n_dev - 1
      !       do j = 1, n_dev - i
      !         resids_mean = resids_mean + resids(i, j)
      !       end do
      !     end do
      !     resids_mean = resids_mean / n_resids
      !     do i = 1, n_dev - 1
      !       do j = 1, n_dev - i
      !         resids(i, j) = resids(i, j) - resids_mean
      !       end do
      !     end do
      !   end if
    end if

    deallocate(indiv_dev_facs)
  end subroutine fit

  ! Compute bootstrap development factors and dispersion parameters.
  subroutine boot_params(triangle, dev_facs_boot, sigmas_boot, ln_shifts_boot, &
    ln_means_boot, ln_sigmas_boot, status)
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(out) :: dev_facs_boot(:), sigmas_boot(:)
    real(c_double), intent(in) :: ln_means_boot(:, :), ln_sigmas_boot(:, :)
    real(c_double), intent(in) :: ln_shifts_boot(:, :)
    integer(c_int), intent(out) :: status

    real(c_double), allocatable :: triangle_boot(:, :)
    integer(c_int) :: i, j, n_rows, n_dev
    real(c_double) :: mean, sd, shape, scale

    real(c_double), allocatable :: resampled_pairs(:, :)
    integer(c_int), allocatable :: pair_indices(:), resampled_pair_indices(:)

    n_dev = size(triangle, 1)
    allocate(triangle_boot(n_dev, n_dev), source=0._c_double)

    status = SUCCESS
    select case (boot_type)
     case (PARAM)
      if(cond) then
        triangle_boot(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle(i, j - 1))
              triangle_boot(i, j) = rnorm_par(rng, i_thread, mean, sd)
              if (triangle_boot(i, j) <= 0) then
                status = FAILURE
                return
              end if

            else if (dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              triangle_boot(i, j) = rgamma_par(rng, i_thread, shape, scale)
              if (triangle_boot(i, j) <= 0) then
                status = FAILURE
                return
              end if
            end if
          end do
        end do

      else
        triangle_boot(:, 1) = triangle(:, 1)
        do j = 2, n_dev
          n_rows = n_dev + 1 - j
          do i = 1, n_rows
            if (dist == NORMAL) then
              mean = dev_facs(j - 1) * triangle_boot(i, j - 1)
              sd = sigmas(j - 1) * sqrt(triangle_boot(i, j - 1))
              triangle_boot(i, j) = rnorm_par(rng, i_thread, mean, sd)
              if (triangle_boot(i, j) <= 0) then
                status = FAILURE
                return
              end if

            else if (dist == GAMMA) then
              shape = (dev_facs(j - 1)**2 * triangle_boot(i, j - 1)) / sigmas(j - 1) **2
              scale = sigmas(j - 1) ** 2 / dev_facs(j - 1)
              triangle_boot(i, j) = rgamma_par(rng, i_thread, shape, scale)
              if (triangle_boot(i, j) <= 0) then
                status = FAILURE
                return
              end if
            end if
          end do
        end do
      end if

      call fit(triangle_boot, dev_facs_boot, sigmas_boot, use_mask=.false., compute_resids=.false.)

     case (RESID)
      if (cond) then
        triangle_boot(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == LOGNORMAL) then
              triangle_boot(i, j + 1) = exp(resampled_resids(i, j) * ln_sigmas(i, j) + &
                ln_means(i, j)) - ln_shifts(i, j)

              triangle_boot(i, j + 1) = dev_facs(j) * triangle(i, j) + &
                sigmas(j) * sqrt(triangle(i, j)) * triangle_boot(i, j + 1)

            else
              mean = dev_facs(j) * triangle(i, j)
              sd = sigmas(j) * sqrt(triangle(i, j))
              triangle_boot(i, j + 1) = mean + sd * resampled_resids(i, j)
            end if
            if (triangle_boot(i, j + 1) <= 0) then
              status = FAILURE
              return
            end if
          end do
        end do

      else
        triangle_boot(:, 1) = triangle(:, 1)
        do j = 1, n_dev - 1
          n_rows = n_dev - j
          do i = 1, n_rows
            if (resids_type == LOGNORMAL) then
              triangle_boot(i, j + 1) = exp(resampled_resids(i, j) * ln_sigmas_boot(i, j) + &
                ln_means_boot(i, j)) - ln_shifts_boot(i, j)
              triangle_boot(i, j + 1) = dev_facs(j) * triangle_boot(i, j) + &
                sigmas(j) * sqrt(triangle_boot(i, j)) * triangle_boot(i, j + 1)

            else
              mean = dev_facs(j) * triangle_boot(i, j)
              sd = sigmas(j) * sqrt(triangle_boot(i, j))
              triangle_boot(i, j + 1) = mean + sd * resampled_resids(i, j)
            end if
            if (triangle_boot(i, j + 1) <= 0) then
              status = FAILURE
              return
            end if
          end do
        end do
      end if

      call fit(triangle_boot, dev_facs_boot, sigmas_boot, use_mask=.false., compute_resids=.false.)

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
  end subroutine boot_params

  ! Simulate predictive distribution of reserve.
  function sim_res(triangle, dev_facs_boot, sigmas_boot, ln_shifts_sim, &
    ln_means_sim, ln_sigmas_sim, status) result(reserve)
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in) :: dev_facs_boot(:), sigmas_boot(:)
    real(c_double), intent(in) :: ln_means_sim(:, :), ln_sigmas_sim(:, :)
    real(c_double), intent(in) :: ln_shifts_sim(:, :)
    integer(c_int), intent(inout) :: status

    integer(c_int) :: n_dev
    integer(c_int) :: i, j, i_diag
    real(c_double) :: mean, sd, shape, scale
    real(c_double) :: reserve
    real(c_double), allocatable :: triangle_sim(:, :), latest(:)

    n_dev = size(triangle, 1)
    allocate(triangle_sim(n_dev, n_dev))

    status = SUCCESS
    triangle_sim = triangle
    if (boot_type == PARAM) then
      if (dist == NORMAL) then
        do i_diag = 1, n_dev - 1
          do i = i_diag + 1, n_dev
            j = n_dev + i_diag + 1 - i

            mean = dev_facs_boot(j - 1) * triangle_sim(i, j - 1)
            sd = sigmas_boot(j - 1) * sqrt(triangle_sim(i, j - 1))
            triangle_sim(i, j) = rnorm_par(rng, i_thread, mean, sd)
            if (triangle_sim(i, j) <= 0) then
              status = FAILURE
              return
            end if
          end do
        end do

      else if (dist == GAMMA) then
        do i_diag = 1, n_dev - 1
          do i = i_diag + 1, n_dev
            j = n_dev + i_diag + 1 - i
            shape = (dev_facs_boot(j - 1)**2 * triangle_sim(i, j - 1)) / sigmas_boot(j - 1) **2
            scale = sigmas_boot(j - 1) ** 2 / dev_facs_boot(j - 1)
            triangle_sim(i, j) = rgamma_par(rng, i_thread, shape, scale)
            if (triangle_sim(i, j) <= 0) then
              status = FAILURE
              return
            end if
          end do
        end do
      end if

    else if (boot_type == RESID) then
      do i_diag = 1, n_dev
        do i = i_diag + 1, n_dev
          j = n_dev + i_diag + 1 - i
          if (resids_type == LOGNORMAL) then
            triangle_sim(i, j) = exp(resampled_resids(i - 1, j - 1) * ln_sigmas_sim(i - 1, j - 1) + &
              ln_means_sim(i - 1, j - 1)) - ln_shifts_sim(i - 1, j - 1) ! triangle_sim temporarily holds value of epsilon

            triangle_sim(i, j) = dev_facs_boot(j - 1) * triangle_sim(i, j - 1) + &
              sigmas_boot(j - 1) * sqrt(triangle_sim(i, j - 1)) * triangle_sim(i, j)

          else
            mean = triangle_sim(i, j - 1) * dev_facs_boot(j - 1)
            sd = sigmas_boot(j - 1) * sqrt(triangle_sim(i, j - 1))
            triangle_sim(i, j) = mean + sd * resampled_resids(i - 1, j - 1)
          end if
          if (triangle_sim(i, j) <= 0 .or. isnan(triangle_sim(i, j))) then
            status = FAILURE
            return
          end if
        end do
      end do
    end if

    allocate(latest(n_dev))
    do j = 1, n_dev
      latest(j) = triangle_sim(n_dev + 1 - j, j)
    end do

    reserve = sum(triangle_sim(:, n_dev)) - sum(latest)
  end function sim_res

  subroutine sample_resids(which, ln_shifts_boot, ln_means_boot, ln_sigmas_boot)
    character(*, kind=c_char), intent(in) :: which

    integer(c_int) :: i, j, k, n_dev, n_resids, n_rows
    integer(c_int) :: rand_idxs(1, 2)
    real(c_double) :: ln_means_boot(:, :), ln_sigmas_boot(:, :)
    real(c_double) :: ln_shifts_boot(:, :)
    integer(c_int), allocatable :: idxs(:, :)

    n_dev = size(resids, 1) + 1
    n_resids = count(resids_mask)

    allocate(idxs(n_resids, 2), source=0)

    k = 1
    do j = 2, n_dev
      n_rows = n_dev + 1 - j
      do i = 1, n_rows
        if (resids_mask(i, j - 1)) then
          idxs(k, :) = [i, j - 1]
          k = k + 1
        end if
      end do
    end do

    if (which == "boot") then
      do j = 2, n_dev
        n_rows = n_dev + 1 - j
        do i = 1, n_rows
          rand_idxs(1, :) = idxs(1 + int(n_resids * runif_par(rng, i_thread)), :)
          resampled_resids(i, j - 1) = resids(rand_idxs(1, 1), rand_idxs(1, 2))
          ln_sigmas_boot(i, j - 1) =  ln_sigmas(rand_idxs(1, 1), rand_idxs(1, 2))
          ln_means_boot(i, j - 1) = ln_means(rand_idxs(1, 1), rand_idxs(1, 2))
          ln_shifts_boot(i, j - 1) = ln_shifts(rand_idxs(1, 1), rand_idxs(1, 2))
        end do
      end do

    else if (which == "sim") then
      do j = 2, n_dev
        n_rows = n_dev + 1 - j
        do i = n_rows + 1, n_dev
          rand_idxs(1, :) = idxs(1 + int(n_resids * runif_par(rng, i_thread)), :)
          resampled_resids(i - 1, j - 1) = resids(rand_idxs(1, 1), rand_idxs(1, 2))
          ln_sigmas_boot(i - 1, j - 1) =  ln_sigmas(rand_idxs(1, 1), rand_idxs(1, 2))
          ln_means_boot(i - 1, j - 1) = ln_means(rand_idxs(1, 1), rand_idxs(1, 2))
          ln_shifts_boot(i - 1, j - 1) = ln_shifts(rand_idxs(1, 1), rand_idxs(1, 2))
        end do
      end do
    end if

  end subroutine sample_resids

  pure real(c_double) function extrapolate_sigma(sigmas, col)
    real(c_double), intent(in) :: sigmas(:)
    integer(c_int), intent(in) :: col

    extrapolate_sigma = sqrt(min(sigmas(col - 1) ** 2, sigmas(col - 2) ** 2, sigmas(col - 1) ** 4 / sigmas(col - 2) ** 2))
  end function extrapolate_sigma

  subroutine boot_cpp(n_dev, triangle, boot_type_, opt, cond_, n_boot, reserve, seed) bind(c, name="mack_boot")
    integer(c_int), intent(in), value :: n_boot, n_dev, opt, boot_type_
    logical(c_bool), intent(in), value :: cond_
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: reserve(n_boot)
    integer(c_int), intent(in), value :: seed

    logical(c_bool) :: print_triangles

    print_triangles = .false.

    allocate(dev_facs(n_dev - 1), source=0._c_double)
    allocate(sigmas(n_dev - 1), source=0._c_double)
    allocate(scale_facs(n_dev, n_dev), source=0._c_double)
    allocate(sigmas_jack(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(ln_shifts(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(ln_means(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(ln_sigmas(n_dev - 1, n_dev - 1), source=0._c_double)
    allocate(resids(n_dev - 1, n_dev - 1), source=0._c_double)

    boot_type = boot_type_
    if (boot_type_ == RESID) then
      resids_type = opt
    else if (boot_type == PARAM) then
      dist = opt
    end if
    cond = cond_

    allocate(obs_mask(n_dev, n_dev))
    obs_mask = .true.

    if (first_call) then
      rng = init_rng(1, seed)
      first_call = .false.
    end if

    i_thread = 0

    call boot(n_dev, triangle, n_boot, reserve, obs_mask, print_triangles)

    deallocate(dev_facs)
    deallocate(sigmas)
    deallocate(scale_facs)
    deallocate(sigmas_jack)
    deallocate(ln_shifts)
    deallocate(ln_means)
    deallocate(ln_sigmas)
    deallocate(resids)
    deallocate(obs_mask)
  end subroutine boot_cpp

  function single_outlier(triangle, dev_facs, sigmas, outlier_rowidx, outlier_colidx, mean_factor, sd_factor, sim_dist)
    integer(c_int), intent(in):: outlier_rowidx, outlier_colidx
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in):: dev_facs(:), sigmas(:)
    integer(c_int), intent(in) :: sim_dist

    real(c_double) :: mean_factor, sd_factor
    real(c_double), allocatable:: single_outlier(:, :)
    real(c_double) :: shape, scale, mean, sd
    integer(c_int) :: n_dev, i, j
    logical :: negative

    if (mean_factor == 1 .and. sd_factor == 1) then
      single_outlier = triangle
      return
    end if

    n_dev = size(triangle, 1)
    allocate(single_outlier(n_dev, n_dev), source=0._c_double)
    single_outlier(:, 1) = triangle(:, 1)

    if (sim_dist == NORMAL) then
      negative = .true.
      do while (negative)
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

        mean = mean_factor * dev_facs(outlier_colidx - 1) * single_outlier(outlier_rowidx, outlier_colidx - 1)
        sd = sd_factor * sigmas(outlier_colidx - 1) * sqrt(single_outlier(outlier_rowidx, outlier_colidx - 1))
        single_outlier(outlier_rowidx, outlier_colidx) = rnorm_par(rng, i_thread, mean, sd)

        if (outlier_colidx < n_dev) then
          do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
            mean = dev_facs(j - 1) * single_outlier(outlier_rowidx, j - 1)
            sd = sigmas(j - 1) * sqrt(single_outlier(outlier_rowidx, j - 1))
            single_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)
          end do
        end if

        if (all(single_outlier >= 0)) negative = .false.
      end do

    else if (sim_dist == GAMMA) then
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

      shape = (mean_factor * dev_facs(outlier_colidx - 1))**2 * single_outlier(outlier_rowidx, outlier_colidx - 1) / &
        (sd_factor * sigmas(outlier_colidx - 1))**2
      scale = (sd_factor * sigmas(outlier_colidx - 1))**2 / (mean_factor * dev_facs(outlier_colidx - 1))
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

  function calendar_outlier(triangle, dev_facs, sigmas, outlier_diagidx, mean_factor, sd_factor, sim_dist)
    integer(c_int), intent(in) :: outlier_diagidx
    real(c_double), intent(in):: mean_factor, sd_factor
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in):: dev_facs(:), sigmas(:)
    integer(c_int), intent(in) :: sim_dist

    integer(c_int) :: i, j, n_dev, n_cols
    real(c_double), allocatable :: calendar_outlier(:, :)
    real(c_double) :: shape, scale, mean, sd
    logical :: negative

    if (mean_factor == 1 .and. sd_factor == 1) then
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

        if (sim_dist == NORMAL) then
          negative = .true.
          do while (negative)
            mean = mean_factor * dev_facs(n_cols - 1) * calendar_outlier(i, n_cols - 1)
            sd = sd_factor * sigmas(n_cols - 1) * sqrt(calendar_outlier(i, n_cols - 1))
            calendar_outlier(i, n_cols) = rnorm_par(rng, i_thread, mean, sd)
            do j = n_cols + 1, n_dev + 1 - i
              mean = dev_facs(j - 1) * calendar_outlier(i, j - 1)
              sd = sigmas(j - 1) * sqrt(calendar_outlier(i, j - 1))
              calendar_outlier(i, j) = rnorm_par(rng, i_thread, mean, sd)
            end do
            
            if (all(calendar_outlier >= 0)) negative = .false.
          end do

        else if (sim_dist == GAMMA) then
          shape = (mean_factor * dev_facs(n_cols - 1))**2 * calendar_outlier(i, n_cols - 1) / (sd_factor * sigmas(n_cols - 1))**2
          scale = (sd_factor * sigmas(n_cols - 1))**2 / (mean_factor * dev_facs(n_cols - 1))
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

  function origin_outlier(triangle, dev_facs, sigmas, outlier_rowidx, mean_factor, sd_factor, sim_dist)
    integer(c_int), intent(in):: outlier_rowidx
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(in) :: mean_factor, sd_factor
    real(c_double), intent(in):: dev_facs(:), sigmas(:)
    integer(c_int), intent(in) :: sim_dist

    real(c_double) :: shape, scale, mean, sd
    real(c_double), allocatable:: origin_outlier(:, :)
    integer(c_int) :: n_dev, j
    logical :: negative

    if (mean_factor == 1 .and. sd_factor == 1) then
      origin_outlier = triangle
      return
    end if

    n_dev = size(triangle, dim=1)
    origin_outlier = triangle

    do j = 2, n_dev + 1 - outlier_rowidx
      if (sim_dist == NORMAL) then
        negative = .true.
        do while (negative)
        mean = mean_factor * dev_facs(j - 1) * origin_outlier(outlier_rowidx, j - 1)
        sd = sd_factor * sigmas(j - 1) * sqrt(origin_outlier(outlier_rowidx, j - 1))
        origin_outlier(outlier_rowidx, j) = rnorm_par(rng, i_thread, mean, sd)

        if (all(origin_outlier >= 0)) negative = .false.
        end do

      else if (sim_dist == GAMMA) then
        shape = (mean_factor * dev_facs(j - 1))**2 * origin_outlier(outlier_rowidx, j - 1) / (sd_factor * sigmas(j - 1))**2
        scale = (sd_factor * sigmas(j - 1))**2 / (mean_factor * dev_facs(j - 1))
        origin_outlier(outlier_rowidx, j) = rgamma_par(rng, i_thread, shape, scale)
      end if
    end do
  end function origin_outlier

end module mod_mack
