module mod_glm

  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use mod_global
  use mod_helpers
  use mod_interface

  implicit none

  integer(c_int) :: boot_type, dist
  logical(c_bool), allocatable :: mask(:, :)
  real(c_double), allocatable :: resids(:, :), triangle_fit(:, :)

  !$omp threadprivate(mask, boot_type, dist, resids, triangle_fit)

contains

  subroutine sim(n_dev, triangle, sim_type, boot_type_, n_conf, m_conf, conf, &
    n_boot, res, show_progress, seed) bind(c, name="glm_sim")
    integer(c_int), intent(in), value :: n_dev, n_boot, n_conf, m_conf
    integer(c_int), intent(in), value :: sim_type, boot_type_
    real(c_double), intent(in) :: conf(n_conf, m_conf)
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(out) :: res(n_boot * n_conf, m_conf + 1)
    logical(c_bool), value, intent(in) :: show_progress
    integer(c_int), value, intent(in) :: seed

    integer(c_int) :: i, j, k, i_sim
    integer(c_int) :: outlier_rowidx, outlier_colidx, outlier_diagidx
    integer(c_int) :: excl_diagidx, excl_rowidx, excl_colidx

    real(c_double), allocatable :: reserve(:)
    real(c_double) :: triangle_sim(n_dev, n_dev)
    real(c_double) :: factor

    integer(c_int) :: n_covs
    real(c_double), allocatable :: betas(:)
    real(c_double) :: disp

    integer(c_int) :: status

    type(c_ptr) ::pgbar

    n_covs = 2*n_dev - 1

    allocate(betas(n_covs))
    call fit(triangle, betas, disp, use_mask = .false., comp = .false., status=status)

    if (first_call) then
      n_threads = init_omp()
      rng = init_rng(n_threads, seed)
      first_call = .false._c_bool
    end if

    pgbar = pgbar_create(n_conf, 1)

    boot_type = boot_type_
    res = 0

    !$omp parallel num_threads(n_threads) shared(conf, res) &
    !$omp& copyin(rng, boot_type, dist, mask) &
    !$omp& private(reserve, outlier_rowidx, outlier_colidx, excl_rowidx, excl_colidx, factor, triangle_sim, i_sim) &
    !$omp& firstprivate(n_boot, n_dev, betas, n_conf, m_conf, triangle, seed)

    allocate(reserve(n_boot), source = 0._c_double)
    allocate(mask(n_dev, n_dev))
    allocate(resids(n_dev, n_dev), source=0._c_double)
    allocate(triangle_fit(n_dev, n_dev), source=0._c_double)
    i_thread = omp_get_thread_num()

    !$omp do schedule(dynamic, 25)
    do i_sim = 1, n_conf
      if (show_progress) call pgbar_incr(pgbar)
      if (sim_type == SINGLE) then
        outlier_rowidx = int(conf(i_sim, 1))
        outlier_colidx = int(conf(i_sim, 2))
        factor = conf(i_sim, 3)
        excl_rowidx = int(conf(i_sim, 4))
        excl_colidx = int(conf(i_sim, 5))
        if (boot_type == PARAM) dist = int(conf(i_sim, 6))

        mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            mask(i, j) = .false.
          end do
        end do
        if (boot_type == RESID) then
          mask(1, n_dev) = .false.
          mask(n_dev, 1) = .false.
        end if
        mask(excl_rowidx, excl_colidx) = .false.

        triangle_sim = triangle
        call single_outlier(triangle_sim, outlier_rowidx, outlier_colidx, factor, betas)
        call boot(n_dev, triangle_sim, n_boot, reserve)

        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_conf) = transpose(spread(conf(i_sim, :), 2, n_boot))
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_conf + 1) = reserve

      else if (sim_type == CALENDAR) then
        outlier_diagidx = int(conf(i_sim, 1))
        factor = conf(i_sim, 2)
        excl_diagidx = int(conf(i_sim, 3))
        if (boot_type == PARAM) dist = int(conf(i_sim, 4))

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
        if (boot_type == PARAM) then
          mask(1, n_dev) = .true.
          mask(n_dev, 1) = .true.
        end if

        triangle_sim = triangle

        call calendar_outlier(triangle_sim, outlier_diagidx, factor, betas)
        call boot(n_dev, triangle_sim, n_boot, reserve)

        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_conf) = transpose(spread(conf(i_sim, :), 2, n_boot))
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_conf + 1) = reserve

      else if (sim_type == ORIGIN) then
        outlier_rowidx = int(conf(i_sim, 1))
        factor = conf(i_sim, 2)
        excl_rowidx = int(conf(i_sim, 3))
        if (boot_type == PARAM) dist = int(conf(i_sim, 4))

        mask = .true.
        do j = 1, n_dev
          do i = n_dev + 2 - j, n_dev
            mask(i, j) = .false.
          end do
        end do
        if (boot_type == PARAM) then
          mask(excl_rowidx, 2:) = .false. ! can't remove first observation
          mask(1, n_dev) = .true.
          mask(n_dev, 1) = .true.
        else
          mask(excl_rowidx, :) = .false.
        end if

        triangle_sim = triangle

        call origin_outlier(triangle_sim, outlier_rowidx, factor, betas)
        call boot(n_dev, triangle_sim, n_boot, reserve)

        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), 1:m_conf) = transpose(spread(conf(i_sim, :), 2, n_boot))
        res(((i_sim - 1)*n_boot + 1):(i_sim*n_boot), m_conf + 1) = reserve
      end if

    end do
    !$omp end do
    deallocate(reserve)
    deallocate(mask)
    deallocate(resids)
    deallocate(triangle_fit)
    !$omp end parallel
  end subroutine sim

  subroutine boot(n_dev, triangle, n_boot, reserve)
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    real(c_double), intent(inout) :: reserve(n_boot)
    integer(c_int), intent(in) :: n_dev, n_boot

    integer(c_int) :: n_pts, n_covs, n_pred, n_obs
    real(c_double), allocatable :: betas(:), X_pred(:, :), y_pred(:)
    real(c_double), allocatable :: triangle_boot(:, :), betas_boot(:), resids_boot(:, :)
    real(c_double), allocatable :: triangle_pred(:, :)

    real(c_double) :: disp, disp_boot
    real(c_double) :: lambda, mean, sd, shape, scale

    real(c_double) :: resid_sim

    real(c_double), allocatable :: flat_resids(:)

    integer(c_int) :: i, j, k, i_boot
    integer(c_int) :: n_excl_resids, n_resids
    integer(c_int) :: status

    n_pts = count(mask)
    n_obs = (n_dev ** 2 + n_dev) / 2
    n_covs = 2*n_dev - 1
    n_pred = n_dev ** 2 - n_obs

    allocate(betas(n_covs), source = 0._c_double)
    allocate(flat_resids(n_pts), source = 0._c_double)

    allocate(resids_boot(n_pts, n_pts), source = 0._c_double)
    allocate(triangle_boot(n_dev, n_dev), source = 0._c_double)
    allocate(triangle_pred(n_dev, n_dev), source = 0._c_double)

    allocate(X_pred(n_pred, n_covs), source = 0._c_double)
    allocate(y_pred(n_pred), source = 0._c_double)
    allocate(betas_boot(n_covs), source = 0._c_double)

    if (boot_type == RESID) then
      call fit(triangle, betas, disp, use_mask = .false., comp = .true., status=status)
      flat_resids = pack(resids, mask)
    else
      call fit(triangle, betas, disp, use_mask = .true., comp = .true., status=status)
    end if

    i_boot = 1
    main_loop: do while (i_boot <= n_boot)
      if (boot_type == RESID) then
        triangle_boot = 0
        resids_boot = 0
        do i = 1, n_dev
          do j = 1, n_dev + 1 - i
            resids_boot(i, j) = flat_resids(1 + int(n_pts * runif_par(rng, i_thread)))
            triangle_boot(i, j) = resids_boot(i, j) * sqrt(triangle_fit(i, j)) + triangle_fit(i, j)
            if (triangle_boot(i, j) <= 0) cycle main_loop
          end do
        end do

      else
        triangle_boot = 0
        do i = 1, n_dev
          do j = 1, n_dev + 1 - i
            if (dist == NORMAL) then
              mean = triangle_fit(i, j)
              sd = sqrt(disp * triangle_fit(i, j))
              triangle_boot(i, j) = rnorm_par(rng, i_thread, mean, sd)
              if (triangle_boot(i, j) < 0) cycle main_loop

            else if (dist == GAMMA) then
              shape = triangle_fit(i, j)**2 / (disp * triangle_fit(i, j))
              scale = (disp * triangle_fit(i, j)) / triangle_fit(i, j)
              triangle_boot(i, j) = rgamma_par(rng, i_thread, shape, scale)

            else if (dist == POISSON) then
              lambda = triangle_fit(i, j) / disp
              triangle_boot(i, j) = disp * rpois_par(rng, i_thread, lambda)
            end if
          end do
        end do

      end if

      call fit(triangle_boot, betas_boot, disp_boot, use_mask = .false., comp = .false., status=status)
      if (status == FAILURE) cycle main_loop

      X_pred = 0.0
      X_pred(:, 1) = 1.0

      k = 1
      do i = 2, n_dev
        do j = n_dev + 1 - i + 1, n_dev
          X_pred(k, i) = 1
          if (j /= 1) X_pred(k, n_dev + j - 1) = 1
          k = k + 1
        end do
      end do

      y_pred = exp(matmul(X_pred, betas_boot))

      k = 1
      do i = 2, n_dev
        do j = n_dev + 1 - i + 1, n_dev
          triangle_pred(i, j) = y_pred(k)
          k = k + 1
        end do
      end do

      if (boot_type == RESID) then
        k = 1
        do i = 2, n_dev
          do j = n_dev + 1 - i + 1, n_dev
            resid_sim = flat_resids(1 + int(n_pts * runif_par(rng, i_thread)))
            triangle_pred(i, j) = triangle_pred(i, j) + resid_sim * sqrt(triangle_pred(i, j))
          end do
        end do
      else
        do i = 1, n_pred
          lambda = y_pred(i)
          if (lambda > 1e7) cycle main_loop
          y_pred(i) = rpois_par(rng, i_thread, lambda)
        end do
      end if

      reserve(i_boot) = sum(y_pred)
      i_boot = i_boot + 1
    end do main_loop
  end subroutine boot

  subroutine boot_cpp(n_dev, triangle, boot_type_, opt, n_boot, reserve, seed) bind(c, name="glm_boot")
    integer(c_int), intent(in), value :: n_boot, n_dev
    real(c_double), intent(in) :: triangle(n_dev, n_dev)
    integer(c_int), intent(in), value :: boot_type_, opt
    real(c_double), intent(out) :: reserve(n_boot)
    integer(c_int), intent(in), value :: seed

    integer(c_int) :: i, j

    boot_type = boot_type_
    if (boot_type_ == PARAM) dist = opt

    allocate(mask(n_dev, n_dev))
    allocate(resids(n_dev, n_dev), source=0._c_double)
    allocate(triangle_fit(n_dev, n_dev), source=0._c_double)

    mask = .true.
    do j = 1, n_dev
      do i = n_dev + 2 - j, n_dev
        mask(i, j) = .false.
      end do
    end do

    if (first_call) then
      rng = init_rng(1, seed)
      first_call = .false.
    end if

    i_thread = 0

    call boot(n_dev, triangle, n_boot, reserve)

    deallocate(mask)
    deallocate(resids)
    deallocate(triangle_fit)
  end subroutine boot_cpp

  subroutine fit(triangle, betas, disp, use_mask, comp, status)
    real(c_double), intent(in) :: triangle(:, :)
    real(c_double), intent(out) :: betas(:)
    real(c_double), intent(out) :: disp
    integer(c_int), intent(out) :: status
    logical, intent(in) :: use_mask, comp
    real(c_double), allocatable :: X(:, :), y(:)
    real(c_double), allocatable :: y_fit(:), X_fit(:, :), y_model(:)

    integer(c_int) :: n_obs, n_covs, n_pts, n_dev, n_fit
    integer(c_int) :: i, j, k, l

    real(c_double) :: diff, eps
    real(c_double), allocatable :: IPIV(:, :), rhs(:), lhs(:, :), mu(:)
    real(c_double), allocatable :: W(:, :), z(:), eta(:),  betas_old(:)
    logical(c_bool), allocatable :: mask_loc(:, :)
    integer(c_int) :: info

    integer(c_int) :: lwork
    real(c_double), allocatable :: work(:, :)

    real(c_double) :: temp

    integer(c_int), parameter :: max_iter = 1e3
    integer(c_int) :: n_iter

    n_dev = size(triangle, dim=1)

    allocate(mask_loc(n_dev, n_dev), source = .true._c_bool)
    if (use_mask) then
      mask_loc = mask
    else
      do j = 1, n_dev
        do i = n_dev + 2 - j, n_dev
          mask_loc(i, j) = .false.
        end do
      end do
    end if

    ! Compute GLM matrix dimensions.
    n_fit = (n_dev**2 + n_dev) / 2
    n_obs = n_fit
    if (use_mask) n_obs = count(mask)
    n_covs = 2*n_dev - 1

    ! Allocate IRWLS variables.
    allocate(X(n_obs, n_covs), source = 0._c_double)
    allocate(y(n_obs), source = 0._c_double)
    allocate(y_model(n_obs), source = 0._c_double)
    allocate(W(n_obs, n_obs), source = 0._c_double)
    allocate(z(n_obs), source = 0._c_double)
    allocate(eta(n_obs), source = 0._c_double)
    allocate(mu(n_obs), source = 0._c_double)
    allocate(X_fit(n_fit, n_covs), source = 0._c_double)
    allocate(y_fit(n_fit), source = 0._c_double)

    ! Allocate LAPACK helper variables.
    allocate(IPIV(n_covs, n_covs), source = 0._c_double)
    allocate(rhs(n_covs), source = 0._c_double)
    allocate(lhs(n_covs, n_covs), source = 0._c_double)
    allocate(betas_old(n_covs), source = 0._c_double)

    lwork = max(1, n_covs * n_obs + max(n_covs * n_obs, 1))
    allocate(work(lwork, lwork), source = 0._c_double)

    ! Set up feature matrix and response vector.
    X = 0
    X(:, 1) = 1._c_double

    X_fit = 0
    X_fit(:, 1) = 1._c_double

    l = 1
    k = 1
    do i = 1, n_dev
      do j = 1, n_dev + 1 - i
        if (mask_loc(i, j)) then
          if (i /= 1) X(k, i) = 1
          if (j /= 1) X(k, n_dev + j - 1) = 1
          y(k) = triangle(i, j)
          k = k + 1
        end if
        if (i /= 1) X_fit(l, i) = 1
        if (j /= 1) X_fit(l, n_dev + j - 1) = 1
        l = l + 1
      end do
    end do

    ! Initialise GLM coefficients.
    betas = spread(0._c_double, 1, n_covs)
    mu = y + 0.1
    eta = log(mu)

    ! Fit Poisson GLM using IRWLS.
    IPIV = 0
    info = 0
    diff = 1E6
    eps = 1E-5
    n_iter = 0
    status = SUCCESS
    do while (diff > eps .and. n_iter < max_iter)
      betas_old = betas

      W = 0
      do i = 1, n_obs
        W(i, i) = sqrt(mu(i))
      end do

      z = eta + y / mu - 1
      lhs = matmul(W, X)
      rhs = matmul(W, z)

      if (any(isnan(rhs)) .or. any(isnan(lhs))) then
        status = FAILURE
        return
      end if

      call dgels('N', n_obs, n_covs, 1, lhs, n_obs, rhs, n_obs, work, lwork, info)

      betas = rhs(1:n_covs)
      diff = norm2(betas - betas_old)

      eta = matmul(X, betas)
      mu = exp(eta)

      if (info /= 0) print *, raise(2)

      if (any(mu <= 0)) then
        status = FAILURE
        return
      end if
      if (any(isnan(mu)) .or. .not. all(ieee_is_finite(mu))) then
        status = FAILURE
        return
      end if

      n_iter = n_iter + 1
    end do

    if (n_iter == max_iter) then
      status = FAILURE
      return
    end if

    y_fit = exp(matmul(X_fit, betas))
    y_model = exp(matmul(X, betas))

    disp = sum((y - y_model)**2 / y_model)  / (n_obs - n_covs)

    if (comp) then
      ! Compute fitted triangle.
      k = 1
      do i = 1, n_dev
        do j = 1, n_dev + 1 - i
          triangle_fit(i, j) = y_fit(k)
          k = k + 1
        end do
      end do

      ! Compute residuals.
      do i = 1, n_dev
        do j = 1, n_dev + 1 - i
          resids(i, j) = (triangle(i, j) - triangle_fit(i, j)) / sqrt(triangle_fit(i, j))
        end do
      end do
    end if
  end subroutine fit

  subroutine single_outlier(triangle, outlier_rowidx, outlier_colidx, factor, betas)
    integer(c_int), intent(in) :: outlier_rowidx
    integer(c_int), intent(in) :: outlier_colidx
    real(c_double), intent(in) :: factor
    real(c_double), intent(in) :: betas(:)
    real(c_double), intent(inout) :: triangle(:, :)

    integer(c_int) :: i, j, n_dev

    real(c_double) :: lambda

    real(c_double) :: temp

    n_dev = size(triangle, 1)

    i = outlier_rowidx
    j = outlier_colidx

    if (j == 1 .and. i == 1) then
      lambda = factor * exp(betas(1))
    else if (j == 1) then
      lambda = factor * exp(betas(1) + betas(i))
    else if (i == 1) then
      lambda = factor * exp(betas(1) + betas(n_dev + j - 1))
    else
      lambda = factor * exp(betas(1) + betas(i) + betas(n_dev + j - 1))
    end if

    triangle(i, j) = rpois_par(rng, i_thread, lambda)
  end subroutine single_outlier

  subroutine calendar_outlier(triangle, outlier_diagidx, factor, betas)
    integer(c_int), intent(in) :: outlier_diagidx
    real(c_double), intent(in) :: factor
    real(c_double), intent(in) :: betas(:)
    real(c_double), intent(inout) :: triangle(:, :)


    real(c_double) :: lambda
    integer(c_int) :: i, j, n_dev

    n_dev = size(triangle, dim=1)

    do i = 1, n_dev
      j = n_dev + 2 - outlier_diagidx - i
      if (j < 1) cycle

      if (j == 1 .and. i == 1) then
        lambda = factor * exp(betas(1))
      else if (j == 1) then
        lambda = factor * exp(betas(1) + betas(i))
      else if (i == 1) then
        lambda = factor * exp(betas(1) + betas(n_dev + j - 1))
      else
        lambda = factor * exp(betas(1) + betas(i) + betas(n_dev + j - 1))
      end if

      triangle(i, j) = rpois_par(rng, i_thread, lambda)
    end do
  end subroutine calendar_outlier

  subroutine origin_outlier(triangle, outlier_rowidx, factor, betas)
    integer(c_int), intent(in) :: outlier_rowidx
    real(c_double), intent(in) :: factor
    real(c_double), intent(in) :: betas(:)
    real(c_double), intent(inout) :: triangle(:, :)

    real(c_double) :: lambda
    integer(c_int) :: i, j, n_dev

    n_dev = size(triangle, dim=1)

    i = outlier_rowidx

    do j = 1, n_dev + 1 - i
      if (j == 1 .and. i == 1) then
        lambda = factor * exp(betas(1))
      else if (j == 1) then
        lambda = factor * exp(betas(1) + betas(i))
      else if (i == 1) then
        lambda = factor * exp(betas(1) + betas(n_dev + j - 1))
      else
        lambda = factor * exp(betas(1) + betas(i) + betas(n_dev + j - 1))
      end if

      triangle(i, j) = rpois_par(rng, i_thread, lambda)
    end do

  end subroutine origin_outlier

end module mod_glm

! lhs = matmul(matmul(transpose(X), W), X)
! rhs = matmul(matmul(transpose(X), W), z)
! call dgesv(n_covs, 1, lhs, n_covs, IPIV, rhs, n_covs, info)
