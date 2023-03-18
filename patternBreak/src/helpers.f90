module helpers

  use iso_c_binding
  use global
  use interface
  use omp_lib

  implicit none

contains

  subroutine single_outlier_glm(triangle, outlier_rowidx, outlier_colidx, factor, betas, rng, i_thread)
    integer(c_int), intent(in) :: outlier_rowidx
    integer(c_int), intent(in) :: outlier_colidx
    real(c_double), intent(in) :: factor
    real(c_double), intent(in) :: betas(:)
    real(c_double), intent(inout) :: triangle(:, :)
    integer(c_int), intent(in) :: i_thread
    type(c_ptr), intent(in) :: rng

    integer(c_int) :: i, j

    real(c_double) :: lambda

    i = outlier_rowidx
    j = outlier_colidx

    if (j == 1 .and. i == 1) then
      lambda = factor * exp(betas(1))
    else if (j == 1) then
      lambda = factor * exp(betas(1) + betas(i - 1))
    else if (i == 1) then
      lambda = factor * exp(betas(1) + betas(j - 1))
    else
      lambda = factor * exp(betas(1) + betas(i - 1) + betas(j - 1))
    end if

    triangle(i, j) = real(rpois_par(rng, i_thread, lambda), kind=c_double)

  end subroutine single_outlier_glm

  subroutine calendar_outlier_glm(triangle, outlier_diagidx, factor, betas, rng, i_thread)
    integer(c_int), intent(in) :: outlier_diagidx
    real(c_double), intent(in) :: factor
    real(c_double), intent(in) :: betas(:)
    real(c_double), intent(inout) :: triangle(:, :)
    integer(c_int), intent(in) :: i_thread
    type(c_ptr), intent(in) :: rng

    real(c_double) :: lambda
    integer(c_int) :: i, j, n_dev

    n_dev = size(triangle, dim=1)

    do i = 1, n_dev
      j = n_dev + 2 - outlier_diagidx - i
      if (j < 1) cycle

      if (j == 1 .and. i == 1) then
        lambda = factor * exp(betas(1))
      else if (j == 1) then
        lambda = factor * exp(betas(1) + betas(i - 1))
      else if (i == 1) then
        lambda = factor * exp(betas(1) + betas(j - 1))
      else
        lambda = factor * exp(betas(1) + betas(i - 1) + betas(j - 1))
      end if

      triangle(i, j) = rpois_par(rng, i_thread, lambda)
    end do

  end subroutine calendar_outlier_glm

  subroutine origin_outlier_glm(triangle, outlier_rowidx, factor, betas, rng, i_thread)
    integer(c_int), intent(in) :: outlier_rowidx
    real(c_double), intent(in) :: factor
    real(c_double), intent(in) :: betas(:)
    real(c_double), intent(inout) :: triangle(:, :)
    integer(c_int), intent(in) :: i_thread
    type(c_ptr), intent(in) :: rng

    real(c_double) :: lambda
    integer(c_int) :: i, j, n_dev

    n_dev = size(triangle, dim=1)

    i = outlier_rowidx

    do j = 1, n_dev + 1 - i

      if (j == 1 .and. i == 1) then
        lambda = factor * exp(betas(1))
      else if (j == 1) then
        lambda = factor * exp(betas(1) + betas(i - 1))
      else if (i == 1) then
        lambda = factor * exp(betas(1) + betas(j - 1))
      else
        lambda = factor * exp(betas(1) + betas(i - 1) + betas(j - 1))
      end if

      triangle(i, j) = rpois_par(rng, i_thread, lambda)

    end do

  end subroutine origin_outlier_glm

  function init_omp() result(n_threads)
    integer(c_size_t) :: n_threads
    character(len=999, kind = c_char) :: OMP_NUM_THREADS
    integer(c_int) :: is_set

    call get_environment_variable("OMP_NUM_THREADS", OMP_NUM_THREADS, status=is_set)
    if (is_set == 0) then
      read(OMP_NUM_THREADS, *) n_threads
    else
      n_threads = nint(0.5 * omp_get_num_procs())
      n_threads = min(omp_get_thread_limit(), n_threads)
      n_threads = min(omp_get_max_threads(), n_threads)
    end if
  end function init_omp

    ! subroutine progress_bar(counter, max, inc)

    !    integer(c_int), intent(in) :: counter
    !    integer(c_int), intent(in) :: max
    !    integer(c_int), intent(in) :: inc

    !    character(len=999) :: buf
    !    character(kind=c_char, len=999) :: progress_str
    !    integer :: cols
    !    real(c_double) :: progress
    !    real(c_double) :: pct_progress

    !    if (mod(counter, inc) == 0) then

    !       progress = real(counter, kind=c_double) / real(max, kind=c_double)
    !       pct_progress = 100*progress

    !       call get_environment_variable("COLUMNS", buf)
    !       read(buf, *) cols

    !       write(progress_str, "('Progress: ', f6.2)") pct_progress
    !       call Rprintf(achar(13) // c_null_char)
    !       call Rprintf(repeat(" ", cols) // c_null_char)
    !       call Rprintf(achar(13) // c_null_char)
    !       call Rprintf(trim(progress_str) // c_null_char )
    !       call Rprintf(achar(13) // c_null_char)
    !       call R_FlushConsole()

    !    end if

    !    if (counter == max) then

    !       call get_environment_variable("COLUMNS", buf)
    !       read(buf, *) cols

    !       call Rprintf(repeat(" ", cols) // c_null_char)

    !    end if

    ! end subroutine progress_bar

  end module helpers
