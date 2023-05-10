module mod_helpers

  use, intrinsic :: iso_c_binding
  use, intrinsic :: omp_lib
  use mod_global
  use mod_interface

  implicit none

contains

  function init_omp() result(n_threads)
    integer(c_int) :: n_threads
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
    ! call omp_set_num_threads(1)
    ! n_threads = 1
  end function init_omp

  

  function validate_rng(n_samples) result(res) bind(c, name="validate_rng_")
    integer(c_int), intent(in), value :: n_samples
    integer(c_int) :: res
    integer(c_int) :: i, j
    real(c_double), allocatable :: samples(:, :)
    real(c_double), allocatable :: flattened(:)
    integer(c_int) :: n_threads, total
    type(c_ptr) :: rng

    n_threads = init_omp()
    rng = init_rng(n_threads, 42)
    allocate(samples(n_samples, n_threads))

    !$omp parallel do num_threads(n_threads) firstprivate(n_threads, n_samples) shared(samples, rng) schedule(static)
    do i = 1, n_threads
      do j = 1, n_samples
        samples(j, i) = runif_par(rng, i)
      end do
    end do
    !$omp end parallel do

    allocate(flattened(size(samples, 1) * size(samples, 2)))
    flattened = pack(samples, mask=.true.)

    res = SUCCESS
    total = size(flattened)
    do j = 1, total
      do i = j + 1, total
        if (flattened(j) == flattened(i)) then
          res = FAILURE
        endif
      end do
    end do

    if (res == SUCCESS) then
      print *, "Succes"
    else
      print *, "Failure"
    end if

  end function validate_rng

  subroutine print_triangle(triangle)
    real(c_double), intent(in) :: triangle(:, :)
    integer :: i, n_dev
    character(len=999, kind=c_char) :: str, str_

    n_dev = size(triangle, 1)

    write(str, '(*(f9.2))') triangle(1, :)
    str = c_new_line // trim(str) // c_new_line // c_carriage_return
    do i = 2, n_dev
      write(str_, '(*(f9.2))') triangle(i, :)
      str = trim(str) // trim(str_) // c_new_line // c_carriage_return
    end do

    call rprint_par(str)
    call r_flush_console()
  end subroutine print_triangle

  subroutine print_vector(vector)
    real(c_double), intent(in) :: vector(:)
    character(len=999, kind=c_char) :: str

    write(str, '(*(f9.2))') vector
    str = c_new_line // trim(str) // c_new_line
    call rprint_par(str)
    call r_flush_console()
  end subroutine print_vector

end module mod_helpers
