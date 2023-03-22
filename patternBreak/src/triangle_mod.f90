module triangle_mod
  use, intrinsic :: iso_c_binding
  use global_mod
  use helpers_mod
  implicit none

  type :: triangle
    integer(c_int), allocatable :: n_dev
    real(c_double), allocatable :: data(:, :)
    real(c_double), allocatable :: indiv_dev_facs(:, :)
    real(c_double), allocatable :: dev_facs(:), sigmas(:)
    real(c_double), allocatable :: resids(:, :)
    real(c_double), allocatable :: scale_facs(:, :), sigmas_jack(:, :)
    real(c_double), allocatable :: log_normal_sigmas(:), log_normal_means(:)
    real(c_double), allocatable :: log_normal_shift(:)
    logical(c_bool), allocatable :: mask(:, :)
    integer(c_int), allocatable :: n_resids
    integer(c_int), allocatable :: resids_type

    contains
      procedure :: fit => fit_triangle
      procedure :: print => print_triangle
      final :: destroy_triangle
  end type triangle

  interface triangle
    procedure :: init_triangle
  end interface

contains

function init_triangle(n_dev, data) result(this)
  integer(c_int), intent(in) :: n_dev
  real(c_double), intent(in) :: data(n_dev, n_dev)
  type(triangle) :: this

  allocate(this%indiv_dev_facs(n_dev, n_dev), source=0._c_double)
  allocate(this%dev_facs(n_dev - 1), source=0._c_double)
  allocate(this%sigmas(n_dev - 1), source=0._c_double)
  allocate(this%scale_facs(n_dev - 1, n_dev - 1), source=0._c_double)
  allocate(this%sigmas_jack(n_dev - 1, n_dev - 1), source=0._c_double)
  allocate(this%log_normal_shift(n_dev - 1), source=0._c_double)
  allocate(this%log_normal_means(n_dev - 1), source=0._c_double)
  allocate(this%log_normal_sigmas(n_dev - 1), source=0._c_double)
  allocate(this%resids(n_dev - 1, n_dev - 1), source=0._c_double)
  allocate(this%mask(n_dev, n_dev), source=.true._c_bool)
  
  this%n_dev = n_dev
  this%data = data
  call this%fit(use_mask = .false., compute_resids = .false.)
  this%n_resids = (n_dev ** 2 - n_dev) / 2
end function init_triangle

subroutine destroy_triangle(this)
  type(triangle) :: this
  deallocate(this%data)
  deallocate(this%indiv_dev_facs)
  deallocate(this%log_normal_shift)
  deallocate(this%log_normal_means)
  deallocate(this%log_normal_sigmas)
  deallocate(this%scale_facs)
  deallocate(this%sigmas_jack)
  deallocate(this%dev_facs)
  deallocate(this%sigmas)
  deallocate(this%resids)
  deallocate(this%mask)
end subroutine destroy_triangle

subroutine print_triangle(this)
  class(triangle) :: this
  integer :: i, j

  do i = 1, this%n_dev
    print "(999(rd, f9.2))", this%data(i, :)
  end do
end subroutine print_triangle

subroutine fit_triangle(this, use_mask, compute_resids)
  class(triangle), intent(inout) :: this
  logical, intent(in) :: compute_resids, use_mask

  integer(c_int) :: i, j, n_rows, n_pts_col, n_resids
  logical(c_bool), allocatable :: col_mask(:)
  real(c_double) :: resids_mean, dev_fac_jack

  if (use_mask) then
    do j = 1, this%n_dev - 1
      n_rows = this%n_dev - j
      do i = 1, n_rows
        this%indiv_dev_facs(i, j) = this%data(i, j + 1) / this%data(i, j)
      end do
      col_mask = this%mask(1:n_rows, j + 1)
      n_pts_col = count(col_mask)
      this%dev_facs(j) = sum(this%data(1:n_rows, j + 1), mask=col_mask) / sum(this%data(1:n_rows, j), mask=col_mask)
      if (n_pts_col >= 2) then
        this%sigmas(j) = sqrt(sum(this%data(1:n_rows, j) * (this%indiv_dev_facs(1:n_rows, j) - &
          this%dev_facs(j)) ** 2, mask=col_mask) / n_pts_col)
      else
        this%sigmas(j) = extrapolate_sigma(this%sigmas, j)
      end if
    end do

  else
    do j = 1, this%n_dev - 1
      n_rows = this%n_dev - j
      do i = 1, n_rows
        this%indiv_dev_facs(i, j) = this%data(i, j + 1) / this%data(i, j)
      end do
      this%dev_facs(j) = sum(this%data(1:n_rows, j + 1)) / sum(this%data(1:n_rows, j))
      if (j < this%n_dev - 1) then
        this%sigmas(j) = sqrt(sum(this%data(1:n_rows, j) * (this%indiv_dev_facs(1:n_rows, j) - &
          this%dev_facs(j)) ** 2) / (n_rows - 1))
      else
        this%sigmas(j) = extrapolate_sigma(this%sigmas, j)
      end if
    end do
  end if

  if (compute_resids) then
    if (this%resids_type /= LOGNORMAL) then
      this%scale_facs = 0
      do j = 1, this%n_dev - 1
        n_rows = this%n_dev - j
        do i = 1, n_rows
          this%scale_facs(i, j) = sqrt(1 - this%data(i, j) / sum(this%data(1:n_rows, j)))
        end do
      end do
      this%scale_facs(1, this%n_dev - 1) = 1
    end if

    n_resids = (this%n_dev ** 2 - this%n_dev) / 2 - 1
    this%resids = 0
    select case (this%resids_type)
     case (STANDARDISED)
      do j = 1, this%n_dev - 1
        n_rows = this%n_dev - j
        do i = 1, n_rows
          this%resids(i, j) = (this%indiv_dev_facs(i, j) - this%dev_facs(j)) * sqrt(this%data(i, j)) / &
            (this%sigmas(j) * this%scale_facs(i, j))
        end do
      end do

     case (MODIFIED)
      do j = 1, this%n_dev - 1
        n_rows = this%n_dev - j
        do i = 1, n_rows
          this%resids(i, j) = (this%indiv_dev_facs(i, j) - this%dev_facs(j)) * sqrt(this%data(i, j)) / this%scale_facs(i, j)
        end do
      end do

     case (STUDENTISED)
      this%sigmas_jack = 0
      allocate(col_mask(this%n_dev - 1))
      do j = 1, this%n_dev - 1
        n_rows = this%n_dev - j
        do i = 1, n_rows
          col_mask = .true.
          col_mask(i) = .false.
          n_pts_col = n_rows - 1
          dev_fac_jack = sum(this%data(1:n_rows, j + 1), mask=col_mask(1:n_rows)) / &
            sum(this%data(1:n_rows, j), mask=col_mask(1:n_rows))
          if (n_pts_col >= 2) then
            this%sigmas_jack(i, j) = sqrt(sum(this%data(1:n_rows, j) * &
              (this%indiv_dev_facs(1:n_rows, j) - dev_fac_jack) ** 2, mask=col_mask(1:n_rows)) / (n_pts_col - 1))
          else
            this%sigmas_jack(i, j) = extrapolate_sigma(this%sigmas_jack(i, :), j)
          end if
          this%resids(i, j) = (this%data(i, j + 1) - this%dev_facs(j) * this%data(i, j)) / &
            (this%sigmas_jack(i, j) * this%scale_facs(i, j) * sqrt(this%data(i, j)))
        end do
      end do
      
     case (LOGNORMAL)
      do j = 1, this%n_dev
        n_rows = this%n_dev - j
        do i = 1, n_rows
          this%log_normal_shift(i) = this%dev_facs(j) * sqrt(this%data(i, j)) / this%sigmas(j)
          this%log_normal_sigmas(i) = sqrt(log(1 + 1 / this%log_normal_shift(i) ** 2))
          this%log_normal_means(i) = log(this%log_normal_shift(i)) - this%log_normal_sigmas(i) ** 2 / 2
          this%resids(i, j) = (this%data(i, j + 1) - this%dev_facs(j) * this%data(i, j)) / &
            (this%sigmas(j) * sqrt(this%data(i, j)))
          this%resids(i, j) = (log(this%resids(i, j) + this%log_normal_shift(i)) - this%log_normal_means(i)) &
            / this%log_normal_sigmas(i)
        end do
      end do
    end select

    this%resids(1, this%n_dev - 1) = 0

    if (this%resids_type /= STUDENTISED) then
      resids_mean = 0
      do i = 1, this%n_dev - 1
        do j = 1, this%n_dev - i
          resids_mean = resids_mean + this%resids(i, j)
        end do
      end do
      resids_mean = resids_mean / n_resids
      do i = 1, this%n_dev - 1
        do j = 1, this%n_dev - i
          this%resids(i, j) = this%resids(i, j) - resids_mean
        end do
      end do
    end if
  end if
end subroutine fit_triangle

end module triangle_mod
