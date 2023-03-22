module sim_mod
  use, intrinsic :: iso_c_binding
  use global_mod
  implicit none
  type :: sim
    integer(c_int) :: n_sim, m_sim
    integer(c_int) :: i_sim = 0
    real(c_double), allocatable :: table(:, :)

    contains
      final :: destroy_sim
  end type

  type, extends(sim) :: mack_sim
    real(c_double), allocatable :: factors(:)
    integer(c_int), allocatable :: boot_types(:)
    integer(c_int), allocatable :: proc_dists(:)
    logical(c_bool), allocatable :: conds(:)
    integer(c_int), allocatable :: resids_types(:)

    real(c_double) :: factor = 0._c_double
    integer(c_int) :: boot_type = 0
    integer(c_int) :: proc_dist = 0
    logical(c_bool) :: cond = .false.
    integer(c_int) :: resids_type = 0

  contains
    final :: destroy_mack
  end type

  type, extends(mack_sim) :: single_sim
    integer(c_int) :: outlier_rowidx, excl_rowidx
    integer(c_int) :: outlier_colidx, excl_colidx

    contains
      procedure :: update => update_single
  end type single_sim

  type, extends(mack_sim) :: calendar_sim
    integer(c_int) :: outlier_diagidx, excl_diagidx

    contains
    procedure :: update => update_calendar
  end type calendar_sim

  type, extends(mack_sim) :: origin_sim
    integer(c_int) :: outlier_rowidx, excl_rowidx

    contains
    procedure :: update => update_origin
  end type origin_sim

  interface single_sim
    procedure :: init_single
  end interface

  interface calendar_sim
    procedure :: init_calendar
  end interface

  interface origin_sim
    procedure :: init_origin
  end interface

contains

  function init_single(n_dev, factors, boot_types, proc_dists, conds, resids_types) result(this)
    integer(c_int), intent(in) :: n_dev
    real(c_double), intent(in) :: factors(:)
    integer(c_int), intent(in) :: boot_types(:)
    integer(c_int), intent(in) :: proc_dists(:)
    logical(c_bool), intent(in) :: conds(:)
    integer(c_int), intent(in) :: resids_types(:)

    integer(c_int) :: n_outliers, n_excl, n_sim_, m_sim_, n_res, m_res, n_sim, m_sim
    integer(c_int) :: n_rows
    integer(c_int) :: n_factors, n_boot_types, n_proc_dists, n_conds, n_resids_types
    integer(c_int) :: i, j, k, i1, i2, i3, i4, i5, i6, i7
    integer(c_int), allocatable :: outliers(:, :), excl(:, :)
    real(c_double), allocatable :: table(:, :)
    type(single_sim) :: this

    n_outliers = (n_dev ** 2 - n_dev) / 2 - 1
    n_excl = n_outliers

    n_factors = size(factors)
    n_boot_types = size(boot_types)
    n_proc_dists = size(proc_dists)
    n_conds = size(conds)
    n_resids_types = size(resids_types)

    n_sim_ = n_outliers * n_excl
    n_sim_ = n_sim_ * n_factors
    n_sim_ = n_sim_ * n_boot_types
    n_sim_ = n_sim_ * n_proc_dists
    n_sim_ = n_sim_ * n_conds
    n_sim_ = n_sim_ * n_resids_types
    m_sim_ = 9

    allocate(outliers(n_outliers, 2), source=0)
    allocate(excl(n_excl, 2), source=0)
    allocate(table(n_sim_, m_sim_), source=0._c_double)

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
            proc_loop: do i5 = 1, n_proc_dists
              cond_loop: do i6 = 1, n_conds
                do i7 = 1, n_resids_types
                  table(k, 1:2) = outliers(i1, :)
                  table(k, 3) = factors(i2)
                  table(k, 4:5) = excl(i3, :)
                  table(k, 6) = boot_types(i4)
                  table(k, 7) = proc_dists(i5)

                  if (table(k, 6) == PAIRS) then
                    table(k, 8) = NONE
                    table(k, 9) = NONE
                    k = k + 1
                    cycle proc_loop

                  else if (table(k, 6) == PARAMETRIC) then
                    table(k, 8) = merge(1._c_double, 0._c_double, conds(i6))
                    table(k, 9) = NONE
                    k = k + 1
                    cycle cond_loop
                  else
                    table(k, 8) = merge(1._c_double, 0._c_double, conds(i6))
                    table(k, 9) = resids_types(i7)
                    k = k + 1
                  end if
                end do
              end do cond_loop
            end do proc_loop
          end do
        end do
      end do
    end do

    this%n_sim = k - 1
    this%m_sim = m_sim_

    allocate(this%table(n_sim, m_sim))
    this%table = table(1:n_sim, :)
    deallocate(table)
    deallocate(outliers)
    deallocate(excl)

    this%factors = factors
    this%boot_types = boot_types
    this%proc_dists = proc_dists
    this%conds = conds
    this%resids_types = resids_types
  end function init_single

  function init_calendar(n_dev, factors, boot_types, proc_dists, conds, resids_types) result(this)
    integer(c_int), intent(in) :: n_dev
    real(c_double), intent(in) :: factors(:)
    integer(c_int), intent(in) :: boot_types(:)
    integer(c_int), intent(in) :: proc_dists(:)
    logical(c_bool), intent(in) :: conds(:)
    integer(c_int), intent(in) :: resids_types(:)

    integer(c_int), allocatable :: outliers(:, :), excl(:, :)
    type(calendar_sim) :: this

    integer(c_int) :: n_outliers, n_excl, n_sim_, m_sim_, n_res, m_res, n_sim, m_sim
    integer(c_int) :: n_rows
    integer(c_int) :: n_factors, n_boot_types, n_proc_dists, n_conds, n_resids_types
    integer(c_int) :: i, j, k, i1, i2, i3, i4, i5, i6, i7
    real(c_double), allocatable :: table(:, :)

    n_outliers = n_dev - 1
    n_excl = n_outliers

    n_sim_ = n_sim_ * n_outliers * n_excl
    n_sim_ = n_sim_ * n_factors
    n_sim_ = n_sim_ * n_boot_types
    n_sim_ = n_sim_ * n_proc_dists
    n_sim_ = n_sim_ * n_conds
    n_sim_ = n_sim_ * n_resids_types
    m_sim_ = 7

    allocate(table(n_sim_, m_sim_), source=0._c_double)

    k = 1
    do i1 = 1, n_outliers
      do i2 = 1, n_factors
        do i3 = 1, n_excl
          do i4 = 1, n_boot_types
            proc_loop: do i5 = 1, n_proc_dists
              cond_loop: do i6 = 1, n_conds
                do i7 = 1, n_resids_types
                  table(k, 1) = i1
                  table(k, 2) = factors(i2)
                  table(k, 3) = i3
                  table(k, 4) = boot_types(i4)
                  table(k, 5) = proc_dists(i5)

                  if (table(k, 4) == PAIRS) then
                    table(k, 6) = NONE
                    table(k, 7) = NONE
                    k = k + 1
                    cycle proc_loop

                  else if (table(k, 4) == PARAMETRIC) then
                    table(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                    table(k, 7) = NONE
                    k = k + 1
                    cycle cond_loop

                  else
                    table(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                    table(k, 7) = resids_types(i7)
                    k = k + 1
                  end if
                end do
              end do cond_loop
            end do proc_loop
          end do
        end do
      end do
    end do

    this%n_sim = k - 1
    this%m_sim = m_sim_

    allocate(this%table(n_sim, m_sim))
    this%table = table(1:n_sim, :)
    deallocate(table)
    deallocate(outliers)
    deallocate(excl)

    this%n_sim = n_sim
    this%m_sim = m_sim
    this%factors = factors
    this%boot_types = boot_types
    this%proc_dists = proc_dists
    this%conds= conds
    this%resids_types = resids_types
  end function init_calendar

  function init_origin(n_dev, factors, boot_types, proc_dists, conds, resids_types) result(this)
    integer(c_int), intent(in) :: n_dev
    real(c_double), intent(in) :: factors(:)
    integer(c_int), intent(in) :: boot_types(:)
    integer(c_int), intent(in) :: proc_dists(:)
    logical(c_bool), intent(in) :: conds(:)
    integer(c_int), intent(in) :: resids_types(:)

    integer(c_int), allocatable :: outliers(:, :), excl(:, :)
    type(origin_sim) :: this

    integer(c_int) :: n_outliers, n_excl, n_sim_, m_sim_, m_res, n_res, n_sim, m_sim
    integer(c_int) :: n_rows
    integer(c_int) :: n_factors, n_boot_types, n_proc_dists, n_conds, n_resids_types
    integer(c_int) :: i, j, k, i1, i2, i3, i4, i5, i6, i7
    real(c_double), allocatable :: table(:, :)

    n_outliers = n_dev
    n_excl = n_outliers

    n_sim_ = n_sim_ * n_outliers * n_excl
    n_sim_ = n_sim_ * n_factors
    n_sim_ = n_sim_ * n_boot_types
    n_sim_ = n_sim_ * n_proc_dists
    n_sim_ = n_sim_ * n_conds
    n_sim_ = n_sim_ * n_resids_types
    m_sim_ = 7

    allocate(table(n_sim_, m_sim_), source=0._c_double)

    k = 1
    do i1 = 1, n_outliers
      do i2 = 1, n_factors
        do i3 = 1, n_excl
          do i4 = 1, n_boot_types
            proc_loop: do i5 = 1, n_proc_dists
              cond_loop: do i6 = 1, n_conds
                do i7 = 1, n_resids_types
                  table(k, 1) = i1
                  table(k, 2) = factors(i2)
                  table(k, 3) = i3
                  table(k, 4) = boot_types(i4)
                  table(k, 5) = proc_dists(i5)

                  if (table(k, 4) == PAIRS) then
                    table(k, 6) = NONE
                    table(k, 7) = NONE
                    k = k + 1
                    cycle proc_loop

                  else if (table(k, 4) == PARAMETRIC) then
                    table(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                    table(k, 7) = NONE
                    k = k + 1
                    cycle cond_loop

                  else
                    table(k, 6) = merge(1._c_double, 0._c_double, conds(i6))
                    table(k, 7) = resids_types(i7)
                    k = k + 1
                  end if
                end do
              end do cond_loop
            end do proc_loop
          end do
        end do
      end do
    end do

  this%n_sim = k - 1
  this%m_sim = m_sim_

  allocate(this%table(n_sim, m_sim))
  this%table = table(1:n_sim, :)
  deallocate(table)
  deallocate(outliers)
  deallocate(excl)

  this%factors = factors
  this%boot_types = boot_types
  this%proc_dists = proc_dists
  this%conds = conds
  this%resids_types = resids_types
end function init_origin

subroutine destroy_mack(this)
  type(mack_sim), intent(inout) :: this
  deallocate(this%factors)
  deallocate(this%boot_types)
  deallocate(this%proc_dists)
  deallocate(this%conds)
  deallocate(this%resids_types)
end subroutine destroy_mack

subroutine destroy_sim(this)
  type(sim), intent(inout) :: this
  deallocate(this%table)
end subroutine destroy_sim

subroutine update_single(this)
  class(single_sim), intent(inout) :: this

  this%i_sim = this%i_sim + 1
  this%outlier_rowidx = int(this%table(this%i_sim, 1))
  this%outlier_colidx = int(this%table(this%i_sim, 2))
  this%excl_rowidx = this%table(this%i_sim, 4)
  this%excl_colidx = this%table(this%i_sim, 5)
  this%factor = this%table(this%i_sim, 3)
  this%boot_type = int(this%table(this%i_sim, 6))
  this%proc_dist = int(this%table(this%i_sim, 7))
  this%cond = int(this%table(this%i_sim, 8))
  this%resids_type = int(this%table(this%i_sim, 9))
end subroutine update_single

subroutine update_calendar(this)
  class(calendar_sim), intent(inout) :: this

  this%i_sim = this%i_sim + 1
  this%boot_type = int(this%table(this%i_sim, 4))
  this%proc_dist = int(this%table(this%i_sim, 5))
  this%cond = int(this%table(this%i_sim, 6))
  this%resids_type = int(this%table(this%i_sim, 7))
  this%factor = this%table(this%i_sim, 2)
  this%outlier_diagidx = int(this%table(this%i_sim, 1))
  this%excl_diagidx = int(this%table(this%i_sim, 3))
end subroutine update_calendar

subroutine update_origin(this)
  class(origin_sim), intent(inout) :: this

  this%i_sim = this%i_sim + 1
  this%boot_type = int(this%table(this%i_sim, 4))
  this%proc_dist = int(this%table(this%i_sim, 5))
  this%cond = int(this%table(this%i_sim, 6))
  this%resids_type = int(this%table(this%i_sim, 7))
  this%factor = this%table(this%i_sim, 2)
  this%outlier_rowidx = int(this%table(this%i_sim, 1))
  this%excl_rowidx = int(this%table(this%i_sim, 3))
end subroutine update_origin

end module sim_mod
