! module helpers

!    use iso_c_binding

!    implicit none

!    integer, parameter :: NORMAL = 1
!    integer, parameter :: GAMMA = 2
!    integer, parameter :: CONDITIONAL = 1
!    integer, parameter :: UNCONDITIONAL = 2
!    integer, parameter :: RAW = 1
!    integer, parameter :: SCALED = 2
!    integer, parameter :: PARAMETRIC = 3

!    include "c_interface.h"

! contains

!    function single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist) result(sim_triangle)

!       integer, intent(in):: outlier_rowidx, outlier_colidx
!       real(c_double), intent(in) :: init_col(:), dev_facs(:), sigmas(:)
!       character(:), allocatable, intent(in) :: dist

!       real(c_double) :: factor
!       real(c_double), allocatable:: sim_triangle(:, :)

!       real(c_double) :: shape, scale
!       real(c_double) :: mean, sd

!       integer :: n_dev, i, j

!       n_dev = size(init_col)

!       allocate(sim_triangle(n_dev, n_dev))
!       sim_triangle = 0

!       sim_triangle(:, 1) = init_col

!       if (dist == NORMAL) then

!          do j = 2, n_dev
!             do i = 1, n_dev + 1 - j
!                if (i == outlier_rowidx) cycle
!                mean = dev_facs(j - 1) * sim_triangle(i, j - 1)
!                sd = sigmas(j - 1) * sqrt(sim_triangle(i, j - 1))
!                sim_triangle(i, j) = rnorm(mean, sd)
!             end do
!          end do

!          if (outlier_colidx > 2) then
!             do j = 2, outlier_colidx - 1
!                mean = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1)
!                sd = sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
!                sim_triangle(outlier_rowidx, j) = rnorm(mean, sd)
!             end do
!          end if

!          mean = factor * dev_facs(outlier_colidx - 1) * sim_triangle(outlier_rowidx, outlier_colidx - 1)
!          sd = sigmas(outlier_colidx - 1) * sqrt(sim_triangle(outlier_rowidx, outlier_colidx - 1))
!          sim_triangle(outlier_rowidx, outlier_colidx) = rnorm(mean, sd)

!          if (outlier_colidx < n_dev) then
!             do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
!                mean = dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1)
!                sd = sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
!                sim_triangle(outlier_rowidx, j) = rnorm(mean, sd)
!             end do
!          end if

!       else if (dist == GAMMA) then

!          do j = 2, n_dev
!             do i = 1, n_dev + 1 - j
!                if (i == outlier_rowidx) cycle
!                shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
!                scale = sigmas(j - 1)**2 / dev_facs(j - 1)
!                sim_triangle(i, j) = rgamma(shape, scale)
!             end do
!          end do

!          if (outlier_colidx > 2) then
!             do j = 2, outlier_colidx - 1
!                shape = dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
!                scale = sigmas(j - 1)**2 / dev_facs(j - 1)
!                sim_triangle(outlier_rowidx, j) = rgamma(shape, scale)
!             end do
!          end if

!          shape = dev_facs(outlier_colidx - 1)**2 * sim_triangle(outlier_rowidx, outlier_colidx - 1) / sigmas(outlier_colidx - 1)**2
!          scale = sigmas(outlier_colidx - 1)**2 / dev_facs(outlier_colidx - 1)
!          sim_triangle(outlier_rowidx, outlier_colidx) = rgamma(shape, scale)

!          if (outlier_colidx < n_dev) then
!             do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
!                shape = dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
!                scale = sigmas(j - 1)**2 / dev_facs(j - 1)
!                sim_triangle(outlier_rowidx, j) = rgamma(shape, scale)
!             end do
!          end if

!       end if

!    end function single_outlier

!    function calendar_outlier(outlier_diagidx, factor, triangle, dev_facs, sigmas, dist) result(sim_triangle)

!       integer, intent(in) :: outlier_diagidx
!       real(c_double), intent(in):: factor
!       real(c_double), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
!       character(:), allocatable, intent(in) :: dist

!       integer :: i, j, n_dev, n_cols
!       real(c_double), allocatable :: sim_triangle(:, :)
!       real(c_double) :: gamma_shape, gamma_rate

!       n_dev = size(triangle, dim=1)

!       allocate(sim_triangle(n_dev, n_dev), source=0._c_double)
!       sim_triangle(:, 1) = triangle(:, 1)

!       do i = 1, n_dev
!          n_cols = n_dev + 2 - outlier_diagidx - i
!          if (n_cols <= 1) then
!             do j = 2, n_dev + 1 - i
!                sim_triangle(i, j) = triangle(i, j)
!             end do
!          else
!             do j = 2, n_cols - 1
!                sim_triangle(i, j) = triangle(i, j)
!             end do

!             if (dist == "normal") then

! #ifdef __INTEL_COMPILER
!                method = VSL_RNG_METHOD_GAUSSIAN_ICDF

!                errcode = vdrnggaussian(method, stream, 1, r, 0._c_double, 1._c_double)

!                sim_triangle(i, n_cols) = factor * dev_facs(n_cols - 1) * sim_triangle(i, n_cols - 1) + sigmas(n_cols - 1) * sqrt(sim_triangle(i, n_cols - 1)) * r(1)

! #elif defined __GFORTRAN__

!                sim_triangle(i, n_cols) = factor * dev_facs(n_cols - 1) * sim_triangle(i, n_cols - 1) + sigmas(n_cols - 1) * sqrt(sim_triangle(i, n_cols - 1)) * norm()

! #endif
!                do j = n_cols + 1, n_dev + 1 - i
! #ifdef __INTEL_COMPILER
!                   errcode = vdrnggaussian(method, stream, 1, r, 0._c_double, 1._c_double)

!                   sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(i, j - 1)) * r(1)

! #elif defined __GFORTRAN__

!                   sim_triangle(i, j) = dev_facs(j - 1) * sim_triangle(i, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(i, j - 1)) * norm()

! #endif
!                end do
!             else if (dist == "gamma") then

!                gamma_shape = factor * dev_facs(n_cols - 1)**2 * sim_triangle(i, n_cols - 1) / sigmas(n_cols - 1)**2
!                gamma_rate = factor * dev_facs(n_cols - 1) / sigmas(n_cols - 1)**2

! #ifdef __INTEL_COMPILER

!                method = VSL_RNG_METHOD_GAMMA_GNORM
!                errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._c_double, 1/gamma_rate)

!                sim_triangle(i, n_cols) = r(1)

! #elif __GFORTRAN__

!                sim_triangle(i, n_cols) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), c_double)

! #endif
!                do j = n_cols + 1, n_dev + 1 - i

!                   gamma_shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
!                   gamma_rate = dev_facs(j - 1) / sigmas(j - 1)**2

! #ifdef __INTEL_COMPILER

!                   errcode = vdrnggamma(method, stream, 1, r, gamma_shape, 0._c_double, 1/gamma_rate)

!                   sim_triangle(i, j) = r(1)

! #elif __GFORTRAN__

!                   sim_triangle(i, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), c_double)

! #endif
!                end do
!             end if
!          end if
!       end do

! #ifdef __INTEL_COMPILER
!       errcode = vsldeletestream(stream)
! #endif

!    end function calendar_outlier

!    function origin_outlier(outlier_rowidx, factor, triangle, dev_facs, sigmas, dist) result(sim_triangle)

!       integer, intent(in):: outlier_rowidx
!       real(c_double), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
!       character(:), allocatable, intent(in) :: dist
!       real(c_double), intent(in) :: factor

!       real(c_double) :: gamma_shape, gamma_rate
!       real(c_double), allocatable:: sim_triangle(:, :)

!       integer :: n_dev, i, j

! #ifdef __INTEL_COMPILER
!       type(vsl_stream_state) :: stream
!       integer(kind=4) :: errcode
!       integer :: brng, method, sd, n
!       real(c_double) :: r(1)

!       brng = VSL_BRNG_MT19937
!       sd = irand()
!       errcode = vslnewstream(stream, brng,  sd)
! #endif

!       n_dev = size(triangle, dim=1)
!       sim_triangle = triangle

!       do j = 2, n_dev + 1 - outlier_rowidx

!          if (dist == "normal") then

! #ifdef __INTEL_COMPILER
!             method = VSL_RNG_METHOD_GAUSSIAN_ICDF

!             errcode = vdrnggaussian(method, stream, 1, r, 0._c_double, 1._c_double)

!             sim_triangle(outlier_rowidx, j) = factor * dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1)) * r(1)

! #elif defined __GFORTRAN__

!             sim_triangle(outlier_rowidx, j) = factor * dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1) + sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1)) * norm()

! #endif

!          else if (dist == "gamma") then

!             gamma_shape = factor * dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
!             gamma_rate = factor * dev_facs(j - 1) / sigmas(j - 1)**2

! #elif __GFORTRAN__

!             sim_triangle(outlier_rowidx, j) = real(gamma(real(gamma_shape), .true.)/real(gamma_rate), c_double)

! #endif
!          end if
!       end do

!    end function origin_outlier

! end module helpers
