module glm

    use, intrinsic :: iso_c_binding
    use constants
    use rng_fort_interface

   implicit none

   ! include "rng_fort_interface.h"

contains

   subroutine glm_boot(triangle, n_dev, n_boot, reserve, excl_resids)

      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: reserve(n_boot)

      integer(c_int), intent(in) :: n_dev, n_boot

      real(c_double), intent(in), optional :: excl_resids(:, :)

      integer(c_int) :: n_pts, n_covs, n_pts_boot
      real(c_double), allocatable :: y(:), X(:, :), b(:), b_boot(:), X_lower(:, :), y_fit_boot(:)
      real(c_double), allocatable :: y_fit(:), triangle_fit(:, :), triangle_boot(:, :), y_boot(:)
      real(c_double), allocatable :: resids(:, :), resids_boot(:, :), flat_resids(:)
      logical(c_bool), allocatable :: resids_mask(:, :)

      integer(c_int) :: i, j, k, i_boot, info

      n_pts = (n_dev**2 + n_dev) / 2
      n_covs = 2*n_dev - 1

      n_pts_boot = n_dev ** 2 - n_pts

      allocate(y(n_pts), X(n_pts, n_covs), X_lower(n_pts_boot, n_covs), b(n_covs), source=0._c_double)

      allocate(resids_boot(n_pts, n_pts), b_boot(n_covs), triangle_boot(n_pts, n_pts), y_fit_boot(n_pts_boot), source = 0._c_double)

      allocate(y_boot(n_pts), y_fit(n_pts), triangle_fit(n_pts, n_pts), source=0._c_double)

      allocate(resids(n_dev, n_dev), flat_resids(n_pts), source=0._c_double)

      allocate(resids_mask(n_dev, n_dev))

      ! Set up feature matrix and response vector.
      X(:, 1) = 1._c_double

      k = 1
      do i = 1, n_dev
         do j = 1, n_dev + 1 - i
            if (i /= 1) X(k, i) = 1
            if (j /= 1) X(k, n_dev + j - 1) = 1
            y(k) = triangle(i, j)
            k = k + 1
         end do
      end do

      y(1) = triangle(1, 1)

      call poisson_fit(X, y, n_pts, n_covs, b)

      y_fit = exp(matmul(X, b))

      k = 1
      do i = 1, n_dev
         do j = 1, n_dev + 1 - i
            triangle_fit(i, j) = y_fit(k)
            k = k + 1
         end do
      end do

      do i = 1, n_dev
         do j = 1, n_dev + 1 - i
            resids(i, j) = (triangle(i, j) - triangle_fit(i, j)) / sqrt(triangle_fit(i, j))
         end do
      end do

      resids_mask = .true.

      flat_resids = pack(resids, resids_mask)

      X_lower(:, 1) = 1._c_double

      k = 1
      do i = 2, n_dev
         do j = n_dev + 1 - i + 1, n_dev
            X_lower(k, i) = 1
            if (j /= 1) X_lower(k, n_dev + j - 1) = 1
            k = k + 1
         end do
      end do

      call GetRNGstate()

      do i_boot = 1, n_boot

         do i = 1, n_dev
            do j = 1, n_dev + 1 - i
               resids_boot(i, j) = flat_resids(1 + int(n_pts * rand()))
               triangle_boot(i, j) = resids_boot(i, j) * sqrt(triangle_fit(i, j)) + triangle_fit(i, j)
            end do
         end do

         k = 1
         do i = 1, n_dev
            do j = 1, n_dev + 1 - i
               y_boot(k) = triangle_boot(i, j)
               k = k + 1
            end do
         end do

         call poisson_fit(X, y_boot, n_pts, n_covs, b_boot)

         y_fit_boot = exp(matmul(X_lower, b_boot))

         do i = 1, n_pts_boot
            y_fit_boot(i) = rpois(y_fit_boot(i))
         end do

         reserve(i_boot) = sum(y_fit_boot)

      end do

      call PutRNGstate()

   end subroutine glm_boot

   subroutine poisson_fit(X, y, N, p, b)

    integer(c_int), intent(in) :: N, p
    real(c_double), intent(in) :: X(N, p), y(N)
    real(c_double), intent(inout) :: b(p)

    real(c_double) :: diff, eps
    real(c_double) :: IPIV(p, p), rhs(p), lhs(p, p)
    real(c_double) :: W(N, N), z(N), eta(N)

    integer(c_int) :: i, info

    ! Initialise coefficients.
    b = spread(0._c_double, 1, p)
    b(1) = log(sum(y) / N)

    ! Apply IRLS to fit Poisson GLM.
    diff = 1E6
    eps = 1E-6
    do while (diff > eps)

       eta = matmul(X, b)

       W = 0

       do i = 1, N
          W(i, i) = exp(eta(i))
       end do

       z = eta + exp(-eta)*y - 1

       lhs = matmul(matmul(transpose(X), W), X)
       rhs = matmul(matmul(transpose(X), W), z)

       call dgesv(p, 1, lhs, p, IPIV, rhs, p, info)

       diff = norm2(b - rhs)

       b = rhs

    end do

 end subroutine poisson_fit

   subroutine glm_boot_centry(n_dev, triangle, n_boot, reserve) bind(C, name='glm_boot_')

      real(c_double), intent(in) :: triangle(n_dev, n_dev)
      real(c_double), intent(inout) :: reserve(n_boot)

      integer(c_int), intent(in), value :: n_dev, n_boot

      call glm_boot(triangle, n_dev, n_boot, reserve)

   end subroutine glm_boot_centry

end module glm
