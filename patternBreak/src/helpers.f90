module helpers

   use iso_c_binding
   use constants
   use rng_fort_interface

   implicit none

contains

   function single_outlier(outlier_rowidx, outlier_colidx, factor, init_col, dev_facs, sigmas, dist) result(triangle)

      integer, intent(in):: outlier_rowidx, outlier_colidx
      real(c_double), intent(in) :: init_col(:), dev_facs(:), sigmas(:)
      integer(c_int), intent(in) :: dist

      real(c_double) :: factor
      real(c_double), allocatable:: triangle(:, :)

      real(c_double) :: shape, scale
      real(c_double) :: mean, sd

      integer :: n_dev, i, j

      n_dev = size(init_col)

      allocate(triangle(n_dev, n_dev), source=0._c_double)
      triangle(:, 1) = init_col

      if (dist == NORMAL) then

         do j = 2, n_dev
            do i = 1, n_dev + 1 - j
               if (i == outlier_rowidx) cycle
               mean = dev_facs(j - 1) * triangle(i, j - 1)
               sd = sigmas(j - 1) * sqrt(triangle(i, j - 1))
               triangle(i, j) = rnorm(mean, sd)
            end do
         end do

         if (outlier_colidx > 2) then
            do j = 2, outlier_colidx - 1
               mean = dev_facs(j - 1) * triangle(outlier_rowidx, j - 1)
               sd = sigmas(j - 1) * sqrt(triangle(outlier_rowidx, j - 1))
               triangle(outlier_rowidx, j) = rnorm(mean, sd)
            end do
         end if

         mean = factor * dev_facs(outlier_colidx - 1) * triangle(outlier_rowidx, outlier_colidx - 1)
         sd = sigmas(outlier_colidx - 1) * sqrt(triangle(outlier_rowidx, outlier_colidx - 1))
         triangle(outlier_rowidx, outlier_colidx) = rnorm(mean, sd)

         if (outlier_colidx < n_dev) then
            do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
               mean = dev_facs(j - 1) * triangle(outlier_rowidx, j - 1)
               sd = sigmas(j - 1) * sqrt(triangle(outlier_rowidx, j - 1))
               triangle(outlier_rowidx, j) = rnorm(mean, sd)
            end do
         end if

      else if (dist == GAMMA) then

         do j = 2, n_dev
            do i = 1, n_dev + 1 - j
               if (i == outlier_rowidx) cycle
               shape = dev_facs(j - 1)**2 * triangle(i, j - 1) / sigmas(j - 1)**2
               scale = sigmas(j - 1)**2 / dev_facs(j - 1)
               triangle(i, j) = rgamma(shape, scale)
            end do
         end do

         if (outlier_colidx > 2) then
            do j = 2, outlier_colidx - 1
               shape = dev_facs(j - 1)**2 * triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
               scale = sigmas(j - 1)**2 / dev_facs(j - 1)
               triangle(outlier_rowidx, j) = rgamma(shape, scale)
            end do
         end if

         shape = dev_facs(outlier_colidx - 1)**2 * triangle(outlier_rowidx, outlier_colidx - 1) / sigmas(outlier_colidx - 1)**2
         scale = sigmas(outlier_colidx - 1)**2 / dev_facs(outlier_colidx - 1)
         triangle(outlier_rowidx, outlier_colidx) = rgamma(shape, scale)

         if (outlier_colidx < n_dev) then
            do j = outlier_colidx + 1, n_dev + 1 - outlier_rowidx
               shape = dev_facs(j - 1)**2 * triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
               scale = sigmas(j - 1)**2 / dev_facs(j - 1)
               triangle(outlier_rowidx, j) = rgamma(shape, scale)
            end do
         end if

      end if

   end function single_outlier

   function calendar_outlier(outlier_diagidx, factor, triangle, dev_facs, sigmas, dist) result(sim_triangle)

      integer, intent(in) :: outlier_diagidx
      real(c_double), intent(in):: factor
      real(c_double), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
      integer(c_int), intent(in) :: dist

      integer :: i, j, n_dev, n_cols
      real(c_double), allocatable :: sim_triangle(:, :)
      real(c_double) :: shape, scale
      real(c_double) :: mean, sd

      n_dev = size(triangle, dim=1)

      allocate(sim_triangle(n_dev, n_dev), source=0._c_double)
      sim_triangle(:, 1) = triangle(:, 1)

      do i = 1, n_dev
         n_cols = n_dev + 2 - outlier_diagidx - i
         if (n_cols <= 1) then
            do j = 2, n_dev + 1 - i
               sim_triangle(i, j) = triangle(i, j)
            end do
         else
            do j = 2, n_cols - 1
               sim_triangle(i, j) = triangle(i, j)
            end do

            if (dist == NORMAL) then

               mean = factor * dev_facs(n_cols - 1) * sim_triangle(i, n_cols - 1)
               sd = sigmas(n_cols - 1) * sqrt(sim_triangle(i, n_cols - 1))

               sim_triangle(i, n_cols) = rnorm(mean, sd)

               do j = n_cols + 1, n_dev + 1 - i
                  mean = dev_facs(j - 1) * sim_triangle(i, j - 1)
                  sd = sigmas(j - 1) * sqrt(sim_triangle(i, j - 1))

                  sim_triangle(i, j) = rnorm(mean, sd)
               end do

            else if (dist == GAMMA) then

               shape = factor * dev_facs(n_cols - 1)**2 * sim_triangle(i, n_cols - 1) / sigmas(n_cols - 1)**2
               scale = sigmas(n_cols - 1)**2 / factor * dev_facs(n_cols - 1)

               sim_triangle(i, n_cols) = rgamma(shape, scale)

               do j = n_cols + 1, n_dev + 1 - i

                  shape = dev_facs(j - 1)**2 * sim_triangle(i, j - 1) / sigmas(j - 1)**2
                  scale = sigmas(j - 1)**2 / dev_facs(j - 1)

                  sim_triangle(i, j) = rgamma(shape, scale)

               end do

            end if

         end if
      end do

   end function calendar_outlier

   function origin_outlier(outlier_rowidx, factor, triangle, dev_facs, sigmas, dist) result(sim_triangle)

      integer, intent(in):: outlier_rowidx
      real(c_double), intent(in) :: triangle(:, :), dev_facs(:), sigmas(:)
      integer(c_int), intent(in) :: dist
      real(c_double), intent(in) :: factor

      real(c_double) :: shape, scale
      real(c_double) :: mean, sd
      real(c_double), allocatable:: sim_triangle(:, :)

      integer(c_int) :: n_dev, i, j

      n_dev = size(triangle, dim=1)
      sim_triangle = triangle

      do j = 2, n_dev + 1 - outlier_rowidx

         if (dist == NORMAL) then

            mean = factor * dev_facs(j - 1) * sim_triangle(outlier_rowidx, j - 1)
            sd = sigmas(j - 1) * sqrt(sim_triangle(outlier_rowidx, j - 1))
            sim_triangle(outlier_rowidx, j) = rnorm(mean, sd)

         else if (dist == GAMMA) then

            shape = factor * dev_facs(j - 1)**2 * sim_triangle(outlier_rowidx, j - 1) / sigmas(j - 1)**2
            scale = sigmas(j - 1)**2 / factor * dev_facs(j - 1)
            sim_triangle(outlier_rowidx, j) = rgamma(shape, scale)

         end if
      end do

   end function origin_outlier

end module helpers
