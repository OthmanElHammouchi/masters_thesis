subroutine reserve_boot(n_dev, triangle, boot_triangle)

  use random, only: norm => random_normal

  implicit none

  integer, intent(in) :: n_dev
  double precision, intent(in) :: triangle(n_dev, n_dev)
  double precision, intent(inout) :: boot_triangle(n_dev, n_dev)

  integer :: i, j, nrows
  double precision :: dev_facs(n_dev - 1), sigmas(n_dev - 1)
  double precision :: indiv_dev_facs(n_dev - 1, n_dev - 1), resids(n_dev - 1, n_dev - 1)

  double precision, allocatable :: flat_resids(:)
  integer :: n_resids

  double precision :: boot_resids(n_dev - 1, n_dev - 1), boot_indiv_dev_facs(n_dev - 1, n_dev - 1), boot_dev_facs(n_dev - 1), &
  boot_sigmas(n_dev - 1)

  do j = 1, n_dev - 1
    
    nrows = n_dev - j

    indiv_dev_facs(1:nrows, j) = triangle(1:nrows, j + 1) / triangle(1:nrows, j)

    dev_facs(j) = sum(triangle(1:nrows, j + 1)) / sum(triangle(1:nrows, j))

    if (j < n_dev - 1) then
      sigmas(j) = sum(triangle(1:nrows, j) * (indiv_dev_facs(1:nrows, j) - dev_facs(j))) / nrows
    else 
      sigmas = sqrt(min(sigmas(j - 1) ** 2, sigmas(j - 2) ** 2, sigmas(j - 1) ** 4 / sigmas(j - 2) ** 2))
    end if

    resids(1:nrows, j) = (indiv_dev_facs(1:nrows, j) - dev_facs(j)) * sqrt(triangle(1:nrows, j)) / sigmas(j)

  end do

  flat_resids = pack(resids, .true.)
  n_resids = size(flat_resids)

  do j = 1, n_dev - 1
    do i = 1, n_dev - j
      boot_resids(i, j) = flat_resids(1 + int(n_resids * rand()))
    end do
  end do

  do j = 1, n_dev - 1

    nrows = n_dev - j

    boot_indiv_dev_facs(1:nrows, j) = dev_facs(j) + boot_resids(1:nrows, j) * sigmas(j) / sqrt(triangle(1:nrows, j))

    boot_dev_facs(j) = sum(triangle(1:nrows, j) * boot_indiv_dev_facs(1:nrows, j)) / sum(triangle(1:nrows, j))

    if (j < n_dev - 1) then
      boot_sigmas(j) = sum(triangle(1:nrows, j) * (boot_indiv_dev_facs(1:nrows, j) - boot_dev_facs(j))) / nrows
    else 
      boot_sigmas = sqrt(min(boot_sigmas(j - 1) ** 2, boot_sigmas(j - 2) ** 2, boot_sigmas(j - 1) ** 4 / boot_sigmas(j - 2) ** 2))
    end if

  end do

  boot_triangle = triangle

  do j = 2, n_dev
    do i = n_dev + 2 - j, n_dev
      boot_triangle(i, j) = (boot_dev_facs(j - 1) * boot_triangle(i, j - 1)) + &
      dble(norm()) * (boot_sigmas(j - 1) * sqrt(boot_triangle(i, j - 1)))
    end do
  end do

end subroutine