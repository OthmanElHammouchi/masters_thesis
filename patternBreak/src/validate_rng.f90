module validate_rng
   use interface
   use helpers
   use constants
   use omp_lib
   use iso_c_binding

   implicit none

contains
   function validate_rng_f(n_samples) result(res) bind(c)
      integer(c_int), intent(in), value :: n_samples
      integer(c_int) :: res
      integer(c_int) :: i, j, k, l
      real(c_double), allocatable :: samples(:, :)
      real(c_double), allocatable :: flattened(:)
      integer(c_int) :: n_threads, i_thread, total
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

end function validate_rng_f
end module validate_rng
