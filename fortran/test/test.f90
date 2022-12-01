program test

    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN

    implicit none

    integer, parameter :: dp = selected_real_kind(15)
    double precision :: nan
    integer, parameter :: n = 7
    double precision :: test_array(n, n), return_array(n, n)
    integer :: i, j

    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)

    test_array = &
    transpose(reshape([3511._dp, 6726._dp, 8992._dp, 10704._dp, 11763._dp, 12350._dp, 12690._dp, &
    4001._dp, 7703._dp, 9981._dp, 11161._dp, 12117._dp, 12746._dp, nan, &
    4355._dp, 8287._dp, 10233._dp, 11755._dp, 12993._dp, nan, nan, &
    4295._dp, 7750._dp, 9773._dp, 11093._dp, nan, nan, nan, &
    4150._dp, 7897._dp, 10217._dp, nan, nan, nan, nan, &
    5102._dp, 9650._dp, nan, nan, nan, nan, nan, &
    6283._dp, nan, nan, nan, nan, nan, nan], [n, n]))

    return_array = nan

    call reserve_boot(n, test_array, return_array)

    do, i = 1, n
        write(*, "(100f10.3)") (return_array(i, j), j = 1, n)
    end do

end program