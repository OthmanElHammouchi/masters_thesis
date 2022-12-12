program test
        use dispmodule
        implicit none
        real :: test_array(3) = [1., 2., 3.]
        call disp(spread(test_array, 2, 3))
end program test
