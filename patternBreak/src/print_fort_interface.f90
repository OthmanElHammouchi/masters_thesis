module print_fort_interface

    use, intrinsic :: iso_c_binding

    implicit none

    interface

        subroutine Rprintf(string) bind(C, name="Rprintf_")
            
            import c_char

            character(c_char), intent(in) :: string(*)

        end subroutine

        subroutine R_FlushConsole() bind(C, name="R_FlushConsole")

        end subroutine R_FlushConsole

    end interface
    
end module print_fort_interface