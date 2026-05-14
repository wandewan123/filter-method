module utils_bc
    implicit none
contains

    pure function get_phi(phi, i, n) result(val)
        integer, intent(in) :: i, n
        real(8), intent(in) :: phi(n)
        real(8) :: val

        if (i < 1 .or. i > n) then
            val = 0.0d0   ! Dirichlet ghost cell
        else
            val = phi(i)
        end if
    end function get_phi

end module utils_bc
