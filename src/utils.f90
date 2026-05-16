module utils
    implicit none
contains

    pure function get_phi(phi, i, n) result(val)
        integer, intent(in) :: i, n
        real(8), intent(in) :: phi(n)
        real(8) :: val

        if (i < 1 .or. i > n) then
            val = 0.0d0
        else
            val = phi(i)
        end if
    end function get_phi

    function trapz(y, dx) result(integral)
        real(8), intent(in) :: y(:)
        real(8), intent(in) :: dx
        real(8)             :: integral
        integer             :: n, i
        
        n        = size(y)
        integral = 0.0d0
        
        if (n > 1) then
            integral = y(1) + y(n)
            do i = 2, n - 1
                integral = integral + 2.0d0 * y(i)
            end do
            integral = integral * (dx / 2.0d0)
        end if
    end function trapz

end module utils
