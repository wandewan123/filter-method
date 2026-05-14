module utils
    implicit none
contains

    function normalize(phi, dx) result(phi_out)
    
        real(8), intent(in) :: phi(:), dx
        real(8)             :: phi_out(size(phi))
        real(8)             :: integral
        integer             :: i, n

        n = size(phi)
        integral = 0.0d0

        do i = 1, n-1
            integral = integral + 0.5d0 * dx * (phi(i)**2 + phi(i+1)**2)
        end do

        if (integral > 0.0d0) then
            phi_out = phi / sqrt(integral)
        else
            phi_out = phi
        end if
    end function normalize

end module utils
