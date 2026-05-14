module utils
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

    ! Normalisasi Fungsi Gelombang dengan Aturan Trapezoid
    function normalize(phi, dx) result(phi_norm)
        real(8), intent(in) :: phi(:)
        real(8), intent(in) :: dx
        real(8)             :: phi_norm(size(phi))
        real(8)             :: integral, norm
        integer             :: n, i

        n = size(phi)
        integral = 0.0d0

        ! Rumus Trapezoid: dx/2 * (f(1)^2 + 2*f(2)^2 + ... + 2*f(n-1)^2 + f(n)^2)
        if (n > 1) then
            integral = phi(1)**2 + phi(n)**2
            do i = 2, n - 1
                integral = integral + 2.0d0 * phi(i)**2
            end do
            integral = integral * (dx / 2.0d0)
        end if

        norm = sqrt(integral)

        ! Hindari pembagian dengan nol jika fungsi kosong
        if (norm > 1.0d-15) then
            phi_norm = phi / norm
        else
            phi_norm = phi
        end if
    end function normalize

end module utils
