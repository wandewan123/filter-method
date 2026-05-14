module solver_module
    implicit none
contains

    subroutine tridag(a, b, c, r, n, u)

    !+──────────────────────────────────────────────────────────────────────────+
    !| SUBROUTINE : TRIDAG                                                       |
    !| TUJUAN     : DIGUNAKAN UNTUK MENYELESAIKAN SISTEM MATRIKS BERBENTUK       |
    !|              KOEFISIEN TRIDIAGONAL                                        |
    !| INPUT      : a (low-diag), b (diag), c(up-diag), r(RHS), n(grid)
    !| OUPUTS     : u(solution)
    !+──────────────────────────────────────────────────────────────────────────+

        integer, intent(in)   :: n
        real(8), intent(in)   :: a(n), b(n), c(n), r(n)
        real(8), intent(out)  :: u(n)
        real(8), dimension(n) :: gam, b_temp
        real(8) :: bet
        integer :: j

        b_temp  = b
        bet     = b_temp(1)
        
        if (bet == 0.0d0) return

        u(1) = r(1) / bet
        do j = 2, n
            gam(j) = c(j - 1) / bet
            bet    = b_temp(j) - (a(j) * gam(j))
            u(j)   = (r(j) - a(j) * u(j - 1)) / bet
        end do

        do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
        end do
    end subroutine tridag
        
end module solver_module
