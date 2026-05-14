module eigenf
  use solver, only   : tridag, pentadag, septadag
  use utils, only    : normalize
  implicit none

contains

  subroutine eigenf2(dx, vpot, psi, n, energy, phi)
    !+──────────────────────────────────────────────────────────────────────────+
    !| SUBROUTINE : PSI_EIGEN_2                                                       |
    !| TUJUAN     : DIGUNAKAN UNTUK MENYELESAIKAN SISTEM MATRIKS BERBENTUK       |
    !|              KOEFISIEN TRIDIAGONAL                                        |
    !| INPUT      : a (low-diag), b (diag), c(up-diag), r(RHS), n(grid)
    !| OUPUTS     : u(solution)
    !+──────────────────────────────────────────────────────────────────────────+

        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, energy, vpot(n), psi(n)
        real(8), intent(out) :: phi(n)

        real(8)              :: dx2, integral
        integer              :: i
        real(8), allocatable :: alfa(:), beta(:)
        allocate(alfa(n), beta(n))
        
        dx2 = dx * dx
        alfa = -0.5d0 / dx2
        
        do i = 1, n
            beta(i) = (1.0d0 / dx2) + vpot(i) - energy
        end do

        call tridag(alfa, beta, alfa, psi, n, phi)

        phi = normalize(phi, dx)

    end subroutine eigenf2

    subroutine eigenf4(dx, vpot, psi, n, energy, phi)
    !+──────────────────────────────────────────────────────────────────────────+
    !| SUBROUTINE : EIGENF4                                                     |
    !| TUJUAN     : MENYELESAIKAN SISTEM MATRIKS BERBENTUK KOEFISIEN            |
    !|              PENTADIAGONAL (5-STENCIL, LAPLACIAN ORDE-4)                  |
    !+──────────────────────────────────────────────────────────────────────────+
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, energy, vpot(n), psi(n)
        real(8), intent(out) :: phi(n)
        
        real(8) :: a(n), b(n), c(n), d(n), e(n), dx2, integral
        integer :: i
        
        dx2 = dx * dx
        
        a = 1.0d0 / (24.0d0 * dx2)               
        b = -2.0d0 / (3.0d0 * dx2)        
        d = -2.0d0 / (3.0d0 * dx2)        
        e = 1.0d0 / (24.0d0 * dx2)        
        
        do i = 1, n
            c(i) = (5.0d0 / (4.0d0 * dx2)) + vpot(i) - energy
        end do

        call pentadag(a, b, c, d, e, psi, n, phi)

        phi = normalize(phi, dx)

    end subroutine eigenf4

    subroutine eigenf6(dx, vpot, psi, n, energy, phi)
    !+──────────────────────────────────────────────────────────────────────────+
    !| SUBROUTINE : EIGENF6                                                     |
    !| TUJUAN     : MENYELESAIKAN SISTEM MATRIKS BERBENTUK KOEFISIEN            |
    !|              HEPTADIAGONAL (7-STENCIL, LAPLACIAN ORDE-6)                  |
    !+──────────────────────────────────────────────────────────────────────────+
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, energy, vpot(n), psi(n)
        real(8), intent(out) :: phi(n)
        
        real(8) :: a3(n), a2(n), a1(n), d0(n), c1(n), c2(n), c3(n)
        real(8) :: dx2, integral
        integer :: i
        
        dx2 = dx * dx
        
        a3 = -1.0d0 / (180.0d0 * dx2)             
        a2 = 3.0d0 / (40.0d0 * dx2)               
        a1 = -3.0d0 / (4.0d0 * dx2)              
        c1 = -3.0d0 / (4.0d0 * dx2)                
        c2 = 3.0d0 / (40.0d0 * dx2)                
        c3 = -1.0d0 / (180.0d0 * dx2)              
        
        do i = 1, n
            d0(i) = (49.0d0 / (36.0d0 * dx2)) + vpot(i) - energy
        end do

        call septadag(a3, a2, a1, d0, c1, c2, c3, psi, n, phi)

        phi = normalize(phi, dx)

    end subroutine eigenf
