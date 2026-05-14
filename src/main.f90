module main_filter
    use solver, only    : tridag, pentadag, septadag
    use eigenf, only    : eigenf2, eigenf4, eigenf6
    use eigene, only    : eigene2, eigene4, eigene6
    use matrix, only    : matrix2, matrix4, matrix6
    implicit none

contains

    subroutine main(x, vpot, guessE, psi0, xi, stencil, tol, maxiter, energy_final, phi_final)
        real(8), intent(in)    :: x(:), vpot(:)
        real(8), intent(in)    :: guessE, xi, tol
        real(8), intent(in)    :: psi0(:)
        integer, intent(in)    :: stencil, maxiter
        real(8), intent(out)   :: energi, phi(size(x))
    
        integer                 :: n, iterasi, i
        real(8)                 :: E_curr, E_next, dx2, residual, norm_res
        real(8)                 :: r(n), r_sq(n), integral_residu
        real(8), allocatable    :: psi_curr(:), phi(:), Hphi(:)
    
        n   = size(x)
        dx2 = (x(2) - x(1)) ** 2
    
        allocate(psi_curr(n), phi(n), H_phi(n))
    
        psi_curr = psi0
        E_curr   = guessE
        iterasi  = 0
        residual = 1.0d0
    
        do while (residual > tol)
            iterasi = iterasi + 1
    
            select case(stencil)
            case(3)
                call eigenf2(dx, vpot, psi, n, E_curr, phi)
                call eigene2(dx, vpot, psi, n, E_next)
                call matrix2(Hphi)
            case(5)
                call eigenf4(dx, vpot, psi, n, E_curr, phi)
                call eigene4(dx, vpot, psi, n, E_next)
                call matrix4(Hphi)
            case(7)
                call eigenf6(dx, vpot, psi, n, E_curr, phi)
                call eigene6(dx, vpot, psi, n, E_next)
                call matrix6(Hphi)
    
            r        = Hphi - (E_next * phi)
            residual = sqrt(sum(r ** 2) * dx) / abs(E_next)
    
            E_curr = ((1.0d0 - xi) * E_next) + (xi * E_current)
            psi    = phi
        
        end do
    
        energy_final = E_next
        phi_final    = psi
        deallocate(psi_curr, phi, H_phi)
        
    end subroutine main_filter

end module filter_interface

    subroutine energi_eigen_fortran(dx, vpot, psi, n, energy_out)
        integer, intent(in) :: n
        real(8), intent(in) :: dx, vpot(n), psi(n)
        real(8), intent(out) :: energy_out
        real(8) :: u(n), psi_u(n), dx2, alfa, seperdx2
        integer :: i

        dx2 = dx*dx
        alfa = -0.5d0/dx2
        seperdx2 = 1.0d0/dx2

        ! Menghitung H * psi (Hamiltonian dikali fungsi gelombang)
        u(1) = (seperdx2 + vpot(1))*psi(1) + alfa*psi(2)
        do i = 2, n-1
            u(i) = (alfa*psi(i-1)) + ((seperdx2 + vpot(i))*psi(i)) + (alfa*psi(i+1))
        end do
        u(n) = (alfa*psi(n-1)) + ((seperdx2 + vpot(n))*psi(n))

        ! Integrasi <psi|H|psi> menggunakan aturan trapesium
        psi_u = psi * u
        energy_out = 0.0d0
        do i = 1, n-1
            energy_out = energy_out + 0.5d0 * dx * (psi_u(i) + psi_u(i+1))
        end do
    end subroutine energi_eigen_fortran

    subroutine main_filter(x, n, vpot, initial_e, dx, xi, tol, energy_final, psi_final)
        integer, intent(in) :: n
        real(8), intent(in) :: x(n), vpot(n), initial_e, dx, xi, tol
        real(8), intent(out) :: energy_final, psi_final(n)
        
        real(8) :: energy, energy_baru, err_e, psi(n), psi_next(n)
        integer :: iterasi

        ! Inisialisasi tebakan awal fungsi gelombang (Gaussian)
        psi = exp(-(x**2)) * 0.25d0
        energy = initial_e
        err_e = 1.0d0
        iterasi = 0

        ! Loop Iterasi (While loop di Python pindah ke sini)
        do while (err_e > tol .and. iterasi < 50000)
            iterasi = iterasi + 1
            
            ! 1. Update Fungsi Gelombang
            call psi_eigen_fortran(dx, vpot, psi, energy, n, psi_next)
            psi = psi_next

            ! 2. Update Energi
            call energi_eigen_fortran(dx, vpot, psi, n, energy_baru)
            
            ! 3. Hitung Error dan Update Energi dengan Mixing (xi)
            err_e = abs((energy_baru - energy) / energy)
            energy = ((1.0d0 - xi) * energy_baru) + (xi * energy)
        end do

        energy_final = energy_baru
        psi_final = psi
    end subroutine main_filter
