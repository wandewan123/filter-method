module main_filter
    use eigenf, only        : eigenf2, eigenf4, eigenf6
    use eigene_mod, only    : eigene
    use matrix, only        : matrix2, matrix4, matrix6
    implicit none

contains

    subroutine solver_filter(x, vpot, guessE, psi0, xi, stencil, tol, maxiter, energy_final, phi_final, status)
    
        real(8), intent(in)    :: x(:), vpot(:), psi0(:)
        real(8), intent(in)    :: guessE, xi, tol
        integer, intent(in)    :: stencil, maxiter
        
        real(8), intent(out)   :: energy_final, phi_final(size(x))
        integer, intent(out)    :: status
        
        integer                 :: n, iterasi, alloc_stat
        real(8)                 :: E_curr, E_next, dx, residual
        real(8), allocatable    :: psi_curr(:), phi(:), Hphi(:), r(:)
    
        n  = size(x)
        dx = x(2) - x(1)
        allocate(psi_curr(n), phi(n), Hphi(n), r(n), stat = alloc_stat)
        if (alloc_stat /= 0) then
            print *, "ERROR: Gagal mengalokasikan memori array!"
            status = -1
            return
        end if
        
        psi_curr = psi0
        E_curr   = guessE
        residual = 1.0d0
        iterasi  = 0
    
        do while (residual > tol .and. iterasi < maxiter)
            iterasi = iterasi + 1
    
            select case(stencil)
            case(3)
                call eigenf2(dx, vpot, psi_curr, n, E_curr, phi)
                call eigene(dx, vpot, phi, n, stencil, E_next)
                call matrix2(dx, vpot, phi, n, Hphi)
            case(5)
                call eigenf4(dx, vpot, psi_curr, n, E_curr, phi)
                call eigene(dx, vpot, phi, n, stencil, E_next)
                call matrix4(dx, vpot, phi, n, Hphi)
            case(7)
                call eigenf6(dx, vpot, psi_curr, n, E_curr, phi)
                call eigene(dx, vpot, phi, n, stencil, E_next)
                call matrix6(dx, vpot, phi, n, Hphi) 
            case default
                print *, "Error: Stencil harus bernilai 3, 5, atau 7."
                status = -2
                deallocate(psi_curr, phi, Hphi, r)
                return
            end select
    
            r = Hphi - (E_next * phi)
            residual = sqrt(sum(r ** 2) * dx) / max(abs(E_next), 1.0d-12)
            E_curr   = ((1.0d0 - xi) * E_next) + (xi * E_curr)
            psi_curr = phi
    
        end do
    
        energy_final = E_next
        phi_final    = phi
        if (residual <= tol) then
            status = 0
        else
            status = 1
        end if
        deallocate(psi_curr, phi, Hphi, r)
    
    end subroutine solver_filter

end module main_filter
