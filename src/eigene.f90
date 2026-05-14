module eigene
    use matrix, only: matrix2, matrix4, matrix6
    implicit none

contains

    ! ==========================================================================
    ! 1. EIGENE2 (3-Stencil / Orde-2)
    ! ==========================================================================
    subroutine eigene2(dx, vpot, psi, n, energy_out)
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, vpot(n), psi(n)
        real(8), intent(out) :: energy_out
        
        real(8) :: H_psi(n), psi_Hpsi(n)
        integer :: i

        ! Menggunakan modul matrix untuk menghitung H * psi
        call matrix2(dx, vpot, psi, n, H_psi)

        ! Integrasi <psi|H|psi> menggunakan aturan trapesium
        psi_Hpsi = psi * H_psi
        energy_out = 0.0d0
        do i = 1, n-1
            energy_out = energy_out + 0.5d0 * dx * (psi_Hpsi(i) + psi_Hpsi(i+1))
        end do
    end subroutine eigene2

    ! ==========================================================================
    ! 2. EIGENE4 (5-Stencil / Orde-4)
    ! ==========================================================================
    subroutine eigene4(dx, vpot, psi, n, energy_out)
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, vpot(n), psi(n)
        real(8), intent(out) :: energy_out
        
        real(8) :: H_psi(n), psi_Hpsi(n)
        integer :: i

        ! Menggunakan modul matrix untuk menghitung H * psi
        call matrix4(dx, vpot, psi, n, H_psi)

        ! Integrasi <psi|H|psi> menggunakan aturan trapesium
        psi_Hpsi = psi * H_psi
        energy_out = 0.0d0
        do i = 1, n-1
            energy_out = energy_out + 0.5d0 * dx * (psi_Hpsi(i) + psi_Hpsi(i+1))
        end do
    end subroutine eigene4

    ! ==========================================================================
    ! 3. EIGENE6 (7-Stencil / Orde-6)
    ! ==========================================================================
    subroutine eigene6(dx, vpot, psi, n, energy_out)
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, vpot(n), psi(n)
        real(8), intent(out) :: energy_out
        
        real(8) :: H_psi(n), psi_Hpsi(n)
        integer :: i

        ! Menggunakan modul matrix untuk menghitung H * psi
        call matrix6(dx, vpot, psi, n, H_psi)

        ! Integrasi <psi|H|psi> menggunakan aturan trapesium
        psi_Hpsi = psi * H_psi
        energy_out = 0.0d0
        do i = 1, n-1
            energy_out = energy_out + 0.5d0 * dx * (psi_Hpsi(i) + psi_Hpsi(i+1))
        end do
    end subroutine eigene6

end module eigene
