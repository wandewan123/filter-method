module eigene
    use matrix, only: matrix2, matrix4, matrix6
    use utils,  only: trapz
    implicit none

contains

subroutine eigene(dx, vpot, psi, n, stencil, energy_out)
    integer, intent(in)  :: n, stencil
    real(8), intent(in)  :: dx, vpot(n), psi(n)
    real(8), intent(out) :: energy_out

    real(8), allocatable :: H_psi(:), psi_Hpsi(:)

    allocate(H_psi(n), psi_Hpsi(n))

    ! Pilih operator
    select case(stencil)
    case(3)
        call matrix2(dx, vpot, psi, n, H_psi)
    case(5)
        call matrix4(dx, vpot, psi, n, H_psi)
    case(7)
        call matrix6(dx, vpot, psi, n, H_psi)
    case default
        error stop "Stencil tidak valid"
    end select

    psi_Hpsi = psi * H_psi

    energy_out = trapz(psi_Hpsi, dx)

    deallocate(H_psi, psi_Hpsi)

end subroutine eigene

end module eigene
