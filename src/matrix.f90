module matrix
    use utils, only: get_phi
    implicit none

contains

! ==========================================================================
! MATRIX2 (Orde-2, 3-point stencil)
! ==========================================================================
subroutine matrix2(dx, vpot, phi, n, Hphi, status)
    integer, intent(in)  :: n
    real(8), intent(in)  :: dx, vpot(n), phi(n)
    real(8), intent(out) :: Hphi(n)
    integer, intent(out) :: status

    real(8) :: dx2, invdx2, alfa
    integer :: i

    status = 0
    if (n < 3) then
        print *, "ERROR [matrix2]: n terlalu kecil untuk stencil 3"
        status = -4 
        return
    end if

    dx2    = dx*dx
    invdx2 = 1.0d0 / dx2
    alfa   = -0.5d0 * invdx2

    do i = 1, n
        Hphi(i) = alfa * get_phi(phi, i-1, n) + &
                  (invdx2 + vpot(i)) * phi(i) + &
                  alfa * get_phi(phi, i+1, n)
    end do

end subroutine matrix2

! ==========================================================================
! MATRIX4 (Orde-4, 5-point stencil)
! ==========================================================================
subroutine matrix4(dx, vpot, phi, n, Hphi, status)
    integer, intent(in)  :: n
    real(8), intent(in)  :: dx, vpot(n), phi(n)
    real(8), intent(out) :: Hphi(n)
    integer, intent(out) :: status
    

    real(8) :: dx2, c0, c1, c2
    integer :: i

    status = 0
    if (n < 5) then
        print *, "ERROR [matrix4]: n terlalu kecil untuk stencil 5"
        status = -4
        return
    end if

    dx2 = dx*dx
    c2 =  1.0d0 / (24.0d0 * dx2)
    c1 = -2.0d0 / (3.0d0 * dx2)
    c0 =  5.0d0 / (4.0d0 * dx2)

    do i = 1, n
        Hphi(i) = c2 * get_phi(phi, i-2, n) + &
                  c1 * get_phi(phi, i-1, n) + &
                  (c0 + vpot(i)) * phi(i)   + &
                  c1 * get_phi(phi, i+1, n) + &
                  c2 * get_phi(phi, i+2, n)
    end do

end subroutine matrix4

! ==========================================================================
! MATRIX6 (Orde-6, 7-point stencil)
! ==========================================================================
subroutine matrix6(dx, vpot, phi, n, Hphi, status)
    integer, intent(in)  :: n
    real(8), intent(in)  :: dx, vpot(n), phi(n)
    real(8), intent(out) :: Hphi(n)
    integer, intent(out) :: status

    real(8) :: dx2, c0, c1, c2, c3
    integer :: i

    status = 0
    if (n < 7) then
        print *, "ERROR [matrix6]: n terlalu kecil untuk stencil 7"
        status = -4
        return
    end if

    dx2 = dx*dx
    c3 = -1.0d0 / (180.0d0 * dx2)
    c2 =  3.0d0 / (40.0d0 * dx2)
    c1 = -3.0d0 / (4.0d0 * dx2)
    c0 =  49.0d0 / (36.0d0 * dx2)

    do i = 1, n
        Hphi(i) = c3 * get_phi(phi, i-3, n) + &
                  c2 * get_phi(phi, i-2, n) + &
                  c1 * get_phi(phi, i-1, n) + &
                  (c0 + vpot(i)) * phi(i)   + &
                  c1 * get_phi(phi, i+1, n) + &
                  c2 * get_phi(phi, i+2, n) + &
                  c3 * get_phi(phi, i+3, n)
    end do

end subroutine matrix6

end module matrix
