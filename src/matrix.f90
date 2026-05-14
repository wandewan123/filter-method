module matrix
    implicit none

contains

    ! ==========================================================================
    ! 1. MATRIX2 (3-Stencil / Orde-2)
    ! ==========================================================================
    subroutine matrix2(dx, vpot, phi, n, Hphi)
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, vpot(n), phi(n)
        real(8), intent(out) :: Hphi(n)
        
        real(8) :: alfa, seperdx2
        integer :: i

        alfa = -0.5d0 / (dx*dx)
        seperdx2 = 1.0d0 / (dx*dx)

        ! Batas Kiri
        Hphi(1) = (seperdx2 + vpot(1)) * phi(1) + alfa * phi(2)
        
        ! Interior
        do i = 2, n-1
            Hphi(i) = (alfa * phi(i-1)) + ((seperdx2 + vpot(i)) * phi(i)) + (alfa * phi(i+1))
        end do
        
        ! Batas Kanan
        Hphi(n) = (alfa * phi(n-1)) + ((seperdx2 + vpot(n)) * phi(n))
    end subroutine matrix2

    ! ==========================================================================
    ! 2. MATRIX4 (5-Stencil / Orde-4)
    ! ==========================================================================
    subroutine matrix4(dx, vpot, phi, n, Hphi)
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, vpot(n), phi(n)
        real(8), intent(out) :: Hphi(n)
        
        real(8) :: c_0, c_1, c_2
        integer :: i

        c_2 = 1.0d0 / (24.0d0 * dx * dx)
        c_1 = -2.0d0 / (3.0d0 * dx * dx)
        c_0 = 5.0d0 / (4.0d0 * dx * dx)

        ! Penanganan Kondisi Batas Tepi Kiri
        Hphi(1) = (c_0 + vpot(1)) * phi(1) + c_1 * phi(2) + c_2 * phi(3)
        Hphi(2) = c_1 * phi(1) + (c_0 + vpot(2)) * phi(2) + c_1 * phi(3) + c_2 * phi(4)
        
        ! Area Interior Utama
        do i = 3, n-2
            Hphi(i) = c_2 * phi(i-2) + c_1 * phi(i-1) + &
                      (c_0 + vpot(i)) * phi(i) + &
                      c_1 * phi(i+1) + c_2 * phi(i+2)
        end do
        
        ! Penanganan Kondisi Batas Tepi Kanan
        Hphi(n-1) = c_2 * phi(n-3) + c_1 * phi(n-2) + (c_0 + vpot(n-1)) * phi(n-1) + c_1 * phi(n)
        Hphi(n)   = c_2 * phi(n-2) + c_1 * phi(n-1) + (c_0 + vpot(n)) * phi(n)
    end subroutine matrix4

    ! ==========================================================================
    ! 3. MATRIX6 (7-Stencil / Orde-6)
    ! ==========================================================================
    subroutine matrix6(dx, vpot, phi, n, Hphi)
        integer, intent(in)  :: n
        real(8), intent(in)  :: dx, vpot(n), phi(n)
        real(8), intent(out) :: Hphi(n)
        
        real(8) :: c_0, c_1, c_2, c_3
        integer :: i

        c_3 = -1.0d0 / (180.0d0 * dx * dx)
        c_2 = 3.0d0 / (40.0d0 * dx * dx)
        c_1 = -3.0d0 / (4.0d0 * dx * dx)
        c_0 = 49.0d0 / (36.0d0 * dx * dx)

        ! Penanganan Kondisi Batas Tepi Kiri
        Hphi(1) = (c_0 + vpot(1)) * phi(1) + c_1 * phi(2) + c_2 * phi(3) + c_3 * phi(4)
        Hphi(2) = c_1 * phi(1) + (c_0 + vpot(2)) * phi(2) + c_1 * phi(3) + c_2 * phi(4) + c_3 * phi(5)
        Hphi(3) = c_2 * phi(1) + c_1 * phi(2) + (c_0 + vpot(3)) * phi(3) + c_1 * phi(4) + c_2 * phi(5) + c_3 * phi(6)

        ! Area Interior Utama
        do i = 4, n-3
            Hphi(i) = c_3 * phi(i-3) + c_2 * phi(i-2) + c_1 * phi(i-1) + &
                      (c_0 + vpot(i)) * phi(i) + &
                      c_1 * phi(i+1) + c_2 * phi(i+2) + c_3 * phi(i+3)
        end do

        ! Penanganan Kondisi Batas Tepi Kanan
        Hphi(n-2) = c_3 * phi(n-5) + c_2 * phi(n-4) + c_1 * phi(n-3) + (c_0 + vpot(n-2)) * phi(n-2) + c_1 * phi(n-1) + c_2 * phi(n)
        Hphi(n-1) = c_3 * phi(n-4) + c_2 * phi(n-3) + c_1 * phi(n-2) + (c_0 + vpot(n-1)) * phi(n-1) + c_1 * phi(n)
        Hphi(n)   = c_3 * phi(n-3) + c_2 * phi(n-2) + c_1 * phi(n-1) + (c_0 + vpot(n)) * phi(n)
    end subroutine matrix6

end module matrix
