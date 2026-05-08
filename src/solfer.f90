module solver_module
    implicit none
contains
    subroutine tridag(a, b, c, r, n, u)
        integer, intent(in)   :: n
        real(8), intent(in)   :: a(n), b(n), c(n), r(n)
        real(8), intent(out)  :: u(n)
        real(8), dimension(n) :: gam, b_temp
        real(8) :: bet
        integer :: j

        b_temp = b
        bet = b_temp(1)
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

    subroutine psi_eigen_fortran(dx, vpot, psi_in, energy, n, psi_out)
    ! Mendefinisikan tipe data dan arah variabel (intent)
    integer, intent(in)  :: n
    real(8), intent(in)  :: dx, energy, vpot(n), psi_in(n)
    real(8), intent(out) :: psi_out(n)
    
    ! Variabel lokal untuk perhitungan internal
    real(8) :: alfa(n), beta(n), dx2, integral
    integer :: i

    dx2 = dx * dx
    
    ! 1. Mengisi array alfa dan beta (Logika Hamiltonia)
    ! Di Fortran, operasi array bisa langsung (vektorisasi)
    alfa = -0.5d0 / dx2
    
    do i = 1, n
        beta(i) = (1.0d0 / dx2) + vpot(i) - energy
    end do

    ! 2. Memanggil subroutine tridag yang sudah kita buat sebelumnya
    ! Kita pakai alfa sebagai sub-diagonal dan super-diagonal
    call tridag(alfa, beta, alfa, psi_in, n, psi_out)

    ! 3. Normalisasi menggunakan aturan trapesium
    integral = 0.0d0
    do i = 1, n-1
        integral = integral + 0.5d0 * dx * (psi_out(i)**2 + psi_out(i+1)**2)
    end do
    
    ! Menghindari pembagian dengan nol jika terjadi error numerik
    if (integral > 0.0d0) then
        psi_out = psi_out / sqrt(integral)
    end if

end subroutine psi_eigen_fortran

        
        
end module solver_module
