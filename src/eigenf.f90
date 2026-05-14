module eigenf
  use matrix, only  : tridag, pentadag, septadag
  implicit none

contains

  subroutine eigenf2(dx, vpot, psi_in, n, E_curr, phi_out)
    !+──────────────────────────────────────────────────────────────────────────+
    !| SUBROUTINE : PSI_EIGEN_2                                                       |
    !| TUJUAN     : DIGUNAKAN UNTUK MENYELESAIKAN SISTEM MATRIKS BERBENTUK       |
    !|              KOEFISIEN TRIDIAGONAL                                        |
    !| INPUT      : a (low-diag), b (diag), c(up-diag), r(RHS), n(grid)
    !| OUPUTS     : u(solution)
    !+──────────────────────────────────────────────────────────────────────────+

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
