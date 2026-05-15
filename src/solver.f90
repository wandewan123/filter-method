module solver
    implicit none
    private
    public :: tridag, pentadag, septadag

    real(8), parameter :: EPS = 1.0d-20

contains

    subroutine tridag(a, b, c, r, n, u)
        integer, intent(in)   :: n
        real(8), intent(in)   :: a(n), b(n), c(n), r(n)
        real(8), intent(out)  :: u(n)
        
        real(8), dimension(n) :: gam
        real(8)               :: bet
        integer               :: j

        bet = b(1)
        if (abs(bet) < EPS) then
            print *, "Error [tridag]: Matriks singular atau diagonal nol!"
            stop
        end if

        u(1) = r(1) / bet
        do j = 2, n
            gam(j) = c(j - 1) / bet
            bet    = b(j) - (a(j) * gam(j))
            if (abs(bet) < EPS) then
                print *, "Error [tridag]: Pembagian dengan nol di indeks: ", j
                stop
            end if
            u(j)   = (r(j) - a(j) * u(j - 1)) / bet
        end do

        do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
        end do
    end subroutine tridag

    subroutine pentadag(a, b, c, d, e, r, n, u)
        integer, intent(in)   :: n
        real(8), intent(in)   :: a(n), b(n), c(n), d(n), e(n), r(n)
        real(8), intent(out)  :: u(n)
        
        real(8), dimension(n) :: c_m, d_m, e_m, r_m
        real(8)               :: xmu
        integer               :: i
    
        c_m(1) = c(1)
        if (abs(c_m(1)) < EPS) stop "Error [pentadag]: Pivot 1 nol!"
        
        d_m(1) = d(1) / c_m(1)
        e_m(1) = e(1) / c_m(1)
        r_m(1) = r(1) / c_m(1)
    
        c_m(2) = c(2) - b(2)*d_m(1)
        if (abs(c_m(2)) < EPS) stop "Error [pentadag]: Pivot 2 nol!"
        
        d_m(2) = (d(2) - b(2)*e_m(1)) / c_m(2)
        e_m(2) = e(2) / c_m(2)
        r_m(2) = (r(2) - b(2)*r_m(1)) / c_m(2)
    
        do i = 3, n
            xmu    = b(i) - a(i)*d_m(i-2)
            c_m(i) = c(i) - a(i)*e_m(i-2) - xmu*d_m(i-1)
            
            if (abs(c_m(i)) < EPS) stop "Error [pentadag]: Pivot nol di loop maju!"
            
            if (i < n)   d_m(i) = (d(i) - xmu*e_m(i-1)) / c_m(i)
            if (i < n-1) e_m(i) = e(i) / c_m(i)
            
            r_m(i) = (r(i) - a(i)*r_m(i-2) - xmu*r_m(i-1)) / c_m(i)
        end do
    
        u(n)   = r_m(n)
        u(n-1) = r_m(n-1) - d_m(n-1)*u(n)
        
        do i = n-2, 1, -1
            u(i) = r_m(i) - d_m(i)*u(i+1) - e_m(i)*u(i+2)
        end do
    end subroutine pentadag

    subroutine septadag(a3, a2, a1, d0, c1, c2, c3, r, n, u)
        integer, intent(in)   :: n
        real(8), intent(in)   :: a3(n), a2(n), a1(n), d0(n), c1(n), c2(n), c3(n), r(n)
        real(8), intent(out)  :: u(n)
        
        ! Salin semua komponen ke array kerja internal agar aman dieliminasi
        real(8), dimension(n) :: d0_m, c1_m, c2_m, c3_m, r_m
        real(8), dimension(n) :: a1_m, a2_m, a3_m
        real(8)               :: factor
        integer               :: i
    
        d0_m = d0; c1_m = c1; c2_m = c2; c3_m = c3; r_m = r
        a1_m = a1; a2_m = a2; a3_m = a3
    
        do i = 1, n
            if (abs(d0_m(i)) < EPS) stop "Error [septadag]: Pivot nol terdeteksi!"
            
            if (i + 1 <= n) then
                factor    = a1_m(i+1) / d0_m(i)
                d0_m(i+1) = d0_m(i+1) - factor * c1_m(i)
                c1_m(i+1) = c1_m(i+1) - factor * c2_m(i)
                c2_m(i+1) = c2_m(i+1) - factor * c3_m(i)
                r_m(i+1)  = r_m(i+1)  - factor * r_m(i)
            end if
            
            if (i + 2 <= n) then
                factor    = a2_m(i+2) / d0_m(i)
                a1_m(i+2) = a1_m(i+2) - factor * c1_m(i) ! MEMPERBAIKI BUG ELIMINASI
                d0_m(i+2) = d0_m(i+2) - factor * c2_m(i)
                c1_m(i+2) = c1_m(i+2) - factor * c3_m(i)
                r_m(i+2)  = r_m(i+2)  - factor * r_m(i)
            end if
            
            if (i + 3 <= n) then
                factor    = a3_m(i+3) / d0_m(i)
                a2_m(i+3) = a2_m(i+3) - factor * c1_m(i) ! MEMPERBAIKI BUG ELIMINASI
                a1_m(i+3) = a1_m(i+3) - factor * c2_m(i) ! MEMPERBAIKI BUG ELIMINASI
                d0_m(i+3) = d0_m(i+3) - factor * c3_m(i)
                r_m(i+3)  = r_m(i+3)  - factor * r_m(i)
            end if
        end do
    
        u(n) = r_m(n) / d0_m(n)
        if (n > 1) u(n-1) = (r_m(n-1) - c1_m(n-1)*u(n)) / d0_m(n-1)
        if (n > 2) u(n-2) = (r_m(n-2) - c1_m(n-2)*u(n-1) - c2_m(n-2)*u(n)) / d0_m(n-2)
        
        do i = n-3, 1, -1
            u(i) = (r_m(i) - c1_m(i)*u(i+1) - c2_m(i)*u(i+2) - c3_m(i)*u(i+3)) / d0_m(i)
        end do
    end subroutine septadag

end module solver
