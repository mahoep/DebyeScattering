module bessel_mod
    implicit none
    
contains
    
    subroutine computeBessel(x, y, Nmx, D1x, D1y, D3x, D3y, A13x, lnA13y)
        use iso_fortran_env, wp => real128  ! quad  
        complex(16), intent(in) :: y
        real(16), intent(in) :: x
        integer, intent(in) :: Nmx
        complex(16), intent(out) :: D1x(0:Nmx), D1y(0:Nmx), D3x(0:Nmx), D3y(0:Nmx), A13x(0:Nmx), lnA13y(0:Nmx)
        integer :: n, j
        complex(16) :: nu
        complex(16), allocatable :: ajx(:), ajy(:), w31(:), D1x0, D1y0
        real(16) :: en
        
        nu = complex(Nmx + 0.5_wp, 0.0_wp)
        
        allocate(ajx(0:Nmx), ajy(0:Nmx), w31(0:Nmx))
        
        ajx(1) = 2.0_wp * nu / x
        ajy(1) = 2.0_wp * nu / y
        
        do j = 2, Nmx
            ajx(j) = (-1.0_wp)**(j-1) * 2.0_wp * (nu + real((j-1),kind=wp)) / x
            ajy(j) = (-1.0_wp)**(j-1) * 2.0_wp * (nu + real((j-1),kind=wp)) / y
        end do

        D1x0 = ajx(Nmx)
        D1y0 = ajy(Nmx)

        do j = 1, Nmx-1
            D1x0 = ajx(Nmx-j) + 1.0_wp / D1x0
            D1y0 = ajy(Nmx-j) + 1.0_wp / D1y0
        end do
        
        D1x(Nmx) = -real(Nmx,kind=wp) / x + D1x0
        D1y(Nmx) = -real(Nmx,kind=wp) / y + D1y0

        do n = 1, Nmx
            en = real(Nmx - n + 1,kind=wp)
            D1x(Nmx-n) = en / x - 1.0_wp / (D1x(Nmx-n+1) + en / x)
            D1y(Nmx-n) = en / y - 1.0_wp / (D1y(Nmx-n+1) + en / y)
        end do
        
        D3x(0) = complex(0.0_wp, 1.0_wp)
        D3y(0) = complex(0.0_wp, 1.0_wp)
        
        w31(0) = -1.0_wp / tan(y) + complex(0.0_wp, 1.0_wp)
        
        do n = 1, Nmx
            en = real(n, kind=wp)
            D3x(n) = -en / x + 1.0_wp / (en / x - D3x(n-1))
            w31(n) = w31(n-1) / ((D3y(n-1) - en / y) * (D1y(n-1) - en / y))
            D3y(n) = D1y(n) + w31(n)
        end do
        
        A13x(0) = 2.0_wp / (1.0_wp - complex(0.0_wp, 1.0_wp) / tan(x))
        lnA13y(0) = 2.0_wp * aimag(y) + log(exp(-2.0_wp * aimag(y)) - &
            cos(2.0_wp * real(y,kind=wp)) + complex(0.0_wp, 1.0_wp) * sin(2.0_wp * real(y,kind=wp)))
        
        do n = 1, Nmx
            A13x(n) = A13x(n-1) * (D1x(n-1) - real(n,kind=wp) / x) / (D3x(n-1) - real(n,kind=wp) / x)
            lnA13y(n) = lnA13y(n-1) + log((D1y(n-1) - real(n,kind=wp) / y) / (D3y(n-1) - real(n,kind=wp) / y))
        end do
        
        deallocate(ajx, ajy, w31)
    end subroutine computeBessel

end module bessel_mod
