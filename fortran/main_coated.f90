program main_coated
    use angular_functions_mod
    use bessel_mod
    use iso_fortran_env, wp => real128  ! quad  
    implicit none
    
    real(16), parameter :: pi_c=4.0_wp*atan(1.0_wp)
    real(16), parameter   :: floorval=1.0e-15_wp
    real(16), allocatable :: theta(:), S1a(:), S1pa(:), XV(:)
    real(16), allocatable :: pi(:,:), tau(:,:)
    complex(16) :: y, m1, m2, alpha, beta
    complex(16), allocatable :: D1x(:), D1y(:), D3x(:), D3y(:), A13x(:), lnA13y(:)
    integer :: n, Nmx, i, nang, Pc, l, layers
    real(16) :: wavelength, k, x, en, ii, check(2)
    complex(16), allocatable :: S1(:), S2(:), S1p(:), S2p(:), Q(:,:), M(:), Qp(:,:)
    complex(16) :: u0, u11, u13, u31, u33, A13y, umR121, umR212, a1, b1, R121(2), R212(2), T1(2)
    complex(16) :: a1p, b1p, ray(2)


    layers = 100
    nang = 1801

    allocate(M(layers+1))
    
    ! M = 1.33_wp
    ! M(1) = 1.35_wp
    ! do i = 1,layers
    !     ii = real(i,kind=wp)
    !     M(i) = 1.45_wp -1.05e-9*(i)**(4.0_wp) + 3.23e-7*(i)**(3.0_wp) - 3.14e-05*(i)**(2.0_wp) - 7.29e-05*(i)
    ! end do
    M(size(M)) = 1.0_wp
    

    allocate(XV(layers))
    do i = 1,layers
        XV(i) = 1.0_wp*i
    end do

    x = XV(size(XV))
    Nmx = int( x + 4.05_wp*(x**(1.0_wp/3.0_wp)) + 2.0_wp )
    wavelength = 532*10**(-9.0_wp)
    k = 2*pi_c / wavelength

    Pc = 1
    allocate(theta(nang), S1(nang), S1a(nang), S1p(nang), S1pa(nang))
    do i = 1, nang
      theta(i) = real((i - 1),kind=wp) * pi_c / real(nang-1,kind=wp)
    end do

    call angularFunctions(theta, Nmx, pi, tau)

    allocate(D1x(0:Nmx), D1y(0:Nmx), D3x(0:Nmx), D3y(0:Nmx), A13x(0:Nmx), lnA13y(0:Nmx))
    

    allocate(Q(size(XV), 2), Qp(size(XV), 2))

    do n  = 1, Nmx
        a1p = 1.0_wp
        b1p = 1.0_wp
        en = real(n, kind=wp)
        Q = 0.0_wp
        do l = 1,size(XV)
            m1 = M(l)
            m2 = M(l+1)
            x = XV(l)
            y = m1/m2 * x

            call computeBessel(x, y, Nmx, D1x, D1y, D3x, D3y, A13x, lnA13y)

            do i = 1,2
                if (i.eq.1) then
                    alpha = m1 /m2
                    beta = 1.0_wp
                else
                    alpha = 1.0_wp
                    beta = m1 /m2
                end if

                u0 = (D1x(n) - D3x(n)) * (D1y(n) - D3y(n))
                u11 = alpha*D1x(n) - beta*D1y(n)
                u31 = alpha*D3x(n) - beta*D1y(n)
                u13 = alpha*D1x(n) - beta*D3y(n)
                u33 = alpha*D3x(n) - beta*D3y(n)

                if (real(lnA13y(n), kind=wp) < 0.0_wp .or. aimag(m1) == 0.0_wp) then
                    A13y = exp(lnA13y(n))
                    T1(i) = m1/m2 * A13x(n) * A13y * u0 / (u33 - A13y*u31)**2.0_wp

                    umR212 = A13x(n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    umR121 = -A13y * u31 / (u33 - A13y*u31)

                    R212(i) = 1.0_wp - A13x(n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    R121(i) = 1.0_wp + A13y * u31 / (u33 - A13y*u31)
                else
                    A13y = exp(-lnA13y(n))
                    T1(i) = m1/m2 * A13x(n) * A13y * u0 / (u33 - A13y*u31)**2.0_wp

                    umR212 = A13x(n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    umR121 = -A13y * u31 / (u33 - A13y*u31)

                    R212(i) = 1.0_wp - A13x(n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    R121(i) = 1.0_wp + A13y * u31 / (u33 - A13y*u31)
                end if
            end do  

            if (l.eq.1) then
                ray = T1 * (R121)**real(Pc-1,kind=wp)

                Q(l,:) = R212 +  T1 / (1.0_wp - R121)
                Qp(l,:) = ray
            else
                ray = (R121 * Q(l-1,:))**real(Pc-1,kind=wp)

                Q(l,:) = R212 +  (T1 * Q(l-1,:)) / (1.0_wp - R121 * Q(l-1,:))
                Qp(l,:) = ray * T1 * Q(l-1,:)
            end if
        end do  

        if (Pc.eq.0) then
            a1p = 0.5_wp - 0.5_wp*R212(1)
            b1p = 0.5_wp - 0.5_wp*R212(2)
        else
            a1p =  -0.5_wp * Qp(size(XV),1)
            b1p =  -0.5_wp * Qp(size(XV),2)
        end if

        a1 = 0.5_wp * ( 1.0_wp - Q(size(XV),1) )
        b1 = 0.5_wp * ( 1.0_wp - Q(size(XV),2) )

        S1 = S1 + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1*pi(n,:) + b1*tau(n,:))
        S2 = S1 + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1*tau(n,:) + b1*pi(n,:))

        S1p = S1p + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1p*pi(n,:) + b1p*tau(n,:))
        ! S2p = S2p + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1p*tau(n,:) + b1p*pi(n,:))

    end do

    S1a = abs(S1)**2.0_wp
    S1pa = abs(S1p)**2.0_wp

    write(*,*) "writing... output.dat"
    open(10,file="output.dat",status="unknown")
    write(10,*)
    write(10,*) "ang (deg)"," S1", "S1(P=Pc)"

    do i=1,nang
      write(10,*) theta(i),S1a(i),S1pa(i)
    enddo
    write(10,*)
    close(10)
    
    ! deallocate(theta, pi, tau, D1x, D1y, D3x, D3y, A13x, lnA13y, S1, ray)

    stop

end program main_coated