program main_coated
    use angular_functions_mod
    use bessel_mod
    use iso_fortran_env, wp => real128  ! quad  
    implicit none
    
    real(16), parameter :: pi_c=4.0_wp*atan(1.0_wp)
    real(16), parameter   :: floorval=1.0e-50_wp
    real(16), allocatable :: theta(:), S1a(:)
    real(16), allocatable :: pi(:,:), tau(:,:)
    complex(16) :: y, m1, m2, alpha, beta
    complex(16), allocatable :: D1x(:), D1y(:), D3x(:), D3y(:), A13x(:), lnA13y(:)
    integer :: n, Nmx, i, nang, P, Pc, pee, l
    real(16) :: wavelength, k, x, en, XV(6), check
    complex(16), allocatable :: S1(:), S2(:), ray(:,:), Q(:,:)
    complex(16) :: u0, u11, u13, u31, u33, A13y, umR121, umR212, a1, b1, R121(2), R212(2), T1(2)
    complex(16) :: a1p, b1p, M(7)
    real(16) :: start, finish

    nang = 1801
    M(1) = (1.35, 0)
    M(2) = (1.3, 0)
    M(3) = (1.25, 0)
    M(4) = (1.2, 0)
    M(5) = (1.15, 0)
    M(6) = (1.1, 0)
    M(7) = (1.0, 0)

    XV(1) = 100
    XV(2) = 200
    XV(3) = 300
    XV(4) = 400
    XV(5) = 500
    XV(6) = 600

    x = XV(size(XV))
    Nmx = int( x + 4.05_wp*(x**(1.0_wp/3.0_wp)) + 2.0_wp )
    wavelength = 532*10**(-9.0_wp)
    k = 2*pi_c / wavelength

    P = 200
    Pc = 0
    allocate(theta(nang), S1(nang), S1a(nang))
    do i = 1, nang
      theta(i) = real((i - 1),kind=wp) * pi_c / real(nang-1,kind=wp)
    end do

    call angularFunctions(theta, Nmx, pi, tau)

    allocate(D1x(0:Nmx), D1y(0:Nmx), D3x(0:Nmx), D3y(0:Nmx), A13x(0:Nmx), lnA13y(0:Nmx))
    

    allocate(ray(2,P), Q(size(XV), 2))

    do n  = 1, Nmx
        a1p = 1.0_wp
        b1p = 1.0_wp
        en = real(n, kind=wp)
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

                if (l.eq.1) then
                    do pee = 1, P
                        ray(i, pee) = T1(i) * (R121(i))**real(pee-1,kind=wp)
                        check = abs(R121(i)**real(pee-1,kind=wp))

                        if  (check.le.floorval) then
                            ray(i,pee+1:P) = 0.0_wp
                            exit
                        endif
                    end do
                    Q(l,i) = R212(i) + sum(ray(i,:))
                else
                    do pee = 1, P
                        ray(i, pee) = (R121(i) * Q(l-1,i))**real(pee-1,kind=wp)
                        check = abs((R121(i) * Q(l-1,i))**real(pee-1,kind=wp))

                        if  (check.le.floorval) then
                            ray(i,pee+1:P) = 0.0_wp
                            exit
                        endif
                    end do
                    Q(l,i) = R212(i) + T1(i)*Q(l-1,i) * sum(ray(i,:))
                end if
            end do  
        end do  

        a1 = 0.5_wp * ( 1.0_wp - Q(size(XV),1) )
        b1 = 0.5_wp * ( 1.0_wp - Q(size(XV),2) )

        S1 = S1 + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1*pi(n,:) + b1*tau(n,:))
        S2 = S1 + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1*tau(n,:) + b1*pi(n,:))

    end do

    S1a = abs(S1)**2.0_wp

    write(*,*) "writing... output.dat"
    open(10,file="output.dat",status="unknown")
    write(10,*)
    write(10,*) "ang (deg)"," S1"

    do i=1,nang
      write(10,*) theta(i),S1a(i)
    enddo
    write(10,*)
    close(10)
    
    deallocate(theta, pi, tau, D1x, D1y, D3x, D3y, A13x, lnA13y, S1, S2, ray)

    stop

end program main_coated