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
    complex(16), allocatable, dimension(:,:) :: D1xL, D1yL, D3xL, D3yL, A13xL, lnA13yL
    integer :: n, Nmx, i, nang, Pc, l, layers
    real(16) :: wavelength, x, en, ii, X_max
    complex(16), allocatable :: S1(:), S2(:), S1p(:), S2p(:), Q(:,:), M(:), Qp(:,:)
    complex(16) :: u0, u11, u13, u31, u33, A13y, umR121, umR212, a1, b1, R121(2), R212(2), T1(2)
    complex(16) :: a1p, b1p, ray(2)
    real :: start, finish

    call cpu_time(start)
    layers = 100 ! specify number of layers, 1= homogeneous case
    nang = 1801 ! angular resolution from 0 to pi.
    X_max = 100.0 ! size parameter of outermost layer
    Pc = 1 ! choice of debye series term to write out, 0,1,2,3, etc

    ! allocate size parameter vector (XV) with following loop. Layers are equally spaced.
    allocate(XV(layers))
    do i = 1,layers
        ii = real(i,kind=wp)
        XV(i) = X_max / real(layers,kind=wp) * ii
    end do

    ! build vector for index of refraction 
    ! the loop has the the index change from 1.45 to 1.35 over a non dimensional radius, ii
    ! M(size(M)) = 1 specifies the index surrounding the sphere (air, vacuum, etc)
    allocate(M(layers+1))
    ! M = 1.45_wp
    ! M(1) = 1.35_wp
    do i = 1,layers
        ii = XV(i) / XV(size(XV)) 
        M(i) = 1.4529_wp + 0.1076_wp*((ii)**(3.0_wp)) - 0.1768_wp*((ii)**(2.0_wp)) - 0.0348_wp*(ii)
    end do
    M(size(M)) = 1.0_wp
    
    ! temporarily store the largest size in x and find the number of iterations to perform, Nmx
    x = XV(size(XV))
    Nmx = int( x + 4.05_wp*(x**(1.0_wp/3.0_wp)) + 2.0_wp )
    wavelength = 532e-9_wp

    ! build array to store angles in radians
    allocate(theta(nang), S1(nang), S1a(nang), S1p(nang), S1pa(nang))
    do i = 1, nang
      theta(i) = real((i - 1),kind=wp) * pi_c / real(nang-1,kind=wp)
    end do

    ! pre compute pi and tau (associated Legendre functions)
    call angularFunctions(theta, Nmx, pi, tau)

    allocate(D1x(0:Nmx), D1y(0:Nmx), D3x(0:Nmx), D3y(0:Nmx), A13x(0:Nmx), lnA13y(0:Nmx))
    allocate(D1xL(layers,0:Nmx), D1yL(layers,0:Nmx), D3xL(layers,0:Nmx), &
        D3yL(layers,0:Nmx), A13xL(layers,0:Nmx), lnA13yL(layers,0:Nmx))

    allocate(Q(size(XV), 2), Qp(size(XV), 2))

    ! pre compute the spherical bessel functions for all index of refraction and iterations, N
    ! with a large number of layers, and a large sphere this will take up a lot of memory
    do l = 1,size(XV)
        m1 = M(l)
        m2 = M(l+1)
        x = XV(l)
        y = m1/m2 * x

        call computeBessel(x, y, Nmx, D1x, D1y, D3x, D3y, A13x, lnA13y)
        D1xL(l,:) = D1x
        D1yL(l,:) = D1y
        D3xL(l,:) = D3x
        D3yL(l,:) = D3y
        A13xL(l,:) = A13x
        lnA13yL(l,:) = lnA13y

    end do

    deallocate(D1x, D1y, D3x, D3y, A13x, lnA13y)

    ! main loop
    do n  = 1, Nmx
        a1p = 1.0_wp
        b1p = 1.0_wp
        en = real(n, kind=wp)
        Q = 0.0_wp
        ! for every n, we must find the contribution from all layers
        do l = 1,size(XV)
            m1 = M(l)
            m2 = M(l+1)
            x = XV(l)
            y = m1/m2 * x

            ! we must find the scattering coefficents for a TE and TM wave (i=1,2)
            do i = 1,2
                if (i.eq.1) then
                    alpha = m1 /m2
                    beta = 1.0_wp
                else
                    alpha = 1.0_wp
                    beta = m1 /m2
                end if

                u0 = (D1xL(l,n) - D3xL(l,n)) * (D1yL(l,n) - D3yL(l,n))
                u11 = alpha*D1xL(l,n) - beta*D1yL(l,n)
                u31 = alpha*D3xL(l,n) - beta*D1yL(l,n)
                u13 = alpha*D1xL(l,n) - beta*D3yL(l,n)
                u33 = alpha*D3xL(l,n) - beta*D3yL(l,n)

                ! depending on the quadrant of the argument of lnA13y we may have to take the negative exponential (else)
                if (real(lnA13yL(l,n), kind=wp) < 0.0_wp .or. aimag(m1) == 0.0_wp) then
                    A13y = exp(lnA13yL(l,n))
                    T1(i) = m1/m2 * A13xL(l,n) * A13y * u0 / (u33 - A13y*u31)**2.0_wp

                    umR212 = A13xL(l,n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    umR121 = -A13y * u31 / (u33 - A13y*u31)

                    R212(i) = 1.0_wp - A13xL(l,n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    R121(i) = 1.0_wp + A13y * u31 / (u33 - A13y*u31)
                else
                    A13y = exp(-lnA13yL(l,n))
                    T1(i) = m1/m2 * A13xL(l,n) * A13y * u0 / (u33 - A13y*u31)**2.0_wp

                    umR212 = A13xL(l,n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    umR121 = -A13y * u31 / (u33 - A13y*u31)

                    R212(i) = 1.0_wp - A13xL(l,n) * (u13 - A13y*u11) / (u33 - A13y*u31)
                    R121(i) = 1.0_wp + A13y * u31 / (u33 - A13y*u31)
                end if
            end do  

            ! the first layer is treated the same a homogeneous sphere
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

        ! if Pc is 0 (diffraction and reflection) we only care about R212 rays
        ! otherwise we use the recurisve Q, matrix to find a,b
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
        S2p = S2p + (2.0_wp*en+1.0_wp)/(en*(en+1.0_wp)) * (a1p*tau(n,:) + b1p*pi(n,:))

    end do

    S1a = abs(S1)**2.0_wp
    S1pa = abs(S1p)**2.0_wp

    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start
    
    write(*,*) "writing... output.dat"
    open(10,file="output.dat",status="unknown")
    write(10,*)
    write(10,*) "ang (deg)"," S1", " S1(P=Pc)"
    do i=1,nang
      write(10,*) theta(i),S1a(i),S1pa(i)
    enddo
    ! do i = 1,layers
    !     write(10,*) M(i)%re
    ! end do
    write(10,*)
    close(10)

    stop

end program main_coated