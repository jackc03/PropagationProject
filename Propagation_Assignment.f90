program FireLaser
    implicit none
    double precision :: lambda  = 8.0d-7
    double precision :: pi      = 3.1415926535d0
    double precision :: tau     = 5.0d-15
    double precision :: xstart  = -4.0d-5
    double precision :: xfinal  = 4.0d-5
    double precision :: tstart  = 0d0
    double precision :: tfinal  = 10d-15
    double precision :: vs      = 3.0d8 !scalar velocity
    double precision :: xp      = 0
    double precision :: deltax, deltat, f, omega, k, t1



    call Propagation()

contains
    subroutine Propagation()
        double precision, allocatable :: t(:), x(:), E0(:), E1(:), E2(:), v(:), frontratio(:)
        integer :: Nx, Nt, n, i
        
        deltax       = (lambda) / 40
        deltat       = deltax / (vs)
        f            = vs / lambda
        omega        = 2 * pi * f
        k            = (2 * pi) / lambda
        Nx           = (xfinal - xstart) / deltax
        Nt           = (tfinal - tstart) / deltat

        allocate(t(Nt),x(Nx),E0(Nx),E1(Nx),E2(Nx),v(Nx),frontratio(Nx))
       
        v            = vs
        v(3000:3200) = vs / 2
        t = 0d0
        x = 0d0
        E0 = 0d0
        E1 = 0d0
        E2 = 0d0
        frontratio = 1d0
        v = vs

        x(1) = xstart
        t(1) = tstart

        do i = 2, Nx, 1
            x(i) = x(i-1) + deltax
        end do
        do i = 2, Nt, 1        
            t(i) = t(i-1) + deltat
        end do
        
        xp = x(Nx/4)


        do i = 1, Nx, 1
            E0(i) = EXP( - ( ( (t(1) - x(i) / vs) + (xp / vs ) )**2 ) / ((tau **2))) * cos((k * x(i)) - (omega * t(1)))
            E1(i) = EXP( - ( ( (t(2) - x(i) / vs) + (xp / vs ) )**2 ) / ((tau **2))) * cos((k * x(i)) - (omega * t(2)))
        end do

        do i = 1, Nx, 1
            frontratio(i) = ((v(i)) **2 *deltat **2 / deltax **2)
        end do

        open(file='EiRefractionFortran.dat', unit=101)
        do i = 1, Nx, 1
            write(101,*) x(i), E0(i)
        end do
        close(101)

        open(file='Ei2RefractionFortran.dat', unit=101)
        do i = 1, Nx, 1
            write(101,*) x(i), E1(i)
        end do
        close(101)


        do n = 1, Nt, 1
            i=1
            E2(i) = frontratio(i) * (E1(i+1) - 2 * E1(i) + E1(Nx)) + 2 * E1(i) - E0(i)
            do i = 2, Nx - 1, 1
                E2(i) = frontratio(i) * (E1(i+1) - 2 * E1(i) + E1(i-1)) + 2 * E1(i) - E0(i)
            end do
            i=Nx
            E2(i) = frontratio(i) * (E1(1) - 2 * E1(i) + E1(i-1)) + 2 * E1(i) - E0(i)
            E0 = E1
            E1 = E2
            !do i = 1, Nx-1, 1
            !    E0(i) = E1(i)
            !    E1(i) = E2(i)
            !end do
                       
        end do

        open(file='EfRefractionFortran.dat', unit=102)
        do i = 1,Nx-1, 1
            write(102,*) x(i), E2(i)
        end do
        close(102)


        
        deallocate(t, x, E0, E1, E2, v, frontratio)



    end subroutine
end program
