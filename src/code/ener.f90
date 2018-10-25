!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 

module ener
    use sysdef
    use config

    integer*8, dimension(:), allocatable :: neb
    real*8, dimension(:), allocatable :: uu, dui

    contains

subroutine ener_ini
    implicit none

    allocate(uu(n), dui(24), neb(24))
    uu = 0.0 ; dui = 0.0

    return
end subroutine ener_ini


subroutine ener_tot
    implicit none
    integer*8 :: i, kk, ixx, iyy, izz, k, m, mo, ixo, iyo

    do i = 1, n
        do kk = 1, 6
            ixx = modulo(xi(i) + nbx(kk) - 1, nx) + 1
            iyy = modulo(yi(i) + nby(kk) - 1, ny) + 1
            izz = modulo(zi(i) + nbz(kk) - 1, nz) + 1
            !print *, ixx, iyy, izz
            k = ss(ixx,iyy,izz)
            uu(i) = uu(i) + ee(ti(i), ti(k))
        end do

        ! next nearest neighbor interactions
        do kk = 1, 6 
            ixx = modulo(xi(i) + mbx(kk) - 1, nx) + 1
            iyy = modulo(yi(i) + mby(kk) - 1, ny) + 1
            k = ss(ixx,iyy,zi(i))
            uu(i) = uu(i) + es(ti(i), ti(k))
        end do
    end do

    return
end subroutine ener_tot


subroutine ener_dif(i, j)
    implicit none
    integer*8 :: i, j, kk, ixx, iyy, izz, k, m, ixo, iyo

    dui = 0.0

    ! position of old i with type ti(j)
    ! NN
    do kk = 1, 6
        ixx = modulo(tix + nbx(kk) - 1, nx) + 1
        iyy = modulo(tiy + nby(kk) - 1, ny) + 1
        izz = modulo(tiz + nbz(kk) - 1, nz) + 1
        k = ss(ixx,iyy,izz)
        neb(kk) = k
        if (k == j) then
            dui(kk) = 0.0
        else
            dui(kk) = ee(ti(j), ti(k)) - ee(ti(i), ti(k))
        end if
    end do

    ! NNN
    do m = 1, 6 
        ixx = modulo(tix + mbx(m) - 1, nx) + 1
        iyy = modulo(tiy + mby(m) - 1, ny) + 1
        k = ss(ixx,iyy,tiz)
        neb(12+m) = k
        dui(12+m) = (es(ti(j), ti(k)) - es(ti(i), ti(k)))
    end do

    ! position of old j with type ti(i)
    do kk = 1, 6
        ixx = modulo(tjx + nbx(kk) - 1, nx) + 1
        iyy = modulo(tjy + nby(kk) - 1, ny) + 1
        izz = modulo(tjz + nbz(kk) - 1, nz) + 1
        k = ss(ixx,iyy,izz)
        neb(6+kk) = k
        if (k == i) then
            dui(6+kk) = 0.0
        else
            dui(6+kk) = ee(ti(i), ti(k)) - ee(ti(j), ti(k))
        end if
    end do

    do m = 1, 6 
        ixx = modulo(tjx + mbx(m) - 1, nx) + 1
        iyy = modulo(tjy + mby(m) - 1, ny) + 1
        k = ss(ixx,iyy,tjz)
        neb(18+m) = k
        dui(18+m) = (es(ti(i), ti(k)) - es(ti(j), ti(k)))
    end do

    return
end subroutine ener_dif

end module ener

