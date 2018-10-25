!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 

module move
    use sysdef
    use config
    use ener
    implicit none
    integer*8 :: itry, jtry

    contains
        
    subroutine move_ini
        implicit none


        return
    end subroutine move_ini

!    subroutine move_swap(i, j)
!        implicit none
!        integer*8 :: idx, idy, idz
!
!        ! select particle location for move (exchange with neighbor)
!        idx = tid/4*(nx/4)
!        idy = modulo(tid, 4)*(ny/4)
!        idz = 0
!        tix = modulo(int(rnd(dummy)*float(nx)) + idx, nx) + 1
!        tiy = modulo(int(rnd(dummy)*float(ny)) + idy, ny) + 1 
!        tiz = modulo(int(rnd(dummy)*float(nz)) + idz, nz) + 1
!        itry = ss(ix, iy, iz)
!
!        ! select direction for move
!        kk = int(rnd(dummy)*6.0) + 1
!        tjx = modulo(ix + nbx(kk), nx) + 1
!        tjy = modulo(iy + nby(kk), ny) + 1
!        tjz = modulo(iz + nbz(kk), nz) + 1
!        jtry = ss(dx, dy, dz)
!
!        return
!    end subroutine move_swap

    subroutine move_accept(i, j)
        implicit none
        integer*8 :: i, j, k, l
        real*8 :: xau, yau, zau

        uu(i) = uu(i) + sum(dui(1:6))
        uu(j) = uu(j) + sum(dui(7:12))

        uu(i) = uu(i) + sum(dui(13:18))
        uu(j) = uu(j) + sum(dui(19:24))
    
        do l = 1, 24
            k = neb(l)
            uu(k) = uu(k) + dui(l)
        end do
    
        xi(i) = tjx
        yi(i) = tjy
        zi(i) = tjz

        xi(j) = tix
        yi(j) = tiy
        zi(j) = tiz

        ss(xi(i), yi(i), zi(i)) = i
        ss(xi(j), yi(j), zi(j)) = j

        return
    end subroutine move_accept

    
end module move
