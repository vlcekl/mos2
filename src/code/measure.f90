!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 

module measure
    use sysdef
    use config
    use ener

    integer*8 :: icount
    integer*8, dimension(:), allocatable :: hsu
    integer*8, dimension(:,:), allocatable :: hz, hu, hs, hsd, hse, hsv
    real*8, dimension(:), allocatable :: shsu
    real*8, dimension(:,:), allocatable :: shz, shu, shs, shsv

    contains

subroutine measure_ini
    implicit none
    integer*8 :: ix, iy, iz, tt
    integer*8, dimension(:,:), allocatable :: gz

    allocate(hsu(128), hsv(2,6), hz(0:ntype, nz), hu(0:ntype,0:ntype))
    allocate(hs(0:ntype,0:ntype), hsd(0:ntype,0:ntype), hse(0:ntype,0:ntype))
    hsu = 0 ; hz = 0 ; hu = 0 ; hs = 0 ; hsd = 0 ; hse = 0
    allocate(shsu(128), shsv(2, 6), shz(0:ntype, nz), shu(0:ntype,0:ntype), shs(0:ntype,0:ntype))
    shsu = 0 ; shz = 0 ; shu = 0 ; shs = 0 ; shsv = 0.0

    icount = 0

    open(3, file='lg.hst', status='unknown')
    write(3,*) '# histograms'
    close(3, status='keep')

    return
end subroutine measure_ini


subroutine measure_do(it)
    implicit none
    integer*8 :: it, ix, iy, iz, ixx, iyy, izz, tt, tj, kk, i, j, dx, dy, mm, m
    integer*8 :: ixo, iyo
    integer*8, dimension(:,:), allocatable :: gz

    ! z-profile
    hz = 0
    do iz = 1, nz
        do iy = 1, ny
            do ix = 1, nx
                tt = ti(ss(ix,iy,iz))
                hz(tt, iz) = hz(tt, iz) + 1
            end do
        end do
    end do

    ! surface statistics
    hsu = 0
    hsv = 0
    iz = zs ! nz-1
    do ix = 1, nx
        do iy = 1, ny
            tt = ti(ss(ix,iy,iz))
            if (tt == 0) cycle
            dx = tt - 1
            do kk = 1, 6
                ixx = modulo(ix + nbx(kk) - 1, nx) + 1
                iyy = modulo(iy + nby(kk) - 1, ny) + 1
                tj = ti(ss(ixx,iyy,iz))
                dx = dx + (tj - 1)*2**kk
            end do
            dx = dx + 1
            hsu(dx) = hsu(dx) + 1
        end do
    end do

    ! pair interaction statistics
    hu = 0
    hs = 0
    hsd = 0
    hse = 0
    do iz = 1, nz
        do iy = 1, ny
            do ix = 1, nx
                i = ss(ix, iy, iz)
                tt = ti(i)
                if (tt == 0) cycle
                do kk = 1, 6
                    ixx = modulo(ix + nbx(kk) - 1, nx) + 1
                    iyy = modulo(iy + nby(kk) - 1, ny) + 1
                    izz = modulo(iz + nbz(kk) - 1, nz) + 1
                    j = ss(ixx,iyy,izz)
                    if (j >= i) then
                        tj = ti(j)
                        if (tt <= tj) then
                            hu(tt, tj) = hu(tt, tj) + 1
                        else
                            hu(tj, tt) = hu(tj, tt) + 1
                        end if
                    end if
                end do

        ! next nearest neighbor interactions depending on the presence of oxygen
                do m = 1, 6 
                    ixx = modulo(xi(i) + mbx(m) - 1, nx) + 1
                    iyy = modulo(yi(i) + mby(m) - 1, ny) + 1
                    j = ss(ixx,iyy,iz)
                    if (j > i) then
                        tj = ti(j)
                        if (tt <= tj) then
                            hs(tt, tj) = hs(tt, tj) + 1
                        else
                            hs(tj, tt) = hs(tj, tt) + 1
                        end if
                    end if
                end do
            end do
        end do
    end do

    open(3, file='lg.hst', access='append', status='old')

    write(3,*) '#', it, 0.5*sum(uu), sum(ee*hu)+ sum(es*hs)

    write(3,*) '# hz'
    do iz = 1, nz
        write(3,*) (hz(tt, iz), tt=0,ntype)
    end do

    write(3,*) '# hsu'
    do tt = 1, 128
        write(3,*) tt, hsu(tt)
    end do
!    write(3,*) '# hsv'
!    do tt = 1, 128
!        write(3,*) tt, hsv(1,tt), hsv(2,tt)
!    end do

    write(3,*) '# hu'
    do tt = 1, ntype
        !if (mov(tt) /= 1) cycle
        do tj = tt, ntype
            write(3,*) tt, tj, hu(tt,tj), hs(tt,tj)
        end do
    end do

    close(3, status='keep')

    icount = icount + 1
    shz = shz + hz
    shu = shu + hu
    shs = shs + hs
    shsu = shsu + hsu
!    shsv = shsv + hsv

    return
end subroutine measure_do

subroutine measure_output
    implicit none
    integer*8 :: ix, iy, iz, tt, tj

    shz = shz/float(icount)
    shu = shu/float(icount)
    shs = shs/float(icount)
    shsu = shsu/float(icount)
    shsv = shsv/float(icount)

    open(3, file='lg.hst', access='append', status='old')

    write(3,"(A6)") 'ENDHST'

    write(3,*) '# hz'
    do iz = 1, nz
        write(3,*) (shz(tt, iz), tt=1,ntype)
    end do

    write(3,*) '# hsu'
    do tt = 1, 128
        write(3,*) tt, shsu(tt)
    end do

    write(3,*) '# hu'
    do tt = 1, ntype
        !if (mov(tt) /= 1) cycle
        do tj = tt, ntype
            write(3,*) tt, tj, shu(tt,tj), shs(tt,tj)
        end do
    end do

    close(3, status='keep')

    return
end subroutine measure_output

end module measure
