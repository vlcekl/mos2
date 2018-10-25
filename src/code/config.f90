!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 

module config
    use rndmod
    use sysdef
    implicit none

    integer*8 :: tix, tiy, tiz, tjx, tjy, tjz
    integer*8 :: n, nx, ny, nz, zs
    integer*8, dimension(:), allocatable :: ti, xi, yi, zi
    integer*8, dimension(:,:,:), allocatable :: ss
    integer*8, dimension(:), allocatable :: nbx, nby, nbz
    integer*8, dimension(:), allocatable :: mbx, mby, mbz, mox, moy
    integer*8, dimension(:,:), allocatable :: ox

    contains

    subroutine config_input(filcfg)
        implicit none
        character*80 :: filcfg
        integer*8 :: i, j, ii, nox
        real*8 :: dummy

        ! read configuration
        open(1, file=filcfg, status='old')
            read(1,*) n
            allocate(ti(0:n), xi(n), yi(n), zi(n))
            ti = 0 ; xi = 0 ; yi = 0 ; zi = 0
 
            read(1,*) nx, ny, nz 
            allocate(ss(nx, ny, nz))
            ss = 0
 
            nt(0) = nx*ny*nz - n
            do i = 1, n
                read(1,*) ti(i), xi(i), yi(i), zi(i)
                ss(xi(i), yi(i), zi(i)) = i
                nt(ti(i)) = nt(ti(i)) + 1
            end do
        close(1, status='keep')

        zs = 1 ! define a surface layer (with special interactions)

        allocate(nbx(6), nby(6), nbz(6))
        nbx = 0 ; nby = 0 ; nbz = 0
        nbx(1) = 2
        nbx(2) = -2
        nbx(3) = 1 
        nby(3) = 1 
        nbx(4) = -1 
        nby(4) = 1
        nbx(5) = 1 
        nby(5) = -1
        nbx(6) = -1 
        nby(6) = -1

        allocate(mbx(6), mby(6), mbz(6))
        mbx = 0 ; mby = 0 ; mbz = 0
        mbx(1) =  3 
        mby(1) =  1
        mbx(2) = -3 
        mby(2) =  1
        mbx(3) =  3
        mby(3) = -1
        mbx(4) = -3
        mby(4) = -1
        mby(5) = 2
        mby(6) = -2

        return
    end subroutine config_input

    subroutine config_output(filcfg)
        implicit none
        character*80 :: filcfg
        integer*8 :: i

        ! write configuration
        open(1, file=filcfg, status='unknown')
            write(1,'(I8)') n
            write(1,'(I3X1I3X1I3)') nx, ny, nz 
            do i = 1, n
                write(1,'(I2X1I3X1I3X1I3)') ti(i), xi(i), yi(i), zi(i)
            end do
        close(1, status='keep')

        return
    end subroutine config_output

end module config
