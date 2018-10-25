!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 

module sysdef
    implicit none
    integer*8 :: ntype
    integer*8, dimension(:), allocatable :: nt, mov
    real*8, dimension(:,:), allocatable :: ee, es, en
    
    contains

    subroutine sysdef_input(filmol)
        implicit none
        character*80 :: filmol, str
        real*8 :: faux
        integer*8 :: i, j, k, npars, iaux

        ! read interaction parameters
        open(1, file=filmol, status='old')
            read(1,*) str, ntype 
            allocate(nt(0:ntype), ee(0:ntype,0:ntype), es(0:ntype,0:ntype), mov(0:ntype), en(0:ntype,0:ntype))
            nt = 0 ; ee = 0.0 ; es = 0.0 ; mov = 0 ; en = 0.0

            read(1,*) str, npars
            do k = 1, npars
                read(1,*) i, j, faux
                ee(i,j) = faux
                ee(j,i) = faux
            end do

            read(1,*) str, npars
            do k = 1, npars
                read(1,*) i, j, faux
                es(i,j) = faux
                es(j,i) = faux
            end do

            read(1,*) str, npars
            do k = 1, npars
                read(1,*) i, iaux ! if movable 1, if not 0
                mov(i) = iaux
            end do

        close(1, status='keep')

        return
    end subroutine sysdef_input

end module sysdef
