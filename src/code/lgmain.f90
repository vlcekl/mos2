!
!  File name: lg.f90
!  Date:      2011/05/15 23:49
!  Author:    Lukas Vlcek
! 


program lg
    use rndmod
    use sysdef
    use config
    use ener
    use move
    use measure
    
    implicit none

    character*80 :: str, filmol, filcfg, fnam
    integer*8 :: i, ix, iy, iz, it, itmax, ierr, try, acc
    integer*8 :: ic, idx, idy, idz, j, kk, tt
    real*8 :: dummy, b, Tr, du, cond
    !OpenMP variables
    integer*4 :: tid, nthreads, tblock
!$    integer*4 :: omp_get_num_threads, omp_get_thread_num

    ! initialize random number generator 
    seed = 113314724325.0

    ! read simulation parameters
    read(*,*) str, Tr     ! reduced temperature
    read(*,*) str, itmax  ! reduced temperature
    read(*,*) str, filmol
    read(*,*) str, filcfg
    b = 1.0/Tr

    call sysdef_input(filmol)
    call config_input(filcfg)
    call ener_ini
    call move_ini

    call ener_tot

    call measure_ini
    call measure_do(0_8)

    print *, 0, Tr, 0.5*sum(uu), 0, 0

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ss,b,n,nx,ny,nz,seed)

    tid = 0
    nthreads = 1
!$  tid = omp_get_thread_num()
!$  nthreads = omp_get_num_threads()
    tblock = nz/nthreads
    !print *, 'threads', tid, nthreads


    do it = 1, itmax
        acc = 0
        try = 0
        do ic = 1, n

            !call move_swap

            ! select particle location for move (exchange with neighbor)
            idx = tid/4*(nx/4)
            idy = modulo(tid, 4)*(ny/4)
            idx = 0
            idy = 0
            idz = 0
            do 
                tix = modulo(int(rnd(dummy)*float(nx)) + idx, nx) + 1
                tiy = modulo(int(rnd(dummy)*float(ny)) + idy, ny) + 1 
                tiz = modulo(int(rnd(dummy)*float(nz)) + idz, nz) + 1
                i = ss(tix, tiy, tiz)
                if (ti(i) > 0) exit
            end do

            ! select direction for move
            kk = int(rnd(dummy)*6.0) + 1
            tjx = modulo(tix + nbx(kk) - 1, nx) + 1
            tjy = modulo(tiy + nby(kk) - 1, ny) + 1
            tjz = modulo(tiz + nbz(kk) - 1, nz) + 1
            j = ss(tjx, tjy, tjz)

            !print *, kk, nbx(kk), nby(kk), nbz(kk)
            !print *, 'i', i, tix, tiy, tiz, ti(i), mov(ti(i))
            !print *, 'j', j, tjx, tjy, tjz, ti(j), mov(ti(j))

            cond = rnd(dummy)

            !if (mov(ti(i))*mov(ti(j)) == 0) print *, '*', mov(ti(i)), mov(ti(j)), ti(i), ti(j), i, j, sum(mov(ti(:)))
            if (mov(ti(i))*mov(ti(j)) == 0) cycle
            !if (mov(ti(i))*mov(ti(j)) /= 0) print *, '-', mov(ti(i)), mov(ti(j)), ti(i), ti(j), i, j, sum(mov(ti(:)))

            try = try + 1

            if (ti(i) == ti(j)) then
                dui = 0.0
                neb = 1
                call move_accept(i, j)
                acc = acc + 1
            else
                call ener_dif(i, j)
                du = sum(dui)
                
                !if (mov(ti(i))*mov(ti(j)) == 0) print *, i, j, ti(i), ti(j), mov(ti(i)), mov(ti(j))
                
                if (exp(-b*du) > cond) then
                    call move_accept(i, j)
                    acc = acc + 1
                end if
            end if

!$OMP BARRIER
        end do
!$OMP FLUSH(ss)
!$OMP MASTER
        if (modulo(it,10) == 0) then
            print *, it, Tr, 0.5*sum(uu), try, acc
            call measure_do(it)
            if (modulo(it,1000) == 0) then
                write(fnam, '(i6)') it/1000
                fnam = adjustl(fnam)
                ierr = len_trim(fnam)
                fnam(ierr+1:ierr+4) = '.dmp'
                call config_output(fnam)
            end if
        end if

!$OMP END MASTER
!$OMP BARRIER

    end do

!$OMP END PARALLEL

    call measure_output

    stop
end program lg
