subroutine global_profiles()

  use mem_basic,       only: rr_v, tair, theta, press, rho
  use mem_grid,        only: mza, zt
  use mem_ijtabs,      only: jtab_w, jtw_prog
  use mem_para,        only: myrank
  use iso_fortran_env, only: r8=>real64

# ifdef OLAM_MPI
  use misc_coms,       only: iparallel
  use mem_para,        only: MPI_COMM_OLAM
  use mpi_f08,         only: MPI_Reduce, MPI_IN_PLACE, MPI_REAL8, MPI_SUM
# endif

  implicit none

  integer, parameter :: nr = 6  ! latitude bands
  integer, parameter :: na = 3  ! 1=all, 2=sea, 3=land
  integer, parameter :: nv = 5  ! # of variables computed

  real, parameter :: glats(nr) = [ 90.1, 90.1, 58.0,  25.0, -25.0, -58.0]
  real, parameter :: glate(nr) = [-90.1, 58.0, 25.0, -25.0, -58.0, -90.1]

  real(r8), allocatable       :: one     (:)
  real(r8), allocatable       :: buffer  (:,:,:,:)
  real(r8), allocatable, save :: arwi_tot(:,:,:)
  logical,               save :: firstime = .true.
  integer,               save :: iout     =  0
  integer                     :: j, iw, n, unit
  character(4)                :: znum
  character(99)               :: fname

  interface dosums
     procedure dosum4
     procedure dosum8
  end interface dosums

  ! First store the sum of the areas to be used to normalize the area sums.
  ! For a parallel run they are only saved on rank 0.

  if (firstime) then

     allocate( arwi_tot(mza,nr,na), source=0._r8 )
     allocate( one     (mza),       source=1._r8 )

     !$omp parallel do private(iw) reduction(+:arwi_tot)
     do j = 1, jtab_w(jtw_prog)%jend ; iw = jtab_w(jtw_prog)%iw(j)

        call dosums( arwi_tot, one, iw )

     enddo
     !$omp end parallel do

#    ifdef OLAM_MPI
     if (iparallel == 1) then
        n = mza*nr*na
        if (myrank == 0) call MPI_Reduce(MPI_IN_PLACE, arwi_tot, n, MPI_REAL8, MPI_SUM, 0, MPI_COMM_OLAM)
        if (myrank /= 0) call MPI_Reduce(arwi_tot,     0._r8,    n, MPI_REAL8, MPI_SUM, 0, MPI_COMM_OLAM)
     endif
#    endif

     if (myrank == 0) arwi_tot(2:,:,:) = 1._r8 / arwi_tot(2:,:,:)
     if (myrank /= 0) deallocate( arwi_tot )

     firstime = .false.

  endif

  ! Now we can compute the area sums

  allocate( buffer(mza,nr,na,nv), source=0._r8 )

  !$omp parallel do private(iw) reduction(+:buffer)
  do j = 1,jtab_w(jtw_prog)%jend ; iw = jtab_w(jtw_prog)%iw(j)

     call dosums( buffer(:,:,:,1), tair (:,iw), iw )
     call dosums( buffer(:,:,:,2), theta(:,iw), iw )
     call dosums( buffer(:,:,:,3), rho  (:,iw), iw )
     call dosums( buffer(:,:,:,4), press(:,iw), iw )
     call dosums( buffer(:,:,:,5), rr_v (:,iw), iw )

  enddo
  !$omp end parallel do

# ifdef OLAM_MPI
  if (iparallel == 1) then
     n = mza*nr*na*nv
     if (myrank == 0) call MPI_Reduce(MPI_IN_PLACE, buffer, n, MPI_REAL8, MPI_SUM, 0, MPI_COMM_OLAM)
     if (myrank /= 0) call MPI_Reduce(buffer,       0._r8,  n, MPI_REAL8, MPI_SUM, 0, MPI_COMM_OLAM)
  endif
# endif

  if (myrank == 0) then

     do n = 1, nv
        buffer(2:,:,:,n) = buffer(2:,:,:,n) * arwi_tot(2:,:,:)
     enddo

     iout = iout + 1
     write(znum,'(I0)') iout

     fname = "profile_" // trim(adjustl(znum)) // ".dat"

     open(newunit=unit, file=fname, form="unformatted", access="stream", status="replace")
     write(unit) mza-1, nr, na, nv
     write(unit) zt(2:)
     write(unit) real( buffer(2:,:,:,:) )
     close(unit)

  endif


contains


  subroutine dosum8( buffer, var8, iw )

    use mem_turb, only: frac_sea
    use mem_grid, only: mza, glatw, lpw, arw0
    import,       only: r8, mza, nr, na, glats, glate

    implicit none (external, type)

    real(r8), intent(inout) :: buffer(mza,nr,na)
    real(r8), intent(in)    :: var8  (mza)
    integer,  intent(in)    :: iw
    integer                 :: n, ka

    ka = lpw(iw)

    do n = 1, nr
       if (glatw(iw) <= glats(n) .and. glatw(iw) > glate(n)) then

          buffer(ka:mza,n,1)    = buffer(ka:mza,n,1) + arw0(iw) * var8(ka:mza)

          if (frac_sea(iw) > .5) then
             buffer(ka:mza,n,2) = buffer(ka:mza,n,2) + arw0(iw) * var8(ka:mza)
          else
             buffer(ka:mza,n,3) = buffer(ka:mza,n,3) + arw0(iw) * var8(ka:mza)
          endif

       endif
    enddo

  end subroutine dosum8


  subroutine dosum4( buffer, var4, iw )

    use mem_turb, only: frac_sea
    use mem_grid, only: mza, glatw, lpw, arw0
    import,       only: r8, mza, nr, na, glats, glate

    implicit none (external, type)

    real(r8), intent(inout) :: buffer(mza,nr,na)
    real,     intent(in)    :: var4  (mza)
    integer,  intent(in)    :: iw
    integer                 :: n, ka

    ka = lpw(iw)

    do n = 1, nr
       if (glatw(iw) <= glats(n) .and. glatw(iw) > glate(n)) then

          buffer(ka:mza,n,1)    = buffer(ka:mza,n,1) + arw0(iw) * var4(ka:mza)

          if (frac_sea(iw) > .5) then
             buffer(ka:mza,n,2) = buffer(ka:mza,n,2) + arw0(iw) * var4(ka:mza)
          else
             buffer(ka:mza,n,3) = buffer(ka:mza,n,3) + arw0(iw) * var4(ka:mza)
          endif

       endif
    enddo

  end subroutine dosum4


end subroutine global_profiles
