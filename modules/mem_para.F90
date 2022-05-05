!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University;
   ! Colorado State University Research Foundation ; ATMET, LLC

   ! This software is free software; you can redistribute it and/or modify it
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version.

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.

   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
   ! (http://www.gnu.org/licenses/gpl.html)
   !----------------------------------------------------------------------------

!===============================================================================

Module mem_para

  implicit none

  integer :: mgroupsize
  integer :: myrank

  integer :: nbytes_int
  integer :: nbytes_real
  integer :: nbytes_real8

  integer                      :: mva_primary
  integer, target, allocatable :: iva_globe_primary(:)
  integer, target, allocatable :: iva_local_primary(:)

  integer                      :: mwa_primary
  integer, target, allocatable :: iwa_globe_primary(:)
  integer, target, allocatable :: iwa_local_primary(:)

  integer                      :: mma_primary
  integer, target, allocatable :: ima_globe_primary(:)
  integer, target, allocatable :: ima_local_primary(:)

  integer                      :: mwsfc_primary
  integer, target, allocatable :: iwsfc_globe_primary(:)
  integer, target, allocatable :: iwsfc_local_primary(:)

  integer                      :: mvsfc_primary
  integer, target, allocatable :: ivsfc_globe_primary(:)
  integer, target, allocatable :: ivsfc_local_primary(:)

  integer                      :: mmsfc_primary
  integer, target, allocatable :: imsfc_globe_primary(:)
  integer, target, allocatable :: imsfc_local_primary(:)

  integer                      :: mland_primary
  integer, target, allocatable :: iland_globe_primary(:)
  integer, target, allocatable :: iland_local_primary(:)

  integer                      :: mlake_primary
  integer, target, allocatable :: ilake_globe_primary(:)
  integer, target, allocatable :: ilake_local_primary(:)

  integer                      :: msea_primary
  integer, target, allocatable :: isea_globe_primary(:)
  integer, target, allocatable :: isea_local_primary(:)

  integer                      :: mwnud_primary
  integer, target, allocatable :: iwnud_globe_primary(:)
  integer, target, allocatable :: iwnud_local_primary(:)

Contains

!===============================================================================

subroutine olam_mpi_init()

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

#ifdef OLAM_MPI

  integer :: ierr, iomp
  !$ integer :: iprov

  ! Initialize MPI and determine process groupsize and myrank

  iomp = 0

  !$ iomp = 1
  !$ call MPI_Init_thread(MPI_THREAD_MULTIPLE, iprov, ierr)
  !$ if (iprov < MPI_THREAD_MULTIPLE) then
  !$    write(*,*) "MPI_THREAD_MULTIPLE is not provided by MPI"
  !$    call olam_stop("Stopping model run")
  !$ endif

  if (iomp == 0) then
     call MPI_Init(ierr)
  endif

  call MPI_Comm_size(MPI_COMM_WORLD,mgroupsize,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

  call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, nbytes_int  , ierr)
  call MPI_Pack_size(1, MPI_REAL   , MPI_COMM_WORLD, nbytes_real , ierr)
  call MPI_Pack_size(1, MPI_REAL8  , MPI_COMM_WORLD, nbytes_real8, ierr)

#else

  mgroupsize = 1
  myrank     = 0

  nbytes_int   = 4
  nbytes_real  = 4
  nbytes_real8 = 8

#endif

end subroutine olam_mpi_init

!===============================================================================

subroutine olam_mpi_finalize()

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  ! Now terminate MPI

#ifdef OLAM_MPI
  integer :: ierr
  call MPI_Finalize(ierr)
#endif

end subroutine olam_mpi_finalize

!===============================================================================

subroutine olam_stop(message)

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

#ifdef OLAM_MPI
  integer :: ierr
#endif

  character(*), intent(in) :: message

#ifdef OLAM_MPI
  write(*,'(A,I0,A)') "Node ", myrank, ":"
  write(*,'(A)') "STOP "//message
  call mpi_abort(MPI_COMM_WORLD,1,ierr)
  stop
#else
  write(*,*) "STOPPING: "//message
  stop
#endif

end subroutine olam_stop

!===============================================================================

subroutine olam_mpi_barrier()

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

#ifdef OLAM_MPI
  integer :: ierr
  call mpi_barrier(MPI_COMM_WORLD, ierr)
#endif

end subroutine olam_mpi_barrier

!===============================================================================

subroutine compute_pario_points()

  use mem_grid,   only: mma, mva, mwa
  use mem_ijtabs, only: itab_v, itab_w, itab_m, itabg_m
  use mem_sfcg,   only: mmsfc, mvsfc, mwsfc, sfcg, itab_msfc, itabg_msfc, &
                        itab_vsfc, itabg_vsfc, itab_wsfc, itabg_wsfc
  use mem_land,   only: itab_land, mland, omland
  use mem_lake,   only: itab_lake, mlake, omlake
  use mem_sea,    only: itab_sea, msea, omsea
  use mem_nudge,  only: mwnud, itab_wnud, nudflag, nudnxp
  use misc_coms,  only: iparallel, mdomain

  implicit none

  integer :: i, ia, ib, ic, id, istart, iland, ilake, isea

  if (iparallel /= 1) then

     mwa_primary = mwa
     allocate(iwa_globe_primary(mwa_primary))
     allocate(iwa_local_primary(mwa_primary))
     do i = 1, mwa
        iwa_globe_primary = i
        iwa_local_primary = i
     enddo

     mva_primary = mva
     allocate(iva_globe_primary(mva_primary))
     allocate(iva_local_primary(mva_primary))
     do i = 1, mva
        iva_globe_primary = i
        iva_local_primary = i
     enddo

     mma_primary = mma
     allocate(ima_globe_primary(mma_primary))
     allocate(ima_local_primary(mma_primary))
     do i = 1, mma
        ima_globe_primary = i
        ima_local_primary = i
     enddo

     mwsfc_primary = mwsfc
     allocate(iwsfc_globe_primary(mwsfc_primary))
     allocate(iwsfc_local_primary(mwsfc_primary))
     do i = 1, mwsfc
        iwsfc_globe_primary = i
        iwsfc_local_primary = i
     enddo

     mvsfc_primary = mvsfc
     allocate(ivsfc_globe_primary(mvsfc_primary))
     allocate(ivsfc_local_primary(mvsfc_primary))
     do i = 1, mvsfc
        ivsfc_globe_primary = i
        ivsfc_local_primary = i
     enddo

     mmsfc_primary = mmsfc
     allocate(imsfc_globe_primary(mmsfc_primary))
     allocate(imsfc_local_primary(mmsfc_primary))
     do i = 1, mmsfc
        imsfc_globe_primary = i
        imsfc_local_primary = i
     enddo

     mland_primary = mland
     allocate(iland_globe_primary(mland_primary))
     allocate(iland_local_primary(mland_primary))
     do i = 1, mland
        iland_globe_primary = i
        iland_local_primary = i
     enddo

     msea_primary = msea
     allocate(isea_globe_primary(msea_primary))
     allocate(isea_local_primary(msea_primary))
     do i = 1, msea
        isea_globe_primary = i
        isea_local_primary = i
     enddo

     mlake_primary = mlake
     allocate(ilake_globe_primary(mlake_primary))
     allocate(ilake_local_primary(mlake_primary))
     do i = 1, mlake
        ilake_globe_primary = i
        ilake_local_primary = i
     enddo

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        mwnud_primary = mwnud
        allocate(iwnud_globe_primary(mwnud_primary))
        allocate(iwnud_local_primary(mwnud_primary))
        do i = 1, mwnud
           iwnud_globe_primary = i
           iwnud_local_primary = i
        enddo
     endif

  else  ! iparallel

     mwa_primary = 0
     do i = 2, mwa
        if (itab_w(i)%irank == myrank) mwa_primary = mwa_primary + 1
     enddo

     mva_primary = 0
     do i = 2, mva
        if (itab_v(i)%irank == myrank) mva_primary = mva_primary + 1
     enddo

     mma_primary = 0
     do i = 2, mma
        if (itab_m(i)%irank == myrank) mma_primary = mma_primary + 1
     enddo

     mmsfc_primary = 0
     do i = 2, mmsfc
        if (itab_msfc(i)%irank == myrank) mmsfc_primary = mmsfc_primary + 1
     enddo

     mvsfc_primary = 0
     do i = 2, mvsfc
        if (itab_vsfc(i)%irank == myrank) mvsfc_primary = mvsfc_primary + 1
     enddo

     mwsfc_primary = 0
     mland_primary = 0
     mlake_primary = 0
     msea_primary  = 0
     do i = 2, mwsfc
        if (itab_wsfc(i)%irank == myrank) then
           mwsfc_primary = mwsfc_primary + 1

           if (sfcg%leaf_class(i) >= 2) then
              mland_primary = mland_primary + 1
           elseif (sfcg%leaf_class(i) == 1) then
              mlake_primary = mlake_primary + 1
           else
              msea_primary = msea_primary + 1
           endif
        endif
     enddo

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        mwnud_primary = 0
        do i = 2, mwnud
           if (itab_wnud(i)%irank == myrank) mwnud_primary = mwnud_primary + 1
        enddo
     endif

     ! additional space for dummy 1st point included with rank 0 only

     if (myrank == 0) then
        mwa_primary = mwa_primary + 1
        mma_primary = mma_primary + 1
        mva_primary = mva_primary + 1
        mmsfc_primary = mmsfc_primary + 1
        mvsfc_primary = mvsfc_primary + 1
        mwsfc_primary = mwsfc_primary + 1
        mland_primary = mland_primary + 1
        mlake_primary = mlake_primary + 1
        msea_primary = msea_primary + 1
        if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
           mwnud_primary = mwnud_primary + 1
        endif
     endif

     ! allocate space for primary global indices

     allocate(iwa_globe_primary(mwa_primary))
     allocate(iwa_local_primary(mwa_primary))

     allocate(iva_globe_primary(mva_primary))
     allocate(iva_local_primary(mva_primary))

     allocate(ima_globe_primary(mma_primary))
     allocate(ima_local_primary(mma_primary))

     allocate(imsfc_globe_primary(mmsfc_primary))
     allocate(imsfc_local_primary(mmsfc_primary))

     allocate(ivsfc_globe_primary(mvsfc_primary))
     allocate(ivsfc_local_primary(mvsfc_primary))

     allocate(iwsfc_globe_primary(mwsfc_primary))
     allocate(iwsfc_local_primary(mwsfc_primary))

     allocate(iland_globe_primary(mland_primary))
     allocate(iland_local_primary(mland_primary))

     allocate(ilake_globe_primary(mlake_primary))
     allocate(ilake_local_primary(mlake_primary))

     allocate(isea_globe_primary(msea_primary))
     allocate(isea_local_primary(msea_primary))

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        allocate(iwnud_globe_primary(mwnud_primary))
        allocate(iwnud_local_primary(mwnud_primary))
     endif

     ! dummy 1st point included with rank 0 only

     if (myrank == 0) then

        iwa_globe_primary(1) = 1
        iwa_local_primary(1) = 1

        iva_globe_primary(1) = 1
        iva_local_primary(1) = 1

        ima_globe_primary(1) = 1
        ima_local_primary(1) = 1

        imsfc_globe_primary(1) = 1
        imsfc_local_primary(1) = 1

        ivsfc_globe_primary(1) = 1
        ivsfc_local_primary(1) = 1

        iwsfc_globe_primary(1) = 1
        iwsfc_local_primary(1) = 1

        iland_globe_primary(1) = 1
        iland_local_primary(1) = 1

        ilake_globe_primary(1) = 1
        ilake_local_primary(1) = 1

        isea_globe_primary(1) = 1
        isea_local_primary(1) = 1

        if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
           iwnud_globe_primary(1) = 1
           iwnud_local_primary(1) = 1
        endif

     endif

     ! set locations of global and local primary points

     if (myrank == 0) then
        istart = 1
     else
        istart = 0
     endif

     ia = istart

     do i = 2, mwa
        if (itab_w(i)%irank == myrank) then
           ia = ia + 1
           iwa_globe_primary(ia) = itab_w(i)%iwglobe
           iwa_local_primary(ia) = i
        endif
     enddo

     ia = istart

     do i = 2, mva
        if (itab_v(i)%irank == myrank) then
           ia = ia + 1
           iva_globe_primary(ia) = itab_v(i)%ivglobe
           iva_local_primary(ia) = i
        endif
     enddo

     ia = istart

     do i = 2, mma
        if (itabg_m( itab_m(i)%imglobe )%irank == myrank) then
           ia = ia + 1
           ima_globe_primary(ia) = itab_m(i)%imglobe
           ima_local_primary(ia) = i
        endif
     enddo

     ia = istart

     do i = 2, mmsfc
        if (itabg_msfc( itab_msfc(i)%imglobe )%irank == myrank) then
           ia = ia + 1
           imsfc_globe_primary(ia) = itab_msfc(i)%imglobe
           imsfc_local_primary(ia) = i
        endif
     enddo

     ia = istart

     do i = 2, mvsfc
        if (itabg_vsfc( itab_vsfc(i)%ivglobe )%irank == myrank) then
           ia = ia + 1
           ivsfc_globe_primary(ia) = itab_vsfc(i)%ivglobe
           ivsfc_local_primary(ia) = i
        endif
     enddo

     ia = istart
     ib = istart
     ic = istart
     id = istart

     do i = 2, mwsfc
        if (itabg_wsfc( itab_wsfc(i)%iwglobe )%irank == myrank) then
           ia = ia + 1
           iwsfc_globe_primary(ia) = itab_wsfc(i)%iwglobe
           iwsfc_local_primary(ia) = i

           if (sfcg%leaf_class(i) >= 2) then
              iland = i - omland
              ib = ib + 1
              iland_globe_primary(ib) = itab_land(iland)%iwglobe
              iland_local_primary(ib) = iland
           elseif (sfcg%leaf_class(i) == 1) then
              ilake = i - omlake
              ic = ic + 1
              ilake_globe_primary(ic) = itab_lake(ilake)%iwglobe
              ilake_local_primary(ic) = ilake
           else
              isea = i - omsea
              id = id + 1
              isea_globe_primary(id) = itab_sea(isea)%iwglobe
              isea_local_primary(id) = isea
           endif
        endif
     enddo

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        ia = istart

        do i = 2, mwnud
           if (itab_wnud(i)%irank == myrank) then
              ia = ia + 1
              iwnud_globe_primary(ia) = itab_wnud(i)%iwnudglobe
              iwnud_local_primary(ia) = i
           endif
        enddo
     endif
  endif

end subroutine compute_pario_points

!===============================================================================

End Module mem_para
