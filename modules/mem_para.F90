Module mem_para

#ifdef OLAM_MPI
  use mpi_f08, only: MPI_Comm
#endif

  implicit none (external, type)

  integer :: mgroupsize
  integer :: myrank

  integer :: nbytes_int
  integer :: nbytes_real
  integer :: nbytes_real8

  ! Define a separate parallel communicator for OLAM, in case we are running
  ! on a subset of the total available processes for a coupled run.

#ifdef OLAM_MPI
  Type(MPI_Comm) :: MPI_COMM_OLAM
#else
  integer        :: MPI_COMM_OLAM = 0
#endif

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
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Comm_size, MPI_Comm_rank,   &
                     MPI_Init, MPI_Init_thread, MPI_THREAD_MULTIPLE, &
                     MPI_Pack_size, MPI_INTEGER, MPI_REAL, MPI_REAL8
  import,      only: MPI_COMM_OLAM
#endif

  import, only: mgroupsize, myrank, nbytes_int, nbytes_real, nbytes_real8

  implicit none (external, type)

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
     call MPI_Init()
  endif

  ! Default MPI communicator set to all processes for now. This may change for
  ! simulations coupled to other models that use some of the available processes.
  MPI_COMM_OLAM = MPI_COMM_WORLD

  call MPI_Comm_size(MPI_COMM_OLAM,mgroupsize,ierr)
  call MPI_Comm_rank(MPI_COMM_OLAM,myrank,ierr)

  call MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_OLAM, nbytes_int  , ierr)
  call MPI_Pack_size(1, MPI_REAL   , MPI_COMM_OLAM, nbytes_real , ierr)
  call MPI_Pack_size(1, MPI_REAL8  , MPI_COMM_OLAM, nbytes_real8, ierr)

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
  use mpi_f08, only: MPI_Finalize
#endif

  import, none

  implicit none (external, type)

  ! Now terminate MPI

#ifdef OLAM_MPI
  call MPI_Finalize()
#endif

end subroutine olam_mpi_finalize

!===============================================================================

subroutine olam_stop(message)

#ifdef OLAM_MPI
  use mpi_f08, only: MPI_Abort
  import,      only: MPI_COMM_OLAM, myrank
#else
  import,      none
#endif

  implicit none (external, type)

  character(*), intent(in) :: message

#ifdef OLAM_MPI
  write(*,'(A,I0,A)') "Node ", myrank, ":"
  write(*,'(A)') "STOP "//message
  call MPI_Abort(MPI_COMM_OLAM, 1)
  stop
#else
  write(*,*) "STOPPING: "//message
  stop
#endif

end subroutine olam_stop

!===============================================================================

subroutine olam_mpi_barrier()

#ifdef OLAM_MPI
  use mpi_f08, only: MPI_Barrier
  import,      only: MPI_COMM_OLAM
#else
  import,     none
#endif

  implicit none (external, type)

#ifdef OLAM_MPI
  call MPI_Barrier(MPI_COMM_OLAM)
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
  import,         only: mva_primary, iva_globe_primary, iva_local_primary, &
                        mwa_primary, iwa_globe_primary, iwa_local_primary, &
                        mma_primary, ima_globe_primary, ima_local_primary, &
                        mwsfc_primary, iwsfc_globe_primary, iwsfc_local_primary, &
                        mvsfc_primary, ivsfc_globe_primary, ivsfc_local_primary, &
                        mmsfc_primary, imsfc_globe_primary, imsfc_local_primary, &
                        mland_primary, iland_globe_primary, iland_local_primary, &
                        mlake_primary, ilake_globe_primary, ilake_local_primary, &
                        msea_primary, isea_globe_primary, isea_local_primary, &
                        mwnud_primary, iwnud_globe_primary, iwnud_local_primary, &
                        myrank

  implicit none (external, type)

  integer :: i, ia, ib, ic, id, istart, iland, ilake, isea

  if (iparallel /= 1) then

     mwa_primary = mwa
     allocate(iwa_globe_primary(mwa_primary))
     allocate(iwa_local_primary(mwa_primary))
     do i = 1, mwa
        iwa_globe_primary(i) = i
        iwa_local_primary(i) = i
     enddo

     mva_primary = mva
     allocate(iva_globe_primary(mva_primary))
     allocate(iva_local_primary(mva_primary))
     do i = 1, mva
        iva_globe_primary(i) = i
        iva_local_primary(i) = i
     enddo

     mma_primary = mma
     allocate(ima_globe_primary(mma_primary))
     allocate(ima_local_primary(mma_primary))
     do i = 1, mma
        ima_globe_primary(i) = i
        ima_local_primary(i) = i
     enddo

     mwsfc_primary = mwsfc
     allocate(iwsfc_globe_primary(mwsfc_primary))
     allocate(iwsfc_local_primary(mwsfc_primary))
     do i = 1, mwsfc
        iwsfc_globe_primary(i) = i
        iwsfc_local_primary(i) = i
     enddo

     mvsfc_primary = mvsfc
     allocate(ivsfc_globe_primary(mvsfc_primary))
     allocate(ivsfc_local_primary(mvsfc_primary))
     do i = 1, mvsfc
        ivsfc_globe_primary(i) = i
        ivsfc_local_primary(i) = i
     enddo

     mmsfc_primary = mmsfc
     allocate(imsfc_globe_primary(mmsfc_primary))
     allocate(imsfc_local_primary(mmsfc_primary))
     do i = 1, mmsfc
        imsfc_globe_primary(i) = i
        imsfc_local_primary(i) = i
     enddo

     mland_primary = mland
     allocate(iland_globe_primary(mland_primary))
     allocate(iland_local_primary(mland_primary))
     do i = 1, mland
        iland_globe_primary(i) = i
        iland_local_primary(i) = i
     enddo

     msea_primary = msea
     allocate(isea_globe_primary(msea_primary))
     allocate(isea_local_primary(msea_primary))
     do i = 1, msea
        isea_globe_primary(i) = i
        isea_local_primary(i) = i
     enddo

     mlake_primary = mlake
     allocate(ilake_globe_primary(mlake_primary))
     allocate(ilake_local_primary(mlake_primary))
     do i = 1, mlake
        ilake_globe_primary(i) = i
        ilake_local_primary(i) = i
     enddo

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        mwnud_primary = mwnud
        allocate(iwnud_globe_primary(mwnud_primary))
        allocate(iwnud_local_primary(mwnud_primary))
        do i = 1, mwnud
           iwnud_globe_primary(i) = i
           iwnud_local_primary(i) = i
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
