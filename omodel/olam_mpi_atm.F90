Module olam_mpi_atm

  integer, parameter :: itagw    = 1
  integer, parameter :: itagv    = 2
  integer, parameter :: itagm    = 3

  integer :: nsends_v, nrecvs_v
  integer :: nsends_w, nrecvs_w
  integer :: nsends_m, nrecvs_m

  Type nodebuffs
     character, allocatable :: buff(:)
     integer,   allocatable :: ipts(:)
     integer                :: nbytes = 0
     integer                :: jend = 0
     integer                :: iremote = -1
  End Type nodebuffs

  type(nodebuffs), allocatable :: send_w(:), recv_w(:)
  type(nodebuffs), allocatable :: send_v(:), recv_v(:)
  type(nodebuffs), allocatable :: send_m(:), recv_m(:)

  integer              :: icurr_v = 1
  integer              :: inext_v = 2
  integer, allocatable :: ireqr_v(:,:)
  integer, allocatable :: ireqs_v(:,:)

  integer              :: icurr_m = 1
  integer              :: inext_m = 2
  integer, allocatable :: ireqr_m(:,:)
  integer, allocatable :: ireqs_m(:,:)

  integer              :: icurr_w = 1
  integer              :: inext_w = 2
  integer, allocatable :: ireqr_w(:,:)
  integer, allocatable :: ireqs_w(:,:)

Contains

!===============================================================================

subroutine olam_mpi_atm_start()

#ifdef OLAM_MPI
  use var_tables, only: nvar_par
  use mem_grid,   only: mza
  use mem_para,   only: nbytes_int, nbytes_real, nbytes_real8
  use mpi
#endif

  implicit none

#ifdef OLAM_MPI
  integer :: ierr, jrecv, jsend, nv
  integer :: nbytes_per_iw
  integer :: nbytes_per_iv
  integer :: nbytes_per_im

  allocate( ireqs_w(nsends_w,2) ) ; ireqs_w = MPI_REQUEST_NULL
  allocate( ireqr_w(nrecvs_w,2) ) ; ireqr_w = MPI_REQUEST_NULL

  allocate( ireqs_v(nsends_v,2) ) ; ireqs_v = MPI_REQUEST_NULL
  allocate( ireqr_v(nrecvs_v,2) ) ; ireqr_v = MPI_REQUEST_NULL

  allocate( ireqs_m(nsends_m,2) ) ; ireqs_m = MPI_REQUEST_NULL
  allocate( ireqr_m(nrecvs_m,2) ) ; ireqr_m = MPI_REQUEST_NULL

  ! Compute W communication buffer sizes

  nv = max(nvar_par, 20)

  nbytes_per_iw =  2 * mza * nbytes_real8 &
                + nv * mza * nbytes_real

  do jrecv = 1, nrecvs_w
     recv_w(jrecv)%nbytes = recv_w(jrecv)%jend * nbytes_per_iw
     allocate( recv_w(jrecv)%buff( recv_w(jrecv)%nbytes ) )
  enddo

  do jsend = 1, nsends_w
     send_w(jsend)%nbytes = send_w(jsend)%jend * nbytes_per_iw
     allocate( send_w(jsend)%buff( send_w(jsend)%nbytes ) )
  enddo

  ! Compute V communication buffer sizes

  nbytes_per_iv = 3       * nbytes_int &
                + 4 * mza * nbytes_real

  do jrecv = 1, nrecvs_v
     recv_v(jrecv)%nbytes  = recv_v(jrecv)%jend * nbytes_per_iv
     allocate( recv_v(jrecv)%buff( recv_v(jrecv)%nbytes ) )
  enddo

  do jsend = 1, nsends_v
     send_v(jsend)%nbytes  = send_v(jsend)%jend * nbytes_per_iv
     allocate( send_v(jsend)%buff( send_v(jsend)%nbytes ) )
  enddo

  ! Compute M communication buffer sizes

  nbytes_per_im = 2 * mza * nbytes_real

  do jrecv = 1, nrecvs_m
     recv_m(jrecv)%nbytes  = recv_m(jrecv)%jend * nbytes_per_im
     allocate( recv_m(jrecv)%buff( recv_m(jrecv)%nbytes ) )
  enddo

  do jsend = 1, nsends_m
     send_m(jsend)%nbytes  = send_m(jsend)%jend * nbytes_per_im
     allocate( send_m(jsend)%buff( send_m(jsend)%nbytes ) )
  enddo

  ! Pre-post MPI non-blocking receives

  do jrecv = 1, nrecvs_v

     call MPI_Irecv(recv_v(jrecv)%buff, recv_v(jrecv)%nbytes, MPI_PACKED, &
                    recv_v(jrecv)%iremote, itagv, MPI_COMM_WORLD,         &
                    ireqr_v(jrecv,inext_v), ierr)
  enddo

  do jrecv = 1, nrecvs_m

     call MPI_Irecv(recv_m(jrecv)%buff, recv_m(jrecv)%nbytes, MPI_PACKED, &
                    recv_m(jrecv)%iremote, itagm, MPI_COMM_WORLD,         &
                    ireqr_m(jrecv,inext_m), ierr)
  enddo

  do jrecv = 1, nrecvs_w

     call MPI_Irecv(recv_w(jrecv)%buff, recv_w(jrecv)%nbytes, MPI_PACKED, &
                    recv_w(jrecv)%iremote, itagw, MPI_COMM_WORLD,         &
                    ireqr_w(jrecv,inext_w), ierr)
  enddo

#endif

end subroutine olam_mpi_atm_start

!===============================================================================

subroutine mpi_send_v(rvara1, rvara2, rvara3, rvara4, &
                      i1dvara1, i1dvara2, i1dvara3)

! Subroutine to perform a parallel MPI send of a "V group" of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_grid, only: mza, mva

  implicit none

  real, optional, intent(in) :: rvara1(mza,mva)
  real, optional, intent(in) :: rvara2(mza,mva)
  real, optional, intent(in) :: rvara3(mza,mva)
  real, optional, intent(in) :: rvara4(mza,mva)

  integer, optional, intent(in) :: i1dvara1(mva)
  integer, optional, intent(in) :: i1dvara2(mva)
  integer, optional, intent(in) :: i1dvara3(mva)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jtmp, jsend
  integer :: j
  integer :: iv

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_v

     ! Make sure previous sends are finished
     call MPI_Waitany(nsends_v, ireqs_v(:,icurr_v), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,ierr,j,iv) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1, send_v(jsend)%jend
        iv = send_v(jsend)%ipts(j)

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Pack(rvara3(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Pack(rvara4(1,iv),mza,MPI_REAL, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Pack(i1dvara1(iv),1,MPI_INTEGER, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Pack(i1dvara2(iv),1,MPI_INTEGER, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Pack(i1dvara3(iv),1,MPI_INTEGER, &
                send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_v(jsend)%buff, ipos, MPI_PACKED, &
                    send_v(jsend)%iremote, itagv, MPI_COMM_WORLD, &
                    ireqs_v(jsend,inext_v), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_v = mod(icurr_v,2) + 1
  inext_v = mod(inext_v,2) + 1

#endif

end subroutine mpi_send_v

!=============================================================================

subroutine mpi_send_m(rvara1, rvara2, i1dvara1)

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real,    optional, contiguous, intent(in) :: rvara1(:,:)
  real,    optional, contiguous, intent(in) :: rvara2(:,:)
  integer, optional, contiguous, intent(in) :: i1dvara1(:)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jtmp, jsend
  integer :: j, im

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_m

     ! Make sure previous send is finished
     call MPI_Waitany(nsends_m, ireqs_m(:,icurr_m), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,ierr,j,im) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1, send_m(jsend)%jend
        im = send_m(jsend)%ipts(j)

        if (present(rvara1)) then
           call MPI_Pack(rvara1(:,im),size(rvara1,1),MPI_REAL, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(:,im),size(rvara2,1),MPI_REAL, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Pack(i1dvara1(im),1,MPI_INTEGER, &
           send_m(jsend)%buff,send_m(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_m(jsend)%buff, ipos, MPI_PACKED, &
                    send_m(jsend)%iremote, itagm, MPI_COMM_WORLD, &
                    ireqs_m(jsend,inext_m), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_m = mod(icurr_m,2) + 1
  inext_m = mod(inext_m,2) + 1

#endif

end subroutine mpi_send_m

!=============================================================================

subroutine mpi_send_w(scalars,                                         &
                      rvara1,   rvara2,   rvara3,   rvara4,   rvara5,  &
                      rvara6,   rvara7,   rvara8,   rvara9,   rvara10, &
                      rvara11,  rvara12,  rvara13,  rvara14,  rvara15, &
                      rvara16,  rvara17,  rvara18,  rvara19,  rvara20, &
                      r1dvara1, r1dvara2, r1dvara3, r1dvara4, r1dvara5,&
                      dvara1,   dvara2,   i1dvara1, i1dvara2, i1dvara3,&
                      svara1,   svara2,   svara3,   svara4,   svara5,  &
                      swvar1,   swvar2,   swvar3,   swvar4,   swvar5   )

! Subroutine to perform a parallel MPI send of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use var_tables, only: nvar_par, vtab_r, nptonv
  use mem_grid,   only: mza, mwa, lsw, nsw_max
  use consts_coms,only: r8

  implicit none

  character(1), optional, intent(in) :: scalars

  real(r8), optional, intent(in) :: dvara1 (mza,mwa)
  real(r8), optional, intent(in) :: dvara2 (mza,mwa)

  real, optional, contiguous, intent(in) :: svara1(:,:)
  real, optional, contiguous, intent(in) :: svara2(:,:)
  real, optional, contiguous, intent(in) :: svara3(:,:)
  real, optional, contiguous, intent(in) :: svara4(:,:)
  real, optional, contiguous, intent(in) :: svara5(:,:)

  real, optional, intent(in) :: rvara1 (mza,mwa)
  real, optional, intent(in) :: rvara2 (mza,mwa)
  real, optional, intent(in) :: rvara3 (mza,mwa)
  real, optional, intent(in) :: rvara4 (mza,mwa)
  real, optional, intent(in) :: rvara5 (mza,mwa)
  real, optional, intent(in) :: rvara6 (mza,mwa)
  real, optional, intent(in) :: rvara7 (mza,mwa)
  real, optional, intent(in) :: rvara8 (mza,mwa)
  real, optional, intent(in) :: rvara9 (mza,mwa)
  real, optional, intent(in) :: rvara10(mza,mwa)
  real, optional, intent(in) :: rvara11(mza,mwa)
  real, optional, intent(in) :: rvara12(mza,mwa)
  real, optional, intent(in) :: rvara13(mza,mwa)
  real, optional, intent(in) :: rvara14(mza,mwa)
  real, optional, intent(in) :: rvara15(mza,mwa)
  real, optional, intent(in) :: rvara16(mza,mwa)
  real, optional, intent(in) :: rvara17(mza,mwa)
  real, optional, intent(in) :: rvara18(mza,mwa)
  real, optional, intent(in) :: rvara19(mza,mwa)
  real, optional, intent(in) :: rvara20(mza,mwa)

  real, optional, intent(in) :: r1dvara1(mwa)
  real, optional, intent(in) :: r1dvara2(mwa)
  real, optional, intent(in) :: r1dvara3(mwa)
  real, optional, intent(in) :: r1dvara4(mwa)
  real, optional, intent(in) :: r1dvara5(mwa)

  integer, optional, intent(in) :: i1dvara1(mwa)
  integer, optional, intent(in) :: i1dvara2(mwa)
  integer, optional, intent(in) :: i1dvara3(mwa)

  real, optional, intent(in) :: swvar1(nsw_max,mwa)
  real, optional, intent(in) :: swvar2(nsw_max,mwa)
  real, optional, intent(in) :: swvar3(nsw_max,mwa)
  real, optional, intent(in) :: swvar4(nsw_max,mwa)
  real, optional, intent(in) :: swvar5(nsw_max,mwa)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: i, ivar
  integer :: j, iw

  ! Before we send anything, post the receives

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_w

     ! Make sure previous sends are finished
     call MPI_Waitany(nsends_w, ireqs_w(:,icurr_w), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,ierr,j,iw,i,ivar) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1, send_w(jsend)%jend
        iw = send_w(jsend)%ipts(j)

        if (present(dvara1)) then
           call MPI_Pack(dvara1(1,iw),mza,MPI_REAL8, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara2)) then
           call MPI_Pack(dvara2(1,iw),mza,MPI_REAL8, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara1)) then
           call MPI_Pack(svara1(:,iw),size(svara1,1),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara2)) then
           call MPI_Pack(svara2(:,iw),size(svara2,1),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara3)) then
           call MPI_Pack(svara3(:,iw),size(svara3,1),MPI_REAL, &
              send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara4)) then
           call MPI_Pack(svara4(:,iw),size(svara4,1),MPI_REAL, &
              send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara5)) then
           call MPI_Pack(svara5(:,iw),size(svara5,1),MPI_REAL, &
              send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara1)) then
           call MPI_Pack(rvara1(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Pack(rvara2(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Pack(rvara3(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Pack(rvara4(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara5)) then
           call MPI_Pack(rvara5(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara6)) then
           call MPI_Pack(rvara6(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara7)) then
           call MPI_Pack(rvara7(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara8)) then
           call MPI_Pack(rvara8(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara9)) then
           call MPI_Pack(rvara9(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara10)) then
           call MPI_Pack(rvara10(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara11)) then
           call MPI_Pack(rvara11(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara12)) then
           call MPI_Pack(rvara12(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara13)) then
           call MPI_Pack(rvara13(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara14)) then
           call MPI_Pack(rvara14(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara15)) then
           call MPI_Pack(rvara15(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara16)) then
           call MPI_Pack(rvara16(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara17)) then
           call MPI_Pack(rvara17(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara18)) then
           call MPI_Pack(rvara18(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara19)) then
           call MPI_Pack(rvara19(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara20)) then
           call MPI_Pack(rvara20(1,iw),mza,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara1)) then
           call MPI_Pack(r1dvara1(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara2)) then
           call MPI_Pack(r1dvara2(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara3)) then
           call MPI_Pack(r1dvara3(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara4)) then
           call MPI_Pack(r1dvara4(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara5)) then
           call MPI_Pack(r1dvara5(iw),1,MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Pack(i1dvara1(iw),1,MPI_INTEGER, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Pack(i1dvara2(iw),1,MPI_INTEGER, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Pack(i1dvara3(iw),1,MPI_INTEGER, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar1)) then
           call MPI_Pack(swvar1(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar2)) then
           call MPI_Pack(swvar2(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar3)) then
           call MPI_Pack(swvar3(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar4)) then
           call MPI_Pack(swvar4(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar5)) then
           call MPI_Pack(swvar5(1,iw),lsw(iw),MPI_REAL, &
                send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(scalars)) then
           do i = 1, nvar_par
              ivar = nptonv(i)

              call MPI_Pack(vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL, &
                   send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

           enddo
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_w(jsend)%buff, ipos, MPI_PACKED, &
                    send_w(jsend)%iremote, itagw, MPI_COMM_WORLD, &
                    ireqs_w(jsend,inext_w), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_w = mod(icurr_w,2) + 1
  inext_w = mod(inext_w,2) + 1

#endif

end subroutine mpi_send_w

!=============================================================================

subroutine mpi_recv_v(rvara1, rvara2, rvara3, rvara4, &
                      i1dvara1, i1dvara2, i1dvara3)

! Subroutine to perform a parallel MPI receive of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_grid, only: mza, mva

  implicit none

  real, optional, intent(inout) :: rvara1(mza,mva)
  real, optional, intent(inout) :: rvara2(mza,mva)
  real, optional, intent(inout) :: rvara3(mza,mva)
  real, optional, intent(inout) :: rvara4(mza,mva)

  integer, optional, intent(in) :: i1dvara1(mva)
  integer, optional, intent(in) :: i1dvara2(mva)
  integer, optional, intent(in) :: i1dvara3(mva)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: iv

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_v

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_v, ireqr_v(:,icurr_v), jrecv, MPI_STATUS_IGNORE, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task firstprivate(jrecv) private(j,ipos,iv,ierr) default(shared)

     ipos = 0

     do j = 1, recv_v(jrecv)%jend
        iv = recv_v(jrecv)%ipts(j)

        if (present(rvara1)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara1(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara2(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara3(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                rvara4(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                i1dvara1(iv),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                i1dvara2(iv),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos, &
                i1dvara3(iv),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

     enddo

     call MPI_Irecv(recv_v(jrecv)%buff, recv_v(jrecv)%nbytes, MPI_PACKED, &
                    recv_v(jrecv)%iremote, itagv, MPI_COMM_WORLD,         &
                    ireqr_v(jrecv,inext_v), ierr                          )

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_v

!=============================================================================

subroutine mpi_recv_m(rvara1, rvara2, i1dvara1)

! Subroutine to perform a parallel MPI receive of a "M group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real,    optional, contiguous, intent(inout) :: rvara1(:,:)
  real,    optional, contiguous, intent(inout) :: rvara2(:,:)
  integer, optional, contiguous, intent(inout) :: i1dvara1(:)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, jtmp
  integer :: j
  integer :: im

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_m

     !  Now, let's wait on our receives

     call MPI_Waitany(nrecvs_m, ireqr_m(:,icurr_m), jrecv, MPI_STATUS_IGNORE, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     !  We got some stuff.  Now unpack it into appropriate space.

     !$omp task firstprivate(jrecv) private(j,ipos,im,ierr) default(shared)

     ipos = 0

     do j = 1, recv_m(jrecv)%jend
        im = recv_m(jrecv)%ipts(j)

        if (present(rvara1)) then
           call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                rvara1(:,im),size(rvara1,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                rvara2(:,im),size(rvara2,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Unpack(recv_m(jrecv)%buff,recv_m(jrecv)%nbytes,ipos, &
                i1dvara1(im),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

     enddo

     call MPI_Irecv(recv_m(jrecv)%buff, recv_m(jrecv)%nbytes, MPI_PACKED, &
                    recv_m(jrecv)%iremote, itagm, MPI_COMM_WORLD,         &
                    ireqr_m(jrecv,inext_m), ierr                          )

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_m

!=============================================================================

subroutine mpi_recv_w(scalars,                                         &
                      rvara1,   rvara2,   rvara3,   rvara4,   rvara5,  &
                      rvara6,   rvara7,   rvara8,   rvara9,   rvara10, &
                      rvara11,  rvara12,  rvara13,  rvara14,  rvara15, &
                      rvara16,  rvara17,  rvara18,  rvara19,  rvara20, &
                      r1dvara1, r1dvara2, r1dvara3, r1dvara4, r1dvara5,&
                      dvara1,   dvara2,   i1dvara1, i1dvara2, i1dvara3,&
                      svara1,   svara2,   svara3,   svara4,   svara5,  &
                      swvar1,   swvar2,   swvar3,   swvar4,   swvar5   )

! Subroutine to perform a parallel MPI receive of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use var_tables, only: vtab_r, nvar_par, nptonv
  use mem_grid,   only: mza, mwa, lsw, nsw_max
  use consts_coms,only: r8

  implicit none

  character(1), optional, intent(in) :: scalars

  real(r8), optional, intent(inout) :: dvara1 (mza,mwa)
  real(r8), optional, intent(inout) :: dvara2 (mza,mwa)

  real, optional, contiguous, intent(inout) :: svara1(:,:)
  real, optional, contiguous, intent(inout) :: svara2(:,:)
  real, optional, contiguous, intent(inout) :: svara3(:,:)
  real, optional, contiguous, intent(inout) :: svara4(:,:)
  real, optional, contiguous, intent(inout) :: svara5(:,:)

  real, optional, intent(inout) :: rvara1 (mza,mwa)
  real, optional, intent(inout) :: rvara2 (mza,mwa)
  real, optional, intent(inout) :: rvara3 (mza,mwa)
  real, optional, intent(inout) :: rvara4 (mza,mwa)
  real, optional, intent(inout) :: rvara5 (mza,mwa)
  real, optional, intent(inout) :: rvara6 (mza,mwa)
  real, optional, intent(inout) :: rvara7 (mza,mwa)
  real, optional, intent(inout) :: rvara8 (mza,mwa)
  real, optional, intent(inout) :: rvara9 (mza,mwa)
  real, optional, intent(inout) :: rvara10(mza,mwa)
  real, optional, intent(inout) :: rvara11(mza,mwa)
  real, optional, intent(inout) :: rvara12(mza,mwa)
  real, optional, intent(inout) :: rvara13(mza,mwa)
  real, optional, intent(inout) :: rvara14(mza,mwa)
  real, optional, intent(inout) :: rvara15(mza,mwa)
  real, optional, intent(inout) :: rvara16(mza,mwa)
  real, optional, intent(inout) :: rvara17(mza,mwa)
  real, optional, intent(inout) :: rvara18(mza,mwa)
  real, optional, intent(inout) :: rvara19(mza,mwa)
  real, optional, intent(inout) :: rvara20(mza,mwa)

  real, optional, intent(inout) :: r1dvara1(mwa)
  real, optional, intent(inout) :: r1dvara2(mwa)
  real, optional, intent(inout) :: r1dvara3(mwa)
  real, optional, intent(inout) :: r1dvara4(mwa)
  real, optional, intent(inout) :: r1dvara5(mwa)

  integer, optional, intent(inout) :: i1dvara1(mwa)
  integer, optional, intent(inout) :: i1dvara2(mwa)
  integer, optional, intent(inout) :: i1dvara3(mwa)

  real, optional, intent(inout) :: swvar1(nsw_max,mwa)
  real, optional, intent(inout) :: swvar2(nsw_max,mwa)
  real, optional, intent(inout) :: swvar3(nsw_max,mwa)
  real, optional, intent(inout) :: swvar4(nsw_max,mwa)
  real, optional, intent(inout) :: swvar5(nsw_max,mwa)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jrecv, ivar, jtmp
  integer :: i, j
  integer :: iw

  !$omp parallel
  !$omp single
  do jtmp = 1, nrecvs_w

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_w, ireqr_w(:,icurr_w), jrecv, MPI_STATUS_IGNORE, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     ! We got some stuff.  Now unpack it into appropriate space.

     !$omp task firstprivate(jrecv) private(j,ipos,iw,ierr,i,ivar) default(shared)

     ipos = 0

     do j = 1, recv_w(jrecv)%jend
        iw = recv_w(jrecv)%ipts(j)

        if (present(dvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                dvara1(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                dvara2(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara1(:,iw),size(svara1,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara2(:,iw),size(svara2,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara3(:,iw),size(svara3,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara4(:,iw),size(svara4,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(svara5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                svara5(:,iw),size(svara5,1),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara1(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara2(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara3(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara4(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara5(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara6)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara6(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara7)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara7(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara8)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara8(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara9)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara9(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara10)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara10(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara11)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara11(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara12)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara12(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara13)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara13(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara14)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara14(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara15)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara15(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara16)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara16(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara17)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara17(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara18)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara18(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara19)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara19(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(rvara20)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                rvara20(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara1(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara2(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara3(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara4(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(r1dvara5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                r1dvara5(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                i1dvara1(iw),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                i1dvara2(iw),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(i1dvara3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                i1dvara3(iw),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar1)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar1(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar2)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar2(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar3)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar3(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar4)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar4(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(swvar5)) then
           call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                swvar5(1,iw),lsw(iw),MPI_REAL,MPI_COMM_WORLD,ierr)
        endif

        if (present(scalars)) then
           do i = 1, nvar_par
              ivar = nptonv(i)

              call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos, &
                   vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

           enddo
        endif

     enddo

     call MPI_Irecv(recv_w(jrecv)%buff, recv_w(jrecv)%nbytes, MPI_PACKED, &
                    recv_w(jrecv)%iremote, itagw, MPI_COMM_WORLD,         &
                    ireqr_w(jrecv,inext_w), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

#endif

end subroutine mpi_recv_w

!=============================================================================

subroutine olam_mpi_atm_stop()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,  only: olam_mpi_barrier
  use misc_coms, only: iparallel

  implicit none

#ifdef OLAM_MPI

  integer :: ierr, jrecv, ii, jsend
  logical :: flags

  if (iparallel == 1) then

     ! First make sure all processes get here so that all sends are finished

     call olam_mpi_barrier()

     ! Cancel any pending communication

     do ii = 1, 2

        do jrecv = 1, nrecvs_v
           if (ireqr_v(jrecv,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_v(jrecv,ii), ierr)
           endif
        enddo

        do jsend = 1, nsends_v
           if (ireqs_v(jsend,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqs_v(jsend,ii), ierr)
           endif
        enddo

        do jrecv = 1, nrecvs_m
           if (ireqr_m(jrecv,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_m(jrecv,ii), ierr)
           endif
        enddo

        do jsend = 1, nsends_m
           if (ireqs_m(jsend,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqs_m(jsend,ii), ierr)
           endif
        enddo

        do jrecv = 1, nrecvs_w
           if (ireqr_w(jrecv,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqr_w(jrecv,ii), ierr)
           endif
        enddo

        do jsend = 1, nsends_w
           if (ireqs_w(jsend,ii) /= MPI_REQUEST_NULL) then
              call MPI_Cancel(ireqs_w(jsend,ii), ierr)
           endif
        enddo

     enddo

     ! Test that all communication requests have been completed or cancelled

     call olam_mpi_barrier()

     call MPI_Testall(nrecvs_v, ireqr_v(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_v, ireqr_v(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_v, ireqs_v(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_v, ireqs_v(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nrecvs_m, ireqr_m(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_m, ireqr_m(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_m, ireqs_m(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_m, ireqs_m(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nrecvs_w, ireqr_w(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_w, ireqr_w(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_w, ireqs_w(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_w, ireqs_w(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call olam_mpi_barrier()

  endif

  ! Deallocate unused arrays

  if (allocated(send_v)) deallocate(send_v)
  if (allocated(recv_v)) deallocate(recv_v)

  if (allocated(send_m)) deallocate(send_m)
  if (allocated(recv_m)) deallocate(recv_m)

  if (allocated(send_w)) deallocate(send_w)
  if (allocated(recv_w)) deallocate(recv_w)

  if (allocated(ireqs_v)) deallocate(ireqs_v)
  if (allocated(ireqr_v)) deallocate(ireqr_v)

  if (allocated(ireqs_m)) deallocate(ireqs_m)
  if (allocated(ireqr_m)) deallocate(ireqr_m)

  if (allocated(ireqs_w)) deallocate(ireqs_w)
  if (allocated(ireqr_w)) deallocate(ireqr_w)

#endif

end subroutine olam_mpi_atm_stop

!==================================================================

End Module olam_mpi_atm
