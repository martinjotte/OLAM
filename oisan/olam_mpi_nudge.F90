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

Module olam_mpi_nud

  integer, parameter :: itagwnud = 21

  integer :: nsends_wnud, nrecvs_wnud

  Type nodebuffs
     character, allocatable :: buff(:)
     integer,   allocatable :: ipts(:)
     integer                :: nbytes = 0
     integer                :: jend = 0
     integer                :: iremote = -1
  End Type nodebuffs

  type(nodebuffs), allocatable :: send_wnud(:), recv_wnud(:)

  integer              :: icurr_wnud = 1
  integer              :: inext_wnud = 2
  integer, allocatable :: ireqr_wnud(:,:)
  integer, allocatable :: ireqs_wnud(:,:)

Contains

!===============================================================================

subroutine olam_mpi_nud_start()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_nudge, only: nudflag, nudnxp
  use misc_coms, only: mdomain

  implicit none

#ifdef OLAM_MPI
  integer :: ierr, jrecv

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then

     allocate( ireqs_wnud(nsends_wnud,2) ) ; ireqs_wnud = MPI_REQUEST_NULL
     allocate( ireqr_wnud(nrecvs_wnud,2) ) ; ireqr_wnud = MPI_REQUEST_NULL

     do jrecv = 1, nrecvs_wnud

        call MPI_Irecv(recv_wnud(jrecv)%buff, recv_wnud(jrecv)%nbytes, MPI_PACKED, &
                       recv_wnud(jrecv)%iremote, itagwnud, MPI_COMM_WORLD,         &
                       ireqr_wnud(jrecv,inext_wnud), ierr)
     enddo
  endif

#endif

end subroutine olam_mpi_nud_start

!=============================================================================

subroutine mpi_send_wnud(dvara1, dvara2, dvara3, dvara4, dvara5, dvara6)

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_grid,   only: mza
  use mem_nudge,  only: mwnud
  use consts_coms,only: r8

  implicit none

  real(r8), optional, intent(in) :: dvara1(mza,mwnud)
  real(r8), optional, intent(in) :: dvara2(mza,mwnud)
  real(r8), optional, intent(in) :: dvara3(mza,mwnud)
  real(r8), optional, intent(in) :: dvara4(mza,mwnud)
  real(r8), optional, intent(in) :: dvara5(mza,mwnud)
  real(r8), optional, intent(in) :: dvara6(mza,mwnud)

#ifdef OLAM_MPI

  integer :: ierr, ipos
  integer :: jsend, jtmp
  integer :: j, iwnud

  !$omp parallel
  !$omp single
  do jtmp = 1, nsends_wnud

     ! Make sure previous sends are finished
     call MPI_Waitany(nsends_wnud, ireqs_wnud(:,icurr_wnud), jsend, MPI_STATUS_IGNORE, ierr)

     if (jsend == MPI_UNDEFINED) jsend = jtmp

     !$omp task private(ipos,ierr,j,iwnud) firstprivate(jsend) default(shared)

     ! Pack the messages into send buffers

     ipos = 0

     do j = 1, send_wnud(jsend)%jend
        iwnud = send_wnud(jsend)%ipts(j)

        if (present(dvara1)) then
           call MPI_Pack(dvara1(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara2)) then
           call MPI_Pack(dvara2(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara3)) then
           call MPI_Pack(dvara3(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara4)) then
           call MPI_Pack(dvara4(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara5)) then
           call MPI_Pack(dvara5(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

        if (present(dvara6)) then
           call MPI_Pack(dvara6(1,iwnud),mza,MPI_REAL8, &
           send_wnud(jsend)%buff,send_wnud(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
        endif

     enddo

     ! Now we can actually go on to sending the stuff

     call MPI_Isend(send_wnud(jsend)%buff, ipos, MPI_PACKED, &
                    send_wnud(jsend)%iremote, itagwnud, MPI_COMM_WORLD, &
                    ireqs_wnud(jsend, inext_wnud), ierr)

     !$omp end task

  enddo
  !$omp end single
  !$omp end parallel

  ! Increment MPI request pointers

  icurr_wnud = mod(icurr_wnud,2) + 1
  inext_wnud = mod(inext_wnud,2) + 1

#endif

end subroutine mpi_send_wnud

!=============================================================================

subroutine mpi_recv_wnud(dvara1, dvara2, dvara3, dvara4, dvara5, dvara6)

! Subroutine to perform a parallel MPI receive of a "WNUD group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_grid,    only: mza
  use mem_nudge,   only: mwnud
  use consts_coms, only: r8

  implicit none

  real(r8), optional, intent(inout) :: dvara1(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara2(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara3(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara4(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara5(mza,mwnud)
  real(r8), optional, intent(inout) :: dvara6(mza,mwnud)

#ifdef OLAM_MPI

  integer  :: ierr, ipos
  integer  :: jrecv, jtmp
  integer  :: j, iwnud
  real(r8) :: vctr1(mza)

  do jtmp = 1, nrecvs_wnud

     ! Now, let's wait on our receives

     call MPI_Waitany(nrecvs_wnud, ireqr_wnud(:,icurr_wnud), jrecv, MPI_STATUS_IGNORE, ierr)

     if (jrecv == MPI_UNDEFINED) jrecv = jtmp

     ! We got some stuff.  Now unpack it into appropriate space.

     ! For now we don't use OpenMP to unpack these buffers because more then one receive
     ! may be added to an individual iwnud point.

     ipos = 0

     do j = 1, recv_wnud(jrecv)%jend
        iwnud = recv_wnud(jrecv)%ipts(j)

        if (present(dvara1)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara1(2:mza,iwnud) = dvara1(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara2)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara2(2:mza,iwnud) = dvara2(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara3)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara3(2:mza,iwnud) = dvara3(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara4)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara4(2:mza,iwnud) = dvara4(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara5)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara5(2:mza,iwnud) = dvara5(2:mza,iwnud) + vctr1(2:mza)
        endif

        if (present(dvara6)) then
           call MPI_Unpack(recv_wnud(jrecv)%buff,recv_wnud(jrecv)%nbytes,ipos, &
                           vctr1,mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

           dvara6(2:mza,iwnud) = dvara6(2:mza,iwnud) + vctr1(2:mza)
        endif

     enddo

     call MPI_Irecv(recv_wnud(jrecv)%buff, recv_wnud(jrecv)%nbytes, MPI_PACKED, &
                    recv_wnud(jrecv)%iremote, itagwnud, MPI_COMM_WORLD,         &
                    ireqr_wnud(jrecv,inext_wnud), ierr)
  enddo

#endif

end subroutine mpi_recv_wnud

!=============================================================================

subroutine olam_mpi_nud_stop()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,  only: olam_mpi_barrier
  use mem_nudge, only: nudflag, nudnxp
  use misc_coms, only: iparallel

  implicit none

#ifdef OLAM_MPI

  integer :: ierr, jrecv, ii, jsend
  logical :: flags

  if (iparallel == 1 .and. nudflag > 0 .and. nudnxp > 0) then

     ! First make sure all processes get here so that all sends are finished

     call olam_mpi_barrier()

     ! Cancel any pending communication

     do ii = 1, 2

        if (nudflag > 0 .and. nudnxp > 0) then

           do jrecv = 1, nrecvs_wnud
              if (ireqr_wnud(jrecv,ii) /= MPI_REQUEST_NULL) then
                 call MPI_Cancel(ireqr_wnud(jrecv,ii), ierr)
              endif
           enddo

           do jsend = 1, nsends_wnud
              if (ireqs_wnud(jsend,ii) /= MPI_REQUEST_NULL) then
                 call MPI_Cancel(ireqs_wnud(jsend,ii), ierr)
              endif
           enddo

        endif

     enddo

     ! Test that all communication requests have been completed or cancelled

     call olam_mpi_barrier()

     call MPI_Testall(nrecvs_wnud, ireqr_wnud(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nrecvs_wnud, ireqr_wnud(:,2), flags, MPI_STATUSES_IGNORE, ierr)

     call MPI_Testall(nsends_wnud, ireqs_wnud(:,1), flags, MPI_STATUSES_IGNORE, ierr)
     call MPI_Testall(nsends_wnud, ireqs_wnud(:,2), flags, MPI_STATUSES_IGNORE, ierr)

  endif

  ! Deallocate unused arrays

  if (allocated(send_wnud )) deallocate(send_wnud)
  if (allocated(recv_wnud )) deallocate(recv_wnud)

  if (allocated(ireqs_wnud)) deallocate(ireqs_wnud)
  if (allocated(ireqr_wnud)) deallocate(ireqr_wnud)

#endif

end subroutine olam_mpi_nud_stop

!=============================================================================

End Module olam_mpi_nud

