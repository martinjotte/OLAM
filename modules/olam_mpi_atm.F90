!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
Module olam_mpi_atm

Contains

subroutine olam_mpi_init()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para,  only: mgroupsize, myrank
  implicit none

#ifdef OLAM_MPI

  integer :: ierr

! Initialize MPI and determine process groupsize and myrank
  
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,mgroupsize,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)

#else

  mgroupsize = 1
  myrank     = 0

#endif

end subroutine olam_mpi_init

!===============================================================================

subroutine alloc_mpi_sndrcv_bufs()

  use mem_para,  only: mgroupsize,                            &
                       send_u, recv_u, send_uf, recv_uf,      &
                       send_v, recv_v, send_vf, recv_vf,      &
                       send_w, recv_w,                        &
                       send_wl, send_wlf, send_ws, send_wsf,  &
                       recv_wl, recv_wlf, recv_ws, recv_wsf

  use misc_coms, only: iparallel, meshtype
  implicit none

  if (iparallel == 0) return

! Allocate send and recv tables

  if (meshtype == 1) then

     allocate(send_u (mgroupsize))
     allocate(recv_u (mgroupsize))
     allocate(send_uf(mgroupsize))
     allocate(recv_uf(mgroupsize))

  else

     allocate(send_v (mgroupsize))
     allocate(recv_v (mgroupsize))
     allocate(send_vf(mgroupsize))
     allocate(recv_vf(mgroupsize))

  endif

  allocate(send_w (mgroupsize))
  allocate(recv_w (mgroupsize))

  allocate(send_wl (mgroupsize))
  allocate(send_wlf(mgroupsize))
  allocate(send_ws (mgroupsize))
  allocate(send_wsf(mgroupsize))

  allocate(recv_wl (mgroupsize))
  allocate(recv_wlf(mgroupsize))
  allocate(recv_ws (mgroupsize))
  allocate(recv_wsf(mgroupsize))

end subroutine alloc_mpi_sndrcv_bufs

!===============================================================================

subroutine olam_mpi_finalize()

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_para, only: send_u, recv_u, send_uf, recv_uf,      &
                      send_v, recv_v, send_vf, recv_vf,      &
                      send_w, recv_w,                        &
                      send_wl, send_wlf, send_ws, send_wsf,  &
                      recv_wl, recv_wlf, recv_ws, recv_wsf

  implicit none

#ifdef OLAM_MPI

  integer :: ierr

  call MPI_Finalize(ierr)
  
  if (allocated(send_u))  deallocate(send_u)
  if (allocated(recv_u))  deallocate(recv_u)
  if (allocated(send_uf)) deallocate(send_uf)
  if (allocated(recv_uf)) deallocate(recv_uf)

  if (allocated(send_v))  deallocate(send_v)
  if (allocated(recv_v))  deallocate(recv_v)
  if (allocated(send_vf)) deallocate(send_vf)
  if (allocated(recv_vf)) deallocate(recv_vf)

  if (allocated(send_w))  deallocate(send_w)
  if (allocated(recv_w))  deallocate(recv_w)

  if (allocated(send_wl))   deallocate(send_wl)
  if (allocated(send_wlf))  deallocate(send_wlf)
  if (allocated(send_ws))   deallocate(send_ws)
  if (allocated(send_wsf))  deallocate(send_wsf)

  if (allocated(recv_wl))  deallocate(recv_wl)
  if (allocated(recv_wlf)) deallocate(recv_wlf)
  if (allocated(recv_ws))  deallocate(recv_ws)
  if (allocated(recv_wsf)) deallocate(recv_wsf)

#endif

end subroutine olam_mpi_finalize

!==================================================================

subroutine olam_alloc_mpi(mza, mrls)

#ifdef OLAM_MPI
  use mpi
#endif

  use mem_ijtabs, only: jtab_u, jtab_v, jtab_w, mloops_u, mloops_v, mloops_w
  use mem_para,   only: myrank, nrecvs_u, nrecvs_v, nrecvs_w,     &
                        nsends_u, nsends_v, nsends_w,             &
                        recv_u, recv_v, recv_w, recv_uf, recv_vf, &
                        send_u, send_v, send_w, send_uf, send_vf
  use misc_coms,  only: io6, meshtype
  use var_tables, only: nvar_par

  implicit none

  integer, intent(in) :: mza
  integer, intent(in) :: mrls

#ifdef OLAM_MPI

  integer :: nbytes_int
  integer :: nbytes_real
  integer :: nbytes_real8

  integer :: nbytes_per_iu
  integer :: nbytes_per_iuf

  integer :: nbytes_per_iv
  integer :: nbytes_per_ivf

  integer :: nbytes_per_iw

  integer :: itag10  =  10
  integer :: itag110 = 110

  integer :: ierr
  integer :: jsend
  integer :: jrecv
  integer :: jtmp

  integer :: nupts, nvpts, nwpts
  integer :: mrl

  integer              :: sbuf(mrls+1)
  integer, allocatable :: rbuf(:,:)

  integer              :: wsbuf(mrls)
  integer, allocatable :: wrbuf(:,:)

! Allocate send buffers

  call MPI_Pack_size(1,MPI_INTEGER,MPI_COMM_WORLD,nbytes_int  ,ierr)
  call MPI_Pack_size(1,MPI_REAL   ,MPI_COMM_WORLD,nbytes_real ,ierr)
  call MPI_Pack_size(1,MPI_REAL8  ,MPI_COMM_WORLD,nbytes_real8,ierr)

  if (meshtype == 1) then

! Post U receives

     allocate(rbuf(mrls+1,nrecvs_u(1)))

     do jrecv = 1,nrecvs_u(1)
        call MPI_IRecv(rbuf(1,jrecv), mrls+1, MPI_INTEGER, recv_u(jrecv)%iremote, &
             itag10, MPI_COMM_WORLD, recv_u(jrecv)%irequest, ierr)
     enddo

  else

! Post V receives

     allocate(rbuf(mrls+1,nrecvs_v(1)))

     do jrecv = 1,nrecvs_v(1)
        call MPI_IRecv(rbuf(1,jrecv), mrls+1, MPI_INTEGER, recv_v(jrecv)%iremote, &
             itag10, MPI_COMM_WORLD, recv_v(jrecv)%irequest, ierr)
     enddo

  endif

! Post W receives

  allocate(wrbuf(mrls,nrecvs_w(1)))

  do jrecv = 1,nrecvs_w(1)

     call MPI_IRecv(wrbuf(1,jrecv), mrls, MPI_INTEGER, recv_w(jrecv)%iremote, &
          itag110, MPI_COMM_WORLD, recv_w(jrecv)%irequest, ierr)

  enddo

  if (meshtype == 1) then

! If triangular grid, determine number of bytes to send per IU column

     nbytes_per_iu = nbytes_int  &
          + mza * 2 * nbytes_real

     nbytes_per_iuf = nbytes_int  &
          + mza * 3 * nbytes_real

! Loop over all U sends for mrl = 1

     do jsend = 1,nsends_u(1)

! Determine size of send_u buffer for mrl = 1

        send_u(jsend)%nbytes = nbytes_int   &
             + nbytes_per_iu * jtab_u(mloops_u+jsend)%jend(1)

        send_uf(jsend)%nbytes = nbytes_int   &
             + nbytes_per_iuf * jtab_u(mloops_u+jsend)%jend(1)

! Allocate buffer

        allocate(send_u(jsend)%buff(send_u(jsend)%nbytes))
        allocate(send_uf(jsend)%buff(send_uf(jsend)%nbytes))

! Send buffer sizes to receive ranks

        sbuf(1) = send_u(jsend)%nbytes
        sbuf(2) = send_uf(jsend)%nbytes

        do mrl = 2,mrls
           sbuf(mrl+1) = jtab_u(mloops_u+jsend)%jend(mrl)

! If at least 1 U point needs to be sent to current remote rank for 
! current mrl, increase nsends_u(mrl) by 1.

           if (jtab_u(mloops_u+jsend)%jend(mrl) > 0) &
                nsends_u(mrl) = nsends_u(mrl) + 1

        enddo

        call MPI_Send(sbuf, mrls+1, MPI_INTEGER, send_u(jsend)%iremote, itag10, &
             MPI_COMM_WORLD, ierr)

     enddo

  else

! If hexagonal grid, determine number of bytes to send per IV column

     nbytes_per_iv = nbytes_int  &
          + mza * 2 * nbytes_real

     nbytes_per_ivf = nbytes_int  &
          + mza * (4 * nbytes_real + nbytes_real8)

! Loop over all V sends for mrl = 1

     do jsend = 1,nsends_v(1)

! Determine size of send_v buffer for mrl = 1

        send_v(jsend)%nbytes = nbytes_int   &
             + nbytes_per_iv * jtab_v(mloops_v+jsend)%jend(1)

        send_vf(jsend)%nbytes = nbytes_int   &
             + nbytes_per_ivf * jtab_v(mloops_v+jsend)%jend(1)

! Allocate buffer

        allocate(send_v(jsend)%buff(send_v(jsend)%nbytes))
        allocate(send_vf(jsend)%buff(send_vf(jsend)%nbytes))

! Send buffer sizes to receive ranks

        sbuf(1) = send_v(jsend)%nbytes
        sbuf(2) = send_vf(jsend)%nbytes

        do mrl = 2,mrls
           sbuf(mrl+1) = jtab_v(mloops_v+jsend)%jend(mrl)

! If at least 1 V point needs to be sent to current remote rank for 
! current mrl, increase nsends_v(mrl) by 1.

           if (jtab_v(mloops_v+jsend)%jend(mrl) > 0) &
                nsends_v(mrl) = nsends_v(mrl) + 1

        enddo

        call MPI_Send(sbuf, mrls+1, MPI_INTEGER, send_v(jsend)%iremote, itag10, &
             MPI_COMM_WORLD,ierr)

     enddo

  endif

! Determine number of bytes to send per IW column

  nbytes_per_iw = nbytes_int                                     &
       + mza * max(3 * nbytes_real8 + 2 * nbytes_real,  &
       nvar_par * nbytes_real)

! Loop over all W sends for mrl = 1

  do jsend = 1,nsends_w(1)

! Determine size of send_w buffer for mrl = 1

     send_w(jsend)%nbytes = nbytes_int                               &
          + nbytes_per_iw * jtab_w(mloops_w+jsend)%jend(1)

! Allocate buffer

     allocate(send_w(jsend)%buff(send_w(jsend)%nbytes))

! Send buffer sizes to receive ranks

     wsbuf(1) = send_w(jsend)%nbytes

     do mrl = 2,mrls
        wsbuf(mrl) = jtab_w(mloops_w+jsend)%jend(mrl)

! If at least 1 W point needs to be sent to current remote rank for 
! current mrl, increase nsends_w(mrl) by 1.

        if (jtab_w(mloops_w+jsend)%jend(mrl) > 0) nsends_w(mrl) = nsends_w(mrl) + 1

     enddo

     call MPI_Send(wsbuf, mrls, MPI_INTEGER, send_w(jsend)%iremote, itag110, &
          MPI_COMM_WORLD, ierr)

  enddo

  if (meshtype == 1) then

! If triangular grid, loop over all U receives for mrl = 1

     do jtmp = 1, nrecvs_u(1)

! Get completed recv_u buffer sizes

        call MPI_Waitany(nrecvs_u(1), recv_u(1:nrecvs_u(1))%irequest, jrecv, &
             MPI_STATUS_IGNORE, ierr)

        recv_u(jrecv)%nbytes  = rbuf(1,jrecv)
        recv_uf(jrecv)%nbytes = rbuf(2,jrecv)

! Allocate recv_u buffers

        allocate(recv_u(jrecv)%buff(recv_u(jrecv)%nbytes))
        allocate(recv_uf(jrecv)%buff(recv_uf(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 U point needs to be 
! received from current remote rank for each mrl, increase nrecvs_u(mrl) by 1.

        do mrl = 2,mrls
           if (rbuf(mrl+1,jrecv) > 0) nrecvs_u(mrl) = nrecvs_u(mrl) + 1
        enddo

     enddo

  else

! If hexagonal grid, loop over all V receives for mrl = 1

     do jtmp = 1,nrecvs_v(1)

! Get recv_v buffer sizes from any node

        call MPI_Waitany(nrecvs_v(1), recv_v(1:nrecvs_v(1))%irequest, jrecv, &
             MPI_STATUS_IGNORE, ierr)

        recv_v(jrecv)%nbytes  = rbuf(1,jrecv)
        recv_vf(jrecv)%nbytes = rbuf(2,jrecv)

! Allocate recv_v buffers

        allocate(recv_v(jrecv)%buff(recv_v(jrecv)%nbytes))
        allocate(recv_vf(jrecv)%buff(recv_vf(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 V point needs to be 
! received from current remote rank for each mrl, increase nrecvs_v(mrl) by 1.

        do mrl = 2,mrls
           if (rbuf(mrl+1,jrecv) > 0) nrecvs_v(mrl) = nrecvs_v(mrl) + 1
        enddo

     enddo

  endif

! Loop over all W receives for mrl = 1

  do jtmp = 1,nrecvs_w(1)

! Get recv_w buffer sizes

     call MPI_Waitany(nrecvs_w(1), recv_w(1:nrecvs_w(1))%irequest, jrecv, &
          MPI_STATUS_IGNORE, ierr)

! Allocate recv_w buffers

     recv_w(jrecv)%nbytes = wrbuf(1,jrecv)
     allocate(recv_w(jrecv)%buff(recv_w(jrecv)%nbytes))

! Loop over all mrl values greater than 1. If at least 1 W point needs to be 
! received from current remote rank for each mrl, increase nrecvs_w(mrl) by 1.

     do mrl = 2,mrls
        if (wrbuf(mrl,jrecv) > 0) nrecvs_w(mrl) = nrecvs_w(mrl) + 1
     enddo

  enddo

#endif

  return
end subroutine olam_alloc_mpi

!===============================================================================

subroutine mpi_send_u(sendgroup,uc0,rpos,rneg)

! Subroutine to perform a parallel MPI send of a "U group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_basic,  only: umc,uc
use mem_para,   only: send_u, recv_u, nsends_u, nrecvs_u
use mem_ijtabs, only: itab_u, jtab_u, mrl_begs, istp, mloops_u

use mem_grid,   only: mza, mua
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: sendgroup

real, optional, intent(in) :: uc0 (mza,mua)
real, optional, intent(in) :: rpos(mza,mua)
real, optional, intent(in) :: rneg(mza,mua)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag1 = 1
integer :: mrl
integer :: j
integer :: iu
integer :: iuglobe

! Set MRL and return if mrl < 1 for this step

if (sendgroup == 'I') then
   mrl = 1
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_u(mrl)

   call MPI_Irecv(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,MPI_PACKED,  &
                  recv_u(jrecv)%iremote,itag1,MPI_COMM_WORLD,          &
                  recv_u(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_u(mrl)

   ipos = 0

   call MPI_Pack(jtab_u(mloops_u+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_u(mloops_u+jsend)%jend(mrl)
      iu = jtab_u(mloops_u+jsend)%iu(j)
      iuglobe = itab_u(iu)%iuglobe
!----------------------------------------------------------------
      call qsub('U',iu)

      call MPI_Pack(iuglobe,1,MPI_INTEGER,  &
         send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'L') then

         call MPI_Pack(uc0(1,iu),mza,MPI_REAL,  &
            send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'R') then

         call MPI_Pack(rpos(1,iu),mza,MPI_REAL,  &
            send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(rneg(1,iu),mza,MPI_REAL,  &
            send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'U' .or. sendgroup == 'I') then
      
         call MPI_Pack(umc(1,iu),mza,MPI_REAL,  &
            send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(uc(1,iu),mza,MPI_REAL,  &
            send_u(jsend)%buff,send_u(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
 
      endif

   enddo
   call rsub('Usend',mloops_u+jsend)

   call MPI_Isend(send_u(jsend)%buff,ipos,MPI_PACKED,          &
                  send_u(jsend)%iremote,itag1,MPI_COMM_WORLD,  &
                  send_u(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_u

!===============================================================================

subroutine mpi_send_v(sendgroup)

! Subroutine to perform a parallel MPI send of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_basic,  only: vmc,vc
use mem_para,   only: send_v, recv_v, nsends_v, nrecvs_v

use mem_ijtabs, only: itab_v, jtab_v, mrl_begs, istp, mloops_v
use mem_grid,   only: mza
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: sendgroup

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag1 = 1
integer :: mrl
integer :: j
integer :: iv
integer :: ivglobe

! Set MRL and return if mrl < 1 for this step

if (sendgroup == 'I') then
   mrl = 1
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_v(mrl)

   call MPI_Irecv(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,MPI_PACKED,  &
                  recv_v(jrecv)%iremote,itag1,MPI_COMM_WORLD,          &
                  recv_v(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_v(mrl)

   ipos = 0

   call MPI_Pack(jtab_v(mloops_v+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_v(mloops_v+jsend)%jend(mrl)
      iv = jtab_v(mloops_v+jsend)%iv(j)
      ivglobe = itab_v(iv)%ivglobe
!----------------------------------------------------------------
      call qsub('V',iv)

      call MPI_Pack(ivglobe,1,MPI_INTEGER,  &
         send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(vmc(1,iv),mza,MPI_REAL,  &
         send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(vc(1,iv),mza,MPI_REAL,  &
         send_v(jsend)%buff,send_v(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   enddo
   call rsub('Vsend',mloops_v+jsend)

   call MPI_Isend(send_v(jsend)%buff,ipos,MPI_PACKED,          &
                  send_v(jsend)%iremote,itag1,MPI_COMM_WORLD,  &
                  send_v(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_v

!===============================================================================

subroutine mpi_send_vf(hcnum_v,hcnum_w,hflux_t,umaru)

! Subroutine to perform a parallel MPI send of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_para,   only: send_vf, recv_vf, nsends_v, nrecvs_v
use mem_ijtabs, only: jtab_v, itab_v, mrl_begs, istp
use mem_grid,   only: mza, mva
use misc_coms,  only: io6
use mem_basic,  only: uc

implicit none

real, intent(in) :: hcnum_v(mza,mva)
real, intent(in) :: hcnum_w(mza,mva)
real, intent(in) :: hflux_t(mza,mva)

real(kind=8), intent(in) :: umaru  (mza,mva)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag2 = 2
integer :: mrl
integer :: j
integer :: iv
integer :: ivglobe

! Return if mrl < 1 for this step

mrl = mrl_begs(istp)
if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_v(mrl)

   call MPI_Irecv(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,MPI_PACKED,  &
                  recv_vf(jrecv)%iremote,itag2,MPI_COMM_WORLD,          &
                  recv_vf(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_v(mrl)

   ipos = 0

   call MPI_Pack(jtab_v(25+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_v(25+jsend)%jend(mrl)
      iv = jtab_v(25+jsend)%iv(j)
      ivglobe = itab_v(iv)%ivglobe
!----------------------------------------------------------------
      call qsub('V',iv)

      call MPI_Pack(ivglobe,1,MPI_INTEGER,  &
         send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(hcnum_v(1,iv),mza,MPI_REAL,  &
         send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(hcnum_w(1,iv),mza,MPI_REAL,  &
         send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(hflux_t(1,iv),mza,MPI_REAL,  &
         send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(umaru(1,iv),mza,MPI_REAL8,  &
         send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      call MPI_Pack(uc(1,iv),mza,MPI_REAL,  &
         send_vf(jsend)%buff,send_vf(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   enddo

   call rsub('Vsendf',25+jsend)

   call MPI_Isend(send_vf(jsend)%buff,ipos,MPI_PACKED,          &
                  send_vf(jsend)%iremote,itag2,MPI_COMM_WORLD,  &
                  send_vf(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_vf

!=============================================================================

subroutine mpi_send_w(sendgroup,thil0,wmc0,scp0,wmarw)

! Subroutine to perform a parallel MPI send of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_basic,  only: wmc,wc,thil,rho,press
use mem_turb,   only: hkm,vkm,vkm_sfc
use var_tables, only: nvar_par, vtab_r, num_scalar, nptonv
use mem_ijtabs, only: jtab_w, itab_w, mrl_begs, mrl_begl, istp, mloops_w
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use micro_coms, only: level
use mem_para,   only: nrecvs_w, nsends_w, recv_w, send_w

implicit none

character(1), intent(in) :: sendgroup

real, optional, intent(in) :: thil0(mza,mwa,1)
real, optional, intent(in) :: wmc0(mza,mwa)
real, optional, intent(in) :: scp0(mza,mwa,num_scalar)

real(kind=8), optional, intent(in) :: wmarw(mza,mwa)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar
integer :: itag3 = 3
integer :: mrl
integer :: i, j
integer :: iw
integer :: iwglobe

! Set MRL and return if mrl < 1 for this step

if (sendgroup == 'I') then
   mrl = 1
elseif (sendgroup == 'S') then
   mrl = mrl_begl(istp)
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

! Before we send anything, post the receives

do jrecv = 1,nrecvs_w(mrl)

   call MPI_Irecv(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,MPI_PACKED,  &
                  recv_w(jrecv)%iremote,itag3,MPI_COMM_WORLD,          &
                  recv_w(jrecv)%irequest,ierr                          )

enddo

! Now we can actually go on to sending the stuff

do jsend = 1,nsends_w(mrl)

   ipos = 0

   call MPI_Pack(jtab_w(mloops_w+jsend)%jend(mrl),1,MPI_INTEGER,  &
      send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,jtab_w(mloops_w+jsend)%jend(mrl)
      iw = jtab_w(mloops_w+jsend)%iw(j)
      iwglobe = itab_w(iw)%iwglobe
!----------------------------------------------------------------
      call qsub('W',iw)

      call MPI_Pack(iwglobe,1,MPI_INTEGER,  &
         send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      if (sendgroup == 'S') then

         do i = 1, nvar_par
            ivar = nptonv(i)

            call MPI_Pack(vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL,  &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         enddo

      elseif (sendgroup == 'L') then
      
         if (present(thil0)) then

            call MPI_Pack(thil0(1,iw,1),mza,MPI_REAL,  &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         elseif (present(wmc0)) then

            call MPI_Pack(wmc0(1,iw),mza,MPI_REAL,  &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         elseif (present(scp0)) then

            do ivar = 1,num_scalar

               call MPI_Pack(scp0(1,iw,ivar),mza,MPI_REAL,  &
                  send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

            enddo
            
         endif

      elseif (sendgroup == 'I') then

         call MPI_Pack(press(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(rho(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(wc(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(thil(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'K') then

         call MPI_Pack(hkm(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(vkm(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(vkm_sfc(iw),1,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

      elseif (sendgroup == 'P') then

         call MPI_Pack(wmc(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(press(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(rho(1,iw),mza,MPI_REAL8,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         if (present(wmarw)) then
            call MPI_Pack(wmarw(1,iw),mza,MPI_REAL8,  &
                 send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
         endif

      elseif (sendgroup == 'T') then

!        call MPI_Pack(wc(1,iw),mza,MPI_REAL,  &
!           send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         call MPI_Pack(thil(1,iw),mza,MPI_REAL,  &
            send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)

         if (level == 3) then
            call MPI_Pack(rho(1,iw),mza,MPI_REAL8,  &
               send_w(jsend)%buff,send_w(jsend)%nbytes,ipos,MPI_COMM_WORLD,ierr)
         endif
 
      endif

   enddo
   call rsub('Wsend',mloops_w+jsend)

   call MPI_Isend(send_w(jsend)%buff,ipos,MPI_PACKED,          &
                  send_w(jsend)%iremote,itag3,MPI_COMM_WORLD,  &
                  send_w(jsend)%irequest,ierr                  )

enddo

#endif

return
end subroutine mpi_send_w

!=============================================================================

subroutine mpi_recv_u(recvgroup,uc0,rpos,rneg)

! Subroutine to perform a parallel MPI receive of a "U group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_basic,  only: umc,uc
use mem_para,   only: send_u, recv_u, nsends_u, nrecvs_u, mgroupsize
use mem_ijtabs, only: itabg_u, mrl_begs, istp, mloops_u
use mem_grid,   only: mza, mua
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: recvgroup

real, optional, intent(inout) :: uc0 (mza,mua)
real, optional, intent(inout) :: rpos(mza,mua)
real, optional, intent(inout) :: rneg(mza,mua)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nupts
integer :: mrl
integer :: j
integer :: iu
integer :: iuglobe

! Set MRL and return if mrl < 1 for this step

if (recvgroup == 'I') then
   mrl = 1
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_u(mrl)

   call MPI_Waitany(nrecvs_u(mrl), recv_u(1:nrecvs_u(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
      nupts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nupts
      call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
         iuglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iu = itabg_u(iuglobe)%iu_myrank
!----------------------------------------------------------------
      call qsub('U',iu)

      if (recvgroup == 'L') then

         call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
            uc0(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'R') then

         call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
            rpos(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
            rneg(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'U' .or. recvgroup == 'I') then

         call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
            umc(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_u(jrecv)%buff,recv_u(jrecv)%nbytes,ipos,  &
            uc(1,iu),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      endif

   enddo
   call rsub('Urecv',mloops_u+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated
call MPI_Waitall(nsends_u(mrl), send_u(1:nsends_u(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_u

!=============================================================================

subroutine mpi_recv_v(recvgroup)

! Subroutine to perform a parallel MPI receive of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_basic,  only: vmc,vc
use mem_para,   only: send_v, recv_v, nsends_v, nrecvs_v, mgroupsize
use mem_ijtabs, only: itabg_v, mrl_begs, istp, mloops_v
use mem_grid,   only: mza
use misc_coms,  only: io6

implicit none

character(1), intent(in) :: recvgroup

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nvpts
integer :: mrl
integer :: j
integer :: iv
integer :: ivglobe

! Set MRL and return if mrl < 1 for this step

if (recvgroup == 'I') then
   mrl = 1
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_v(mrl)

   call MPI_Waitany(nrecvs_v(mrl), recv_v(1:nrecvs_v(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos,  &
      nvpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nvpts
      call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos,  &
         ivglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iv = itabg_v(ivglobe)%iv_myrank
!----------------------------------------------------------------
      call qsub('V',iv)

      call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos,  &
         vmc(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_v(jrecv)%buff,recv_v(jrecv)%nbytes,ipos,  &
         vc(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

   enddo
   call rsub('Vrecv',mloops_v+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated
call MPI_Waitall(nsends_v(mrl), send_v(1:nsends_v(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_v

!=============================================================================

subroutine mpi_recv_vf(hcnum_v,hcnum_w,hflux_t,umaru)

! Subroutine to perform a parallel MPI receive of a "V group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_para,   only: send_vf, recv_vf, nsends_v, nrecvs_v, mgroupsize
use mem_ijtabs, only: mrl_begs, itabg_v, istp
use mem_grid,   only: mza, mva
use misc_coms,  only: io6
use mem_basic,  only: uc

implicit none

real, intent(inout) :: hcnum_v(mza,mva)
real, intent(inout) :: hcnum_w(mza,mva)
real, intent(inout) :: hflux_t(mza,mva)

real(kind=8), intent(inout) :: umaru(mza,mva)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nvpts
integer :: mrl
integer :: j
integer :: iv
integer :: ivglobe

! Return if mrl < 1 for this step

mrl = mrl_begs(istp)
if (mrl < 1) return

!  Make sure sends are all finished and de-allocated

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_v(mrl)

   call MPI_Waitany(nrecvs_v(mrl), recv_vf(1:nrecvs_v(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
      nvpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nvpts
      call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
         ivglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iv = itabg_v(ivglobe)%iv_myrank
!----------------------------------------------------------------
      call qsub('V',iv)

      call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
         hcnum_v(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
         hcnum_w(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
         hflux_t(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
         umaru(1,iv),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

      call MPI_Unpack(recv_vf(jrecv)%buff,recv_vf(jrecv)%nbytes,ipos,  &
         uc(1,iv),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

   enddo
   call rsub('Vrecv',25+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated
call MPI_Waitall(nsends_v(mrl), send_vf(1:nsends_v(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_vf

!=============================================================================

subroutine mpi_recv_w(recvgroup,thil0,wmc0,scp0,wmarw)

! Subroutine to perform a parallel MPI receive of a "W group" or "S group"
! of field variables

#ifdef OLAM_MPI
  use mpi
#endif

use mem_basic,  only: wmc,wc,thil,rho,press
use mem_turb,   only: hkm,vkm,vkm_sfc
use var_tables, only: vtab_r, nvar_par, num_scalar, nptonv
use mem_para,   only: nrecvs_w, nsends_w, recv_w, send_w, mgroupsize
use mem_ijtabs, only: itabg_w, mrl_begs, mrl_begl, istp, mloops_w
use mem_grid,   only: mza, mwa
use misc_coms,  only: io6
use micro_coms, only: level

implicit none

character(1), intent(in) :: recvgroup

real, optional, intent(inout) :: thil0(mza,mwa,1)
real, optional, intent(inout) :: wmc0(mza,mwa)
real, optional, intent(inout) :: scp0(mza,mwa,num_scalar)

real(kind=8), optional, intent(inout) :: wmarw(mza,mwa)

#ifdef OLAM_MPI

integer :: ierr,ipos
integer :: jrecv,jsend,ivar,jtmp
integer :: nwpts
integer :: mrl
integer :: i, j
integer :: iw
integer :: iwglobe

! Set MRL and return if mrl < 1 for this step

if (recvgroup == 'I') then
   mrl = 1
elseif (recvgroup == 'S') then
   mrl = mrl_begl(istp)
else
   mrl = mrl_begs(istp)
endif

if (mrl < 1) return

!  Now, let's wait on our receives

do jtmp = 1,nrecvs_w(mrl)

   call MPI_Waitany(nrecvs_w(mrl), recv_w(1:nrecvs_w(mrl))%irequest, jrecv, &
        MPI_STATUS_IGNORE, ierr)

!  We got all our stuff.  Now unpack it into appropriate space.

   ipos = 0

   call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
      nwpts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

   call psub()
!----------------------------------------------------------------
   do j = 1,nwpts
      call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
         iwglobe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

      iw = itabg_w(iwglobe)%iw_myrank
!----------------------------------------------------------------
      call qsub('W',iw)

      if (recvgroup == 'S') then

         do i = 1, nvar_par
            ivar = nptonv(i)

            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
               vtab_r(ivar)%rvar2_p(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         enddo

      elseif (recvgroup == 'L') then

         if (present(thil0)) then

            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
               thil0(1,iw,1),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         elseif (present(wmc0)) then

            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
               wmc0(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         elseif (present(scp0)) then
         
            do ivar = 1,num_scalar

               call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
                  scp0(1,iw,ivar),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

            enddo

         endif

      elseif (recvgroup == 'I') then

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            press(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            rho(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            wc(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            thil(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'K') then

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            hkm(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            vkm(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            vkm_sfc(iw),1,MPI_REAL,MPI_COMM_WORLD,ierr)

      elseif (recvgroup == 'P') then
            
         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            wmc(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            press(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            rho(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)

         if (present(wmarw)) then
            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
                 wmarw(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
         endif

      elseif (recvgroup == 'T') then

!        call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
!           wc(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
            thil(1,iw),mza,MPI_REAL,MPI_COMM_WORLD,ierr)

         if (level == 3) then
            call MPI_Unpack(recv_w(jrecv)%buff,recv_w(jrecv)%nbytes,ipos,  &
               rho(1,iw),mza,MPI_REAL8,MPI_COMM_WORLD,ierr)
         endif

      endif

   enddo
   call rsub('Wrecv',mloops_w+jrecv)

enddo

! Make sure all of our sends are finished and de-allocated

call MPI_Waitall(nsends_w(mrl), send_w(1:nsends_w(mrl))%irequest, &
     MPI_STATUSES_IGNORE, ierr)

#endif

return
end subroutine mpi_recv_w

!=============================================================================

subroutine olam_stop(message)

#ifdef OLAM_MPI
  use mpi
#endif

use mem_para,  only: myrank

implicit none

#ifdef OLAM_MPI
integer :: ierr
#endif

character(*), intent(in) :: message

#ifdef MPI  
  write(*,'(A,I0,A)') "Node ", myrank, ":"
  write(*,'(A)') "STOP "//message
  call mpi_abort(MPI_COMM_WORLD,1,ier)
  stop
#else
  write(*,*) "STOPPING: "//message
  stop
#endif

end subroutine olam_stop

End Module olam_mpi_atm
