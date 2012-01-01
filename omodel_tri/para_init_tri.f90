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
subroutine para_init()

use misc_coms,  only: io6, meshtype

use mem_ijtabs, only: itab_m,      itab_u,      itab_w,      &
                      itab_m_vars, itab_u_vars, itab_w_vars, &
                      itabg_m,     itabg_u,     itabg_w,     &
                      itab_m_pd,   itab_u_pd,   itab_w_pd,   &
                      alloc_itabs, mrls

use mem_grid,   only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, &
                      alloc_gridz, alloc_xyzem, alloc_xyzew, &
                      alloc_grid1, alloc_grid2


use mem_para,   only: mgroupsize, myrank, &
                      send_u, recv_u, send_v, recv_v, send_w, recv_w, &
                      nsends_u, nsends_v, nsends_w, &
                      nrecvs_u, nrecvs_v, nrecvs_w

use mem_sflux,  only: nseaflux,  mseaflux,  seaflux,  seafluxg,  &
                      nlandflux, mlandflux, landflux, landfluxg, &
                      flux_vars

use sea_coms,   only: nws

use leaf_coms,  only: nwl, isfcl

use mem_sea,    only: itabg_ws

use mem_leaf,   only: itabg_wl

implicit none

integer :: j,k,imn,iun,ivn,iwn
integer :: im,iu,iv,iw
integer :: itopm,iup,ivp,iwp
integer :: isf,ilf,iws,iwl
integer :: npoly

integer :: im_myrank = 1 ! Counter for M points to be included on this rank
integer :: iu_myrank = 1 ! Counter for U points to be included on this rank
integer :: iv_myrank = 1 ! Counter for V points to be included on this rank
integer :: iw_myrank = 1 ! Counter for W points to be included on this rank

! Automatic arrays

logical :: myrankflag_m(nma) ! Flag for M points existing on this rank
logical :: myrankflag_u(nua) ! Flag for U points existing on this rank
logical :: myrankflag_w(nwa) ! Flag for W points existing on this rank

logical :: seaflag(nws)
logical :: landflag(nwl)

! Temporary datatypes

type(flux_vars), allocatable :: landflux_temp(:)
type(flux_vars), allocatable ::  seaflux_temp(:)

integer :: ierr

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_u(mrls)) ; nsends_u(1:mrls) = 0
allocate (nsends_w(mrls)) ; nsends_w(1:mrls) = 0
allocate (nrecvs_u(mrls)) ; nrecvs_u(1:mrls) = 0
allocate (nrecvs_w(mrls)) ; nrecvs_w(1:mrls) = 0

! Initialize myrank flag arrays to .false.

if (isfcl == 1) then
   landflag(:) = .false.
   seaflag (:) = .false.
endif

myrankflag_m(:) = .false.
myrankflag_u(:) = .false.
myrankflag_w(:) = .false.

! Loop over all U points, and for each whose assigned irank is equal to myrank,
! flag all U and W points in its computational stencil for inclusion on this
! rank, excluding IUP and IWP.

do iu = 2,nua

   if (itabg_u(iu)%irank == myrank) then

      myrankflag_u(iu) = .true.

      myrankflag_u( itab_u_pd(iu)%iu(1:12) ) = .true.
      myrankflag_w( itab_u_pd(iu)%iw(1:6)  ) = .true.

   endif
enddo

! Loop over all W points, and for each whose assigned irank is equal to myrank,
! flag all U and W points in its computational stencil for inclusion on this
! rank, excluding IUP and IWP.

do iw = 2,nwa

   if (itabg_w(iw)%irank == myrank) then

      myrankflag_w(iw) = .true.

! The standard computational stencil of mem_ijtabs

      myrankflag_w( itab_w_pd(iw)%iw(1:9) ) = .true.
      myrankflag_u( itab_w_pd(iw)%iu(1:9) ) = .true.

! Special for Zalesak monotonic advection (each W point needs the 9 
! U points from each of the 3 bordering triangles available)

      do j=1,3
         iwn = itab_w_pd(iw)%iw(j)
         myrankflag_u( itab_w_pd(iwn)%iu(1:9) ) = .true.
      enddo

   endif
enddo

! Loop over all U points, and for each that has been flagged for inclusion
! on this rank, flag both its M points for inclusion on this rank.
! Count U points also.

do iu = 2,nua

   if (myrankflag_u(iu)) then

      myrankflag_m( itab_u_pd(iu)%im(1:2) ) = .true.
      iu_myrank = iu_myrank + 1

   endif
enddo

! Loop over all M and W points and count the ones that have been flagged
! for inclusion on this rank.

do im = 2,nma
   if (myrankflag_m(im)) then
      im_myrank = im_myrank + 1
   endif
enddo

do iw = 2,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1
   endif
enddo

! Set mma, mua, mwa values for this rank

mma = im_myrank
mua = iu_myrank
mva = mua
mwa = iw_myrank


! Allocate grid structure variables

call alloc_gridz()
call alloc_itabs(meshtype,mma,mua,mva,mwa)
call alloc_xyzem(mma)
call alloc_xyzew(mwa)
call alloc_grid1(meshtype, mma, mua, mva, mwa)
call alloc_grid2(meshtype, mma, mua, mva, mwa)

! Reset point counts to 1

im_myrank = 1
iu_myrank = 1
iw_myrank = 1

! Store new myrank M, U, W indices in itabg data structures

do im = 1,nma
   if (myrankflag_m(im)) then
      im_myrank = im_myrank + 1
      
      itabg_m(im)%im_myrank = im_myrank      
   endif
enddo

do iu = 1,nua
   if (myrankflag_u(iu)) then
      iu_myrank = iu_myrank + 1
      
      itabg_u(iu)%iu_myrank = iu_myrank
      
! Fill itabg_u(iu)%iu_myrank value for IUP point

      iup = itab_u_pd(iu)%iup
      itabg_u(iup)%iu_myrank = iu_myrank
   endif
enddo

do iw = 1,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1

      itabg_w(iw)%iw_myrank = iw_myrank      

! Fill itabg_w(iw)%iw_myrank value for IWP point

      iwp = itab_w_pd(iw)%iwp
      itabg_w(iwp)%iw_myrank = iw_myrank
   endif
enddo


! Defining global index of local points (will be used on gridfile_read)

do im = 1, nma
   if (myrankflag_m(im)) itab_m(itabg_m(im)%im_myrank)%imglobe = im
enddo

do iu = 1, nua
   if (myrankflag_u(iu)) itab_u(itabg_u(iu)%iu_myrank)%iuglobe = iu
enddo

do iw = 1, nwa
   if (myrankflag_w(iw)) itab_w(itabg_w(iw)%iw_myrank)%iwglobe = iw
enddo

! Reading the local grid structure

call gridfile_read()

!!!! ITAB_? POPULATION

! put itab_m on place

do im = 1,nma
   if (myrankflag_m(im)) then
      npoly = itab_m_pd(im)%npoly
   
      im_myrank = itabg_m(im)%im_myrank

! Reset IM neighbor indices to 1

      itab_m(im_myrank)%itopm = 1
      itab_m(im_myrank)%iw(1:npoly) = 1
      itab_m(im_myrank)%iu(1:npoly) = 1

! Global indices of neighbors of IM
! Set indices of neighbors of IM that are present on this rank

      itopm = itab_m_pd(im)%itopm
      if (myrankflag_m(itopm)) itab_m(im_myrank)%itopm = itabg_m(itopm)%im_myrank

      do j = 1,npoly
         iu = itab_m_pd(im)%iu(j)
         iw = itab_m_pd(im)%iw(j)
      
         if (myrankflag_u(iu)) itab_m(im_myrank)%iu(j) = itabg_u(iu)%iu_myrank
         if (myrankflag_w(iw)) itab_m(im_myrank)%iw(j) = itabg_w(iw)%iw_myrank
      enddo
   endif
enddo

! put itab_u on place

do iu = 1,nua
   if (myrankflag_u(iu)) then
   
! Look up global IUP index and subdomain IU index

      iup = itab_u_pd(iu)%iup
      iu_myrank = itabg_u(iu)%iu_myrank

! Next, redefine some individual itab_u members

      itab_u(iu_myrank)%irank   = itabg_u(iu)%irank

! Check if this U point is primary on a remote rank

      if (itab_u(iu_myrank)%irank /= myrank) then      

! Turn off some loop flags for these points (some will be turned back on later)

         call uloops('n',iu_myrank,-7,-8,-12,-13,-16,-21,-22,-23, 0, 0)

! Turn off LBC copy (n/a for global domain) if IUP point is on remote node 

         if (.not. myrankflag_u(iup))  &
            call uloops('n',iu_myrank,-9,-18,0,0,0,0,0,0,0,0)

      endif

! Reset IU neighbor indices to 1

      itab_u(iu_myrank)%iup = 1

      itab_u(iu_myrank)%im(1:2)  = 1
      itab_u(iu_myrank)%iu(1:12) = 1
      itab_u(iu_myrank)%iw(1:6)  = 1

! Set indices of neighbors of IU that are present on this rank

      if (myrankflag_u(iup)) itab_u(iu_myrank)%iup = itabg_u(iup)%iu_myrank

      do j = 1,2
         imn = itab_u_pd(iu)%im(j)
         if (myrankflag_m(imn)) itab_u(iu_myrank)%im(j) = itabg_m(imn)%im_myrank
      enddo

      do j = 1,12
         iun = itab_u_pd(iu)%iu(j)
         if (myrankflag_u(iun)) itab_u(iu_myrank)%iu(j) = itabg_u(iun)%iu_myrank
      enddo

      do j = 1,6
         iwn = itab_u_pd(iu)%iw(j)
         if (myrankflag_w(iwn)) itab_u(iu_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
      enddo
   endif
enddo

! put itab_w on place

do iw = 1,nwa
   if (myrankflag_w(iw)) then

! Look up global IWP index and subdomain IW index
   
      iwp = itab_w_pd(iw)%iwp
      iw_myrank = itabg_w(iw)%iw_myrank

! Next, redefine some individual itab_w members

      itab_w(iw_myrank)%irank   = itabg_w(iw)%irank

! Check if this W point is primary on a remote rank

      if (itab_w(iw_myrank)%irank /= myrank) then

! Turn off some loop flags for these points (some will be turned back on later)
! Loop flags 18 and 21 will be over all IW points (primary and border) on each node

         call wloops('n',iw_myrank,-12,-13,-15,-16,-17,-19,-20,-21,-26,-27)
         call wloops('n',iw_myrank,-28,-29,-30,-34,  0,  0,  0,  0,  0,  0)

! Turn off LBC copy (n/a for global domain) if IWP point is on remote node

         if (.not. myrankflag_w(iwp))  &
            call wloops('n',iw_myrank,-22,-24,-31,-32,-35,0,0,0,0,0)

      endif

! Reset IW neighbor indices to 1

      itab_w(iw_myrank)%iwp = 1

      itab_w(iw_myrank)%im(1:3) = 1
      itab_w(iw_myrank)%iu(1:9) = 1
      itab_w(iw_myrank)%iw(1:9) = 1

! Set indices of neighbors of IW that are present on this rank

      if (myrankflag_w(iwp)) itab_w(iw_myrank)%iwp = itabg_w(iwp)%iw_myrank

      do j = 1,3
         imn = itab_w_pd(iw)%im(j)
         if (myrankflag_m(imn)) itab_w(iw_myrank)%im(j) = itabg_m(imn)%im_myrank
      enddo
      
      do j = 1,9
         iun = itab_w_pd(iw)%iu(j)
         if (myrankflag_u(iun)) itab_w(iw_myrank)%iu(j) = itabg_u(iun)%iu_myrank
      enddo
      
      do j = 1,9
         iwn = itab_w_pd(iw)%iw(j)
         if (myrankflag_w(iwn)) itab_w(iw_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
      enddo

   endif
enddo

! Turn back on some loop flags needed with respect to the local IW point

do iw = 1,nwa
   if (itabg_w(iw)%irank == myrank) then

! Set uloop flag 22 for the 3 nearest U neighbors if IW is primary
! on this rank.

      do j=1,3
         iun = itab_w_pd(iw)%iu(j)
         call uloops('n',itabg_u(iun)%iu_myrank,22,0,0,0,0,0,0,0,0,0)
      enddo

! Set uloop flag 21 for the 9 nearest U neighbors if IW is primary
! on this rank (all the U points of the 3 neighboring triangles)

      do j=1,9
         iun = itab_w_pd(iw)%iu(j)
         call uloops('n',itabg_u(iun)%iu_myrank,21,0,0,0,0,0,0,0,0,0)
      enddo

! Set wloop flag 28 for the 3 nearest W neighbors if IW is primary
! on this rank

      do j=1,3
         iwn = itab_w_pd(iw)%iw(j)
         call wloops('n',itabg_w(iwn)%iw_myrank,28,0,0,0,0,0,0,0,0,0)
      enddo
     
   endif
enddo

! Loop over all U points and for each that is primary on a remote rank, 
! access all U and W points in its stencil.

do iu = 2,nua
   if (itabg_u(iu)%irank /= myrank) then
   
! IU point is primary on remote rank.  

! If IU point is in memory of myrank, its value must be received from 
! remote rank.  Add that remote rank to receive table.

      if (myrankflag_u(iu)) call recv_table_u(itabg_u(iu)%irank)

! Add to send table any U or W point that is primary on myrank and is 
! in the stencil of IU.  

      iup  = itab_u_pd(iu)%iup 
      if (itabg_u(iup)%irank == myrank) call send_table_u(iup,itabg_u(iu)%irank)

      do j = 1,12
         iun = itab_u_pd(iu)%iu(j) 
         if (itabg_u(iun)%irank == myrank) call send_table_u(iun,itabg_u(iu)%irank)
      enddo         

      do j = 1,6
         iwn = itab_u_pd(iu)%iw(j)
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_u(iu)%irank)
      enddo         

   endif
enddo

! Loop over all W points and for each that is primary on a remote rank,
! access all U and W points in its stencil.

do iw = 2,nwa
   if (itabg_w(iw)%irank /= myrank) then

! IW point is primary on remote rank.  

! If IW point is in memory of myrank, it must be received from remote rank.
! Add that remote rank to receive table.

      if (myrankflag_w(iw)) call recv_table_w(itabg_w(iw)%irank)

! Add to send table any U or W point that is primary on myrank and is 
! in the stencil of IW.  

      iwp = itab_w_pd(iw)%iwp
      if (itabg_w(iwp)%irank == myrank) call send_table_w(iwp,itabg_w(iw)%irank)

      do j = 1,9
         iun = itab_w_pd(iw)%iu(j) 
         if (itabg_u(iun)%irank == myrank) call send_table_u(iun,itabg_w(iw)%irank)
      enddo

      do j = 1,9
         iwn = itab_w_pd(iw)%iw(j) 
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_w(iw)%irank)
      enddo

! Special for Zalesak monotonic advection:

      do j = 1,3
         iwn = itab_w_pd(iw)%iw(j) 
         do k=1,9
            iun = itab_w_pd(iwn)%iu(k)
            if (itabg_u(iun)%irank == myrank) then
               call send_table_u(iun,itabg_w(iw)%irank)
            endif
         enddo
      enddo

   endif
enddo

!!!! END OF ITAB_? POPULATION




!!!!!!!!!! ISSO DEVE SER DELETADO
if (isfcl == 1) then
   call move_alloc(landflux, landflux_temp)
   call move_alloc( seaflux,  seaflux_temp)
endif
!!!!!!!!!! ISSO DEVE SER DELETADO

! Check whether LAND/SEA models are used

if (isfcl == 1) then

! Copy SEAFLUX values

   mseaflux = 1

   do isf = 2,nseaflux
      iw  = seaflux_temp(isf)%iw
      iws = seaflux_temp(isf)%iwls

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_ws(iws)%irank == myrank) then

         mseaflux = mseaflux + 1
         seaflag(iws) = .true.
      endif
   enddo

   allocate (seaflux(mseaflux))

   mseaflux = 1

   seaflux(1)%ifglobe = 1
   seaflux(1)%iw = 1
   seaflux(1)%iwls = 1

   do isf = 2,nseaflux
      iw  = seaflux_temp(isf)%iw   ! full-domain index
      iws = seaflux_temp(isf)%iwls ! full-domain index

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_ws(iws)%irank == myrank) then

         mseaflux = mseaflux + 1

         seaflux(mseaflux) = seaflux_temp(isf) ! retain global indices

         seafluxg(isf)%isf_myrank = mseaflux
      endif
   enddo

   ! Set the rank of the seaflux cell to the rank of the atm cell
   ! (only needed for the parcombine step)

   do isf = 1,mseaflux
      iw = seaflux(isf)%iw
      seaflux(isf)%iwrank = itabg_w(iw)%irank
   enddo
   
! Copy LANDLUX values

   mlandflux = 1

   do ilf = 2,nlandflux
      iw  = landflux_temp(ilf)%iw
      iwl = landflux_temp(ilf)%iwls

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_wl(iwl)%irank == myrank) then

         mlandflux = mlandflux + 1
         landflag(iwl) = .true.
      endif
   enddo

   allocate (landflux(mlandflux))

   mlandflux = 1

   landflux(1)%ifglobe = 1
   landflux(1)%iw = 1
   landflux(1)%iwls = 1

   do ilf = 2,nlandflux
      iw  = landflux_temp(ilf)%iw   ! full-domain index
      iwl = landflux_temp(ilf)%iwls ! full-domain index

      if (itabg_w (iw )%irank == myrank .or.  &
          itabg_wl(iwl)%irank == myrank) then

         mlandflux = mlandflux + 1

         landflux(mlandflux) = landflux_temp(ilf) ! retain global indices

         landfluxg(ilf)%ilf_myrank = mlandflux
      endif
   enddo

   ! Set the rank of the landlflux cell to the rank of the atm cell
   ! (only needed for the parcombine step)

   do ilf = 1,mlandflux
      iw = landflux(ilf)%iw
      landflux(ilf)%iwrank = itabg_w(iw)%irank
   enddo

   call para_init_sea ( seaflag,  seaflux_temp)
   call para_init_land(landflag, landflux_temp)

endif

! Deallocate temporary data structures and arrays

deallocate (landflux_temp, seaflux_temp)

! Deallocate para_decomp _pd arrays

deallocate (itab_m_pd, itab_u_pd, itab_w_pd)

return
end subroutine para_init

!===============================================================================

subroutine recv_table_u(iremote)

use mem_para,  only: nrecvs_u, recv_u
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

do jrecv=1,nrecvs_u(1)
   if (recv_u(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_u(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_u(1).

if (jrecv > nrecvs_u(1)) nrecvs_u(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_u(jrecv)%iremote = iremote

return
end subroutine recv_table_u

!===============================================================================

subroutine recv_table_w(iremote)

use mem_para,  only: nrecvs_w, recv_w
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote_w is already in table of ranks to receive from

do jrecv = 1,nrecvs_w(1)
   if (recv_w(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_w(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_w(1).

if (jrecv > nrecvs_w(1)) nrecvs_w(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_w(jrecv)%iremote = iremote

return
end subroutine recv_table_w

!===============================================================================

subroutine send_table_u(iu,iremote)

use mem_ijtabs, only: itab_u, itabg_u, mloops_u
use mem_para,   only: nsends_u, send_u, mgroupsize
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iu
integer, intent(in) :: iremote

integer :: jsend
integer :: iu_myrank

! Check whether iremote_u is already in table of ranks to send to

do jsend=1,nsends_u(1)
   if (send_u(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_u, jsend represents a rank not yet entered in the
! table, so increase nsends_u.

if (jsend > nsends_u(1)) nsends_u(1) = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iu_myrank = itabg_u(iu)%iu_myrank

itab_u(iu_myrank)%loop(mloops_u+jsend) = .true.
send_u(jsend)%iremote = iremote

return
end subroutine send_table_u

!===============================================================================

subroutine send_table_w(iw,iremote)

use mem_ijtabs, only: itab_w, itabg_w, mloops_w
use mem_para,   only: nsends_w, send_w
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iw
integer, intent(in) :: iremote

integer :: jsend
integer :: iw_myrank

! Check whether iremote_w is already in table of ranks to send to

do jsend=1,nsends_w(1)
   if (send_w(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_w, jsend represents a rank not yet entered in the
! table, so increase nsends_w.

if (jsend > nsends_w(1)) nsends_w = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iw_myrank = itabg_w(iw)%iw_myrank

itab_w(iw_myrank)%loop(mloops_w+jsend) = .true.
send_w(jsend)%iremote = iremote

return
end subroutine send_table_w

