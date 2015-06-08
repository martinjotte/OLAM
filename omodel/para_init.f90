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
subroutine para_init()

use misc_coms,  only: io6, mdomain, isubdomain

use mem_ijtabs, only: itab_m,      itab_v,      itab_w,      &
                      itab_m_vars, itab_v_vars, itab_w_vars, &
                      itabg_m,     itabg_v,     itabg_w,     &
                      itab_m_pd,   itab_v_pd,   itab_w_pd,   &
                      alloc_itabs, mrls, &
                      jtm_vadj, jtm_lbcp, &
                      jtv_init, jtv_prog, jtv_wadj, jtv_wstn, jtv_lbcp, &
                      jtw_init, jtw_prog, jtw_wadj, jtw_wstn, jtw_lbcp

use mem_grid,   only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, &
                      alloc_gridz, alloc_xyzem, alloc_xyzew, &
                      alloc_grid1, alloc_grid2

use mem_para,   only: mgroupsize, myrank,                        &
                      nsends_v, nsends_w, nsends_m, nsends_wnud, &
                      nrecvs_v, nrecvs_w, nrecvs_m, nrecvs_wnud

use sea_coms,   only: nws, mws

use leaf_coms,  only: nwl, mwl, isfcl

use mem_sea,    only: itab_ws

use mem_leaf,   only: itab_wl

use mem_nudge,  only: nudflag, nudnxp, nwnud, mwnud, itab_wnud, itabg_wnud, &
                      alloc_nudge1

implicit none

integer :: j,imn,ivn,iwn,jnud
integer :: im,iv,iw,iw1,iw2,iwnud,iwnud1,iwnud2,iwnud3
integer :: imp,ivp,iwp
integer :: iws,iwl,ipass,nland,nsea
integer :: npoly
integer :: wadj_flag

integer :: im_myrank ! Counter for M points to be included on this rank
integer :: iv_myrank ! Counter for V points to be included on this rank
integer :: iw_myrank ! Counter for W points to be included on this rank
integer :: iwnud_myrank ! Counter for WNUD points to be included on this rank

! Automatic arrays

logical :: myrankflag_m(nma) ! Flag for M points existing on this rank
logical :: myrankflag_v(nva) ! Flag for V points existing on this rank
logical :: myrankflag_w(nwa) ! Flag for W points existing on this rank
logical :: myrankflag_wnud(nwnud)  ! Flag that ITABW(IW)%IWNUD(1:3) exist on
                                   ! this rank (for IW primary on this rank)
logical :: myrankflag_wnud1(nwnud) ! Flag that ITABW(IW)%IWNUD(1) exists on
                                   ! this rank (for IW primary on this rank)

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_v(mrls)) ; nsends_v(1:mrls) = 0
allocate (nsends_w(mrls)) ; nsends_w(1:mrls) = 0
allocate (nsends_m(mrls)) ; nsends_m(1:mrls) = 0

nsends_wnud = 0

allocate (nrecvs_v(mrls)) ; nrecvs_v(1:mrls) = 0
allocate (nrecvs_w(mrls)) ; nrecvs_w(1:mrls) = 0
allocate (nrecvs_m(mrls)) ; nrecvs_m(1:mrls) = 0

nrecvs_wnud = 0

! Initialize myrank flag arrays to .false.

myrankflag_m(:) = .false.
myrankflag_v(:) = .false.
myrankflag_w(:) = .false.
myrankflag_wnud(:) = .false.
myrankflag_wnud1(:) = .false.

! Loop over all V points, and for each whose assigned irank is equal to myrank,
! flag all M, V, and W points in its computational stencil for inclusion on this
! rank, excluding IVP.

do iv = 2,nva
   if (itabg_v(iv)%irank == myrank) then

      myrankflag_v(iv) = .true.

      myrankflag_m( itab_v_pd(iv)%im(1:6) ) = .true.
      myrankflag_v( itab_v_pd(iv)%iv(1:4) ) = .true.
      myrankflag_w( itab_v_pd(iv)%iw(1:4) ) = .true.

   endif
enddo

! Loop over all W points, and for each whose assigned irank is equal to myrank,
! flag all M, V, W and WNUD points in its computational stencil for
! inclusion on this rank, excluding IWP.

do iw = 2,nwa
   if (itabg_w(iw)%irank == myrank) then

      myrankflag_w(iw) = .true.

! The standard computational stencil of mem_ijtabs

      myrankflag_m( itab_w_pd(iw)%im(1:7) ) = .true.
      myrankflag_v( itab_w_pd(iw)%iv(1:7) ) = .true.
      myrankflag_w( itab_w_pd(iw)%iw(1:7) ) = .true.

      if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
         myrankflag_wnud ( itab_w_pd(iw)%iwnud(1:3) ) = .true.
         myrankflag_wnud1( itab_w_pd(iw)%iwnud(1)   ) = .true.
      endif

   endif
enddo

! Ignore the "dummy" 1st point in case it was activated

myrankflag_m(1) = .false.
myrankflag_v(1) = .false.
myrankflag_w(1) = .false.
myrankflag_wnud(1) = .false.
myrankflag_wnud1(1) = .false.

! Loop over all M, V, W, and WNUD points and count the ones that
! have been flagged for inclusion on this rank.

im_myrank = 1
iv_myrank = 1
iw_myrank = 1
iwnud_myrank = 1

do im = 2,nma
   if (myrankflag_m(im)) then
      im_myrank = im_myrank + 1
   endif
enddo

do iv = 2,nva
   if (myrankflag_v(iv)) then
      iv_myrank = iv_myrank + 1
   endif
enddo

do iw = 2,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1
   endif
enddo

if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
   do iwnud = 2,nwnud
      if (myrankflag_wnud(iwnud)) then
         iwnud_myrank = iwnud_myrank + 1
      endif
   enddo

   mwnud = iwnud_myrank
endif

! Set mma, mva, mwa values for this rank

mma = im_myrank
mva = iv_myrank
mua = mva
mwa = iw_myrank

! Allocate grid structure variables

call alloc_gridz()
call alloc_itabs(mma, mva, mwa, 1)
call alloc_xyzem(mma)
call alloc_xyzew(mwa)
call alloc_grid1(mma, mva, mwa)
call alloc_grid2(mma, mva, mwa)

if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) call alloc_nudge1(mwnud)

! Reset point counts to 1

im_myrank = 1
iv_myrank = 1
iw_myrank = 1
iwnud_myrank = 1

! Store new myrank M, V, W, and WNUD indices in itabg data structures

do im = 1,nma
   if (myrankflag_m(im)) then
      im_myrank = im_myrank + 1
      
      itabg_m(im)%im_myrank = im_myrank      
   endif
enddo

do iv = 1,nva
   if (myrankflag_v(iv)) then
      iv_myrank = iv_myrank + 1
      
      itabg_v(iv)%iv_myrank = iv_myrank
   endif
enddo

do iw = 1,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1

      itabg_w(iw)%iw_myrank = iw_myrank      
   endif
enddo

if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
   do iwnud = 1,nwnud
      if (myrankflag_wnud(iwnud)) then
         iwnud_myrank = iwnud_myrank + 1

         itabg_wnud(iwnud)%iwnud_myrank = iwnud_myrank      
      endif
   enddo
endif

! Defining global index of local points (will be used on gridfile_read)

do im = 1, nma
   if (myrankflag_m(im)) itab_m(itabg_m(im)%im_myrank)%imglobe = im
enddo

do iv = 1, nva
   if (myrankflag_v(iv)) itab_v(itabg_v(iv)%iv_myrank)%ivglobe = iv
enddo

do iw = 1, nwa
   if (myrankflag_w(iw)) itab_w(itabg_w(iw)%iw_myrank)%iwglobe = iw
enddo

if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
   do iwnud = 1, nwnud
      if (myrankflag_wnud(iwnud)) &
         itab_wnud(itabg_wnud(iwnud)%iwnud_myrank)%iwnudglobe = iwnud
   enddo
endif

! Read the grid structure for all points in local parallel subdomain
! for this rank (or for all points in domain if run is sequential)

call gridfile_read()

! itab_m, itab_v, and itab_w data structures that exist in local subdomain
! memory were filled in subroutine gridfile_read, but with member values
! for the full domain.  Those member values that depend on local subdomain
! are reset next.

! Loop over all M points in global domain

do im = 1,nma
   if (myrankflag_m(im)) then ! M point is in memory on local subdomain

      im_myrank = itabg_m(im)%im_myrank ! Local index of M point

! Set indices of neighbors of IM that are present on this rank

      npoly = itab_m_pd(im)%npoly

      do j = 1,npoly
         ivn = itab_m_pd(im)%iv(j) ! Global index
         iwn = itab_m_pd(im)%iw(j) ! Global index

         ! Get the other im neighbor in this direction

         if ( itab_v_pd(ivn)%im(1) == im ) then
            imn = itab_v_pd(ivn)%im(2)
         else
            imn = itab_v_pd(ivn)%im(1)
         endif

! If IVN point is in memory on local parallel subdomain, assign its local
! index to IM stencil; otherwise set index to 1

         if (myrankflag_v(ivn)) then
            itab_m(im_myrank)%iv(j) = itabg_v(ivn)%iv_myrank
         else
            itab_m(im_myrank)%iv(j) = 1
         endif

! If IWN point is in memory on local parallel subdomain, assign its local
! index to IM stencil; otherwise set index to 1

         if (myrankflag_w(iwn)) then
            itab_m(im_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
         else
            itab_m(im_myrank)%iw(j) = 1
         endif
      enddo

! Turn off jtm_vadj loop flag if M point is not primary on this subdomain

      if (itabg_m(im)%irank /= myrank) then
         call mloopf('n', im_myrank, -jtm_vadj, 0, 0, 0, 0, 0)
      endif

! If IMP point exists on local parallel subdomain, assign its local
! index to IM stencil; otherwise set index to 1 and turn off lbcp flag

      imp = itab_m_pd(im)%imp ! Global index

      if (myrankflag_m(imp)) then
         itab_m(im_myrank)%imp = itabg_m(imp)%im_myrank ! Local index
      else
         itab_m(im_myrank)%imp = 1
         call mloopf('n', im_myrank, -jtm_lbcp, 0, 0, 0, 0, 0)
      endif

   endif
enddo

! Loop over all V points in global domain

do iv = 1,nva
   if (myrankflag_v(iv)) then ! V point is in memory on local subdomain

      iv_myrank = itabg_v(iv)%iv_myrank ! Local index of V point

      itab_v(iv_myrank)%irank = itabg_v(iv)%irank

! Set local indices of M points in V stencil

      do j = 1,6
         imn = itab_v_pd(iv)%im(j) ! Global index

         if (myrankflag_m(imn)) then
            itab_v(iv_myrank)%im(j) = itabg_m(imn)%im_myrank
         else
            itab_v(iv_myrank)%im(j) = 1
         endif
      enddo

! Set local indices of V points in V stencil

      do j = 1,4
         ivn = itab_v_pd(iv)%iv(j) ! Global index

         if (myrankflag_v(ivn)) then
            itab_v(iv_myrank)%iv(j) = itabg_v(ivn)%iv_myrank
         else
            itab_v(iv_myrank)%iv(j) = 1
         endif
      enddo

! Set local indices of W points in V stencil

      do j = 1,4
         iwn = itab_v_pd(iv)%iw(j) ! Global index

         if (myrankflag_w(iwn)) then
            itab_v(iv_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
         else
            itab_v(iv_myrank)%iw(j)  = 1
         endif
      enddo

! If IV point is not primary on local parallel subdomain,
! turn off jtv_init and jtv_prog flags

      if (itabg_v(iv)%irank /= myrank) then
         call vloopf('n',iv_myrank, -jtv_init, -jtv_prog, 0, 0, 0, 0)
      endif

! If adjacent W points are both not primary on local parallel subdomain,
! turn off jtv_wadj and jtv_wstn flags

      iw1 = itab_v_pd(iv)%iw(1) ! Global index
      iw2 = itab_v_pd(iv)%iw(2) ! Global index

      if (itabg_w(iw1)%irank /= myrank .and. itabg_w(iw2)%irank /= myrank) then
         call vloopf('n',iv_myrank, -jtv_wadj, -jtv_wstn, 0, 0, 0, 0)
      endif
      
! If IVP point exists on local parallel subdomain, assign its local
! index to IV stencil; otherwise set index to 1 and turn off lbcp flag

      ivp = itab_v_pd(iv)%ivp ! Global index

      if (myrankflag_v(ivp)) then
         itab_v(iv_myrank)%ivp = itabg_v(ivp)%iv_myrank
      else
         itab_v(iv_myrank)%ivp = 1
         call vloopf('n',iv_myrank, -jtv_lbcp, 0, 0, 0, 0, 0)
      endif

   endif
enddo

! Loop over all W points in global domain

do iw = 1,nwa
   if (myrankflag_w(iw)) then ! W point is in memory on local subdomain

      iw_myrank = itabg_w(iw)%iw_myrank ! Local index of W point

      itab_w(iw_myrank)%irank = itabg_w(iw)%irank

! Set local indices of M, V, W points in W stencil

      wadj_flag = 0
      if (itabg_w(iw)%irank == myrank) wadj_flag = 1

      do j = 1,7
         imn = itab_w_pd(iw)%im(j) ! Global index
         ivn = itab_w_pd(iw)%iv(j) ! Global index
         iwn = itab_w_pd(iw)%iw(j) ! Global index

         if (myrankflag_m(imn)) then
            itab_w(iw_myrank)%im(j) = itabg_m(imn)%im_myrank
         else
            itab_w(iw_myrank)%im(j) = 1
         endif

         if (myrankflag_v(ivn)) then
            itab_w(iw_myrank)%iv(j) = itabg_v(ivn)%iv_myrank
         else
            itab_w(iw_myrank)%iv(j) = 1
         endif

         if (myrankflag_w(iwn)) then
            itab_w(iw_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
         else
            itab_w(iw_myrank)%iw(j) = 1
         endif

         if (itabg_w(iwn)%irank == myrank) wadj_flag = 1
      enddo

! Set local indices of WNUD points in W stencil

      if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
         do jnud = 1,3
            iwnud = itab_w_pd(iw)%iwnud(jnud) ! Global index
            itab_w(iw_myrank)%iwnud(jnud) = itabg_wnud(iwnud)%iwnud_myrank
         enddo
      endif

! If IW point is not primary on local parallel subdomain,
! turn off jtw_init and jtw_prog flags

      if (itabg_w(iw)%irank /= myrank) then
         call wloopf('n',iw_myrank, -jtw_init, -jtw_prog, 0, 0, 0, 0)
      endif

! If none of the adjacent W points is primary on local parallel subdomain,
! turn off jtw_wadj and jtw_wstn flags

      if (wadj_flag == 0) then
         call wloopf('n',iw_myrank, -jtw_wadj, -jtw_wstn, 0, 0, 0, 0)
      endif

! If IWP point exists on local parallel subdomain, assign its local
! index to IW stencil; otherwise set index to 1 and turn off lbcp flag

      iwp = itab_w_pd(iw)%iwp ! Global index

      if (myrankflag_w(iwp)) then
         itab_w(iw_myrank)%iwp = itabg_w(iwp)%iw_myrank
      else
         itab_w(iw_myrank)%iwp = 1
         call wloopf('n',iw_myrank, -jtw_lbcp, 0, 0, 0, 0, 0)
      endif

   endif
enddo

! Loop over all WNUD points in global domain

do iwnud = 1,nwnud
   if (myrankflag_wnud(iwnud)) then ! WNUD point is in memory on local subdomain

      iwnud_myrank = itabg_wnud(iwnud)%iwnud_myrank ! Local index of WNUD point

      itab_wnud(iwnud_myrank)%irank = itabg_wnud(iwnud)%irank
   endif
enddo

! Loop over all V points and for each that is primary on a remote rank, 
! access all V, W and M points in its stencil.

do iv = 2,nva
   if (itabg_v(iv)%irank /= myrank) then

! IV point is primary on remote rank.  

! If IV point is in memory of myrank, its value must be received from 
! remote rank.  Add that remote rank to receive table.

      ivp  = itab_v_pd(iv)%ivp

      if (myrankflag_v(iv)) then
         if (itabg_v(ivp)%irank /= myrank) call recv_table_v(itabg_v(ivp)%irank)
      endif

! Add to send table any V or W point that is primary on myrank and is 
! in the stencil of IV.

      if (itabg_v(ivp)%irank == myrank) call send_table_v(ivp,itabg_v(iv)%irank)

      do j = 1,4
         ivn = itab_v_pd(iv)%iv(j)
         ivn = itab_v_pd(ivn)%ivp
         if (itabg_v(ivn)%irank == myrank) call send_table_v(ivn,itabg_v(iv)%irank)
      enddo

      do j = 1,4
         iwn = itab_v_pd(iv)%iw(j)
         iwn = itab_w_pd(iwn)%iwp
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_v(iv)%irank)
      enddo

      do j = 1,6
         imn = itab_v_pd(iv)%im(j)
         imn = itab_m_pd(imn)%imp
         if (itabg_m(imn)%irank == myrank) call send_table_m(imn,itabg_v(iv)%irank)
      enddo

   endif

! If IV point in in memory of myrank but its IVP is not, add the remote rank
! of the IVP point to the receive table

   ivp = itab_v_pd(iv)%ivp
   if (itabg_v(ivp)%irank /= myrank) then
      if (myrankflag_v(iv)) call recv_table_v(itabg_v(ivp)%irank)
   endif

enddo

! Loop over all W points and for each that is primary on a remote rank,
! access all V and W points in its stencil.

do iw = 2,nwa

   if (itabg_w(iw)%irank /= myrank) then

! IW point is primary on remote rank.  

! If IW point is in memory of myrank, it must be received from remote rank.
! Add that remote rank to receive table.

      iwp = itab_w_pd(iw)%iwp

      if (myrankflag_w(iw)) then
         if (itabg_w(iwp)%irank /= myrank) call recv_table_w(itabg_w(iwp)%irank)
      endif

! Add to send table any V or W point that is primary on myrank and is 
! in the stencil of IW.  

      if (itabg_w(iwp)%irank == myrank) call send_table_w(iwp,itabg_w(iw)%irank)

      do j = 1,7
         ivn = itab_w_pd(iw)%iv(j)
         ivn = itab_v_pd(ivn)%ivp
         if (itabg_v(ivn)%irank == myrank) call send_table_v(ivn,itabg_w(iw)%irank)
      enddo

      do j = 1,7
         iwn = itab_w_pd(iw)%iw(j)
         iwn = itab_w_pd(iwn)%iwp
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_w(iw)%irank)
      enddo

      do j = 1,7
         imn = itab_w_pd(iw)%im(j)
         imn = itab_m_pd(imn)%imp
         if (itabg_m(imn)%irank == myrank) call send_table_m(imn,itabg_w(iw)%irank)
      enddo

! Although a primary rank is assigned to each nudging point, the assignment is
! used only to designate which rank writes values for that nudging point to
! history files.  Parallel MPI communication, on the other hand, must take into
! account that computation of nudging point values is based on sums over model
! grid points on multiple ranks.  Each of the ranks that contribute to a given
! nudging point accumulate a partial sum and send the partial sum to all other
! ranks that need that nudging point value.  After receiving the partial sums,
! each rank completes the sum (a small duplication of effort) over the partial
! sums to obtain its own copy of the full sum.  Thus, for building send and
! recv tables, we check itab_w_pd(iw)%iwnud(1:3) and itab_w_pd(iw)%iwnud(1)
! instead of itab_wnud()%irank.

      if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
         iwnud1 = itab_w_pd(iw)%iwnud(1)
         iwnud2 = itab_w_pd(iw)%iwnud(2)
         iwnud3 = itab_w_pd(iw)%iwnud(3)

         if (myrankflag_wnud(iwnud1)) call recv_table_wnud(itabg_w(iw)%irank)

         if (myrankflag_wnud1(iwnud1)) call send_table_wnud(iwnud1,itabg_w(iw)%irank)
         if (myrankflag_wnud1(iwnud2)) call send_table_wnud(iwnud2,itabg_w(iw)%irank)
         if (myrankflag_wnud1(iwnud3)) call send_table_wnud(iwnud3,itabg_w(iw)%irank)
      endif

   endif

! If IW point in in memory of myrank but its IWP is not, add the remote rank
! of the IWP point to the receive table

   iwp = itab_w_pd(iw)%iwp
   if (itabg_w(iwp)%irank /= myrank) then
      if (myrankflag_w(iw)) call recv_table_w(itabg_w(iwp)%irank)
   endif

enddo

! Loop over all M points and for each that is primary on a remote rank,
! access all points in its stencil.

do im = 2, nma
   if (itabg_m(im)%irank /= myrank) then

! IM point is primary on remote rank.  

! If IM point is in memory of myrank, it must be received from remote rank.
! Add that remote rank to receive table.

      imp = itab_m_pd(im)%imp

      if (myrankflag_m(im)) then
         if (itabg_m(imp)%irank /= myrank) call recv_table_m(itabg_m(imp)%irank)
      endif

! Add to send table any V or M point that is primary on myrank and is 
! in the stencil of IM.  

      if (itabg_m(imp)%irank == myrank) call send_table_m(imp,itabg_m(im)%irank)

   endif

! If IM point in in memory of myrank but its IMP is not, add the remote rank
! of the IMP point to the receive table

   imp = itab_m_pd(im)%imp
   if (itabg_m(imp)%irank /= myrank) then
      if (myrankflag_m(im)) call recv_table_m(itabg_m(imp)%irank)
   endif

enddo

do iw = 2, nwa
   if (myrankflag_w(iw)) then
      iwp = itab_w_pd(iw)%iwp
      if (.not. myrankflag_w(iwp)) then
         if (itabg_w(iwp)%iw_myrank_iwp == -1) then
            ! never been set, we will MPI copy to this point
            itabg_w(iwp)%iw_myrank_iwp = itabg_w(iw)%iw_myrank
         else
            ! We are already MPI copying this point, so do an LBC copy instead
            iw_myrank  = itabg_w(iw)%iw_myrank
            itab_w(iw_myrank)%iwp = itabg_w(iwp)%iw_myrank_iwp
            call wloopf('n', iw_myrank, jtw_lbcp, 0, 0, 0, 0, 0)
         endif
      endif
   endif
enddo

do iv = 2, nva
   if (myrankflag_v(iv)) then
      ivp = itab_v_pd(iv)%ivp
      if (.not. myrankflag_v(ivp)) then
         if (itabg_v(ivp)%iv_myrank_ivp == -1) then
            ! never been set, we will MPI copy to this point
            itabg_v(ivp)%iv_myrank_ivp = itabg_v(iv)%iv_myrank
         else
            ! We are already MPI copying this point, so do an LBC copy instead
            iv_myrank  = itabg_v(iv)%iv_myrank
            itab_v(iv_myrank)%ivp = itabg_v(ivp)%iv_myrank_ivp
            call vloopf('n', iv_myrank, jtv_lbcp, 0, 0, 0, 0, 0)
         endif
      endif
   endif
enddo

do im = 2, nma
   if (myrankflag_m(im)) then
      imp = itab_m_pd(im)%imp
      if (.not. myrankflag_m(imp)) then
         if (itabg_m(imp)%im_myrank_imp == -1) then
            ! never been set, we will MPI copy to this point
            itabg_m(imp)%im_myrank_imp = itabg_m(im)%im_myrank
         else
            ! We are already MPI copying this point, so do an LBC copy instead
            im_myrank  = itabg_m(im)%im_myrank
            itab_m(im_myrank)%imp = itabg_m(imp)%im_myrank_imp
            call mloopf('n', im_myrank, jtm_lbcp, 0, 0, 0, 0, 0)
         endif
      endif
   endif
enddo

! Check whether LAND/SEA models are used

if (isfcl == 1) then

   call para_init_sea ()
   call para_init_land()

! Do two passes through the following code to build lists of attached land
! and sea cells for each IW column

   do ipass = 1,2

! Set nland and nsea counters to zero for all atmosphere IW columns

      itab_w(1:mwa)%nland = 0
      itab_w(1:mwa)%nsea = 0

! Loop over all land cells and get atmosphere column iw index

      do iwl = 2,mwl
         iw = itab_wl(iwl)%iw
         if (isubdomain == 1) then
            iw = itabg_w(iw)%iw_myrank
         endif

! Increment land cell counter for IW column

         itab_w(iw)%nland = itab_w(iw)%nland + 1

! If second pass, enter IWL land cell index in itab_w(iw)%iland array

         if (ipass == 2) itab_w(iw)%iland(itab_w(iw)%nland) = iwl
      enddo

! Loop over all sea cells and get atmosphere column iw index

      do iws = 2,mws
         iw = itab_ws(iws)%iw
         if (isubdomain == 1) then
            iw = itabg_w(iw)%iw_myrank
         endif

! Increment sea cell counter for IW column

         itab_w(iw)%nsea = itab_w(iw)%nsea + 1

! If second pass, enter IWS sea cell index in itab_w(iw)%isea array

         if (ipass == 2) itab_w(iw)%isea(itab_w(iw)%nsea) = iws
      enddo

! If first pass, allocate iland and isea members of itab_w(iw)

      if (ipass == 1) then
         do iw = 2,mwa
            nland = itab_w(iw)%nland
            nsea  = itab_w(iw)%nsea
            allocate(itab_w(iw)%iland(max(1,nland)))
            allocate(itab_w(iw)%isea(max(1,nsea)))
         enddo
      endif

   enddo  ! ipass

endif  ! isfcl = 1

call compute_primary_points()

! Deallocate para_decomp _pd arrays

deallocate (itab_m_pd, itab_v_pd, itab_w_pd)

end subroutine para_init

!===============================================================================

subroutine recv_table_v(iremote)

use mem_para,  only: nrecvs_v, recv_v
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

do jrecv=1,nrecvs_v(1)
   if (recv_v(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_v(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_v(1).

if (jrecv > nrecvs_v(1)) nrecvs_v(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_v(jrecv)%iremote = iremote

return
end subroutine recv_table_v

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

subroutine recv_table_m(iremote)

use mem_para,  only: nrecvs_m, recv_m
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote is already in table of ranks to receive from

do jrecv = 1, nrecvs_m(1)
   if (recv_m(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_v(1), jrecv represents a rank not yet entered in the
! table, so increase nrecvs_v(1).

if (jrecv > nrecvs_m(1)) nrecvs_m(1) = jrecv

! Enter remote rank in recv-remote-rank table.

recv_m(jrecv)%iremote = iremote

end subroutine recv_table_m

!===============================================================================

subroutine recv_table_wnud(iremote)

use mem_para,  only: nrecvs_wnud, recv_wnud
use misc_coms, only: io6

implicit none

integer, intent(in) :: iremote

integer :: jrecv

! Check whether iremote_w is already in table of ranks to receive from

do jrecv = 1,nrecvs_wnud
   if (recv_wnud(jrecv)%iremote == iremote) exit
enddo

! If jrecv exceeds nrecvs_wnud, jrecv represents a rank not yet entered in the
! table, so increase nrecvs_wnud.

if (jrecv > nrecvs_wnud) nrecvs_wnud = jrecv

! Enter remote rank in recv-remote-rank table.

recv_wnud(jrecv)%iremote = iremote

return
end subroutine recv_table_wnud

!===============================================================================

subroutine send_table_v(iv,iremote)

use mem_ijtabs, only: itab_v, itabg_v, mloops
use mem_para,   only: nsends_v, send_v, mgroupsize
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iv
integer, intent(in) :: iremote

integer :: jsend
integer :: iv_myrank

! Check whether iremote_v is already in table of ranks to send to

do jsend=1,nsends_v(1)
   if (send_v(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_v, jsend represents a rank not yet entered in the
! table, so increase nsends_v.

if (jsend > nsends_v(1)) nsends_v(1) = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iv_myrank = itabg_v(iv)%iv_myrank

itab_v(iv_myrank)%loop(mloops+jsend) = .true.
send_v(jsend)%iremote = iremote

return
end subroutine send_table_v

!===============================================================================

subroutine send_table_w(iw,iremote)

use mem_ijtabs, only: itab_w, itabg_w, mloops
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

itab_w(iw_myrank)%loop(mloops+jsend) = .true.
send_w(jsend)%iremote = iremote

return
end subroutine send_table_w

!===============================================================================

subroutine send_table_m(im,iremote)

use mem_ijtabs, only: itab_m, itabg_m, mloops
use mem_para,   only: nsends_m, send_m
use misc_coms,  only: io6

implicit none

integer, intent(in) :: im
integer, intent(in) :: iremote

integer :: jsend
integer :: im_myrank

! Check whether iremote_m is already in table of ranks to send to

do jsend = 1, nsends_m(1)
   if (send_m(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_m, jsend represents a rank not yet entered in the
! table, so increase nsends_m.

if (jsend > nsends_m(1)) nsends_m(1) = jsend

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

im_myrank = itabg_m(im)%im_myrank

itab_m(im_myrank)%loop(mloops+jsend) = .true.
send_m(jsend)%iremote = iremote

end subroutine send_table_m

!===============================================================================

subroutine send_table_wnud(iwnud,iremote)

use mem_nudge, only: nloops_wnud, itab_wnud, itabg_wnud
use mem_para,  only: nsends_wnud, send_wnud
use misc_coms, only: io6

implicit none

integer, intent(in) :: iwnud
integer, intent(in) :: iremote

integer :: jsend
integer :: iwnud_myrank

! Check whether iremote is already in table of ranks to send to

do jsend=1,nsends_wnud
   if (send_wnud(jsend)%iremote == iremote) exit
enddo

! If jsend exceeds nsends_wnud, jsend represents a rank not yet entered in the
! table, so increase nsends_w.

if (jsend > nsends_wnud) then
    nsends_wnud = jsend
    write(io6,*) 'increasing nsends_wnud to ',nsends_wnud
    if (nsends_wnud > nloops_wnud) stop 'stop: nsends_wnud > nloops_wnud '
endif

! Enter point in send-point table, and enter remote rank in send-remote-rank table.

iwnud_myrank = itabg_wnud(iwnud)%iwnud_myrank

itab_wnud(iwnud_myrank)%loop(jsend) = .true.
send_wnud(jsend)%iremote = iremote

return
end subroutine send_table_wnud

!===============================================================================

subroutine compute_primary_points()

  use mem_grid,   only: mma, mva, mwa
  use mem_ijtabs, only: itab_v, itab_w, itab_m, itabg_m
  use mem_nudge,  only: mwnud, itab_wnud
  use misc_coms,  only: io6
  use leaf_coms,  only: mwl
  use mem_leaf,   only: itab_wl
  use sea_coms,   only: mws
  use mem_sea,    only: itab_ws
  use mem_para,   only: mva_primary, iva_globe_primary, iva_local_primary, &
                        mwa_primary, iwa_globe_primary, iwa_local_primary, &
                        mma_primary, ima_globe_primary, ima_local_primary, &
                        mwl_primary, iwl_globe_primary, iwl_local_primary, &
                        mws_primary, iws_globe_primary, iws_local_primary, &
                        mwnud_primary, iwnud_globe_primary, iwnud_local_primary, &
                        myrank

  implicit none

  integer :: i, ia, istart

  do i = 2, mwa
     if (itab_w(i)%irank == myrank) mwa_primary = mwa_primary + 1
  enddo

  do i = 2, mva
     if (itab_v(i)%irank == myrank) mva_primary = mva_primary + 1
  enddo

  do i = 2, mma
     if (itabg_m( itab_m(i)%imglobe )%irank == myrank) mma_primary = mma_primary + 1
  enddo

  do i = 2, mwl
     mwl_primary = mwl_primary + 1
  enddo

  do i = 2, mws
     mws_primary = mws_primary + 1
  enddo

  do i = 2, mwnud
     if (itab_wnud(i)%irank == myrank) mwnud_primary = mwnud_primary + 1
  enddo

  ! additional space for dummy 1st point included with rank 0 only

  if (myrank == 0) then
     mwa_primary = mwa_primary + 1 
     mma_primary = mma_primary + 1
     mva_primary = mva_primary + 1
     mwl_primary = mwl_primary + 1
     mws_primary = mws_primary + 1
     mwnud_primary = mwnud_primary + 1 
  endif

  ! allocate space for primary global indices

  allocate(iwa_globe_primary(mwa_primary))
  allocate(iwa_local_primary(mwa_primary))

  allocate(iva_globe_primary(mva_primary))
  allocate(iva_local_primary(mva_primary))

  allocate(ima_globe_primary(mma_primary))
  allocate(ima_local_primary(mma_primary))

  allocate(iwl_globe_primary(mwl_primary))
  allocate(iwl_local_primary(mwl_primary))

  allocate(iws_globe_primary(mws_primary))
  allocate(iws_local_primary(mws_primary))

  allocate(iwnud_globe_primary(mwnud_primary))
  allocate(iwnud_local_primary(mwnud_primary))

  ! dummy 1st point included with rank 0 only

  if (myrank == 0) then

     iwa_globe_primary(1) = 1
     iwa_local_primary(1) = 1

     iva_globe_primary(1) = 1
     iva_local_primary(1) = 1

     ima_globe_primary(1) = 1
     ima_local_primary(1) = 1

     iwl_globe_primary(1) = 1
     iwl_local_primary(1) = 1

     iws_globe_primary(1) = 1
     iws_local_primary(1) = 1

     iwnud_globe_primary(1) = 1
     iwnud_local_primary(1) = 1

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

  if (ia /= mwa_primary) stop "error computing number of primary points1"

  ia = istart

  do i = 2, mva
     if (itab_v(i)%irank == myrank) then
        ia = ia + 1
        iva_globe_primary(ia) = itab_v(i)%ivglobe
        iva_local_primary(ia) = i
     endif
  enddo

  if (ia /= mva_primary) stop "error computing number of primary points2"

  ia = istart

  do i = 2, mma
     if (itabg_m( itab_m(i)%imglobe )%irank == myrank) then
        ia = ia + 1
        ima_globe_primary(ia) = itab_m(i)%imglobe
        ima_local_primary(ia) = i
     endif
  enddo

  if (ia /= mma_primary) stop "error computing number of primary points3"

  ia = istart

  do i = 2, mwl
     ia = ia + 1
     iwl_globe_primary(ia) = itab_wl(i)%iwglobe
     iwl_local_primary(ia) = i
  enddo

  if (ia /= mwl_primary) stop "error computing number of primary points4"

  ia = istart

  do i = 2, mws
     ia = ia + 1
     iws_globe_primary(ia) = itab_ws(i)%iwglobe
     iws_local_primary(ia) = i
  enddo

  if (ia /= mws_primary) stop "error computing number of primary points5"

  ia = istart

  do i = 2, mwnud
     if (itab_wnud(i)%irank == myrank) then
        ia = ia + 1
        iwnud_globe_primary(ia) = itab_wnud(i)%iwnudglobe
        iwnud_local_primary(ia) = i
     endif
  enddo

  if (ia /= mwnud_primary) stop "error computing number of primary points6"

!!!! temporary checks !!!!!!!!!!!!!!!
do i = 2, mva_primary
   if (iva_globe_primary(i) < iva_globe_primary(i-1)) then
      write(io6,*) 'error: VA is out of order!!!!'
      stop
   endif
enddo

do i = 2, mwa_primary
   if (iwa_globe_primary(i) < iwa_globe_primary(i-1)) then
      write(io6,*) 'error: WA is out of order!!!!'
      stop
   endif
enddo

do i = 2, mma_primary
   if (ima_globe_primary(i) < ima_globe_primary(i-1)) then
      write(io6,*) 'error: MA is out of order!!!!'
      stop
   endif
enddo

do i = 2, mwl_primary
   if (iwl_globe_primary(i) < iwl_globe_primary(i-1)) then
      write(io6,*) 'error: WL is out of order!!!!'
      stop
   endif
enddo

do i = 2, mws_primary
   if (iws_globe_primary(i) < iws_globe_primary(i-1)) then
      write(io6,*) 'error: WS is out of order!!!!'
      stop
   endif
enddo

do i = 2, mwnud_primary
   if (iwnud_globe_primary(i) < iwnud_globe_primary(i-1)) then
      write(io6,*) 'error: WNUD is out of order!!!!'
      stop
   endif
enddo

!!!!!!! temporary checks 2

do i = 2, mva
   if (itab_v(i)%ivglobe < itab_v(i-1)%ivglobe) then
      write(io6,*) 'error: VAGLOBE is out of order!!!!'
      stop
   endif
enddo

do i = 2, mwa
   if (itab_w(i)%iwglobe < itab_w(i-1)%iwglobe) then
      write(io6,*) 'error: WAGLOBE is out of order!!!!'
      stop
   endif
enddo

do i = 2, mma
   if (itab_m(i)%imglobe < itab_m(i-1)%imglobe) then
      write(io6,*) 'error: MAGLOBE is out of order!!!!'
      stop
   endif
enddo

do i = 2, mwl
   if (itab_wl(i)%iwglobe < itab_wl(i-1)%iwglobe) then
      write(io6,*) 'error: WLGLOBE is out of order!!!!'
      stop
   endif
enddo

do i = 2, mws
   if (itab_ws(i)%iwglobe < itab_ws(i-1)%iwglobe) then
      write(io6,*) 'error: WSGLOBE is out of order!!!!'
      stop
   endif
enddo

do i = 2, mwnud
   if (itab_wnud(i)%iwnudglobe < itab_wnud(i-1)%iwnudglobe) then
      write(io6,*) 'error: WNUDGLOBE is out of order!!!!'
      stop
   endif
enddo

end subroutine compute_primary_points
