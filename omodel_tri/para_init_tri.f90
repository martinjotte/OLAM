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
                      alloc_itabs, mrls

use mem_grid,   only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, &
                      lpm, lpu, lcu, lpw, lsw, &
                      xem, yem, zem, xeu, yeu, zeu, &
                      xev, yev, zev, xew, yew, zew, &
                      unx, uny, unz, vnx, vny, vnz, wnx, wny, wnz, &
                      dnu, dniu, dnv, dniv, arm0, arw0, topm, topw, &
                      glatm, glonm, glatu, glonu, glatv, glonv, glatw, glonw, &
                      aru, arv, arw, volui, volvi, volwi, volt, volti, &
                      alloc_xyzem, alloc_xyzew, alloc_grid1, alloc_grid2

use mem_para,   only: mgroupsize, myrank, &
                      send_u, recv_u, send_v, recv_v, send_w, recv_w, &
                      nsends_u, nsends_v, nsends_w, &
                      nrecvs_u, nrecvs_v, nrecvs_w, &
                      send_uf, recv_uf, send_vf, recv_vf

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

integer :: lpm_temp(nma)
integer :: lpu_temp(nua),lcu_temp(nua)
integer :: lpw_temp(nwa),lsw_temp(nwa)

real :: topm_temp(nma), glatm_temp(nma), glonm_temp(nma), arm0_temp(nma)

real :: xem_temp(nma),yem_temp(nma),zem_temp(nma)
real :: xeu_temp(nua),yeu_temp(nua),zeu_temp(nua)
real :: xew_temp(nwa),yew_temp(nwa),zew_temp(nwa)
real :: unx_temp(nua),uny_temp(nua),unz_temp(nua)
real :: vnx_temp(nua),vny_temp(nua),vnz_temp(nua)
real :: wnx_temp(nwa),wny_temp(nwa),wnz_temp(nwa)

real :: glatw_temp(nwa),glonw_temp(nwa),arw0_temp(nwa),topw_temp(nwa)
real :: dnu_temp(nua),dniu_temp(nua)
real :: dnv_temp(nua),dniv_temp(nua)

real :: aru_temp(nza,nua),arw_temp(nza,nwa)
real :: glatu_temp(nua),glonu_temp(nua)

real :: volui_temp(nza,nua),volwi_temp(nza,nwa)

real(kind=8) :: volt_temp(nza,nwa),volti_temp(nza,nwa)

! Temporary datatypes

type (itab_m_vars), allocatable :: ltab_m(:)
type (itab_u_vars), allocatable :: ltab_u(:)
type (itab_w_vars), allocatable :: ltab_w(:)

type(flux_vars), allocatable :: landflux_temp(:)
type(flux_vars), allocatable ::  seaflux_temp(:)

! Move data to temporary data structures, nullifying the old datatype

call move_alloc(itab_m, ltab_m)
call move_alloc(itab_u, ltab_u)
call move_alloc(itab_w, ltab_w)

if (isfcl == 1) then
   call move_alloc(landflux, landflux_temp)
   call move_alloc( seaflux,  seaflux_temp)
endif

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

! Copy grid coordinates to temporary arrays

do im = 1,nma
   lpm_temp(im)  = lpm(im)

   topm_temp(im) = topm(im)

   xem_temp(im)  = xem(im)
   yem_temp(im)  = yem(im)
   zem_temp(im)  = zem(im)

   glatm_temp(im) = glatm(im)
   glonm_temp(im) = glonm(im)

   arm0_temp(im) = arm0(im)
enddo

do iu = 1,nua
   lpu_temp(iu) = lpu(iu)
   lcu_temp(iu) = lcu(iu)

   xeu_temp(iu) = xeu(iu)
   yeu_temp(iu) = yeu(iu)
   zeu_temp(iu) = zeu(iu)
   
   unx_temp(iu) = unx(iu)
   uny_temp(iu) = uny(iu)
   unz_temp(iu) = unz(iu)

   vnx_temp(iu) = vnx(iu)
   vny_temp(iu) = vny(iu)
   vnz_temp(iu) = vnz(iu)

   dnu_temp(iu)  = dnu(iu)
   dniu_temp(iu) = dniu(iu)
   
   dnv_temp(iu)  = dnv(iu)
   dniv_temp(iu) = dniv(iu)
   
   glatu_temp(iu) = glatu(iu)
   glonu_temp(iu) = glonu(iu)

   do k = 1,nza
      aru_temp(k,iu) = aru(k,iu)
      volui_temp(k,iu) = volui(k,iu)
   enddo
enddo

do iw = 1,nwa
   topw_temp(iw) = topw(iw)

   lpw_temp(iw) = lpw(iw)
   lsw_temp(iw) = lsw(iw)

   xew_temp(iw) = xew(iw)
   yew_temp(iw) = yew(iw)
   zew_temp(iw) = zew(iw)

   wnx_temp(iw) = wnx(iw)
   wny_temp(iw) = wny(iw)
   wnz_temp(iw) = wnz(iw)

   glatw_temp(iw) = glatw(iw)
   glonw_temp(iw) = glonw(iw)

   arw0_temp(iw) = arw0(iw)
   
   do k = 1,nza
      arw_temp(k,iw) = arw(k,iw)
      volwi_temp(k,iw) = volwi(k,iw)
      volt_temp(k,iw) = volt(k,iw)
      volti_temp(k,iw) = volti(k,iw)
   enddo
enddo

! Deallocate grid coordinate arrays

deallocate (lpm, lpu, lcu, lpw, lsw)
deallocate (xem, yem, zem, topm, glatm, glonm, arm0)
deallocate (xeu, yeu, zeu, unx, uny, unz, vnx, vny, vnz, dnu, dniu, dnv, dniv)
deallocate (aru, volui, glatu, glonu)
deallocate (xew, yew, zew, wnx, wny, wnz, glatw, glonw, arw0, topw)
deallocate (arw, volwi, volt, volti)

! Loop over all U points, and for each whose assigned irank is equal to myrank,
! flag all U and W points in its computational stencil for inclusion on this
! rank, excluding IUP and IWP.

do iu = 2,nua

   if (itabg_u(iu)%irank == myrank) then

      myrankflag_u(iu) = .true.

      myrankflag_u( ltab_u(iu)%iu(1:12) ) = .true.
      myrankflag_w( ltab_u(iu)%iw(1:6)  ) = .true.

   endif
enddo

! Loop over all W points, and for each whose assigned irank is equal to myrank,
! flag all U and W points in its computational stencil for inclusion on this
! rank, excluding IUP and IWP.

do iw = 2,nwa

   if (itabg_w(iw)%irank == myrank) then

      myrankflag_w(iw) = .true.

! The standard computational stencil of mem_ijtabs

      myrankflag_w( ltab_w(iw)%iw(1:9) ) = .true.
      myrankflag_u( ltab_w(iw)%iu(1:9) ) = .true.

! Special for Zalesak monotonic advection (each W point needs the 9 
! U points from each of the 3 bordering triangles available)

      do j=1,3
         iwn = ltab_w(iw)%iw(j)
         myrankflag_u( ltab_w(iwn)%iu(1:9) ) = .true.
      enddo

   endif
enddo

! Loop over all U points, and for each that has been flagged for inclusion
! on this rank, flag both its M points for inclusion on this rank.
! Count U points also.

do iu = 2,nua

   if (myrankflag_u(iu)) then

      myrankflag_m( ltab_u(iu)%im(1:2) ) = .true.
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

! Allocate itab data structures and main grid coordinate arrays

call alloc_itabs(meshtype,mma,mua,mva,mwa)
call alloc_xyzem(mma)
call alloc_xyzew(mwa)
call alloc_grid1(meshtype)
call alloc_grid2(meshtype)

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

      iup = ltab_u(iu)%iup
      itabg_u(iup)%iu_myrank = iu_myrank
   endif
enddo

do iw = 1,nwa
   if (myrankflag_w(iw)) then
      iw_myrank = iw_myrank + 1

      itabg_w(iw)%iw_myrank = iw_myrank      

! Fill itabg_w(iw)%iw_myrank value for IWP point

      iwp = ltab_w(iw)%iwp
      itabg_w(iwp)%iw_myrank = iw_myrank
   endif
enddo

! M point memory copy

do im = 1,nma
   if (myrankflag_m(im)) then
      npoly = ltab_m(im)%npoly
   
      im_myrank = itabg_m(im)%im_myrank

! First, copy entire data type from global index to subdomain index

      itab_m(im_myrank) = ltab_m(im)

! Next, redefine some individual itab_m members

      itab_m(im_myrank)%imglobe = im

! Reset IM neighbor indices to 1

      itab_m(im_myrank)%itopm = 1
      itab_m(im_myrank)%iw(1:npoly) = 1
      itab_m(im_myrank)%iu(1:npoly) = 1

! Global indices of neighbors of IM
! Set indices of neighbors of IM that are present on this rank

      itopm = ltab_m(im)%itopm
      if (myrankflag_m(itopm)) itab_m(im_myrank)%itopm = itabg_m(itopm)%im_myrank

      do j = 1,npoly
         iu = ltab_m(im)%iu(j)
         iw = ltab_m(im)%iw(j)
      
         if (myrankflag_u(iu)) itab_m(im_myrank)%iu(j) = itabg_u(iu)%iu_myrank
         if (myrankflag_w(iw)) itab_m(im_myrank)%iw(j) = itabg_w(iw)%iw_myrank
      enddo
      
! Copy M point grid values

      lpm(im_myrank) = lpm_temp(im)

      xem(im_myrank) = xem_temp(im)
      yem(im_myrank) = yem_temp(im)
      zem(im_myrank) = zem_temp(im)

      topm(im_myrank) = topm_temp(im)

      glatm(im_myrank) = glatm_temp(im)
      glonm(im_myrank) = glonm_temp(im)

      arm0(im_myrank) = arm0_temp(im)
   endif
enddo

! Loop over all U points and select those that are in memory of this subdomain

do iu = 1,nua
   if (myrankflag_u(iu)) then
   
! Look up global IUP index and subdomain IU index

      iup = ltab_u(iu)%iup
      iu_myrank = itabg_u(iu)%iu_myrank

! First, copy entire data type from global index to subdomain index

      itab_u(iu_myrank) = ltab_u(iu)

! Next, redefine some individual itab_u members

      itab_u(iu_myrank)%irank   = itabg_u(iu)%irank
      itab_u(iu_myrank)%iuglobe = iu

! Check if this U point is primary on a remote rank

      if (itab_u(iu_myrank)%irank /= myrank) then      

! Turn off some loop flags for these points (some will be turned back on later)

         call uloops('n',iu_myrank,-7,-8,-12,-13,-16,-21,-22,-23, 0, 0)

! Turn off LBC copy (n/a for global domain) if IUP point is on remote node 

         iup = ltab_u(iu)%iup

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
         imn = ltab_u(iu)%im(j)
         if (myrankflag_m(imn)) itab_u(iu_myrank)%im(j) = itabg_m(imn)%im_myrank
      enddo

      do j = 1,12
         iun = ltab_u(iu)%iu(j)
         if (myrankflag_u(iun)) itab_u(iu_myrank)%iu(j) = itabg_u(iun)%iu_myrank
      enddo

      do j = 1,6
         iwn = ltab_u(iu)%iw(j)
         if (myrankflag_w(iwn)) itab_u(iu_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
      enddo

! Copy U point grid values

      lpu(iu_myrank) = lpu_temp(iu)
      lcu(iu_myrank) = lcu_temp(iu)

      xeu(iu_myrank) = xeu_temp(iu)
      yeu(iu_myrank) = yeu_temp(iu)
      zeu(iu_myrank) = zeu_temp(iu)

      unx(iu_myrank) = unx_temp(iu)
      uny(iu_myrank) = uny_temp(iu)
      unz(iu_myrank) = unz_temp(iu)

      vnx(iu_myrank) = vnx_temp(iu)
      vny(iu_myrank) = vny_temp(iu)
      vnz(iu_myrank) = vnz_temp(iu)

      dnu(iu_myrank) = dnu_temp(iu)
      dniu(iu_myrank) = dniu_temp(iu)
      
      dnv(iu_myrank) = dnv_temp(iu)
      dniv(iu_myrank) = dniv_temp(iu)
      
      glatu(iu_myrank) = glatu_temp(iu)
      glonu(iu_myrank) = glonu_temp(iu)

      do k = 1,nza
         aru  (k,iu_myrank) = aru_temp (k,iu)
         volui(k,iu_myrank) = volui_temp(k,iu)
      enddo
   endif
enddo

! Loop over all W points and select those that are in memory of this subdomain

do iw = 1,nwa
   if (myrankflag_w(iw)) then

! Look up global IWP index and subdomain IW index
   
      iwp = ltab_w(iw)%iwp
      iw_myrank = itabg_w(iw)%iw_myrank

! First, copy entire data type from global index to subdomain index

      itab_w(iw_myrank) = ltab_w(iw)

! Next, redefine some individual itab_w members

      itab_w(iw_myrank)%irank   = itabg_w(iw)%irank
      itab_w(iw_myrank)%iwglobe = iw

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
         imn = ltab_w(iw)%im(j)
         if (myrankflag_m(imn)) itab_w(iw_myrank)%im(j) = itabg_m(imn)%im_myrank
      enddo
      
      do j = 1,9
         iun = ltab_w(iw)%iu(j)
         if (myrankflag_u(iun)) itab_w(iw_myrank)%iu(j) = itabg_u(iun)%iu_myrank
      enddo
      
      do j = 1,9
         iwn = ltab_w(iw)%iw(j)
         if (myrankflag_w(iwn)) itab_w(iw_myrank)%iw(j) = itabg_w(iwn)%iw_myrank
      enddo
      
! Copy W point grid values

      lpw(iw_myrank) = lpw_temp(iw)
      lsw(iw_myrank) = lsw_temp(iw)

      topw(iw_myrank) = topw_temp(iw)

      xew(iw_myrank) = xew_temp(iw)
      yew(iw_myrank) = yew_temp(iw)
      zew(iw_myrank) = zew_temp(iw)

      wnx(iw_myrank) = wnx_temp(iw)
      wny(iw_myrank) = wny_temp(iw)
      wnz(iw_myrank) = wnz_temp(iw)

      glatw(iw_myrank) = glatw_temp(iw)
      glonw(iw_myrank) = glonw_temp(iw)

      arw0(iw_myrank) = arw0_temp(iw)

      do k = 1,nza
         arw  (k,iw_myrank) = arw_temp  (k,iw)
         volwi(k,iw_myrank) = volwi_temp(k,iw)
         volt (k,iw_myrank) = volt_temp (k,iw)
         volti(k,iw_myrank) = volti_temp(k,iw)
      enddo
      
   endif
enddo

! Turn back on some loop flags needed with respect to the local IW point

do iw = 1,nwa
   if (itabg_w(iw)%irank == myrank) then

! Set uloop flag 22 for the 3 nearest U neighbors if IW is primary
! on this rank.

      do j=1,3
         iun = ltab_w(iw)%iu(j)
         call uloops('n',itabg_u(iun)%iu_myrank,22,0,0,0,0,0,0,0,0,0)
      enddo

! Set uloop flag 21 for the 9 nearest U neighbors if IW is primary
! on this rank (all the U points of the 3 neighboring triangles)

      do j=1,9
         iun = ltab_w(iw)%iu(j)
         call uloops('n',itabg_u(iun)%iu_myrank,21,0,0,0,0,0,0,0,0,0)
      enddo

! Set wloop flag 28 for the 3 nearest W neighbors if IW is primary
! on this rank

      do j=1,3
         iwn = ltab_w(iw)%iw(j)
         call wloops('n',itabg_w(iwn)%iw_myrank,28,0,0,0,0,0,0,0,0,0)
      enddo
     
   endif
enddo

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
   landflux(1)%iw   = 1
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

      iup  = ltab_u(iu)%iup 
      if (itabg_u(iup)%irank == myrank) call send_table_u(iup,itabg_u(iu)%irank)

      do j = 1,12
         iun = ltab_u(iu)%iu(j) 
         if (itabg_u(iun)%irank == myrank) call send_table_u(iun,itabg_u(iu)%irank)
      enddo         

      do j = 1,6
         iwn = ltab_u(iu)%iw(j)
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

      iwp = ltab_w(iw)%iwp
      if (itabg_w(iwp)%irank == myrank) call send_table_w(iwp,itabg_w(iw)%irank)

      do j = 1,9
         iun = ltab_w(iw)%iu(j) 
         if (itabg_u(iun)%irank == myrank) call send_table_u(iun,itabg_w(iw)%irank)
      enddo

      do j = 1,9
         iwn = ltab_w(iw)%iw(j) 
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_w(iw)%irank)
      enddo

! Special for Zalesak monotonic advection:

      do j = 1,3
         iwn = ltab_w(iw)%iw(j) 
         do k=1,9
            iun = ltab_w(iwn)%iu(k)
            if (itabg_u(iun)%irank == myrank) then
               call send_table_u(iun,itabg_w(iw)%irank)
            endif
         enddo
      enddo

   endif
enddo

! Deallocate temporary data structures and arrays

deallocate (ltab_m,ltab_u,ltab_w)
deallocate (landflux_temp, seaflux_temp)

return
end subroutine para_init

!===============================================================================

subroutine recv_table_u(iremote)

use mem_para,  only: nrecvs_u, recv_u, recv_uf
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
recv_uf(jrecv)%iremote = iremote

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
use mem_para,   only: nsends_u, send_u, send_uf, mgroupsize
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
send_uf(jsend)%iremote = iremote

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

