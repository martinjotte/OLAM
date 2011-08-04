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

use mem_ijtabs, only: itab_m,      itab_v,      itab_w,      &
                      itab_m_vars, itab_v_vars, itab_w_vars, &
                      itabg_m,     itabg_v,     itabg_w,     &
                      alloc_itabs, mrls

use mem_grid,   only: nza, nma, nua, nva, nwa, mma, mua, mva, mwa, &
                      lpm, lpv, lcv, lpw, lsw, &
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
logical :: myrankflag_v(nva) ! Flag for V points existing on this rank
logical :: myrankflag_w(nwa) ! Flag for W points existing on this rank

logical :: seaflag(nws)
logical :: landflag(nwl)

integer :: lpm_temp(nma)
integer :: lpv_temp(nva),lcv_temp(nva)
integer :: lpw_temp(nwa),lsw_temp(nwa)

real :: topm_temp(nma), glatm_temp(nma), glonm_temp(nma), arm0_temp(nma)

real :: xem_temp(nma),yem_temp(nma),zem_temp(nma)
real :: xev_temp(nva),yev_temp(nva),zev_temp(nva)
real :: xew_temp(nwa),yew_temp(nwa),zew_temp(nwa)
real :: unx_temp(nva),uny_temp(nva),unz_temp(nva)
real :: vnx_temp(nva),vny_temp(nva),vnz_temp(nva)
real :: wnx_temp(nwa),wny_temp(nwa),wnz_temp(nwa)

real :: glatw_temp(nwa),glonw_temp(nwa),arw0_temp(nwa),topw_temp(nwa)
real :: dnu_temp(nva),dniu_temp(nva)
real :: dnv_temp(nva),dniv_temp(nva)

real :: aru_temp(nza,nva),arv_temp(nza,nva),arw_temp(nza,nwa)
real :: glatv_temp(nva),glonv_temp(nva)

real :: volvi_temp(nza,nva),volwi_temp(nza,nwa)

real(kind=8) :: volt_temp(nza,nwa),volti_temp(nza,nwa)

! Temporary datatypes

type (itab_m_vars), allocatable :: ltab_m(:)
type (itab_v_vars), allocatable :: ltab_v(:)
type (itab_w_vars), allocatable :: ltab_w(:)

type(flux_vars), allocatable :: landflux_temp(:)
type(flux_vars), allocatable ::  seaflux_temp(:)

! Move data to temporary data structures, nullifying the old datatype

call move_alloc(itab_m, ltab_m)
call move_alloc(itab_v, ltab_v)
call move_alloc(itab_w, ltab_w)

if (isfcl == 1) then
   call move_alloc(landflux, landflux_temp)
   call move_alloc( seaflux,  seaflux_temp)
endif

! Allocate send & recv counter arrays and initialize to zero

allocate (nsends_v(mrls)) ; nsends_v(1:mrls) = 0
allocate (nsends_w(mrls)) ; nsends_w(1:mrls) = 0
allocate (nrecvs_v(mrls)) ; nrecvs_v(1:mrls) = 0
allocate (nrecvs_w(mrls)) ; nrecvs_w(1:mrls) = 0

! Initialize myrank flag arrays to .false.

if (isfcl == 1) then
   landflag(:) = .false.
   seaflag (:) = .false.
endif

myrankflag_m(:) = .false.
myrankflag_v(:) = .false.
myrankflag_w(:) = .false.

! Copy grid coordinates to temporary arrays

do im = 1,nma
   lpm_temp(im) = lpm(im)

   topm_temp(im) = topm(im)

   xem_temp(im)  = xem(im)
   yem_temp(im)  = yem(im)
   zem_temp(im)  = zem(im)

   glatm_temp(im) = glatm(im)
   glonm_temp(im) = glonm(im)

   arm0_temp(im) = arm0(im)
enddo

do iv = 1,nva
   lpv_temp(iv) = lpv(iv)
   lcv_temp(iv) = lcv(iv)

   xev_temp(iv) = xev(iv)
   yev_temp(iv) = yev(iv)
   zev_temp(iv) = zev(iv)
   
   unx_temp(iv) = unx(iv)
   uny_temp(iv) = uny(iv)
   unz_temp(iv) = unz(iv)

   vnx_temp(iv) = vnx(iv)
   vny_temp(iv) = vny(iv)
   vnz_temp(iv) = vnz(iv)

   dnu_temp(iv)  = dnu(iv)
   dniu_temp(iv) = dniu(iv)

   dnv_temp(iv)  = dnv(iv)
   dniv_temp(iv) = dniv(iv)
   
   glatv_temp(iv) = glatv(iv)
   glonv_temp(iv) = glonv(iv)

   do k = 1,nza
      aru_temp(k,iv) = aru(k,iv)
      arv_temp(k,iv) = arv(k,iv)
      volvi_temp(k,iv) = volvi(k,iv)
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

deallocate (lpm, lpv, lcv, lpw, lsw)
deallocate (xem, yem, zem, topm, glatm, glonm, arm0)
deallocate (xev, yev, zev, unx, uny, unz, vnx, vny, vnz, dnu, dnv, dniu, dniv)
deallocate (aru, arv, volvi, glatv, glonv)
deallocate (xew, yew, zew, wnx, wny, wnz, glatw, glonw, arw0, topw)
deallocate (arw, volwi, volt, volti)

! Loop over all V points, and for each whose assigned irank is equal to myrank,
! flag all V and W points in its computational stencil for inclusion on this
! rank, excluding IVP and IWP.

do iv = 2,nva

   if (itabg_v(iv)%irank == myrank) then

      myrankflag_v(iv) = .true.

      myrankflag_v( ltab_v(iv)%iv(1:16) ) = .true.
      myrankflag_w( ltab_v(iv)%iw(1:4)  ) = .true.

   endif
enddo

! Loop over all V points, and for each that has been flagged for inclusion
! on this rank, flag its M points for inclusion on this rank.
! Count V points also.

do iv = 2,nva

   if (myrankflag_v(iv)) then

      do j = 1,6
         imn = ltab_v(iv)%im(j) 
         myrankflag_m(imn) = .true.
      enddo

      iv_myrank = iv_myrank + 1

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

! Set mma, mva, mwa values for this rank

mma = im_myrank
mva = iv_myrank
mua = mva
mwa = iw_myrank

! Allocate itab data structures and main grid coordinate arrays

call alloc_itabs(meshtype,mma,mua,mva,mwa)
call alloc_xyzem(mma)
call alloc_xyzew(mwa)
call alloc_grid1(meshtype)
call alloc_grid2(meshtype)

! Reset point counts to 1

im_myrank = 1
iv_myrank = 1
iw_myrank = 1

! Store new myrank M, V, W indices in itabg data structures

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
      
! Fill itabg_v(iv)%iv_myrank value for IVP point

      ivp = ltab_v(iv)%ivp
      itabg_v(ivp)%iv_myrank = iv_myrank
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
      itab_m(im_myrank)%iv(1:npoly) = 1

! Global indices of neighbors of IM
! Set indices of neighbors of IM that are present on this rank

      itopm = ltab_m(im)%itopm
      if (myrankflag_m(itopm)) itab_m(im_myrank)%itopm = itabg_m(itopm)%im_myrank

      do j = 1,npoly
         iv = ltab_m(im)%iv(j)
         iw = ltab_m(im)%iw(j)
      
         if (myrankflag_v(iv)) itab_m(im_myrank)%iv(j) = itabg_v(iv)%iv_myrank
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

! Loop over all V points and select those that are in memory of this subdomain

do iv = 1,nva
   if (myrankflag_v(iv)) then

! Look up global IVP index and subdomain IV index

      ivp = ltab_v(iv)%ivp
      iv_myrank = itabg_v(iv)%iv_myrank

! First, copy entire data type from global index to subdomain index

      itab_v(iv_myrank) = ltab_v(iv)

! Next, redefine some individual itab_v members

      itab_v(iv_myrank)%irank   = itabg_v(iv)%irank
      itab_v(iv_myrank)%ivglobe = iv

! Check if this V point is primary on a remote rank

      if (itab_v(iv_myrank)%irank /= myrank) then      

! Turn off some loop flags for these points

         call vloops('n',iv_myrank,-7,-8,-12,-15,-16,-18,  0, 0, 0, 0)

! Turn off LBC copy (n/a for global domain) if IVP point is on remote node

         if (.not. myrankflag_v(ivp))  &
            call vloops('n',iv_myrank,-9,-18,0,0,0,0,0,0,0)

      endif

! Reset IV neighbor indices to 1

      itab_v(iv_myrank)%ivp = 1

      itab_v(iv_myrank)%im(1:6)  = 1
      itab_v(iv_myrank)%iv(1:16) = 1
      itab_v(iv_myrank)%iw(1:4)  = 1


! Set indices of neighbors of IV that are present on this rank

      if (myrankflag_v(ivp)) itab_v(iv_myrank)%ivp = itabg_v(ivp)%iv_myrank

      do j = 1,6
         imn = ltab_v(iv)%im(j)
         if (myrankflag_m(imn)) itab_v(iv_myrank)%im(j) = itabg_m(imn)%im_myrank
      enddo

      do j = 1,16
         ivn = ltab_v(iv)%iv(j)
         if (myrankflag_v(ivn)) itab_v(iv_myrank)%iv(j) = itabg_v(ivn)%iv_myrank
      enddo

      do j = 1,4
         iwn = ltab_v(iv)%iw(j)
         if (myrankflag_w(iwn)) itab_v(iv_myrank)%iw(j) = itabg_w(ivn)%iw_myrank
      enddo

! Copy V point grid values

      lpv(iv_myrank) = lpv_temp(iv)
      lcv(iv_myrank) = lcv_temp(iv)

      xev(iv_myrank) = xev_temp(iv)
      yev(iv_myrank) = yev_temp(iv)
      zev(iv_myrank) = zev_temp(iv)

      unx(iv_myrank) = unx_temp(iv)
      uny(iv_myrank) = uny_temp(iv)
      unz(iv_myrank) = unz_temp(iv)

      vnx(iv_myrank) = vnx_temp(iv)
      vny(iv_myrank) = vny_temp(iv)
      vnz(iv_myrank) = vnz_temp(iv)

      dnu (iv_myrank) = dnu_temp (iv)
      dniu(iv_myrank) = dniu_temp(iv)

      dnv (iv_myrank) = dnv_temp (iv)
      dniv(iv_myrank) = dniv_temp(iv)
      
      glatv(iv_myrank) = glatv_temp(iv)
      glonv(iv_myrank) = glonv_temp(iv)

      do k = 1,nza
         aru  (k,iv_myrank) = aru_temp (k,iv)
         arv  (k,iv_myrank) = arv_temp (k,iv)
         volvi(k,iv_myrank) = volvi_temp(k,iv)
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

! Turn off some loop flags for these points

         call wloops('n',iw_myrank,-12,-13,-15,-16,-17,-19,-20,-23,-26,-27)
         call wloops('n',iw_myrank,-29,-30,-34,  0,  0,  0,  0,  0,  0,  0)

! Turn off LBC copy if IWP point is on remote node

         if (.not. myrankflag_w(iwp))  &
            call wloops('n',iw_myrank,-9,-22,-24,-31,-35,0,0,0,0,0)

      endif

! Reset IW neighbor indices to 1

      itab_w(iw_myrank)%iwp = 1

      itab_w(iw_myrank)%im(1:7) = 1
      itab_w(iw_myrank)%iv(1:7) = 1
      itab_w(iw_myrank)%iw(1:7) = 1

! Set indices of neighbors of IW that are present on this rank

      if (myrankflag_w(iwp)) itab_w(iw_myrank)%iwp = itabg_w(iwp)%iw_myrank

      do j = 1,7
         imn = ltab_w(iw)%im(j)
         if (myrankflag_m(imn)) itab_w(iw_myrank)%im(j) = itabg_m(imn)%im_myrank
      enddo
      
      do j = 1,7
         ivn = ltab_w(iw)%iv(j)
         if (myrankflag_v(ivn)) itab_w(iw_myrank)%iv(j) = itabg_v(ivn)%iv_myrank
      enddo
      
      do j = 1,7
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

! Loop over all V points and for each that is primary on a remote rank, 
! access all V and W points in its stencil.

do iv = 2,nva
   if (itabg_v(iv)%irank /= myrank) then

! IV point is primary on remote rank.  

! If IV point is in memory of myrank, its value must be received from 
! remote rank.  Add that remote rank to receive table.

      if (myrankflag_v(iv)) call recv_table_v(itabg_v(iv)%irank)

! Add to send table any V or W point that is primary on myrank and is 
! in the stencil of IU.  

      ivp  = ltab_v(iv)%ivp 
      if (itabg_v(ivp)%irank == myrank) call send_table_v(ivp,itabg_v(iv)%irank)

      do j = 1,16
         ivn = ltab_v(iu)%iv(j) 
         if (itabg_v(ivn)%irank == myrank) call send_table_v(ivn,itabg_v(iv)%irank)
      enddo         

      do j = 1,4
         iwn = ltab_v(iv)%iw(j) 
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_v(iv)%irank)
      enddo         

   endif
enddo

! Loop over all W points and for each that is primary on a remote rank,
! access all V and W points in its stencil.

do iw = 2,nwa
   if (itabg_w(iw)%irank /= myrank) then

! IW point is primary on remote rank.  

! If IW point is in memory of myrank, it must be received from remote rank.
! Add that remote rank to receive table.

      if (myrankflag_w(iw)) call recv_table_w(itabg_w(iw)%irank)

! Add to send table any V or W point that is primary on myrank and is 
! in the stencil of IW.  

      iwp = ltab_w(iw)%iwp
      if (itabg_w(iwp)%irank == myrank) call send_table_w(iwp,itabg_w(iw)%irank)

      do j = 1,7
         ivn = ltab_w(iw)%iv(j) 
         if (itabg_v(ivn)%irank == myrank) call send_table_v(ivn,itabg_w(iw)%irank)
      enddo         

      do j = 1,7
         iwn = ltab_w(iw)%iw(j) 
         if (itabg_w(iwn)%irank == myrank) call send_table_w(iwn,itabg_w(iw)%irank)
      enddo         

   endif
enddo

! Deallocate temporary data structures and arrays

deallocate (ltab_m,ltab_v,ltab_w)
deallocate (landflux_temp, seaflux_temp)

return
end subroutine para_init

!===============================================================================

subroutine recv_table_v(iremote)

use mem_para,  only: nrecvs_v, recv_v, recv_vf
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
recv_vf(jrecv)%iremote = iremote

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

subroutine send_table_v(iv,iremote)

use mem_ijtabs, only: itab_v, itabg_v, mloops_v
use mem_para,   only: nsends_v, send_v, send_vf, mgroupsize
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

itab_v(iv_myrank)%loop(mloops_v+jsend) = .true.
send_v(jsend)%iremote = iremote
send_vf(jsend)%iremote = iremote

return
end subroutine send_table_v

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

