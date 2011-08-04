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
subroutine para_decomp()

! Decompose global grid into multiple subdomains for parallel computation

use mem_para,   only: mgroupsize, myrank
use misc_coms,  only: io6, meshtype

use mem_ijtabs, only: itab_u, itab_v, itab_w, itabg_m, itabg_u, itabg_v, itabg_w
use mem_grid,   only: nma, nua, nva, nwa, xem, yem, zem

use leaf_coms,  only: nml, nul, nwl, isfcl
use mem_leaf,   only: itab_ul, itab_wl, itabg_ml, itabg_ul, itabg_wl, land

use sea_coms,   only: nms, nus, nws
use mem_sea,    only: itab_us, itab_ws, itabg_ms, itabg_us, itabg_ws, sea

use mem_sflux,  only: nseaflux, nlandflux, seafluxg, landfluxg, landflux, seaflux

implicit none

integer :: im,iu,iv,iw
integer :: im1,im2,im3,iw1,iw2
integer :: igp,jgp
integer :: i, j, ii, jj, iil, jjl, iis, jjs, npoly

integer :: iter,ibin
integer :: ngroups
integer :: numtot, numcut, numcent

real :: xewm (nwa),yewm (nwa),zewm (nwa)
real :: xewml(nwl),yewml(nwl),zewml(nwl)
real :: xewms(nws),yewms(nws),zewms(nws)

real :: xmin,ymin,zmin
real :: xmax,ymax,zmax
real :: xrange, yrange, zrange
real :: cmin, cmax, cmin0, cut

integer :: igsize(mgroupsize)
integer :: nwg   (mgroupsize)
integer :: nwgl  (mgroupsize)
integer :: nwgs  (mgroupsize)

integer :: iwtemp (nwa-1), jwtemp (nwa-1)
integer :: iwltemp(nwl-1), jwltemp(nwl-1)
integer :: iwstemp(nws-1), jwstemp(nws-1)

integer :: num(1002)

real :: val (nwa)
real :: vall(nwl)
real :: vals(nws)

Type grp_var
   integer, allocatable :: iw(:)
   integer, allocatable :: iwl(:)
   integer, allocatable :: iws(:)
End type

type (grp_var) :: grp(mgroupsize)

integer :: iwl_atm_ranks(nwl,mgroupsize)
integer :: niwl_atm(nwl)
integer :: iwl

integer :: iws_atm_ranks(nws,mgroupsize)
integer :: niws_atm(nws)
integer :: iws

integer :: nuv_per_node(0:mgroupsize-1)
integer :: iwcr, iwor

iwl_atm_ranks = -1
niwl_atm = 0

iws_atm_ranks = -1
niws_atm = 0

nuv_per_node = 0

! Allocate permanent itabg data structures

allocate (itabg_m(nma))
allocate (itabg_w(nwa))

if (meshtype == 1) then
   allocate (itabg_u(nua))
else
   allocate (itabg_v(nva))
endif

if (isfcl == 1) then
   allocate (itabg_ml(nml))
   allocate (itabg_ul(nul))
   allocate (itabg_wl(nwl))

   allocate (itabg_ms(nms))
   allocate (itabg_us(nus))
   allocate (itabg_ws(nws))

   allocate (seafluxg(nseaflux))
   allocate (landfluxg(nlandflux))
endif

! Define temp variables for W

do iw = 2,nwa
   xewm(iw) = -1.e9
   yewm(iw) = -1.e9
   zewm(iw) = -1.e9

   do j = 1,itab_w(iw)%npoly
      im = itab_w(iw)%im(j)
      
      if (xewm(iw) < xem(im)) xewm(iw) = xem(im)
      if (yewm(iw) < yem(im)) yewm(iw) = yem(im)
      if (zewm(iw) < zem(im)) zewm(iw) = zem(im)
   enddo
enddo

! Allocate and fill grp%iw, grp%iwl, and grp%iws for group 1

allocate (grp(1)%iw(nwa-1))

do iw = 2,nwa
   grp(1)%iw(iw-1) = iw
enddo

! Check if LEAF and SEA models are being used

if (isfcl == 1) then

! Define temp variables for WL

   do iw = 2,nwl
      xewml(iw) = -1.e9
      yewml(iw) = -1.e9
      zewml(iw) = -1.e9

      do j = 1,itab_wl(iw)%npoly
         im = itab_wl(iw)%im(j)
      
         if (xewml(iw) < land%xem(im)) xewml(iw) = land%xem(im)
         if (yewml(iw) < land%yem(im)) yewml(iw) = land%yem(im)
         if (zewml(iw) < land%zem(im)) zewml(iw) = land%zem(im)
      enddo
   enddo

! Define temp variables for WS

   do iw = 2,nws
      xewms(iw) = -1.e9
      yewms(iw) = -1.e9
      zewms(iw) = -1.e9

      do j = 1,itab_ws(iw)%npoly
         im = itab_ws(iw)%im(j)
      
         if (xewms(iw) < sea%xem(im)) xewms(iw) = sea%xem(im)
         if (yewms(iw) < sea%yem(im)) yewms(iw) = sea%yem(im)
         if (zewms(iw) < sea%zem(im)) zewms(iw) = sea%zem(im)
      enddo
   enddo

! Allocate and fill grp%iwl, and grp%iws for group 1

   allocate (grp(1)%iwl(nwl-1))
   allocate (grp(1)%iws(nws-1))

   do iw = 2,nwl
      grp(1)%iwl(iw-1) = iw
   enddo

   do iw = 2,nws
      grp(1)%iws(iw-1) = iw
   enddo

endif

ngroups = 1
jgp = 1
igsize(1) = mgroupsize

nwg (1) = nwa - 1
nwgl(1) = nwl - 1
nwgs(1) = nws - 1

do while (ngroups < mgroupsize)

   do igp = 1,ngroups

      if (igsize(igp) > 1) then

         jgp = jgp + 1

         xmin = 1.e9
         ymin = 1.e9
         zmin = 1.e9

         xmax = -1.e9
         ymax = -1.e9
         zmax = -1.e9

! Find max and min x,y,z of current group

         do i = 1,nwg(igp)
            iw = grp(igp)%iw(i)

            if (xmin > xewm(iw)) xmin = xewm(iw)
            if (ymin > yewm(iw)) ymin = yewm(iw)
            if (zmin > zewm(iw)) zmin = zewm(iw)

            if (xmax < xewm(iw)) xmax = xewm(iw)
            if (ymax < yewm(iw)) ymax = yewm(iw)
            if (zmax < zewm(iw)) zmax = zewm(iw)
         enddo

! Determine whether to cut in x, y, or z direction

         if (1.01 * (zmax - zmin) > xmax - xmin  .and.  &
             1.01 * (zmax - zmin) > ymax - ymin) then

! ATM cells - z direction

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               val(iw) = zewm(iw)
            enddo

            cmin = zmin
            cmax = zmax

           if (isfcl == 1) then

! LAND cells - z direction

               do i = 1,nwgl(igp)
                  iw = grp(igp)%iwl(i)
                  vall(iw) = zewml(iw)
               enddo

! SEA cells - z direction

               do i = 1,nwgs(igp)
                  iw = grp(igp)%iws(i)
                  vals(iw) = zewms(iw)
               enddo
               
            endif

         elseif (1.01 * (xmax - xmin) > ymax - ymin) then

! ATM cells - x direction

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               val(iw) = xewm(iw)
            enddo

            cmin = xmin
            cmax = xmax

            if (isfcl == 1) then

! LAND cells - x direction

               do i = 1,nwgl(igp)
                  iw = grp(igp)%iwl(i)
                  vall(iw) = xewml(iw)
               enddo

! SEA cells - x direction

               do i = 1,nwgs(igp)
                  iw = grp(igp)%iws(i)
                  vals(iw) = xewms(iw)
               enddo
               
            endif

         else

! ATM cells - y direction

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               val(iw) = yewm(iw)
            enddo

            cmin = ymin
            cmax = ymax

            if (isfcl == 1) then

! LAND cells - y direction

               do i = 1,nwgl(igp)
                  iw = grp(igp)%iwl(i)
                  vall(iw) = yewml(iw)
               enddo

! SEA cells - y direction

               do i = 1,nwgs(igp)
                  iw = grp(igp)%iws(i)
                  vals(iw) = yewms(iw)
               enddo
               
            endif

         endif

! Determine cut value, iterating 3 times

         do iter = 1,3

            num(:) = 0
            numcent = 0

            do i = 1,nwg(igp)
               iw = grp(igp)%iw(i)
               if (val(iw) <= cmin) then
                  ibin = 1
               elseif (val(iw) >= cmax) then
                  ibin = 1002
               else
                  ibin = int(1000. * (val(iw) - cmin) / (cmax - cmin)) + 2
               endif
               num(ibin) = num(ibin) + 1

! Sum points that are beyond geometric center of group

               if (ibin >= 502) numcent = numcent + 1
            enddo

! Set igsize(jgp) based on numcent

            if (iter == 1) then
               igsize(jgp) = (numcent * igsize(igp)) / nwg(igp)
               igsize(jgp) = max(1,min(igsize(igp)-1,igsize(jgp)))

               igsize(igp) = igsize(igp) - igsize(jgp)
               numcut = (nwg(igp) * igsize(igp)) / (igsize(igp) + igsize(jgp))
            endif

! Sum number in each bin until reaching half of total

            numtot = num(1)
            ibin = 1
            do while (numtot + num(ibin+1) <= numcut)
               numtot = numtot + num(ibin+1)
               ibin = ibin + 1
            enddo
            
            cmin0 = cmin + (cmax - cmin) * .001 * (real(ibin) - 1.1)
            cmax = cmin + (cmax - cmin) * .001 * (real(ibin) +  .1)
            cmin = cmin0
         enddo
         cut = .5 * (cmin + cmax)

! Transfer a number of IW points from igp to jgp

         jj = 0
         ii = 0
         do i = 1,nwg(igp)
            iw = grp(igp)%iw(i)

            if (val(iw) > cut) then
               jj = jj + 1
               jwtemp(jj) = iw
            else
               ii = ii + 1
               iwtemp(ii) = iw
            endif
         enddo

         nwg(igp) = ii
         nwg(jgp) = jj

! Deallocate 1 old group and allocate 2 new groups

         deallocate(grp(igp)%iw)
         allocate(grp(igp)%iw(ii))
         allocate(grp(jgp)%iw(jj))

! Fill 2 new groups of ATM cellsfrom temporary arrays

         do j = 1,jj
            grp(jgp)%iw(j)=jwtemp(j)
         enddo

         do i = 1,ii
            grp(igp)%iw(i)=iwtemp(i)
         enddo

         if (isfcl == 1) then

! Transfer a number of IWL points from igp to jgp

            jjl = 0
            iil = 0
            do i = 1,nwgl(igp)
               iw = grp(igp)%iwl(i)

               if (vall(iw) > cut) then
                  jjl = jjl + 1
                  jwltemp(jjl) = iw
               else
                  iil = iil + 1
                  iwltemp(iil) = iw
               endif
            enddo

            nwgl(igp) = iil
            nwgl(jgp) = jjl

! Deallocate 1 old group and allocate 2 new groups

            deallocate(grp(igp)%iwl)
            allocate(grp(igp)%iwl(iil))
            allocate(grp(jgp)%iwl(jjl))

! Fill 2 new groups of LAND cells from temporary arrays

            do j = 1,jjl
               grp(jgp)%iwl(j)=jwltemp(j)
            enddo

            do i = 1,iil
               grp(igp)%iwl(i)=iwltemp(i)
            enddo

! Transfer a number of IWS points from igp to jgp

            jjs = 0
            iis = 0
            do i = 1,nwgs(igp)
               iw = grp(igp)%iws(i)

               if (vals(iw) > cut) then
                  jjs = jjs + 1
                  jwstemp(jjs) = iw
               else
                  iis = iis + 1
                  iwstemp(iis) = iw
               endif
            enddo

            nwgs(igp) = iis
            nwgs(jgp) = jjs

! Deallocate 1 old group and allocate 2 new groups

            deallocate(grp(igp)%iws)
            allocate(grp(igp)%iws(iis))
            allocate(grp(jgp)%iws(jjs))

! Fill 2 new groups of SEA cells from temporary arrays

            do j = 1,jjs
               grp(jgp)%iws(j)=jwstemp(j)
            enddo

            do i = 1,iis
               grp(igp)%iws(i)=iwstemp(i)
            enddo

         endif

      endif
   enddo

   ngroups = jgp

enddo

! Fill irank for each IW point from group array

do igp = 1,ngroups

! ATM cells

   do i = 1,nwg(igp)
      iw = grp(igp)%iw(i)
      itabg_w(iw)%irank = igp - 1
   enddo

   if (isfcl == 1) then

! LAND cells

      do i = 1,nwgl(igp)
         iw = grp(igp)%iwl(i)
         itabg_wl(iw)%irank = igp - 1
      enddo

! SEA cells

      do i = 1,nwgs(igp)
         iw = grp(igp)%iws(i)
         itabg_ws(iw)%irank = igp - 1
      enddo

   endif

enddo

! Temporary:
! If all of the atm cells above a land/sea cell are on the same node,
! put the sfc cell on that same node if it isn't already

do i = 1, nlandflux
   iw  = landflux(i)%iw
   iwl = landflux(i)%iwls
   niwl_atm(iwl) = niwl_atm(iwl) + 1
   j = niwl_atm(iwl)
   iwl_atm_ranks(iwl,j) = itabg_w(iw)%irank
enddo

do i = 1, nseaflux
   iw  = seaflux(i)%iw
   iws = seaflux(i)%iwls
   niws_atm(iws) = niws_atm(iws) + 1
   j = niws_atm(iws)
   iws_atm_ranks(iws,j) = itabg_w(iw)%irank
enddo

do iwl = 1, nwl
   if (niwl_atm(iwl) > 0) then
      j = niwl_atm(iwl)
      if ( all( iwl_atm_ranks(iwl,1:j) == iwl_atm_ranks(iwl,1) )) then
         if (iwl_atm_ranks(iwl,1) /= itabg_wl(iwl)%irank) then
            itabg_wl(iwl)%irank = iwl_atm_ranks(iwl,1)
         endif
      endif
   endif
enddo

do iws = 1, nws
   if (niws_atm(iws) > 0) then
      j = niws_atm(iws)
      if ( all( iws_atm_ranks(iws,1:j) == iws_atm_ranks(iws,1) )) then
         if (iws_atm_ranks(iws,1) /= itabg_ws(iws)%irank) then
            itabg_ws(iws)%irank = iws_atm_ranks(iws,1)
         endif
      endif
   endif
enddo

! Loop over all U/V points and assign its rank to the higher rank of its
! two IW neighbors

! ATM cells

if (meshtype == 1) then

   do iu = 2,nua
      iw1 = itab_u(iu)%iw(1)
      iw2 = itab_u(iu)%iw(2)

      itabg_u(iu)%irank = itabg_w(iw2)%irank
      ! itabg_u(iu)%irank = max(itabg_w(iw1)%irank,itabg_w(iw2)%irank)

      nuv_per_node( itabg_u(iu)%irank ) = nuv_per_node( itabg_u(iu)%irank ) + 1
   enddo

   do iu = 2,nua
      iw1 = itab_u(iu)%iw(1)
      iw2 = itab_u(iu)%iw(2)

      if (itabg_w(iw1)%irank == itabg_w(iw2)%irank) cycle

      iwcr = itabg_u(iu)%irank
      if (itabg_w(iw1)%irank /= iwcr) then
         iwor = itabg_w(iw1)%irank
      else
         iwor = itabg_w(iw2)%irank
      endif

      if (nuv_per_node(iwor) < nuv_per_node(iwcr)-1) then
         itabg_u(iu)%irank = iwor
         nuv_per_node( iwor ) = nuv_per_node( iwor ) + 1
         nuv_per_node( iwcr ) = nuv_per_node( iwcr ) - 1
      endif
   enddo

else

   do iv = 2,nva
      iw1 = itab_v(iv)%iw(1)
      iw2 = itab_v(iv)%iw(2)

      itabg_v(iv)%irank = itabg_w(iw2)%irank
      ! itabg_v(iv)%irank = max(itabg_w(iw1)%irank,itabg_w(iw2)%irank)
      
      nuv_per_node( itabg_v(iv)%irank ) = nuv_per_node( itabg_v(iv)%irank ) + 1
   enddo

   do iv = 2,nva
      iw1 = itab_v(iv)%iw(1)
      iw2 = itab_v(iv)%iw(2)

      if (itabg_w(iw1)%irank == itabg_w(iw2)%irank) cycle

      iwcr = itabg_v(iv)%irank
      if (itabg_w(iw1)%irank /= iwcr) then
         iwor = itabg_w(iw1)%irank
      else
         iwor = itabg_w(iw2)%irank
      endif
  
      if (nuv_per_node(iwor) < nuv_per_node(iwcr)-1) then
         itabg_v(iv)%irank = iwor
         nuv_per_node( iwor ) = nuv_per_node( iwor ) + 1
         nuv_per_node( iwcr ) = nuv_per_node( iwcr ) - 1
      endif
   enddo
   
endif

if (isfcl == 1) then

! LAND cells

   do iu = 2,nul
      iw1 = itab_ul(iu)%iw(1)
      iw2 = itab_ul(iu)%iw(2)

! iw1 or iw2 may be zero at edge of land grid

      if (iw1 < 2) then
         itabg_ul(iu)%irank = itabg_wl(iw2)%irank
      elseif (iw2 < 2) then
         itabg_ul(iu)%irank = itabg_wl(iw1)%irank
      else
         itabg_ul(iu)%irank = max(itabg_wl(iw1)%irank,itabg_wl(iw2)%irank)
      endif
   enddo

! SEA cells

   do iu = 2,nus
      iw1 = itab_us(iu)%iw(1)
      iw2 = itab_us(iu)%iw(2)

! iw1 or iw2 may be zero at edge of sea grid

      if (iw1 < 2) then
         itabg_us(iu)%irank = itabg_ws(iw2)%irank
      elseif (iw2 < 2) then
         itabg_us(iu)%irank = itabg_ws(iw1)%irank
      else
         itabg_us(iu)%irank = max(itabg_ws(iw1)%irank,itabg_ws(iw2)%irank)
      endif
      
   enddo

endif

return
end subroutine para_decomp
