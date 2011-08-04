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
subroutine tileslab_horiz_mp(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mza, mma, mwa, lpw, zm, zt, xem, yem, zem, &
                      xeu, yeu, zeu, xev, yev, zev, xew, yew, zew
use mem_ijtabs, only: itab_m, itab_u
use misc_coms,  only: io6, meshtype

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: kt, k
integer :: ici
integer :: ng
integer :: j, jn, jnn, im, iw, npoly, iv1, iv2, im1, im2
integer :: notavail, navail

real :: hpt, vpt
real :: fldval
real :: xeu_cross, yeu_cross, zeu_cross

real :: htpn(7),vtpn(7)
real :: hqpn(4),vqpn(4)

integer :: ktf(mwa),km(mma)
real :: wtbot(mma),wttop(mma)

! Find cell K indices on the given plot surface

call horizplot_k(iplt,mma,ktf,km,wtbot,wttop)

! If field is 3d, first plot underground points with underground color

if (action == 'T' .and. op%dimens == '3') then
   call plot_underground_w(iplt,ktf)
endif

do im = 2,mma

! Check for points to be skipped over

   if (.not. itab_m(im)%loop(1)) cycle  ! For now, skip pts that don't read in topm

! Get tile plot coordinates.

   call oplot_transform(iplt,xem(im),yem(im),zem(im),hpt,vpt)

   npoly = itab_m(im)%npoly

! Initialize navail counter

   navail = 0

   do j = 1,npoly

      iw = itab_m(im)%iw(j)

! Skip this M point if current IW point index < 2 (which occurs at lateral boundary
! of limited-area domain or parallel subdomain)      

      if (iw < 2) go to 9

      call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),htpn(j),vtpn(j))

! Avoid wrap-around for lat-lon plot

      if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))

! Jump out of loop if any cell corner is on other side of earth

      if (htpn(j) > 1.e11) go to 9

! If any IW point for this IM point is above ground and below model top,
! set notavail flag to zero

      if (ktf(iw) == 0) navail = navail + 1
   enddo

   if (navail == 0) cycle 

! Jump out of loop if entire cell is outside plot window. 

   if ( all(htpn(1:npoly) < op%xmin) .or. &
        all(htpn(1:npoly) > op%xmax) .or. &
        all(vtpn(1:npoly) < op%ymin) .or. &
        all(vtpn(1:npoly) > op%ymax) ) cycle

! Get field value

   call oplot_lib(km(im),im,'VALUE',op%fldname(iplt),wtbot(im),wttop(im), &
                  fldval,notavail)

   if (notavail > 0) cycle 

! If ALL surrounding T cells are available (e.g., above ground), plot entire
! M cell as single polygon.

   if (navail == npoly) then

      call celltile(iplt,im,npoly,htpn,vtpn,hpt,vpt,fldval,action)

! If any surrounding T cell is unavailable, plot M cell by avalable sectors

   else

      hqpn(1) = hpt
      vqpn(1) = vpt

      do j = 1,npoly

         iw = itab_m(im)%iw(j)

         if (ktf(iw) == 0) then
            jn = j + 1
            if (jn > npoly) jn = 1
            jnn = jn + 1
            if (jnn > npoly) jnn = 1

            hqpn(3) = htpn(j)
            vqpn(3) = vtpn(j)

            if (meshtype == 1) then

! Specific way to get IV1 and IV2 since ordering of W and U/V neighbors of M 
! is not identical for both grid systems

               iv1 = itab_m(im)%iu(j)
               iv2 = itab_m(im)%iu(jn)
               
! Evaluate point of intersection between U and V edges

               im1 = itab_u(iv1)%im(1)
               im2 = itab_u(iv1)%im(2)

               xeu_cross = xem(im1) + (xem(im2) - xem(im1)) * itab_u(iv1)%crossmm
               yeu_cross = yem(im1) + (yem(im2) - yem(im1)) * itab_u(iv1)%crossmm
               zeu_cross = zem(im1) + (zem(im2) - zem(im1)) * itab_u(iv1)%crossmm
               
               call oplot_transform(iplt,xeu_cross,yeu_cross,zeu_cross, &
                                    hqpn(2),vqpn(2))

               im1 = itab_u(iv2)%im(1)
               im2 = itab_u(iv2)%im(2)

               xeu_cross = xem(im1) + (xem(im2) - xem(im1)) * itab_u(iv2)%crossmm
               yeu_cross = yem(im1) + (yem(im2) - yem(im1)) * itab_u(iv2)%crossmm
               zeu_cross = zem(im1) + (zem(im2) - zem(im1)) * itab_u(iv2)%crossmm
               
               call oplot_transform(iplt,xeu_cross,yeu_cross,zeu_cross, &
                                    hqpn(4),vqpn(4))

            else

! Specific way to get IV1 and IV2 since ordering of W and U/V neighbors of M 
! is not identical for both grid systems

               iv1 = itab_m(im)%iv(jn)
               iv2 = itab_m(im)%iv(jnn)
               
               call oplot_transform(iplt,xev(iv1),yev(iv1),zev(iv1),hqpn(2),vqpn(2))
               call oplot_transform(iplt,xev(iv2),yev(iv2),zev(iv2),hqpn(4),vqpn(4))
            endif

            call celltile(iplt,im,4,hqpn,vqpn,hpt,vpt,fldval,action)

         endif

      enddo

   endif

   9 continue     

enddo

return
end subroutine tileslab_horiz_mp

!===============================================================================

subroutine tileslab_horiz_tw(iplt,action)

use oplot_coms, only: op, xepc, yepc, zepc
use mem_grid,   only: mza, mwa, zm, zt, xew, yew, zew, xem, yem, zem, lpw
use mem_ijtabs, only: itab_w
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: npoly, j, im
integer :: kt, k
integer :: iw
integer :: iv1,iv2
integer :: ici
integer :: ng
integer :: iok
integer :: notavail

real :: hpt, vpt
real :: fldval

real :: topo1, topo2

real :: htpn(7), vtpn(7)

integer :: ktf(mwa),kw(mwa)
real :: wtbot(mwa),wttop(mwa)

! Find cell K indices on the given plot surface

call horizplot_k(iplt,mwa,ktf,kw,wtbot,wttop)

! If field is 3d, first plot underground points with underground color

if (action == 'T' .and. op%dimens == '3') then
   call plot_underground_w(iplt,ktf)
endif

do iw = 2,mwa

! Skip this point if it is underground

   if (ktf(iw) /= 0) cycle

! Get tile plot coordinates.  

   call oplot_transform(iplt,xew(iw),yew(iw),zew(iw),hpt,vpt)

   npoly = itab_w(iw)%npoly

   do j = 1,npoly

      im = itab_w(iw)%im(j)
      
      call oplot_transform(iplt,xem(im),yem(im),zem(im),htpn(j),vtpn(j))

! Avoid wrap-around for lat-lon plot

      if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))

! Jump out of loop if cell corner is on other side of earth

      if (htpn(j) > 1.e11) go to 9

   enddo

! Jump out of loop if entire cell is outside plot window

   if ( all(htpn(1:npoly) < op%xmin) .or. all(htpn(1:npoly) > op%xmax) .or.  &
        all(vtpn(1:npoly) < op%ymin) .or. all(vtpn(1:npoly) > op%ymax) ) cycle

! Get cell value and plot if 'available'

   call oplot_lib(kw(iw),iw,'VALUE',op%fldname(iplt),wtbot(iw),wttop(iw), &
                  fldval,notavail)

   if (notavail > 0) cycle 

   call celltile(iplt,iw,npoly,htpn,vtpn,hpt,vpt,fldval,action)

! Plot cone circle if so specified

   if (op%pltcone(iplt) == 'C') then
   
      call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
      
      if (iok == 1) then
      
         call oplot_transform(iplt,xepc(1),yepc(1),zepc(1),htpn(1),vtpn(1))
         call oplot_transform(iplt,xepc(2),yepc(2),zepc(2),htpn(2),vtpn(2))
   
         call o_frstpt(htpn(1),vtpn(1))
         call o_vector(htpn(2),vtpn(2))

      endif
 
   endif
   
   9 continue
      
enddo

return
end subroutine tileslab_horiz_tw

!===============================================================================

subroutine tileslab_horiz_vn(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mza, mva, mwa, zm, xeu, yeu, zeu, xev, yev, zev, &
                      xem, yem, zem, xew, yew, zew, lpw
use mem_ijtabs, only: itab_u, itab_v
use misc_coms,  only: io6, meshtype

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: kt, k
integer :: iv
integer :: iw1, iw2
integer :: itpn
integer :: ici
integer :: ng
integer :: im1, im2
integer :: notavail

real :: fldval
real :: hpt, vpt
real :: x1, x2, x3, x4

real :: htpn(4), vtpn(4)

real :: htpn2,htpn4
real :: vtpn2,vtpn4

integer :: ktf(mwa),kv(mva)
real :: wtbot(mva),wttop(mva)

! Find cell K indices on the given plot surface

call horizplot_k(iplt,mva,ktf,kv,wtbot,wttop)

! If field is 3d, first plot underground points with underground color

if (action == 'T' .and. op%dimens == '3') then
   call plot_underground_w(iplt,ktf)
endif

do iv = 2,mva

! Transform tile plot X and Y coordinates.  

   if (meshtype == 1) then
   
      call oplot_transform(iplt,xeu(iv),yeu(iv),zeu(iv),hpt,vpt)

      im1 = itab_u(iv)%im(1)
      im2 = itab_u(iv)%im(2)
      iw1 = itab_u(iv)%iw(1)
      iw2 = itab_u(iv)%iw(2)

   else

      call oplot_transform(iplt,xev(iv),yev(iv),zev(iv),hpt,vpt)

      im1 = itab_v(iv)%im(1)
      im2 = itab_v(iv)%im(2)
      iw1 = itab_v(iv)%iw(1)
      iw2 = itab_v(iv)%iw(2)

   endif

! Skip this V point if both its W neighbors are below ground

   if (ktf(iw1) /= 0 .and. ktf(iw2) /= 0) cycle

   call oplot_transform(iplt,xem(im1),yem(im1),zem(im1),htpn(1),vtpn(1))
   call oplot_transform(iplt,xem(im2),yem(im2),zem(im2),htpn(3),vtpn(3))

   if (iw2 > 1) then
      call oplot_transform(iplt,xew(iw2),yew(iw2),zew(iw2),htpn(2),vtpn(2))
   else
      htpn(2) = htpn(3)
      vtpn(2) = vtpn(3)
   endif

   if (iw1 > 1) then
      call oplot_transform(iplt,xew(iw1),yew(iw1),zew(iw1),htpn(4),vtpn(4))
   else
      htpn(4) = htpn(3)
      vtpn(4) = vtpn(3)
   endif

!  Avoid wrap-around for lat-lon plot

   if (op%projectn(iplt) == 'L') then
      do itpn = 1,4
         call ll_unwrap(hpt,htpn(itpn))
      enddo
   endif

! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:4) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

   if (all(htpn(1:4) < op%xmin) .or. all(htpn(1:4) > op%xmax) .or.  &
       all(vtpn(1:4) < op%ymin) .or. all(vtpn(1:4) > op%ymax)) cycle

! Save copy of 2nd and 4th htpn,vtpn pts

   htpn2 = htpn(2)
   htpn4 = htpn(4)
   vtpn2 = vtpn(2)
   vtpn4 = vtpn(4)
   
! Get cell value and 'available' flag

   call oplot_lib(kv(iv),iv,'VALUE',op%fldname(iplt),wtbot(iv),wttop(iv), &
                  fldval,notavail)    

   if (notavail > 0) cycle 

! If both W neighbors are above ground, plot full cell

   if (ktf(iw1) == 0 .and. ktf(iw2) == 0) then

      call celltile(iplt,iv,4,htpn,vtpn,hpt,vpt,fldval,action)
      
! Else, check if IW1 is above ground

   elseif (ktf(iw1) == 0) then

! Plot IW1 half of cell with tile color

      htpn(4) = htpn4
      vtpn(4) = vtpn4

      htpn(2) = htpn(3)
      vtpn(2) = vtpn(3)

      call celltile(iplt,iv,4,htpn,vtpn,hpt,vpt,fldval,action)
      
   else

! Plot IW2 half of cell with tile color

      htpn(2) = htpn2
      vtpn(2) = vtpn2

      htpn(4) = htpn(3)
      vtpn(4) = vtpn(3)

      call celltile(iplt,iv,4,htpn,vtpn,hpt,vpt,fldval,action)

   endif

enddo

return
end subroutine tileslab_horiz_vn

!===============================================================================

subroutine tileslab_horiz_s(iplt,action)

use oplot_coms, only: op
use mem_sea,    only: sea, itab_ws
use sea_coms,   only: mws
use misc_coms,  only: io6
use max_dims,   only: maxnlspoly

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: npoly
integer :: itpn
integer :: ici
integer :: ng
integer :: ipm, jpm
integer :: notavail
integer :: iws
integer :: nspoly, jms, ims

real :: hpt, vpt
real :: fldval

real :: htpn(maxnlspoly), vtpn(maxnlspoly)
real :: wtbot = 1., wttop = 0.

do iws = 2,mws

   nspoly = itab_ws(iws)%npoly

! Initialize hpt and vpt to zero

   hpt = 0.
   vpt = 0.

! Get tile plot coordinates.  

   do jms = 1,nspoly

      ims = itab_ws(iws)%im(jms)
      
      call oplot_transform(iplt          &
                          ,sea%xem(ims)  &
                          ,sea%yem(ims)  &
                          ,sea%zem(ims)  &
                          ,htpn(jms)     &
                          ,vtpn(jms)     )

      hpt = hpt + htpn(jms)
      vpt = vpt + vtpn(jms)
   enddo

   hpt = hpt / real(nspoly)
   vpt = vpt / real(nspoly)

! Avoid wrap-around for lat-lon plot

   if (op%projectn(iplt) == 'L') then
      do jms = 1,nspoly
         call ll_unwrap(hpt,htpn(jms))
      enddo
   endif
   
! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:nspoly) > 1.e11)) cycle

! Jump out of loop if any is cell corner is outside plot window. 

   if ( all(htpn(1:nspoly) < op%xmin) .or. all(htpn(1:nspoly) > op%xmax) .or. &
        all(vtpn(1:nspoly) < op%ymin) .or. all(vtpn(1:nspoly) > op%ymax) ) cycle

! Plot cell

   call oplot_lib(1,iws,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)    
   if (notavail > 0) cycle 
   call celltile(iplt,iws,nspoly,htpn,vtpn,hpt,vpt,fldval,action)

enddo

return
end subroutine tileslab_horiz_s

!===============================================================================

subroutine tileslab_horiz_l(iplt,action)

use oplot_coms, only: op
use mem_leaf,   only: land, itab_wl
use leaf_coms,  only: mwl, nzg, nzs
use misc_coms,  only: io6
use max_dims,   only: maxnlspoly

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: npoly
integer :: ici
integer :: ng
integer :: ipm
integer :: k
integer :: notavail
integer :: iwl
integer :: nlpoly, jml, iml

real :: hpt, vpt
real :: fldval
real :: wtbot = 1., wttop = 0.

real :: htpn(maxnlspoly), vtpn(maxnlspoly)

! Find K level to plot if field is 3d

if (op%dimens == '3G') then
   k = min(nzg,max(1,nint(op%slabloc(iplt))))
elseif (op%dimens == '3S') then
   k = min(nzs,max(1,nint(op%slabloc(iplt))))
else
   k = 1
endif

do iwl = 2,mwl

   nlpoly = itab_wl(iwl)%npoly

! Initialize hpt and vpt to zero

   hpt = 0.
   vpt = 0.

! Get tile plot coordinates.  

   do jml = 1,nlpoly

      iml = itab_wl(iwl)%im(jml)

      call oplot_transform(iplt           &
                          ,land%xem(iml)  &
                          ,land%yem(iml)  &
                          ,land%zem(iml)  &
                          ,htpn(jml)      &
                          ,vtpn(jml)      )

      hpt = hpt + htpn(jml)
      vpt = vpt + vtpn(jml)

   enddo

   hpt = hpt / real(nlpoly)
   vpt = vpt / real(nlpoly)
   
! Avoid wrap-around for lat-lon plots

   if (op%projectn(iplt) == 'L') then
      do jml = 1,nlpoly
         call ll_unwrap(hpt,htpn(jml))
      enddo
   endif
   
! Jump out of loop if any cell corner is on other side of earth

   if (any(htpn(1:nlpoly) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

   if ( all(htpn(1:nlpoly) < op%xmin) .or. all(htpn(1:nlpoly) > op%xmax) .or.  &
        all(vtpn(1:nlpoly) < op%ymin) .or. all(vtpn(1:nlpoly) > op%ymax) ) cycle

! Plot cell

   call oplot_lib(k,iwl,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)    
   
   call celltile(iplt,iwl,nlpoly,htpn,vtpn,hpt,vpt,fldval,action)

enddo

return
end subroutine tileslab_horiz_l

!===============================================================================

subroutine tileslab_horiz_fs(iplt,action)

use oplot_coms, only: op
use mem_sea,    only: sea, itab_ws
use leaf_coms,  only: mwl, nzg, nzs
use misc_coms,  only: io6
use mem_sflux,  only: mseaflux,seaflux,xemsfpat,yemsfpat,zemsfpat,nsfpatm

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ici
integer :: ng
integer :: ipm
integer :: k
integer :: notavail
integer :: iws
integer :: j
integer :: ims, jms
integer :: isf,isfpat,jsfpats,ksfpat,jsfpatm
integer :: isfprev

real :: hpt, vpt
real :: fldval

real :: htpn(5), vtpn(5)
real :: wtbot = 1., wttop = 0.

k = 1
isfprev = 0

do isf = 2,mseaflux

   jsfpats = seaflux(isf)%jpats
   isfpat  = seaflux(isf)%ipat

   call oplot_transform(iplt              &
                       ,seaflux(isf)%xef  &
                       ,seaflux(isf)%yef  &
                       ,seaflux(isf)%zef  &
                       ,hpt               &
                       ,vpt               )

   do ksfpat = isfpat,isfpat + jsfpats - 1
   
      jsfpatm = nsfpatm(ksfpat)

      do j = 1,jsfpatm
         call oplot_transform(iplt               &
                             ,xemsfpat(j,ksfpat)  &
                             ,yemsfpat(j,ksfpat)  &
                             ,zemsfpat(j,ksfpat)  &
                             ,htpn(j)            &
                             ,vtpn(j)            )

! Avoid wrap-around for lat-lon plots

         if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))

      enddo

! Jump out of loop if any cell corner is on other side of earth

      if (any(htpn(1:jsfpatm) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

      if ( all(htpn(1:jsfpatm) < op%xmin) .or.  &
           all(htpn(1:jsfpatm) > op%xmax) .or.  &
           all(vtpn(1:jsfpatm) < op%ymin) .or.  &
           all(vtpn(1:jsfpatm) > op%ymax) ) cycle

      if (isfprev /= isf) then
         isfprev = isf
         
         call oplot_lib(k,isf,'VALUE',op%fldname(iplt),wtbot,wttop, &
                        fldval,notavail)   
      endif
      
      call celltile(iplt,isf,jsfpatm,htpn,vtpn,hpt,vpt,fldval,action)
      if (action == 'P') cycle

   enddo

enddo

return
end subroutine tileslab_horiz_fs

!===============================================================================

subroutine tileslab_horiz_fl(iplt,action)

use oplot_coms, only: op
use leaf_coms,  only: mwl, nzg, nzs
use misc_coms,  only: io6
use mem_sflux,  only: mlandflux,landflux,xemlfpat,yemlfpat,zemlfpat,nlfpatm

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: ici
integer :: ng
integer :: ipm
integer :: k
integer :: notavail
integer :: iwl
integer :: j
integer :: iml, jml
integer :: ilf,ilfpat,jlfpats,klfpat,jlfpatm
integer :: ilfprev

real :: hpt,  vpt
real :: fldval

real :: htpn(5), vtpn(5)
real :: wtbot = 1., wttop = 0.

k = 1
ilfprev = 0

do ilf = 2,mlandflux

   jlfpats = landflux(ilf)%jpats
   ilfpat  = landflux(ilf)%ipat
   
   call oplot_transform(iplt               &
                       ,landflux(ilf)%xef  &
                       ,landflux(ilf)%yef  &
                       ,landflux(ilf)%zef  &
                       ,hpt                &
                       ,vpt                )


   do klfpat = ilfpat,ilfpat + jlfpats - 1

      jlfpatm = nlfpatm(klfpat)

      do j = 1,jlfpatm
         call oplot_transform(iplt               &
                             ,xemlfpat(j,klfpat)  &
                             ,yemlfpat(j,klfpat)  &
                             ,zemlfpat(j,klfpat)  &
                             ,htpn(j)            &
                             ,vtpn(j)            )

          if (op%projectn(iplt) == 'L') call ll_unwrap(hpt,htpn(j))
      enddo

! Jump out of loop if any cell corner is on other side of earth

      if (any(htpn(1:jlfpatm) > 1.e11)) cycle

! Jump out of loop if entire cell is outside plot window. 

      if ( all(htpn(1:jlfpatm) < op%xmin) .or.  &
           all(htpn(1:jlfpatm) > op%xmax) .or.  &
           all(vtpn(1:jlfpatm) < op%ymin) .or.  &
           all(vtpn(1:jlfpatm) > op%ymax) ) cycle

      if (ilfprev /= ilf) then
         ilfprev = ilf
         
         call oplot_lib(k,ilf,'VALUE',op%fldname(iplt),wtbot,wttop, &
                        fldval,notavail)
      endif

      call celltile(iplt,ilf,jlfpatm,htpn,vtpn,hpt,vpt,fldval,action)
      if (action == 'P') cycle

   enddo

enddo

return
end subroutine tileslab_horiz_fl

!===============================================================================

subroutine tileslab_vert_t(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, zm, zt, lpw
use misc_coms,  only: io6, meshtype

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k
integer :: iw
integer :: iv1,iv2
integer :: iok
integer :: notavail

real :: hpt,  vpt
real :: fldval

real :: htpn(4), vtpn(4)
real :: topo1,topo2
real :: wtbot = 1., wttop = 0.

! Loop over W points

do iw = 2,mwa

! Get horizontal plot coordinates for this W point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
   elseif (op%projectn(iplt)(1:1) == 'V') then

! temporary IF statement - need to make xyplot_w version for hexagons   
      if (meshtype == 1) then
         call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)
      endif
! end temporary      

   endif

   hpt = .5 * (htpn(1) + htpn(2))

! Jump out of loop if this W point does not intersect plot cone

   if (iok /= 1) cycle

! Jump out of loop if either cell side is outside plot window. 

   if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
       htpn(2) < op%xmin .or. htpn(2) > op%xmax) cycle
   
   do k = 2,mza-1

! Get T-cell vertical coordinates

      vtpn(1) = zm(k-1)
      vtpn(2) = vtpn(1)
      vtpn(3) = zm(k)
      vtpn(4) = vtpn(3)
      vpt = zt(k)

! Check if cell is above ground

      if (k >= lpw(iw)) then 

! Cell is above ground

         call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                        fldval,notavail)
         if (notavail > 0) cycle 

         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

      else
           
! Cell is underground - Tile-plot it with underground color
  
         call fillpolyg(4,htpn(1),vtpn(1),op%icigrnd)
      endif

   enddo

enddo
   
return
end subroutine tileslab_vert_t

!===============================================================================

subroutine tileslab_vert_v(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, zm, zt, lpw
use mem_ijtabs, only: itab_u, itab_v
use misc_coms,  only: io6, meshtype

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k
integer :: iw
integer :: ng
integer :: im1,im2
integer :: iw1,iw2
integer :: iv1,iv2
integer :: iok
integer :: notavail

real :: hpt,vpt
real :: fldval

real :: hcpn(4),fldvals(2)
real :: htpn(4),vtpn(4)

real :: hpt1,hpt2
real :: topo1,topo2
real :: wtbot = 1., wttop = 0.

! First plot underground T cells with underground color

call plot_underground_w(iplt,(/0/))

! Loop over W points

do iw = 2,mwa

! Get horizontal plot coordinates for IW point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,hcpn)
   elseif (op%projectn(iplt)(1:1) == 'V') then
      call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,hcpn) ! need to fix for hex??
   endif

! Jump out of loop if this IW point does not intersect plot cone

   if (iok /= 1) cycle

! Avoid wrap_around

   if (op%projectn(iplt) == 'C') then
      call ll_unwrap(hcpn(1),hcpn(2))
   endif

   hpt = .5 * (hcpn(1) + hcpn(2))

! Skip current IV point if either cell side is outside plot window. 

   if (hcpn(1) < op%xmin .or. hcpn(1) > op%xmax .or.  &
       hcpn(2) < op%xmin .or. hcpn(2) > op%xmax) cycle
   
   hcpn(3) = hcpn(2)
   hcpn(4) = hcpn(1)
   
   do k = lpw(iw),mza-1  ! Loop is over T levels

! Skip this K point if either upper or lower cell center is outside plot window. 

      if (zm(k-1) < op%ymin .or. zm(k) > op%ymax) cycle
   
! Get vertical coordinates

      vtpn(1) = zm(k-1)
      vtpn(2) = vtpn(1)
      vtpn(3) = zm(k)
      vtpn(4) = vtpn(3)
      vpt = zt(k)

! Get value for left half of cell
   
      call oplot_lib(k,iv1,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldval,notavail)
      if (notavail > 0) cycle 
         
      htpn(1) = hcpn(1)
      htpn(2) = .5 * (hcpn(1) + hcpn(2))
      htpn(3) = htpn(2)
      htpn(4) = htpn(1)

      call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

! Get value for right half of cell
   
      call oplot_lib(k,iv2,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldval,notavail)
      if (notavail > 0) cycle 
         
      htpn(1) = .5 * (hcpn(1) + hcpn(2))
      htpn(2) = hcpn(2)
      htpn(3) = htpn(2)
      htpn(4) = htpn(1)

      call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

   enddo

enddo

return
end subroutine tileslab_vert_v

!===============================================================================

subroutine tileslab_vert_w(iplt,action)

use oplot_coms, only: op
use mem_grid,   only: mwa, mza, zm, zt, lpw
use misc_coms,  only: io6, meshtype

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k
integer :: iw
integer :: ng
integer :: iv1,iv2
integer :: iok
integer :: notavail

real :: hpt,vpt
real :: fldval

real :: htpn(4),vtpn(4)
real :: topo1,topo2
real :: wtbot = 1., wttop = 0.

! Loop over W points

do iw = 2,mwa

! Get horizontal plot coordinates for this W point

   if (op%projectn(iplt) == 'C') then
      call coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)
   elseif (op%projectn(iplt) == 'V') then

! temporary IF statement - need to make xyplot_w version for hexagons   
      if (meshtype == 1) then
         call xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)
      endif
! end temporary 

   endif

   hpt = .5 * (htpn(1) + htpn(2))

! Jump out of loop if this W point does not intersect plot cone

   if (iok /= 1) cycle

! Jump out of loop if either cell side is outside plot window. 

   if (htpn(1) < op%xmin .or. htpn(1) > op%xmax .or.  &
       htpn(2) < op%xmin .or. htpn(2) > op%xmax) cycle
   
   do k = 1,mza-1
      vpt = zm(k)
      
! Check if lower T cell is above ground

      if (k >= lpw(iw)) then 

! Yes - plot full cell with tile color

         vtpn(1) = zt(k)
         vtpn(2) = vtpn(1)
         vtpn(3) = min(zt(k+1),zm(mza-1))
         vtpn(4) = vtpn(3)
         
         call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                        fldval,notavail)
         if (notavail > 0) cycle 
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)

! Else, check if only upper T cell is above ground

      elseif (k+1 >= lpw(iw)) then 

! Yes - plot upper half W cell with tile color

         vtpn(1) = zm(k)
         vtpn(2) = vtpn(1)
         vtpn(3) = zt(k+1)
         vtpn(4) = vtpn(3)

         call oplot_lib(k,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
                        fldval,notavail)
         if (notavail > 0) cycle 
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)
         
      else
      
! At this point, upper T cell is below ground - Plot with underground color

         vtpn(1) = zm(k)
         vtpn(2) = vtpn(1)
         vtpn(3) = zm(k+1)
         vtpn(4) = vtpn(3)

         call fillpolyg(4,htpn(1),vtpn(1),op%icigrnd)

      endif

   enddo

! Plot top half W cell with tile color

!   vtpn(1) = zt(mza-1)
!   vtpn(2) = vtpn(1)
!   vtpn(3) = zm(mza-1)
!   vtpn(4) = vtpn(3)
!   vpt = zm(mza-1)

!   call oplot_lib(mza-1,iw,'VALUE',op%fldname(iplt),wtbot,wttop, &
!                  fldval,notavail)
!   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval)
      
enddo
   
return
end subroutine tileslab_vert_w

!===============================================================================

subroutine tileslab_vert_l(iplt,action)

use oplot_coms, only: op
use leaf_coms,  only: mwl, nzg, nzs
use mem_leaf,   only: land
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k
integer :: iw
integer :: ng
integer :: iv1, iv2
integer :: iu1, iu2
integer :: ip
integer :: itpn
integer :: notavail

real :: hpt, vpt
real :: fldval
real :: patwidth
real :: botk
real :: delzk

real :: htpn(4), vtpn(4)
real :: wtbot = 1., wttop = 0.

!!!!!!!! VERTICAL XSECTION NEEDS WORK

do ip = 2,mwl
!   iw = land%iw(ip)

! Skip iw column if not intersected by slabloc or if outside window bounds


! Get horizontal plot coordinates for cells in this column

   op%stagpt = 'LA'  ! Get land cell fractional area
   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)

!   patwidth = (xu(iu2) - xu(iu1)) * fldval

!   htpn(1) = vt2da(iw)
!   htpn(2) = htpn(1) + patwidth
!   htpn(3) = htpn(2)
!   htpn(4) = htpn(1)

!   vt2da(iw) = htpn(2)

!   hpt = .5 * (htpn(1) + htpn(2))

   botk = - real(nzg+nzs+4)   ! level of bottom of bottom soil layer (negative)
   delzk = op%ymin / botk  ! This is a positive number

   do k = 1,nzg  ! Loop over soil layers

! Get vertical coordinates for soil layers

      vtpn(1) = delzk * (float(k-1) + botk) 
      vtpn(2) = vtpn(1)
      vtpn(3) = delzk * (float(k)   + botk)
      vtpn(4) = vtpn(3)
      vpt = .5 * (vtpn(1) + vtpn(3))

! plot soil layers

      op%stagpt = 'L'  ! Get soil value
      call oplot_lib(k,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                     fldval,notavail)
      call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

      if (op%pltgrid(iplt) == 'G') then
         call o_frstpt(htpn(4),vtpn(4))
         do itpn = 1,4
            call o_vector(htpn(itpn),vtpn(itpn))
         enddo
      endif

   enddo

   do k = 1,nzs  ! Loop over sfcwater layers

! check for existence of sfcwater in current level k

      call oplot_lib(k,ip,'VALUE','SFCWATER_MASS',wtbot,wttop, &
                     fldval,notavail)

      if (fldval > 0.) then
      
! Get vertical coordinates for sfcwater layer k

         vtpn(1) = delzk * (float(nzg+k-1) + .5 + botk) 
         vtpn(2) = vtpn(1)
         vtpn(3) = delzk * (float(nzg+k)   + .5 + botk) 
         vtpn(4) = vtpn(3)
         vpt = .5 * (vtpn(1) + vtpn(3))

! plot sfcwater layer k

         op%stagpt = 'LW'  ! Get sfcwater value
         call oplot_lib(k,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                        fldval,notavail)
         call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

         if (op%pltgrid(iplt) == 'G') then
            call o_frstpt(htpn(4),vtpn(4))
            do itpn = 1,4
               call o_vector(htpn(itpn),vtpn(itpn))
            enddo
         endif
         
      endif
      
   enddo

! plot vegetation layer

   vtpn(1) = delzk * (float(nzg+nzs) + 1. + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg+nzs+1) + 1. + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'LV'  ! Get vegetation value
   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

! plot canopy air layer

   vtpn(1) = delzk * (float(nzg+nzs+1) + 1.5 + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg+nzs+2) + 1.5 + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'LC'  ! Get vegetation value
   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

enddo
   
return
end subroutine tileslab_vert_l

!===============================================================================

subroutine tileslab_vert_s(iplt,action)

use oplot_coms, only: op
use leaf_coms,  only: nzg, nzs
use sea_coms,   only: mws
use mem_sea,    only: sea
use misc_coms,  only: io6

implicit none

integer,      intent(in) :: iplt
character(1), intent(in) :: action

integer :: k
integer :: iw
integer :: ng
integer :: iv1, iv2
integer :: iu1, iu2
integer :: ip
integer :: itpn
integer :: notavail

real :: hpt, vpt
real :: fldval
real :: patwidth
real :: botk
real :: delzk

real :: htpn(4), vtpn(4)
real :: wtbot = 1., wttop = 0.

!!!!!!!! VERTICAL XSECTION NEEDS WORK

do ip = 2,mws
!   iw = sea%iw(ip)

! Skip iw column if not intersected by slabloc or if outside window bounds

! Get horizontal plot coordinates for cells in this column

   op%stagpt = 'SA'  ! Get sea cell fractional area
   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)

   botk = - real(nzg+nzs+4)   ! level of bottom of bottom soil layer (negative)
   delzk = op%ymin / botk  ! This is a positive number

! plot (top) sea layer

   vtpn(1) = delzk * (float(nzg-1) + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg) + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'S'  ! Get vegetation value
   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

! plot canopy air layer

   vtpn(1) = delzk * (float(nzg+nzs+1) + 1.5 + botk) 
   vtpn(2) = vtpn(1)
   vtpn(3) = delzk * (float(nzg+nzs+2) + 1.5 + botk) 
   vtpn(4) = vtpn(3)
   vpt = .5 * (vtpn(1) + vtpn(3))

   op%stagpt = 'SC'  ! Get canopy air value
   call oplot_lib(0,ip,'VALUE',op%fldname(iplt),wtbot,wttop, &
                  fldval,notavail)
   call celltile(iplt,0,4,htpn,vtpn,hpt,vpt,fldval,action)  ! Checks 'TILE'

   if (op%pltgrid(iplt) == 'G') then
      call o_frstpt(htpn(4),vtpn(4))
      do itpn = 1,4
         call o_vector(htpn(itpn),vtpn(itpn))
      enddo
   endif

enddo
   
return
end subroutine tileslab_vert_s

!===============================================================================

subroutine celltile(iplt,i,npoly,htpn,vtpn,hpt,vpt,fldval,action)

use oplot_coms, only: op
use plotcolors, only: clrtab
use misc_coms,  only: io6

implicit none

integer, intent(in) :: iplt
integer, intent(in) :: i
integer, intent(in) :: npoly

real, intent(in) :: htpn(*)
real, intent(in) :: vtpn(*)
real, intent(in) :: hpt
real, intent(in) :: vpt
real, intent(in) :: fldval

character(1), intent(in) :: action

real    :: fldval1
integer :: icolor
integer :: itab
integer :: ival

if (action == 'T') then

   itab = op%icolortab(iplt)

! Cyclic treatment of color palette (used with integer-type data)

   fldval1 = fldval
   if (clrtab(itab)%ifmt(1) == 20)  &
      fldval1 = mod(fldval-1.,real(clrtab(itab)%nvals-2)) + 1.

! Extract contour color from color table

   ival = 1
   do while (fldval1 > clrtab(itab)%vals(ival) .and.  &
               ival < clrtab(itab)%nvals             )
      ival = ival + 1
   enddo
   icolor = clrtab(itab)%ipal(ival)
   
   call fillpolyg(npoly,htpn,vtpn,icolor)

   if (op%fldval_min > fldval) op%fldval_min = fldval
   if (op%fldval_max < fldval) op%fldval_max = fldval

endif

if (action == 'P') then
   if ( (hpt > op%xmin) .and. (hpt < op%xmax) .and. &
        (vpt > op%ymin) .and. (vpt < op%ymax) ) then
      call oplot_prtvalue(fldval,hpt,vpt,op%vsprd,.7*op%psiz,op%icolortab(iplt))
   endif
endif

return
end subroutine celltile

!===============================================================================

subroutine coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)

use oplot_coms,  only: op, xepc, yepc, zepc
use mem_grid,    only: xem, yem, zem, topm
use mem_ijtabs,  only: itab_w
use consts_coms, only: pio180, erad, piu180, pi2
use misc_coms,   only: io6, meshtype

implicit none

integer, intent(in) :: iw
integer, intent(out) :: iv1,iv2
integer, intent(out) :: iok
real, intent(out) :: htpn(4)
real, intent(out) :: topo1,topo2

integer :: npoly,j,jm1,jm2,jm3,jmin,jv,im,im1,im2,iv

real :: unxec3,unyec3,unzec3
real :: xcent,ycent,zcent
real :: xwincent,ywincent,zwincent
real :: radcone
real :: raxis
real :: vecx,vecy,vecz
real :: vecleftx,veclefty,vecleftz
real :: valx,valy
real :: angm(7),angm1,angm2
real :: angmin,angmax
real :: ang0
real :: wt1,wt2
real :: xrad1,yrad1,zrad1
real :: xrad2,yrad2,zrad2
real :: ange1,ange2
real :: sinplat,cosplat
real :: sinvaz,cosvaz
real :: sinconang,cosconang
real :: sinconlat,cosconlat
real :: sinconlon,cosconlon

real :: xem1,xem2,yem1,yem2,zem1,zem2
real :: topm1,topm2

sinplat = sin(op%plat3 * pio180)
cosplat = cos(op%plat3 * pio180)

! cone is viewed from INSIDE, so cone axis is in opposite direction from viewazim 

sinvaz = -sin(op%viewazim * pio180)
cosvaz = -cos(op%viewazim * pio180)

sinconang = sin(op%coneang * pio180)
cosconang = cos(op%coneang * pio180)

! Use formulas (5-5) and (5-6) of USGS Map Projections Manual to get lat/lon
! of plot-cone axis given lat/lon of plot center, plot-cone angle, and view azimuth
! (Avoid special case of where cone center is +/- 90 deg longitude from plot center)

op%conelat = asin(sinplat * cosconang + cosplat * sinconang * cosvaz)

if (abs(cosplat * cosconang - sinplat * sinconang * cosvaz) > 1.e-6) then

   op%conelon = op%plon3 * pio180                                   &
              + atan2(sinconang * sinvaz,                           &
               (cosplat * cosconang - sinplat * sinconang * cosvaz))

elseif (op%viewazim < 180.) then
   op%conelon = (op%plon3 - 90.) * pio180
else
   op%conelon = (op%plon3 + 90.) * pio180
endif

sinconlat = sin(op%conelat)
cosconlat = cos(op%conelat)

sinconlon = sin(op%conelon)
cosconlon = cos(op%conelon)

! Earth components of unit vector outward along cone axis

unxec3 = cosconlat * cosconlon
unyec3 = cosconlat * sinconlon
unzec3 = sinconlat

! Intersection of cone and earth is a circle - find earth coords of circle center

xcent = unxec3 * erad * cosconang 
ycent = unyec3 * erad * cosconang 
zcent = unzec3 * erad * cosconang 

! Cone radius

radcone = erad * sinconang

! Earth coordinates of plot window center

zwincent = erad * sinplat
raxis = erad * cosplat
xwincent = raxis * cos(op%plon3 * pio180) 
ywincent = raxis * sin(op%plon3 * pio180) 

! Initialize ANGMIN and ANGMAX

angmin = 400.
angmax = -400.

npoly = itab_w(iw)%npoly

! Loop over neighbor M points for this IW cell

do j = 1,npoly
   im = itab_w(iw)%im(j)

! Compute angle of M point with plot cone center using dot products 
! (Adequate precision may require cross product for cone angles 
! close to 0 or 180)

   angm(j) = piu180  &
            * acos((xem(im)*unxec3 + yem(im)*unyec3 + zem(im)*unzec3) / erad)

   if (angmin > angm(j)) then
      angmin = angm(j)
      jmin = j
   endif

   if (angmax < angm(j)) then
      angmax = angm(j)
   endif

enddo

! Return with iok = 0 if IW column is not intersected by plot slab

iok = 0

if (angmin > op%coneang .or. angmax < op%coneang) return

iok = 1  ! Since we got here, IW column is intersected by plot cone

! Fill arrays of values at M points in cyclic order around IW column, beginning
! with M point that is closest to cone axis

do j = 1,npoly
   jm1 = j + jmin - 1
   jm2 = jm1 + 1
   jm3 = jm2 + 1

   if (jm1 > npoly) jm1 = jm1 - npoly
   if (jm2 > npoly) jm2 = jm2 - npoly
   if (jm3 > npoly) jm3 = jm3 - npoly
   
   im1 = itab_w(iw)%im(jm1)
   im2 = itab_w(iw)%im(jm2)

   if (meshtype == 1) then
      jv = jm3
      iv  = itab_w(iw)%iu(jv)
   else
      jv = jm2
      iv  = itab_w(iw)%iv(jv)
   endif
      
   xem1 = xem(im1)
   yem1 = yem(im1)
   zem1 = zem(im1)

   xem2 = xem(im2)
   yem2 = yem(im2)
   zem2 = zem(im2)
   
   topm1 = topm(im1)
   topm2 = topm(im2)
   
   angm1 = angm(jm1)
   angm2 = angm(jm2)

! Find two points of intersection between current IW polygon and cone

   if (angm1 <= op%coneang .and. angm2 >= op%coneang) then

! This interval touches cone

      if (angm1 == angm2) then

! This interval is tangent to cone      

         if (jm1 == 1 .or. jm2 == 1) then

! This interval is on minimum side of polygon

            xepc(1) = xem1
            yepc(1) = yem1
            zepc(1) = zem1
            topo1   = topm1
            iv1     = iv

            xepc(2) = xem2
            yepc(2) = yem2
            zepc(2) = zem2
            topo2   = topm2
            iv2     = iv

         else

! This interval is on maximum side of polygon

            xepc(1) = xem2
            yepc(1) = yem2
            zepc(1) = zem2
            topo1   = topm2
            iv1     = iv

            xepc(2) = xem1
            yepc(2) = yem1
            zepc(2) = zem1
            topo2   = topm1
            iv2     = iv

         endif

         exit

      else
         
! This interval touches cone at 1 point

         wt2 = (op%coneang - angm1) / (angm2 - angm1)
         wt1 = 1. - wt2

         xepc(2) = wt1 * xem1  + wt2 * xem2
         yepc(2) = wt1 * yem1  + wt2 * yem2
         zepc(2) = wt1 * zem1  + wt2 * zem2
         topo2   = wt1 * topm1 + wt2 * topm2
         iv2     = iv

      endif
      
   elseif (angm1 > op%coneang .and. angm2 <= op%coneang) then

! This interval touches cone at 1 point

         wt2 = (op%coneang - angm2) / (angm1 - angm2)
         wt1 = 1. - wt2

         xepc(1) = wt1 * xem2  + wt2 * xem1
         yepc(1) = wt1 * yem2  + wt2 * yem1
         zepc(1) = wt1 * zem2  + wt2 * zem1
         topo1   = wt1 * topm2 + wt2 * topm1
         iv1     = iv

   endif

enddo

! Transform horizontal point coordinates

! Components of vector from circle center to plot window center

vecx = xwincent - xcent
vecy = ywincent - ycent
vecz = zwincent - zcent

! Components of vector 90 degrees to the left (in azimuth) from preceding vector

vecleftx = unyec3 * vecz - unzec3 * vecy
veclefty = unzec3 * vecx - unxec3 * vecz
vecleftz = unxec3 * vecy - unyec3 * vecx

! Compute dot product between vector from circle center to plot center
! and vector from circle center to current point

valx = vecx * (xepc(1) - xcent)  &
     + vecy * (yepc(1) - ycent)  &
     + vecz * (zepc(1) - zcent)

! Repeat with 90-left vector

valy = vecleftx * (xepc(1) - xcent)  &
     + veclefty * (yepc(1) - ycent)  &
     + vecleftz * (zepc(1) - zcent)

ange1 = atan2(-valy,valx)  ! Angle increases clockwise

! Repeat dot product for second point

valx = vecx * (xepc(2) - xcent)  &
     + vecy * (yepc(2) - ycent)  &
     + vecz * (zepc(2) - zcent)

valy = vecleftx * (xepc(2) - xcent)  &
     + veclefty * (yepc(2) - ycent)  &
     + vecleftz * (zepc(2) - zcent)

ange2 = atan2(-valy,valx)  ! Angle increases clockwise

! Avoid wrap_around

if (ange2 < ange1) ange2 = ange2 + pi2

! Scale angles to htpn coordinates (in meters along cone circle)

htpn(1) = ange1 * radcone 
htpn(2) = ange2 * radcone
htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine coneplot_w

!===============================================================================

subroutine xyplot_w(iplt,iw,iv1,iv2,topo1,topo2,iok,htpn)

! This subroutine finds points of intersection between vertical plot slab and
! triangular (prism) grid cells

use oplot_coms,  only: op
use mem_grid,    only: xem, yem, topm
use mem_ijtabs,  only: itab_w
use consts_coms, only: pio180
use oname_coms,  only: nl
use misc_coms,  only: io6

implicit none

integer, intent(in)  :: iplt
integer, intent(in)  :: iw
integer, intent(out) :: iv1,iv2
integer, intent(out) :: iok

real, intent(out) :: htpn(4)
real, intent(out) :: topo1,topo2

integer :: im1,im2,im3
integer :: iu1,iu2,iu3
integer :: iuc1,iuc2,iuc3

real :: s1,s2,s3
real :: sc1,sc2,sc3
real :: smin,smax
real :: wt1,wt2
real :: x1,x2
real :: y1,y2
real :: sinvaz,cosvaz

real :: xemc(3),yemc(3)  ! Coords of 3 M points
real :: topmc(3) ! Topo height of 3 M points

im1 = itab_w(iw)%im(1)
im2 = itab_w(iw)%im(2)
im3 = itab_w(iw)%im(3)

iu1 = itab_w(iw)%iu(1)
iu2 = itab_w(iw)%iu(2)
iu3 = itab_w(iw)%iu(3)

sinvaz = sin((90. - op%viewazim) * pio180)
cosvaz = cos((90. - op%viewazim) * pio180)

! Location of 3 M points along line perpendicular to plot slab

s1 = (xem(im1) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
   + (yem(im1) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

s2 = (xem(im2) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
   + (yem(im2) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

s3 = (xem(im3) - nl%plotspecs(iplt)%plotcoord1) * cosvaz  &
   + (yem(im3) - nl%plotspecs(iplt)%plotcoord2) * sinvaz

smin = min(s1,s2,s3)
smax = max(s1,s2,s3)

! Return with iok = 0 if iw column is not intersected by plot slab

iok = 0

if (smin > op%slabloc(iplt) .or. smax < op%slabloc(iplt)) return

iok = 1  ! Since we got here, IW column is intersected by plot cone

! Find IM point that has lowest S coordinate and copy to temporary M point #1
! Fill other 2 points in cyclic order

if (s1 <= s2 .and. s1 <= s3) then  ! m1 has lowest S value

   xemc(1)  = xem(im1)
   yemc(1)  = yem(im1)
   topmc(1) = topm(im1)
   
   xemc(2)  = xem(im2)
   yemc(2)  = yem(im2)
   topmc(2) = topm(im2)

   xemc(3)  = xem(im3)
   yemc(3)  = yem(im3)
   topmc(3) = topm(im3)

   sc1 = s1
   sc2 = s2
   sc3 = s3

   iuc1 = iu1
   iuc2 = iu2
   iuc3 = iu3

elseif (s2 <= s1 .and. s2 <= s3) then  ! m2 has lowest S value

   xemc(1)  = xem(im2)
   yemc(1)  = yem(im2)
   topmc(1) = topm(im2)

   xemc(2)  = xem(im3)
   yemc(2)  = yem(im3)
   topmc(2) = topm(im3)

   xemc(3)  = xem(im1)
   yemc(3)  = yem(im1)
   topmc(3) = topm(im1)

   sc1 = s2
   sc2 = s3
   sc3 = s1

   iuc1 = iu2
   iuc2 = iu3
   iuc3 = iu1

elseif (s3 <= s1 .and. s3 <= s2) then  ! m3 has lowest S value

   xemc(1)  = xem(im3)
   yemc(1)  = yem(im3)
   topmc(1) = topm(im3)

   xemc(2)  = xem(im1)
   yemc(2)  = yem(im1)
   topmc(2) = topm(im1)

   xemc(3)  = xem(im2)
   yemc(3)  = yem(im2)
   topmc(3) = topm(im2)

   sc1 = s3
   sc2 = s1
   sc3 = s2

   iuc1 = iu3
   iuc2 = iu1
   iuc3 = iu2

endif
   
! Find two points of intersection between current IW triangle and X slab

if (sc2 > op%slabloc(iplt)) then

   wt2 = (op%slabloc(iplt) - sc1) / (sc2 - sc1)
   wt1 = 1. - wt2

   x2    = wt1 * xemc(1)  + wt2 * xemc(2)  ! x coord of htpn(2)
   y2    = wt1 * yemc(1)  + wt2 * yemc(2)  ! y coord of htpn(2)
   topo2 = wt1 * topmc(1) + wt2 * topmc(2) ! topo height(2)

   htpn(2) = (x2 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
           - (y2 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
           
   iv2 = iuc3
   
   if (sc3 > op%slabloc(iplt)) then

      wt2 = (op%slabloc(iplt) - sc1) / (sc3 - sc1)
      wt1 = 1. - wt2

      x1    = wt1 * xemc(1)  + wt2 * xemc(3)  ! x coord of htpn(1)
      y1    = wt1 * yemc(1)  + wt2 * yemc(3)  ! y coord of htpn(1)
      topo1 = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

      htpn(1) = (x1 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
              - (y1 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
   
      iv1 = iuc2
      
   else

      wt2 = (op%slabloc(iplt) - sc3) / (sc2 - sc3)
      wt1 = 1. - wt2

      x1    = wt1 * xemc(3)  + wt2 * xemc(2)  ! x coord of htpn(1)
      y1    = wt1 * yemc(3)  + wt2 * yemc(2)  ! y coord of htpn(1)
      topo1 = wt1 * topmc(3) + wt2 * topmc(2) ! topo height(1)

      htpn(1) = (x1 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
              - (y1 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
   
      iv1 = iuc1
      
   endif

else

   wt2 = (op%slabloc(iplt) - sc2) / (sc3 - sc2)
   wt1 = 1. - wt2

   x2    = wt1 * xemc(2)  + wt2 * xemc(3)  ! x coord of htpn(2)
   y2    = wt1 * yemc(2)  + wt2 * yemc(3)  ! y coord of htpn(2)
   topo2 = wt1 * topmc(2) + wt2 * topmc(3) ! topo height(2)

   htpn(2) = (x2 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
           - (y2 - nl%plotspecs(iplt)%plotcoord2) * cosvaz
   
   iv2 = iuc1

   wt2 = (op%slabloc(iplt) - sc1) / (sc3 - sc1)
   wt1 = 1. - wt2

   x1    = wt1 * xemc(1)  + wt2 * xemc(3)  ! x coord of htpn(1)
   y1    = wt1 * yemc(1)  + wt2 * yemc(3)  ! y coord of htpn(1)
   topo1 = wt1 * topmc(1) + wt2 * topmc(3) ! topo height(1)

   htpn(1) = (x1 - nl%plotspecs(iplt)%plotcoord1) * sinvaz  &
           - (y1 - nl%plotspecs(iplt)%plotcoord2) * cosvaz

   iv1 = iuc2

endif

htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine xyplot_w


