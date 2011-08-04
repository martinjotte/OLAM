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
subroutine spawn_nest()

! This subroutine adds nested grid regions at the beginning of a simulation.
! Later will make modified version to add nested grid region(s) during a
!   simulation.

use mem_ijtabs,  only: itab_md, itab_ud, itab_wd, ltab_md, ltab_ud, ltab_wd,  &
                       nest_ud, nest_wd, nloops_m, nloops_u, nloops_w, mrls,  &
                       alloc_itabsd
use mem_grid,    only: nma, nua, nwa, xem, yem, zem,  &
                       alloc_xyzem, nrows, mrows
use misc_coms,   only: io6, ngrids, mdomain, nxp, ngrdll, grdrad, grdlat, grdlon
use consts_coms, only: pio180, erad, pi1, pi2
use max_dims,    only: maxgrds, maxngrdll

implicit none

integer :: iu,iw,iregion,im,iw1,iw2,im1,im2,im3 &
   ,iu1,iu2,iu3,ndiv,iu1o,iu2o,iu3o,iu1o_iw1,iu2o_iw1,iu3o_iw1  &
   ,iu4,iu5,iu6,iw3,ngr,mrlo,mrloo,nwaa,j,npoly,nw
   
integer :: nma0,nua0,nwa0

real :: focdist,xf1,xf2,yf1,yf2,widfac,xw,yw,xm1,ym1,xm2,ym2,xm3,ym3
real :: expansion

integer :: nside    ! Number of sides of nested grid polygon
integer :: nper(16) ! Value of perimeter side counter at end of each side
integer, allocatable :: imper(:) ! Ouside IW index at each perimeter index
integer, allocatable :: iuper(:) ! Boundary IU index at each perimeter index
integer :: kma,kua,kwa   ! New M,U,W indices while constructing nest perimeter
integer :: ngrp   ! Number of perimeter groups
integer :: npts = 10000  ! Sufficient array space for # of CM pts along FM perimeter

integer, allocatable :: jm(:,:),ju(:,:),iurow_pent(:),igsize(:),nwdivg(:)

real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)

integer :: nper2 = 0 ! Actual # of perimeter pts
integer :: nconcave
integer :: nccv, jj, ipt, jpt, inflag

integer, allocatable :: npolyper(:) ! npoly at each perimeter M pt
integer, allocatable :: nwdivper(:) ! # divided W pts at at each perimeter M pt
integer, allocatable :: nearpent(:) ! flag = 1 if adjacent outside M pt is poly5

!-------------------------------------------------------------------------------
! New option introduced April 2010:  Define nested grid region as comprising the 
! W points adjacent to all M points that are within a specified distance of one
! or more connected line segments.

real, external :: linesegdist

real :: seglat, seglon ! lat/lon of segment midpoint
real :: xs(2),ys(2)    ! PS coordinates of segment endpoints
!-------------------------------------------------------------------------------

! Make duplicate of current grid dimensions

nma0 = nma
nua0 = nua
nwa0 = nwa 

do ngr = 2,ngrids  ! Loop over nested grids

!-------------------------------------------------------------------------------
! Set default number of coarse-grid transition rows (MROWS) across which the
! transition from coarse-to-fine resolution is made (allowed values are 1:5).

! Set number of transition rows to be used; may put this in OLAMIN eventually
! Set default value of NCONCAVE (allowed values are 1:2).
! These could be made grid-dependent if desired.

   mrows = 3
   nconcave = 1

! NCONCAVE = 1 eliminates only sharp concavities and activates new
!    refined mesh algorithm 

! NCONCAVE = 2 eliminates both sharp and gentler concavities and uses
!    older refined mesh algorithm.  MROWS /= 2 may only be used with
!    this older algorithm.

! If NCONCAVE is set to 1, MROWS will (must) be set to 2

   if (nconcave == 1) mrows = 2
!-------------------------------------------------------------------------------

! Allocate temporary tables

   allocate (nest_ud(nua), nest_wd(nwa))  ! Nest relations
   allocate (xem_temp(nma), yem_temp(nma), zem_temp(nma))
   allocate (iurow_pent(nua))

   allocate (jm(nrows+1,npts), ju(nrows,npts))
   allocate (imper(npts), iuper(npts), igsize(npts), nearpent(npts), nwdivg(npts))

   allocate (npolyper(npts), nwdivper(npts))

   call pent_urows(iurow_pent)

! Copy ITAB information and M-point coordinates to temporary tables

   call move_alloc(itab_md, ltab_md)
   call move_alloc(itab_ud, ltab_ud)
   call move_alloc(itab_wd, ltab_wd)

   do im = 1,nma
      xem_temp(im) = xem(im)
      yem_temp(im) = yem(im)
      zem_temp(im) = zem(im)
   enddo

! Deallocate main tables

   deallocate (xem, yem, zem)      

! Locate and flag all W triangles to be subdivided
! Define 3 new W indices and 3 new U indices for each W - attach to W.

   do im = 2,nma

! Check whether M location is within elliptical current nested region

      inflag = 0

      do ipt = 1,ngrdll(ngr)
         jpt = min(ipt+1,ngrdll(ngr)) ! jpt used in case one wants to define
                                          ! only a single endpoint

! Transform segment endpoints to PS space using mean lat/lon of each segment

! (If multiple segments are used, none should be excessively long in order
! to avoid large PS transformation discontinuities at segment endpoints.)

         seglat = .5 * (grdlat(ngr,ipt) + grdlat(ngr,jpt))
         seglon = .5 * (grdlon(ngr,ipt) + grdlon(ngr,jpt))

         call ll_xy (grdlat(ngr,ipt),grdlon(ngr,ipt), &
            seglat,seglon,xs(1),ys(1))

         call ll_xy (grdlat(ngr,jpt),grdlon(ngr,jpt), &
            seglat,seglon,xs(2),ys(2))

! Transform current M point to PS space using mean lat/lon of each segment

         call e_ps(xem_temp(im),yem_temp(im),zem_temp(im), &
            seglat,seglon,xm1,ym1)

! If current M point is close enough to line segment, flag it for inclusion
! in refined grid interior

         if (linesegdist(xm1,ym1,xs(1),ys(1),xs(2),ys(2)) &
            < grdrad(ngr)) inflag = 1

      enddo
      
      if (inflag == 1) then

! W cells surrounding M are to be subdivided

         npoly = ltab_md(im)%npoly
         
         do j = 1,npoly
         
            iw = ltab_md(im)%iw(j)
            
            if (nest_wd(iw)%iw(1) == 0) then

!        write(io6,*) 'subdividing W cell ',iw

               nest_wd(iw)%iu(1) = nua0 + 1
               nest_wd(iw)%iu(2) = nua0 + 2
               nest_wd(iw)%iu(3) = nua0 + 3

               nest_wd(iw)%iw(1) = nwa0 + 1
               nest_wd(iw)%iw(2) = nwa0 + 2
               nest_wd(iw)%iw(3) = nwa0 + 3

               nua0 = nua0 + 3
               nwa0 = nwa0 + 3

            endif
      
         enddo

      endif
      
   enddo

! Add W points to nested grid region in order to eliminate concavities
! (or to eliminate sharp concavities).  This requires iterative procedure

   nwaa = 0  ! Counter of already-existing W points - initialize to zero to
             ! force at least one pass through the following DO WHILE loop

   do while (nwa0 > nwaa)

      nwaa = nwa0

      do im = 2,nma

         npoly = ltab_md(im)%npoly
         nw = 0  ! Initialize counter for subdivided W points around this M point

         do j = 1,npoly
            iw = ltab_md(im)%iw(j)
            
            if (nest_wd(iw)%iw(1) > 0) then
               nw = nw + 1  ! Count up subdivided W points around this M point
            endif
         enddo

! Check npoly and nw for illegal values

         if (nw > 0 .and. nw < npoly .and. &
              (npoly == 5 .or. npoly - nw <= nconcave)) then

! This M point is a concavity so activate remaining unactivated W points around it

            do j = 1,npoly
               iw = ltab_md(im)%iw(j)

               if (nest_wd(iw)%iw(1) == 0) then
                  nest_wd(iw)%iu(1) = nua0 + 1
                  nest_wd(iw)%iu(2) = nua0 + 2
                  nest_wd(iw)%iu(3) = nua0 + 3

                  nest_wd(iw)%iw(1) = nwa0 + 1
                  nest_wd(iw)%iw(2) = nwa0 + 2
                  nest_wd(iw)%iw(3) = nwa0 + 3

                  nua0 = nua0 + 3
                  nwa0 = nwa0 + 3
                  
                  write(io6,*) 'Activating W point ',iw,' to prevent concavity'
               
               endif
            enddo

         endif

      enddo

      if (nconcave == 1) then

         call perim_map2(npts,nper2,imper,iuper,npolyper,nwdivper,nearpent)

! Add points for refinement if 2 consecutive concave points occur

         nccv = 0
         do j = 1,nper2
            if (nwdivper(j) > 3) then
               nccv = nccv + 1
            else
               nccv = 0
            endif
         
            if (nccv > 1) then

               do jj = j-2,j
                  iw1 = ltab_ud(iuper(jj))%iw(1)
                  iw2 = ltab_ud(iuper(jj))%iw(2)
      
                  if (nest_wd(iw1)%iw(3) == 0) then
                     nest_wd(iw1)%iu(1) = nua0 + 1
                     nest_wd(iw1)%iu(2) = nua0 + 2
                     nest_wd(iw1)%iu(3) = nua0 + 3

                     nest_wd(iw1)%iw(1) = nwa0 + 1
                     nest_wd(iw1)%iw(2) = nwa0 + 2
                     nest_wd(iw1)%iw(3) = nwa0 + 3
                  else
                     nest_wd(iw2)%iu(1) = nua0 + 1
                     nest_wd(iw2)%iu(2) = nua0 + 2
                     nest_wd(iw2)%iu(3) = nua0 + 3

                     nest_wd(iw2)%iw(1) = nwa0 + 1
                     nest_wd(iw2)%iw(2) = nwa0 + 2
                     nest_wd(iw2)%iw(3) = nwa0 + 3
                  endif

                  nua0 = nua0 + 3
                  nwa0 = nwa0 + 3

               enddo

            endif

         enddo

      endif

   enddo

! Print perimeter map

   print*, ' '
   do j = 1,nper2
      write(6,'(a,6i7)') 'perim ',j,imper(j),iuper(j),npolyper(j), &
                                  nwdivper(j),nearpent(j)
   enddo
   print*, ' '

! Nested region should be fully expanded now without any concavities, 
! 5-edge vertices, or consecutive weak concavities.

! Map out perimeter of new mesh refined region (returning information that
! depends on nconcave value)

   if (nconcave == 2) then
      call perim_map(nper,nside,npts,imper,iuper)
   endif

! Define new vertex index for midpoint of each original U edge that is adjacent
! to an original triangle that is being subdivided.  Attach new vertex
! index to old U edge.  Also, define new U index for second half of U.
! Attach to U.  [Make adjacent to U%m2.]

   do iu = 2,nua

      ! Check whether this U is adjacent to a W that is being subdivided

      iw1 = ltab_ud(iu)%iw(1)
      iw2 = ltab_ud(iu)%iw(2)

      if (nest_wd(iw1)%iw(3) > 0 .or. nest_wd(iw2)%iw(3) > 0) then

         nest_ud(iu)%im = nma0 + 1
         nest_ud(iu)%iu = nua0 + 1

         nma0 = nma0 + 1
         nua0 = nua0 + 1

      endif
   enddo

! Save current values of nma0, nua0, nwa0 prior to adding boundary points

   kma = nma0
   kua = nua0
   kwa = nwa0

! Form groups of original M and U indices on perimeter and increase nma0, nua0, 
! and nwa0 according to what groups will require

   if (nconcave == 2) then

      call perim_add(jm,ju,npts,nside,nma0,nua0,nwa0, &
                     nper,ngrp,igsize,imper,iuper,iurow_pent)

   else
   
      call perim_add2(npts,nper2,imper,iuper,nwdivper,nearpent,iurow_pent, &
           nma0,nua0,nwa0,ngrp,igsize,jm,ju,nwdivg)

   endif

! Allocate main tables to expanded size
! Initialize all neighbor indices to zero

   call alloc_itabsd(nma0,nua0,nwa0)

   call alloc_xyzem(nma0)

! Memory copy to main tables

   do im = 1,nma
      itab_md(im)%loop(1:nloops_m) = ltab_md(im)%loop(1:nloops_m)
      itab_md(im)%itopm = ltab_md(im)%itopm
      xem(im) = xem_temp(im)
      yem(im) = yem_temp(im)
      zem(im) = zem_temp(im)
   enddo

   itab_ud(1:nua) = ltab_ud(1:nua)
   itab_wd(1:nwa) = ltab_wd(1:nwa)
  
! Average coordinates to new M points

   do iu = 2,nua
      if (nest_ud(iu)%im > 0) then
         im  = nest_ud(iu)%im
         im1 = itab_ud(iu)%im(1)
         im2 = itab_ud(iu)%im(2)

         xem(im) = .5 * (xem(im1) + xem(im2))
         yem(im) = .5 * (yem(im1) + yem(im2))
         zem(im) = .5 * (zem(im1) + zem(im2))
      endif
   enddo

! Contruct tables for new fully subdivided triangles

   mrloo = 0  ! Initialize check variable for uniform mrlw over current nested grid

   do iw = 2,nwa

! Check if IW is fully subdivided cell

      if (nest_wd(iw)%iw(3) > 0) then

! This is fully subdivided W cell

! Mapping of original undivided triangles

         iu1o = ltab_wd(iw)%iu(1)
         iu2o = ltab_wd(iw)%iu(2)
         iu3o = ltab_wd(iw)%iu(3)
         mrlo = ltab_wd(iw)%mrlw

! Check of mrlw value for current nested grid

         if (mrloo == 0) mrloo = mrlo  ! Set to first nonzero mrlw encountered
         if (mrlo /= mrloo) then
            write(io6,*) 'Current nested grid ',ngr
            write(io6,*) 'crosses pre-existing grid boundary.'
            write(io6,*) 'iw = ',iw
            write(io6,*) 'stopping model'
            stop 'stop - nested grid out of bounds'
         endif

         iu1o_iw1 = ltab_ud(iu1o)%iw(1)
         iu2o_iw1 = ltab_ud(iu2o)%iw(1)
         iu3o_iw1 = ltab_ud(iu3o)%iw(1)

! Mapping of new divided triangles

         iu1 = nest_wd(iw)%iu(1)
         iu2 = nest_wd(iw)%iu(2)
         iu3 = nest_wd(iw)%iu(3)

         itab_ud(iu1)%iup = iu1
         itab_ud(iu2)%iup = iu2
         itab_ud(iu3)%iup = iu3
         
         call udloops('f',iu1, 1, 4, 7, 8,11,12,13,14,16,20)
         call udloops('n',iu1,21,22,23, 0, 0, 0, 0, 0, 0, 0)

         call udloops('f',iu2, 1, 4, 7, 8,11,12,13,14,16,20)
         call udloops('n',iu2,21,22,23, 0, 0, 0, 0, 0, 0, 0)

         call udloops('f',iu3, 1, 4, 7, 8,11,12,13,14,16,20)
         call udloops('n',iu3,21,22,23, 0, 0, 0, 0, 0, 0, 0)

         iu4 = nest_ud(iu1o)%iu       
         iu5 = nest_ud(iu2o)%iu       
         iu6 = nest_ud(iu3o)%iu

         itab_ud(iu4)%loop(1:nloops_u) = itab_ud(iu1o)%loop(1:nloops_u)
         itab_ud(iu5)%loop(1:nloops_u) = itab_ud(iu2o)%loop(1:nloops_u)
         itab_ud(iu6)%loop(1:nloops_u) = itab_ud(iu3o)%loop(1:nloops_u)

         itab_ud(iu4)%iup = iu4
         itab_ud(iu5)%iup = iu5
         itab_ud(iu6)%iup = iu6

         iw1 = nest_wd(iw)%iw(1)        
         iw2 = nest_wd(iw)%iw(2)        
         iw3 = nest_wd(iw)%iw(3) 

! Fill tables with new values

         itab_wd(iw)%iu(1) = iu1       
         itab_wd(iw)%iu(2) = iu2       
         itab_wd(iw)%iu(3) = iu3   

         itab_wd(iw1)%iu(1) = iu1    
         itab_wd(iw2)%iu(1) = iu2    
         itab_wd(iw3)%iu(1) = iu3

         itab_wd(iw)%mrlw  = mrlo + 1
         itab_wd(iw1)%mrlw = mrlo + 1
         itab_wd(iw2)%mrlw = mrlo + 1
         itab_wd(iw3)%mrlw = mrlo + 1

         itab_wd(iw1)%mrlw_orig = mrlo + 1
         itab_wd(iw2)%mrlw_orig = mrlo + 1
         itab_wd(iw3)%mrlw_orig = mrlo + 1

         itab_ud(iu1o)%im(2) = nest_ud(iu1o)%im
         itab_ud(iu2o)%im(2) = nest_ud(iu2o)%im
         itab_ud(iu3o)%im(2) = nest_ud(iu3o)%im

         itab_ud(iu4)%im(1) = nest_ud(iu1o)%im
         itab_ud(iu4)%im(2) = ltab_ud(iu1o)%im(2)
         itab_ud(iu5)%im(1) = nest_ud(iu2o)%im
         itab_ud(iu5)%im(2) = ltab_ud(iu2o)%im(2)
         itab_ud(iu6)%im(1) = nest_ud(iu3o)%im
         itab_ud(iu6)%im(2) = ltab_ud(iu3o)%im(2)

         if (iw == iu1o_iw1) then
            itab_wd(iw3)%iu(2) = iu1o
            itab_wd(iw2)%iu(3) = iu4

            itab_ud(iu1)%im(1) = nest_ud(iu2o)%im            
            itab_ud(iu1)%im(2) = nest_ud(iu3o)%im            
            itab_ud(iu1)%iw(1) = iw1
            itab_ud(iu1)%iw(2) = iw

            itab_ud(iu1o)%iw(1) = iw3
            itab_ud(iu4)%iw(1) = iw2
         else
            itab_wd(iw3)%iu(2) = iu4
            itab_wd(iw2)%iu(3) = iu1o

            itab_ud(iu1)%im(1) = nest_ud(iu3o)%im            
            itab_ud(iu1)%im(2) = nest_ud(iu2o)%im            
            itab_ud(iu1)%iw(1) = iw
            itab_ud(iu1)%iw(2) = iw1

            itab_ud(iu1o)%iw(2) = iw2
            itab_ud(iu4)%iw(2) = iw3
         endif    

         if (iw == iu2o_iw1) then
            itab_wd(iw1)%iu(2) = iu2o
            itab_wd(iw3)%iu(3) = iu5

            itab_ud(iu2)%im(1) = nest_ud(iu3o)%im            
            itab_ud(iu2)%im(2) = nest_ud(iu1o)%im            
            itab_ud(iu2)%iw(1) = iw2
            itab_ud(iu2)%iw(2) = iw

            itab_ud(iu2o)%iw(1) = iw1
            itab_ud(iu5)%iw(1) = iw3
         else
            itab_wd(iw1)%iu(2) = iu5
            itab_wd(iw3)%iu(3) = iu2o

            itab_ud(iu2)%im(1) = nest_ud(iu1o)%im            
            itab_ud(iu2)%im(2) = nest_ud(iu3o)%im            
            itab_ud(iu2)%iw(1) = iw
            itab_ud(iu2)%iw(2) = iw2
            
            itab_ud(iu2o)%iw(2) = iw3
            itab_ud(iu5)%iw(2) = iw1
         endif    

         if (iw == iu3o_iw1) then
            itab_wd(iw2)%iu(2) = iu3o
            itab_wd(iw1)%iu(3) = iu6

            itab_ud(iu3)%im(1) = nest_ud(iu1o)%im            
            itab_ud(iu3)%im(2) = nest_ud(iu2o)%im            
            itab_ud(iu3)%iw(1) = iw3
            itab_ud(iu3)%iw(2) = iw

            itab_ud(iu3o)%iw(1) = iw2
            itab_ud(iu6)%iw(1) = iw1
         else
            itab_wd(iw2)%iu(2) = iu6
            itab_wd(iw1)%iu(3) = iu3o

            itab_ud(iu3)%im(1) = nest_ud(iu2o)%im            
            itab_ud(iu3)%im(2) = nest_ud(iu1o)%im            
            itab_ud(iu3)%iw(1) = iw
            itab_ud(iu3)%iw(2) = iw3

            itab_ud(iu3o)%iw(2) = iw1
            itab_ud(iu6)%iw(2) = iw2
         endif    

      endif

   enddo    ! end of iw loop

! Fill transition zone

   if (nconcave == 2) then
      call perim_fill(ngr,mrloo,kma,kua,kwa,jm,ju,npts,ngrp,igsize)
   else
      call perim_fill2(ngr,mrloo,kma,kua,kwa,jm,ju,npts,ngrp,igsize,nwdivg)
   endif

! Fill itabs loop tables for newly spawned points (U pts already done above)

   do im = nma+1,nma0
      itab_md(im)%itopm = im
      call mdloops('f',im,1,0,1,0)
   enddo

   do iw = nwa+1,nwa0
      itab_wd(iw)%iwp = iw
      call wdloops('f',iw, 1, 3, 5, 6, 7, 8,11,12,13,14)
      call wdloops('n',iw,15,16,17,18,19,20,21,23,25,26)
      call wdloops('n',iw,27,28,29,30,33,34, 0, 0 ,0 ,0)
   enddo

! Copy new counter values

   nma = nma0
   nua = nua0
   nwa = nwa0

   deallocate (ltab_md,ltab_ud,ltab_wd)
   deallocate (nest_ud,nest_wd)
   deallocate (xem_temp,yem_temp,zem_temp)

   deallocate (jm,ju)
   deallocate (imper,iuper,iurow_pent,igsize,nearpent,nwdivg)
   deallocate (npolyper,nwdivper)

   call tri_neighbors()

!! Call subroutine to ID W cells just outside and just inside current NGR
! border.  This is permanent ID, used in spring dynamics even when new
! grids are added.

   call perim_mrow()

! This is the place to do spring dynamics

   call spring_dynamics(ngr)

   write(io6,'(/,a,i2)') 'Finished spawning grid number ',ngr
   write(io6,'(a,i8)')   ' nma = ',nma
   write(io6,'(a,i8)')   ' nua = ',nua
   write(io6,'(a,i8)')   ' nwa = ',nwa

enddo   ! end of ngr loop

! Set mrls equal to maximum mrlw value

do iw = 2,nwa
   if (mrls < itab_wd(iw)%mrlw) mrls = itab_wd(iw)%mrlw
enddo

! Push all M point coordinates out to earth radius if not already there

do im = 2,nma
   expansion = erad / sqrt(xem(im) ** 2  &
                         + yem(im) ** 2  &
                         + zem(im) ** 2  )

   xem(im) = xem(im) * expansion
   yem(im) = yem(im) * expansion
   zem(im) = zem(im) * expansion
enddo

return
end subroutine spawn_nest

!==============================================================================

subroutine pent_urows(iurow_pent)

use mem_ijtabs, only: itab_md,itab_ud,itab_wd
use mem_grid,   only: nua, nwa, impent, mrows
use misc_coms,  only: io6

implicit none

integer, intent(out) :: iurow_pent(nua)

integer :: ipent,im,j,iw,irow,jrow,iw1,iw2,iw3,iwrow,iu

! automatic arrays

integer :: iwrow_temp1(nwa)
integer :: iwrow_temp2(nwa)

! Initialize temporary arrays to zero before loop over pentagon points

iurow_pent(1:nua) = 0
iwrow_temp1(1:nwa) = 0
iwrow_temp2(1:nwa) = 0

! Loop over all 12 pentagon M points

do ipent = 1,12

   im = impent(ipent)
! Loop over W points that surround current M point; set row flag to 1

   do j = 1,5
      iw = itab_md(im)%iw(j)
      iwrow_temp1(iw) = 1
      iwrow_temp2(iw) = 1
   enddo
   
enddo
   
! Advance outward and flag each row

do irow = 1,2*mrows-1
   jrow = mod(irow,2)

   do iw = 2,nwa

      if (iwrow_temp1(iw) == 0) then

! If IW is adjacent to any other IW cell with nonzero mrow, 
! set mrow for IW cell 

         iw1 = itab_wd(iw)%iw(1)
         iw2 = itab_wd(iw)%iw(2)
         iw3 = itab_wd(iw)%iw(3)

! Check for positive mrow values

         iwrow = max(iwrow_temp1(iw1)  &
                    ,iwrow_temp1(iw2)  &
                    ,iwrow_temp1(iw3))

         if (iwrow > 0) iwrow_temp2(iw) = iwrow + jrow

      endif

   enddo

   do iw = 2,nwa
      iwrow_temp1(iw) = iwrow_temp2(iw)
   enddo

enddo

! Loop over all U points and flag those between unequal nonzero iwrow values

do iu = 2,nua
   iw1 = itab_ud(iu)%iw(1)
   iw2 = itab_ud(iu)%iw(2)

   if (iwrow_temp1(iw1) > 0 .and. iwrow_temp1(iw2) > 0 .and.  &
       iwrow_temp1(iw1) /= iwrow_temp1(iw2)) then
          
       iurow_pent(iu) = min(iwrow_temp1(iw1),iwrow_temp1(iw2))
   endif
enddo

return
end subroutine pent_urows

!==============================================================================

subroutine perim_map(nper,nside,npts,imper,iuper)

! Perim_map maps the perimeter points of a nested grid that is being spawned

use mem_grid,  only: nma
use misc_coms, only: io6

implicit none

integer, intent(out) :: nper(16) ! Value of perimeter side counter at end of each side
integer, intent(out) :: nside    ! Number of sides of nested grid polygon
integer, intent(in)  :: npts
integer, intent(out) :: imper(npts) ! Boundary IM index at each perimeter index
integer, intent(out) :: iuper(npts) ! Boundary IU index at each perimeter index

integer :: imstart  ! IM index of starting M point (at a corner of ngr perimeter)
integer :: im       ! dummy im index
integer :: nwdiv    ! number of subdivided W pts adjacent to current M pt
integer :: ima      ! Current M pt in counterclockwise path around perimeter
integer :: imb      ! Next M pt in counterclockwise path around perimeter
integer :: iua      ! Current U pt in counterclockwise path around perimeter
integer :: iper     ! Current perimeter point counter
integer :: iside    ! Current perimeter side counter

! Set IM starting point to zero

imstart = 0

! Loop over all ORIGINAL M points and find the first one that has exactly 2 
! fully-divided W neighbors in the grid being spawned.

do im = 2,nma   
   call count_wdiv(im,nwdiv)
   if (nwdiv == 2) then
      imstart = im
      exit
   endif
enddo

! Bug check:  make sure that imstart is not zero now.

if (imstart == 0) then
   write(io6,*) 'imstart is zero - stopping model'
   stop 'stop imstart'
endif

! March around NGR boundary in a counterclockwise direction beginning at IMSTART

ima = imstart
imper(1) = ima
iper = 0
iside = 1

nper(:) = 0

do  ! Loop over all original M points on ngr perimeter 

! Find next M and U points (imb, iua) on perimeter and outside W point (iwout)

   call perim_ngr(ima, imb, iua)

! Increment perimeter point counter and store point indices

   nper(iside) = nper(iside) + 1
   iper = iper + 1

   iuper(iper) = iua
   imper(iper+1) = imb

! Check if imb is a corner point.

   call count_wdiv(imb,nwdiv)

   if (nwdiv == 2) then

! imb is corner point.  Store current iper index in nper array.

      nper(iside) = iper

! Check if imb equals istart.  If it does, exit loop

      if (imb == imstart) then

! imb equals istart.  Store total number of sides in nside.  Exit Do loop
      
         nside = iside
         exit
      endif
   
! imb does not equal istart.  Advance to next side

      iside = iside + 1
         
   endif

   ima = imb

enddo

return
end subroutine perim_map

!==============================================================================

subroutine perim_map2(npts,nper2,imper,iuper,npolyper,nwdivper,nearpent)

! Perim_map maps the perimeter points of a nested grid that is being spawned

use mem_grid,  only: nma
use misc_coms, only: io6
use mem_ijtabs, only: ltab_md, nest_wd, ltab_ud

implicit none

integer, intent(in)  :: npts  ! Estimated # of perimeter pts (for array dims only)
integer, intent(out) :: nper2 ! Actual # of perimeter pts

integer, intent(out) :: imper(npts) ! Boundary IM index at each perimeter index
integer, intent(out) :: iuper(npts) ! Boundary IU index at each perimeter index
integer, intent(out) :: npolyper(npts) ! npoly at each perimeter M pt
integer, intent(out) :: nwdivper(npts) ! # divided W pts at at each perimeter M pt
integer, intent(out) :: nearpent(npts) ! flag = 1 if adjacent outside M pt is poly5

integer :: imstart  ! IM index of starting M point (at a corner of ngr perimeter)
integer :: im       ! dummy im index
integer :: nwdiv    ! number of subdivided W pts adjacent to current M pt
integer :: ima      ! Current M pt in counterclockwise path around perimeter
integer :: imb      ! Next M pt in counterclockwise path around perimeter
integer :: iua      ! Current U pt in counterclockwise path around perimeter

integer :: j, iw, npoly, iu, im1, im2, iw1, iw2

! Set IM starting point to zero

imstart = 0

! Loop over all ORIGINAL M points and find the first one that is on boundary
! of refined area being spawned

do im = 2,nma   

   npoly = ltab_md(im)%npoly
   nwdiv = 0

! Loop over all W points that are adjacent to this M point

   do j = 1,npoly
      iw = ltab_md(im)%iw(j)

      if (nest_wd(iw)%iw(3) > 0) then
         nwdiv = nwdiv + 1
      endif
   enddo

   if (nwdiv > 0 .and. nwdiv < npoly) then
      imstart = im
      exit
   endif

enddo

! Bug check: make sure that imstart is not zero now.

if (imstart == 0) then
   write(io6,*) 'imstart is zero - stopping model'
   stop 'stop imstart'
endif

! March around NGR boundary in a counterclockwise direction beginning at IMSTART

ima = imstart
nper2 = 1

do  ! Loop over all original M points on ngr perimeter 

! Find next M and U points (imb, iua) on perimeter and outside W point (iwout)

   npoly = ltab_md(ima)%npoly
   nwdiv = 0
   nearpent(nper2) = 0

! Loop over all W points that are adjacent to 'ma' point

   do j = 1,npoly
      iw = ltab_md(ima)%iw(j)

      if (nest_wd(iw)%iw(3) > 0) then
         nwdiv = nwdiv + 1
      endif

      iu = ltab_md(ima)%iu(j)

      im1 = ltab_ud(iu)%im(1)
      im2 = ltab_ud(iu)%im(2)

      iw1 = ltab_ud(iu)%iw(1)
      iw2 = ltab_ud(iu)%iw(2)

      if (nest_wd(iw1)%iw(3) == 0 .and. nest_wd(iw2)%iw(3) == 0) then

         if (ima == im1 .and. ltab_md(im2)%npoly == 5) nearpent(nper2) = 1
         if (ima == im2 .and. ltab_md(im1)%npoly == 5) nearpent(nper2) = 1

      endif

   enddo

   call perim_ngr(ima, imb, iua)

   imper   (nper2) = ima
   iuper   (nper2) = iua
   npolyper(nper2) = npoly
   nwdivper(nper2) = nwdiv

! Check if imb equals istart.  If it does, exit loop

   if (imb == imstart) exit

   nper2 = nper2 + 1
   ima = imb

enddo

return
end subroutine perim_map2

!==============================================================================

subroutine count_wdiv(im,nwdiv)

! Subroutine count_wdiv is to be used during the process of spawning a nested grid,
! after temporary arrays nest_ud and nest_wd have been filled for interior points.

! Given any M point IM, determine how many fully-divided W points are adjacent to it.

use mem_ijtabs, only: ltab_md, nest_wd
use misc_coms,  only: io6

implicit none

integer, intent(in) :: im
integer, intent(out) :: nwdiv

integer :: j, iw

nwdiv = 0

! Loop over all W points that are adjacent to this M point

do j = 1,ltab_md(im)%npoly
   iw = ltab_md(im)%iw(j)

   if (nest_wd(iw)%iw(3) > 0) then
      nwdiv = nwdiv + 1
   endif
enddo

return
end subroutine count_wdiv

!==============================================================================

subroutine perim_ngr(imstart, imnext, iunext)

! Subroutine perim_ngr is to be used during the process of spawning a nested grid,
! after temporary arrays nest_ud and nest_wd have been filled for interior points.

! Given any M point on the nested grid boundary that existed before the spawn 
! process began, and proceeding along the boundary in a counterclockwise direction,
! find the adjacent U and M points that existed before the spawn process began.

use mem_ijtabs, only: ltab_md, ltab_ud, ltab_wd, nest_wd
use misc_coms,  only: io6

implicit none

integer, intent(in) :: imstart  ! starting original M point

integer, intent(out) :: imnext  ! next original M pt (counterclockwise)
integer, intent(out) :: iunext  ! next original U pt (counterclockwise)

integer :: j, im1, im2, iw1, iw2, iu

imnext = 0
iunext = 0

! Loop over all U points connected to current M point

do j = 1,ltab_md(imstart)%npoly
   iu = ltab_md(imstart)%iu(j)
   
   im1 = ltab_ud(iu)%im(1)
   im2 = ltab_ud(iu)%im(2)

   iw1 = ltab_ud(iu)%iw(1)
   iw2 = ltab_ud(iu)%iw(2)

! Current IU point is either along or inside NGR boundary

   if (im1 == imstart .and. nest_wd(iw1)%iw(3) == 0 .and. &
                            nest_wd(iw2)%iw(3) > 0) then

! Current IU point is on boundary and in clockwise direction,
! and next M point is im2.          

      iunext = iu
      imnext = im2

      exit

   elseif (im2 == imstart .and. nest_wd(iw2)%iw(3) == 0 .and. &
                                nest_wd(iw1)%iw(3) > 0) then

! Current IU point is on boundary and in clockwise direction,
! and next M point is im1.          

      iunext = iu
      imnext = im1

      exit

   endif

enddo

if (iunext < 2) then
   write(io6,*) 'iunext < 2 in subroutine perim_ngr - stopping model'
   stop 'perim_ngr1'
elseif (imnext < 2) then
   write(io6,*) 'imnext < 2 in subroutine perim_ngr - stopping model'
   stop 'perim_ngr2'
endif

return
end subroutine perim_ngr

!==============================================================================

subroutine perim_add(jm,ju,npts,nside,nma0,nua0,nwa0  &
                    ,nper,ngrp,igsize,imper,iuper,iurow_pent)

use mem_grid,   only: nrows, mrows, nua
use misc_coms,  only: io6

implicit none

integer, intent(inout) :: npts,nside,nma0,nua0,nwa0
integer, intent(out) :: jm(nrows+1,npts)
integer, intent(out) :: ju(nrows,npts)
integer, intent(out) :: ngrp
integer, intent(out) :: igsize(npts)
integer, intent(in)  :: imper(npts)
integer, intent(in)  :: iuper(npts)
integer, intent(inout) :: nper(16)
integer, intent(in) :: iurow_pent(nua)

integer :: iper,ig,iside,ileft,iwid
integer :: igs   ! counter for group number on current side

integer :: iu,jside

! Loop over all original U points on perimeter of new grid

iper = 1
igs = 0  ! group size
iside = 1

do while (iside <= nside)

   do while (iper <= nper(iside))

      iu = iuper(iper)

! Check if current U point is influenced by pentagon M point.
! If so, increase group size.

      if (iurow_pent(iu) > 0) then
         igs = igs + 1
         
! Check if group size equals iurow_pent(iu)

         if (igs == iurow_pent(iu)) then         
         
! This group of points is influenced by pentagon M point.  Place them
! Into a new "side" of the nest.  Check if there are any remaining points
! on current ORIGINAL side.  Shift remaining iper points to new sides
! (a shift of 1 or 2 sides).
            
            if (iper == nper(iside)) then
               do jside = nside,iside,-1
                  nper(jside+1) = nper(jside)
               enddo
               nside = nside + 1
            else
               do jside = nside,iside,-1
                  nper(jside+2) = nper(jside)
               enddo
               nside = nside + 2
               nper(iside + 1) = iper
            endif
            
            nper(iside) = iper - igs
            igs = 0

         endif 
         
      endif

      iper = iper + 1

   enddo
   
   iside = iside + 1
   
enddo

! Determine width of transition zone at each point on perimeter

iper = 1
ig = 0  ! group number

do iside = 1,nside

   igs = 0

   do while (iper <= nper(iside))

      ig = ig + 1
      igs = igs + 1

      ileft = nper(iside) - iper + 1

! FIRST METHOD: START WITH MAX AND DECREASE AS END OF SIDE IS APPROACHED

      if     (ileft == 16 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 12 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 11 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 9 .and. mrows == 4) then
         iwid = 3
      elseif (ileft == 8 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 7 .and. mrows == 5) then
         iwid = 4
      elseif (ileft == 6 .and. mrows >= 4) then
         iwid = 3
      elseif (ileft == 5 .and. mrows == 4) then
         iwid = 3
      elseif (ileft == 4 .and. mrows == 3) then
         iwid = 2
      else
         iwid = min(ileft,mrows)
      endif

! SECOND METHOD: USE 2 AT BOTH ENDS

!      if (igs == 1) then
!         iwid = min(ileft,2)
!      else
!         iwid = min(mrows,max(2,ileft-2),ileft)
!      endif

      igsize(ig) = iwid

                    jm(1,ig) = imper(iper)
                    jm(2,ig) = imper(iper+1)
      if (iwid > 1) jm(3,ig) = imper(iper+2)
      if (iwid > 2) jm(4,ig) = imper(iper+3)
      if (iwid > 3) jm(5,ig) = imper(iper+4)
      if (iwid > 4) jm(6,ig) = imper(iper+5)

                    ju(1,ig) = iuper(iper)
      if (iwid > 1) ju(2,ig) = iuper(iper+1)
      if (iwid > 2) ju(3,ig) = iuper(iper+2)       
      if (iwid > 3) ju(4,ig) = iuper(iper+3)      
      if (iwid > 4) ju(5,ig) = iuper(iper+4)

      iper = iper + iwid

! Add required number of W,U,M points for this group

      nwa0 = nwa0 + iwid**2
      nua0 = nua0 + iwid * (3 * iwid - 1) / 2
      nma0 = nma0 + iwid * (iwid - 1) / 2

   enddo

enddo

ngrp = ig

return
end subroutine perim_add

!==============================================================================

subroutine perim_add2(npts,nper2,imper,iuper,nwdivper,nearpent,iurow_pent, &
     nma0,nua0,nwa0,ngrp,igsize,jm,ju,nwdivg)

use mem_grid,   only: nrows, mrows, nua
use misc_coms,  only: io6

implicit none

integer, intent(in) :: npts,nper2
integer, intent(in) :: imper(npts)
integer, intent(in) :: iuper(npts)
integer, intent(in) :: nwdivper(npts)
integer, intent(in) :: nearpent(npts)
integer, intent(in) :: iurow_pent(nua)

integer, intent(inout) :: nma0,nua0,nwa0
integer, intent(out) :: ngrp
integer, intent(out) :: igsize(npts)
integer, intent(out) :: jm(nrows+1,npts)
integer, intent(out) :: ju(nrows,npts)
integer, intent(out) :: nwdivg(npts)

logical :: ready
logical :: ufilled(nper2) ! automatic array

integer :: iper,iper2,iperm,iperp,ig,ileft,iwid
integer :: igs   ! counter for group number

integer :: iu

ufilled(1:nper2) = .false.
ig = 0  ! group number

!---------------------------------------------------------------------------
! FIRST PASS: PROCESS ALL ISOLATED CONVEX POINTS (UNLESS ADJACENT TO 5-NODE)
!---------------------------------------------------------------------------

! Loop over all original M/U points on perimeter of new grid

do iper = 1,nper2
   iperm = iper - 1
   iperp = iper + 1

   if (iper == 1)     iperm = nper2
   if (iper == nper2) iperp = 1

   if (nwdivper(iper ) == 2 .and. &
       nwdivper(iperm) /= 2 .and. &
       nwdivper(iperp) /= 2 .and. &
       nearpent(iper ) /= 1) then

! Skip IPER point if either adjacent U is already filled

   if (ufilled(iperm) .or. ufilled(iper)) cycle

      ig = ig + 1

      igsize(ig) = 2

      jm(1,ig) = imper(iperm)
      jm(2,ig) = imper(iper )
      jm(3,ig) = imper(iperp)

      ju(1,ig) = iuper(iperm)
      ju(2,ig) = iuper(iper )

      nwdivg(ig) = nwdivper(iper)

      ufilled(iperm) = .true.
      ufilled(iper ) = .true.

! Add required number of W,U,M points for this group

      nma0 = nma0 + 1
      nua0 = nua0 + 5
      nwa0 = nwa0 + 4

   endif

enddo

!---------------------------------------------------------------------------
! SECOND PASS: PROCESS ALL CONCAVE POINTS
!---------------------------------------------------------------------------

! Loop over all original M/U points on perimeter of new grid

do iper = 1,nper2
   iperm = iper - 1
   iperp = iper + 1

   if (iper == 1)     iperm = nper2
   if (iper == nper2) iperp = 1

   if (nwdivper(iper) == 4) then

! Skip IPER point if either adjacent U is already filled

      if (ufilled(iperm) .or. ufilled(iper)) cycle

      ig = ig + 1

      igsize(ig) = 2

      jm(1,ig) = imper(iperm)
      jm(2,ig) = imper(iper )
      jm(3,ig) = imper(iperp)

      ju(1,ig) = iuper(iperm)
      ju(2,ig) = iuper(iper )

      nwdivg(ig) = nwdivper(iper)

      ufilled(iperm) = .true.
      ufilled(iper ) = .true.

! Add required number of W,U,M points for this group

      nua0 = nua0 + 2
      nwa0 = nwa0 + 2

   endif

enddo

!---------------------------------------------------------------------------
! THIRD PASS: PROCESS ALL REMAINING POINTS (CENT1 & CENT2)
!---------------------------------------------------------------------------

! Loop over all original M/U points on perimeter of new grid.  (March twice
! around perimeter to allow a delay in processing until immediately after
! reaching the first of the convex or concave points that were filled above.)

ready = .false.

do iper2 = 1,nper2*2
   iper = iper2
   if (iper2 > nper2) iper = iper - nper2

   iperm = iper - 1
   iperp = iper + 1

   if (iper == 1)     iperm = nper2
   if (iper == nper2) iperp = 1

   if (ufilled(iperm) .and. ufilled(iper)) ready = .true.
   
   if (.not. ready) cycle

! Move to next M point if IPERM point is already filled

   if (ufilled(iperm)) cycle

! If IPER U point is already filled or if iper M point is next to a 5-node point,
! fill IPERM point with 'cent1'

   if (ufilled(iper) .or. nearpent(iper) == 1) then

      ig = ig + 1

      igsize(ig) = 1

      jm(1,ig) = imper(iperm)
      jm(2,ig) = imper(iper )

      ju(1,ig) = iuper(iperm)

      nwdivg(ig) = nwdivper(iper)

      ufilled(iperm) = .true.

! Add required number of W,U,M points for this group

      nua0 = nua0 + 1
      nwa0 = nwa0 + 1

! If we got here, U points at IPERM and IPER are both unfilled and IPER M point
! is not next to a 5-node point.  Fill with 'cent2'.

   else
      
      ig = ig + 1

      igsize(ig) = 2

      jm(1,ig) = imper(iperm)
      jm(2,ig) = imper(iper )
      jm(3,ig) = imper(iperp)

      ju(1,ig) = iuper(iperm)
      ju(2,ig) = iuper(iper )

      nwdivg(ig) = nwdivper(iper)

      ufilled(iperm) = .true.
      ufilled(iper ) = .true.

! Add required number of W,U,M points for this group

      nma0 = nma0 + 1
      nua0 = nua0 + 5
      nwa0 = nwa0 + 4

   endif
   
enddo

ngrp = ig

return
end subroutine perim_add2

!===========================================================================

subroutine perim_mrow()

use mem_ijtabs, only: itab_wd, itab_ud
use mem_grid,   only: nwa
use misc_coms,  only: io6

implicit none

integer :: iper, iu, iw, iw1, iw2, iw3, iter, irow, jrow, mrow, mrowh

integer :: mrow_temp(nwa),mrowh_temp(nwa)

! Loop over all W points

do iw = 2,nwa

! Initialize all mrow values to zero.

   itab_wd(iw)%mrow = 0
   mrow_temp(iw)   = 0
   mrowh_temp(iw)  = 0

! Set values on nested grid border to +/- 1.

   iw1 = itab_wd(iw)%iw(1)
   iw2 = itab_wd(iw)%iw(2)
   iw3 = itab_wd(iw)%iw(3)
   
   if     (itab_wd(iw)%mrlw < itab_wd(iw1)%mrlw .or.  &
           itab_wd(iw)%mrlw < itab_wd(iw2)%mrlw .or.  &
           itab_wd(iw)%mrlw < itab_wd(iw3)%mrlw) then

      itab_wd(iw)%mrow = 1
      itab_wd(iw)%mrowh = 1
      mrow_temp(iw) = 1
      mrowh_temp(iw) = 1

   elseif (itab_wd(iw)%mrlw > itab_wd(iw1)%mrlw .or.  &
           itab_wd(iw)%mrlw > itab_wd(iw2)%mrlw .or.  &
           itab_wd(iw)%mrlw > itab_wd(iw3)%mrlw) then

      itab_wd(iw)%mrow = -1
      itab_wd(iw)%mrowh = -1
      mrow_temp(iw) = -1
      mrowh_temp(iw) = -1

   endif

enddo

do irow = 2,10  ! First row already done above
   jrow = mod(irow,2)

   do iw = 2,nwa

      if (itab_wd(iw)%mrow == 0) then

! If IW is adjacent to any other IW cell with nonzero mrow, 
! set mrow for IW cell 

         iw1 = itab_wd(iw)%iw(1)
         iw2 = itab_wd(iw)%iw(2)
         iw3 = itab_wd(iw)%iw(3)

! Check for positive mrow & mrowh values

         mrow = max(itab_wd(iw1)%mrow  &
                   ,itab_wd(iw2)%mrow  &
                   ,itab_wd(iw3)%mrow)

         mrowh = max(itab_wd(iw1)%mrowh  &
                   ,itab_wd(iw2)%mrowh  &
                   ,itab_wd(iw3)%mrowh)

         if (mrow > 0)  mrow_temp (iw) = mrow + jrow
         if (mrowh > 0) mrowh_temp(iw) = mrowh + 1

! Check for negative mrow & mrowh values

         mrow = min(itab_wd(iw1)%mrow  &
                   ,itab_wd(iw2)%mrow  &
                   ,itab_wd(iw3)%mrow)

         mrowh = min(itab_wd(iw1)%mrowh  &
                    ,itab_wd(iw2)%mrowh  &
                    ,itab_wd(iw3)%mrowh)

         if (mrow < 0)  mrow_temp (iw) = mrow - jrow
         if (mrowh < 0) mrowh_temp(iw) = mrowh - 1

      endif

   enddo

   do iw = 2,nwa
      itab_wd(iw)%mrow  = mrow_temp(iw)
      itab_wd(iw)%mrowh = mrowh_temp(iw)
   enddo

enddo

return
end subroutine perim_mrow

!===========================================================================

real function linesegdist(x0,y0,x1,y1,x2,y2)

! Determine distance on 2D Cartesian plane between point (x0,y0) and line
! segment [(x1,y1),(x2,y2)].

implicit none

real, intent(in) :: x0,y0,x1,y1,x2,y2

real :: dsq01,dsq02,dsq12,dline

! Distance squared from (x0,y0) to (x1,y1)

dsq01 = (x1-x0)**2 + (y1-y0)**2

! Distance squared from (x0,y0) to (x2,y2)

dsq02 = (x2-x0)**2 + (y2-y0)**2

! Distance squared from (x1,y1) to (x2,y2)

dsq12 = (x2-x1)**2 + (y2-y1)**2

if (dsq12 < 1.e-6) then
   linesegdist = sqrt(dsq01)
elseif (dsq01 > dsq12 + dsq02) then
   linesegdist = sqrt(dsq02)
elseif (dsq02 > dsq12 + dsq01) then
   linesegdist = sqrt(dsq01)
else 
   dline = abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1)) / sqrt(dsq12)
   linesegdist = dline
endif

return
end function linesegdist





