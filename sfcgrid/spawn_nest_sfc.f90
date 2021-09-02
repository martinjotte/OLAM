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

subroutine spawn_nest_sfc()

  ! This subroutine adds nested grid regions at the beginning of a simulation.

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, &
                          nest_ud_vars, nest_wd_vars, alloc_itabsd, &
                          nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd

  use mem_sfcg,     only: nsfcgrids, nsfcgrdll, sfcgrdrad, sfcgrdlat, &
                          sfcgrdlon, nxp_sfc

  use misc_coms,    only: io6, mdomain, runtype, ngrids
  use consts_coms,  only: pio180, erad, pi1, pi2
  use oname_coms,   only: nl

  implicit none

  type (itab_md_vars), allocatable :: ltab_md(:)
  type (itab_ud_vars), allocatable :: ltab_ud(:)
  type (itab_wd_vars), allocatable :: ltab_wd(:)

  type (nest_ud_vars), allocatable :: nest_ud(:)
  type (nest_wd_vars), allocatable :: nest_wd(:)

  integer :: iu,iw,im,iw1,iw2,im1,im2
  integer :: iu1,iu2,iu3,iu1o,iu2o,iu3o,iu1o_iw1,iu2o_iw1,iu3o_iw1
  integer :: iu4,iu5,iu6,iw3,ngr,nsfcgr,mrlo,mrloo,nwdd,j,npoly,nw

  integer :: nmd0,nud0,nwd0

  real :: expansion

  integer, allocatable :: imper(:) ! Ouside IW index at each perimeter index
  integer, allocatable :: iuper(:) ! Boundary IU index at each perimeter index
  integer :: kma,kua,kwa   ! New M,U,W indices while constructing nest perimeter
  integer :: ngrp   ! Number of perimeter groups
  integer, parameter :: npts = 10000  ! Sufficient array space for # of CM pts along FM perimeter

  integer, allocatable :: jm(:,:),ju(:,:),igsize(:),nwdivg(:)

  real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)

  integer, allocatable :: lista(:), listb(:), jdone(:,:)

  integer :: nper2 = 0 ! Actual # of perimeter pts
  integer :: nccv, jj, minside
  integer :: iper, jm2, ju2, jw2

  integer :: imbeg, ipent, nlista, nlistb, immmm, ndone, ilistb
  integer :: impen, imcent, imnear
  integer :: mlist(6)

  integer, allocatable :: npolyper(:) ! npoly at each perimeter M pt
  integer, allocatable :: nwdivper(:) ! # divided W pts at at each perimeter M pt
  integer, allocatable :: nearpent(:) ! flag = 1 if adjacent outside M pt is poly5

  real :: xp1, xp2, xq1, xq2
  real :: yp1, yp2, yq1, yq2
  integer :: iskip

  ! Make duplicate of current grid dimensions

  nmd0 = nmd
  nud0 = nud
  nwd0 = nwd

  do nsfcgr = 1, nsfcgrids  ! Loop over nested grids
     ngr = nsfcgr + ngrids

     ! Allocate temporary tables

     allocate (nest_ud(nud), nest_wd(nwd))  ! Nest relations
     allocate (jm(3,npts), ju(2,npts))
     allocate (imper(npts), iuper(npts), igsize(npts), nearpent(npts), nwdivg(npts))

     allocate (npolyper(npts), nwdivper(npts))

     ! Copy ITAB information and M-point coordinates to temporary tables

     call move_alloc(itab_md, ltab_md)
     call move_alloc(itab_ud, ltab_ud)
     call move_alloc(itab_wd, ltab_wd)

     call move_alloc(xemd, xem_temp)
     call move_alloc(yemd, yem_temp)
     call move_alloc(zemd, zem_temp)

     ! Locate and flag all W triangles to be subdivided
     ! Define 3 new W indices and 3 new U indices for each W - attach to W.

     do im = 2,nmd

        ! Check whether M location is within specified region for NGR refinement

        call ngr_area(nsfcgr,minside,xem_temp(im),yem_temp(im),zem_temp(im), &
                      nsfcgrdll, sfcgrdrad, sfcgrdlat, sfcgrdlon)

        if (minside == 1) then

           ! W cells surrounding M are to be subdivided

           npoly = ltab_md(im)%npoly

           do j = 1,npoly

              iw = ltab_md(im)%iw(j)

              if (nest_wd(iw)%iw(3) == 0) then

                 nest_wd(iw)%iu(1) = nud0 + 1
                 nest_wd(iw)%iu(2) = nud0 + 2
                 nest_wd(iw)%iu(3) = nud0 + 3

                 nest_wd(iw)%iw(1) = nwd0 + 1
                 nest_wd(iw)%iw(2) = nwd0 + 2
                 nest_wd(iw)%iw(3) = nwd0 + 3

                 nud0 = nud0 + 3
                 nwd0 = nwd0 + 3

              endif

           enddo

        endif

     enddo ! im

     ! Add W points to nested grid region in order to eliminate concavities
     ! (or to eliminate sharp concavities).  This requires iterative procedure

     nwdd = 0  ! Counter of already-existing W points - initialize to zero to
               ! force at least one pass through the following DO WHILE loop

     do while (nwd0 > nwdd)

        nwdd = nwd0

        do im = 2,nmd

           npoly = ltab_md(im)%npoly
           nw = 0  ! Initialize counter for subdivided W points around this M point

           do j = 1,npoly
              iw = ltab_md(im)%iw(j)

              if (nest_wd(iw)%iw(3) > 0) then
                 nw = nw + 1  ! Count up subdivided W points around this M point
              endif
           enddo

           ! Check npoly and nw for illegal values

           if (nw == 0 .or. nw == npoly) cycle
           ! if (npoly == 6 .and. nw < 5) cycle
           ! MJO: This fills in around any pentagon point that happens to be on
           ! a mesh boundary. I will relax this criteria for now:
           if (nw < npoly - 1) cycle

           ! If we got here, then either of the following is true at current M pt:
           !    (1) npoly = 5 and nw > 0 and nw < npoly
           !    (2) nw == 5 (a sharp concavity)
           ! Thus, we must fill in all around current M point

           ! This M point is a concavity so activate remaining unactivated W points around it

           do j = 1,npoly
              iw = ltab_md(im)%iw(j)

              if (nest_wd(iw)%iw(3) == 0) then
                 nest_wd(iw)%iu(1) = nud0 + 1
                 nest_wd(iw)%iu(2) = nud0 + 2
                 nest_wd(iw)%iu(3) = nud0 + 3

                 nest_wd(iw)%iw(1) = nwd0 + 1
                 nest_wd(iw)%iw(2) = nwd0 + 2
                 nest_wd(iw)%iw(3) = nwd0 + 3

                 nud0 = nud0 + 3
                 nwd0 = nwd0 + 3

!                write(io6,*) 'Activating W point ',iw,' to prevent concavity'

              endif
           enddo

        enddo

        call perim_map2(npts, nmd, nud, nwd, nper2, imper, iuper, npolyper, &
                        nwdivper, nearpent, ltab_md, ltab_ud, ltab_wd, nest_wd)

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
                    nest_wd(iw1)%iu(1) = nud0 + 1
                    nest_wd(iw1)%iu(2) = nud0 + 2
                    nest_wd(iw1)%iu(3) = nud0 + 3

                    nest_wd(iw1)%iw(1) = nwd0 + 1
                    nest_wd(iw1)%iw(2) = nwd0 + 2
                    nest_wd(iw1)%iw(3) = nwd0 + 3
                 else
                    nest_wd(iw2)%iu(1) = nud0 + 1
                    nest_wd(iw2)%iu(2) = nud0 + 2
                    nest_wd(iw2)%iu(3) = nud0 + 3

                    nest_wd(iw2)%iw(1) = nwd0 + 1
                    nest_wd(iw2)%iw(2) = nwd0 + 2
                    nest_wd(iw2)%iw(3) = nwd0 + 3
                 endif

                 nud0 = nud0 + 3
                 nwd0 = nwd0 + 3

              enddo

           endif

        enddo

     enddo ! (nwd0 > nwdd)

     ! Print perimeter map
     !
     ! print*, ' '
     ! do j = 1,nper2
     !    write(6,'(a,8i7)') 'perim ',ngr,nper2,j,imper(j),iuper(j),npolyper(j), &
     !                                nwdivper(j),nearpent(j)
     ! enddo
     ! print*, ' '

     ! Nested region should be fully expanded now without any concavities,
     ! 5-edge vertices, or consecutive weak concavities.

     ! Define new vertex index for midpoint of each original U edge that is adjacent
     ! to an original triangle that is being subdivided.  Attach new vertex
     ! index to old U edge.  Also, define new U index for second half of U, and
     ! attach to U.  [Make adjacent to U%m2.]

     do iu = 2,nud

        ! Check whether this U is adjacent to a W that is being subdivided

        iw1 = ltab_ud(iu)%iw(1)
        iw2 = ltab_ud(iu)%iw(2)

        if (nest_wd(iw1)%iw(3) > 0 .or. nest_wd(iw2)%iw(3) > 0) then

           nest_ud(iu)%im = nmd0 + 1
           nest_ud(iu)%iu = nud0 + 1

           nmd0 = nmd0 + 1
           nud0 = nud0 + 1

        endif
     enddo

     ! Save current values of nmd0, nud0, nwd0 prior to adding boundary points

     kma = nmd0
     kua = nud0
     kwa = nwd0

     ! Form groups of original M and U indices on perimeter and increase nmd0, nud0,
     ! and nwd0 according to what groups will require

     call perim_add1(npts,nper2,imper,iuper,nwdivper,nearpent, &
                     nmd0,nud0,nwd0,ngrp,igsize,jm,ju,nwdivg)

     ! Allocate main tables to expanded size
     ! Initialize all neighbor indices to zero

     call alloc_itabsd(nmd0,nud0,nwd0)

     ! Memory copy to main tables

     do im = 1,nmd
        itab_md(im)%imp       = ltab_md(im)%imp
        itab_md(im)%mrlm      = ltab_md(im)%mrlm
        itab_md(im)%mrlm_orig = ltab_md(im)%mrlm_orig
        itab_md(im)%ngr       = ltab_md(im)%ngr

        xemd(im) = xem_temp(im)
        yemd(im) = yem_temp(im)
        zemd(im) = zem_temp(im)
     enddo

     itab_ud(1:nud) = ltab_ud(1:nud)
     itab_wd(1:nwd) = ltab_wd(1:nwd)

     ! Loop over U points and check for those flagged for subdivision

     mrlo = 0

     do iu = 2,nud
        if (nest_ud(iu)%im > 1) then
           im  = nest_ud(iu)%im
           im1 = itab_ud(iu)%im(1)
           im2 = itab_ud(iu)%im(2)

           ! Save mrl value of parent grid for new grid being spawned (saved within IF
           ! block because itab_md(im1)%mrlm gets changed later in DO loop)

           if (mrlo == 0) then
              mrlo = itab_md(im1)%mrlm
           endif

           ! Average coordinates to new M points

           xemd(im) = .5 * (xemd(im1) + xemd(im2))
           yemd(im) = .5 * (yemd(im1) + yemd(im2))
           zemd(im) = .5 * (zemd(im1) + zemd(im2))

           ! Set itab_md()%mrlm = mrlo + 1 at end and mid points of subdivided U edge

           itab_md(im1)%mrlm = mrlo + 1
           itab_md(im2)%mrlm = mrlo + 1
           itab_md(im )%mrlm = mrlo + 1

           itab_md(im )%mrlm_orig = mrlo + 1

           itab_md(im1)%ngr = ngr
           itab_md(im2)%ngr = ngr
           itab_md(im )%ngr = ngr
        endif
     enddo

     ! Contruct tables for new fully subdivided triangles

     mrloo = 0  ! Initialize check variable for uniform mrlw over current nested grid

     do iw = 2,nwd

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
              write(io6,*) 'mrlo,mrloo ',mrlo,mrloo
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

           iu4 = nest_ud(iu1o)%iu
           iu5 = nest_ud(iu2o)%iu
           iu6 = nest_ud(iu3o)%iu

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

           itab_wd(iw)%ngr  = ngr
           itab_wd(iw1)%ngr = ngr
           itab_wd(iw2)%ngr = ngr
           itab_wd(iw3)%ngr = ngr

           if (nest_ud(iu1o)%im > 1) then
              itab_ud(iu1o)%im(2) = nest_ud(iu1o)%im
              itab_ud(iu4)%im(1)  = nest_ud(iu1o)%im
              itab_ud(iu4)%im(2)  = ltab_ud(iu1o)%im(2)
           endif

           if (nest_ud(iu2o)%im > 1) then
              itab_ud(iu2o)%im(2) = nest_ud(iu2o)%im
              itab_ud(iu5)%im(1)  = nest_ud(iu2o)%im
              itab_ud(iu5)%im(2)  = ltab_ud(iu2o)%im(2)
           endif

           if (nest_ud(iu3o)%im > 1) then
              itab_ud(iu3o)%im(2) = nest_ud(iu3o)%im
              itab_ud(iu6)%im(1)  = nest_ud(iu3o)%im
              itab_ud(iu6)%im(2)  = ltab_ud(iu3o)%im(2)
           endif

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

     call perim_fill1(ngr,mrloo,nud0,nwd0,nud,nwd,kma,kua,kwa,&
                      jm,ju,npts,ngrp,igsize,nwdivg, &
                      itab_ud, itab_wd, ltab_ud, ltab_wd, nest_ud, nest_wd)

     ! Copy new counter values

     nmd = nmd0
     nud = nud0
     nwd = nwd0

     deallocate (ltab_md,ltab_ud,ltab_wd)
     deallocate (nest_ud,nest_wd)
     deallocate (xem_temp,yem_temp,zem_temp)

     deallocate (jm,ju)
     deallocate (imper,iuper,igsize,nearpent,nwdivg)
     deallocate (npolyper,nwdivper)

     ! Plot grid lines

     if (.false.) then

        call o_reopnwk()
        call plotback()

        call oplot_set(1)

        do iu = 2,nud
           im1 = itab_ud(iu)%im(1)
           im2 = itab_ud(iu)%im(2)

           call oplot_transform(1,xemd(im1),yemd(im1),zemd(im1),xp1,yp1)
           call oplot_transform(1,xemd(im2),yemd(im2),zemd(im2),xp2,yp2)

           call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

           if (iskip == 1) cycle

           call o_frstpt (xq1,yq1)
           call o_vector (xq2,yq2)
        enddo

        call o_frame()
        call o_clswk()

     endif

     call tri_neighbors(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

     ! Call subroutine to ID W cells just outside and just inside current NGR
     ! border.  This is permanent ID, used in spring dynamics even when new
     ! grids are added.

     call perim_mrow(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

     ! This is the place to do spring dynamics

     call spring_dynamics(1, 0, ngr, nxp_sfc, nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd)

     write(io6,'(/,a,i2)') 'Finished spawning surface grid number ',ngr
     write(io6,'(a,i8)')   ' nmd = ',nmd
     write(io6,'(a,i8)')   ' nud = ',nud
     write(io6,'(a,i8)')   ' nwd = ',nwd

     if (runtype == 'MAKEGRID_PLOT' .and. nsfcgr >= nl%sfcgridplot_base) then

        if (nsfcgr == nl%sfcgridplot_base) then
           call o_reopnwk()
           call plotback()
           call oplot_set_makegrid(2,mrlo)
        endif

        ! Set plot line color (red) and thickness
        call o_gsplci(1)
        call o_gsfaci(1)
        call o_gstxci(1)
        call o_gslwsc(2.5)
        call o_sflush()

        do iu = 2, nud
           iw1 = itab_ud(iu)%iw(1)
           iw2 = itab_ud(iu)%iw(2)

           if ( ( all(itab_md( itab_wd(iw1)%im(1:3) )%ngr == ngr) .and. &
                  any(itab_md( itab_wd(iw2)%im(1:3) )%ngr /= ngr) ) .or. &
                ( all(itab_md( itab_wd(iw2)%im(1:3) )%ngr == ngr) .and. &
                  any(itab_md( itab_wd(iw1)%im(1:3) )%ngr /= ngr) ) ) then

              im1 = itab_ud(iu)%im(1)
              im2 = itab_ud(iu)%im(2)

              call oplot_transform(1,xemd(im1),yemd(im1),zemd(im1),xp1,yp1)
              call oplot_transform(1,xemd(im2),yemd(im2),zemd(im2),xp2,yp2)

              call trunc_segment(xp1,xp2,yp1,yp2,xq1,xq2,yq1,yq2,iskip)

              if (iskip == 1) cycle

              call o_frstpt (xq1,yq1)
              call o_vector (xq2,yq2)

           endif
        enddo

     endif

  enddo   ! end of ngr loop

  if (runtype == 'MAKEGRID_PLOT') then
     call o_sflush
     call mkmap_makegrid()
     call o_frame()
     call o_clswk()
  endif

end subroutine spawn_nest_sfc

!==============================================================================

subroutine perim_add1(npts,nper2,imper,iuper,nwdivper,nearpent, &
                      nmd0,nud0,nwd0,ngrp,igsize,jm,ju,nwdivg)

  implicit none

  integer, intent(in) :: npts,nper2
  integer, intent(in) :: imper(npts)
  integer, intent(in) :: iuper(npts)
  integer, intent(in) :: nwdivper(npts)
  integer, intent(in) :: nearpent(npts)

  integer, intent(inout) :: nmd0,nud0,nwd0
  integer, intent(out) :: ngrp
  integer, intent(out) :: igsize(npts)
  integer, intent(out) :: jm(3,npts)
  integer, intent(out) :: ju(2,npts)
  integer, intent(out) :: nwdivg(npts)

  logical :: ufilled(nper2) ! automatic array

  integer :: iper,iper2,iperm,iperp,ig,ileft,iwid
  integer :: igs   ! counter for group number

  integer :: iu

  ufilled(1:nper2) = .false.
  ig = 0  ! group number

  !---------------------------------------------------------------------------
  ! FIRST PASS: PROCESS ALL CONCAVE POINTS
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

        nud0 = nud0 + 2
        nwd0 = nwd0 + 2

     endif

  enddo

  !---------------------------------------------------------------------------
  ! SECOND PASS: PROCESS ALL REMAINING POINTS (CENT1)
  !---------------------------------------------------------------------------

  ! Loop over all original M/U points on perimeter of new grid.

  do iper = 1,nper2
     iperm = iper - 1

     if (iper == 1) iperm = nper2

     ! Move to next M point if IPERM point is already filled

     if (ufilled(iperm)) cycle

     ! Fill IPERM point with 'cent1'

     ig = ig + 1

     igsize(ig) = 1

     jm(1,ig) = imper(iperm)
     jm(2,ig) = imper(iper )

     ju(1,ig) = iuper(iperm)

     nwdivg(ig) = nwdivper(iper)

     ufilled(iperm) = .true.

! Add required number of W,U,M points for this group

     nud0 = nud0 + 1
     nwd0 = nwd0 + 1

  enddo

  ngrp = ig

end subroutine perim_add1

!===============================================================================

subroutine perim_fill1(ngr,mrloo,nud0,nwd0,nud,nwd,kma,kua,kwa, &
                       jm,ju,npts,ngrp,igsize,nwdivg, &
                       itab_ud, itab_wd, ltab_ud, ltab_wd, nest_ud, nest_wd)

  use mem_delaunay, only: itab_ud_vars, itab_wd_vars, nest_ud_vars, nest_wd_vars

  implicit none

  integer, intent(in)    :: ngr,mrloo
  integer, intent(in) :: nud0, nwd0, nud, nwd
  integer, intent(inout) :: kma,kua,kwa
  integer, intent(in)    :: npts
  integer, intent(inout) :: jm(3,npts)
  integer, intent(inout) :: ju(2,npts)
  integer, intent(in)    :: ngrp
  integer, intent(in)   :: igsize(npts)
  integer, intent(in)    :: nwdivg(npts)

  type (itab_ud_vars), intent(inout) :: itab_ud(nud), ltab_ud(nud)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwd), ltab_wd(nwd)
  type (nest_ud_vars), intent(inout) :: nest_ud(nud)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwd)

  integer :: jmo(3),juo(2)  ! temp storage for output row

  integer :: irow, ig, iu, iw  ! last 2 special

  ! Determine width of transition zone at each point on perimeter

  do ig = 1,ngrp

     if (igsize(ig) == 1) then

        call perim_fill_cent1(ngr,mrloo,nud0,nwd0,nud,nwd,kma,kua,kwa, &
                              jm(1,ig),ju(1,ig),jm(2,ig), &
                              itab_ud, itab_wd, ltab_ud, ltab_wd, nest_ud, nest_wd)

     else ! nwdivg(ig) = 4

        call perim_fill_concave2(ngr,mrloo,nud0,nwd0,nud,nwd,kma,kua,kwa, &
                                 jm(1,ig),ju(1,ig),jm(2,ig),ju(2,ig),jm(3,ig), &
                                 itab_ud, itab_wd, ltab_ud, ltab_wd, nest_ud, nest_wd)
     endif

  enddo

end subroutine perim_fill1

!===========================================================================

subroutine perim_fill_cent1(ngr,mrloo,nud0,nwd0,nud,nwd,kma,kua,kwa, &
                            jm1,ju1,jm2, &
                            itab_ud, itab_wd, ltab_ud, ltab_wd, nest_ud, nest_wd)

  use mem_delaunay, only: itab_ud_vars, itab_wd_vars, nest_ud_vars, nest_wd_vars

  implicit none

  integer, intent(in) :: ngr,mrloo
  integer, intent(in) :: nud0, nwd0, nud, nwd
  integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

  integer, intent(in) :: jm1,jm2  ! Original border M indices
  integer, intent(in) :: ju1      ! Original border U index

  type (itab_ud_vars), intent(inout) :: itab_ud(nud0), ltab_ud(nud)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwd0), ltab_wd(nwd)
  type (nest_ud_vars), intent(inout) :: nest_ud(nud)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwd)

  integer :: ju4,ju7
  integer :: iw1,iw2
  integer :: iu1,iu2,iu3,iu4,iu5
  integer :: im1,im2,im3,im4

  ! Assign new I indices for all M, W points, and for U points that need no checks

  im1 = jm1
  im2 = nest_ud(ju1)%im
  im3 = jm2

  iu4 = kua + 1  ! Newly added point
  iw2 = kwa + 1  ! Newly added point

  ! Increment indices for newly added points

  kua = kua + 1
  kwa = kwa + 1

  ! Determine orientation of positive ju1 direction

  if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
     iw1 = ltab_ud(ju1)%iw(1)
     iu1 = ju1
     iu2 = nest_ud(ju1)%iu

     itab_ud(iu1)%iw(1) = iw1
     itab_ud(iu2)%iw(1) = iw2
  else                              ! Positive ju1 points outward
     iw1 = ltab_ud(ju1)%iw(2)
     iu1 = nest_ud(ju1)%iu
     iu2 = ju1

     itab_ud(iu1)%iw(2) = iw1
     itab_ud(iu2)%iw(2) = iw2
  endif

  if (ju1 == ltab_wd(iw1)%iu(1)) then
     iu3 = ltab_wd(iw1)%iu(2)
     iu5 = ltab_wd(iw1)%iu(3)
  elseif (ju1 == ltab_wd(iw1)%iu(2)) then
     iu3 = ltab_wd(iw1)%iu(3)
     iu5 = ltab_wd(iw1)%iu(1)
  else
     iu3 = ltab_wd(iw1)%iu(1)
     iu5 = ltab_wd(iw1)%iu(2)
  endif

  ! Fill remaining neighbor indices for U and W points

  ! nothing needed for iu3

  if (iw1 == itab_ud(iu5)%iw(1)) then
     im4 = itab_ud(iu5)%im(2)
     itab_ud(iu5)%iw(1) = iw2
  else
     im4 = itab_ud(iu5)%im(1)
     itab_ud(iu5)%iw(2) = iw2
  endif

  itab_ud(iu4)%im(1) = im2
  itab_ud(iu4)%im(2) = im4
  itab_ud(iu4)%iw(1) = iw1
  itab_ud(iu4)%iw(2) = iw2

  itab_wd(iw1)%iu(1) = iu1
  itab_wd(iw1)%iu(2) = iu3
  itab_wd(iw1)%iu(3) = iu4

  itab_wd(iw2)%iu(1) = iu2
  itab_wd(iw2)%iu(2) = iu4
  itab_wd(iw2)%iu(3) = iu5

  ! Fill mrl for new W points

  if (itab_wd(iw1)%mrlw /= mrloo) go to 5

  itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
  itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig

  ! NEW for spawn_nest_sfc

  itab_wd(iw1)%ngr = ngr
  itab_wd(iw2)%ngr = ngr

  return

  5 continue

  write(6,*) ''
  write(6, '(1x,A)')    'Error in subroutine perim_fill_cent1.'
  write(6, '(1x,A,I0)') 'Current nested grid ', ngr
  write(6, '(1x,A)')    'crosses pre-existing grid boundary.'
  write(6, '(1x,A,I0)') 'iw1 = ',iw1
  stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_cent1

!===========================================================================

subroutine perim_fill_concave2(ngr,mrloo,nud0,nwd0,nud,nwd,kma,kua,kwa, &
                               jm1,ju1,jm2,ju2,jm3, &
                               itab_ud, itab_wd, ltab_ud, ltab_wd, nest_ud, nest_wd)

  use mem_delaunay, only: itab_ud_vars, itab_wd_vars, nest_ud_vars, nest_wd_vars

  implicit none

  integer, intent(in)    :: ngr,mrloo
  integer, intent(in)    :: nud0, nwd0, nud, nwd
  integer, intent(inout) :: kma, kua, kwa   ! Index values for latest added points

  integer, intent(in) :: jm1,jm2,jm3   ! Original M indices

  integer, intent(in) :: ju1,ju2       ! Original U indices

  type (itab_ud_vars), intent(inout) :: itab_ud(nud0), ltab_ud(nud)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwd0), ltab_wd(nwd)
  type (nest_ud_vars), intent(inout) :: nest_ud(nud)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwd)

  integer :: jm4
  integer :: ju3,ju4,ju5
  integer :: jw1,jw2

  integer :: im1,im2,im3,im4,im5,im6             ! Temporary new M indices
  integer :: iw1,iw2,iw3,iw4                     ! Temporary new W indices
  integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9 ! Temporary new U indices

  ! Newly added edges

  iu6 = kua + 1
  iu8 = kua + 2

  ! Newly added areas

  iw2 = kwa + 1
  iw3 = kwa + 2

  ! Increment indices for newly added points

  kua = kua + 2
  kwa = kwa + 2

  ! Determine orientation of positive ju1 direction to get iw1 (=jw1)

  if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points into fine mesh
     jw1 = ltab_ud(ju1)%iw(1)
     jw2 = ltab_ud(ju1)%iw(4)

     ju3 = ltab_ud(ju1)%iu(1)
     ju4 = ltab_ud(ju1)%iu(2)
     ju5 = ltab_ud(ju1)%iu(7)

     iu1 = ju1
     iu2 = nest_ud(ju1)%iu

     itab_ud(iu1)%iw(1) = jw1  ! not needed?
     itab_ud(iu2)%iw(1) = iw3
  else                              ! Positive ju1 points out from fine mesh
     jw1 = ltab_ud(ju1)%iw(2)
     jw2 = ltab_ud(ju1)%iw(5)

     ju3 = ltab_ud(ju1)%iu(4)
     ju4 = ltab_ud(ju1)%iu(3)
     ju5 = ltab_ud(ju1)%iu(10)

     iu1 = nest_ud(ju1)%iu
     iu2 = ju1

     itab_ud(iu1)%iw(2) = jw1
     itab_ud(iu2)%iw(2) = iw3
  endif

  !write(io6,'(a,20i7)') 'pf11 ',jm1,jm2,jm3,ju1,ju2, &
  !   ltab_ud(ju1)%im(1),ltab_ud(ju1)%im(2),ltab_ud(ju2)%im(1),ltab_ud(ju2)%im(2), &
  !   ltab_ud(ju1)%iw(1),ltab_ud(ju1)%iw(2),ltab_ud(ju2)%iw(1),ltab_ud(ju2)%iw(2), &
  !   nest_ud(ju1)%iu

  ! Determine orientation of positive ju2 direction to get iw4 (=jw2)

  if (jm2 == ltab_ud(ju2)%im(1)) then  ! Positive ju2 points inward
     iu3 = ju2
     iu4 = nest_ud(ju2)%iu

     itab_ud(iu3)%iw(1) = iw3
     itab_ud(iu4)%iw(1) = jw2
  else                              ! Positive ju2 points outward
     iu3 = nest_ud(ju2)%iu
     iu4 = ju2

     itab_ud(iu3)%iw(2) = iw3  ! not needed?
     itab_ud(iu4)%iw(2) = jw2
  endif

  ! Determine neighbors of ju4

  if (jm2 == ltab_ud(ju4)%im(1)) then
     jm4 = ltab_ud(ju4)%im(2)
  else
     jm4 = ltab_ud(ju4)%im(1)
  endif

  ! Assign remaining I indices that already existed as J indices

  ! (iu1:iu4 already assigned above)

  ! Assign new I indices for all M, W points, and for U points that need no checks

  im1 = jm1
  im2 = nest_ud(ju1)%im
  im3 = jm2
  im4 = nest_ud(ju2)%im
  im5 = jm3
  im6 = jm4

  iu5 = ju3
  iu7 = ju4
  iu9 = ju5

  iw1 = jw1
  iw4 = jw2

  !print*, ' '
  !write(io6,'(a,20i7)') 'c2a.1 ',im1,im2,im3,im4,im5,im6
  !write(io6,'(a,20i7)') 'c2a.2 ',iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9
  !write(io6,'(a,20i7)') 'c2a.3 ',iw1,iw2,iw3,iw4
  !write(io6,'(a,20i7)') 'c2a.4  ',itab_ud(iu1)%im(1),itab_ud(iu1)%im(2), &
  !                                itab_ud(iu1)%iw(1),itab_ud(iu1)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.5  ',itab_ud(iu2)%im(1),itab_ud(iu2)%im(2), &
  !                                itab_ud(iu2)%iw(1),itab_ud(iu2)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.6  ',itab_ud(iu3)%im(1),itab_ud(iu3)%im(2), &
  !                                itab_ud(iu3)%iw(1),itab_ud(iu3)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.7  ',itab_ud(iu4)%im(1),itab_ud(iu4)%im(2), &
  !                                itab_ud(iu4)%iw(1),itab_ud(iu4)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.8  ',itab_ud(iu5)%im(1),itab_ud(iu5)%im(2), &
  !                                itab_ud(iu5)%iw(1),itab_ud(iu5)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.9  ',itab_ud(iu6)%im(1),itab_ud(iu6)%im(2), &
  !                                itab_ud(iu6)%iw(1),itab_ud(iu6)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.10 ',itab_ud(iu7)%im(1),itab_ud(iu7)%im(2), &
  !                                itab_ud(iu7)%iw(1),itab_ud(iu7)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.11 ',itab_ud(iu8)%im(1),itab_ud(iu8)%im(2), &
  !                                itab_ud(iu8)%iw(1),itab_ud(iu8)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.12 ',itab_ud(iu9)%im(1),itab_ud(iu9)%im(2), &
  !                                itab_ud(iu9)%iw(1),itab_ud(iu9)%iw(2)
  !write(io6,'(a,20i7)') 'c2a.13 ',itab_wd(iw1)%iu(1),itab_wd(iw1)%iu(2), &
  !                                itab_wd(iw1)%iu(3)
  !write(io6,'(a,20i7)') 'c2a.14 ',itab_wd(iw2)%iu(1),itab_wd(iw2)%iu(2), &
  !                                itab_wd(iw2)%iu(3)
  !write(io6,'(a,20i7)') 'c2a.15 ',itab_wd(iw3)%iu(1),itab_wd(iw3)%iu(2), &
  !                                itab_wd(iw3)%iu(3)
  !write(io6,'(a,20i7)') 'c2a.16 ',itab_wd(iw4)%iu(1),itab_wd(iw4)%iu(2), &
  !                                itab_wd(iw4)%iu(3)
  !print*, ' '



  ! Fill remaining neighbor indices for U and W points

  ! nothing needed for iu5

  itab_ud(iu6)%im(1) = im2
  itab_ud(iu6)%im(2) = im6
  itab_ud(iu6)%iw(1) = iw1
  itab_ud(iu6)%iw(2) = iw2

  itab_ud(iu7)%im(1) = im2
  itab_ud(iu7)%im(2) = im4
  itab_ud(iu7)%iw(1) = iw2
  itab_ud(iu7)%iw(2) = iw3

  itab_ud(iu8)%im(1) = im4
  itab_ud(iu8)%im(2) = im6
  itab_ud(iu8)%iw(1) = iw2
  itab_ud(iu8)%iw(2) = iw4

  ! nothing needed for iu9

  itab_wd(iw1)%iu(1) = iu1
  itab_wd(iw1)%iu(2) = iu5
  itab_wd(iw1)%iu(3) = iu6

  itab_wd(iw2)%iu(1) = iu6
  itab_wd(iw2)%iu(2) = iu8
  itab_wd(iw2)%iu(3) = iu7

  itab_wd(iw3)%iu(1) = iu2
  itab_wd(iw3)%iu(2) = iu7
  itab_wd(iw3)%iu(3) = iu3

  itab_wd(iw4)%iu(1) = iu4
  itab_wd(iw4)%iu(2) = iu8
  itab_wd(iw4)%iu(3) = iu9

  ! Fill mrl for new W points

  if (itab_wd(iw1)%mrlw /= mrloo) go to 5
  if (itab_wd(iw4)%mrlw /= mrloo) go to 5

  itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
  itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig

  itab_wd(iw3)%mrlw      = itab_wd(iw1)%mrlw
  itab_wd(iw3)%mrlw_orig = itab_wd(iw1)%mrlw_orig

  ! NEW for spawn_nest_sfc

  itab_wd(iw1)%ngr = ngr
  itab_wd(iw2)%ngr = ngr
  itab_wd(iw3)%ngr = ngr
  itab_wd(iw4)%ngr = ngr

  return

  5 continue

  write(6,*) 'In subroutine perim_fill_concave2, current nested grid ',ngr
  write(6,*) 'crosses pre-existing grid boundary.'
  write(6,*) 'iw1 = ',iw1
  write(6,*) 'stopping model'
  stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_concave2

