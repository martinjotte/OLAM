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

subroutine spawn_nest()

! This subroutine adds nested grid regions at the beginning of a simulation.
! Later will make modified version to add nested grid region(s) during a
! simulation.

  use mem_ijtabs,   only: mloops, mrls, &
                          jtm_grid, jtu_grid, jtv_grid, jtw_grid, &
                          jtm_init, jtu_init, jtv_init, jtw_init, &
                          jtm_prog, jtu_prog, jtv_prog, jtw_prog, &
                          jtm_wadj, jtu_wadj, jtv_wadj, jtw_wadj, &
                          jtm_wstn, jtu_wstn, jtv_wstn, jtw_wstn, &
                          jtm_lbcp, jtu_lbcp, jtv_lbcp, jtw_lbcp, &
                          jtm_vadj, jtu_wall, jtv_wall, jtw_vadj

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, &
                          nest_ud_vars, nest_wd_vars, alloc_itabsd, &
                          itab_md, itab_ud, itab_wd, copy_tri_grid, &
                          xemd, yemd, zemd, nmd, nud, nwd

  use mem_grid,     only: impent, nrows, mrows

  use misc_coms,    only: io6, ngrids, mdomain, nxp, ngrdll, grdrad, &
                          grdlat, grdlon

  use consts_coms,  only: pio180, erad, pi1, pi2
  use mem_sfcg,     only: nsfcgrid_root

  implicit none

  type (itab_md_vars), allocatable :: ltab_md(:)
  type (itab_ud_vars), allocatable :: ltab_ud(:)
  type (itab_wd_vars), allocatable :: ltab_wd(:)

  type (nest_ud_vars), allocatable :: nest_ud(:)
  type (nest_wd_vars), allocatable :: nest_wd(:)

  integer :: iu,iw,im,iw1,iw2,im1,im2
  integer :: iu1,iu2,iu3,iu1o,iu2o,iu3o,iu1o_iw1,iu2o_iw1,iu3o_iw1
  integer :: iu4,iu5,iu6,iw3,ngr,mrlo,mrloo,nwda,j,npoly,nw

  integer :: nmd0,nud0,nwd0

  integer, allocatable :: imper(:) ! Ouside IW index at each perimeter index
  integer, allocatable :: iuper(:) ! Boundary IU index at each perimeter index
  integer :: kma,kua,kwa   ! New M,U,W indices while constructing nest perimeter

  integer, parameter :: npts = 10000  ! Sufficient array space for # of CM pts along FM perimeter

  integer, allocatable :: jm(:,:),ju(:,:),igsize(:),nwdivg(:)

  real, allocatable :: xem_temp(:),yem_temp(:),zem_temp(:)

  integer, allocatable :: lista(:), listb(:), jdone(:,:)

  integer :: nper2 = 0 ! Actual # of perimeter pts
  integer :: jj, minside
  integer :: iper, jm2, ju2, jw2

  integer :: imbeg, ipent, nlista, nlistb, immmm, ndone, ilistb
  integer :: impen, imcent, imnear
  integer :: mlist(6)
  real :: reg, xeg, yeg, zeg, dist, distmin

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

  ! This routine spawns nest using the Method C refined mesh algorithm.  It forces
  ! straight segments of the refined mesh "boundary" to have lengths that are multiples
  ! of 3 triangle edges on the unrefined side.  Method C also fits exactly 3 transition
  ! rows across a gap that was originally 2 coarse grid rows wide, centered on
  ! the aforementioned refined mesh "boundary".

  mrows = 3

  do ngr = 2, ngrids  ! Loop over nested grids

     write(io6,'(/,a,i0)') 'Spawning grid number ',ngr

     ! Allocate temporary tables

     allocate (nest_ud(nud), nest_wd(nwd))  ! Nest relations
!    allocate (xem_temp(nmd), yem_temp(nmd), zem_temp(nmd))

     allocate (jm(nrows+1,npts), ju(nrows,npts))
     allocate (imper(npts), iuper(npts), igsize(npts), nearpent(npts), nwdivg(npts))

     allocate (npolyper(npts), nwdivper(npts))

     ! Copy ITAB information and M-point coordinates to temporary tables

     call move_alloc(itab_md, ltab_md)
     call move_alloc(itab_ud, ltab_ud)
     call move_alloc(itab_wd, ltab_wd)

     call move_alloc(xemd, xem_temp)
     call move_alloc(yemd, yem_temp)
     call move_alloc(zemd, zem_temp)

     allocate (lista(nmd), listb(nmd), jdone(6,nmd))

     ! Find closest M point to first specified NGR center point.

     ! Get earth coordinates for [grdlat(ngr,1),grdlon(ngr,1)]

     if (mdomain < 2) then
        zeg = erad * sin(grdlat(ngr,1) * pio180)
        reg = erad * cos(grdlat(ngr,1) * pio180)
        xeg = reg  * cos(grdlon(ngr,1) * pio180)
        yeg = reg  * sin(grdlon(ngr,1) * pio180)
     else
        zeg = 0.
        xeg = grdlon(ngr,1)
        yeg = grdlat(ngr,1)
     endif

     ! Initialize distance

     distmin = 1.e12

     ! Loop over all M points in domain

     do im = 2,nmd
        dist = sqrt((xem_temp(im) - xeg) ** 2 &
                  + (yem_temp(im) - yeg) ** 2 &
                  + (zem_temp(im) - zeg) ** 2)

        if (distmin > dist) then
           distmin = dist
           imcent = im
        endif
     enddo

     ! If using global grid (MDOMAIN < 2), search over 12 original ipent points and
     ! determine if any lies inside current NGR refinement area.  If not, determine
     ! whether any ipent point is close to NGR refinement area AND has the same
     ! mrlm value as the imcent point.

     imbeg = 0
     impen = 0

     if (mdomain < 2) then

        do ipent = 1,12
           im = impent(ipent)

           call ngr_area(ngr,minside,xem_temp(im),yem_temp(im),zem_temp(im), &
                         ngrdll, grdrad, grdlat, grdlon)

           if (minside == 1) then
              imbeg = im
              exit
           elseif (minside == 2 .and. &
              ltab_md(im)%mrlm == ltab_md(imcent)%mrlm) then

              impen = im
           endif
        enddo

     endif

     ! If imbeg is zero but impen is not, then find nearest M point to impen
     ! that is inside NGR refinement area.

     if (imbeg == 0 .and. impen > 0) then

        ! Initialize distance

        distmin = 1.e12

        ! Loop over all M points in domain

        do im = 2,nmd

           ! Check whether M location is within specified region for NGR refinement

           call ngr_area(ngr,minside,xem_temp(im),yem_temp(im),zem_temp(im), &
                         ngrdll, grdrad, grdlat, grdlon)

           if (minside == 1) then

              dist = sqrt((xem_temp(im) - xem_temp(impen)) ** 2 &
                        + (yem_temp(im) - yem_temp(impen)) ** 2 &
                        + (zem_temp(im) - zem_temp(impen)) ** 2)

              if (distmin > dist) then
                 distmin = dist
                 imnear = im
              endif

           endif

        enddo

        ! Now, march three M points at a time from impen toward imnear until finding
        ! an M point that IS inside the NGR refinement area.

        im = impen

        searchloop: do

           jdone(:,im) = 0  ! (6,nmd)
           mlist(1:6) = 0

           call thirdm(nmd, nud, im, jdone, mlist, ltab_md, ltab_ud)

           ! Search through THIRDM neighbors of IM that are in current mlist

           distmin = 1.e12

           do j = 1,6
              immmm = mlist(j)
              if (immmm > 1) then

                 call ngr_area(ngr,minside,xem_temp(immmm),yem_temp(immmm), &
                                                           zem_temp(immmm), &
                               ngrdll, grdrad, grdlat, grdlon)

                 if (minside == 1) then
                    imbeg = immmm
                    exit searchloop
                 endif

                 dist = sqrt((xem_temp(immmm) - xem_temp(imnear)) ** 2 &
                           + (yem_temp(immmm) - yem_temp(imnear)) ** 2 &
                           + (zem_temp(immmm) - zem_temp(imnear)) ** 2)

                 if (distmin > dist) then
                    distmin = dist
                    im = immmm
                 endif

              endif ! immmm > 1

           enddo ! j,immmm

        enddo searchloop

     endif   ! (imbeg == 0 .and. impen > 0)

     ! If imbeg is still zero, then use imcent as starting point

     if (imbeg == 0) imbeg = imcent

     ! Now that starting point IMBEG has been determined, build full list of
     ! M points that are inside specified NGR refinement area

     ! Initialize quantities for search

     jdone(:,:) = 0

     nlista = 1
     lista(nlista) = imbeg

     nlistb = 0

     ! Loop over points in LISTA, as long as any exist

     do while (nlista > 0)

     ! Initialize to zero the array of up to 6 potential THIRDM neighbors of IM

        mlist(1:6) = 0

        ! Search for THIRDM neighbors of IM (the last value in LISTA), rejecting
        ! paths that have already been done

        im = lista(nlista)

        call thirdm(nmd, nud, im, jdone, mlist, ltab_md, ltab_ud)

        ! Now that IM's THIRDM neighbors have been found, remove IM from LISTA and add
        ! it to LISTB.

        lista(nlista) = 0
        nlista = nlista - 1

        nlistb = nlistb + 1
        listb(nlistb) = im

        ! Search through THIRDM neighbors of IM that are in current mlist

        do j = 1,6
           immmm = mlist(j)
           if (immmm > 1) then

              ! Loop over jdone values of IMMMM point.  If two or more of these edges have already
              ! been traversed, then this IMMMM point has already been added to LISTA.
              ! Thus, do not add it again.  (One edge was traversed in latest call to thirdm, but that
              ! traverse has not yet been considered for lista.)

              ndone = 0
              do jj = 1,6
                 if (jdone(jj,immmm) == 1) ndone = ndone + 1
              enddo

              if (ndone < 2) then

                 ! Check whether IMMMM point is inside NGR refinement area

                 call ngr_area(ngr,minside,xem_temp(immmm),yem_temp(immmm), &
                                                           zem_temp(immmm), &
                               ngrdll, grdrad, grdlat, grdlon)

                 if (minside == 1) then

                    ! IMMMM point is inside NGR refinement area; add it to LISTA

                    nlista = nlista + 1
                    lista(nlista) = immmm

                 endif ! minside == 1

              endif ! ndone < 2

           endif ! immmm > 1

        enddo ! j,immmm

     enddo ! nlista > 0

     ! Now, listb contains the full list of M point indices, numbering nlistb.

     ! Flag W points that are within radius_3 polygons of M points in listb.
     ! (We can use a nest_wd()%iw(3) as a flag.)

     nest_wd(:)%iw(3) = 0

     ! Loop over IM points that are in listb

     do ilistb = 1,nlistb
        im = listb(ilistb)

        call fill_rad3(nmd, nwd, im, ltab_md, ltab_wd, nest_wd)
     enddo

     deallocate (lista,listb,jdone)

     ! Add W points to nested grid region in order to eliminate concavities
     ! (or to eliminate sharp concavities).  This requires iterative procedure

     nwda = 0  ! Counter of already-existing W points - initialize to zero to
               ! force at least one pass through the following DO WHILE loop

     do while (nwd0 > nwda)

        nwda = nwd0

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

           call fill_rad3(nmd, nwd, im, ltab_md, ltab_wd, nest_wd)

           nwd0 = nwd0 + 1 ! just to keep DO WHILE going as long as necessary

        enddo

        call perim_map2(npts, nmd, nud, nwd, nper2, imper, iuper, npolyper, &
                        nwdivper, nearpent, ltab_md, ltab_ud, ltab_wd, nest_wd)

     enddo ! (nwd0 > nwda)

     ! Print perimeter map
     !  print*, ' '
     !  do j = 1,nper2
     !     write(6,'(a,8i7)') 'perim ',ngr,nper2,j,imper(j),iuper(j),npolyper(j), &
     !                                 nwdivper(j),nearpent(j)
     !  enddo
     !  print*, ' '

     ! Nested region should be fully expanded now without any concavities,
     ! 5-edge vertices, or consecutive weak concavities.

     ! Reset subdivide flag to -1 for W triangle adjacent to center segment
     ! of each set of 3 segments.  This will suppress subdivision of W triangle
     ! and also of 3 adjacent U edges.  Doing this will cause remaining interior
     ! subdivisions to exactly match memory requirement.

     do iper = 1,nper2-2,3
        jm2 = imper(iper+1)
        ju2 = iuper(iper+1)

        if (jm2 == ltab_ud(ju2)%im(1)) then
           jw2 = ltab_ud(ju2)%iw(2)
        else
           jw2 = ltab_ud(ju2)%iw(1)
        endif

        nest_wd(jw2)%iw(3) = -1
     enddo

     ! Reset nwd0 counter for actual count

     nwd0 = nwd

     ! Loop over all W points, counting up those with nest_wd()%iw(3) still flagged.
     ! Reset all nest_wd members according to current count

     do iw = 2,nwd

        if (nest_wd(iw)%iw(3) > 0) then

!          write(io6,*) 'subdividing W cell ',iw

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

     ! Define new vertex index for midpoint of each original U edge that is adjacent
     ! to an original triangle that is being subdivided, unless it is also adjacent
     ! to an original triangle with subdivide flag = -1.  Attach new vertex
     ! index to old U edge.  Also, define new U index for second half of U, and
     ! attach to U.  [Make adjacent to U%m2.]

     do iu = 2,nud

        ! Check whether this U is adjacent to a W that is being subdivided

        iw1 = ltab_ud(iu)%iw(1)
        iw2 = ltab_ud(iu)%iw(2)

        if (nest_wd(iw1)%iw(3) > 0 .or. nest_wd(iw2)%iw(3) > 0) then

           if (nest_wd(iw1)%iw(3) < 0 .or. nest_wd(iw2)%iw(3) < 0) then

              nest_ud(iu)%im = 1
              nest_ud(iu)%iu = iu

           else

              nest_ud(iu)%im = nmd0 + 1
              nest_ud(iu)%iu = nud0 + 1

              nmd0 = nmd0 + 1
              nud0 = nud0 + 1

           endif

        endif
     enddo

     ! Save current values of nmd0, nud0, nwd0 prior to adding boundary points

     kma = nmd0
     kua = nud0
     kwa = nwd0

     ! Allocate main tables to expanded size
     ! Initialize all neighbor indices to zero

     call alloc_itabsd(nmd0,nud0,nwd0)

     ! Memory copy to main tables

     do im = 1,nmd
        itab_md(im)%loop(1:mloops) = ltab_md(im)%loop(1:mloops)
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

     call perim_fill3(ngr, nmd, nud, nwd, mrlo, nper2, imper, iuper, ltab_ud, nest_ud, nest_wd)

     ! Fill itabs loop tables for newly spawned points

     do im = nmd+1,nmd0
        itab_md(im)%imp = im
        call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_wstn, 0)
     enddo

     do iu = nud+1,nud0
        itab_ud(iu)%iup = iu
        call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
     enddo

     do iw = nwd+1,nwd0
        itab_wd(iw)%iwp = iw
        call wdloopf('f',iw, jtw_grid, jtw_vadj, 0, 0, 0, 0)
     enddo

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

     call spring_dynamics(3, 1, ngr, nxp, nmd, nud, nwd, xemd, yemd, zemd, &
                          itab_md, itab_ud, itab_wd)

     write(io6,'(/,a,i0)') 'Finished spawning grid number ',ngr
     write(io6,'(a,i0)')   ' nma = ',nmd
     write(io6,'(a,i0)')   ' nua = ',nud
     write(io6,'(a,i0)')   ' nwa = ',nwd

     ! If independent refinement of surface grid will be done and uses current NGR
     ! refinement of ATM grid as its root, copy ATM grid quantities to surface grid

     if (nsfcgrid_root == ngr) then
        call copy_tri_grid()
     endif

  enddo   ! end of ngr loop

  ! Set mrls equal to maximum mrlw value

  do iw = 2,nwd
     if (mrls < itab_wd(iw)%mrlw) mrls = itab_wd(iw)%mrlw
  enddo

end subroutine spawn_nest

!==============================================================================

subroutine perim_fill3(ngr, nma, nua, nwa, mrlo, nper, imper, iuper, ltab_ud, nest_ud, nest_wd)

  use mem_delaunay, only: itab_ud_vars, nest_ud_vars, nest_wd_vars, &
                          itab_md, itab_ud, itab_wd, xemd, yemd, zemd
  implicit none

  integer, intent(in) :: ngr, nma, nua, nwa, mrlo, nper, imper(nper), iuper(nper)

  type (itab_ud_vars), intent(inout) :: ltab_ud(nua)
  type (nest_ud_vars), intent(inout) :: nest_ud(nua)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwa)

  integer :: jm1, jm2, jm3, ju1, ju2, ju3, ku1, ku2, ku3

  integer :: im5,im12,im13,im17,im18,im19,im20,im24
  integer :: im16,im21,im22,im23,im25,im26

  integer :: iu15,iu16,iu25,iu26,iu33,iu34,iu35,iu41,iu42,iu43
  integer :: iu44,iu45,iu48,iu49,iu50,iu51
  integer :: iu46,iu53

  integer :: iw6,iw7,iw8,iw9,iw19,iw20,iw21,iw27,iw29,iw31
  integer :: iw26,iw28,iw30,iw32

  integer :: iper

  do iper = 1,nper,3

     ! Inventory of CM mesh points

     jm1 = imper(iper)
     jm2 = imper(iper+1)
     jm3 = imper(iper+2)

     ju1 = iuper(iper)
     ju2 = iuper(iper+1)
     ju3 = iuper(iper+2)

     if (jm1 == ltab_ud(ju1)%im(1)) then
        iu41 = ju1
        iu42 = nest_ud(ju1)%iu
        iu46 = ltab_ud(ju1)%iu(5)
        iw26 = ltab_ud(ju1)%iw(3)
        iw27 = ltab_ud(ju1)%iw(1)
     else
        iu41 = nest_ud(ju1)%iu
        iu42 = ju1
        iu46 = ltab_ud(ju1)%iu(12)
        iw26 = ltab_ud(ju1)%iw(6)
        iw27 = ltab_ud(ju1)%iw(2)
     endif

     if (jm2 == ltab_ud(ju2)%im(1)) then

        iu49 = ltab_ud(ju2)%iu(1)
        iu50 = ltab_ud(ju2)%iu(2)
        iu34 = ltab_ud(ju2)%iu(3)
        iu35 = ltab_ud(ju2)%iu(4)
        iu48 = ltab_ud(ju2)%iu(5)
        iu51 = ltab_ud(ju2)%iu(8)

        iw6 = ltab_ud(ju2)%iw(5)
        iw9 = ltab_ud(ju2)%iw(6)
        iw29 = ltab_ud(ju2)%iw(1)
        iw20 = ltab_ud(ju2)%iw(2)
        iw28 = ltab_ud(ju2)%iw(3)
        iw30 = ltab_ud(ju2)%iw(4)

     else

        iu49 = ltab_ud(ju2)%iu(4)
        iu50 = ltab_ud(ju2)%iu(3)
        iu34 = ltab_ud(ju2)%iu(2)
        iu35 = ltab_ud(ju2)%iu(1)
        iu48 = ltab_ud(ju2)%iu(12)
        iu51 = ltab_ud(ju2)%iu(9)

        iw6 = ltab_ud(ju2)%iw(4)
        iw9 = ltab_ud(ju2)%iw(3)
        iw29 = ltab_ud(ju2)%iw(2)
        iw20 = ltab_ud(ju2)%iw(1)
        iw28 = ltab_ud(ju2)%iw(6)
        iw30 = ltab_ud(ju2)%iw(5)

     endif

     if (jm3 == ltab_ud(ju3)%im(1)) then
        im21 = ltab_ud(ju3)%im(2)
        iu44 = ju3
        iu45 = nest_ud(ju3)%iu
        iu53 = ltab_ud(ju3)%iu(8)
        iw31 = ltab_ud(ju3)%iw(1)
        iw32 = ltab_ud(ju3)%iw(4)
     else
        im21 = ltab_ud(ju3)%im(1)
        iu44 = nest_ud(ju3)%iu
        iu45 = ju3
        iu53 = ltab_ud(ju3)%iu(9)
        iw31 = ltab_ud(ju3)%iw(2)
        iw32 = ltab_ud(ju3)%iw(5)
     endif

     im16 = jm1
     im17 = nest_ud(ju1)%im
     im18 = jm2
     im19 = jm3
     im20 = nest_ud(ju3)%im
     iu43 = ju2

     ku1 = nest_wd(iw6)%iu(1)
     ku2 = nest_wd(iw6)%iu(2)
     ku3 = nest_wd(iw6)%iu(3)

     if (itab_ud(ku1)%im(1) > 1 .and. itab_ud(ku1)%im(2) > 1) then
        iu25 = ku2
        iu15 = ku3
     elseif (itab_ud(ku2)%im(1) > 1 .and. itab_ud(ku2)%im(2) > 1) then
        iu25 = ku3
        iu15 = ku1
     elseif (itab_ud(ku3)%im(1) > 1 .and. itab_ud(ku3)%im(2) > 1) then
        iu25 = ku1
        iu15 = ku2
     endif

     if (itab_ud(iu15)%iw(1) == iw6) then
        iw7 = itab_ud(iu15)%iw(2)
     else
        iw7 = itab_ud(iu15)%iw(1)
     endif

     if (itab_ud(iu25)%iw(1) == iw6) then
        iw19 = itab_ud(iu25)%iw(2)
        im12 = itab_ud(iu25)%im(2) ! to reposition xyzemd(12)
     else
        iw19 = itab_ud(iu25)%iw(1)
        im12 = itab_ud(iu25)%im(1) ! to reposition xyzemd(12)
     endif

     ku1 = nest_wd(iw9)%iu(1)
     ku2 = nest_wd(iw9)%iu(2)
     ku3 = nest_wd(iw9)%iu(3)

     if (itab_ud(ku1)%im(1) > 1 .and. itab_ud(ku1)%im(2) > 1) then
        iu16 = ku2
        iu26 = ku3
     elseif (itab_ud(ku2)%im(1) > 1 .and. itab_ud(ku2)%im(2) > 1) then
        iu16 = ku3
        iu26 = ku1
     elseif (itab_ud(ku3)%im(1) > 1 .and. itab_ud(ku3)%im(2) > 1) then
        iu16 = ku1
        iu26 = ku2
     endif

     if (itab_ud(iu16)%iw(1) == iw9) then
        iw8 = itab_ud(iu16)%iw(2)
     else
        iw8 = itab_ud(iu16)%iw(1)
     endif

     if (itab_ud(iu26)%iw(1) == iw9) then
        iw21 = itab_ud(iu26)%iw(2)
        im13 = itab_ud(iu26)%im(1) ! to reposition xyzemd(13)
     else
        iw21 = itab_ud(iu26)%iw(1)
        im13 = itab_ud(iu26)%im(2) ! to reposition xyzemd(13)
     endif

     if (im16 == itab_ud(iu46)%im(1)) then
        im22 = itab_ud(iu46)%im(2)
     else
        im22 = itab_ud(iu46)%im(1)
     endif

     if (im18 == itab_ud(iu48)%im(1)) then
        im23 = itab_ud(iu48)%im(2)
     else
        im23 = itab_ud(iu48)%im(1)
     endif

     if (im18 == itab_ud(iu49)%im(1)) then
        im24 = itab_ud(iu49)%im(2)
     else
        im24 = itab_ud(iu49)%im(1)
     endif

     if (im19 == itab_ud(iu51)%im(1)) then
        im25 = itab_ud(iu51)%im(2)
     else
        im25 = itab_ud(iu51)%im(1)
     endif

     if (im21 == itab_ud(iu53)%im(1)) then
        im26 = itab_ud(iu53)%im(2)
     else
        im26 = itab_ud(iu53)%im(1)
     endif

     ! Fill neighbor indices:

     if (itab_ud(iu15)%im(1) == 1) then
        itab_ud(iu15)%im(1) = im18
     else
        itab_ud(iu15)%im(2) = im18
     endif

     if (itab_ud(iu16)%im(1) == 1) then
        itab_ud(iu16)%im(1) = im18
     else
        itab_ud(iu16)%im(2) = im18
     endif

     if (itab_ud(iu25)%im(1) == 1) then
        itab_ud(iu25)%im(1) = im18
     else
        itab_ud(iu25)%im(2) = im18
     endif

     if (itab_ud(iu26)%im(1) == 1) then
        itab_ud(iu26)%im(1) = im18
     else
        itab_ud(iu26)%im(2) = im18
     endif

     if (itab_ud(iu34)%im(1) == im18) then
        itab_ud(iu34)%iw(1) = iw8
        itab_ud(iu34)%iw(2) = iw7
        im5 = itab_ud(iu34)%im(2)
     else
        itab_ud(iu34)%iw(1) = iw7
        itab_ud(iu34)%iw(2) = iw8
        im5 = itab_ud(iu34)%im(1)
     endif

     if (itab_ud(iu35)%im(1) == im19) then
        itab_ud(iu35)%iw(2) = iw19
        itab_ud(iu35)%iw(1) = iw21
        itab_ud(iu35)%im(2) = im18
     else
        itab_ud(iu35)%iw(1) = iw19
        itab_ud(iu35)%iw(2) = iw21
        itab_ud(iu35)%im(1) = im18
     endif

     if (itab_ud(iu41)%im(2) == im17) then
        itab_ud(iu41)%iw(1) = iw27
     else
        itab_ud(iu41)%iw(2) = iw27
     endif

     if (itab_ud(iu42)%im(1) == im17) then
        itab_ud(iu42)%im(2) = im19
        itab_ud(iu42)%iw(1) = iw20
     else
        itab_ud(iu42)%im(1) = im19
        itab_ud(iu42)%iw(2) = iw20
     endif

     if (itab_ud(iu43)%im(2) == im19) then
        itab_ud(iu43)%im(1) = im24
     else
        itab_ud(iu43)%im(2) = im24
     endif

     if (itab_ud(iu44)%im(1) == im19) then
        itab_ud(iu44)%iw(1) = iw29
     else
        itab_ud(iu44)%iw(2) = iw29
     endif

     if (itab_ud(iu45)%im(1) == im20) then
        itab_ud(iu45)%iw(1) = iw31
     else
        itab_ud(iu45)%iw(2) = iw31
     endif

     if (itab_ud(iu48)%iw(2) == iw27) then
        itab_ud(iu48)%im(2) = im17
     else
        itab_ud(iu48)%im(1) = im17
     endif

     if (itab_ud(iu49)%im(2) == im24) then
        itab_ud(iu49)%im(1) = im17
        itab_ud(iu49)%iw(2) = iw20
     else
        itab_ud(iu49)%im(2) = im17
        itab_ud(iu49)%iw(1) = iw20
     endif

     if (itab_ud(iu50)%im(1) == im24) then
        itab_ud(iu50)%im(2) = im20
     else
        itab_ud(iu50)%im(1) = im20
     endif

     if (itab_ud(iu51)%iw(2) == iw31) then
        itab_ud(iu51)%im(1) = im20
     else
        itab_ud(iu51)%im(2) = im20
     endif

     if (itab_wd(iw8)%iu(1) == iu16) then
        itab_wd(iw8)%iu(3) = iu34
     elseif (itab_wd(iw8)%iu(2) == iu16) then
        itab_wd(iw8)%iu(1) = iu34
     else
        itab_wd(iw8)%iu(2) = iu34
     endif

     if (itab_wd(iw19)%iu(1) == iu25) then
        iu33 = itab_wd(iw19)%iu(2)
        itab_wd(iw19)%iu(3) = iu35
     elseif (itab_wd(iw19)%iu(2) == iu25) then
        iu33 = itab_wd(iw19)%iu(3)
        itab_wd(iw19)%iu(1) = iu35
     else
        iu33 = itab_wd(iw19)%iu(1)
        itab_wd(iw19)%iu(2) = iu35
     endif

     ! special case

     if (itab_ud(iu33)%iw(2) == iw19) then
        itab_ud(iu33)%im(2) = im19
     else
        itab_ud(iu33)%im(1) = im19
     endif

     if (itab_wd(iw20)%iu(1) == iu43) then
        itab_wd(iw20)%iu(2) = iu42
        itab_wd(iw20)%iu(3) = iu49
     elseif (itab_wd(iw20)%iu(2) == iu43) then
        itab_wd(iw20)%iu(3) = iu42
        itab_wd(iw20)%iu(1) = iu49
     else
        itab_wd(iw20)%iu(1) = iu42
        itab_wd(iw20)%iu(2) = iu49
     endif

     if (itab_wd(iw27)%iu(1) == iu48) then
        itab_wd(iw27)%iu(2) = iu41
     elseif (itab_wd(iw27)%iu(2) == iu48) then
        itab_wd(iw27)%iu(3) = iu41
     else
        itab_wd(iw27)%iu(1) = iu41
     endif

     if (itab_wd(iw29)%iu(1) == iu50) then
        itab_wd(iw29)%iu(2) = iu44
        itab_wd(iw29)%iu(3) = iu43
     elseif (itab_wd(iw29)%iu(2) == iu50) then
        itab_wd(iw29)%iu(3) = iu44
        itab_wd(iw29)%iu(1) = iu43
     else
        itab_wd(iw29)%iu(1) = iu44
        itab_wd(iw29)%iu(2) = iu43
     endif

     if (itab_wd(iw31)%iu(1) == iu51) then
        itab_wd(iw31)%iu(3) = iu45
     elseif (itab_wd(iw31)%iu(2) == iu51) then
        itab_wd(iw31)%iu(1) = iu45
     else
        itab_wd(iw31)%iu(2) = iu45
     endif

     itab_md(im22)%ngr = ngr
     itab_md(im23)%ngr = ngr
     itab_md(im24)%ngr = ngr
     itab_md(im25)%ngr = ngr
     itab_md(im26)%ngr = ngr

     itab_wd(iw20)%ngr = ngr
     itab_wd(iw26)%ngr = ngr
     itab_wd(iw27)%ngr = ngr
     itab_wd(iw28)%ngr = ngr
     itab_wd(iw29)%ngr = ngr
     itab_wd(iw30)%ngr = ngr
     itab_wd(iw31)%ngr = ngr
     itab_wd(iw32)%ngr = ngr

     ! NEW M locations

     itab_md(im17)%mrlm_orig = itab_md(im18)%mrlm_orig
     itab_md(im20)%mrlm_orig = itab_md(im19)%mrlm_orig
     itab_md(im18)%mrlm_orig = mrlo + 1
     itab_md(im19)%mrlm_orig = mrlo + 1

     xemd(im19) = .5 * (xemd(im24) + xemd(im5))
     yemd(im19) = .5 * (yemd(im24) + yemd(im5))
     zemd(im19) = .5 * (zemd(im24) + zemd(im5))

     xemd(im18) = .5 * (xemd(im19) + xemd(im5))
     yemd(im18) = .5 * (yemd(im19) + yemd(im5))
     zemd(im18) = .5 * (zemd(im19) + zemd(im5))

     xemd(im17) = .75 * xemd(im17) + .25 * xemd(im19)
     yemd(im17) = .75 * yemd(im17) + .25 * yemd(im19)
     zemd(im17) = .75 * zemd(im17) + .25 * zemd(im19)

     xemd(im20) = .75 * xemd(im20) + .25 * xemd(im19)
     yemd(im20) = .75 * yemd(im20) + .25 * yemd(im19)
     zemd(im20) = .75 * zemd(im20) + .25 * zemd(im19)

     xemd(im12) = .833 * xemd(im12) + .167 * xemd(im18)
     yemd(im12) = .833 * yemd(im12) + .167 * yemd(im18)
     zemd(im12) = .833 * zemd(im12) + .167 * zemd(im18)

     xemd(im13) = .833 * xemd(im13) + .167 * xemd(im18)
     yemd(im13) = .833 * yemd(im13) + .167 * yemd(im18)
     zemd(im13) = .833 * zemd(im13) + .167 * zemd(im18)

  enddo ! iper

end subroutine perim_fill3

!===========================================================================

subroutine pent_urows(iurow_pent)

  use mem_delaunay, only: itab_md, itab_ud, itab_wd, nud, nwd
  use mem_grid,     only: impent, mrows

  implicit none

  integer, intent(inout) :: iurow_pent(nud)

  integer :: ipent,im,j,iw,irow,jrow,iw1,iw2,iw3,iwrow,iu

  ! automatic arrays

  integer :: iwrow_temp1(nwd)
  integer :: iwrow_temp2(nwd)

  ! Initialize temporary arrays to zero before loop over pentagon points

  iurow_pent (1:nud) = 0
  iwrow_temp1(1:nwd) = 0
  iwrow_temp2(1:nwd) = 0

  ! Loop over all 12 pentagon M points

  do ipent = 1,12

     im = impent(ipent)

     ! Loop over W points that surround current M point; set row flag to 1

     do j = 1, 5
        iw = itab_md(im)%iw(j)
        iwrow_temp1(iw) = 1
        iwrow_temp2(iw) = 1
     enddo

  enddo

  ! Advance outward and flag each row

  do irow = 1, 2*mrows-1
     jrow = mod(irow,2)

     do iw = 2, nwd

        if (iwrow_temp1(iw) == 0) then

           ! If IW is adjacent to any other IW cell with nonzero mrow,
           ! set mrow for IW cell

           iw1 = itab_wd(iw)%iw(1)
           iw2 = itab_wd(iw)%iw(2)
           iw3 = itab_wd(iw)%iw(3)

           ! Check for positive mrow values

           iwrow = max(iwrow_temp1(iw1), &
                       iwrow_temp1(iw2), &
                       iwrow_temp1(iw3))

           if (iwrow > 0) iwrow_temp2(iw) = iwrow + jrow

        endif

     enddo

     do iw = 2, nwd
        iwrow_temp1(iw) = iwrow_temp2(iw)
     enddo

  enddo

  ! Loop over all U points and flag those between unequal nonzero iwrow values

  do iu = 2, nud
     iw1 = itab_ud(iu)%iw(1)
     iw2 = itab_ud(iu)%iw(2)

     if (iwrow_temp1(iw1) > 0 .and. iwrow_temp1(iw2) > 0 .and.  &
         iwrow_temp1(iw1) /= iwrow_temp1(iw2)) then

         iurow_pent(iu) = min(iwrow_temp1(iw1),iwrow_temp1(iw2))
     endif
  enddo

end subroutine pent_urows

!==============================================================================

subroutine perim_map2(npts, nma, nua, nwa, nper2, imper, iuper, npolyper, &
                      nwdivper, nearpent, ltab_md, ltab_ud, ltab_wd, nest_wd)

  ! Perim_map maps the perimeter points of a nested grid that is being spawned

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, nest_wd_vars

  implicit none

  integer, intent(in)  :: npts  ! Estimated # of perimeter pts (for array dims only)
  integer, intent(in)  :: nma, nua, nwa
  integer, intent(out) :: nper2 ! Actual # of perimeter pts

  integer, intent(out) :: imper(npts) ! Boundary IM index at each perimeter index
  integer, intent(out) :: iuper(npts) ! Boundary IU index at each perimeter index
  integer, intent(out) :: npolyper(npts) ! npoly at each perimeter M pt
  integer, intent(out) :: nwdivper(npts) ! # divided W pts at at each perimeter M pt
  integer, intent(out) :: nearpent(npts) ! flag = 1 if adjacent outside M pt is poly5

  type (itab_md_vars), intent(inout) :: ltab_md(nma)
  type (itab_ud_vars), intent(inout) :: ltab_ud(nua)
  type (itab_wd_vars), intent(inout) :: ltab_wd(nwa)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwa)

  integer :: imstart  ! IM index of starting M point (at a corner of ngr perimeter)
  integer :: im       ! dummy im index
  integer :: nwdiv    ! number of subdivided W pts adjacent to current M pt
  integer :: ima      ! Current M pt in counterclockwise path around perimeter
  integer :: imb      ! Next M pt in counterclockwise path around perimeter
  integer :: iua      ! Current U pt in counterclockwise path around perimeter

  integer :: j, iw, npoly, iu, im1, im2, iw1, iw2

  ! Set IM starting point to zero

  imstart = 0

  ! Loop over all ORIGINAL M points and find the first one that is at
  ! convex corner on boundary of refined area being spawned

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

     if (nwdiv == 2) then
        imstart = im
        exit
     endif

  enddo

  ! Bug check: make sure that imstart is not zero now.

  if (imstart == 0) then
     write(*,*) 'imstart is zero - stopping model'
     stop 'stop imstart perim_map2'
  endif

  ! March around NGR boundary in a counterclockwise direction beginning at IMSTART

  ima = imstart
  nper2 = 1

  do  ! Loop over all original M points on ngr perimeter

     ! Find next M and U points (imb, iua) on perimeter and outside W point (iwout)

     npoly = ltab_md(ima)%npoly
     nwdiv = 0
     nearpent(nper2) = 0

     ! Loop over all W points that are adjacent to ima point

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

     call perim_ngr(nma, nua, nwa, ima, imb, iua, &
                    ltab_md, ltab_ud, ltab_wd, nest_wd)

     imper   (nper2) = ima
     iuper   (nper2) = iua
     npolyper(nper2) = npoly
     nwdivper(nper2) = nwdiv

     ! Check if imb equals istart.  If it does, exit loop

     if (imb == imstart) exit

     nper2 = nper2 + 1
     ima = imb

  enddo

end subroutine perim_map2

!==============================================================================

subroutine perim_ngr(nma, nua, nwa, imstart, imnext, iunext, &
                     ltab_md, ltab_ud, ltab_wd, nest_wd)

  ! Subroutine perim_ngr is to be used during the process of spawning a nested grid,
  ! after temporary arrays nest_ud and nest_wd have been filled for interior points.

  ! Given any M point on the nested grid boundary that existed before the spawn
  ! process began, and proceeding along the boundary in a counterclockwise direction,
  ! find the adjacent U and M points that existed before the spawn process began.

  use mem_delaunay, only: itab_md_vars, itab_ud_vars, itab_wd_vars, nest_wd_vars

  implicit none

  integer, intent(in)  :: nma, nua, nwa
  integer, intent(in) :: imstart  ! starting original M point

  integer, intent(out) :: imnext  ! next original M pt (counterclockwise)
  integer, intent(out) :: iunext  ! next original U pt (counterclockwise)

  type (itab_md_vars), intent(inout) :: ltab_md(nma)
  type (itab_ud_vars), intent(inout) :: ltab_ud(nua)
  type (itab_wd_vars), intent(inout) :: ltab_wd(nwa)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwa)

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
     write(*,*) 'iunext < 2 in subroutine perim_ngr - stopping model'
     stop 'perim_ngr1'
  elseif (imnext < 2) then
     write(*,*) 'imnext < 2 in subroutine perim_ngr - stopping model'
     stop 'perim_ngr2'
  endif

end subroutine perim_ngr

!===========================================================================

  subroutine perim_mrow(nma, nua, nwa, itab_md, itab_ud, itab_wd)

  use mem_delaunay, only: itab_wd_vars, itab_ud_vars, itab_md_vars

  implicit none

  integer, intent(in) :: nma, nua, nwa

  type (itab_md_vars), intent(inout) :: itab_md(nma)
  type (itab_ud_vars), intent(inout) :: itab_ud(nua)
  type (itab_wd_vars), intent(inout) :: itab_wd(nwa)

  integer :: iw, iw1, iw2, iw3, im1, im2, im3
  integer :: irow, jrow, mrow

  integer :: mrow_temp(nwa)

  ! Loop over all W points

  do iw = 2,nwa

     ! Initialize all mrow values to zero.

     itab_wd(iw)%mrow  = 0

     mrow_temp (iw)    = 0

     ! Set mrow values on nested grid border to +/- 1.

     iw1 = itab_wd(iw)%iw(1)
     iw2 = itab_wd(iw)%iw(2)
     iw3 = itab_wd(iw)%iw(3)

     im1 = itab_wd(iw)%im(1)
     im2 = itab_wd(iw)%im(2)
     im3 = itab_wd(iw)%im(3)

     if     (itab_wd(iw)%mrlw < itab_wd(iw1)%mrlw .or. &
             itab_wd(iw)%mrlw < itab_wd(iw2)%mrlw .or. &
             itab_wd(iw)%mrlw < itab_wd(iw3)%mrlw) then

        itab_wd(iw)%mrow = 1
        mrow_temp(iw) = 1

     elseif (itab_wd(iw)%mrlw > itab_wd(iw1)%mrlw .or. &
             itab_wd(iw)%mrlw > itab_wd(iw2)%mrlw .or. &
             itab_wd(iw)%mrlw > itab_wd(iw3)%mrlw) then

        itab_wd(iw)%mrow = -1
        mrow_temp(iw) = -1

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

           ! Check for positive mrow values

           mrow = max(itab_wd(iw1)%mrow &
                     ,itab_wd(iw2)%mrow &
                     ,itab_wd(iw3)%mrow)

           if (mrow > 0) mrow_temp (iw) = mrow + jrow

           ! Check for negative mrow values

           mrow = min(itab_wd(iw1)%mrow &
                     ,itab_wd(iw2)%mrow &
                     ,itab_wd(iw3)%mrow)

           if (mrow < 0) mrow_temp (iw) = mrow - jrow

        endif

     enddo

     do iw = 2,nwa
        itab_wd(iw)%mrow  = mrow_temp(iw)
     enddo

  enddo

end subroutine perim_mrow

!==============================================================================

subroutine ngr_area(ngr, inside, x, y, z, ngrdll, grdrad, grdlat, grdlon)

  ! Subroutine ngr_area checks whether a point located at coordinates (x,y,z) 
  ! is within specified region for NGR refinement

  use max_dims,  only: maxgrds, maxngrdll
  use misc_coms, only: mdomain

  implicit none

  integer, intent(in ) :: ngr
  integer, intent(out) :: inside
  real   , intent(in ) :: x,y,z
  integer, intent(in ) :: ngrdll(maxgrds)
  real,    intent(in ) :: grdrad(maxgrds,maxngrdll)
  real,    intent(in ) :: grdlat(maxgrds,maxngrdll)
  real,    intent(in ) :: grdlon(maxgrds,maxngrdll)

  !-------------------------------------------------------------------------------
  ! New option introduced April 2010:  Define nested grid region as comprising the 
  ! W points adjacent to all M points that are within a specified distance of one
  ! or more connected line segments.
  integer :: ipt, jpt
  real :: seglat, seglon ! lat/lon of segment midpoint
  real :: xs(2),ys(2)    ! PS coordinates of segment endpoints
  real :: xm1, ym1, dist
  real :: t, gradius
  !-------------------------------------------------------------------------------

  inside = 0

  if (ngrdll(ngr) == 1) then

     if (mdomain < 2) then

        ! Transform (x,y,z) location to PS space using the lat/lon of the refinement point

        call e_ps(x,y,z,grdlat(ngr,1),grdlon(ngr,1),xm1,ym1)

        ! If using Cartesian geometry, use direct mapping in Cartesian plane

     else

        xm1 = x - grdlon(ngr,1)
        ym1 = y - grdlat(ngr,1)

     endif

     ! If (x,y,z) location is close enough to line segment, flag it for inclusion
     ! in refined grid interior

     dist = sqrt(xm1 * xm1 + ym1 * ym1)

     if (dist < grdrad(ngr,1)) then
        inside = 1
     elseif (inside == 0 .and. dist < grdrad(ngr,1) * 1.2 ) then
        inside = 2  ! larger radius when searching for ipent points
     endif

  else

     searchloop: do ipt = 1,ngrdll(ngr)-1
        jpt = ipt+1

        ! If using spherical earth geometry (earth-Cartesian coordinates)...

        if (mdomain < 2) then

        ! Transform segment endpoints to PS space using mean lat/lon of each segment

        ! (If multiple segments are used, none should be excessively long in order
        ! to avoid large PS transformation discontinuities at segment endpoints.)

           seglat = .5 * (grdlat(ngr,ipt) + grdlat(ngr,jpt))
           seglon = .5 * (grdlon(ngr,ipt) + grdlon(ngr,jpt))

           ! Correct seglon if segment crosses 180 W

           if (abs(grdlon(ngr,ipt) - grdlon(ngr,jpt)) > 180.) then
              if (seglon <= 0.) then
                 seglon = seglon + 180.
              else
                 seglon = seglon - 180.
              endif
           endif

           call ll_xy (grdlat(ngr,ipt),grdlon(ngr,ipt), &
                seglat,seglon,xs(1),ys(1))

           call ll_xy (grdlat(ngr,jpt),grdlon(ngr,jpt), &
                seglat,seglon,xs(2),ys(2))

           ! Transform (x,y,z) location to PS space using mean lat/lon of each segment

           call e_ps(x,y,z,seglat,seglon,xm1,ym1)

           ! If using Cartesian geometry, use direct mapping in Cartesian plane

        else

           xs(1) = grdlon(ngr,ipt)
           ys(1) = grdlat(ngr,ipt)
           xs(2) = grdlon(ngr,jpt)
           ys(2) = grdlat(ngr,jpt)

           xm1 = x
           ym1 = y

        endif

        ! If (x,y,z) location is close enough to line segment, flag it for inclusion
        ! in refined grid interior

        call linesegdist2(xm1,ym1,xs(1),ys(1),xs(2),ys(2),dist,t)

        gradius = (1.0 - t) * grdrad(ngr,ipt) + t * grdrad(ngr,jpt)

        if (dist < gradius) then
           inside = 1
           exit searchloop
        elseif (inside == 0 .and. dist < gradius * 1.2 ) then
           inside = 2  ! larger radius when searching for ipent points
        endif

     enddo searchloop

  endif

end subroutine ngr_area

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

end function linesegdist

!===========================================================================

subroutine linesegdist2(x0,y0,x1,y1,x2,y2,dist,t)

  ! Determine distance on 2D Cartesian plane between point (x0,y0) and line
  ! segment [(x1,y1),(x2,y2)], and the closest point on the line segment

  implicit none

  real, intent(in ) :: x0,y0,x1,y1,x2,y2
  real, intent(out) :: dist,t
  real              :: dx, dy, xp, yp

  dx = x2 - x1
  dy = y2 - y1

  xp = x0 - x1
  yp = y0 - y1

  t = (xp*dx + yp*dy) / (dx*dx + dy*dy)
  t = max(0., min(1., t))

  dist = sqrt( (xp - t * dx)**2 + (yp - t * dy)**2 )

end subroutine linesegdist2

!==============================================================================

subroutine thirdm(nma, nua, im, jdone, mlist, ltab_md, ltab_ud)

  ! Find set of 5 or 6 M points that are 3 edges away from current IM point along
  ! "straight" paths, i.e., along a path that enters and exits both intervening M
  ! points along opposite edges.  Both intervening M points will always have 6 edges,
  ! given the conditions under which this subroutine is called.

  use mem_delaunay, only: itab_md_vars, itab_ud_vars

  implicit none

  integer, intent(in) :: nma, nua, im
  integer, intent(inout) :: jdone(6,nma)
  integer, intent(out) :: mlist(6)

  type (itab_md_vars), intent(inout) :: ltab_md(nma)
  type (itab_ud_vars), intent(inout) :: ltab_ud(nua)

  integer :: npoly, j, jj, jjop, iu, iuu, iuuu, imm, immm, immmm

  npoly = ltab_md(im)%npoly

  ! Loop over the 5 or 6 edges that connect to current IM point

  do j = 1,npoly

     ! Skip current edge if it has already been traced

     if (jdone(j,im) == 1) cycle

     ! Mark current edge as being traced

     jdone(j,im) = 1

     ! IU is index of current edge

     iu = ltab_md(im)%iu(j)

     ! IMM is M point at far end of IU

     if (im == ltab_ud(iu)%im(1)) then
        imm = ltab_ud(iu)%im(2)
     else
        imm = ltab_ud(iu)%im(1)
     endif

     ! Search over 6 edges that connect to IMM; jjop is opposite edge

     do jj = 1,6
        jjop = mod(jj+2,6) + 1

        ! If current edge is equal to IU, set IUU to opposite edge

        if (ltab_md(imm)%iu(jj) == iu) then
           iuu = ltab_md(imm)%iu(jjop)
           exit
        endif
     enddo

     ! IMMM is M point at far end of IUU

     if (imm == ltab_ud(iuu)%im(1)) then
        immm = ltab_ud(iuu)%im(2)
     else
        immm = ltab_ud(iuu)%im(1)
     endif

     ! Search over 6 edges that connect to IMMM; jjop is opposite edge

     do jj = 1,6
        jjop = mod(jj+2,6) + 1

        ! If current edge is equal to IUU, set IUUU to opposite edge

        if (ltab_md(immm)%iu(jj) == iuu) then
           iuuu = ltab_md(immm)%iu(jjop)
           exit
        endif
     enddo

     ! IMMMM is M point at far end of IUUU

     if (immm == ltab_ud(iuuu)%im(1)) then
        immmm = ltab_ud(iuuu)%im(2)
     else
        immmm = ltab_ud(iuuu)%im(1)
     endif

     mlist(j) = immmm

     ! Search over 6 edges that connect to IMMMM

     do jj = 1,6

        ! If current edge is equal to IUUU, mark it as done

        if (ltab_md(immmm)%iu(jj) == iuuu) then
           jdone(jj,immmm) = 1
           exit
        endif
     enddo

  enddo ! end of j loop

end subroutine thirdm

!==============================================================================

subroutine fill_rad3(nma, nwa, im, ltab_md, ltab_wd, nest_wd)

  use mem_delaunay, only: itab_md_vars, itab_wd_vars, nest_wd_vars

  implicit none

  integer, intent(in) :: nma, nwa, im

  type (itab_md_vars), intent(inout) :: ltab_md(nma)
  type (itab_wd_vars), intent(inout) :: ltab_wd(nwa)
  type (nest_wd_vars), intent(inout) :: nest_wd(nwa)

  integer :: npoly, j, iw, jj, imx, iwx, iwy, im1, im2, im3

  npoly = ltab_md(im)%npoly

  ! Loop over all IW neighbors of current IM point

  do j = 1,npoly
     iw = ltab_md(im)%iw(j)

     ! Flag IW point for NGR refinement

     nest_wd(iw)%iw(3) = 1

     ! Find indices of adjacent M and W points for this IW sector of IM

     if (im == ltab_wd(iw)%im(1)) then
        imx = ltab_wd(iw)%im(2)
        iwx = ltab_wd(iw)%iw(4)
        iwy = ltab_wd(iw)%iw(5)
     elseif (im == ltab_wd(iw)%im(2)) then
        imx = ltab_wd(iw)%im(3)
        iwx = ltab_wd(iw)%iw(6)
        iwy = ltab_wd(iw)%iw(7)
     else
        imx = ltab_wd(iw)%im(1)
        iwx = ltab_wd(iw)%iw(8)
        iwy = ltab_wd(iw)%iw(9)
     endif

     ! Find 3 distant M points for this IW sector of IM

     if (imx == ltab_wd(iwx)%im(1)) then
        im1 = ltab_wd(iwx)%im(2)
        im2 = ltab_wd(iwx)%im(3)
     elseif (imx == ltab_wd(iwx)%im(2)) then
        im1 = ltab_wd(iwx)%im(3)
        im2 = ltab_wd(iwx)%im(1)
     else
        im1 = ltab_wd(iwx)%im(1)
        im2 = ltab_wd(iwx)%im(2)
     endif

     if (im2 == ltab_wd(iwy)%im(1)) then
        im3 = ltab_wd(iwy)%im(2)
     elseif (im2 == ltab_wd(iwy)%im(2)) then
        im3 = ltab_wd(iwy)%im(3)
     else
        im3 = ltab_wd(iwy)%im(1)
     endif

     ! Flag all immediate W neighbors of IM1, IM2, and IM3 (all have npoly = 6)

     do jj = 1,6
        nest_wd(ltab_md(im1)%iw(jj))%iw(3) = 1
        nest_wd(ltab_md(im2)%iw(jj))%iw(3) = 1
        nest_wd(ltab_md(im3)%iw(jj))%iw(3) = 1
     enddo

  enddo ! j,iw

end subroutine fill_rad3
