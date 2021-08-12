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

subroutine icosahedron(nxp0)

  use mem_ijtabs,   only: jtm_grid, jtu_grid, jtv_grid, jtw_grid, &
                          jtm_init, jtu_init, jtv_init, jtw_init, &
                          jtm_prog, jtu_prog, jtv_prog, jtw_prog, &
                          jtm_wadj, jtu_wadj, jtv_wadj, jtw_wadj, &
                          jtm_wstn, jtu_wstn, jtv_wstn, jtw_wstn, &
                          jtm_lbcp, jtu_lbcp, jtv_lbcp, jtw_lbcp, &
                          jtm_vadj, jtu_wall, jtv_wall, jtw_vadj

  use mem_delaunay, only: itab_md, itab_ud, itab_wd, alloc_itabsd, &
                          xemd, yemd, zemd, nmd, nud, nwd

  use mem_grid,     only: impent
  use consts_coms,  only: pi2, erad, erador5

  implicit none

  integer, intent(in) :: nxp0

  real, parameter :: pwrd = 0.9  ! 0.9 is close to making uniform-sized triangles
! real, parameter :: pwrd = 1.0  ! 1.0 is original value

  integer :: ibigd,i,j,idiamond,im_left,iu0,iu1,iu2,iu3,iu4,iw1,iw2,im &
     ,idiamond_top,im_top,im_right,im_bot,nn10,idiamond_right,idiamond_bot &
     ,iu,iw
  integer :: id
 
  real :: wts,wtn,wtw,wte,expansion,anglen,anglew,anglee,angles,wtw0,wte0,sumwt

  integer, parameter :: ibigd_ne(10) = [ 6,7,8,9,10,7,8,9,10,6 ]
  integer, parameter :: ibigd_se(10) = [ 2,3,4,5, 1,2,3,4, 5,1 ]

  real, dimension(10) :: xed_s,xed_n,xed_w,xed_e, &
                         yed_s,yed_n,yed_w,yed_e, &
                         zed_s,zed_n,zed_w,zed_e

  ! Define triangles, edges, and vertices for icosahedral faces and subdivisions

  ! For now, use nxp0 to divide each face

  nn10 = nxp0 * nxp0 * 10

  ! ADD 1 to total number of points needed

  nmd =     nn10 + 2 + 1  ! ADDING 1 for reference point (index = 1)
  nud = 3 * nn10     + 1  ! ADDING 1 for reference point (index = 1)
  nwd = 2 * nn10     + 1  ! ADDING 1 for reference point (index = 1)

  ! Allocate memory for itabs and M earth coords
  ! Initialize all neighbor indices to zero

  call alloc_itabsd(nmd,nud,nwd)

  do im = 2,nmd
     itab_md(im)%imp = im
     itab_md(im)%mrlm = 1
     itab_md(im)%mrlm_orig = 1
     itab_md(im)%ngr = 1
     call mdloopf('f',im, jtm_grid, jtm_init, jtm_prog, jtm_wadj, jtm_wstn, 0)
  enddo

  do iu = 2,nud
     itab_ud(iu)%iup = iu
     call udloopf('f',iu, jtu_grid, jtu_init, jtu_prog, jtu_wadj, jtu_wstn, 0)
  enddo

  do iw = 2,nwd
     itab_wd(iw)%iwp = iw
     call wdloopf('f',iw, jtw_grid, jtw_vadj, 0, 0, 0, 0)
  enddo

  ! Fill big diamond corner coordinates

  do id = 1,5

     anglen = .2 * (id-1) * pi2
     anglew = anglen - .1 * pi2
     anglee = anglen + .1 * pi2

     zed_s(id) = -erad
     xed_s(id) = 0.
     yed_s(id) = 0.

     zed_n(id) = erador5
     xed_n(id) = erador5 * 2. * cos(anglen)
     yed_n(id) = erador5 * 2. * sin(anglen)

     zed_w(id) = -erador5
     xed_w(id) = erador5 * 2. * cos(anglew)
     yed_w(id) = erador5 * 2. * sin(anglew)

     zed_e(id) = -erador5
     xed_e(id) = erador5 * 2. * cos(anglee)
     yed_e(id) = erador5 * 2. * sin(anglee)

  enddo

  do id = 6,10

     angles = .2 * (id-6) * pi2 + .1 * pi2
     anglew = angles - .1 * pi2
     anglee = angles + .1 * pi2

     zed_s(id) = -erador5
     xed_s(id) = erador5 * 2. * cos(angles)
     yed_s(id) = erador5 * 2. * sin(angles)

     zed_n(id) = erad
     xed_n(id) = 0.
     yed_n(id) = 0.

     zed_w(id) = erador5
     xed_w(id) = erador5 * 2. * cos(anglew)
     yed_w(id) = erador5 * 2. * sin(anglew)

     zed_e(id) = erador5
     xed_e(id) = erador5 * 2. * cos(anglee)
     yed_e(id) = erador5 * 2. * sin(anglee)

  enddo

  ! Store IM index of south-pole and north-pole pentagonal points

  impent(1) = 2
  impent(12) = nmd

  do ibigd = 1,10

     do j = 1,nxp0
        do i = 1,nxp0

           idiamond = (ibigd - 1) * nxp0 * nxp0 &
                    + (j - 1)     * nxp0        &
                    +  i

  ! Indices that are "attached" to this diamond

           im_left  = idiamond + 2

           iu0 = 3 * idiamond
           iu1 = 3 * idiamond - 1
           iu3 = 3 * idiamond + 1

           iw1 = 2 * idiamond
           iw2 = 2 * idiamond + 1

  ! Store IM index of 10 out of 12 pentagonal points

           if (i == 1 .and. j == nxp0) impent(ibigd+1) = im_left

  ! Indices that are "attached" to another diamond

           if (ibigd < 6) then   ! Southern 5 diamonds

              ! Top diamond indices

              if (i < nxp0) then
                 idiamond_top = idiamond + 1
              else
                 idiamond_top = (ibigd_ne(ibigd) - 1) * nxp0 * nxp0 &
                              + (j - 1)               * nxp0        &
                              +  1
              endif

              im_top   = idiamond_top + 2
              iu4 = 3 * idiamond_top - 1   ! (it's the iu1 for id_top)

              ! Right diamond indices

              if (j > 1 .and. i < nxp0) then
                 idiamond_right = idiamond - nxp0 + 1
              elseif (j == 1) then
                 idiamond_right = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                                + (i - 1)             * nxp0        &
                                +  1
                 iu2 = 3 * idiamond_right - 1 ! (it's the iu1 for id_right)
              else            ! case for i = nxp0 and j > 1
                 idiamond_right = (ibigd_ne(ibigd)-1) * nxp0 * nxp0 &
                                + (j - 2)             * nxp0        &
                                +  1
              endif

              im_right = idiamond_right + 2

              ! Bottom diamond indices

              if (j > 1) then
                 idiamond_bot = idiamond - nxp0
                 iu2 = 3 * idiamond_bot + 1 ! (it's the iu3 for id_bot)
              else
                 idiamond_bot = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                              + (i - 2)             * nxp0        &
                              +  1
              endif
              im_bot = idiamond_bot + 2

              if (i == 1 .and. j == 1) im_bot = 2

           else                  ! Northern 5 diamonds

             ! Top diamond indices

              if (i < nxp0) then
                 idiamond_top = idiamond + 1
                 iu4 = 3 * idiamond_top - 1   ! (it's the iu1 for id_top)
              else
                 idiamond_top = (ibigd_ne(ibigd)-1) * nxp0 * nxp0 &
                              + (nxp0 - 1)          * nxp0        &
                              +  j + 1
              endif

              im_top = idiamond_top + 2

              ! Right diamond indices

              if (j > 1 .and. i < nxp0) then
                 idiamond_right = idiamond - nxp0 + 1
              elseif (j == 1 .and. i < nxp0) then
                 idiamond_right = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                                + (nxp0 - 1)          * nxp0        &
                                + i + 1
              else            ! case for i = nxp0
                 idiamond_right = (ibigd_ne(ibigd)-1) * nxp0 * nxp0 &
                                + (nxp0 - 1)          * nxp0        &
                                +  j
                 iu4 = 3 * idiamond_right + 1 ! (it's the iu3 for id_right)
              endif

              im_right = idiamond_right + 2

              ! Bottom diamond indices

              if (j > 1) then
                 idiamond_bot = idiamond - nxp0
              else
                 idiamond_bot = (ibigd_se(ibigd)-1) * nxp0 * nxp0 &
                              + (nxp0 - 1)          * nxp0        &
                              + i
              endif

              im_bot = idiamond_bot + 2
              iu2 = 3 * idiamond_bot + 1 ! (it's the iu3 for id_bot)

              if (i == nxp0 .and. j == nxp0) &
                 im_top = 10 * nxp0 * nxp0 + 3

           endif

           call fill_diamond(im_left,im_right,im_top,im_bot &
              ,iu0,iu1,iu2,iu3,iu4,iw1,iw2)

           ! M point (xemd,yemd,zemd) coordinates

           if (i + j <= nxp0) then
              wts  = max(0.,min(1.,real(nxp0 + 1 - i - j) / real(nxp0)))
              wtn  = 0.
              wtw0 = max(0.,min(1.,real(j) / real(i + j - 1)))
              wte0 = 1. - wtw0
           else
              wts  = 0.
              wtn  = max(0.,min(1.,real(i + j - nxp0 - 1) / real(nxp0)))
              wte0 = max(0.,min(1.,real(nxp0 - j) &
                   / real(2 * nxp0 + 1 - i - j)))
              wtw0 = 1. - wte0
           endif

           ! Experimental adjustment in spacing
           ! Compute sum of weights raised to pwrd

           wtw = (1. - wts - wtn) * wtw0
           wte = (1. - wts - wtn) * wte0
           sumwt = wts**pwrd + wtn**pwrd + wtw**pwrd + wte**pwrd

           wts = wts**pwrd / sumwt
           wtn = wtn**pwrd / sumwt
           wtw = wtw**pwrd / sumwt
           wte = wte**pwrd / sumwt

           xemd(im_left) = wts * xed_s(ibigd) &
                         + wtn * xed_n(ibigd) &
                         + wtw * xed_w(ibigd) &
                         + wte * xed_e(ibigd)

           yemd(im_left) = wts * yed_s(ibigd) &
                         + wtn * yed_n(ibigd) &
                         + wtw * yed_w(ibigd) &
                         + wte * yed_e(ibigd)

           zemd(im_left) = wts * zed_s(ibigd) &
                         + wtn * zed_n(ibigd) &
                         + wtw * zed_w(ibigd) &
                         + wte * zed_e(ibigd)

           ! Push M point coordinates out to earth radius

           expansion = erad / sqrt( xemd(im_left) ** 2 &
                                  + yemd(im_left) ** 2 &
                                  + zemd(im_left) ** 2 )

           xemd(im_left) = xemd(im_left) * expansion
           yemd(im_left) = yemd(im_left) * expansion
           zemd(im_left) = zemd(im_left) * expansion

        enddo  ! end i loop
     enddo     ! end j loop

  enddo        ! end idbig loop

  xemd(2) = 0.
  yemd(2) = 0.
  zemd(2) = -erad

  xemd(nmd) = 0.
  yemd(nmd) = 0.
  zemd(nmd) = erad

  ! call twist(nxp0)

  call tri_neighbors(nmd, nud, nwd, itab_md, itab_ud, itab_wd)

  ! This is the place to do spring dynamics

  call spring_dynamics(3, 1, 1, nxp0, nmd, nud, nwd, xemd, yemd, zemd, &
                       itab_md, itab_ud, itab_wd)

end subroutine icosahedron

!===============================================================================

subroutine fill_diamond(im_left,im_right,im_top,im_bot,  &
                        iu0,iu1,iu2,iu3,iu4,iw1,iw2)

  use mem_delaunay, only: itab_ud, itab_wd

  implicit none

  integer, intent(in) :: im_left,im_right,im_top,im_bot
  integer, intent(in) :: iu0,iu1,iu2,iu3,iu4,iw1,iw2

  itab_ud(iu0)%im(1) = im_left
  itab_ud(iu0)%im(2) = im_right
  itab_ud(iu0)%iw(1) = iw1
  itab_ud(iu0)%iw(2) = iw2
  itab_ud(iu0)%mrlu = 1

  itab_ud(iu1)%im(1) = im_left
  itab_ud(iu1)%im(2) = im_bot
  itab_ud(iu1)%iw(2) = iw1
  itab_ud(iu1)%mrlu = 1

  itab_ud(iu2)%iw(1) = iw1

  itab_ud(iu3)%im(1) = im_top
  itab_ud(iu3)%im(2) = im_left
  itab_ud(iu3)%iw(2) = iw2
  itab_ud(iu3)%mrlu = 1

  itab_ud(iu4)%iw(1) = iw2

  itab_wd(iw1)%iu(1) = iu0
  itab_wd(iw1)%iu(2) = iu1
  itab_wd(iw1)%iu(3) = iu2
  itab_wd(iw1)%mrlw = 1
  itab_wd(iw1)%mrlw_orig = 1
  itab_wd(iw1)%ngr = 1

  itab_wd(iw2)%iu(1) = iu0
  itab_wd(iw2)%iu(2) = iu4
  itab_wd(iw2)%iu(3) = iu3
  itab_wd(iw2)%mrlw = 1
  itab_wd(iw2)%mrlw_orig = 1
  itab_wd(iw2)%ngr = 1

end subroutine fill_diamond

!===============================================================================

subroutine twist(nxp)

  use mem_delaunay, only: itab_ud, itab_wd, xemd, yemd, zemd, nmd, nud
  use consts_coms,  only: pi1, pi2

  implicit none

  integer, intent(in) :: nxp

  integer :: im,im1,im2,im1_east,im2_east
  integer :: iu,iu1,iu2,iu3,iu_east,iueq,iueq_east,iueq_fill
  integer :: iw

  integer :: nxpo2  ! half of nxp (nxp must be even if twist is done)

  integer :: iueq_tab(10*nxp)
  integer :: iweq_tab(10*nxp)

  real :: radm  ! Distance of M point from Earth axis [m]
  real :: angm  ! longitude (radians) of M point before shift is applied

  nxpo2 = nxp / 2

  ! Initialize iueq_tab and iweq_tab to zero

  iueq_tab(:) = 0
  iweq_tab(:) = 0

  ! Shift all M points in southern hemisphere (excluding equator)
  ! 36 degrees to the east

  do im = 2,nmd
     if (zemd(im) < -100.) then

        radm = sqrt(xemd(im) ** 2 + yemd(im) ** 2)
        angm = atan2(yemd(im),xemd(im))

        xemd(im) = radm * cos(angm + .2 * pi1)
        yemd(im) = radm * sin(angm + .2 * pi1)

     endif
  enddo

  ! Loop over all U points and search for those that are on the equator

  do iu = 1,nud
     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

     if (abs(zemd(im1)) < 100. .and. abs(zemd(im2)) < 100.) then

     ! Compute special index for equatorial U points that increases with longitude.
     ! This uses knowledge of icosahedron subroutine that im2 is always east of iu.
     ! IUEQ skips approximately every other integer and therefore needs to be collapsed.

        iueq = int(10 * nxp * (atan2(yemd(im2),xemd(im2)) + pi1) / pi2) + 1

        ! Fill table of iu and iw1 indices for each iueq

        iueq_tab(iueq) = iu
        iweq_tab(iueq) = itab_ud(iu)%iw(1)

     endif
  enddo

  ! Collapse IUEQ tables in order to use consecutive integer indices

  iueq_fill = 0

  do iueq = 1,nxp * 10

     if (iueq_tab(iueq) > 1) then
        iueq_fill = iueq_fill + 1

        iueq_tab(iueq_fill) = iueq_tab(iueq)
        iweq_tab(iueq_fill) = iweq_tab(iueq)
     endif

  enddo

  ! Loop over equatorial U points

  do iueq = 1,nxp * 5

  ! Find equatorial index of U point that is 36 degrees to the east

     iueq_east = iueq + nxpo2
     if (iueq_east > nxp * 5) iueq_east = iueq_east - nxp * 5

  ! Get U point indices for iueq and iueq_east points

     iu      = iueq_tab(iueq)
     iu_east = iueq_tab(iueq_east)

  ! Get index for adjacent W point to the south of iu

     iw = iweq_tab(iueq)

  ! Get indices for M endpoints of iu and iu_east

     im1 = itab_ud(iu)%im(1)
     im2 = itab_ud(iu)%im(2)

     im1_east = itab_ud(iu_east)%im(1)
     im2_east = itab_ud(iu_east)%im(2)

  ! Get U neighbor indices for adjacent W point to the south

     iu1 = itab_wd(iw)%iu(1)
     iu2 = itab_wd(iw)%iu(2)
     iu3 = itab_wd(iw)%iu(3)

  ! For U points on equator: shift their south W point neighbor index

     itab_ud(iu_east)%iw(1) = iw

  ! For W points bordering equator on south: shift their equatorial U point
  ! neighbor index

     if     (iu1 == iu) then
        itab_wd(iw)%iu(1) = iu_east
     elseif (iu2 == iu) then
        itab_wd(iw)%iu(2) = iu_east
     elseif (iu3 == iu) then
        itab_wd(iw)%iu(3) = iu_east
     endif

  ! For U points bordering equator on south: shift their equatorial M point
  ! neighbor index

     if (iu1 /= iu) then
        if     (itab_ud(iu1)%im(1) == im1) then
                itab_ud(iu1)%im(1) =  im1_east
        elseif (itab_ud(iu1)%im(2) == im1) then
                itab_ud(iu1)%im(2) =  im1_east
        elseif (itab_ud(iu1)%im(1) == im2) then
                itab_ud(iu1)%im(1) =  im2_east
        elseif (itab_ud(iu1)%im(2) == im2) then
                itab_ud(iu1)%im(2) =  im2_east
        endif
     endif

     if (iu2 /= iu) then
        if     (itab_ud(iu2)%im(1) == im1) then
                itab_ud(iu2)%im(1) =  im1_east
        elseif (itab_ud(iu2)%im(2) == im1) then
                itab_ud(iu2)%im(2) =  im1_east
        elseif (itab_ud(iu2)%im(1) == im2) then
                itab_ud(iu2)%im(1) =  im2_east
        elseif (itab_ud(iu2)%im(2) == im2) then
                itab_ud(iu2)%im(2) =  im2_east
        endif
     endif

     if (iu3 /= iu) then
        if     (itab_ud(iu3)%im(1) == im1) then
                itab_ud(iu3)%im(1) =  im1_east
        elseif (itab_ud(iu3)%im(2) == im1) then
                itab_ud(iu3)%im(2) =  im1_east
        elseif (itab_ud(iu3)%im(1) == im2) then
                itab_ud(iu3)%im(1) =  im2_east
        elseif (itab_ud(iu3)%im(2) == im2) then
                itab_ud(iu3)%im(2) =  im2_east
        endif
     endif

  enddo

end subroutine twist
