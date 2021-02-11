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
subroutine voronoi_sfc()

  use mem_sfcg,    only: nsfcgrids, nsfcgrid_root, nmd, nud, nwd, &
                         itab_md, itab_ud, itab_wd, xemd, yemd, zemd, &
                         nmsfc, nvsfc, nwsfc, &
                         itab_msfc, itab_vsfc, itab_wsfc, sfcg, alloc_sfcgrid1

  use misc_coms,   only: io6, mdomain
  use consts_coms, only: erad, piu180
  use oname_coms,  only: nl

  implicit none

  integer :: iw1, iw2, iw3, im, iv, iw
  integer :: imd,iud,iwd,j,iwn
  real    :: expansion, raxis
  real    :: xebc,yebc,zebc
  real    :: glatbc,glonbc
  real    :: x1,x2,x3,y1,y2,y3
  real    :: dx12,dx13,dx23
  real    :: s1,s2,s3
  real    :: xcc,ycc

  ! Interchange grid dimensions

  nmsfc = nwd
  nvsfc = nud
  nwsfc = nmd

  ! Allocate Voronoi set of arrays

  call alloc_sfcgrid1(nmsfc, nvsfc, nwsfc)

  ! Transfer information from Delaunay to Voronoi arrays

  do iw = 1,nwsfc
     imd = iw

     sfcg%xew(iw) = xemd(imd)
     sfcg%yew(iw) = yemd(imd)
     sfcg%zew(iw) = zemd(imd)
  enddo

  ! Deallocate xemd, yemd, zemd

  deallocate (xemd,yemd,zemd)

  ! Compute XEM,YEM,ZEM location as circumcentric coordinates of 3 W points.
  ! This establishes W cell as voronoi.

  do im = 2,nmsfc
     iwd = im

     if (any(itab_wd(iwd)%im(1:3) < 2)) cycle

     ! Indices of 3 M points surrounding WD point

     iw1 = itab_wd(iwd)%im(1)
     iw2 = itab_wd(iwd)%im(2)
     iw3 = itab_wd(iwd)%im(3)

     ! First, compute barycenter of 3 W points

     xebc = (sfcg%xew(iw1) + sfcg%xew(iw2) + sfcg%xew(iw3)) / 3.
     yebc = (sfcg%yew(iw1) + sfcg%yew(iw2) + sfcg%yew(iw3)) / 3.
     zebc = (sfcg%zew(iw1) + sfcg%zew(iw2) + sfcg%zew(iw3)) / 3.

     if (mdomain <= 1) then

        ! If mdomain <= 1, push M point coordinates out to earth radius

        expansion = erad / sqrt(xebc ** 2  &
                              + yebc ** 2  &
                              + zebc ** 2  )

         xebc = xebc * expansion
         yebc = yebc * expansion
         zebc = zebc * expansion

        ! Get latitude and longitude of barycentric point

        raxis = sqrt(xebc ** 2 + yebc ** 2)

        glatbc = atan2(zebc,raxis) * piu180
        glonbc = atan2(yebc,xebc) * piu180

        ! Transform 3 W points to PS coordinates

        call e_ps(sfcg%xew(iw1),sfcg%yew(iw1),sfcg%zew(iw1),glatbc,glonbc,x1,y1)
        call e_ps(sfcg%xew(iw2),sfcg%yew(iw2),sfcg%zew(iw2),glatbc,glonbc,x2,y2)
        call e_ps(sfcg%xew(iw3),sfcg%yew(iw3),sfcg%zew(iw3),glatbc,glonbc,x3,y3)

     else

        ! For Cartesian domain, use given planar X,Y coordinates

        x1 = sfcg%xew(iw1) - xebc
        x2 = sfcg%xew(iw2) - xebc
        x3 = sfcg%xew(iw3) - xebc

        y1 = sfcg%yew(iw1) - yebc
        y2 = sfcg%yew(iw2) - yebc
        y3 = sfcg%yew(iw3) - yebc

     endif

     ! Compute intermediate quanties

     dx12 = x2 - x1
     dx13 = x3 - x1
     dx23 = x3 - x2

     s1 = x1**2 + y1**2
     s2 = x2**2 + y2**2
     s3 = x3**2 + y3**2

     ! Algebraic solution for circumcenter Y coordinate

     ycc = .5 * (dx13 * s2 - dx12 * s3 - dx23 * s1) &
              / (dx13 * y2 - dx12 * y3 - dx23 * y1)

     ! Algebraic solution for circumcenter X coordinate

     if (abs(dx12) > abs(dx13)) then
        xcc = (s2 - s1 - ycc * 2. * (y2 - y1)) / (2. * dx12)
     else
        xcc = (s3 - s1 - ycc * 2. * (y3 - y1)) / (2. * dx13)
     endif

     ! For global domain, transform circumcenter from PS to earth coordinates

     if (mdomain <= 1) then
        call ps_e(sfcg%xem(im),sfcg%yem(im),sfcg%zem(im),glatbc,glonbc,xcc,ycc)
     else
        sfcg%xem(im) = xcc + xebc
        sfcg%yem(im) = ycc + yebc
     endif

  enddo

  ! Loop over V points

  do iv = 2,nvsfc
     iud = iv

     itab_vsfc(iv)%imn(1:2)  = itab_ud(iud)%iw(1:2)
     itab_vsfc(iv)%iwn(1:2)  = itab_ud(iud)%im(1:2)
  enddo

  ! Loop over WSFC points

  do iw = 2,nwsfc
     imd = iw

     if (nsfcgrid_root < 0 .or. &
        (nsfcgrids > 0 .and. itab_md(imd)%ngr > abs(nsfcgrid_root))) then

        ! Any surface grid cells that are refined independently of the
        ! atmospheric grid are flagged to remain voronoi cells (i.e., to not be
        ! subsequently subdivided according to atmospheric grid cut cells in
        ! subroutine makesfc3.  If nsfcgrid_root < 0, then ALL surface grid
        ! are flagged to remain voronoi cells.  These cells will not
        ! necessarily coincide exactly with atmospheric grid columns and may
        ! overlap more than one atmospheric grid column.  A value of
        ! ivoronoi = 3 is assigned here to represent this status.
        ! (May want to revise conditions of IF STATEMENT in the future.)

        itab_wsfc(iw)%ivoronoi = 3

     else

        ! All other cells (for which ivoronoi = 0) should exactly match an 
        ! atmospheric grid column (until subdivision, if any, is done in
        ! subroutine makesfc3).  Accordingly, the surface grid cell's 
        ! atmosphere column index is assigned here.  (A check for exact
        ! coincidence is done in subroutine makesfc3.)

        itab_wsfc(iw)%ivoronoi = 0
        itab_wsfc(iw)%nwatm    = 1
        itab_wsfc(iw)%iwatm(1) = iw

     endif

     itab_wsfc(iw)%npoly   = itab_md(imd)%npoly

     ! Loop over IM/IV neighbors of IW

     do j = 1,itab_wsfc(iw)%npoly
        im = itab_md(imd)%iw(j)
        iwd = im
        iv = itab_md(imd)%iu(j)

        iw1 = itab_vsfc(iv)%iwn(1)
        iw2 = itab_vsfc(iv)%iwn(2)

        itab_wsfc(iw)%imn(j) = im
        itab_wsfc(iw)%ivn(j) = iv

        if (iw1 == iw) then
           itab_wsfc(iw)%iwn(j) = iw2
        else
           itab_wsfc(iw)%iwn(j) = iw1
        endif
     enddo

  enddo

  ! Loop over MSFC points

  do im = 2,nmsfc
     iwd = im

     itab_msfc(im)%ivn(1:3)   = itab_wd(iwd)%iu(1:3)
     itab_msfc(im)%iwn(1:3)   = itab_wd(iwd)%im(1:3)
  enddo

  ! Loop over all WSFC points and find those that are Voronoi cells

  do iw = 2,nwsfc
     if (itab_wsfc(iw)%ivoronoi == 3) then

        ! If any IWN neighbor of this Voronoi SFC cell is not itself a Voronoi
        ! cell, set its ivoronoi value to 2.  This flags the MAKESFC process
        ! (1) to avoid the topocut procedure of cutting those cells by
        ! levels of the atmosphere grid, since some of its M vertices may have
        ! moved in spring dynamics, and (2) to avoid filling the IVN and IWN
        ! members of itab_wsfc, since they may not all be present.

        do j = 1,itab_wsfc(iw)%npoly
           iwn = itab_wsfc(iw)%iwn(j)
           if (itab_wsfc(iwn)%ivoronoi < 2) itab_wsfc(iwn)%ivoronoi = 2
        enddo
     endif
  enddo

  deallocate(itab_md,itab_ud,itab_wd)

end subroutine voronoi_sfc

!===============================================================================

subroutine grid_geometry_hex_sfc()

  use mem_sfcg,   only: nmsfc, nvsfc, nwsfc, itab_vsfc, itab_wsfc, sfcg
  use misc_coms,   only: io6, mdomain, nxp
  use consts_coms, only: erad, erad2, piu180
  use oplot_coms,  only: op
  use oname_coms,  only: nl
  use mem_para,    only: myrank

  implicit none

  integer       :: im,iv,iw,ivn
  integer       :: im1,im2
  integer       :: iw1,iw2
  integer       :: j,npoly
  real          :: raxis, expansion
  real          :: xm1,xm2,ym1,ym2
  real          :: xw1,xw2,yw1,yw2
  real          :: xq1, yq1, xq2, yq2, psiz, vsprd
  integer       :: iskip
  character(10) :: string

  ! Loop over all M points and compute their latitude and longitude

  !$omp parallel
  !$omp do private(raxis)
  do im = 2,nmsfc
     if (mdomain <= 1) then
        raxis = sqrt(sfcg%xem(im) ** 2 + sfcg%yem(im) ** 2)  ! dist from earth axis
        sfcg%glatm(im) = atan2(sfcg%zem(im),raxis)        * piu180
        sfcg%glonm(im) = atan2(sfcg%yem(im),sfcg%xem(im)) * piu180
     else
        sfcg%glatm(im) = 0.  ! want it this way?
        sfcg%glonm(im) = 0.  ! want it this way?
     endif
  enddo
  !$omp end do

  ! Loop over all V points

  !$omp do private(im1,im2,iw1,iw2,expansion)
  do iv = 2,nvsfc

     ! M-point indices of two end points of V segment

     im1 = itab_vsfc(iv)%imn(1)
     im2 = itab_vsfc(iv)%imn(2)

     ! W-point indices on either side of V segment

     iw1 = itab_vsfc(iv)%iwn(1)
     iw2 = itab_vsfc(iv)%iwn(2)

     ! V point is midway between W points of Voronoi cells

     sfcg%xev(iv) = .5 * (sfcg%xew(iw1) + sfcg%xew(iw2))
     sfcg%yev(iv) = .5 * (sfcg%yew(iw1) + sfcg%yew(iw2))
     sfcg%zev(iv) = .5 * (sfcg%zew(iw1) + sfcg%zew(iw2))

     ! If mdomain <= 1, push V point coordinates out to earth radius

     if (mdomain <= 1) then

        expansion = erad / sqrt(sfcg%xev(iv) ** 2 &
                              + sfcg%yev(iv) ** 2 &
                              + sfcg%zev(iv) ** 2 )

        sfcg%xev(iv) = sfcg%xev(iv) * expansion
        sfcg%yev(iv) = sfcg%yev(iv) * expansion
        sfcg%zev(iv) = sfcg%zev(iv) * expansion

     endif

     ! Normal distance across U face

     sfcg%dnu(iv) = sqrt( (sfcg%xem(im1) - sfcg%xem(im2))**2 &
                       + (sfcg%yem(im1) - sfcg%yem(im2))**2 &
                       + (sfcg%zem(im1) - sfcg%zem(im2))**2 )
   ! sfcg%dnu(iv) = erad2 * asin(sfcg%dnu(iv) / erad2) ! geodesic arc length
     sfcg%dniu(iv) = 1. / sfcg%dnu(iv)

     ! Normal distance across V face

     sfcg%dnv(iv) = sqrt( (sfcg%xew(iw1) - sfcg%xew(iw2))**2 &
                       + (sfcg%yew(iw1) - sfcg%yew(iw2))**2 &
                       + (sfcg%zew(iw1) - sfcg%zew(iw2))**2 )
   ! sfcg%dnv(iv) = erad2 * asin(sfcg%dnv(iv) / erad2) ! geodesic arc length
     sfcg%dniv(iv) = 1. / sfcg%dnv(iv)

  enddo
  !$omp end do
  !$omp end parallel

  !$omp parallel do private(raxis,npoly,j,ivn)
  do iw = 2,nwsfc

     ! Fill outward unit vector components and latitude and longitude of W point

     if (mdomain <= 1) then

        raxis = sqrt(sfcg%xew(iw) ** 2 + sfcg%yew(iw) ** 2)


        sfcg%glatw(iw) = atan2(sfcg%zew(iw),raxis)        * piu180
        sfcg%glonw(iw) = atan2(sfcg%yew(iw),sfcg%xew(iw)) * piu180

     else

        sfcg%glatw(iw) = 0. ! want it this way?
        sfcg%glonw(iw) = 0. ! want it this way?

     endif

     ! Number of polygon edges/vertices

     npoly = itab_wsfc(iw)%npoly
     sfcg%area(iw) = 0.

     ! Loop over all polygon edges

     do j = 1,npoly
        ivn = itab_wsfc(iw)%ivn(j)
        sfcg%area(iw) = sfcg%area(iw) + 0.25 * sfcg%dnv(ivn) * sfcg%dnu(ivn)
     enddo

  enddo
  !$omp end parallel do

  ! Plot grid lines

  if (.false.) then

     if (myrank /= 0) return

     call o_reopnwk()
     call plotback()
     call oplot_set(1)
     psiz = .035 / real(nxp) ! not good with nested grids
     vsprd = .10 * sqrt(sfcg%area(nwsfc))

     do iv = 2, nvsfc

        im1 = itab_vsfc(iv)%imn(1)
        im2 = itab_vsfc(iv)%imn(2)

        iw1 = itab_vsfc(iv)%iwn(1)
        iw2 = itab_vsfc(iv)%iwn(2)

        call oplot_transform(1,sfcg%xem(im1),sfcg%yem(im1),sfcg%zem(im1),xm1,ym1)
        call oplot_transform(1,sfcg%xem(im2),sfcg%yem(im2),sfcg%zem(im2),xm2,ym2)
        call oplot_transform(1,sfcg%xew(iw2),sfcg%yew(iw2),sfcg%zew(iw2),xw2,yw2)
        call oplot_transform(1,sfcg%xew(iw1),sfcg%yew(iw1),sfcg%zew(iw1),xw1,yw1)

        call trunc_segment(xm1,xm2,ym1,ym2,xq1,xq2,yq1,yq2,iskip)

        if (iskip == 1) cycle

        call o_frstpt (xq1,yq1)
        call o_vector (xq2,yq2)

        if ( xm1 < op%xmin .or.  &
             xm1 > op%xmax .or.  &
             ym1 < op%ymin .or.  &
             ym1 > op%ymax ) cycle

        write(string,'(I0)') im1
        call o_plchlq (xm1,ym1,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') im2
        call o_plchlq (xm2,ym2,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iw1
        call o_plchlq (xw1,yw1,trim(adjustl(string)),psiz,0.,0.)

        write(string,'(I0)') iw2
        call o_plchlq (xw2,yw2,trim(adjustl(string)),psiz,0.,0.)

     enddo  ! IV

     call o_frame()
     call o_clswk()

  endif

end subroutine grid_geometry_hex_sfc

!===============================================================================

subroutine sfc_atm_hex_overlay(iwsfc)

  use mem_grid,   only: nza, nwa, arw0, xem, yem, zem, xew, yew, zew, &
                        zm, topm, topw, arw, volt, nma, dzt, &
                        glatm, glonm, glatw, glonw
  use mem_ijtabs,  only: itab_w, mrls
  use mem_sfcg,   only: itab_wsfc, sfcg
  use consts_coms, only: r8

  implicit none

  integer, parameter :: npmax = 7  ! Heptagons are max polygon for SINGLE GRID LEVEL
                                   ! in atm polygon cell 

  integer, parameter :: nqmax = 7  ! Land cells can also be up to heptagons at this stage

  integer, intent(in) :: iwsfc

  integer :: iw, npoly, nsfcpoly, jmsfc, imsfc, jm, im, j, nwatm, idum

  real :: xp0, yp0, xq0, yq0, dum

  real(r8) :: xw, yw, alpha
  real(r8) :: xp(npmax),yp(npmax)
  real(r8) :: xq(nqmax),yq(nqmax)

  real :: area

  real :: xesfcmin, yesfcmin, zesfcmin
  real :: xesfcmax, yesfcmax, zesfcmax

  real, save, allocatable :: xeamin(:),yeamin(:),zeamin(:)
  real, save, allocatable :: xeamax(:),yeamax(:),zeamax(:)

  logical, save :: firstcall = .true.

  if (firstcall) then
     firstcall = .false.

     allocate (xeamin(nwa)); xeamin = 1.e9
     allocate (yeamin(nwa)); yeamin = 1.e9
     allocate (zeamin(nwa)); zeamin = 1.e9

     allocate (xeamax(nwa)); xeamax = -1.e9
     allocate (yeamax(nwa)); yeamax = -1.e9
     allocate (zeamax(nwa)); zeamax = -1.e9

     do iw = 2,nwa
        npoly = itab_w(iw)%npoly

        do jm = 1,npoly
           im = itab_w(iw)%im(jm)

           if (xeamin(iw) > xem(im)) xeamin(iw) = xem(im)
           if (yeamin(iw) > yem(im)) yeamin(iw) = yem(im)
           if (zeamin(iw) > zem(im)) zeamin(iw) = zem(im)

           if (xeamax(iw) < xem(im)) xeamax(iw) = xem(im)
           if (yeamax(iw) < yem(im)) yeamax(iw) = yem(im)
           if (zeamax(iw) < zem(im)) zeamax(iw) = zem(im)
        enddo
     enddo

  endif

  xesfcmin = 1.e9
  yesfcmin = 1.e9
  zesfcmin = 1.e9

  xesfcmax = -1.e9
  yesfcmax = -1.e9
  zesfcmax = -1.e9

  nwatm = 0
  nsfcpoly = itab_wsfc(iwsfc)%npoly

  do jmsfc = 1,nsfcpoly
     imsfc = itab_wsfc(iwsfc)%imn(jmsfc)

     if (xesfcmin > sfcg%xem(imsfc)) xesfcmin = sfcg%xem(imsfc)
     if (yesfcmin > sfcg%yem(imsfc)) yesfcmin = sfcg%yem(imsfc)
     if (zesfcmin > sfcg%zem(imsfc)) zesfcmin = sfcg%zem(imsfc)

     if (xesfcmax < sfcg%xem(imsfc)) xesfcmax = sfcg%xem(imsfc)
     if (yesfcmax < sfcg%yem(imsfc)) yesfcmax = sfcg%yem(imsfc)
     if (zesfcmax < sfcg%zem(imsfc)) zesfcmax = sfcg%zem(imsfc)
  enddo

  ! Loop over all neighbor M points of this iwsfc

  do jmsfc = 1,nsfcpoly
     imsfc = itab_wsfc(iwsfc)%imn(jmsfc)

     ! Evaluate x,y coordinates of LAND cell M points on gnomic
     ! plane tangent at iwsfc

     call e_gn(sfcg%xem(imsfc), sfcg%yem(imsfc), sfcg%zem(imsfc), &
               sfcg%glatw(iwsfc), sfcg%glonw(iwsfc), xq0, yq0)

     xq(jmsfc) = real(xq0,r8)
     yq(jmsfc) = real(yq0,r8)
  enddo

  ! Loop over all ATM IW cells

  do iw = 1,nwa

     ! Skip interaction using non-overlap check

     if (abs(xew(iw) + sfcg%xew(iwsfc)) < 12.e6) then  ! Skip if both near X pole
        if (xeamin(iw) > xesfcmax + 10.) cycle
        if (xesfcmin > xeamax(iw) + 10.) cycle
     endif

     if (abs(yew(iw) + sfcg%yew(iwsfc)) < 12.e6) then  ! Skip if both near Y pole
        if (yeamin(iw) > yesfcmax + 10.) cycle
        if (yesfcmin > yeamax(iw) + 10.) cycle
     endif

     if (abs(zew(iw) + sfcg%zew(iwsfc)) < 12.e6) then  ! Skip if both near Z pole
        if (zeamin(iw) > zesfcmax + 10.) cycle
        if (zesfcmin > zeamax(iw) + 10.) cycle
     endif

     ! Evaluate x,y coordinates of current W point on gnomic plane
     ! tangent at iwsfc

     call e_gn(xew(iw), yew(iw), zew(iw), &
               sfcg%glatw(iwsfc), sfcg%glonw(iwsfc), xp0, yp0)

     xw = real(xp0,r8)
     yw = real(yp0,r8)

     ! Loop over all neighbor M points of this IW

     npoly = itab_w(iw)%npoly

     do jm = 1,npoly
        im = itab_w(iw)%im(jm)

        ! Evaluate x,y coordinates of current M point on gnomic
        ! plane tangent at iwsfc

        call e_gn(xem(im), yem(im), zem(im), &
                  sfcg%glatw(iwsfc), sfcg%glonw(iwsfc), xp0, yp0)

        xp(jm) = real(xp0,r8)
        yp(jm) = real(yp0,r8)
     enddo

     ! Evaluate possible overlap of ATM and SURFACE polygons

     call polygon_overlap2(npoly,nsfcpoly,xp,yp,xq,yq,area)

      ! Skip iw if overlap is zero

     if (area < 1.0e-7) cycle

     ! If overlap is positive, even though it may be very small, check whether
     ! this SFCG cell overlaps with the IW or any of the IM points of this ATM
     ! cell.  Where overlap is found, set topo height of ATM cell point from
     ! TOPW of SURFACE cell (unless ATM topo has already been set).

     if (topw(iw) < -1.e3) then
        call inout_check(nsfcpoly,xq,yq,xw,yw,alpha)
        if (alpha > 1.0_r8) topw(iw) = sfcg%topw(iwsfc)
     endif

     do jm = 1,npoly
        im = itab_w(iw)%im(jm)

        if (topm(im) < -1.e3) then
           call inout_check(nsfcpoly,xq,yq,xp(jm),yp(jm),alpha)
           if (alpha > 1.0_r8) topm(im) = sfcg%topw(iwsfc)
        endif
     enddo

     ! If overlap area is less than 1.e-5 of sfcg cell area, this overlap will
     ! not be counted.  

     if (area < 1.0e-5 * sfcg%area(iwsfc)) cycle

     ! This iwsfc SURFACE cell overlaps with IW ATM cell.

     nwatm = nwatm + 1

     itab_wsfc(iwsfc)%nwatm = nwatm
     itab_wsfc(iwsfc)%iwatm(nwatm) = iw
     itab_wsfc(iwsfc)%arc  (nwatm) = area

  enddo  ! iw

  ! Order multiple overlap areas so that largest is first

  if (nwatm > 1) then
     do j = 2,nwatm
        if (itab_wsfc(iwsfc)%arc(j)  > itab_wsfc(iwsfc)%arc(1)) then
           idum                      = itab_wsfc(iwsfc)%iwatm(j)
           itab_wsfc(iwsfc)%iwatm(j) = itab_wsfc(iwsfc)%iwatm(1)
           itab_wsfc(iwsfc)%iwatm(1) = idum

           dum                     = itab_wsfc(iwsfc)%arc(j)
           itab_wsfc(iwsfc)%arc(j) = itab_wsfc(iwsfc)%arc(1)
           itab_wsfc(iwsfc)%arc(1) = dum
        endif
     enddo
  endif

end subroutine sfc_atm_hex_overlay

!============================================================================

subroutine polygon_overlap2(np,nq,xp,yp,xq,yq,area)

  ! Given x,y coordinates of the vertices of polygons p and q, compute the area
  ! of overlap between the polygons using a sweepline algorithm.

  ! Method adapted from:
  ! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: np,nq ! Number of vertices in p and q
  real(r8), intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
  real(r8), intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices

  real, intent(out) :: area                  ! area of overlap of p and q

  real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

  integer :: nev      ! # of events
  integer :: nsect    ! # of intersections between strip centerline and p,q edges
  real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
  real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges 
  real(r8) :: xbot(np+nq) ! x-coords where strip bottom line intersects p and q edges 
  real(r8) :: xtop(np+nq) ! x-coords where strip top line intersects p and q edges 
  real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values

  integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,ia,ib,is
  real(r8) :: p0,q0,dx,dy,dxtrap

  real(r8) :: alpha    ! angle swept by ray from point to polygon perimeter during circuit

  ! Find vertices in p that are not outside q

  nev = 0
  area = 0.

  do ip = 1,np
     call inout_check(nq,xq,yq,xp(ip),yp(ip),alpha)

     if (abs(alpha) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yp(ip)
     endif
  enddo

  ! Find vertices in q that are not outside p

  do iq = 1,nq
     call inout_check(np,xp,yp,xq(iq),yq(iq),alpha)

     if (abs(alpha) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yq(iq)
     endif
  enddo

  ! Find intersecting edges of polygons p and q

  do ipa = 1,np
     ipb = ipa + 1
     if (ipa == np) ipb = 1

     do iqa = 1,nq
        iqb = iqa + 1
        if (iqa == nq) iqb = 1

        call intersect(xp(ipa),yp(ipa),xp(ipb),yp(ipb)  &
                      ,xq(iqa),yq(iqa),xq(iqb),yq(iqb),p0,q0,iflag)

        if (iflag == -1) cycle
        if (p0 < -0.000000001_r8 .or. p0 > 1.000000001_r8) cycle
        if (q0 < -0.000000001_r8 .or. q0 > 1.000000001_r8) cycle

  ! Line segments pa-pb and qa-qb intersect; find y-coordinate of intersection

        nev = nev + 1
        yev(nev) = (1.0_r8 - p0) * yp(ipa) + p0 * yp(ipb)
     enddo
  enddo

  ! Sort event list to obtain increasing order of yev

  call fltsort(nev,yev)

  ! Loop over event points

  do iev = 1,nev-1

     ! dy = width of strip between current and next event

     dy = yev(iev+1) - yev(iev)

     ! Reject overlap event if dy is less than 0.1 meter (threshold ok for single precision)
   
     if (dy < 0.1_r8) cycle

     ! ymid = y-coordinate of strip centerline.
     ! Initialize dx (centerline length sum) to 0.

     ymid = yev(iev) + 0.5_r8 * dy
     dx = 0.0_r8

     ! Find x-coordinate of intersections of strip centerline with edges of p and q

     nsect = 0

     do ia = 1,np
        ib = ia + 1
        if (ia == np) ib = 1

        if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

        nsect = nsect + 1
        xmid(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (ymid      -yp(ia)) / (yp(ib) - yp(ia))
        xbot(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev)  -yp(ia)) / (yp(ib) - yp(ia))
        xtop(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev+1)-yp(ia)) / (yp(ib) - yp(ia))

     enddo

     do ia = 1,nq
        ib = ia + 1
        if (ia == nq) ib = 1

        if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

        nsect = nsect + 1
        xmid(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (ymid      -yq(ia)) / (yq(ib) - yq(ia))
        xbot(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev)  -yq(ia)) / (yq(ib) - yq(ia))
        xtop(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev+1)-yq(ia)) / (yq(ib) - yq(ia))
     enddo

     ! Sort xmid values into increasing order

     call fltsort3(nsect,xmid,xbot,xtop)

     ! see if the segment is inside both polygons

     do is = 1,nsect - 1
        xcent = 0.5_r8 * (xmid(is) + xmid(is+1))

        if (xcent == xmid(is)) cycle          ! if segment length is 0

        call inout_check(np,xp,yp,xcent,ymid,alpha)
        if (abs(alpha) < 1.0_r8) cycle

        call inout_check(nq,xq,yq,xcent,ymid,alpha)
        if (abs(alpha) < 1.0_r8) cycle

        dxtrap = xmid(is+1) - xmid(is)

        dx = dx + dxtrap
     enddo

     area = area + dx * dy

  enddo

end subroutine polygon_overlap2

!============================================================================

subroutine sfc_grid_topocut(iwsfc,nwsfc,nvsfc,mwsfc,mmsfc,miw,miv,mim)

  use mem_grid,    only: nza, zm, topm, topw, &
                         xem, yem, zem, xev, yev, zev, xew, yew, zew, &
                         dnv, dnu, dniv, dniu
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_ijtabs,  only: itab_w
  use consts_coms, only: erad, piu180
  use misc_coms,   only: mdomain

  implicit none

  integer, intent(in)    :: iwsfc, nwsfc, nvsfc
  integer, intent(inout) :: mwsfc, mmsfc
  integer, intent(inout) :: miw(nwsfc), miv(nvsfc,nza), mim(nwsfc,7,nza)

  integer :: npoly, imsfc, kw, jv, jm1, jm2, iv, im1, im2, tm1, km1
  integer :: t1, t2, t3, ipat, iwsfc_fill
  integer :: jp1, jp2

  real :: z1,z2,z3,x1,x2,x3,y1,y2,y3
  real :: zm1, xm1, ym1, xm2, ym2, xem1, yem1, zem1
  real :: topc, xec, yec, zec, area
  real :: zmpat(5,nza), xmpat(5,nza), ympat(5,nza)
  integer :: tmpat(5,nza)  ! Flag denoting the edge or vertex where a new M 
                           ! point is defined
  integer :: kmpat(5,nza)
  integer :: kwpat(nza),lpoly(nza),npat
  real :: efw, efm1, efm2
  real :: pnx, pny, pnz, dvm1, dvm2, topv
  real :: xc, yc
  real :: expansion, raxis

  npoly = itab_wsfc(iwsfc)%npoly

  ! Loop over all polygon edges

  do jv = 1,npoly
     jm1 = jv - 1
     if (jm1 == 0) jm1 = npoly
     jm2 = jv

     iv  = itab_wsfc(iwsfc)%ivn(jv)
     im1 = itab_wsfc(iwsfc)%imn(jm1)
     im2 = itab_wsfc(iwsfc)%imn(jm2)

     ! Each triangular sector of a hexagonal (or pentagonal or heptagonal) grid
     ! cell, formed by the IW point and two adjacent IM points, has topographic
     ! height specified on each of its vertices, so the triangle is treated as
     ! a planar topographic surface.  Find lines on this surface where it is
     ! intersected by each zm level of the atmospheric grid, and construct
     ! polygonal areas bounded by the lines of intersection and triangular
     ! sector edges (using the algorithm in cont3f).  Each of these polygons
     ! defines a separate cell of the land/sea grid, and its area is later
     ! summed (in subroutine ctrlvols) with other cells to compute arw(k,iw)
     ! and volt(k,iw) for the atmosphere cell.

     ! Expand earth relative surface coordinates out to terrain height
     ! for computing unit normals for sloping surface cells

     efw  = (erad + sfcg%topw(iwsfc)) / erad
     efm1 = (erad + sfcg%topm(im1))  / erad
     efm2 = (erad + sfcg%topm(im2))  / erad

     call unit_normal ( sfcg%xew(iwsfc)*efw, sfcg%yew(iwsfc)*efw, sfcg%zew(iwsfc)*efw, &
                        sfcg%xem(im1)*efm1, sfcg%yem(im1)*efm1, sfcg%zem(im1)*efm1, &
                        sfcg%xem(im2)*efm2, sfcg%yem(im2)*efm2, sfcg%zem(im2)*efm2, &
                        pnx, pny, pnz )

     ! For cont3f algorithm, use local coordinate system with IW point at
     ! origin and IV edge being a line of constant positive y.

     dvm1 = sqrt((sfcg%xev(iv) - sfcg%xem(im1))**2 &
          +      (sfcg%yev(iv) - sfcg%yem(im1))**2 &
          +      (sfcg%zev(iv) - sfcg%zem(im1))**2)

     dvm2 = sqrt((sfcg%xev(iv) - sfcg%xem(im2))**2 &
          +      (sfcg%yev(iv) - sfcg%yem(im2))**2 &
          +      (sfcg%zev(iv) - sfcg%zem(im2))**2)

     topv = (dvm1 * sfcg%topm(im2) + dvm2 * sfcg%topm(im1)) * sfcg%dniu(iv)

     x1 = 0.
     x2 = dvm1
     x3 = -dvm2
     y1 = 0.
     y2 = 0.5 * sfcg%dnv(iv)
     y3 = 0.5 * sfcg%dnv(iv)
     z1 = sfcg%topw(iwsfc)
     z2 = sfcg%topm(im1)
     z3 = sfcg%topm(im2)
     t1 = 1
     t2 = 2
     t3 = 3

     if     (z1 >= z2 .and. z2 >= z3) then
        call cont3sfc2(z1,z2,z3,x1,x2,x3,y1,y2,y3,t1,t2,t3, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
     elseif (z1 >= z3 .and. z3 >= z2) then
        call cont3sfc2(z1,z3,z2,x1,x3,x2,y1,y3,y2,t1,t3,t2, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
        call reverse_polygon(zmpat,xmpat,ympat,tmpat,kmpat,lpoly,npat)
     elseif (z2 >= z1 .and. z1 >= z3) then
        call cont3sfc2(z2,z1,z3,x2,x1,x3,y2,y1,y3,t2,t1,t3, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
        call reverse_polygon(zmpat,xmpat,ympat,tmpat,kmpat,lpoly,npat)
     elseif (z2 >= z3 .and. z3 >= z1) then
        call cont3sfc2(z2,z3,z1,x2,x3,x1,y2,y3,y1,t2,t3,t1, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
     elseif (z3 >= z1 .and. z1 >= z2) then
        call cont3sfc2(z3,z1,z2,x3,x1,x2,y3,y1,y2,t3,t1,t2, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
     elseif (z3 >= z2 .and. z2 >= z1) then
        call cont3sfc2(z3,z2,z1,x3,x2,x1,y3,y2,y1,t3,t2,t1, &
                       zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)
        call reverse_polygon(zmpat,xmpat,ympat,tmpat,kmpat,lpoly,npat)
     endif

     ! Loop over new patches that have been found in this sector

     do ipat = 1,npat

        ! Increment number of surface cells and set iwsfc_fill to that number, 
        ! except when filling last patch of last polygon sector, when instead
        ! iwsfc_fill is set to the parent Voronoi cell index iwsfc

        if (jv < npoly .or. ipat < npat) then
           mwsfc = mwsfc + 1
           iwsfc_fill = mwsfc
        else
           iwsfc_fill = iwsfc
        endif

        ! Compute polygon area and centroid for current ipat

        area = 0.
        xc = 0.
        yc = 0.
        kw = kwpat(ipat)

        do jp1 = 1,lpoly(ipat)
           jp2 = jp1 + 1
           if (jp2 > lpoly(ipat)) jp2 = 1

           zm1 = zmpat(jp1,ipat)
           xm1 = xmpat(jp1,ipat)
           ym1 = ympat(jp1,ipat)
           tm1 = tmpat(jp1,ipat)
           km1 = kmpat(jp1,ipat)

           xm2 = xmpat(jp2,ipat)
           ym2 = ympat(jp2,ipat)

           area = area + 0.5 * (xm1 * ym2 - xm2 * ym1)

           xc = xc + (xm1 + xm2) * (xm1 * ym2 - xm2 * ym1)
           yc = yc + (ym1 + ym2) * (xm1 * ym2 - xm2 * ym1)

           ! Transform M points to earth coordinates

           xem1 = sfcg%xew(iwsfc) &
                + xm1 * (sfcg%xem(im1) - sfcg%xem(im2 )) * dniu(iv) &
                + ym1 * (sfcg%xev(iv ) - sfcg%xew(iwsfc)) * dniv(iv) * 2.0
           yem1 = sfcg%yew(iwsfc) &
                + xm1 * (sfcg%yem(im1) - sfcg%yem(im2 )) * dniu(iv) &
                + ym1 * (sfcg%yev(iv ) - sfcg%yew(iwsfc)) * dniv(iv) * 2.0
           zem1 = sfcg%zew(iwsfc) &
                + xm1 * (sfcg%zem(im1) - sfcg%zem(im2 )) * dniu(iv) &
                + ym1 * (sfcg%zev(iv ) - sfcg%zew(iwsfc)) * dniv(iv) * 2.0

           if (mdomain < 2) then
              expansion = erad / sqrt(xem1**2 + yem1**2 + zem1**2)

              xem1 = xem1 * expansion
              yem1 = yem1 * expansion
              zem1 = zem1 * expansion
           endif

           ! Using information extracted from cont3sfc2, determine if current M
           ! point of this patch already exists as an M point for other patches
           ! or needs to be added to surface grid

           if (tm1 == 1) then
              if (miw(iwsfc) == 0) then
                 mmsfc = mmsfc + 1
                 miw(iwsfc) = mmsfc
              endif
              imsfc = miw(iwsfc)
           elseif (tm1 == 2) then
              imsfc = im1
           elseif (tm1 == 3) then
              imsfc = im2
           elseif (tm1 == 10) then
              if (miv(iv,km1) == 0) then
                 mmsfc = mmsfc + 1
                 miv(iv,km1) = mmsfc
              endif
              imsfc = miv(iv,km1)
           elseif (tm1 == 20) then
              if (mim(iwsfc,jm2,km1) == 0) then
                 mmsfc = mmsfc + 1
                 mim(iwsfc,jm2,km1) = mmsfc
              endif
              imsfc = mim(iwsfc,jm2,km1)
           elseif (tm1 == 30) then
              if (mim(iwsfc,jm1,km1) == 0) then
                 mmsfc = mmsfc + 1
                 mim(iwsfc,jm1,km1) = mmsfc
              endif
              imsfc = mim(iwsfc,jm1,km1)
           endif

           ! Enter current M point into itab_wsfc table

           itab_wsfc(iwsfc_fill)%imn(jp1) = imsfc

           ! If M point is being added, populate it with spatial information

           if (imsfc == mmsfc) then
              sfcg%xem(imsfc)  = xem1
              sfcg%yem(imsfc)  = yem1
              sfcg%zem(imsfc)  = zem1
              sfcg%topm(imsfc) = zm1

              raxis = sqrt(xem1**2 + yem1**2)

              sfcg%glatm(imsfc) = atan2(zem1,raxis) * piu180
              sfcg%glonm(imsfc) = atan2(yem1,xem1) * piu180
           endif

        enddo ! jp1

        xc = xc / (6. * area)
        yc = yc / (6. * area)

        ! Interpolate earth coordinates and topo height to centroid

        xec = sfcg%xew(iwsfc) &
            + xc * (sfcg%xem(im1) - sfcg%xem(im2 )) * dniu(iv) &
            + yc * (sfcg%xev(iv ) - sfcg%xew(iwsfc)) * dniv(iv) * 2.0
        yec = sfcg%yew(iwsfc) &
            + xc * (sfcg%yem(im1) - sfcg%yem(im2 )) * dniu(iv) &
            + yc * (sfcg%yev(iv ) - sfcg%yew(iwsfc)) * dniv(iv) * 2.0
        zec = sfcg%zew(iwsfc) &
            + xc * (sfcg%zem(im1) - sfcg%zem(im2 )) * dniu(iv) &
            + yc * (sfcg%zev(iv ) - sfcg%zew(iwsfc)) * dniv(iv) * 2.0

        topc = sfcg%topw(iwsfc) &
             + xc * (sfcg%topm(im1) - sfcg%topm(im2 )) * dniu(iv) &
             + yc * (    topv      - sfcg%topw(iwsfc)) * dniv(iv) * 2.0

        if (mdomain < 2) then

           expansion = erad / sqrt(xec**2 + yec**2 + zec**2)

           xec = xec * expansion
           yec = yec * expansion
           zec = zec * expansion

           ! Compute latitude and longitude of centroid; catalogue sfc cell variables

           raxis = sqrt(xec**2 + yec**2)

           sfcg%glatw(iwsfc_fill) = atan2(zec,raxis) * piu180
           sfcg%glonw(iwsfc_fill) = atan2(yec,xec) * piu180

        else

           sfcg%glatw(iwsfc_fill) = 0.  ! want it this way?
           sfcg%glonw(iwsfc_fill) = 0.  ! want it this way?

        endif

        sfcg%area(iwsfc_fill) = area

        sfcg%xew(iwsfc_fill) = xec
        sfcg%yew(iwsfc_fill) = yec
        sfcg%zew(iwsfc_fill) = zec

        sfcg%topw(iwsfc_fill) = topc

        sfcg%wnx(iwsfc_fill) = pnx
        sfcg%wny(iwsfc_fill) = pny
        sfcg%wnz(iwsfc_fill) = pnz

        itab_wsfc(iwsfc_fill)%ivoronoi = 0
        itab_wsfc(iwsfc_fill)%nwatm = 1
        itab_wsfc(iwsfc_fill)%iwatm(1) = iwsfc
        itab_wsfc(iwsfc_fill)%kwatm(1) = kwpat(ipat)
        itab_wsfc(iwsfc_fill)%npoly    = lpoly(ipat)
        itab_wsfc(iwsfc_fill)%arc(1)   = area

     enddo  ! ipat

  enddo   ! jv

end subroutine sfc_grid_topocut

!===============================================================================

subroutine cont3sfc2(z1,z2,z3,x1,x2,x3,y1,y2,y3,t1,t2,t3,zmpat,xmpat,ympat,tmpat,kmpat,kwpat,lpoly,npat)

  use mem_grid, only: nza, zm

  implicit none

  real,    intent(in   ) :: z1,z2,z3,x1,x2,x3,y1,y2,y3
  integer, intent(in   ) :: t1,t2,t3
  real,    intent(inout) :: zmpat(5,nza), xmpat(5,nza), ympat(5,nza)
  integer, intent(inout) :: tmpat(5,nza), kmpat(5,nza)
  integer, intent(out  ) :: kwpat(nza),lpoly(nza),npat

  integer :: to1, to2, km, iflag
  real :: zo1, zo2, xo1, xo2, yo1, yo2, contlev ! "old" values of xmpat1,xmpat2,ympat1,ympat2

  ! Determine lowest zm level that exceeds z3, the lowest triangle vertex

  km = 1
  do while (z3 >= zm(km) .and. km < nza-1)
     km = km + 1
  enddo
  contlev = zm(km)
  npat = 1
  kmpat(:,:) = 0

  ! If all 3 points are in the same contour interval

  if (contlev >= z1 .or. contlev <= z3) then

     zmpat(1,npat) = z1
     xmpat(1,npat) = x1
     ympat(1,npat) = y1
     tmpat(1,npat) = t1

     zmpat(2,npat) = z2
     xmpat(2,npat) = x2
     ympat(2,npat) = y2
     tmpat(2,npat) = t2

     zmpat(3,npat) = z3
     xmpat(3,npat) = x3
     ympat(3,npat) = y3
     tmpat(3,npat) = t3

     kwpat(npat) = km

     lpoly(npat) = 3

     return
  endif

  ! Draw all contour lines of value less than z1

  iflag = 0

  do while (contlev < z1)

     if (iflag > 0) then
        zo1 = zmpat(1,npat-1)
        xo1 = xmpat(1,npat-1)
        yo1 = ympat(1,npat-1)
        to1 = tmpat(1,npat-1)

        zo2 = zmpat(2,npat-1)
        xo2 = xmpat(2,npat-1)
        yo2 = ympat(2,npat-1)
        to2 = tmpat(2,npat-1)
     endif

     if ( abs(z1-z3) > 1.e-25 ) then
        zmpat(1,npat) = contlev
        xmpat(1,npat) = x3 + (x1-x3) * (contlev-z3) / (z1-z3)
        ympat(1,npat) = y3 + (y1-y3) * (contlev-z3) / (z1-z3)
        tmpat(1,npat) = t2*10
        kmpat(1,npat) = km
     else
        zmpat(1,npat) = z3
        xmpat(1,npat) = x3
        ympat(1,npat) = y3
        tmpat(1,npat) = t3
     endif

     if (contlev < z2) then

        if ( abs(z2-z3) > 1.e-25 ) then
           zmpat(2,npat) = contlev
           xmpat(2,npat) = x3 + (x2-x3) * (contlev-z3) / (z2-z3)
           ympat(2,npat) = y3 + (y2-y3) * (contlev-z3) / (z2-z3)
           tmpat(2,npat) = t1*10
           kmpat(2,npat) = km
        else
           zmpat(2,npat) = z3
           xmpat(2,npat) = x3
           ympat(2,npat) = y3
           tmpat(2,npat) = t3
        endif

        if (iflag == 0) then  ! lowest contour interval: z3 is a node
           zmpat(3,npat) = z3
           xmpat(3,npat) = x3   
           ympat(3,npat) = y3
           tmpat(3,npat) = t3

           kwpat(npat) = km
           lpoly(npat) = 3
        else
           zmpat(3,npat) = zo2
           xmpat(3,npat) = xo2
           ympat(3,npat) = yo2
           tmpat(3,npat) = to2
           kmpat(3,npat) = km - 1

           zmpat(4,npat) = zo1
           xmpat(4,npat) = xo1
           ympat(4,npat) = yo1
           tmpat(4,npat) = to1
           kmpat(4,npat) = km - 1

           kwpat(npat) = km
           lpoly(npat) = 4
        endif

        if (km == nza-1) then  ! highest model level: z1 and z2 are nodes
                               ! (should not actually happen)
           zmpat(3,npat) = z2
           xmpat(3,npat) = x2
           ympat(3,npat) = y2
           tmpat(3,npat) = t2

           zmpat(4,npat) = z1         
           xmpat(4,npat) = x1         
           ympat(4,npat) = y1
           tmpat(4,npat) = t1

           kwpat(npat) = km
           lpoly(npat) = 4

           return
        elseif (zm(km+1) >= z1) then  ! highest contour interval:
                                      ! z1 and z2 are nodes
           zmpat(1,npat+1) = zmpat(2,npat)
           xmpat(1,npat+1) = xmpat(2,npat)
           ympat(1,npat+1) = ympat(2,npat)
           tmpat(1,npat+1) = tmpat(2,npat)
           kmpat(1,npat+1) = km

           zmpat(2,npat+1) = zmpat(1,npat)
           xmpat(2,npat+1) = xmpat(1,npat)
           ympat(2,npat+1) = ympat(1,npat)
           tmpat(2,npat+1) = tmpat(1,npat)
           kmpat(2,npat+1) = km

           zmpat(3,npat+1) = z1
           xmpat(3,npat+1) = x1
           ympat(3,npat+1) = y1
           tmpat(3,npat+1) = t1

           zmpat(4,npat+1) = z2         
           xmpat(4,npat+1) = x2         
           ympat(4,npat+1) = y2
           tmpat(4,npat+1) = t2

           kwpat(npat+1) = km + 1
           lpoly(npat+1) = 4

           npat = npat + 1

           return
        endif

        iflag = 1

     else   

        if ( abs(z1-z2) > 1.e-25 ) then
           zmpat(2,npat) = contlev
           xmpat(2,npat) = x2 + (x1-x2) * (contlev-z2) / (z1-z2)
           ympat(2,npat) = y2 + (y1-y2) * (contlev-z2) / (z1-z2)
           tmpat(2,npat) = t3*10
           kmpat(2,npat) = km
        else
           zmpat(2,npat) = z2
           xmpat(2,npat) = x2
           ympat(2,npat) = y2
           tmpat(2,npat) = t2
        endif

        if (iflag == 0) then  ! lowest contour color: z2 and z3 are nodes
           zmpat(3,npat) = z2
           xmpat(3,npat) = x2
           ympat(3,npat) = y2
           tmpat(3,npat) = t2

           zmpat(4,npat) = z3        
           xmpat(4,npat) = x3        
           ympat(4,npat) = y3         
           tmpat(4,npat) = t3

           kwpat(npat) = km
           lpoly(npat) = 4
        elseif (iflag == 1) then ! switching from contlev < z2 to contlev > z2
           zmpat(3,npat) = z2
           xmpat(3,npat) = x2
           ympat(3,npat) = y2
           tmpat(3,npat) = t2

           zmpat(4,npat) = zo2
           xmpat(4,npat) = xo2
           ympat(4,npat) = yo2
           tmpat(4,npat) = to2
           kmpat(4,npat) = km - 1

           zmpat(5,npat) = zo1         
           xmpat(5,npat) = xo1         
           ympat(5,npat) = yo1
           tmpat(5,npat) = to1
           kmpat(5,npat) = km - 1

           kwpat(npat) = km
           lpoly(npat) = 5
        else                     ! keeping to contlev > z2
           zmpat(3,npat) = zo2
           xmpat(3,npat) = xo2
           ympat(3,npat) = yo2
           tmpat(3,npat) = to2
           kmpat(3,npat) = km - 1

           zmpat(4,npat) = zo1
           xmpat(4,npat) = xo1
           ympat(4,npat) = yo1
           tmpat(4,npat) = to1
           kmpat(4,npat) = km - 1

           kwpat(npat) = km
           lpoly(npat) = 4
        endif

        if (km == nza-1) then  ! highest model level: z1 is a node 
           zmpat(3,npat) = z1  ! (should not actually happen)
           xmpat(3,npat) = x1  ! (should not actually happen)
           ympat(3,npat) = y1
           tmpat(3,npat) = t1

           kwpat(npat) = km
           lpoly(npat) = 3

           return
        elseif (zm(km+1) >= z1) then  ! highest contour interval: z1 is a node
           zmpat(1,npat+1) = zmpat(2,npat)
           xmpat(1,npat+1) = xmpat(2,npat)
           ympat(1,npat+1) = ympat(2,npat)
           tmpat(1,npat+1) = tmpat(2,npat)
           kmpat(1,npat+1) = km

           zmpat(2,npat+1) = zmpat(1,npat)
           xmpat(2,npat+1) = xmpat(1,npat)
           ympat(2,npat+1) = ympat(1,npat)
           tmpat(2,npat+1) = tmpat(1,npat)
           kmpat(2,npat+1) = km

           zmpat(3,npat+1) = z1
           xmpat(3,npat+1) = x1
           ympat(3,npat+1) = y1
           tmpat(3,npat+1) = t1

           kwpat(npat+1) = km + 1
           lpoly(npat+1) = 3

           npat = npat + 1

           return
        endif

        iflag = 2

     endif

     npat = npat + 1
     km = km + 1
     contlev = zm(km)

  enddo

end subroutine cont3sfc2

!==========================================================================

subroutine reverse_polygon(zmpat, xmpat, ympat, tmpat, kmpat, lpoly, npat)

  use mem_grid, only: nza

  implicit none

  real,    intent(inout) :: zmpat(5,nza), xmpat(5,nza), ympat(5,nza)
  integer, intent(inout) :: tmpat(5,nza), kmpat(5,nza)
  integer, intent(in   ) :: npat, lpoly(nza)

  integer :: ipat, j, jc, istore
  real :: store

  do ipat = 1,npat
     do j = 1,lpoly(ipat) / 2
        jc = lpoly(ipat) + 1 - j

        store = zmpat(j,ipat)
        zmpat(j,ipat) = zmpat(jc,ipat)
        zmpat(jc,ipat) = store

        store = xmpat(j,ipat)
        xmpat(j,ipat) = xmpat(jc,ipat)
        xmpat(jc,ipat) = store

        store = ympat(j,ipat)
        ympat(j,ipat) = ympat(jc,ipat)
        ympat(jc,ipat) = store

        istore = tmpat(j,ipat)
        tmpat(j,ipat) = tmpat(jc,ipat)
        tmpat(jc,ipat) = istore

        istore = kmpat(j,ipat)
        kmpat(j,ipat) = kmpat(jc,ipat)
        kmpat(jc,ipat) = istore
     enddo
  enddo

end subroutine reverse_polygon

