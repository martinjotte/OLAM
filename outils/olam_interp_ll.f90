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

subroutine interp_htw_ll(npts,iws_loc,wts_loc,nlevin,nlevout,field,field_ll)

  use mem_grid,     only: mwa
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w

  implicit none

  integer, intent(in)    :: npts, nlevin, nlevout
  integer, intent(in)    :: iws_loc(npts,3)
  real,    intent(in)    :: wts_loc(npts,3)
  real,    intent(inout) :: field(nlevin,mwa)
  real,    intent(inout) :: field_ll(npts,nlevout)

  integer :: ipt, kin, kout

  if (iparallel == 1) then
     call mpi_send_w(1,svara=field)
     call mpi_recv_w(1,svara=field)
  endif

  do ipt = 1, npts
     do kout = 1, nlevout
        kin = kout + nlevin - nlevout
        field_ll(ipt,kout) = wts_loc(ipt,1) * field(kin,iws_loc(ipt,1)) &
                           + wts_loc(ipt,2) * field(kin,iws_loc(ipt,2)) &
                           + wts_loc(ipt,3) * field(kin,iws_loc(ipt,3))
     enddo
  enddo

end subroutine interp_htw_ll

!================================================================================

subroutine find_3iws_ll(nlon,nlat,alon,alat,iws_ll,wts_ll)

  use mem_grid,   only: glatw, glonw, mwa, nwa, xew, yew, zew, dnv
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
  use consts_coms,only: pio180, r8, erad

  implicit none

  integer, intent(in)    :: nlon, nlat
  real,    intent(in)    :: alon(nlon), alat(nlat)
  integer, intent(inout) :: iws_ll(nlon,nlat,3)
  real,    intent(inout) :: wts_ll(nlon,nlat,3)

  integer :: ilat, ilon, j, j1, j2, jw, iw, iwn, np, jnext, npoly

  real :: xwn(7), ywn(7), dot00(7), dot01(7), dot11(7), denomi(7), dn(7)
  real :: qx, qy, raxis, dnvmax, dist, distn
  real :: dot02,dot12,u,v

  real :: xea(nlon,nlat), yea(nlon,nlat), zea(nlat)
  real :: dxe, dye, dze

  real :: coswlon, sinwlon
  real :: coswlat, sinwlat

  real, parameter :: fuzz = 0.001

  ! Compute and store earth coordinates (xea,yea,zea) of each lat/lon point

  do ilat = 1, nlat
     zea(ilat) = erad * sin(alat(ilat) * pio180)
     raxis     = erad * cos(alat(ilat) * pio180)

     do ilon = 1, nlon
        xea(ilon,ilat) = raxis * sin(alon(ilon) * pio180)
        yea(ilon,ilat) = raxis * cos(alon(ilon) * pio180)
     enddo
  enddo

  ! Loop over all prognostic W points

  do jw = 1, jtab_w(jtw_prog)%jend(1)
     iw = jtab_w(jtw_prog)%iw(jw)

     sinwlat = sin(glatw(iw) * pio180)
     coswlat = cos(glatw(iw) * pio180)
     sinwlon = sin(glonw(iw) * pio180)
     coswlon = cos(glonw(iw) * pio180)

     ! Find max distance to neighbor W points

     npoly = itab_w(iw)%npoly
     dnvmax = maxval(dnv(itab_w(iw)%iv(1:npoly)))

     ! Loop over neighbor W points

     do j = 1,npoly
        iwn = itab_w(iw)%iw(j)

        ! Transform neighbor W points to PS coordinates tangent at IW point

        dxe = xew(iwn) - xew(iw)
        dye = yew(iwn) - yew(iw)
        dze = zew(iwn) - zew(iw)

        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xwn(j),ywn(j))
     enddo
           
     ! Loop over each pair of consecutive neighbor W points, and set up
     ! triangle-check coefficients that depend only on W points

     do j1 = 1,npoly
        j2 = j1 + 1
        if (j1 == npoly) j2 = 1

        dot00(j1) = xwn(j1) * xwn(j1) + ywn(j1) * ywn(j1)
        dot01(j1) = xwn(j1) * xwn(j2) + ywn(j1) * ywn(j2)
        dot11(j1) = xwn(j2) * xwn(j2) + ywn(j2) * ywn(j2)

        denomi(j1) = 1. / (dot00(j1) * dot11(j1) - dot01(j1) * dot01(j1))
        dn(j1)     = 1. / (xwn(j2) * ywn(j1) - xwn(j1) * ywn(j2))
     enddo
           
     ! Loop over all lat/lon points and determine which are closer to current
     ! W point than to any other W point on globe.  It is sufficient to show
     ! that a lat/lon point is closer to current W point than to any neighbor W
     ! point, AND that it is closer to current W point than most distant
     ! neighbor W point is.

     do ilat = 1, nlat
        lonloop: do ilon = 1, nlon

           ! Skip this lat/lon point if its nearest IW point has already been found

           if (iws_ll(ilon,ilat,1) > 1) cycle

           ! Skip this lat/lon point if any of its earth coordinates differs from
           ! lat/lon point by more than dnvmax

           if (abs(xea(ilon,ilat) - xew(iw)) > dnvmax) cycle
           if (abs(yea(ilon,ilat) - yew(iw)) > dnvmax) cycle
           if (abs(zea(ilat)      - zew(iw)) > dnvmax) cycle

           ! Compute distance between IW point and lat/lon point

           dist = sqrt((xea(ilon,ilat)-xew(iw))**2 &
                     + (yea(ilon,ilat)-yew(iw))**2 &
                     + (zea(ilat)     -zew(iw))**2)

           ! Skip this lat/lon point if it is farther from IW than most distant
           ! neighbor IWN is

           if (dist > dnvmax) cycle

           ! Loop over neighbor W points

           do j = 1,npoly
              iwn = itab_w(iw)%iw(j)

              ! Compute distance between lat/lon point and IWN point

              distn = sqrt((xea(ilon,ilat)-xew(iwn))**2 &
                         + (yea(ilon,ilat)-yew(iwn))**2 &
                         + (zea(ilat)     -zew(iwn))**2)

              ! If lat/lon point is closer to IWN point than to IW point, move
              ! on to next lat/lon point.  Bias is used to reduce chance of
              ! lat/lon point being rejected by all IW points in domain; this 
              ! might lead to a few lat/lon values being interpolated on
              ! multiple MPI subdomains, but this is sorted out later.

              if (distn < 0.999999 * dist) cycle lonloop
           enddo

           ! If this point was reached, current lat/lon point is inside IW cell.
           ! Store IW index for the lat/lon point and transform lat/lon point
           ! to PS coordinates tangent at IW point.

           dxe = xea(ilon,ilat) - xew(iw)
           dye = yea(ilon,ilat) - yew(iw)
           dze = zea(ilat)      - zew(iw)

           call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,qx,qy)

           ! Loop over each pair of consecutive neighbor W points 

           do j1 = 1,npoly
              j2 = j1 + 1
              if (j1 == npoly) j2 = 1

              ! Set up triangle-check coefficients that depend on lat/lon points

              dot02 = xwn(j1) * qx + ywn(j1) * qy
              dot12 = xwn(j2) * qx + ywn(j2) * qy

              u = (dot11(j1) * dot02 - dot01(j1) * dot12) * denomi(j1)
              v = (dot00(j1) * dot12 - dot01(j1) * dot02) * denomi(j1)

              if (u > -fuzz .and. v > -fuzz .and. u + v < 1.0 + fuzz) then

                 ! lat/lon point is inside current triangle; store indices of
                 ! both IWN neighbors

                 iws_ll(ilon,ilat,1) = iw
                 iws_ll(ilon,ilat,2) = itab_w(iw)%iw(j1)
                 iws_ll(ilon,ilat,3) = itab_w(iw)%iw(j2)

                 ! Compute and store 3 interpolation weights

                 wts_ll(ilon,ilat,2) = dn(j1) * (-ywn(j2) * qx + xwn(j2) * qy)
                 wts_ll(ilon,ilat,3) = dn(j1) * ( ywn(j1) * qx - xwn(j1) * qy)
                 wts_ll(ilon,ilat,1) = 1. - wts_ll(ilon,ilat,2) - wts_ll(ilon,ilat,3)

                 exit
              endif
           enddo  ! j1 loop

        enddo lonloop

     enddo ! ilat loop

  enddo  ! iw loop

end subroutine find_3iws_ll
