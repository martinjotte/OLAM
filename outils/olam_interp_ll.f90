subroutine interp_htw_ll(npts,iws_loc,wts_loc,nlevin,nlevout,field,field_ll)

  use mem_grid,     only: mwa, mza, lpw
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: jtab_w, jtw_prog

  implicit none

  integer, intent(in)    :: npts, nlevin, nlevout
  integer, intent(in)    :: iws_loc(npts,3)
  real,    intent(in)    :: wts_loc(npts,3)
  real,    intent(inout) :: field(nlevin,mwa)
  real,    intent(inout) :: field_ll(npts,nlevout)

  integer :: ipt, kin, kout, j, iw, ka, k

  if (nlevin == mza) then
     !$omp parallel do private (iw,ka,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        ka = lpw(iw)
        do k = ka-1, 1, -1
           field(k,iw) = field(ka,iw)
        enddo
     enddo
     !$omp end parallel do
  endif

  if (iparallel == 1) then
     call mpi_send_w(svara1=field)
     call mpi_recv_w(svara1=field)
  endif

  !$omp parallel do private (kout,kin)
  do ipt = 1, npts
     do kout = 1, nlevout
        kin = kout + nlevin - nlevout
        field_ll(ipt,kout) = wts_loc(ipt,1) * field(kin,iws_loc(ipt,1)) &
                           + wts_loc(ipt,2) * field(kin,iws_loc(ipt,2)) &
                           + wts_loc(ipt,3) * field(kin,iws_loc(ipt,3))
     enddo
  enddo
  !$omp end parallel do

end subroutine interp_htw_ll

!================================================================================

subroutine find_3iws_ll(nlon,nlat,alon,alat,iws_ll,wts_ll)

  use mem_grid,   only: xew, yew, zew, dnv
  use mem_ijtabs, only: itab_w, jtab_w, jtw_prog
  use consts_coms,only: pio180, r8, erad, eradi

  implicit none

  integer, intent(in)    :: nlon, nlat
  real,    intent(in)    :: alon(nlon), alat(nlat)
  integer, intent(inout) :: iws_ll(nlon,nlat,3)
  real,    intent(inout) :: wts_ll(nlon,nlat,3)

  integer :: ilat, ilon, j, j1, j2, jw, iw, iwn, npoly

  real :: xwn(7), ywn(7), dot00(7), dot01(7), dot11(7), denomi(7), dn(7)
  real :: qx, qy, dnvmax, dist, distn
  real :: dot02,dot12,u,v

  real :: xea(nlon,nlat), yea(nlon,nlat), zea(nlat)
  real :: dxe, dye, dze, rads

  real :: coswlon, sinwlon
  real :: coswlat, sinwlat

  real :: cosalon(nlon), sinalon(nlon)
  real :: cosalat(nlat), sinalat(nlat)

  real :: raxis, raxisi

  real, parameter :: fuzz = 0.001

  ! Compute and store earth coordinates (xea,yea,zea) of each lat/lon point

  do ilat = 1, nlat
     rads          = alat(ilat) * pio180
     cosalat(ilat) = cos(rads)
     sinalat(ilat) = sin(rads)
  enddo

  do ilon = 1, nlon
     rads          = alon(ilon) * pio180
     cosalon(ilon) = cos(rads)
     sinalon(ilon) = sin(rads)
  enddo

  do ilat = 1, nlat
     zea(ilat) = erad * sinalat(ilat)
     raxis     = erad * cosalat(ilat)
     do ilon = 1, nlon
        xea(ilon,ilat) = raxis * cosalon(ilon)
        yea(ilon,ilat) = raxis * sinalon(ilon)
     enddo
  enddo

  ! Loop over all prognostic W points

  do jw = 1, jtab_w(jtw_prog)%jend
     iw = jtab_w(jtw_prog)%iw(jw)

     raxis  = sqrt( xew(iw)**2 + yew(iw)**2 )

     sinwlat = zew(iw) * eradi
     coswlat = raxis   * eradi

     ! For points less than 100 m from Earth's polar axis, make arbitrary
     ! assumption that longitude = 0 deg.  This is just to settle on a PS
     ! planar coordinate system in which to do the algebra.

     if (raxis >= 1.e2) then
        raxisi = 1.0 / raxis
        sinwlon = yew(iw) * raxisi
        coswlon = xew(iw) * raxisi
     else
        sinwlon = 0.
        coswlon = 1.
     endif

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
        dn    (j1) = 1. / (xwn(j2) * ywn(j1) - xwn(j1) * ywn(j2))
     enddo

     ! Loop over all lat/lon points and determine which are closer to current
     ! W point than to any other W point on globe.  It is sufficient to show
     ! that a lat/lon point is closer to current W point than to any neighbor W
     ! point, AND that it is closer to current W point than most distant
     ! neighbor W point is.

     do ilat = 1, nlat
        if (abs(zea(ilat) - zew(iw)) > dnvmax) cycle

        lonloop: do ilon = 1, nlon

           ! Skip this lat/lon point if its nearest IW point has already been found

           if (iws_ll(ilon,ilat,1) > 1) cycle

           ! Skip this lat/lon point if any of its earth coordinates differs from
           ! lat/lon point by more than dnvmax

           if (abs(xea(ilon,ilat) - xew(iw)) > dnvmax) cycle
           if (abs(yea(ilon,ilat) - yew(iw)) > dnvmax) cycle

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
