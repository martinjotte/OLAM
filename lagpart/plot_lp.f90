subroutine plot_lp(iplt)

  use oplot_coms, only: op
  use plotcolors, only: clrtab
  use mem_grid,   only: glatw, glonw
  use mem_lp,     only: nlp, atp

  implicit none

  integer, intent(in) :: iplt

  integer :: ilp
  integer :: iw, kk, icolor, notavail, ival, itab

  real :: hp, vp, dp, wttop, fldval, bsize

  real :: delh, delv
  character(len=6) :: title

!print*, 'in plot_lp ', iplt, nlp

  kk = 2     ! dummy value
  wttop = 0. ! dummy value

  call oplot_set(iplt)

  call gslwsc(3.)    ! line thickness

  ! Loop over particles

  do ilp = 1,nlp

     ! Skip this particle if it is not active

     if (atp(ilp)%ncluster == 0) cycle

!print*, 'plp0 ', iplt, ilp, nlp, op%projectn(iplt)

     ! Get plot coordinates of current particle

     if ( op%projectn(iplt) == 'L' .or.  &
          op%projectn(iplt) == 'P' .or.  &
          op%projectn(iplt) == 'G' .or.  &
          op%projectn(iplt) == 'O' .or.  &
          op%projectn(iplt) == 'Z') then  ! For projection onto horizontal surface

        iw = atp(ilp)%iw
        call oplot_transform(iplt,atp(ilp)%xep0,atp(ilp)%yep0,atp(ilp)%zep0,glonw(iw),glatw(iw),hp,vp)
        dp = 1.0  ! dummy value

!write(6,'(a,2i7,9e12.2)') 'plp1 ', iplt, ilp, atp(ilp)%xep0, atp(ilp)%yep0, atp(ilp)%zep0, hp, vp, dp

     else  ! For projection onto vertical surface

        call coneplot_tri_lp(iplt, ilp, atp(ilp)%xep0, atp(ilp)%yep0, atp(ilp)%zep0, hp, dp)
        vp = atp(ilp)%zp

!write(6,'(a,2i7,9e12.2)') 'plp3 ', iplt, ilp, atp(ilp)%xep0, atp(ilp)%yep0, atp(ilp)%zep0, hp, vp, dp

     endif

     if ( hp < op%xmin .or. hp > op%xmax .or.  &
          vp < op%ymin .or. vp > op%ymax ) cycle

     ! Get particle value (dp uses wtbot position in oplot_lib call)

     ! Get particle value

     call oplot_lib(kk, ilp, 'VALUE', op%fldname(iplt), dp, wttop, fldval, notavail)

     if (notavail > 0) cycle

     if (trim(op%fldname(iplt)) == 'LP') then

        ! Plot LP locations with black color if not plotting any LP values

        icolor = 10

     else

        ! Extract contour color from color table if plotting LP value

        ival = 1
        itab = op%icolortab(iplt)
        do while (fldval > clrtab(itab)%vals(ival) .and. &
                  ival < clrtab(itab)%nvals              )
           ival = ival + 1
        enddo
        icolor = clrtab(itab)%ipal(ival)

     endif

     call o_sflush
     call o_gsplci(icolor)
     call o_gstxci(icolor)

     if (.true.) then

        ! draw individual LPs

        delh = 0.0002 * (op%xmax - op%xmin)
        delv = 0.0002 * (op%ymax - op%ymin)
        call o_frstpt (hp-delh,vp-delv)
        call o_vector (hp+delh,vp-delv)
        call o_vector (hp+delh,vp+delv)
        call o_vector (hp-delh,vp+delv)
        call o_vector (hp-delh,vp-delv)

        ! Plot particle index (use 'K' for this in OLAMIN/PLTSPEC2 ?)

        write (title, '(i6)') ilp
        bsize = .0005 * (op%hp2 - op%hp1)
        call o_plchhq(hp,vp+delv*4.,trim(adjustl(title)),bsize,0.,0.)
     else

        ! Fill time series arrays of successive LP locations and draw vectors to show trajectory.

     endif

  enddo

end subroutine plot_lp

!===============================================================================

subroutine coneplot_tri_lp(iplt, ilp, xep0, yep0, zep0, hp, dp)

  ! Version of coneplot that transforms earth coordinates of a Lagrangian particle
  ! to plot coordinates (lateral and perpendicular to plot window)

  use oplot_coms,  only: op
  use consts_coms, only: pio180, erad, eradi, piu180, pi2
  use misc_coms,   only: mdomain

  implicit none

  integer, intent(in) :: iplt, ilp
  real, intent(in) :: xep0,yep0,zep0
  real, intent(out) :: hp, dp

  real :: unxec3,unyec3,unzec3
  real :: xcent,ycent,zcent
  real :: xwincent,ywincent,zwincent
  real :: radcone
  real :: raxis
  real :: vecx,vecy,vecz
  real :: vecleftx,veclefty,vecleftz
  real :: valx,valy
  real :: angep
  real :: sinplat,cosplat
  real :: sinvaz,cosvaz
  real :: sinconang,cosconang
  real :: sinconlat,cosconlat
  real :: sinconlon,cosconlon
  real :: sp

  ! If Cartesian domain (plon3 and plat3 are x,y coords in meters in this case)

  if (mdomain >= 2) then

     sinvaz = sin((90. - op%viewazim) * pio180)
     cosvaz = cos((90. - op%viewazim) * pio180)

     sp = (xep0 - op%plon3) * cosvaz + (yep0 - op%plat3) * sinvaz

     dp = sp - op%slabloc(iplt)

  else

     ! cone is viewed from INSIDE, so cone axis is in opposite direction from viewazim

     sinvaz = -sin(op%viewazim * pio180)
     cosvaz = -cos(op%viewazim * pio180)

     sinplat = sin(op%plat3 * pio180)
     cosplat = cos(op%plat3 * pio180)

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

     ! Compute angle of particle (ray from earth center) with plot cone center using dot products
     ! (Adequate precision may require cross product for cone angles
     ! close to 0 or 180)

     sp = piu180 * acos((xep0 * unxec3 + yep0 * unyec3 + zep0 * unzec3) * eradi)

     ! Compute distance (m) of particle from cone surface (negative means inside cone
     ! and therefore closer to viewer at cone axis)

     dp = erad * pio180 * (sp - op%coneang)

  endif

  ! Transform horizontal point coordinates to distances along cone

  if (mdomain >= 2) then

     ! For Cartesian grid

     hp = (xep0 - op%plon3) * sinvaz - (yep0 - op%plat3) * cosvaz

  else

     ! On sphere...

     ! Components of vector from circle center to plot window center

     vecx = xwincent - xcent
     vecy = ywincent - ycent
     vecz = zwincent - zcent

     ! Components of vector 90 degrees to the left (in azimuth) from preceding vector

     vecleftx = unyec3 * vecz - unzec3 * vecy
     veclefty = unzec3 * vecx - unxec3 * vecz
     vecleftz = unxec3 * vecy - unyec3 * vecx

     ! Compute dot product between vector from circle center to plot center
     ! and vector from circle center to LP

     valx = vecx * (xep0 - xcent)  &
          + vecy * (yep0 - ycent)  &
          + vecz * (zep0 - zcent)

     ! Repeat with 90-left vector

     valy = vecleftx * (xep0 - xcent)  &
          + veclefty * (yep0 - ycent)  &
          + vecleftz * (zep0 - zcent)

     angep = atan2(-valy,valx)  ! Angle increases clockwise

     ! Scale angle to htpn coordinates (in meters along cone circle)

     hp = angep * radcone

  endif

end subroutine coneplot_tri_lp

