subroutine coneplot_tri(iplt,iw,xe1,ye1,ze1,xe2,ye2,ze2,xe3,ye3,ze3, &
                        wta1,wta2,wta3,wtb1,wtb2,wtb3,iok,htpn)

! Version of coneplot that checks whether cone surface intersects a specified
! surface triangle, and if it does, also obtains horizontal coordinates of both
! intersection points

use oplot_coms,  only: op
use consts_coms, only: pio180, erad, eradi, piu180, pi2
use misc_coms,   only: mdomain

implicit none

integer, intent(in) :: iplt, iw
real, intent(in) :: xe1,ye1,ze1,xe2,ye2,ze2,xe3,ye3,ze3
real, intent(out) :: wta1,wta2,wta3,wtb1,wtb2,wtb3
integer, intent(out) :: iok
real, intent(out) :: htpn(4)

real :: unxec3,unyec3,unzec3
real :: xcent,ycent,zcent
real :: xwincent,ywincent,zwincent
real :: radcone
real :: raxis
real :: vecx,vecy,vecz
real :: vecleftx,veclefty,vecleftz
real :: valx,valy
real :: s1, s2, s3
real :: angea,angeb
real :: sinplat,cosplat
real :: sinvaz,cosvaz
real :: sinconang,cosconang
real :: sinconlat,cosconlat
real :: sinconlon,cosconlon
real :: xea, yea, zea, xeb, yeb, zeb
real :: slabloc

! If Cartesian domain (plon3 and plat3 are x,y coords in meters in this case)

if (mdomain >= 2) then

   sinvaz = sin((90. - op%viewazim) * pio180)
   cosvaz = cos((90. - op%viewazim) * pio180)

   s1 = (xe1 - op%plon3) * cosvaz + (ye1 - op%plat3) * sinvaz
   s2 = (xe2 - op%plon3) * cosvaz + (ye2 - op%plat3) * sinvaz
   s3 = (xe3 - op%plon3) * cosvaz + (ye3 - op%plat3) * sinvaz

   slabloc = op%slabloc(iplt)

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

! Compute angle of triangle vertices with plot cone center using dot products
! (Adequate precision may require cross product for cone angles
! close to 0 or 180)

   s1 = piu180 * acos((xe1 * unxec3 + ye1 * unyec3 + ze1 * unzec3) * eradi)
   s2 = piu180 * acos((xe2 * unxec3 + ye2 * unyec3 + ze2 * unzec3) * eradi)
   s3 = piu180 * acos((xe3 * unxec3 + ye3 * unyec3 + ze3 * unzec3) * eradi)

   slabloc = op%coneang

endif

! Initialize all weights to zero for points A and B

wta1 = 0.
wta2 = 0.
wta3 = 0.
wtb1 = 0.
wtb2 = 0.
wtb3 = 0.

! Return with iok = 0 if triangle column is not intersected by plot slab
! or touches slab at only 1 vertex

iok = 0

if (min(s1,s2,s3) > slabloc .or. max(s1,s2,s3) < slabloc) return

if (s1 >= slabloc .and. s2 > slabloc .and. s3 > slabloc) return
if (s2 >= slabloc .and. s3 > slabloc .and. s1 > slabloc) return
if (s3 >= slabloc .and. s1 > slabloc .and. s2 > slabloc) return

if (s1 <= slabloc .and. s2 < slabloc .and. s3 < slabloc) return
if (s2 <= slabloc .and. s3 < slabloc .and. s1 < slabloc) return
if (s3 <= slabloc .and. s1 < slabloc .and. s2 < slabloc) return

iok = 1  ! Since we got here, IW column is intersected by plot cone or slab

! Check first interval

if (s1 <= slabloc .and. s2 >= slabloc) then

! Interval 1:2 touches cone

   if (s1 == s2) then

! Interval 1:2 is tangent to cone

      if (s3 > slabloc) then

! Interval 1:2 is on minimum side of triangle

         wta1 = 1.
         wtb2 = 1.
         go to 10

      else

! Interval 1:2 is on maximum side of triangle

         wta2 = 1.
         wtb1 = 1.
         go to 10

      endif

   else

! Interval 1:2 touches cone at 1 point

      wtb2 = (slabloc - s1) / (s2 - s1)
      wtb1 = 1. - wtb2

   endif

elseif (s1 >= slabloc .and. s2 <= slabloc) then

! Interval 1:2 touches cone at 1 point

   wta1 = (slabloc - s2) / (s1 - s2)
   wta2 = 1. - wta1

endif

! Check second interval

if (s2 <= slabloc .and. s3 >= slabloc) then

! Interval 2:3 touches cone

   if (s2 == s3) then

! Interval 2:3 is tangent to cone

      if (s1 > slabloc) then

! Interval 2:3 is on minimum side of triangle

         wta2 = 1.
         wtb3 = 1.
         go to 10

      else

! Interval 2:3 is on maximum side of triangle

         wta3 = 1.
         wtb2 = 1.
         go to 10

      endif

   else

! Interval 2:3 touches cone at 1 point

      wtb3 = (slabloc - s2) / (s3 - s2)
      wtb2 = 1. - wtb3

   endif

elseif (s2 >= slabloc .and. s3 <= slabloc) then

! Interval 2:3 touches cone at 1 point

   wta2 = (slabloc - s3) / (s2 - s3)
   wta3 = 1. - wta2

endif

! Check third interval

if (s3 <= slabloc .and. s1 >= slabloc) then

! Interval 3:1 touches cone

   if (s3 == s1) then

! Interval 3:1 is tangent to cone

      if (s2 > slabloc) then

! Interval 3:1 is on minimum side of triangle

         wta3 = 1.
         wtb1 = 1.
         go to 10

      else

! Interval 3:1 is on maximum side of triangle

         wta1 = 1.
         wtb3 = 1.
         go to 10

      endif

   else

! Interval 3:1 touches cone at 1 point

      wtb1 = (slabloc - s3) / (s1 - s3)
      wtb3 = 1. - wtb1

   endif

elseif (s3 >= slabloc .and. s1 <= slabloc) then

! Interval 3:1 touches cone at 1 point

   wta3 = (slabloc - s1) / (s3 - s1)
   wta1 = 1. - wta3

endif

10 continue

xea = wta1 * xe1 + wta2 * xe2 + wta3 * xe3
yea = wta1 * ye1 + wta2 * ye2 + wta3 * ye3
zea = wta1 * ze1 + wta2 * ze2 + wta3 * ze3

xeb = wtb1 * xe1 + wtb2 * xe2 + wtb3 * xe3
yeb = wtb1 * ye1 + wtb2 * ye2 + wtb3 * ye3
zeb = wtb1 * ze1 + wtb2 * ze2 + wtb3 * ze3

! Transform horizontal point coordinates to distances along cone

if (mdomain >= 2) then

! For Cartesian grid

   htpn(1) = (xea - op%plon3) * sinvaz - (yea - op%plat3) * cosvaz
   htpn(2) = (xeb - op%plon3) * sinvaz - (yeb - op%plat3) * cosvaz

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
! and vector from circle center to current point

   valx = vecx * (xea - xcent)  &
        + vecy * (yea - ycent)  &
        + vecz * (zea - zcent)

! Repeat with 90-left vector

   valy = vecleftx * (xea - xcent)  &
        + veclefty * (yea - ycent)  &
        + vecleftz * (zea - zcent)

   angea = atan2(-valy,valx)  ! Angle increases clockwise

! Repeat dot product for second point

   valx = vecx * (xeb - xcent)  &
        + vecy * (yeb - ycent)  &
        + vecz * (zeb - zcent)

   valy = vecleftx * (xeb - xcent)  &
        + veclefty * (yeb - ycent)  &
        + vecleftz * (zeb - zcent)

   angeb = atan2(-valy,valx)  ! Angle increases clockwise

! Avoid wrap_around

   if (angeb < angea) angeb = angeb + pi2

! Scale angles to htpn coordinates (in meters along cone circle)

   htpn(1) = angea * radcone
   htpn(2) = angeb * radcone

endif

htpn(3) = htpn(2)
htpn(4) = htpn(1)

return
end subroutine coneplot_tri

!===============================================================================

subroutine coneplot_w(iw,iv1,iv2,topo1,topo2,iok,htpn)

use oplot_coms,  only: op, xepc, yepc, zepc
use mem_grid,    only: xem, yem, zem, topm
use mem_ijtabs,  only: itab_w
use consts_coms, only: pio180, erad, piu180, pi2, pi1

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
real :: wt1,wt2
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

   jv = jm2
   iv  = itab_w(iw)%iv(jv)

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

if (ange2 + pi1 < ange1) ange2 = ange2 + pi2

! Scale angles to htpn coordinates (in meters along cone circle)

htpn(1) = ange1 * radcone
htpn(2) = ange2 * radcone
htpn(3) = htpn(2)
htpn(4) = htpn(1)

end subroutine coneplot_w
