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
subroutine ll_xy (qlat,qlon,polelat,polelon,x,y)

implicit none

real, intent(in) :: qlat
real, intent(in) :: qlon
real, intent(in) :: polelat
real, intent(in) :: polelon

real, intent(out) :: x
real, intent(out) :: y

real, parameter :: erad = 6.367e6
real, parameter :: erad2 = 1.2734e7
real, parameter :: pi180 = 3.14159265 / 180.

real :: sinplat
real :: cosplat
real :: sinplon
real :: cosplon
real :: sinqlat
real :: cosqlat
real :: sinqlon
real :: cosqlon
real :: x3p
real :: y3p
real :: z3p
real :: z3q
real :: x3q
real :: y3q
real :: xq
real :: yq
real :: zq
real :: t
       
! This subroutine computes cartesian coordinates (x,y) in a polar stereographic
! projection whose pole point is located at geographic latitude-longitude
! coordinates (polelat,polelon) of a point located at geographic 
! latitude-longitude coordinates (qlat,qlon).

! Evaluate sine and cosine of latitude and longitude of pole point p and
! input point q.

sinplat = sin(polelat * pi180)
cosplat = cos(polelat * pi180)
sinplon = sin(polelon * pi180)
cosplon = cos(polelon * pi180)

sinqlat = sin(qlat * pi180)
cosqlat = cos(qlat * pi180)
sinqlon = sin(qlon * pi180)
cosqlon = cos(qlon * pi180)

! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime
! meridian, and the y axis is the equator and 90 E.

! For the pole point, these are:

x3p = erad * cosplat * cosplon
y3p = erad * cosplat * sinplon
z3p = erad * sinplat

! For the given lat,lon point, these are:

z3q = erad * sinqlat
x3q = erad * cosqlat * cosqlon
y3q = erad * cosqlat * sinqlon

! Transform q point from (x3,y3,z3) coordinates in the above system to
! polar stereographic coordinates (x,y,z):

xq = - sinplon * (x3q-x3p) + cosplon * (y3q-y3p)
yq =   cosplat * (z3q-z3p)  &
     - sinplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )
zq =   sinplat * (z3q-z3p)  &
     + cosplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )

! Parametric equation for line from antipodal point at (0,0,-2 erad) to
! point q has the following parameter (t) value on the polar stereographic
! plane:

t = erad2 / (erad2 + zq)

! This gives the following x and y coordinates for the projection of point q
! onto the polar stereographic plane:

x = xq * t
y = yq * t

return
end subroutine ll_xy

!============================================================================

subroutine xy_ll (qlat,qlon,polelat,polelon,x,y)
implicit none

real, intent(out) :: qlat
real, intent(out) :: qlon

real, intent(in) :: polelat
real, intent(in) :: polelon
real, intent(in) :: x
real, intent(in) :: y

real, parameter :: erad = 6.367e6
real, parameter :: erad2 = 1.2734e7
real, parameter :: pi180 = 3.14159265 / 180.

real :: sinplat
real :: cosplat
real :: sinplon
real :: cosplon
real :: x3p
real :: y3p
real :: z3p
real :: z3q
real :: x3q
real :: y3q
real :: xq
real :: yq
real :: zq
real :: t
real :: d
real :: alpha
real :: r3q

! This subroutine computes latitude-longitude coordinates (qlat,qlon) at
! a point located at cartesian coordinates (x,y) in a polar stereographic
! projection whose pole point is located at geographic latitude-longitude
! coordinates (polelat,polelon).

! Evaluate sine and cosine of latitude and longitude of pole point p.

sinplat = sin(polelat * pi180)
cosplat = cos(polelat * pi180)
sinplon = sin(polelon * pi180)
cosplon = cos(polelon * pi180)

! Compute (x3,y3,z3) coordinates of the pole point where the origin is the
! center of the earth, the z axis is the north pole, the x axis is the
! equator and prime meridian, and the y axis is the equator and 90 E.

x3p = erad * cosplat * cosplon
y3p = erad * cosplat * sinplon
z3p = erad * sinplat

! Compute distance d from given point R on the polar stereographic plane
! to the pole point P:

d = sqrt (x ** 2 + y ** 2)

! Compute angle QCP where C is the center of the Earth.  This is twice
! angle QAP where A is the antipodal point.  Angle QAP is the same as
! angle RAP:

alpha = 2. * atan2(d,erad2)

! Compute zq, the height of Q relative to the polar stereographic plane:

zq = erad * (cos(alpha) - 1.)

! Compute the parameter t which is the the distance ratio AQ:AR

t = (erad2 + zq) / erad2

! Compute xq and yq, the x and y coordinates of Q in polar stereographic space:

xq = t * x
yq = t * y

! Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):

x3q = x3p - xq * sinplon - yq * cosplon * sinplat  &
    + zq * cosplat * cosplon
y3q = y3p + xq * cosplon - yq * sinplon * sinplat  &
    + zq * cosplat * sinplon
z3q = z3p + yq * cosplat + zq * sinplat

! Compute the latitude and longitude of Q:

qlon = atan2(y3q,x3q) / pi180
r3q = sqrt(x3q ** 2 + y3q ** 2)
qlat = atan2(z3q,r3q) / pi180

return
end subroutine xy_ll

!============================================================================

subroutine e_ps(xeq,yeq,zeq,polelat,polelon,x,y)

use consts_coms, only: pio180, erad, erad2

implicit none

real, intent(in) :: xeq
real, intent(in) :: yeq
real, intent(in) :: zeq
real, intent(in) :: polelat
real, intent(in) :: polelon

real, intent(out) :: x
real, intent(out) :: y

real :: sinplat
real :: cosplat
real :: sinplon
real :: cosplon
real :: xep
real :: yep
real :: zep
real :: xq
real :: yq
real :: zq
real :: t

! This subroutine computes coordinates (x,y) of point q projected onto polar
! stereographic plane whose pole point is located at geographic coordinates
! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
! "earth cartesian space", where the origin is the center of the earth, 
! the z axis is the north pole, the x axis is the equator and prime meridian, 
! and the y axis is the equator and 90 E..
       
! Evaluate sine and cosine of latitude and longitude of pole point p

sinplat = sin(polelat * pio180)
cosplat = cos(polelat * pio180)
sinplon = sin(polelon * pio180)
cosplon = cos(polelon * pio180)

! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

xep = erad * cosplat * cosplon
yep = erad * cosplat * sinplon
zep = erad * sinplat

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! polar stereographic plane with origin at the pole point, the z axis pointing 
! radially outward from the center of the earth, and the y axis pointing 
! northward along the local earth meridian from the pole point.

xq =                                  - sinplon * (xeq-xep) + cosplon * (yeq-yep)
yq =   cosplat * (zeq-zep) - sinplat * (cosplon * (xeq-xep) + sinplon * (yeq-yep))
zq =   sinplat * (zeq-zep) + cosplat * (cosplon * (xeq-xep) + sinplon * (yeq-yep))

! Parametric equation for line from antipodal point at (0,0,-2 erad) in 3D 
! coordinates of polar stereographic plane to point q has the following 
! parameter (t) value on the polar stereographic plane (zq <= 0):

t = erad2 / (erad2 + zq)

! This gives the following x and y coordinates for the projection of point q
! onto the polar stereographic plane:

x = xq * t
y = yq * t

return
end subroutine e_ps

!============================================================================

subroutine e_or(xeq,yeq,zeq,polelat,polelon,x,y)

use consts_coms, only: pio180, erad

implicit none

real, intent(in) :: xeq
real, intent(in) :: yeq
real, intent(in) :: zeq
real, intent(in) :: polelat
real, intent(in) :: polelon

real, intent(out) :: x
real, intent(out) :: y

real :: sinplat
real :: cosplat
real :: sinplon
real :: cosplon
real :: xep
real :: yep
real :: zep
real :: xq
real :: yq
real :: zq
real :: t

! This subroutine computes coordinates (x,y) of point q projected onto 
! orthographic plane whose pole point is located at geographic coordinates
! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
! "earth cartesian space", where the origin is the center of the earth, 
! the z axis is the north pole, the x axis is the equator and prime meridian, 
! and the y axis is the equator and 90 E..
       
! Evaluate sine and cosine of latitude and longitude of pole point p

sinplat = sin(polelat * pio180)
cosplat = cos(polelat * pio180)
sinplon = sin(polelon * pio180)
cosplon = cos(polelon * pio180)

! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

xep = erad * cosplat * cosplon
yep = erad * cosplat * sinplon
zep = erad * sinplat

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! polar stereographic plane with origin at the pole point, the z axis pointing 
! radially outward from the center of the earth, and the y axis pointing 
! northward along the local earth meridian from the pole point.

xq = - sinplon * (xeq-xep) + cosplon * (yeq-yep)
yq =   cosplat * (zeq-zep) - sinplat * (cosplon * (xeq-xep) + sinplon * (yeq-yep))
zq =   sinplat * (zeq-zep) + cosplat * (cosplon * (xeq-xep) + sinplon * (yeq-yep))

if (zq < -erad) then   ! point is on other side of earth - do not plot
   x = 1.e12  ! This will serve as a flag to not plot point
   y = 1.e12  ! This will serve as a flag to not plot point
else
   x = xq
   y = yq
endif

return
end subroutine e_or

!============================================================================

subroutine e_ec(xeq,yeq,zeq,rlon,rlat)

use consts_coms, only: piu180

implicit none

real, intent(in) :: xeq
real, intent(in) :: yeq
real, intent(in) :: zeq

real, intent(out) :: rlon
real, intent(out) :: rlat

real :: rax

! This subroutine computes coordinates (rlon,rlat) of point q
! projected onto an equidistant cylindrical surface.  Input coordinates 
! (xeq,yeq,zeq) of q are defined in "earth cartesian space", where the origin 
! is the center of the earth, the z axis is the north pole, the x axis is the 
! equator and prime meridian, and the y axis is the equator and 90 E.
       
rax = sqrt(xeq ** 2 + yeq ** 2)  ! distance of q from earth axis
rlat = atan2(zeq,rax) * piu180
rlon = atan2(yeq,xeq) * piu180

return
end subroutine e_ec

!============================================================================

subroutine ps_e(xeq,yeq,zeq,polelat,polelon,x,y)

use consts_coms, only: pio180, erad, erad2

implicit none

real, intent(out) :: xeq
real, intent(out) :: yeq
real, intent(out) :: zeq

real, intent(in) :: polelat
real, intent(in) :: polelon
real, intent(in) :: x
real, intent(in) :: y

real :: sinplat
real :: cosplat
real :: sinplon
real :: cosplon
real :: xep
real :: yep
real :: zep
real :: xq
real :: yq
real :: zq
real :: t
real :: alpha

! Given coordinates (x,y) of point q projected onto polar stereographic plane
! whose pole point is located at geographic coordinates (polelat,polelon), 
! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
! surface defined in "earth cartesian space", where the origin is the center
! of the earth, the z axis is the north pole, the x axis is the equator and 
! prime meridian, and the y axis is the equator and 90 E..

! Evaluate sine and cosine of latitude and longitude of pole point p

sinplat = sin(polelat * pio180)
cosplat = cos(polelat * pio180)
sinplon = sin(polelon * pio180)
cosplon = cos(polelon * pio180)

! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

xep = erad * cosplat * cosplon
yep = erad * cosplat * sinplon
zep = erad * sinplat

! Parametric equation for line from antipodal point at (0,0,-2 erad) in 3D 
! coordinates of polar stereographic plane to point (x,y) on the polar stereographic
! plane has the following parameter (t) value on the earth's surface (zq <= 0):

alpha = 2. * atan2(sqrt(x**2 + y**2),erad2)

t = .5 * (1. + cos(alpha))

! This gives the following xq, yq, zq coordinates relative to the polar
! stereographic plane for the projection of point q from the polar stereographic
! plane to the earth surface:

xq = x * t
yq = y * t
zq = erad2 * (t - 1.)

! Transform q point located on the earth's surface from ps coordinates (xq,yq,zq)
! to earth coordinates (xe,ye,ze).  The polar stereographic plane has its origin
! at the pole point, with the z axis pointing radially outward from the center
! of the earth, and the y axis pointing northward along the local earth meridian
! from the pole point.

xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
yeq = yep + cosplon * xq - sinplon * ( sinplat * yq + cosplat * zq)
zeq = zep                            + cosplat * yq + sinplat * zq

return
end subroutine ps_e


