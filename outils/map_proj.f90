subroutine ll_xy (qlat,qlon,polelat,polelon,x,y)

use consts_coms, only: erad, erad2, pio180

implicit none

real, intent(in) :: qlat
real, intent(in) :: qlon
real, intent(in) :: polelat
real, intent(in) :: polelon

real, intent(out) :: x
real, intent(out) :: y

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

sinplat = sin(polelat * pio180)
cosplat = cos(polelat * pio180)
sinplon = sin(polelon * pio180)
cosplon = cos(polelon * pio180)

sinqlat = sin(qlat * pio180)
cosqlat = cos(qlat * pio180)
sinqlon = sin(qlon * pio180)
cosqlon = cos(qlon * pio180)

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

end subroutine ll_xy

!============================================================================

subroutine xy_ll (qlat,qlon,polelat,polelon,x,y)

use consts_coms, only: erad, erad2, pio180

implicit none

real, intent(out) :: qlat
real, intent(out) :: qlon

real, intent(in) :: polelat
real, intent(in) :: polelon
real, intent(in) :: x
real, intent(in) :: y

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

sinplat = sin(polelat * pio180)
cosplat = cos(polelat * pio180)
sinplon = sin(polelon * pio180)
cosplon = cos(polelon * pio180)

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

qlon = atan2(y3q,x3q) / pio180
r3q = sqrt(x3q ** 2 + y3q ** 2)
qlat = atan2(z3q,r3q) / pio180

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

end subroutine e_ps

!============================================================================

subroutine e_gn(xeq,yeq,zeq,polelat,polelon,x,y)

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
! gnomonic tangent plane whose pole point is located at geographic coordinates
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
! gnomonic tangent plane with origin at the pole point, the z axis pointing
! radially outward from the center of the earth, and the y axis pointing
! northward along the local earth meridian from the pole point.

xq =                                  - sinplon * (xeq-xep) + cosplon * (yeq-yep)
yq =   cosplat * (zeq-zep) - sinplat * (cosplon * (xeq-xep) + sinplon * (yeq-yep))
zq =   sinplat * (zeq-zep) + cosplat * (cosplon * (xeq-xep) + sinplon * (yeq-yep))

if (zq < -erad) then   ! point is on other side of earth - do not plot
   x = 1.e12  ! This will serve as a flag to not plot point
   y = 1.e12  ! This will serve as a flag to not plot point
else

   ! Parametric equation for line from earth center point at (0,0,-erad) in 3D
   ! coordinates of gnomonic tangent plane to point q has the following
   ! parameter (t) value on the gnomonic tangent plane (zq <= 0):

   t = erad / (erad + zq)

   ! This gives the following x and y coordinates for the projection of point q
   ! onto the gnomonic tangent plane:

   x = xq * t
   y = yq * t

endif

end subroutine e_gn

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

end subroutine e_or

!============================================================================

subroutine de_or(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

  use consts_coms, only: pio180, erad

  implicit none

  real, intent(in) :: dxe
  real, intent(in) :: dye
  real, intent(in) :: dze

  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon

  real, intent(out) :: x
  real, intent(out) :: y

  real :: xq
  real :: yq
  real :: zq

! This subroutine computes coordinates (x,y) of point q projected onto
! orthographic plane whose pole point is located at geographic coordinates
! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
! "earth cartesian space", where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime meridian,
! and the y axis is the equator and 90 E..

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! polar stereographic plane with origin at the pole point, the z axis pointing
! radially outward from the center of the earth, and the y axis pointing
! northward along the local earth meridian from the pole point.

  xq = - sinplon * dxe + cosplon * dye
  yq =   cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)
  zq =   sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)

  if (zq < -erad) then   ! point is on other side of earth - do not plot
     x = 1.e12           ! This will serve as a flag to not plot point
     y = 1.e12           ! This will serve as a flag to not plot point
  else
     x = xq
     y = yq
  endif

end subroutine de_or

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

end subroutine e_ec

!============================================================================

subroutine ec_e(rlon,rlat,xeq,yeq,zeq)

  use consts_coms, only: pio180, erad

  implicit none

  real, intent(out) :: xeq
  real, intent(out) :: yeq
  real, intent(out) :: zeq

  real, intent( in) :: rlon
  real, intent( in) :: rlat

  real :: sinrlat
  real :: cosrlat
  real :: sinrlon
  real :: cosrlon

  ! This subroutine computes "earth cartesian space" coordinates
  ! (xeq,yeq,zeq) or a point located at (rlon,rlat).

  sinrlat = sin(rlat * pio180)
  cosrlat = cos(rlat * pio180)
  sinrlon = sin(rlon * pio180)
  cosrlon = cos(rlon * pio180)

  xeq = cosrlon * cosrlat * erad
  yeq = sinrlon * cosrlat * erad
  zeq =           sinrlat * erad

end subroutine ec_e

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
yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
zeq = zep                            + cosplat * yq + sinplat * zq

end subroutine ps_e

!============================================================================

subroutine gn_e(xeq,yeq,zeq,polelat,polelon,x,y)

use consts_coms, only: pio180, erad, eradsq

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
!real :: alpha

! Given coordinates (x,y) of point q projected onto gnomonic tangent plane
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

! Parametric equation for line from earth center point at (0,0,-erad) in 3D
! coordinates of gnomonic tangent plane to point (x,y) on the tangent gnomonic
! plane has the following parameter (t) value on the earth's surface (zq <= 0):

! alpha = atan2(sqrt(x**2 + y**2),erad)
! t     = cos(alpha)
t     = erad / sqrt(x*x + y*y + eradsq)

! This gives the following xq, yq, zq coordinates relative to the polar
! stereographic plane for the projection of point q from the polar stereographic
! plane to the earth surface:

xq = x * t
yq = y * t
zq = erad * (t - 1.)

! Transform q point located on the earth's surface from ps coordinates (xq,yq,zq)
! to earth coordinates (xe,ye,ze).  The polar stereographic plane has its origin
! at the pole point, with the z axis pointing radially outward from the center
! of the earth, and the y axis pointing northward along the local earth meridian
! from the pole point.

xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
zeq = zep                            + cosplat * yq + sinplat * zq

end subroutine gn_e

!============================================================================

subroutine de_ps(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

use consts_coms, only: pio180, erad2

implicit none

real, intent(in) :: dxe
real, intent(in) :: dye
real, intent(in) :: dze

real, intent(in) :: cosplat
real, intent(in) :: sinplat
real, intent(in) :: cosplon
real, intent(in) :: sinplon

real, intent(out) :: x
real, intent(out) :: y

real :: xq
real :: yq
real :: zq
real :: t

! This subroutine computes coordinates (x,y) of point q projected onto polar
! stereographic plane given the sines and cosines of the pole point
! located at geographic coordinates (polelat,polelon).
! Input vector (dxe,dye,dze) is the distance of point q from the pole point
! in "earth cartesian space", where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime meridian,
! and the y axis is the equator and 90 E..

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! polar stereographic plane with origin at the pole point, the z axis pointing
! radially outward from the center of the earth, and the y axis pointing
! northward along the local earth meridian from the pole point.

xq =                           - sinplon * dxe + cosplon * dye
yq =  cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)
zq =  sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)

! Parametric equation for line from antipodal point at (0,0,-2 erad) in 3D
! coordinates of polar stereographic plane to point q has the following
! parameter (t) value on the polar stereographic plane (zq <= 0):

t = erad2 / (erad2 + zq)

! This gives the following x and y coordinates for the projection of point q
! onto the polar stereographic plane:

x = xq * t
y = yq * t

end subroutine de_ps

!============================================================================

subroutine de_ps_mult(n,dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

use consts_coms, only: pio180, erad2

implicit none

integer, intent(in) :: n

real, intent(in) :: dxe(n)
real, intent(in) :: dye(n)
real, intent(in) :: dze(n)

real, intent(in) :: cosplat
real, intent(in) :: sinplat
real, intent(in) :: cosplon
real, intent(in) :: sinplon

real, intent(out) :: x(n)
real, intent(out) :: y(n)

real :: xq
real :: yq
real :: zq
real :: t
integer :: i

! This subroutine computes coordinates (x,y) of point q projected onto polar
! stereographic plane given the sines and cosines of the pole point
! located at geographic coordinates (polelat,polelon).
! Input vector (dxe,dye,dze) is the distance of point q from the pole point
! in "earth cartesian space", where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime meridian,
! and the y axis is the equator and 90 E..

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! polar stereographic plane with origin at the pole point, the z axis pointing
! radially outward from the center of the earth, and the y axis pointing
! northward along the local earth meridian from the pole point.

do i = 1, n

   xq =                              - sinplon * dxe(i) + cosplon * dye(i)
   yq =  cosplat * dze(i) - sinplat * (cosplon * dxe(i) + sinplon * dye(i))
   zq =  sinplat * dze(i) + cosplat * (cosplon * dxe(i) + sinplon * dye(i))

! Parametric equation for line from antipodal point at (0,0,-2 erad) in 3D
! coordinates of polar stereographic plane to point q has the following
! parameter (t) value on the polar stereographic plane (zq <= 0):

   t = erad2 / (erad2 + zq)

! This gives the following x and y coordinates for the projection of point q
! onto the polar stereographic plane:

   x(i) = xq * t
   y(i) = yq * t

enddo

end subroutine de_ps_mult

!============================================================================

subroutine ll_xy_mult (n,qlat,qlon,cosplat,sinplat,cosplon,sinplon,x3p,y3p,z3p,x,y)

  use consts_coms, only: erad, erad2, pio180

  implicit none

  integer, intent(in) :: n

  real, intent(in) :: qlat(n)
  real, intent(in) :: qlon(n)
  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon
  real, intent(in) :: x3p
  real, intent(in) :: y3p
  real, intent(in) :: z3p

  real, intent(out) :: x(n)
  real, intent(out) :: y(n)

  real :: sinqlat
  real :: cosqlat
  real :: sinqlon
  real :: cosqlon
  real :: z3q
  real :: x3q
  real :: y3q
  real :: xq
  real :: yq
  real :: zq
  real :: t

  integer :: i

! This subroutine computes cartesian coordinates (x,y) in a polar stereographic
! projection whose pole point is located at geographic latitude-longitude
! coordinates (polelat,polelon) of a point located at geographic
! latitude-longitude coordinates (qlat,qlon).

! Evaluate sine and cosine of latitude and longitude of pole point p and
! input point q.

  do i = 1, n

     sinqlat = sin(qlat(i) * pio180)
     cosqlat = cos(qlat(i) * pio180)
     sinqlon = sin(qlon(i) * pio180)
     cosqlon = cos(qlon(i) * pio180)

! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime
! meridian, and the y axis is the equator and 90 E.
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

     x(i) = xq * t
     y(i) = yq * t

  enddo

end subroutine ll_xy_mult

!============================================================================

subroutine ps_de(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

  use consts_coms, only: erad2, erad2sq

  implicit none

  real, intent(out) :: dxe
  real, intent(out) :: dye
  real, intent(out) :: dze

  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon

  real, intent(in) :: x
  real, intent(in) :: y

  real :: xq
  real :: yq
  real :: zq
  real :: t
! real :: alpha

! Given coordinates (x,y) of point q projected onto polar stereographic plane
! whose pole point is located at geographic coordinates (polelat,polelon),
! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
! surface defined in "earth cartesian space", where the origin is the center
! of the earth, the z axis is the north pole, the x axis is the equator and
! prime meridian, and the y axis is the equator and 90 E..

! Parametric equation for line from antipodal point at (0,0,-2 erad) in 3D
! coordinates of polar stereographic plane to point (x,y) on the polar stereographic
! plane has the following parameter (t) value on the earth's surface (zq <= 0):

! alpha = 2. * atan2(sqrt(x**2 + y**2),erad2)
! t     = .5 * (1. + cos(alpha))
  t     = erad2sq / (x*x + y*y + erad2sq)

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

  dxe = -sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
  dye =  cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
  dze =  cosplat * yq + sinplat * zq

end subroutine ps_de

!============================================================================

subroutine de_gn(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

  use consts_coms, only: pio180, erad

  implicit none

  real, intent(in) :: dxe
  real, intent(in) :: dye
  real, intent(in) :: dze

  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon

  real, intent(out) :: x
  real, intent(out) :: y

  real :: xq
  real :: yq
  real :: zq
  real :: t

! This subroutine computes coordinates (x,y) of point q projected onto
! gnomonic tangent plane whose pole point is located at geographic coordinates
! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
! "earth cartesian space", where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime meridian,
! and the y axis is the equator and 90 E..

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! gnomonic tangent plane with origin at the pole point, the z axis pointing
! radially outward from the center of the earth, and the y axis pointing
! northward along the local earth meridian from the pole point.

  xq =                          - sinplon * dxe + cosplon * dye
  yq = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)
  zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)

  if (zq < -.999 * erad) then  ! point is on other side of earth - do not plot

     x = 1.e12                 ! This will serve as a flag to not plot point
     y = 1.e12                 ! This will serve as a flag to not plot point

  else

! Parametric equation for line from earth center point at (0,0,-erad) in 3D
! coordinates of gnomonic tangent plane to point q has the following
! parameter (t) value on the gnomonic tangent plane (zq <= 0):

     t = erad / (erad + zq)

! This gives the following x and y coordinates for the projection of point q
! onto the gnomonic tangent plane:

     x = xq * t
     y = yq * t

  endif

end subroutine de_gn

!============================================================================

subroutine de_gn_mult(n,dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

  use consts_coms, only: pio180, erad

  implicit none

  integer, intent(in) :: n

  real, intent(in) :: dxe(n)
  real, intent(in) :: dye(n)
  real, intent(in) :: dze(n)

  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon

  real, intent(out) :: x(n)
  real, intent(out) :: y(n)

  real :: xq
  real :: yq
  real :: zq
  real :: t
  integer :: i

! This subroutine computes coordinates (x,y) of point q projected onto
! gnomonic tangent plane whose pole point is located at geographic coordinates
! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
! "earth cartesian space", where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime meridian,
! and the y axis is the equator and 90 E..

! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
! gnomonic tangent plane with origin at the pole point, the z axis pointing
! radially outward from the center of the earth, and the y axis pointing
! northward along the local earth meridian from the pole point.

  do i = 1, n

     xq =                             - sinplon * dxe(i) + cosplon * dye(i)
     yq = cosplat * dze(i) - sinplat * (cosplon * dxe(i) + sinplon * dye(i))
     zq = sinplat * dze(i) + cosplat * (cosplon * dxe(i) + sinplon * dye(i))

! Parametric equation for line from earth center point at (0,0,-erad) in 3D
! coordinates of gnomonic tangent plane to point q has the following
! parameter (t) value on the gnomonic tangent plane (zq <= 0):

     t = erad / max(erad + zq, 1.) ! Guard if point is on other side of earth

! This gives the following x and y coordinates for the projection of point q
! onto the gnomonic tangent plane:

     x(i) = min(xq * t, 1.e12)   ! 1.e12 will serve as a flag to not plot point
     y(i) = min(yq * t, 1.e12)   ! 1.e12 will serve as a flag to not plot point

  enddo

end subroutine de_gn_mult

!============================================================================

subroutine gn_de(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

  use consts_coms, only: erad, eradsq

  implicit none

  real, intent(out) :: dxe
  real, intent(out) :: dye
  real, intent(out) :: dze

  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon

  real, intent(in) :: x
  real, intent(in) :: y

  real :: xq
  real :: yq
  real :: zq
  real :: t
! real :: alpha

! Given coordinates (x,y) of point q projected onto gnomonic tangent plane
! whose pole point is located at geographic coordinates (polelat,polelon),
! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
! surface defined in "earth cartesian space", where the origin is the center
! of the earth, the z axis is the north pole, the x axis is the equator and
! prime meridian, and the y axis is the equator and 90 E..

! Parametric equation for line from earth center point at (0,0,-erad) in 3D
! coordinates of gnomonic tangent plane to point (x,y) on the gnomonic tangent
! plane has the following parameter (t) value on the earth's surface (zq <= 0):

! alpha = atan2(sqrt(x**2 + y**2),erad)
! t     = cos(alpha)
  t     = erad / sqrt(x*x + y*y + eradsq)

! This gives the following xq, yq, zq coordinates relative to the polar
! stereographic plane for the projection of point q from the polar stereographic
! plane to the earth surface:

  xq = x * t
  yq = y * t
  zq = erad * (t - 1.)

! Transform q point located on the earth's surface from ps coordinates (xq,yq,zq)
! to earth coordinates (xe,ye,ze).  The polar stereographic plane has its origin
! at the pole point, with the z axis pointing radially outward from the center
! of the earth, and the y axis pointing northward along the local earth meridian
! from the pole point.

  dxe = - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
  dye =   cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
  dze =                            + cosplat * yq + sinplat * zq

end subroutine gn_de

!============================================================================

subroutine ll_gn_mult (n,qlat,qlon,cosplat,sinplat,cosplon,sinplon,x3p,y3p,z3p,x,y)

  use consts_coms, only: erad, erad2, pio180

  implicit none

  integer, intent(in) :: n

  real, intent(in) :: qlat(n)
  real, intent(in) :: qlon(n)
  real, intent(in) :: cosplat
  real, intent(in) :: sinplat
  real, intent(in) :: cosplon
  real, intent(in) :: sinplon
  real, intent(in) :: x3p
  real, intent(in) :: y3p
  real, intent(in) :: z3p

  real, intent(out) :: x(n)
  real, intent(out) :: y(n)

  real :: sinqlat
  real :: cosqlat
  real :: sinqlon
  real :: cosqlon
  real :: dxe
  real :: dye
  real :: dze
  real :: xq
  real :: yq
  real :: zq
  real :: t

  integer :: i

! This subroutine computes cartesian coordinates (x,y) in a polar stereographic
! projection whose pole point is located at geographic latitude-longitude
! coordinates (polelat,polelon) of a point located at geographic
! latitude-longitude coordinates (qlat,qlon).

! Evaluate sine and cosine of latitude and longitude of pole point p and
! input point q.

  do i = 1, n

     sinqlat = sin(qlat(i) * pio180)
     cosqlat = cos(qlat(i) * pio180)
     sinqlon = sin(qlon(i) * pio180)
     cosqlon = cos(qlon(i) * pio180)

! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
! the z axis is the north pole, the x axis is the equator and prime
! meridian, and the y axis is the equator and 90 E.
! For the given lat,lon point, these are:

     dze = erad * sinqlat           - z3p
     dxe = erad * cosqlat * cosqlon - x3p
     dye = erad * cosqlat * sinqlon - y3p

! Transform q point from (x3,y3,z3) coordinates in the above system to
! polar stereographic coordinates (x,y,z):

     xq =                          - sinplon * dxe + cosplon * dye
     yq = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)
     zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)

! Parametric equation for line from earth center point at (0,0,-erad) in 3D
! coordinates of gnomonic tangent plane to point q has the following
! parameter (t) value on the gnomonic tangent plane (zq <= 0):

     t = erad / max(erad + zq, 1.) ! Guard if point is on other side of earth

! This gives the following x and y coordinates for the projection of point q
! onto the gnomonic tangent plane:

     x(i) = xq * t
     y(i) = yq * t

  enddo

end subroutine ll_gn_mult

!============================================================================

subroutine get_sincos_latlon(coslon,sinlon,coslat,sinlat,xe,ye,ze)

  use consts_coms, only: eradi

  implicit none

  real, intent(in)  :: xe, ye, ze   ! Earth cartesion distances
  real, intent(out) :: coslon, sinlon, coslat, sinlat
  real              :: ra

  ra = sqrt( xe**2 + ye**2 )

  sinlat = ze * eradi
  coslat = ra * eradi

  ! For points less than 100 m from Earth's polar axis, make arbitrary
  ! assumption that longitude = 0 deg.  This is just to settle on a PS
  ! planar coordinate system in which to do the algebra.

  if (ra >= 1.e2) then
     sinlon = ye / ra
     coslon = xe / ra
  else
     sinlon = 0.
     coslon = 1.
  endif

end subroutine get_sincos_latlon

!============================================================================
