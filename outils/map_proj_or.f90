module map_proj_or

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Routines for converting TO an orthographic projection:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! These routines compute coordinates (x,y) of point q projected onto a
  ! orthographic tangent plane whose pole point is located at geographic coordinates
  ! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
  ! "earth cartesian space", where the origin is the center of the earth,
  ! the z axis is the north pole, the x axis is the equator and prime meridian,
  ! and the y axis is the equator and 90 E.
  interface ec_or
     module procedure ec_or_sing
     module procedure ec_or_1d
     module procedure ec_or_2d
  end interface ec_or

  ! These routine computes coordinates (x,y) of point q located at (glat,glon)
  ! projected onto an orthographic tangent plane whose pole point is located at
  ! geographic coordinates (polelat,polelon).
  interface ll_or
     module procedure ll_or_sing
     module procedure ll_or_1d
     module procedure ll_or_2d
  end interface ll_or

  ! These routines compute coordinates (x,y) of point q projected onto a
  ! orthographic tangent plane given the sines and cosines of the lat and lon
  ! of the projection's pole point and the vector distance of point q
  ! (in earth cartesian coordinates) from the pole point.
  interface de_or
     module procedure de_or_sing
     module procedure de_or_1d
     module procedure de_or_2d
  end interface de_or


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Routines for converting FROM an orthographic projection:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Given coordinates (x,y) of point q projected onto an orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! these routines compute coordinates (xeq,yeq,zeq) of q on earth spherical
  ! surface defined in "earth cartesian space"
  interface or_ec
     module procedure or_ec_sing
     module procedure or_ec_1d
     module procedure or_ec_2d
     module procedure or_ec_2d_xygrid
  end interface or_ec

  ! Given coordinates (x,y) of point q projected onto an orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! these routines compute geographic coordinate (qlat,qlon) of q
  interface or_ll
     module procedure or_ll_sing
     module procedure or_ll_1d
     module procedure or_ll_2d
     module procedure or_ll_2d_xygrid
  end interface or_ll

  ! Given coordinates (x,y) of point q projected onto an orthographic tangent plane
  ! and the sines and cosines of the lat and lon of the tangent plane's pole
  ! point, compute the vector distance of q (in earth cartesian space) from the
  ! pole point.
  interface or_de
     module procedure or_de_sing
     module procedure or_de_1d
     module procedure or_de_2d
     module procedure or_de_2d_xygrid
  end interface or_de


contains

!============================================================================

subroutine ec_or_sing(xeq,yeq,zeq,polelat,polelon,x,y)

  use consts_coms, only: pio180, erad
  implicit none

  real, intent(in)  :: xeq
  real, intent(in)  :: yeq
  real, intent(in)  :: zeq
  real, intent(in)  :: polelat
  real, intent(in)  :: polelon
  real, intent(out) :: x
  real, intent(out) :: y

  real :: sinplat, cosplat
  real :: sinplon, cosplon
  real :: xep, yep, zep
  real :: dxe, dye, dze
  real :: zq

  ! This subroutine computes coordinates (x,y) of point q projected onto
  ! orthographic tangent plane whose pole point is located at geographic coordinates
  ! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
  ! "earth cartesian space", where the origin is the center of the earth,
  ! the z axis is the north pole, the x axis is the equator and prime meridian,
  ! and the y axis is the equator and 90 E.

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
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  dxe = xeq - xep
  dye = yeq - yep
  dze = zeq - zep

  x =                          - sinplon * dxe + cosplon * dye
  y = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

  ! Orthographic projection can only represent one hemisphere; set x and y
  ! to a large value to indicate the projection is out of bounds

  zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
  if (zq < -.999*erad) then
     x = 1.e30
     y = 1.e30
  endif

end subroutine ec_or_sing

!============================================================================

subroutine ec_or_1d(xeq,yeq,zeq,polelat,polelon,x,y,n)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: xeq(n)
  real,    intent(in)  :: yeq(n)
  real,    intent(in)  :: zeq(n)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(out) :: x(n)
  real,    intent(out) :: y(n)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, yep, zep
  real    :: dxe, dye, dze
  real    :: zq
  integer :: i

  ! This subroutine computes coordinates (x,y) of point q projected onto
  ! orthographic tangent plane whose pole point is located at geographic coordinates
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
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  do i = 1, n

     dxe = xeq(i) - xep
     dye = yeq(i) - yep
     dze = zeq(i) - zep

     x(i) =                          - sinplon * dxe + cosplon * dye
     y(i) = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

     ! Orthographic projection can only represent one hemisphere; set x and y
     ! to a large value to indicate the projection is out of bounds

     zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
     if (zq < -.999*erad) then
        x(i) = 1.e30
        y(i) = 1.e30
     endif

  enddo

end subroutine ec_or_1d

!============================================================================

subroutine ec_or_2d(xeq,yeq,zeq,polelat,polelon,x,y,n1,n2)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(in)  :: xeq(n1,n2)
  real,    intent(in)  :: yeq(n1,n2)
  real,    intent(in)  :: zeq(n1,n2)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(out) :: x(n1,n2)
  real,    intent(out) :: y(n1,n2)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, yep, zep
  real    :: dxe, dye, dze
  real    :: zq
  integer :: i, j

  ! This subroutine computes coordinates (x,y) of point q projected onto
  ! orthographic tangent plane whose pole point is located at geographic coordinates
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
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  !$omp parallel do private(i,dxe,dye,dze,zq)
  do j = 1, n2
     do i = 1, n1

        dxe = xeq(i,j) - xep
        dye = yeq(i,j) - yep
        dze = zeq(i,j) - zep

        x(i,j) =                          - sinplon * dxe + cosplon * dye
        y(i,j) = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

        ! Orthographic projection can only represent one hemisphere; set x and y
        ! to a large value to indicate the projection is out of bounds

        zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
        if (zq < -.999*erad) then
           x(i,j) = 1.e30
           y(i,j) = 1.e30
        endif

     enddo
  enddo
  !$omp end parallel do

end subroutine ec_or_2d

!============================================================================

subroutine ll_or_sing(qlat,qlon,polelat,polelon,x,y)

  use consts_coms, only: pio180, erad
  implicit none

  real, intent(in)  :: qlat
  real, intent(in)  :: qlon
  real, intent(in)  :: polelat
  real, intent(in)  :: polelon
  real, intent(out) :: x
  real, intent(out) :: y

  real :: sinplat, sinqlat
  real :: cosplat, cosqlat
  real :: sinplon, sinqlon
  real :: cosplon, cosqlon
  real :: xep, xeq, dxe
  real :: yep, yeq, dye
  real :: zep, zeq, dze
  real :: zq

  ! This subroutine computes coordinates (x,y) of a point (qlat,qlon) projected
  ! onto an orthographic tangent plane whose pole point is located at geographic
  ! coordinates (polelat,polelon).

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

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  ! Compute (xeq,yeq,zeq) coordinates of the input point in "earth cartesian space"

  xeq = erad * cosqlat * cosqlon
  yeq = erad * cosqlat * sinqlon
  zeq = erad * sinqlat

  ! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  dxe = xeq - xep
  dye = yeq - yep
  dze = zeq - zep

  x =                          - sinplon * dxe + cosplon * dye
  y = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

  ! Orthographic projection can only represent one hemisphere; set x and y
  ! to a large value to indicate the projection is out of bounds

  zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
  if (zq < -.999*erad) then
     x = 1.e30
     y = 1.e30
  endif

end subroutine ll_or_sing

!============================================================================

subroutine ll_or_1d(qlat,qlon,polelat,polelon,x,y,n)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: qlat(n)
  real,    intent(in)  :: qlon(n)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(out) :: x(n)
  real,    intent(out) :: y(n)

  real    :: sinplat, sinqlat
  real    :: cosplat, cosqlat
  real    :: sinplon, sinqlon
  real    :: cosplon, cosqlon
  real    :: xep, xeq, dxe
  real    :: yep, yeq, dye
  real    :: zep, zeq, dze
  real    :: zq
  integer :: i

  ! This subroutine computes coordinates (x,y) of a point (qlat,qlon) projected
  ! onto an orthographic tangent plane whose pole point is located at geographic
  ! coordinates (polelat,polelon).

  ! Evaluate sine and cosine of latitude and longitude of pole point p

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  do i = 1, n

     ! Evaluate sine and cosine of latitude and longitude of input point q

     sinqlat = sin(qlat(i) * pio180)
     cosqlat = cos(qlat(i) * pio180)
     sinqlon = sin(qlon(i) * pio180)
     cosqlon = cos(qlon(i) * pio180)

     ! Compute (xeq,yeq,zeq) coordinates of the input point in "earth cartesian space"

     xeq = erad * cosqlat * cosqlon
     yeq = erad * cosqlat * sinqlon
     zeq = erad * sinqlat

     ! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
     ! orthographic tangent plane with origin at the pole point, the z axis pointing
     ! radially outward from the center of the earth, and the y axis pointing
     ! northward along the local earth meridian from the pole point.

     dxe = xeq - xep
     dye = yeq - yep
     dze = zeq - zep

     x(i) =                          - sinplon * dxe + cosplon * dye
     y(i) = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

     ! Orthographic projection can only represent one hemisphere; set x and y
     ! to a large value to indicate the projection is out of bounds

     zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
     if (zq < -.999*erad) then
        x(i) = 1.e30
        y(i) = 1.e30
     endif

  enddo

end subroutine ll_or_1d

!============================================================================

subroutine ll_or_2d(qlat,qlon,polelat,polelon,x,y,n1,n2)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(in)  :: qlat(n1,n2)
  real,    intent(in)  :: qlon(n1,n2)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(out) :: x(n1,n2)
  real,    intent(out) :: y(n1,n2)

  real    :: sinplat, sinqlat
  real    :: cosplat, cosqlat
  real    :: sinplon, sinqlon
  real    :: cosplon, cosqlon
  real    :: xep, xeq, dxe
  real    :: yep, yeq, dye
  real    :: zep, zeq, dze
  real    :: zq
  integer :: i, j

  ! This subroutine computes coordinates (x,y) of a point (qlat,qlon) projected
  ! onto an orthographic tangent plane whose pole point is located at geographic
  ! coordinates (polelat,polelon).

  ! Evaluate sine and cosine of latitude and longitude of pole point p

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  !$omp parallel do private(i,sinqlat,cosqlat,sinqlon,cosqlon,xeq,yeq,zeq, &
  !$omp                     dxe,dye,dze,zq)
  do j = 1, n2
     do i = 1, n1

        ! Evaluate sine and cosine of latitude and longitude of input point q

        sinqlat = sin(qlat(i,j) * pio180)
        cosqlat = cos(qlat(i,j) * pio180)
        sinqlon = sin(qlon(i,j) * pio180)
        cosqlon = cos(qlon(i,j) * pio180)

        ! Compute (xeq,yeq,zeq) coordinates of the input point in "earth cartesian space"

        xeq = erad * cosqlat * cosqlon
        yeq = erad * cosqlat * sinqlon
        zeq = erad * sinqlat

        ! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
        ! orthographic tangent plane with origin at the pole point, the z axis pointing
        ! radially outward from the center of the earth, and the y axis pointing
        ! northward along the local earth meridian from the pole point.

        dxe = xeq - xep
        dye = yeq - yep
        dze = zeq - zep

        x(i,j) =                          - sinplon * dxe + cosplon * dye
        y(i,j) = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

        ! Orthographic projection can only represent one hemisphere; set x and y
        ! to a large value to indicate the projection is out of bounds

        zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
        if (zq < -.999*erad) then
           x(i,j) = 1.e30
           y(i,j) = 1.e30
        endif

     enddo
  enddo
  !$omp end parallel do

end subroutine ll_or_2d

!============================================================================

subroutine de_or_sing(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y)

  use consts_coms, only: pio180, erad
  implicit none

  real, intent(in)  :: dxe
  real, intent(in)  :: dye
  real, intent(in)  :: dze
  real, intent(in)  :: cosplat
  real, intent(in)  :: sinplat
  real, intent(in)  :: cosplon
  real, intent(in)  :: sinplon
  real, intent(out) :: x
  real, intent(out) :: y

  real :: zq

  ! This subroutine computes coordinates (x,y) of point q projected onto
  ! orthographic tangent plane whose pole point is located at geographic coordinates
  ! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
  ! "earth cartesian space", where the origin is the center of the earth,
  ! the z axis is the north pole, the x axis is the equator and prime meridian,
  ! and the y axis is the equator and 90 E..

  ! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  x =                          - sinplon * dxe + cosplon * dye
  y = cosplat * dze - sinplat * (cosplon * dxe + sinplon * dye)

  ! Orthographic projection can only represent one hemisphere; set x and y
  ! to a large value to indicate the projection is out of bounds

  zq = sinplat * dze + cosplat * (cosplon * dxe + sinplon * dye)
  if (zq < -.999*erad) then
     x = 1.e30
     y = 1.e30
  endif

end subroutine de_or_sing

!============================================================================

subroutine de_or_1d(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y,n)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: dxe(n)
  real,    intent(in)  :: dye(n)
  real,    intent(in)  :: dze(n)
  real,    intent(in)  :: cosplat
  real,    intent(in)  :: sinplat
  real,    intent(in)  :: cosplon
  real,    intent(in)  :: sinplon
  real,    intent(out) :: x(n)
  real,    intent(out) :: y(n)

  real    :: zq
  integer :: i

  ! This subroutine computes coordinates (x,y) of point q projected onto
  ! orthographic tangent plane whose pole point is located at geographic coordinates
  ! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
  ! "earth cartesian space", where the origin is the center of the earth,
  ! the z axis is the north pole, the x axis is the equator and prime meridian,
  ! and the y axis is the equator and 90 E..

  ! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  do i = 1, n

     x(i) =                             - sinplon * dxe(i) + cosplon * dye(i)
     y(i) = cosplat * dze(i) - sinplat * (cosplon * dxe(i) + sinplon * dye(i))

     ! Orthographic projection can only represent one hemisphere; set x and y
     ! to a large value to indicate the projection is out of bounds

     zq = sinplat * dze(i) + cosplat * (cosplon * dxe(i) + sinplon * dye(i))
     if (zq < -.999*erad) then
        x(i) = 1.e30
        y(i) = 1.e30
     endif

  enddo

end subroutine de_or_1d

!============================================================================

subroutine de_or_2d(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y,n1,n2)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n1,n2
  real,    intent(in)  :: dxe(n1,n2)
  real,    intent(in)  :: dye(n1,n2)
  real,    intent(in)  :: dze(n1,n2)
  real,    intent(in)  :: cosplat
  real,    intent(in)  :: sinplat
  real,    intent(in)  :: cosplon
  real,    intent(in)  :: sinplon
  real,    intent(out) :: x(n1,n2)
  real,    intent(out) :: y(n1,n2)

  real    :: zq
  integer :: i, j

  ! This subroutine computes coordinates (x,y) of point q projected onto
  ! orthographic tangent plane whose pole point is located at geographic coordinates
  ! (polelat,polelon).  Input coordinates (xeq,yeq,zeq) of q are defined in
  ! "earth cartesian space", where the origin is the center of the earth,
  ! the z axis is the north pole, the x axis is the equator and prime meridian,
  ! and the y axis is the equator and 90 E..

  ! Transform q point from (xe,ye,ze) coordinates to 3D coordinates relative to
  ! orthographic tangent plane with origin at the pole point, the z axis pointing
  ! radially outward from the center of the earth, and the y axis pointing
  ! northward along the local earth meridian from the pole point.

  !$omp parallel do private(i,zq)
  do j = 1, n2
     do i = 1, n1

        x(i,j) =                               - sinplon * dxe(i,j) + cosplon * dye(i,j)
        y(i,j) = cosplat * dze(i,j) - sinplat * (cosplon * dxe(i,j) + sinplon * dye(i,j))

        ! Orthographic projection can only represent one hemisphere; set x and y
        ! to a large value to indicate the projection is out of bounds

        zq = sinplat * dze(i,j) + cosplat * (cosplon * dxe(i,j) + sinplon * dye(i,j))
        if (zq < -.999*erad) then
           x(i,j) = 1.e30
           y(i,j) = 1.e30
        endif

     enddo
  enddo
  !$omp end parallel do

end subroutine de_or_2d

!============================================================================

subroutine or_ec_sing(xeq,yeq,zeq,polelat,polelon,x,y,iskip)

  use consts_coms, only: pio180, erad, eradsq
  implicit none

  real, intent(out) :: xeq
  real, intent(out) :: yeq
  real, intent(out) :: zeq
  real, intent(in)  :: polelat
  real, intent(in)  :: polelon
  real, intent(in)  :: x
  real, intent(in)  :: y

  logical, optional, intent(out) :: iskip

  real :: sinplat, cosplat
  real :: sinplon, cosplon
  real :: xep, yep, zep
  real :: xq, yq, zq

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
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

  ! This gives the following xq, yq, zq coordinates relative to the
  ! orthographic plane for the projection of point q from the orthographic
  ! plane to the earth surface:

  xq = x
  yq = y
  zq = sqrt( max(0., eradsq - x*x - y*y) ) - erad

  if (present(iskip)) then
     iskip = ( x*x + y*y > .999 * eradsq )
  endif

  ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
  ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
  ! at the pole point, with the z axis pointing radially outward from the center
  ! of the earth, and the y axis pointing northward along the local earth meridian
  ! from the pole point.

  xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
  yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
  zeq = zep                            + cosplat * yq + sinplat * zq

end subroutine or_ec_sing

!============================================================================

subroutine or_ec_1d(xeq,yeq,zeq,polelat,polelon,x,y,n,iskip)

  use consts_coms, only: pio180, erad, eradsq
  implicit none

  integer, intent(in)  :: n
  real,    intent(out) :: xeq(n)
  real,    intent(out) :: yeq(n)
  real,    intent(out) :: zeq(n)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(in)  :: x(n)
  real,    intent(in)  :: y(n)

  logical, optional, intent(out) :: iskip(n)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, yep, zep
  real    :: xq, yq, zq
  integer :: i

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
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

  do i = 1, n

     ! This gives the following xq, yq, zq coordinates relative to the
     ! orthographic plane for the projection of point q from the orthographic
     ! plane to the earth surface:

     xq = x(i)
     yq = y(i)
     zq = sqrt( max(0., eradsq - x(i)*x(i) - y(i)*y(i)) ) - erad

     if (present(iskip)) then
        iskip(i) = ( x(i)*x(i) + y(i)*y(i) > .999 * eradsq )
     endif

     ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
     ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
     ! at the pole point, with the z axis pointing radially outward from the center
     ! of the earth, and the y axis pointing northward along the local earth meridian
     ! from the pole point.

     xeq(i) = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
     yeq(i) = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
     zeq(i) = zep                            + cosplat * yq + sinplat * zq

  enddo

end subroutine or_ec_1d

!============================================================================

subroutine or_ec_2d(xeq,yeq,zeq,polelat,polelon,x,y,n1,n2,iskip)

  use consts_coms, only: pio180, erad, eradsq
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(out) :: xeq(n1,n2)
  real,    intent(out) :: yeq(n1,n2)
  real,    intent(out) :: zeq(n1,n2)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(in)  :: x(n1,n2)
  real,    intent(in)  :: y(n1,n2)

  logical, optional, intent(out) :: iskip(n1,n2)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, yep, zep
  real    :: xq, yq, zq
  integer :: i, j

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
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

  !$omp parallel do private(i,xq,yq,zq)
  do j = 1, n2
     do i = 1, n1

        ! This gives the following xq, yq, zq coordinates relative to the
        ! orthographic plane for the projection of point q from the orthographic
        ! plane to the earth surface:

        xq = x(i,j)
        yq = y(i,j)
        zq = sqrt( max(0., eradsq - x(i,j)*x(i,j) - y(i,j)*y(i,j)) ) - erad

        if (present(iskip)) then
           iskip(i,j) = ( x(i,j)*x(i,j) + y(i,j)*y(i,j) > .999 * eradsq )
        endif

        ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
        ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
        ! at the pole point, with the z axis pointing radially outward from the center
        ! of the earth, and the y axis pointing northward along the local earth meridian
        ! from the pole point.

        xeq(i,j) = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
        yeq(i,j) = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
        zeq(i,j) = zep                            + cosplat * yq + sinplat * zq

     enddo
  enddo
  !$omp end parallel do

end subroutine or_ec_2d

!============================================================================

subroutine or_ec_2d_xygrid(xeq,yeq,zeq,polelat,polelon,x,y,n1,n2,iskip)

  use consts_coms, only: pio180, erad, eradsq
  implicit none

  integer, intent(in)    :: n1, n2
  real,    intent(inout) :: xeq(n1,n2)
  real,    intent(inout) :: yeq(n1,n2)
  real,    intent(inout) :: zeq(n1,n2)
  real,    intent(in)    :: polelat
  real,    intent(in)    :: polelon
  real,    intent(in)    :: x(n1)
  real,    intent(in)    :: y(n2)

  logical, optional, intent(out) :: iskip(n1,n2)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, yep, zep
  real    :: xq, yq, zq
  real    :: eeyy
  integer :: i, j

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
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

  !$omp parallel do private(eeyy,i,xq,yq,zq)
  do j = 1, n2

     eeyy = eradsq - y(j) * y(j)
     yq   = y(j)

     do i = 1, n1

        ! This gives the following xq, yq, zq coordinates relative to the
        ! orthographic plane for the projection of point q from the orthographic
        ! plane to the earth surface:

        xq = x(i)
        zq = sqrt( max(0., eeyy - x(i)*x(i)) ) - erad

        if (present(iskip)) then
           iskip(i,j) = ( x(i)*x(i) + y(j)*y(j) > .999 * eradsq )
        endif

        ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
        ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
        ! at the pole point, with the z axis pointing radially outward from the center
        ! of the earth, and the y axis pointing northward along the local earth meridian
        ! from the pole point.

        xeq(i,j) = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
        yeq(i,j) = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
        zeq(i,j) = zep                            + cosplat * yq + sinplat * zq

     enddo
  enddo
  !$omp end parallel do

end subroutine or_ec_2d_xygrid

!============================================================================

subroutine or_ll_sing(qlat,qlon,polelat,polelon,x,y,iskip)

  use consts_coms, only: pio180, erad, eradsq, piu180
  implicit none

  real, intent(out) :: qlat
  real, intent(out) :: qlon
  real, intent(in)  :: polelat
  real, intent(in)  :: polelon
  real, intent(in)  :: x
  real, intent(in)  :: y

  logical, optional, intent(out) :: iskip

  real :: sinplat, cosplat
  real :: sinplon, cosplon
  real :: xep, xeq
  real :: yep, yeq
  real :: zep, zeq
  real :: xq, yq, zq
  real :: req

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes geographic coordinates (qlat,qlon) of q.

  ! Evaluate sine and cosine of latitude and longitude of pole point p

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  ! This gives the following xq, yq, zq coordinates relative to the
  ! orthographic plane for the projection of point q from the orthographic
  ! plane to the earth surface:

  xq = x
  yq = y
  zq = sqrt( max(0., eradsq - x*x - y*y) ) - erad

  if (present(iskip)) then
     iskip = ( x*x + y*y > .999 * eradsq )
  endif

  ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
  ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
  ! at the pole point, with the z axis pointing radially outward from the center
  ! of the earth, and the y axis pointing northward along the local earth meridian
  ! from the pole point.

  xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
  yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
  zeq = zep                            + cosplat * yq + sinplat * zq

  ! Compute the latitude and longitude of Q:

  req = sqrt(xeq**2 + yeq**2)
  qlon = atan2(yeq,xeq) * piu180
  qlat = atan2(zeq,req) * piu180

end subroutine or_ll_sing

!============================================================================

subroutine or_ll_1d(qlat,qlon,polelat,polelon,x,y,n,iskip)

  use consts_coms, only: pio180, erad, eradsq, piu180
  implicit none

  integer, intent(in)  :: n
  real,    intent(out) :: qlat(n)
  real,    intent(out) :: qlon(n)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(in)  :: x(n)
  real,    intent(in)  :: y(n)

  logical, optional, intent(out) :: iskip(n)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, xeq
  real    :: yep, yeq
  real    :: zep, zeq
  real    :: xq, yq, zq
  real    :: req
  integer :: i

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes geographic coordinates (qlat,qlon) of q.

  ! Evaluate sine and cosine of latitude and longitude of pole point p

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  do i = 1, n

     ! This gives the following xq, yq, zq coordinates relative to the
     ! orthographic plane for the projection of point q from the orthographic
     ! plane to the earth surface:

     xq = x(i)
     yq = y(i)
     zq = sqrt( max(0., eradsq - x(i)*x(i) - y(i)*y(i)) ) - erad

     if (present(iskip)) then
        iskip(i) = ( x(i)*x(i) + y(i)*y(i) > .999 * eradsq )
     endif

     ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
     ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
     ! at the pole point, with the z axis pointing radially outward from the center
     ! of the earth, and the y axis pointing northward along the local earth meridian
     ! from the pole point.

     xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
     yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
     zeq = zep                            + cosplat * yq + sinplat * zq

     ! Compute the latitude and longitude of Q:

     req    = sqrt(xeq**2 + yeq**2)
     qlon(i) = atan2(yeq,xeq) * piu180
     qlat(i) = atan2(zeq,req) * piu180

  enddo

end subroutine or_ll_1d

!============================================================================

subroutine or_ll_2d(qlat,qlon,polelat,polelon,x,y,n1,n2,iskip)

  use consts_coms, only: pio180, erad, eradsq, piu180
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(out) :: qlat(n1,n2)
  real,    intent(out) :: qlon(n1,n2)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(in)  :: x(n1,n2)
  real,    intent(in)  :: y(n1,n2)

  logical, optional, intent(out) :: iskip(n1,n2)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, xeq
  real    :: yep, yeq
  real    :: zep, zeq
  real    :: xq, yq, zq
  real    :: req
  integer :: i, j

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes geographic coordinates (qlat,qlon) of q.

  ! Evaluate sine and cosine of latitude and longitude of pole point p

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  !$omp parallel do private(i,xq,yq,zq,xeq,yeq,zeq,req)
  do j = 1, n2
     do i = 1, n1

        ! This gives the following xq, yq, zq coordinates relative to the
        ! orthographic plane for the projection of point q from the orthographic
        ! plane to the earth surface:

        xq = x(i,j)
        yq = y(i,j)
        zq = sqrt( max(0., eradsq - x(i,j)*x(i,j) - y(i,j)*y(i,j)) ) - erad

        if (present(iskip)) then
           iskip(i,j) = ( x(i,j)*x(i,j) + y(i,j)*y(i,j) > .999 * eradsq )
        endif

        ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
        ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
        ! at the pole point, with the z axis pointing radially outward from the center
        ! of the earth, and the y axis pointing northward along the local earth meridian
        ! from the pole point.

        xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
        yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
        zeq = zep                            + cosplat * yq + sinplat * zq

        ! Compute the latitude and longitude of Q:

        req       = sqrt(xeq**2 + yeq**2)
        qlon(i,j) = atan2(yeq,xeq) * piu180
        qlat(i,j) = atan2(zeq,req) * piu180

     enddo
  enddo
  !$omp end parallel do

end subroutine or_ll_2d

!============================================================================

subroutine or_ll_2d_xygrid(qlat,qlon,polelat,polelon,x,y,n1,n2,iskip)

  use consts_coms, only: pio180, erad, eradsq, piu180
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(out) :: qlat(n1,n2)
  real,    intent(out) :: qlon(n1,n2)
  real,    intent(in)  :: polelat
  real,    intent(in)  :: polelon
  real,    intent(in)  :: x(n1)
  real,    intent(in)  :: y(n2)

  logical, optional, intent(out) :: iskip(n1,n2)

  real    :: sinplat, cosplat
  real    :: sinplon, cosplon
  real    :: xep, xeq
  real    :: yep, yeq
  real    :: zep, zeq
  real    :: xq, yq, zq
  real    :: req
  real    :: eeyy
  integer :: i, j

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes geographic coordinates (qlat,qlon) of q.

  ! Evaluate sine and cosine of latitude and longitude of pole point p

  sinplat = sin(polelat * pio180)
  cosplat = cos(polelat * pio180)
  sinplon = sin(polelon * pio180)
  cosplon = cos(polelon * pio180)

  ! Compute (xep,yep,zep) coordinates of the pole point in "earth cartesian space"

  xep = erad * cosplat * cosplon
  yep = erad * cosplat * sinplon
  zep = erad * sinplat

  !$omp parallel do private(eeyy,i,xq,yq,zq,xeq,yeq,zeq,req)
  do j = 1, n2

     eeyy = eradsq - y(j) * y(j)
     yq   = y(j)

     do i = 1, n1

        ! This gives the following xq, yq, zq coordinates relative to the
        ! orthographic plane for the projection of point q from the orthographic
        ! plane to the earth surface:

        xq = x(i)
        zq = sqrt( max(0., eeyy - x(i)*x(i)) ) - erad

        if (present(iskip)) then
           iskip(i,j) = ( x(i)*x(i) + y(j)*y(j) > .999 * eradsq )
        endif

        ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
        ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
        ! at the pole point, with the z axis pointing radially outward from the center
        ! of the earth, and the y axis pointing northward along the local earth meridian
        ! from the pole point.

        xeq = xep - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
        yeq = yep + cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
        zeq = zep                            + cosplat * yq + sinplat * zq

        ! Compute the latitude and longitude of Q:

        req       = sqrt(xeq**2 + yeq**2)
        qlon(i,j) = atan2(yeq,xeq) * piu180
        qlat(i,j) = atan2(zeq,req) * piu180

     enddo
  enddo
  !$omp end parallel do

end subroutine or_ll_2d_xygrid

!============================================================================

subroutine or_de_sing(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y,iskip)

  use consts_coms, only: erad, eradsq
  implicit none

  real, intent(out) :: dxe
  real, intent(out) :: dye
  real, intent(out) :: dze
  real, intent(in)  :: cosplat
  real, intent(in)  :: sinplat
  real, intent(in)  :: cosplon
  real, intent(in)  :: sinplon
  real, intent(in)  :: x
  real, intent(in)  :: y

  logical, optional, intent(out) :: iskip

  real :: xq, yq, zq

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
  ! surface defined in "earth cartesian space", where the origin is the center
  ! of the earth, the z axis is the north pole, the x axis is the equator and
  ! prime meridian, and the y axis is the equator and 90 E..

  ! This gives the following xq, yq, zq coordinates relative to the
  ! orthographic plane for the projection of point q from the orthographic
  ! plane to the earth surface:

  xq = x
  yq = y
  zq = sqrt( max(0., eradsq - x*x - y*y) ) - erad

  if (present(iskip)) then
     iskip = ( x*x + y*y > .999 * eradsq )
  endif

  ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
  ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
  ! at the pole point, with the z axis pointing radially outward from the center
  ! of the earth, and the y axis pointing northward along the local earth meridian
  ! from the pole point.

  dxe = - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
  dye =   cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
  dze =                            + cosplat * yq + sinplat * zq

end subroutine or_de_sing

!============================================================================

subroutine or_de_1d(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y,n,iskip)

  use consts_coms, only: erad, eradsq
  implicit none

  integer, intent(in)  :: n
  real,    intent(out) :: dxe(n)
  real,    intent(out) :: dye(n)
  real,    intent(out) :: dze(n)
  real,    intent(in)  :: cosplat
  real,    intent(in)  :: sinplat
  real,    intent(in)  :: cosplon
  real,    intent(in)  :: sinplon
  real,    intent(in)  :: x(n)
  real,    intent(in)  :: y(n)

  logical, optional, intent(out) :: iskip(n)

  real    :: xq, yq, zq
  integer :: i

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
  ! surface defined in "earth cartesian space", where the origin is the center
  ! of the earth, the z axis is the north pole, the x axis is the equator and
  ! prime meridian, and the y axis is the equator and 90 E..

  do i = 1, n

     ! This gives the following xq, yq, zq coordinates relative to the
     ! orthographic plane for the projection of point q from the orthographic
     ! plane to the earth surface:

     xq = x(i)
     yq = y(i)
     zq = sqrt( max(0., eradsq - x(i)*x(i) - y(i)*y(i)) ) - erad

     if (present(iskip)) then
        iskip(i) = ( x(i)*x(i) + y(i)*y(i) > .999 * eradsq )
     endif

     ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
     ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
     ! at the pole point, with the z axis pointing radially outward from the center
     ! of the earth, and the y axis pointing northward along the local earth meridian
     ! from the pole point.

     dxe(i) = - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
     dye(i) =   cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
     dze(i) =                            + cosplat * yq + sinplat * zq

  enddo

end subroutine or_de_1d

!============================================================================

subroutine or_de_2d(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y,n1,n2,iskip)

  use consts_coms, only: erad, eradsq
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(out) :: dxe(n1,n2)
  real,    intent(out) :: dye(n1,n2)
  real,    intent(out) :: dze(n1,n2)
  real,    intent(in)  :: cosplat
  real,    intent(in)  :: sinplat
  real,    intent(in)  :: cosplon
  real,    intent(in)  :: sinplon
  real,    intent(in)  :: x(n1,n2)
  real,    intent(in)  :: y(n1,n2)

  logical, optional, intent(out) :: iskip(n1,n2)

  real    :: xq, yq, zq
  integer :: i, j

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
  ! surface defined in "earth cartesian space", where the origin is the center
  ! of the earth, the z axis is the north pole, the x axis is the equator and
  ! prime meridian, and the y axis is the equator and 90 E..

  !$omp parallel do private(i,xq,yq,zq)
  do j = 1, n2
     do i = 1, n1

        ! This gives the following xq, yq, zq coordinates relative to the
        ! orthographic plane for the projection of point q from the orthographic
        ! plane to the earth surface:

        xq = x(i,j)
        yq = y(i,j)
        zq = sqrt( max(0., eradsq - x(i,j)*x(i,j) - y(i,j)*y(i,j)) ) - erad

        if (present(iskip)) then
           iskip(i,j) = ( x(i,j)*x(i,j) + y(i,j)*y(i,j) > .999 * eradsq )
        endif

        ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
        ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
        ! at the pole point, with the z axis pointing radially outward from the center
        ! of the earth, and the y axis pointing northward along the local earth meridian
        ! from the pole point.

        dxe(i,j) = - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
        dye(i,j) =   cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
        dze(i,j) =                            + cosplat * yq + sinplat * zq

     enddo
  enddo
  !$omp end parallel do

end subroutine or_de_2d

!============================================================================

subroutine or_de_2d_xygrid(dxe,dye,dze,cosplat,sinplat,cosplon,sinplon,x,y,n1,n2,iskip)

  use consts_coms, only: erad, eradsq
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(out) :: dxe(n1,n2)
  real,    intent(out) :: dye(n1,n2)
  real,    intent(out) :: dze(n1,n2)
  real,    intent(in)  :: cosplat
  real,    intent(in)  :: sinplat
  real,    intent(in)  :: cosplon
  real,    intent(in)  :: sinplon
  real,    intent(in)  :: x(n1)
  real,    intent(in)  :: y(n2)

  logical, optional, intent(out) :: iskip(n1,n2)

  real    :: xq, yq, zq
  real    :: eeyy
  integer :: i, j

  ! Given coordinates (x,y) of point q projected onto orthographic tangent plane
  ! whose pole point is located at geographic coordinates (polelat,polelon),
  ! this subroutine computes coordinates (xeq,yeq,zeq) of q on earth spherical
  ! surface defined in "earth cartesian space", where the origin is the center
  ! of the earth, the z axis is the north pole, the x axis is the equator and
  ! prime meridian, and the y axis is the equator and 90 E..

  !$omp parallel do private(eeyy,i,xq,yq,zq)
  do j = 1, n2

     eeyy = eradsq - y(j) * y(j)
     yq   = y(j)

     do i = 1, n1

        ! This gives the following xq, yq, zq coordinates relative to the
        ! orthographic plane for the projection of point q from the orthographic
        ! plane to the earth surface:

        xq = x(i)
        zq = sqrt( max(0., eeyy - x(i)*x(i)) ) - erad

        if (present(iskip)) then
           iskip(i,j) = ( x(i)*x(i) + y(j)*y(j) > .999 * eradsq )
        endif

        ! Transform q point located on the earth's surface from or coordinates (xq,yq,zq)
        ! to earth coordinates (xe,ye,ze).  The orthographic plane has its origin
        ! at the pole point, with the z axis pointing radially outward from the center
        ! of the earth, and the y axis pointing northward along the local earth meridian
        ! from the pole point.

        dxe(i,j) = - sinplon * xq + cosplon * (-sinplat * yq + cosplat * zq)
        dye(i,j) =   cosplon * xq - sinplon * ( sinplat * yq - cosplat * zq)
        dze(i,j) =                            + cosplat * yq + sinplat * zq

     enddo
  enddo
  !$omp end parallel do

end subroutine or_de_2d_xygrid

!============================================================================

end module map_proj_or
