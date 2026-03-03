module map_proj

  use map_proj_gn
  use map_proj_or
  use map_proj_ps

  ! Convert an earth cartesion (xe,ye,ze) point to geographic lat/lon
  interface ec_ll
     module procedure ec_ll_sing
     module procedure ec_ll_1d
     module procedure ec_ll_2d
  end interface ec_ll

  ! Convert a geographic lat/lon point to earth cartesion (xe,ye,ze) ccordinates
  interface ll_ec
     module procedure ll_ec_sing
     module procedure ll_ec_1d
     module procedure ll_ec_2d
  end interface ll_ec

  ! Compute the sines/cosines of geographic lat/lon from earth cartesian
  ! (xe,ye,ze) coordinates
  interface get_cossin_lonlat
     module procedure get_cossin_lonlat_sing
     module procedure get_cossin_lonlat_1d
     module procedure get_cossin_lonlat_2d
  end interface get_cossin_lonlat

contains

!============================================================================

subroutine ec_ll_sing(xeq,yeq,zeq,rlon,rlat)

  use consts_coms, only: piu180
  implicit none

  real, intent(in)  :: xeq
  real, intent(in)  :: yeq
  real, intent(in)  :: zeq
  real, intent(out) :: rlon
  real, intent(out) :: rlat
  real              :: rax

  ! This subroutine computes coordinates (rlon,rlat) of point q
  ! projected onto an equidistant cylindrical surface.  Input coordinates
  ! (xeq,yeq,zeq) of q are defined in "earth cartesian space", where the origin
  ! is the center of the earth, the z axis is the north pole, the x axis is the
  ! equator and prime meridian, and the y axis is the equator and 90 E.

  rax  = sqrt( xeq**2 + yeq**2 )  ! distance of q from earth axis
  rlat = atan2(zeq,rax) * piu180
  rlon = atan2(yeq,xeq) * piu180

end subroutine ec_ll_sing

!============================================================================

subroutine ec_ll_1d(xeq,yeq,zeq,rlon,rlat,n)

  use consts_coms, only: piu180
  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: xeq (n)
  real,    intent(in)  :: yeq (n)
  real,    intent(in)  :: zeq (n)
  real,    intent(out) :: rlon(n)
  real,    intent(out) :: rlat(n)
  real                 :: rax
  integer              :: i

  ! This subroutine computes coordinates (rlon,rlat) of point q
  ! projected onto an equidistant cylindrical surface.  Input coordinates
  ! (xeq,yeq,zeq) of q are defined in "earth cartesian space", where the origin
  ! is the center of the earth, the z axis is the north pole, the x axis is the
  ! equator and prime meridian, and the y axis is the equator and 90 E.

  do i = 1, n
     rax     = sqrt( xeq(i)**2 + yeq(i)**2 )  ! distance of q from earth axis
     rlat(i) = atan2(zeq(i),rax)    * piu180
     rlon(i) = atan2(yeq(i),xeq(i)) * piu180
  enddo

end subroutine ec_ll_1d

!============================================================================

subroutine ec_ll_2d(xeq,yeq,zeq,rlon,rlat,n1,n2)

  use consts_coms, only: piu180
  implicit none

  integer, intent(in)  :: n1, n2
  real,    intent(in)  :: xeq (n1,n2)
  real,    intent(in)  :: yeq (n1,n2)
  real,    intent(in)  :: zeq (n1,n2)
  real,    intent(out) :: rlon(n1,n2)
  real,    intent(out) :: rlat(n1,n2)
  real                 :: rax
  integer              :: i, j

  ! This subroutine computes coordinates (rlon,rlat) of point q
  ! projected onto an equidistant cylindrical surface.  Input coordinates
  ! (xeq,yeq,zeq) of q are defined in "earth cartesian space", where the origin
  ! is the center of the earth, the z axis is the north pole, the x axis is the
  ! equator and prime meridian, and the y axis is the equator and 90 E.

  !$omp parallel do private(i,rax)
  do j = 1, n2
     do i = 1, n1
        rax       = sqrt( xeq(i,j)**2 + yeq(i,j)**2 )  ! distance of q from earth axis
        rlat(i,j) = atan2(zeq(i,j),rax)      * piu180
        rlon(i,j) = atan2(yeq(i,j),xeq(i,j)) * piu180
     enddo
  enddo
  !$omp end parallel do

end subroutine ec_ll_2d

!============================================================================

subroutine ll_ec_sing(rlon,rlat,xeq,yeq,zeq)

  use consts_coms, only: pio180, erad
  implicit none

  real, intent(out) :: xeq
  real, intent(out) :: yeq
  real, intent(out) :: zeq
  real, intent(in)  :: rlon
  real, intent(in)  :: rlat
  real              :: sinrlat
  real              :: cosrlat
  real              :: sinrlon
  real              :: cosrlon

  ! This subroutine computes "earth cartesian space" coordinates
  ! (xeq,yeq,zeq) or a point located at (rlon,rlat).

  sinrlat = sin(rlat * pio180)
  cosrlat = cos(rlat * pio180)
  sinrlon = sin(rlon * pio180)
  cosrlon = cos(rlon * pio180)

  xeq = cosrlon * cosrlat * erad
  yeq = sinrlon * cosrlat * erad
  zeq =           sinrlat * erad

end subroutine ll_ec_sing

!============================================================================

subroutine ll_ec_1d(rlon,rlat,xeq,yeq,zeq,n)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: rlon(n)
  real,    intent(in)  :: rlat(n)
  real,    intent(out) :: xeq (n)
  real,    intent(out) :: yeq (n)
  real,    intent(out) :: zeq (n)
  real                 :: sinrlat
  real                 :: cosrlat
  real                 :: sinrlon
  real                 :: cosrlon
  integer              :: i

  ! This subroutine computes "earth cartesian space" coordinates
  ! (xeq,yeq,zeq) or a point located at (rlon,rlat).

  do i = 1, n
     sinrlat = sin(rlat(i) * pio180)
     cosrlat = cos(rlat(i) * pio180)
     sinrlon = sin(rlon(i) * pio180)
     cosrlon = cos(rlon(i) * pio180)

     xeq(i) = cosrlon * cosrlat * erad
     yeq(i) = sinrlon * cosrlat * erad
     zeq(i) =           sinrlat * erad
  enddo

end subroutine ll_ec_1d

!============================================================================

subroutine ll_ec_2d(rlon,rlat,xeq,yeq,zeq,n1,n2)

  use consts_coms, only: pio180, erad
  implicit none

  integer, intent(in)  :: n1,n2
  real,    intent(in)  :: rlon(n1,n2)
  real,    intent(in)  :: rlat(n1,n2)
  real,    intent(out) :: xeq (n1,n2)
  real,    intent(out) :: yeq (n1,n2)
  real,    intent(out) :: zeq (n1,n2)
  real                 :: sinrlat
  real                 :: cosrlat
  real                 :: sinrlon
  real                 :: cosrlon
  integer              :: i, j

  ! This subroutine computes "earth cartesian space" coordinates
  ! (xeq,yeq,zeq) or a point located at (rlon,rlat).

  !$omp parallel do private(i,sinrlat,cosrlat,sinrlon,cosrlon)
  do j = 1, n2
     do i = 1, n1
        sinrlat = sin(rlat(i,j) * pio180)
        cosrlat = cos(rlat(i,j) * pio180)
        sinrlon = sin(rlon(i,j) * pio180)
        cosrlon = cos(rlon(i,j) * pio180)

        xeq(i,j) = cosrlon * cosrlat * erad
        yeq(i,j) = sinrlon * cosrlat * erad
        zeq(i,j) =           sinrlat * erad
     enddo
  enddo
  !$omp end parallel do

end subroutine ll_ec_2d

!============================================================================

subroutine get_cossin_lonlat_sing(coslon,sinlon,coslat,sinlat,xe,ye,ze)

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

end subroutine get_cossin_lonlat_sing

!============================================================================

subroutine get_cossin_lonlat_1d(coslon,sinlon,coslat,sinlat,xe,ye,ze,n)

  use consts_coms, only: eradi
  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: xe(n), ye(n), ze(n)  ! Earth cartesion distances
  real,    intent(out) :: coslon(n), sinlon(n), coslat(n), sinlat(n)
  real                 :: ra
  integer              :: i

  do i = 1, n

     ra = sqrt( xe(i)**2 + ye(i)**2 )

     sinlat(i) = ze(i) * eradi
     coslat(i) = ra   * eradi

     sinlon(i) = ye(i) / max(ra,1.e2)
     coslon(i) = xe(i) / max(ra,1.e2)

     ! For points less than 100 m from Earth's polar axis, make arbitrary
     ! assumption that longitude = 0 deg.  This is just to settle on a PS
     ! planar coordinate system in which to do the algebra.

     if (ra < 1.e2) then
        sinlon(i) = 0.
        coslon(i) = 1.
     endif

  enddo

end subroutine get_cossin_lonlat_1d

!============================================================================

subroutine get_cossin_lonlat_2d(coslon,sinlon,coslat,sinlat,xe,ye,ze,n1,n2)

  use consts_coms, only: eradi
  implicit none

  integer, intent(in)  :: n1,n2
  real,    intent(in)  :: xe(n1,n2), ye(n1,n2), ze(n1,n2)  ! Earth cartesion distances
  real,    intent(out) :: coslon(n1,n2), sinlon(n1,n2), coslat(n1,n2), sinlat(n1,n2)
  real                 :: ra
  integer              :: i, j

  !$omp parallel do private(i,ra)
  do j = 1, n2
     do i = 1, n1

        ra = sqrt( xe(i,j)**2 + ye(i,j)**2 )

        sinlat(i,j) = ze(i,j) * eradi
        coslat(i,j) = ra      * eradi

        sinlon(i,j) = ye(i,j) / max(ra,1.e2)
        coslon(i,j) = xe(i,j) / max(ra,1.e2)

        ! For points less than 100 m from Earth's polar axis, make arbitrary
        ! assumption that longitude = 0 deg.  This is just to settle on a PS
        ! planar coordinate system in which to do the algebra.

        if (ra < 1.e2) then
           sinlon(i,j) = 0.
           coslon(i,j) = 1.
        endif

     enddo
  enddo
  !$omp end parallel do

end subroutine get_cossin_lonlat_2d

!============================================================================

end module map_proj
