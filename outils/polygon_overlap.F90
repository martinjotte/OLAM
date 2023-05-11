subroutine polygon_overlap(np,nq,xp,yp,xq,yq,arp,arq,area)

  ! Given x,y coordinates of the vertices of polygons p and q, compute the area
  ! of overlap between the polygons using a sweepline algorithm.

  ! Method adapted from:
  ! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

  use consts_coms, only: r8

  implicit none

  integer,  intent(in ) :: np,nq ! Number of vertices in p and q
  real(r8), intent(in ) :: xp(np),yp(np) ! x,y coordinates of p vertices
  real(r8), intent(in ) :: xq(nq),yq(nq) ! x,y coordinates of q vertices
  real,     intent(in ) :: arp, arq
  real,     intent(out) :: area        ! area of overlap of p and q

  real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

  integer  :: nev         ! # of events
  integer  :: nsect       ! # of intersections between strip centerline and p,q edges
  real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
  real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges
  real(r8) :: xbot(np+nq) ! x-coords where strip bottom line intersects p and q edges
  real(r8) :: xtop(np+nq) ! x-coords where strip top line intersects p and q edges
  real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values
  real(r8) :: alphap(np)  ! if any perimeter points of p are on/inside q
  real(r8) :: alphaq(nq)  ! if any perimeter points of q are on/inside p
  real(r8) :: alpha

  integer  :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,ia,ib,is
  real(r8) :: p0,q0,dx,dy,dxtrap,dyi

  nev = 0

  ! Find vertices in p that are not outside q

  do ip = 1,np
     call inout_check(nq,xq,yq,xp(ip),yp(ip),alphap(ip))

     if (abs(alphap(ip)) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yp(ip)
     endif
  enddo

  if ( all( alphap(1:np) > 3._r8 ) ) then
     ! Polygon p vertices are entirely within or on the border of q.
     ! For our regulary-shaped polygons with non-concave faces we can exit here
     ! and assume p is entirely within q.
     area = arp
     return
  endif

  ! Find vertices in q that are not outside p

  do iq = 1,nq
     call inout_check(np,xp,yp,xq(iq),yq(iq),alphaq(iq))

     if (abs(alphaq(iq)) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yq(iq)
     endif
  enddo

  if ( all( alphaq(1:nq) > 3._r8 ) ) then
     ! Polygon q vertices are entirely within or on the border of p.
     ! For our regulary-shaped polygons with non-concave faces we can exit here
     ! and assume q is entirely within p.
     area = arq
     return
  endif

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

  area = 0.0

  ! Loop over event points

  do iev = 1,nev-1

     ! dy = width of strip between current and next event

     dy = yev(iev+1) - yev(iev)

     ! Reject overlap event if dy is less than 0.1 meter (threshold ok for single precision)

     if (dy < 0.1_r8) cycle

     ! ymid = y-coordinate of strip centerline.
     ! Initialize dx (centerline length sum) to 0.

     ymid = yev(iev) + 0.5_r8 * dy
     dx   = 0.0_r8

     ! Find x-coordinate of intersections of strip centerline with edges of p and q

     nsect = 0

     do ia = 1,np
        ib = ia + 1
        if (ia == np) ib = 1

        if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

        nsect = nsect + 1
        dyi = 1._r8 / (yp(ib) - yp(ia))

        xmid(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (ymid      -yp(ia)) * dyi
        xbot(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev)  -yp(ia)) * dyi
        xtop(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev+1)-yp(ia)) * dyi
     enddo

     do ia = 1,nq
        ib = ia + 1
        if (ia == nq) ib = 1

        if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

        nsect = nsect + 1
        dyi = 1._r8 / (yq(ib) - yq(ia))

        xmid(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (ymid      -yq(ia)) * dyi
        xbot(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev)  -yq(ia)) * dyi
        xtop(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev+1)-yq(ia)) * dyi
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

end subroutine polygon_overlap

!============================================================================

  subroutine polygon_overlap2(np,nq,xp,yp,xq,yq,arp,arq,area,alphap,alphaq)

  ! Given x,y coordinates of the vertices of polygons p and q, compute the area
  ! of overlap between the polygons using a sweepline algorithm.

  ! Method adapted from:
  ! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

  use consts_coms, only: r8

  implicit none

  integer,  intent(in) :: np,nq ! Number of vertices in p and q
  real(r8), intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
  real(r8), intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices
  real,     intent(in) :: arp, arq

  real,     intent(out) :: area        ! area of overlap of p and q
  real(r8), intent(out) :: alphap(np)  ! if any perimeter points of p are on/inside q
  real(r8), intent(out) :: alphaq(nq)  ! if any perimeter points of q are on/inside p

  real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

  integer :: nev      ! # of events
  integer :: nsect    ! # of intersections between strip centerline and p,q edges
  real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
  real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges
  real(r8) :: xbot(np+nq) ! x-coords where strip bottom line intersects p and q edges
  real(r8) :: xtop(np+nq) ! x-coords where strip top line intersects p and q edges
  real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values
  real(r8) :: alpha

  integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,ia,ib,is
  real(r8) :: p0,q0,dx,dy,dxtrap,dyi

  nev = 0

  ! Find vertices in p that are not outside q

  do ip = 1,np
     call inout_check(nq,xq,yq,xp(ip),yp(ip),alphap(ip))

     if (abs(alphap(ip)) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yp(ip)
     endif
  enddo

  if ( all( alphap(1:np) > 6._r8 ) ) then
     ! Polygon p vertices are entirely within (and not on the border) of q.
     ! For our regulary-shaped polygons with non-concave faces we can exit here
     ! and assume p is entirely within q.
     ! If any p vertices are on the border of q, we still need to check for
     ! q overlaps.
     alphaq(1:nq) = 0.0
     area         = arp
     return
  endif

  ! Find vertices in q that are not outside p

  do iq = 1,nq
     call inout_check(np,xp,yp,xq(iq),yq(iq),alphaq(iq))

     if (abs(alphaq(iq)) > 0.2_r8) then
        nev = nev + 1
        yev(nev) = yq(iq)
     endif
  enddo

  if ( all( alphaq(1:nq) > 3._r8 ) ) then
     ! Polygon q vertices are entirely within or on the border of p.
     ! For our regulary-shaped polygons with non-concave faces we can exit here
     ! and assume q is entirely within p.
     area = arq
     return
  endif

  if ( all( alphap(1:np) > 3._r8 ) ) then
     area = arp
     return
  endif

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

  area = 0.0

  ! Loop over event points

  do iev = 1,nev-1

     ! dy = width of strip between current and next event

     dy = yev(iev+1) - yev(iev)

     ! Reject overlap event if dy is less than 0.1 meter (threshold ok for single precision)

     if (dy < 0.1_r8) cycle

     ! ymid = y-coordinate of strip centerline.
     ! Initialize dx (centerline length sum) to 0.

     ymid = yev(iev) + 0.5_r8 * dy
     dx   = 0.0_r8

     ! Find x-coordinate of intersections of strip centerline with edges of p and q

     nsect = 0

     do ia = 1,np
        ib = ia + 1
        if (ia == np) ib = 1

        if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

        nsect = nsect + 1
        dyi = 1._r8 / (yp(ib) - yp(ia))

        xmid(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (ymid      -yp(ia)) * dyi
        xbot(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev)  -yp(ia)) * dyi
        xtop(nsect) = xp(ia) &
                    + (xp(ib) - xp(ia)) * (yev(iev+1)-yp(ia)) * dyi
     enddo

     do ia = 1,nq
        ib = ia + 1
        if (ia == nq) ib = 1

        if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

        nsect = nsect + 1
        dyi = 1._r8 / (yq(ib) - yq(ia))

        xmid(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (ymid      -yq(ia)) * dyi
        xbot(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev)  -yq(ia)) * dyi
        xtop(nsect) = xq(ia) &
                    + (xq(ib) - xq(ia)) * (yev(iev+1)-yq(ia)) * dyi
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

subroutine inout_check(n,x,y,x0,y0,th1)

  ! Given planar Cartesian coordinates of the vertices of a simple closed polygon,
  ! x(1),y(1),...,x(n),y(n), and of an additional point, x0,y0, determine
  ! whether the point is inside, outside, or on the border of the polygon.

  ! Set:
  ! th1 = 0 if outside
  ! th1 = pi if on the border
  ! th1 = 2 * pi if inside

  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: n  ! number of points in polygon

  real(r8), intent(in)  :: x(n), y(n) ! x,y coordinates of polygon points
  real(r8), intent(in)  :: x0,y0      ! x,y coordinates of additional point
  real(r8), intent(out) :: th1        ! output value

  ! Local variables

  integer  :: i, ip1
  real(r8) :: theta
  real(r8) :: x1(n), y1(n), xh(n), dd(n)

  real(r8), parameter :: pi1 = 3.1415926535898_r8, pi2 = 2.0_r8 * pi1

  do i = 1, n
     x1(i) = x(i) - x0
     y1(i) = y(i) - y0
     dd(i) = x1(i)*x1(i) + y1(i)*y1(i)
  enddo

  ! Check whether x0,y0 is exactly on a vertex (to within 1.0 m)

  if ( any( dd < 1.0_r8 ) ) then
     th1 = pi1
     return
  endif

#ifdef USE_MASS

  ! Use faster atan2 on IBM Power
  call vatan2(xh, y1, x1, n)
  do i = 1, n
     if (xh(i) < 0.0_r8) xh(i) = xh(i) + pi2
  enddo

#else

  do i = 1, n
     xh(i) = atan2(y1(i),x1(i))
     if (xh(i) < 0.0_r8) xh(i) = xh(i) + pi2
  enddo

#endif

  th1 = 0.0_r8

  do i = 1, n

     if (i == n) then
        ip1 = 1
     else
        ip1 = i + 1
     endif

     theta = xh(ip1) - xh(i)
     if (theta < -pi1) theta = theta + pi2

     if ((abs(abs(theta) - pi1) < 0.000000001_r8)) then
        th1 = pi1
        return
     endif

     if (theta > pi1) theta = theta - pi2
     th1 = th1 + theta
  enddo

  th1 = abs(th1)

end subroutine inout_check

!============================================================================

subroutine intersect(xpa,ypa,xpb,ypb,xqa,yqa,xqb,yqb,p0,q0,iflag)

  ! Given x,y coordinates of points pa and pb on line p and points qa and qb on
  ! line q, find parameteric values p0 and q0 where p and q intersect.
  ! If no intersection, set iflag = -1.

  use consts_coms, only: r8

  implicit none

  real(r8), intent(in) :: xpa,ypa,xpb,ypb  ! x,y coordinates of pa and pb
  real(r8), intent(in) :: xqa,yqa,xqb,yqb  ! x,y coordinates of qa and qb

  integer,  intent(out) :: iflag
  real(r8), intent(out) :: p0,q0

  real(r8) :: pabxqab  ! cross product of vectors pa-pb and qa-qb
  real(r8) :: pabxqabi

  iflag = -1

  pabxqab = (xpb - xpa) * (yqb - yqa) - (ypb - ypa) * (xqb - xqa)

  if (pabxqab /= 0.0_r8) then  ! Infinite lines intersect
     iflag = 0
     if (pabxqab < 0.0_r8) iflag = 1

     pabxqabi = 1._r8 / pabxqab

     p0 = ((yqb - yqa) * (xqa - xpa) + (yqa - ypa) * (xqa - xqb)) * pabxqabi
     q0 = ((ypb - ypa) * (xqa - xpa) + (yqa - ypa) * (xpa - xpb)) * pabxqabi
  endif

end subroutine intersect

!============================================================================

subroutine fltsort(n,f)

  ! Sort n floating point numbers f into ascending order

  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: n
  real(r8), intent(inout) :: f(n)

  integer :: i,j

  real(r8) :: f0

  do i = 1,n-1
     do j = i+1,n
        if (f(j) < f(i)) then
           f0 = f(i)
           f(i) = f(j)
           f(j) = f0
        endif
     enddo
  enddo

end subroutine fltsort

!============================================================================

subroutine fltsort2(n,f1,f2)

  ! Sort n floating point numbers in each of f1 and f2 into ascending order by f1

  use consts_coms, only: r8

  implicit none

  integer, intent(in) :: n
  real(r8), intent(inout) :: f1(n),f2(n)

  integer :: i,j

  real(r8) :: f0

  do i = 1,n-1
     do j = i+1,n
        if (f1(j) < f1(i)) then
           f0 = f1(i)
           f1(i) = f1(j)
           f1(j) = f0

           f0 = f2(i)
           f2(i) = f2(j)
           f2(j) = f0
        endif
     enddo
  enddo

end subroutine fltsort2

!============================================================================

subroutine fltsort3(n,f1,f2,f3)

  ! Sort n floating point numbers in each of f1, f2, and f3 into ascending
  ! order by f1

  use consts_coms, only: r8

  implicit none

  integer,  intent(in) :: n
  real(r8), intent(inout) :: f1(n),f2(n),f3(n)

  integer  :: i,j
  real(r8) :: f0

  do i = 1,n-1
     do j = i+1,n
        if (f1(j) < f1(i)) then
           f0 = f1(i)
           f1(i) = f1(j)
           f1(j) = f0

           f0 = f2(i)
           f2(i) = f2(j)
           f2(j) = f0

           f0 = f3(i)
           f3(i) = f3(j)
           f3(j) = f0
        endif
     enddo
  enddo

end subroutine fltsort3
