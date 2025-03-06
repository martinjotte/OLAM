module polygon_lib

  interface inout_check
     module procedure inout_check_r4
     module procedure inout_check_r8
  end interface inout_check

  interface polygon_overlap
     module procedure polygon_overlap_r4
     module procedure polygon_overlap_r8
  end interface polygon_overlap

  interface area_simple_polygon
     module procedure area_simple_polygon_r4
     module procedure area_simple_polygon_r8
  end interface area_simple_polygon

  private
  public :: polygon_overlap, inout_check, area_simple_polygon, &
            inout_check_nonconvex

contains

!============================================================================

subroutine polygon_overlap_r4(np,nq,xp,yp,xq,yq,arp,arq,area,alpha1,alpha2)

  ! Given x,y coordinates of the vertices of polygons p and q, compute the area
  ! of overlap between the polygons using a sweepline algorithm.

  ! Method adapted from:
  ! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  use               sortlib,         only: insertion_sort

  implicit none

  integer, intent(in)  :: np,nq         ! Number of vertices in p and q
  real,    intent(in)  :: xp(np),yp(np) ! x,y coordinates of p vertices
  real,    intent(in)  :: xq(nq),yq(nq) ! x,y coordinates of q vertices
  real,    intent(in)  :: arp, arq      ! areas of polygon p and q
  real,    intent(out) :: area          ! area of overlap of p and q

  integer :: nev               ! # of events
  real    :: yev(np+nq+np*nq)  ! y-coordinates of event

  integer :: nsect       ! # of intersections between strip centerline and p,q edges
  real    :: ymid        ! y-coord of centerline of strip between consecutive events
  real    :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges
  real    :: xcent       ! x-coord of midpoint of segment between xmid values

  integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,ia,ib,is
  real    :: p0,q0,dx,dy,dxtrap,pabxqab,pabxqabi,xev
  real    :: xqmax,xqmin,yqmax,yqmin
  real    :: xpmax,xpmin,ypmax,ypmin
  real    :: xppmax(np),xppmin(np),yppmax(np),yppmin(np)
  real    :: xqqmax(nq),xqqmin(nq),yqqmax(nq),yqqmin(nq)

  logical :: alphap(np)  ! if any perimeter points of p are on/inside q
  logical :: alphaq(nq)  ! if any perimeter points of q are on/inside p
  logical :: alpha

  real    :: dxp(np), dyp(np), dxp_dyp(np)
  real    :: dxq(nq), dyq(nq), dxq_dyq(nq)
  real    :: xpqmax, xpqmin
  real(r8):: ar8

  logical, optional, intent(in) :: alpha1(np), alpha2(nq)

  ! Note: if alpha1 is present it indicates if the polynomial p vertices are
  ! in/on q, and if alpha2 is present it indicates if the polynomial q vertices
  ! are in/on p.
  ! This is useful if we are checking the overlaps between an olam cell and an
  ! input grid, since we just need to pre-compute if a corner point is within an
  ! olam cell once rather than four times for each adjacent analysis box.

  area = 0.0

  ! First check if we may have overlaps

  xqmax = maxval(xq)
  xpmin = minval(xp)
  if (xqmax < xpmin) return

  xqmin = minval(xq)
  xpmax = maxval(xp)
  if (xqmin > xpmax) return

  yqmax = maxval(yq)
  ypmin = minval(yp)
  if (yqmax < ypmin) return

  yqmin = minval(yq)
  ypmax = maxval(yp)
  if (yqmin > ypmax) return

  ! Find vertices in p that are not outside q

  if (.not. present(alpha1)) then

     do iqa = 1, nq
        iqb = iqa + 1
        if (iqa == nq) iqb = 1
        dxq(iqa) = xq(iqb) - xq(iqa)
        dyq(iqa) = yq(iqb) - yq(iqa)
     enddo

     alphap = .false.

     do ip = 1, np

        if ( xp(ip) <= xqmax .and. xp(ip) >= xqmin .and. &
             yp(ip) <= yqmax .and. yp(ip) >= yqmin ) then

           call inout_check_r4(nq,xq,yq,xp(ip),yp(ip),alphap(ip),dx=dxq,dy=dyq)

        endif
     enddo

  else

     alphap = alpha1

  endif

  ! If Polygon p vertices are entirely within q, then for simple convex
  ! polygons we can exit here and assume all of p is entirely within q.

  if ( all( alphap ) ) then
     area = arp
     return
  endif

  ! Find vertices in q that are not outside p

  if (.not. present(alpha2)) then

     do ipa = 1, np
        ipb = ipa + 1
        if (ipa == np) ipb = 1
        dxp(ipa) = xp(ipb) - xp(ipa)
        dyp(ipa) = yp(ipb) - yp(ipa)
     enddo

     alphaq = .false.

     do iq = 1, nq

        if ( xq(iq) <= xpmax .and. xq(iq) >= xpmin .and. &
             yq(iq) <= ypmax .and. yq(iq) >= ypmin ) then

           call inout_check_r4(np,xp,yp,xq(iq),yq(iq),alphaq(iq),dx=dxp,dy=dyp)

        endif
     enddo

  else

     alphaq = alpha2

  endif

  ! If Polygon q vertices are entirely within p, then for simple convex
  ! polygons we can exit here and assume all of q is entirely within p.

  if ( all( alphaq ) ) then
     area = arq
     return
  endif

  ! If we are here either p or q is not entirely in the other polygon.
  ! Start searching for overlap events.

  nev = 0

  do ip = 1, np
     if ( alphap(ip) ) then
        nev = nev + 1
        yev(nev) = yp(ip)
     endif
  enddo

  do iq = 1, nq
     if ( alphaq(iq) ) then
        nev = nev + 1
        yev(nev) = yq(iq)
     endif
  enddo

  do ipa = 1,np
     ipb = ipa + 1
     if (ipa == np) ipb = 1

     if (present(alpha2)) then
        dxp(ipa) = xp(ipb) - xp(ipa)
        dyp(ipa) = yp(ipb) - yp(ipa)
     endif

     xppmin(ipa) = min(xp(ipa),xp(ipb))
     xppmax(ipa) = max(xp(ipa),xp(ipb))

     yppmin(ipa) = min(yp(ipa),yp(ipb))
     yppmax(ipa) = max(yp(ipa),yp(ipb))
  enddo

  do iqa = 1, nq
     iqb = iqa + 1
     if (iqa == nq) iqb = 1

     if (present(alpha1)) then
        dxq(iqa) = xq(iqb) - xq(iqa)
        dyq(iqa) = yq(iqb) - yq(iqa)
     endif

     xqqmin(iqa) = min(xq(iqa),xq(iqb))
     xqqmax(iqa) = max(xq(iqa),xq(iqb))

     yqqmin(iqa) = min(yq(iqa),yq(iqb))
     yqqmax(iqa) = max(yq(iqa),yq(iqb))
  enddo

  ! Find intersecting edges of polygons p and q

  do ipa = 1, np

     if ( xppmin(ipa) > xqmax ) cycle
     if ( xppmax(ipa) < xqmin ) cycle

     if ( yppmin(ipa) > yqmax ) cycle
     if ( yppmax(ipa) < yqmin ) cycle

     do iqa = 1, nq

        if ( xppmin(ipa) > xqqmax(iqa) ) cycle
        if ( xppmax(ipa) < xqqmin(iqa) ) cycle

        if ( yppmin(ipa) > yqqmax(iqa) ) cycle
        if ( yppmax(ipa) < yqqmin(iqa) ) cycle

        pabxqab = dxp(ipa) * dyq(iqa) - dyp(ipa) * dxq(iqa)

        ! Lines nearly parallel
        if (abs(pabxqab) < 1.e-15) cycle

        pabxqabi = 1.0 / pabxqab

        p0 = ( dyq(iqa) * (xq(iqa) - xp(ipa)) &
             - dxq(iqa) * (yq(iqa) - yp(ipa)) ) * pabxqabi

        if (p0 < -0.0000001 .or. p0 > 1.0000001) cycle

        q0 = ( dyp(ipa) * (xq(iqa) - xp(ipa)) &
             - dxp(ipa) * (yq(iqa) - yp(ipa)) ) * pabxqabi

        if (q0 < -0.0000001 .or. q0 > 1.0000001) cycle

        ! Line segments pa-pb and qa-qb intersect;
        ! find y-coordinate of intersection

        nev = nev + 1
        yev(nev) = yp(ipa) + p0 * dyp(ipa)

     enddo
  enddo

  ! Exit here if there are no more than one overlap events

  if (nev <= 1) return

  ! Sort event list to obtain increasing order of yev

  call insertion_sort(yev(1:nev))

  ar8 = 0.0_r8

  do ip = 1, np
     dxp_dyp(ip) = dxp(ip) / sign( max( abs(dyp(ip)), 1.e-15), dyp(ip) )
  enddo

  do iq = 1, nq
     dxq_dyq(iq) = dxq(iq) / sign( max( abs(dyq(iq)), 1.e-15), dyq(iq) )
  enddo

  xpqmax = min(xpmax,xqmax)
  xpqmin = max(xpmin,xqmin)

  ! Loop over event points

  do iev = 1, nev-1

     ! dy = width of strip between current and next event

     dy = yev(iev+1) - yev(iev)

     ! Reject overlap event if dy is less than 0.1 meter
     ! (threshold ok for single precision)

     if (dy < 0.1) cycle

     ! ymid = y-coordinate of strip centerline.
     ! Initialize dx (centerline length sum) to 0.

     ymid = yev(iev) + 0.5 * dy

     ! Find x-coordinate of intersections of strip centerline with edges of p and q

     nsect = 0

     do ia = 1, np
        if (ymid < yppmin(ia) .or. ymid > yppmax(ia)) cycle

        xev = xp(ia) + dxp_dyp(ia) * (ymid - yp(ia))

        if (xev > xpqmax .or. xev < xpqmin) cycle

        nsect = nsect + 1
        xmid(nsect) = xev
     enddo

     do ia = 1, nq
        if (ymid < yqqmin(ia) .or. ymid > yqqmax(ia)) cycle

        xev = xq(ia) + dxq_dyq(ia) * (ymid - yq(ia))

        if (xev > xpqmax .or. xev < xpqmin) cycle

        nsect = nsect + 1
        xmid(nsect) = xev
     enddo

     ! Skip this event if there are no more than one segment

     if (nsect <= 1) cycle

     ! Sort xmid values into increasing order

     call insertion_sort(xmid(1:nsect))

     ! check if the segment is inside both polygons

     dx = 0.0

     do is = 1, nsect - 1

        dxtrap = xmid(is+1) - xmid(is)

        ! Reject overlap event if dx segment is less than 0.1 meter
        ! (threshold ok for single precision)

        if (dxtrap < 0.1) cycle

        xcent = xmid(is) + 0.5 * dxtrap

        call inout_check_r4(np,xp,yp,xcent,ymid,alpha,dx=dxp,dy=dyp)
        if (.not. alpha) cycle

        call inout_check_r4(nq,xq,yq,xcent,ymid,alpha,dx=dxq,dy=dyq)
        if (.not. alpha) cycle

        dx = dx + dxtrap
     enddo

     ar8 = ar8 + dx * dy

  enddo

  area = ar8

end subroutine polygon_overlap_r4

!===============================================================================

subroutine polygon_overlap_r8(np,nq,xp,yp,xq,yq,arp,arq,area,alpha1,alpha2)

  ! Given x,y coordinates of the vertices of polygons p and q, compute the area
  ! of overlap between the polygons using a sweepline algorithm.

  ! Method adapted from:
  ! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  use               sortlib,         only: insertion_sort

  implicit none

  integer,  intent(in)  :: np,nq         ! Number of vertices in p and q
  real(r8), intent(in)  :: xp(np),yp(np) ! x,y coordinates of p vertices
  real(r8), intent(in)  :: xq(nq),yq(nq) ! x,y coordinates of q vertices
  real,     intent(in)  :: arp, arq      ! areas of polygon p and q
  real,     intent(out) :: area          ! area of overlap of p and q

  integer  :: nev               ! # of events
  real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

  integer  :: nsect       ! # of intersections between strip centerline and p,q edges
  real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
  real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges
  real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values

  integer  :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,ia,ib,is
  real(r8) :: p0,q0,dx,dy,dxtrap,pabxqab,pabxqabi,xev
  real(r8) :: xqmax,xqmin,yqmax,yqmin
  real(r8) :: xpmax,xpmin,ypmax,ypmin
  real(r8) :: xppmax(np),xppmin(np),yppmax(np),yppmin(np)
  real(r8) :: xqqmax(nq),xqqmin(nq),yqqmax(nq),yqqmin(nq)

  logical  :: alphap(np)  ! if any perimeter points of p are on/inside q
  logical  :: alphaq(nq)  ! if any perimeter points of q are on/inside p
  logical  :: alpha

  real(r8) :: dxp(np), dyp(np), dxp_dyp(np)
  real(r8) :: dxq(nq), dyq(nq), dxq_dyq(nq)
  real(r8) :: xpqmax, xpqmin
  real(r8) :: ar8

  logical, optional, intent(in) :: alpha1(np), alpha2(nq)

  ! Note: if alpha1 is present it indicates if the polynomial p vertices are
  ! in/on q, and if alpha2 is present it indicates if the polynomial q vertices
  ! are in/on p.
  ! This is useful if we are checking the overlaps between an olam cell and an
  ! input grid, since we just need to pre-compute if a corner point is within an
  ! olam cell once rather than four times for each adjacent analysis box.

  area = 0.0

  ! First check if we may have overlaps

  xqmax = maxval(xq)
  xpmin = minval(xp)
  if (xqmax < xpmin) return

  xqmin = minval(xq)
  xpmax = maxval(xp)
  if (xqmin > xpmax) return

  yqmax = maxval(yq)
  ypmin = minval(yp)
  if (yqmax < ypmin) return

  yqmin = minval(yq)
  ypmax = maxval(yp)
  if (yqmin > ypmax) return

  ! Find vertices in p that are not outside q

  if (.not. present(alpha1)) then

     do iqa = 1, nq
        iqb = iqa + 1
        if (iqa == nq) iqb = 1
        dxq(iqa) = xq(iqb) - xq(iqa)
        dyq(iqa) = yq(iqb) - yq(iqa)
     enddo

     alphap = .false.

     do ip = 1, np

        if ( xp(ip) <= xqmax .and. xp(ip) >= xqmin .and. &
             yp(ip) <= yqmax .and. yp(ip) >= yqmin ) then

           call inout_check_r8(nq,xq,yq,xp(ip),yp(ip),alphap(ip),dx=dxq,dy=dyq)

        endif
     enddo

  else

     alphap = alpha1

  endif

  ! If Polygon p vertices are entirely within q, then for simple convex
  ! polygons we can exit here and assume all of p is entirely within q.

  if ( all( alphap ) ) then
     area = arp
     return
  endif

  ! Find vertices in q that are not outside p

  if (.not. present(alpha2)) then

     do ipa = 1, np
        ipb = ipa + 1
        if (ipa == np) ipb = 1
        dxp(ipa) = xp(ipb) - xp(ipa)
        dyp(ipa) = yp(ipb) - yp(ipa)
     enddo

     alphaq = .false.

     do iq = 1, nq

        if ( xq(iq) <= xpmax .and. xq(iq) >= xpmin .and. &
             yq(iq) <= ypmax .and. yq(iq) >= ypmin ) then

           call inout_check_r8(np,xp,yp,xq(iq),yq(iq),alphaq(iq),dx=dxp,dy=dyp)

        endif
     enddo

  else

     alphaq = alpha2

  endif

  ! If Polygon q vertices are entirely within p, then for simple convex
  ! polygons we can exit here and assume all of q is entirely within p.

  if ( all( alphaq ) ) then
     area = arq
     return
  endif

  ! If we are here either p or q is not entirely in the other polygon.
  ! Start searching for overlap events.

  nev = 0

  do ip = 1, np
     if ( alphap(ip) ) then
        nev = nev + 1
        yev(nev) = yp(ip)
     endif
  enddo

  do iq = 1, nq
     if ( alphaq(iq) ) then
        nev = nev + 1
        yev(nev) = yq(iq)
     endif
  enddo

  do ipa = 1,np
     ipb = ipa + 1
     if (ipa == np) ipb = 1

     if (present(alpha2)) then
        dxp(ipa) = xp(ipb) - xp(ipa)
        dyp(ipa) = yp(ipb) - yp(ipa)
     endif

     xppmin(ipa) = min(xp(ipa),xp(ipb))
     xppmax(ipa) = max(xp(ipa),xp(ipb))

     yppmin(ipa) = min(yp(ipa),yp(ipb))
     yppmax(ipa) = max(yp(ipa),yp(ipb))
  enddo

  do iqa = 1, nq
     iqb = iqa + 1
     if (iqa == nq) iqb = 1

     if (present(alpha1)) then
        dxq(iqa) = xq(iqb) - xq(iqa)
        dyq(iqa) = yq(iqb) - yq(iqa)
     endif

     xqqmin(iqa) = min(xq(iqa),xq(iqb))
     xqqmax(iqa) = max(xq(iqa),xq(iqb))

     yqqmin(iqa) = min(yq(iqa),yq(iqb))
     yqqmax(iqa) = max(yq(iqa),yq(iqb))
  enddo

  ! Find intersecting edges of polygons p and q

  do ipa = 1, np

     if ( xppmin(ipa) > xqmax ) cycle
     if ( xppmax(ipa) < xqmin ) cycle

     if ( yppmin(ipa) > yqmax ) cycle
     if ( yppmax(ipa) < yqmin ) cycle

     do iqa = 1,nq

        if ( xppmin(ipa) > xqqmax(iqa) ) cycle
        if ( xppmax(ipa) < xqqmin(iqa) ) cycle

        if ( yppmin(ipa) > yqqmax(iqa) ) cycle
        if ( yppmax(ipa) < yqqmin(iqa) ) cycle

        ! Lines nearly parallel
        pabxqab = dxp(ipa) * dyq(iqa) - dyp(ipa) * dxq(iqa)

        if (abs(pabxqab) < 1.e-20_r8) cycle

        pabxqabi = 1._r8 / pabxqab

        p0 = ( dyq(iqa) * (xq(iqa) - xp(ipa)) &
             - dxq(iqa) * (yq(iqa) - yp(ipa)) ) * pabxqabi

        if (p0 < -0.000000001_r8 .or. p0 > 1.000000001_r8) cycle

        q0 = ( dyp(ipa) * (xq(iqa) - xp(ipa)) &
             - dxp(ipa) * (yq(iqa) - yp(ipa)) ) * pabxqabi

        if (q0 < -0.000000001_r8 .or. q0 > 1.000000001_r8) cycle

        ! Line segments pa-pb and qa-qb intersect;
        ! find y-coordinate of intersection

        nev = nev + 1
        yev(nev) = yp(ipa) + p0 * dyp(ipa)

     enddo
  enddo

  ! Exit here if there are no more than one overlap events

  if (nev <= 1) return

  ! Sort event list to obtain increasing order of yev

  call insertion_sort(yev(1:nev))

  ar8 = 0.0_r8

  do ip = 1, np
     dxp_dyp(ip) = dxp(ip) / sign( max( abs(dyp(ip)), 1.e-20_r8), dyp(ip) )
  enddo

  do iq = 1, nq
     dxq_dyq(iq) = dxq(iq) / sign( max( abs(dyq(iq)), 1.e-20_r8), dyq(iq) )
  enddo

  xpqmax = min(xpmax,xqmax)
  xpqmin = max(xpmin,xqmin)

  ! Loop over event points

  do iev = 1, nev-1

     ! dy = width of strip between current and next event

     dy = yev(iev+1) - yev(iev)

     ! Reject overlap event if dy is less than 0.1 meter
     ! (threshold ok for single precision)

     if (dy < 0.1_r8) cycle

     ! ymid = y-coordinate of strip centerline.
     ! Initialize dx (centerline length sum) to 0.

     ymid = yev(iev) + 0.5_r8 * dy

     ! Find x-coordinate of intersections of strip centerline with edges of p and q

     nsect = 0

     do ia = 1, np
        if (ymid < yppmin(ia) .or. ymid > yppmax(ia)) cycle

        xev = xp(ia) + dxp_dyp(ia) * (ymid - yp(ia))

        if (xev > xpqmax .or. xev < xpqmin) cycle

        nsect = nsect + 1
        xmid(nsect) = xev
     enddo

     do ia = 1, nq
        if (ymid < yqqmin(ia) .or. ymid > yqqmax(ia)) cycle

        xev = xq(ia) + dxq_dyq(ia) * (ymid - yq(ia))

        if (xev > xpqmax .or. xev < xpqmin) cycle

        nsect = nsect + 1
        xmid(nsect) = xev
     enddo

     ! Skip this event if there are no more than one segment

     if (nsect <= 1) cycle

     ! Sort xmid values into increasing order

     call insertion_sort(xmid(1:nsect))

     ! check if the segment is inside both polygons

     dx = 0.0_r8

     do is = 1, nsect - 1

        dxtrap = xmid(is+1) - xmid(is)

        ! Reject overlap event if dx segment is less than 0.1 meter
        ! (threshold ok for single precision)

        if (dxtrap < 0.1_r8) cycle

        xcent = xmid(is) + 0.5_r8 * dxtrap

        call inout_check_r8(np,xp,yp,xcent,ymid,alpha,dx=dxp,dy=dyp)
        if (.not. alpha) cycle

        call inout_check_r8(nq,xq,yq,xcent,ymid,alpha,dx=dxq,dy=dyq)
        if (.not. alpha) cycle

        dx = dx + dxtrap
     enddo

     ar8 = ar8 + dx * dy

  enddo

  area = ar8

end subroutine polygon_overlap_r8

!===============================================================================

subroutine inout_check_r4(n,x,y,x0,y0,icheck,dx,dy)

  ! Given planar Cartesian coordinates of the vertices of a simple closed
  ! CONVEX polygon, x(1),y(1),...,x(n),y(n), and of an additional point, x0,y0,
  ! determine if the point is inside, outside, or on the border of the polygon.

  implicit none

  integer,  intent(in)  :: n          ! number of polygon vertices
  real,     intent(in)  :: x(n), y(n) ! x,y coordinates of polygon vertices
  real,     intent(in)  :: x0, y0     ! x,y coordinate of point to check
  logical,  intent(out) :: icheck     ! result (.true. if point in or on polygon)

  real, optional, intent(in) :: dx(n), dy(n)

  real, parameter :: tol = 1.e-4

  real    :: xprod(n)
  integer :: ia, ib

  ! For a convex polygon (all interior angles less than 180 degrees), a point
  ! lies inside the polygon if it is on the same side of all the line segments
  ! making up the edges of the polygon

  if (present(dx) .and. present(dy)) then

     do ia = 1, n
        xprod(ia) = (y0 - y(ia)) * dx(ia)  &
                  - (x0 - x(ia)) * dy(ia)
     enddo

  else

     do ia = 1, n
        ib = ia + 1
        if (ia == n) ib = 1

        xprod(ia) = (y0 - y(ia)) * (x(ib) - x(ia)) &
                  - (x0 - x(ia)) * (y(ib) - y(ia))
     enddo

  endif

  ! Check if point is to the right of all line segments

  icheck = all( xprod > -tol )

  ! If not, check if point is to the left of all line segments

  if (.not. icheck) then
     icheck = all( xprod < tol )
  endif

end subroutine inout_check_r4

!===============================================================================

subroutine inout_check_r8(n,x,y,x0,y0,icheck,dx,dy)

  ! Given planar Cartesian coordinates of the vertices of a simple closed
  ! CONVEX polygon, x(1),y(1),...,x(n),y(n), and of an additional point, x0,y0,
  ! determine if the point is inside, outside, or on the border of the polygon.

  use, intrinsic :: iso_fortran_env, only: r8=>real64

  implicit none

  integer,  intent(in)  :: n          ! number of polygon vertices
  real(r8), intent(in)  :: x(n), y(n) ! x,y coordinates of polygon vertices
  real(r8), intent(in)  :: x0, y0     ! x,y coordinate of point to check
  logical , intent(out) :: icheck     ! result (.true. if point in or on polygon)

  real(r8), optional, intent(in) :: dx(n), dy(n)

  real(r8), parameter :: tol = 1.e-4_r8

  real(r8) :: xprod(n)
  integer  :: ia, ib

  ! For a convex polygon (all interior angles less than 180 degrees), a point
  ! lies inside the polygon if it is on the same side of all the line segments
  ! making up the edges of the polygon

  if (present(dx) .and. present(dy)) then

     do ia = 1, n
        xprod(ia) = (y0 - y(ia)) * dx(ia)  &
                  - (x0 - x(ia)) * dy(ia)
     enddo

  else

     do ia = 1, n
        ib = ia + 1
        if (ia == n) ib = 1

        xprod(ia) = (y0 - y(ia)) * (x(ib) - x(ia)) &
                  - (x0 - x(ia)) * (y(ib) - y(ia))
     enddo

  endif

  ! Check if point is to the right of all line segments

  icheck = all( xprod > -tol )

  ! If not, check if point is to the left of all line segments

  if (.not. icheck) then
     icheck = all( xprod < tol )
  endif

end subroutine inout_check_r8

!============================================================================

subroutine inout_check_nonconvex(n,x,y,x0,y0,th1)

  ! Given planar Cartesian coordinates of the vertices of a simple closed polygon,
  ! x(1),y(1),...,x(n),y(n), and of an additional point, x0,y0, determine
  ! whether the point is inside, outside, or on the border of the polygon.

  ! Set:
  ! th1 = 0 if outside
  ! th1 = pi if on the border
  ! th1 = 2 * pi if inside

  use, intrinsic :: iso_fortran_env, only: r8=>real64

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

  do i = 1, n
     xh(i) = atan2(y1(i),x1(i))
     if (xh(i) < 0.0_r8) xh(i) = xh(i) + pi2
  enddo

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

end subroutine inout_check_nonconvex

!============================================================================

subroutine intersect(xpa,ypa,xpb,ypb,xqa,yqa,xqb,yqb,p0,q0,iflag)

  ! Given x,y coordinates of points pa and pb on line p and points qa and qb on
  ! line q, find parameteric values p0 and q0 where p and q intersect.
  ! If no intersection, set iflag = -1.

  use, intrinsic :: iso_fortran_env, only: r8=>real64

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

subroutine area_simple_polygon_r4(n,x,y,a)

  ! Computes the area of any simple convex polygon from
  ! the coordinates of its vertices

  use, intrinsic :: iso_fortran_env, only: r8=>real64

  implicit none

  integer, intent(in)  :: n    ! number of vertices
  real,    intent(in)  :: x(n) ! x coordinate of vertices
  real,    intent(in)  :: y(n) ! y coordinate of vertices
  real,    intent(out) :: a    ! area

  real(r8)             :: area
  integer              :: i, ii

  area = real( x(n) * y(1), r8 ) &
       - real( y(n) * x(i), r8 )

  do i = 1, n-1
     ii = i + 1
     area = area + real( x(i) * y(ii), r8 ) &
                 - real( y(i) * x(ii), r8 )
  enddo

  a = abs( 0.5 * real( area ) )

end subroutine area_simple_polygon_r4

!============================================================================

subroutine area_simple_polygon_r8(n,x,y,a)

  ! Computes the area of any simple convex polygon from
  ! the coordinates of its vertices

  use, intrinsic :: iso_fortran_env, only: r8=>real64

  implicit none

  integer,  intent(in)  :: n    ! number of vertices
  real(r8), intent(in)  :: x(n) ! x coordinate of vertices
  real(r8), intent(in)  :: y(n) ! y coordinate of vertices
  real(r8), intent(out) :: a    ! area
  integer               :: i, ii
  real(r8)              :: area ! area

  area = x(n) * y(1) - y(n) * x(1)

  do i = 1, n-1
     ii = i + 1
     area = area + x(i) * y(ii) - y(i) * x(ii)
  enddo

  a = abs( 0.5_r8 * area )

end subroutine area_simple_polygon_r8

end module polygon_lib
