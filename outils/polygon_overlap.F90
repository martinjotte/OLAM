!============================================================================

subroutine polygon_overlap(iwls,np,nq,xp,yp,xq,yq,area,ntrap,xtrap,ytrap,traparea)

! Given x,y coordinates of the vertices of polygons p and q, compute the area
! of overlap between the polygons using a sweepline algorithm.

! Method adapted from:
! Zerzan, J., 1989, Computers & Geosciences, Vol. 15, No. 7, pp. 1109-1114.

use consts_coms, only: r8

implicit none

integer, intent(in) :: np,nq,iwls ! Number of vertices in p and q
real(r8), intent(in) :: xp(np),yp(np) ! x,y coordinates of p vertices
real(r8), intent(in) :: xq(nq),yq(nq) ! x,y coordinates of q vertices

integer, intent(out) :: ntrap              ! number of trapezoids
real, intent(out) :: xtrap(4,np+nq+np*nq)  ! trapezoid x coordinates
real, intent(out) :: ytrap(4,np+nq+np*nq)  ! trapezoid y coordinates
real, intent(out) :: traparea(np+nq+np*nq) ! trapezoid area

real, intent(out) :: area                  ! area of overlap of p and q

real(r8) :: yev(np+nq+np*nq)  ! y-coordinates of event

integer :: nev      ! # of events
integer :: nsect    ! # of intersections between strip centerline and p,q edges
real(r8) :: ymid        ! y-coord of centerline of strip between consecutive events
real(r8) :: xmid(np+nq) ! x-coords where strip centerline intersects p and q edges 
real(r8) :: xbot(np+nq) ! x-coords where strip bottom line intersects p and q edges 
real(r8) :: xtop(np+nq) ! x-coords where strip top line intersects p and q edges 
real(r8) :: xcent       ! x-coord of midpoint of segment between xmid values

integer :: ip,iq,ipa,ipb,iqa,iqb,iflag,iev,i,ia,ib,is
real(r8) :: p0,q0,dx,dy,dxtrap

real(r8) :: alpha    ! angle swept by ray from point to polygon perimeter during circuit

! Find vertices in p that are not outside q

nev = 0
ntrap = 0
area = 0.

do ip = 1,np
   call inout_check(nq,xq,yq,xp(ip),yp(ip),alpha)

   if (abs(alpha) > 0.2_r8) then
      nev = nev + 1
      yev(nev) = yp(ip)
   endif
enddo

! Find vertices in q that are not outside p

do iq = 1,nq
   call inout_check(np,xp,yp,xq(iq),yq(iq),alpha)

   if (abs(alpha) > 0.2_r8) then
      nev = nev + 1
      yev(nev) = yq(iq)
   endif
enddo

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

if (nev == 0) then
!   print*, 'no overlap'
   return
endif

! Sort event list to obtain increasing order of yev

call fltsort(nev,yev)

! Loop over event points

do iev = 1,nev-1

! dy = width of strip between current and next event

   dy = yev(iev+1) - yev(iev)
   
! Reject overlap event if dy is less than 0.1 meter (threshold ok for single precision)
   
   if (dy < 0.1_r8) cycle

! ymid = y-coordinate of strip centerline.
! Initialize dx (centerline length sum) to 0.

   ymid = yev(iev) + 0.5_r8 * dy
   dx = 0.0_r8

! Find x-coordinate of intersections of strip centerline with edges of p and q

   nsect = 0

   do ia = 1,np
      ib = ia + 1
      if (ia == np) ib = 1

      if (ymid < min(yp(ia),yp(ib)) .or. ymid > max(yp(ia),yp(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xp(ia) &
                  + (xp(ib) - xp(ia)) * (ymid      -yp(ia)) / (yp(ib) - yp(ia))
      xbot(nsect) = xp(ia) &
                  + (xp(ib) - xp(ia)) * (yev(iev)  -yp(ia)) / (yp(ib) - yp(ia))
      xtop(nsect) = xp(ia) &
                  + (xp(ib) - xp(ia)) * (yev(iev+1)-yp(ia)) / (yp(ib) - yp(ia))

   enddo

   do ia = 1,nq
      ib = ia + 1
      if (ia == nq) ib = 1

      if (ymid < min(yq(ia),yq(ib)) .or. ymid > max(yq(ia),yq(ib))) cycle

      nsect = nsect + 1
      xmid(nsect) = xq(ia) &
                  + (xq(ib) - xq(ia)) * (ymid      -yq(ia)) / (yq(ib) - yq(ia))
      xbot(nsect) = xq(ia) &
                  + (xq(ib) - xq(ia)) * (yev(iev)  -yq(ia)) / (yq(ib) - yq(ia))
      xtop(nsect) = xq(ia) &
                  + (xq(ib) - xq(ia)) * (yev(iev+1)-yq(ia)) / (yq(ib) - yq(ia))

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

! Find x at 4 corners of current trapezoidal region

      ntrap = ntrap + 1
      
      xtrap(1,ntrap) = xbot(is)
      xtrap(2,ntrap) = xbot(is+1)
      xtrap(3,ntrap) = xtop(is+1)
      xtrap(4,ntrap) = xtop(is)
      
      ytrap(1,ntrap) = yev(iev) 
      ytrap(2,ntrap) = yev(iev) 
      ytrap(3,ntrap) = yev(iev+1) 
      ytrap(4,ntrap) = yev(iev+1) 

      traparea(ntrap) = dxtrap * dy
   enddo

   area = area + dx * dy

enddo

return
end subroutine polygon_overlap

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

real(r8), intent(in) :: x(n),y(n)  ! x,y coordinates of polygon points      
real(r8), intent(in) :: x0,y0      ! x,y coordinates of additional point
real(r8), intent(out) :: th1       ! output value

! Local variables

integer :: i
real(r8) :: theta
real(r8) :: xh0, xh1, xh2
real(r8) :: edgemax
real(r8) :: x1(n), y1(n)

real(r8), parameter :: pi1 = 3.1415926535898_r8, pi2 = 2.0_r8 * pi1

do i = 1,n
   x1(i) = x(i) - x0
   y1(i) = y(i) - y0

! Check whether x0,y0 is exactly on a vertex (to within 1.e-6 m)
!   if ((abs(x1(i)) < 0.000001_r8) .and. (abs(y1(i)) < 0.000001_r8)) then

! Check whether x0,y0 is exactly on a vertex (to within 1.0 m)

   if ((abs(x1(i)) < 1.0_r8) .and. (abs(y1(i)) < 1.0_r8)) then
      th1 = pi1
      return
   endif
enddo

th1 = 0.0_r8

xh2 = atan2(y1(1),x1(1))
if (xh2 < 0.0_r8) xh2 = xh2 + pi2

xh0 = xh2

do i = 1,n
   
   if (i == n) then
      xh1 = xh0
   else
      xh1 = atan2(y1(i+1),x1(i+1))
      if (xh1 < 0.0_r8) xh1 = xh1 + pi2
   endif

   theta = xh1 - xh2
   if (theta < -pi1) theta = theta + pi2

   if ((abs(abs(theta) - pi1) < 0.000000001_r8)) then
      th1 = pi1
      return
   endif

   if (theta > pi1) theta = theta - pi2
   th1 = th1 + theta
   
   xh2 = xh1
enddo

th1 = abs(th1)

return
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

integer, intent(out) :: iflag
real(r8), intent(out) :: p0,q0

real(r8) :: pabxqab  ! cross product of vectors pa-pb and qa-qb 

iflag = -1

pabxqab = (xpb - xpa) * (yqb - yqa) - (ypb - ypa) * (xqb - xqa)

if (pabxqab /= 0.0_r8) then  ! Infinite lines intersect
   iflag = 0
   if (pabxqab < 0.0_r8) iflag = 1

   p0 = ((yqb - yqa) * (xqa - xpa) + (yqa - ypa) * (xqa - xqb)) / pabxqab
   q0 = ((ypb - ypa) * (xqa - xpa) + (yqa - ypa) * (xpa - xpb)) / pabxqab
endif

!print*, 'pabxqab ',pabxqab,p0,q0,iflag

return
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

return
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

return
end subroutine fltsort2

!============================================================================

subroutine fltsort3(n,f1,f2,f3)

! Sort n floating point numbers in each of f1, f2, and f3 into ascending order by f1

use consts_coms, only: r8

implicit none

integer, intent(in) :: n
real(r8), intent(inout) :: f1(n),f2(n),f3(n)

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

         f0 = f3(i)
         f3(i) = f3(j)
         f3(j) = f0
      endif
   enddo
enddo

return
end subroutine fltsort3

!===============================================================================

subroutine matrix8_3x3(a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3,x1,x2,x3)

use consts_coms, only: r8

implicit none

real(r8), intent(in) :: a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3
real(r8), intent(out) :: x1,x2,x3

integer :: i
real(r8), dimension(4,3) :: abr

abr(1,1) = a11
abr(2,1) = a21
abr(3,1) = a31
abr(4,1) = b1

abr(1,2) = a12
abr(2,2) = a22
abr(3,2) = a32
abr(4,2) = b2

abr(1,3) = a13
abr(2,3) = a23
abr(3,3) = a33
abr(4,3) = b3

! Interchange rows if necessary so that first row has 
! largest (magnitude) element of first column

if (abs(abr(1,2)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange8(abr(i,1),abr(i,2))
   enddo
endif

if (abs(abr(1,3)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange8(abr(i,1),abr(i,3))
   enddo
endif

! Add -abr(1,2)/abr(1,1) times first row to second row and
! add -abr(1,3)/abr(1,1) times first row to third row.

do i = 2,4
   abr(i,2) = abr(i,2) - abr(1,2)/abr(1,1)*abr(i,1)
   abr(i,3) = abr(i,3) - abr(1,3)/abr(1,1)*abr(i,1)
enddo

! Interchange rows 2 and 3 if necessary so that second row
! has larger (magnitude) element of second column

if (abs(abr(2,3)) > abs(abr(2,2))) then
   do i = 2,4
      call rchange8(abr(i,2),abr(i,3))
   enddo
endif

! Add -abr(2,3)/abr(2,2) times second row to third row.

do i = 3,4
   abr(i,3) = abr(i,3) - abr(2,3)/abr(2,2)*abr(i,2)
enddo

! Back substitution

x3 = abr(4,3) / abr(3,3)
x2 = (abr(4,2) - abr(3,2) * x3) / abr(2,2)
x1 = (abr(4,1) - abr(2,1) * x2 - abr(3,1) * x3) / abr(1,1)

return
end subroutine matrix8_3x3

!===============================================================================

  subroutine rchange8(r1,r2)

  use consts_coms, only: r8

  implicit none

  real(r8), intent(inout) :: r1,r2
  real(r8) :: c

  c = r1
  r1 = r2
  r2 = c

  return
  end subroutine rchange8

