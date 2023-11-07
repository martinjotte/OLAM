subroutine hintrp_ee(na,vctra,eleva,nb,vctrb,elevb)

  implicit none

  ! Interpolates vctra values at eleva heights linearly by height vctrb
  ! values at elevb heights.  For any elevb heights outside the range of
  ! eleva heights, vctrb values are extrapolated from the first or last two
  ! vctra values.

  integer, intent(in) :: na
  integer, intent(in) :: nb

  real, intent(in) :: vctra(na)
  real, intent(in) :: eleva(na)
  real, intent(in) :: elevb(nb)

  real, intent(out) :: vctrb(nb)

  integer :: ka
  integer :: kb
  real    :: grada

  ka = 1
  do kb = 1,nb

   do while(ka < na-1 .and. elevb(kb) > eleva(ka+1))
      ka = ka + 1
   enddo

   grada = (vctra(ka+1) - vctra(ka)) / (eleva(ka+1) - eleva(ka))
   vctrb(kb) = vctra(ka) + grada * (elevb(kb) - eleva(ka))

enddo

end subroutine hintrp_ee

!===============================================================================

subroutine hintrp_cc(na,vctra,eleva,nb,vctrb,elevb)

  implicit none

  ! Interpolates vctra values at eleva heights linearly by height vctrb
  ! values at elevb heights.  For any elevb heights outside the range of
  ! eleva heights, vctrb values are held constant with height and set equal
  ! to the first or last vctra value.

  integer, intent(in) :: na
  integer, intent(in) :: nb

  real, intent(in) :: vctra(na)
  real, intent(in) :: eleva(na)
  real, intent(in) :: elevb(nb)

  real, intent(out) :: vctrb(nb)

  integer :: ka
  integer :: kb
  real    :: grada

  ka = 1
  do kb = 1,nb

     do while(ka < na-1 .and. elevb(kb) > eleva(ka+1))
        ka = ka + 1
     enddo

     if (elevb(kb) < eleva(1)) then
        vctrb(kb) = vctra(1)
     elseif (elevb(kb) > eleva(na)) then
        vctrb(kb) = vctra(na)
     else
        grada = (vctra(ka+1) - vctra(ka)) / (eleva(ka+1) - eleva(ka))
        vctrb(kb) = vctra(ka) + grada * (elevb(kb) - eleva(ka))
     endif

enddo

end subroutine hintrp_cc

!===============================================================================

subroutine pintrp_ee(na,vctra,pressa,nb,vctrb,pressb)

  implicit none

  integer, intent(in) :: na
  integer, intent(in) :: nb

  real, intent(in) :: vctra(na)
  real, intent(in) :: pressa(na)
  real, intent(in) :: pressb(nb)

  real, intent(out) :: vctrb(nb)

  integer :: ka
  integer :: kb
  real    :: grada

  ka = 1
  do kb = 1, nb

     do while(ka < na-1 .and. pressb(kb) < pressa(ka+1))
        ka = ka + 1
     enddo

     grada = (vctra(ka+1) - vctra(ka)) / (pressa(ka+1) - pressa(ka))
     vctrb(kb) = vctra(ka) + grada * (pressb(kb) - pressa(ka))

  enddo

end subroutine pintrp_ee

!===============================================================================

subroutine pintrp_cc(na,vctra,pressa,nb,vctrb,pressb)
  implicit none

  integer, intent(in)  :: na
  real,    intent(in)  :: vctra (na)
  real,    intent(in)  :: pressa(na)

  integer, intent(in)  :: nb
  real,    intent(out) :: vctrb (nb)
  real,    intent(in)  :: pressb(nb)

  integer :: ka
  integer :: kb
  real    :: grada

  ka = 1
  do kb = 1,nb

     do while(ka < na-1 .and. pressb(kb) < pressa(ka+1))
        ka = ka + 1
     enddo

     if (pressb(kb) >= pressa(1)) then
        vctrb(kb) = vctra(1)
     elseif (pressb(kb) <= pressa(na)) then
        vctrb(kb) = vctra(na)
     else
        grada = (vctra(ka+1) - vctra(ka)) / (pressa(ka+1) - pressa(ka))
        vctrb(kb) = vctra(ka) + grada * (pressb(kb) - pressa(ka))
     endif

  enddo

end subroutine pintrp_cc

!===============================================================================

subroutine spline1(n,xdat,ydat,yppdat)

  implicit none

  integer, intent(in)  :: n
  real,    intent(in)  :: xdat(n), ydat(n)
  real,    intent(out) :: yppdat(n)

  integer :: i,j
  real    :: fx,fy
  real    :: scr(n)

  yppdat(1) = 0.
  yppdat(n) = 0.
  scr   (1) = 0.

  do i = 2, n-1
     fx = (xdat(i)-xdat(i-1)) / (xdat(i+1)-xdat(i-1))
     fy = fx * yppdat(i-1) + 2.
     yppdat(i) = (fx-1.) / fy
     scr(i) = (6. * ( (ydat(i+1)-ydat(i)) / (xdat(i+1)-xdat(i))   &
                    - (ydat(i)-ydat(i-1)) / (xdat(i)-xdat(i-1)) ) &
                  / (xdat(i+1)-xdat(i-1)) - fx * scr(i-1)) / fy
  enddo

  do j = n-1,2,-1
     yppdat(j) = yppdat(j) * yppdat(j+1) + scr(j)
  enddo

end subroutine spline1

!===============================================================================

subroutine spline2(n,xdat,ydat,yppdat,x,y)

  implicit none

  integer, intent(in) :: n
  real,    intent(in) :: xdat(n), ydat(n), yppdat(n)
  real,    intent(in) :: x
  real,    intent(out) :: y

  integer :: j,jhi,jlo
  real    :: a,b,h

  jlo = 1
  jhi = n

  do while (jhi-jlo > 1)
     j = (jhi+jlo) / 2
     if (xdat(j) > x) then
        jhi = j
     else
        jlo = j
     endif
  enddo

  if (jhi-jlo /= 1) stop 'bad input in spline2'

  h = xdat(jhi) - xdat(jlo)
  a = (xdat(jhi)-x) / h
  b = (x-xdat(jlo)) / h

  y = a * ydat(jlo) + b * ydat(jhi)  &
    + ((a**3-a) * yppdat(jlo) + (b**3-b) * yppdat(jhi)) * h**2 / 6.

end subroutine spline2

!===============================================================================

subroutine spline2_vec(n,m,xdat,ydat,yppdat,x,y)

  implicit none

  integer, intent(in ) :: n, m
  real,    intent(in ) :: xdat(n)
  real,    intent(in ) :: ydat(m,n), yppdat(m,n)
  real,    intent(in ) :: x
  real,    intent(out) :: y(m)

  integer :: j,jhi,jlo
  real    :: a,b,h,c,d,e

  jlo = 1
  jhi = n

  do while (jhi-jlo > 1)
     j = (jhi+jlo) / 2
     if (xdat(j) > x) then
        jhi = j
     else
        jlo = j
     endif
  enddo

  if (jhi-jlo /= 1) stop 'bad input in spline2_vec'

  h = xdat(jhi) - xdat(jlo)
  e = h**2 / 6.
  a = (xdat(jhi)-x) / h
  b = (x-xdat(jlo)) / h
  c = (a**3 - a) * e
  d = (b**3 - b) * e

  do j = 1, m
     y(j) = a * ydat  (j,jlo) + b * ydat  (j,jhi) &
          + c * yppdat(j,jlo) + d * yppdat(j,jhi)
  enddo

end subroutine spline2_vec

!===============================================================================

subroutine gdtost_orig(a,ix,iy,stax,stay,staval)

  implicit none

  integer :: ix,iy
  real    :: a(ix,iy),r(4),scr(4),stax,stay,staval

  ! SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
  ! FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
  ! GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
  ! IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
  ! LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
  ! AND STATION COLUMN.
  ! VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.

  integer :: iy1,iy2,ix1,ix2,ii,i,jj,j
  real    :: fiym2,fixm2,yy,xx

  iy1=int(stay)-1
  iy2=iy1+3
  ix1=int(stax)-1
  ix2=ix1+3
  staval=1e30
  fiym2=float(iy1)-1
  fixm2=float(ix1)-1
  ii=0
  do 100 i=ix1,ix2
     ii=ii+1
     if(i.ge.1.and.i.le.ix) go to 101
     scr(ii)=1e30
     go to 100
101  jj=0
     do 111 j=iy1,iy2
        jj=jj+1
        if(j.ge.1.and.j.le.iy) go to 112
        r(jj)=1e30
        go to 111
112     r(jj)=a(i,j)
111  continue
     yy=stay-fiym2
     call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
100 continue
  xx=stax-fixm2
  call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)

end subroutine gdtost_orig

!===============================================================================

subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)

  implicit none

  real :: x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy
  real :: wt1,wt2,yz22,yz23,yz24,yz11,yz12,yz13,yoo
  integer :: istend

  yyy=1.e30
  if ( x2.gt.1.e19.or.x3.gt.1.e19 .or.  &
       y2.gt.1.e19.or.y3.gt.1.e19 ) return
  wt1=(xxx-x3)/(x2-x3)
  wt2=1.0-wt1
  istend=0
  if(y4.lt.1.e19.and.x4.lt.1.e19) go to 410
  yz22=wt1
  yz23=wt2
  yz24=0.0
  istend= 1
410 if(y1.lt.1.e19.and.x1.lt.1.e19) go to 430
  yz11=0.0
  yz12=wt1
  yz13=wt2
  if(istend.eq.1)go to 480
  go to 450
430 yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
  yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
  yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
  if(istend.eq.  1    ) go to 470
450 yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
  yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
  yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
470 yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
  go to 490
480 yyy=wt1*y2+wt2*y3
490 yoo=yyy

end subroutine binom

!===============================================================================

