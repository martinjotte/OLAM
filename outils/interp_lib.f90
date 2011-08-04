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
subroutine gdtost(a,ix,iy,stax,stay,staval)
implicit none
integer :: ix,iy
real :: a(ix,iy),r(4),scr(4),stax,stay,staval

!     SUBROUTINE TO RETURN STATIONS BACK-INTERPOLATED VALUES(STAVAL)
!     FROM UNIFORM GRID POINTS USING OVERLAPPING-QUADRATICS.
!     GRIDDED VALUES OF INPUT ARRAY A DIMENSIONED A(IX,IY),WHERE
!     IX=GRID POINTS IN X, IY = GRID POINTS IN Y .  STATION
!     LOCATION GIVEN IN TERMS OF GRID RELATIVE STATION X (STAX)
!     AND STATION COLUMN.
!     VALUES GREATER THAN 1.0E30 INDICATE MISSING DATA.

integer :: iy1,iy2,ix1,ix2,ii,i,jj,j
real :: fiym2,fixm2,yy,xx

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
101   jj=0
do 111 j=iy1,iy2
jj=jj+1
if(j.ge.1.and.j.le.iy) go to 112
r(jj)=1e30
go to 111
112   r(jj)=a(i,j)
111   continue
yy=stay-fiym2
call binom(1.,2.,3.,4.,r(1),r(2),r(3),r(4),yy,scr(ii))
100   continue
xx=stax-fixm2
call binom(1.,2.,3.,4.,scr(1),scr(2),scr(3),scr(4),xx,staval)
return
end

!===============================================================================

subroutine binom(x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy)

implicit none

real :: x1,x2,x3,x4,y1,y2,y3,y4,xxx,yyy
real :: wt1,wt2,yz22,yz23,yz24,yz11,yz12,yz13,yoo
integer :: istend

 yyy=1e30
 if(x2.gt.1.e19.or.x3.gt.1.e19.or.  &
   y2.gt.1.e19.or.y3.gt.1.e19)return
wt1=(xxx-x3)/(x2-x3)
wt2=1.0-wt1
istend=0
if(y4.lt.1.e19.and.x4.lt.1.e19) go to 410
yz22=wt1
yz23=wt2
yz24=0.0
istend= 1
410   if(y1.lt.1.e19.and.x1.lt.1.e19) go to 430
yz11=0.0
yz12=wt1
yz13=wt2
if(istend.eq.1)go to 480
go to 450
430   yz11=(xxx-x2)*(xxx-x3)/((x1-x2)*(x1-x3))
yz12=(xxx-x1)*(xxx-x3)/((x2-x1)*(x2-x3))
yz13=(xxx-x1)*(xxx-x2)/((x3-x1)*(x3-x2))
if(istend.eq.  1    ) go to 470
450   yz22=(xxx-x3)*(xxx-x4)/((x2-x3)*(x2-x4))
yz23=(xxx-x2)*(xxx-x4)/((x3-x2)*(x3-x4))
yz24=(xxx-x2)*(xxx-x3)/((x4-x2)*(x4-x3))
470   yyy=wt1*(yz11*y1+yz12*y2+yz13*y3)+wt2*(yz22*y2+yz23*y3+yz24*y4)
 go to 490
480      yyy=wt1*y2+wt2*y3
490   yoo=yyy
return
end

!===============================================================================

subroutine htint(nzz1,vctra,eleva,nzz2,vctrb,elevb)
implicit none
integer :: nzz1,nzz2
real :: vctra(nzz1),vctrb(nzz2),eleva(nzz1),elevb(nzz2)

integer :: l,k,kk
real :: wt

l=1
do 20 k=1,nzz2
30 continue
if(elevb(k).lt.eleva(1))go to 35
if(elevb(k).ge.eleva(l).and.elevb(k).le.eleva(l+1))go to 35
if(elevb(k).gt.eleva(nzz1))go to 36
l=l+1
if(l.eq.nzz1) then
  print *,'htint:nzz1',nzz1
  do kk=1,l
    print*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
  enddo
  stop 'htint'
endif
go to 30
35 continue
wt=(elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
vctrb(k)=vctra(l)+(vctra(l+1)-vctra(l))*wt
go to 20
36 continue
wt=(elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
vctrb(k)=vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
20 continue

return
end

!===============================================================================

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

real :: grada

ka = 1
do kb = 1,nb
   do while(ka < na-1 .and. elevb(kb) > eleva(ka+1))
      ka = ka + 1
   enddo
   grada = (vctra(ka+1) - vctra(ka)) / (eleva(ka+1) - eleva(ka))
   vctrb(kb) = vctra(ka) + grada * (elevb(kb) - eleva(ka))
enddo
return
end

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

real :: grada

ka = 1
do kb = 1,nb
   do while(ka < na-1 .and. elevb(kb) > eleva(ka+1))
      ka = ka + 1
   enddo
   grada = (vctra(ka+1) - vctra(ka)) / (eleva(ka+1) - eleva(ka))

   if (elevb(kb) < eleva(1)) then
      vctrb(kb) = vctra(1)
   elseif (elevb(kb) > eleva(na)) then
      vctrb(kb) = vctra(na)
   else
      vctrb(kb) = vctra(ka) + grada * (elevb(kb) - eleva(ka))
   endif
enddo
return
end

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

real :: grada

ka = 1
do kb = 1,nb
   do while(ka < na-1 .and. pressb(kb) < pressa(ka+1))
      ka = ka + 1
   enddo
   grada = (vctra(ka+1) - vctra(ka)) / (pressa(ka+1) - pressa(ka))
   vctrb(kb) = vctra(ka) + grada * (pressb(kb) - pressa(ka))
enddo
return
end

!===============================================================================

subroutine spline1(n,xdat,ydat,yppdat)

implicit none

integer, intent(in) :: n
real, intent(in) :: xdat(n),ydat(n)
real, intent(out) :: yppdat(n)

integer i,j
real fx,fy
real, allocatable, dimension(:) :: scr

allocate(scr(n))

yppdat(1) = 0.
yppdat(n) = 0.
scr(1) = 0.

do i = 2,n-1
   fx = (xdat(i)-xdat(i-1)) / (xdat(i+1)-xdat(i-1))
   fy = fx * yppdat(i-1) + 2.
   yppdat(i) = (fx-1.) / fy
   scr(i) = (6.*((ydat(i+1)-ydat(i)) / (xdat(i+1)-xdat(i))  &
      - (ydat(i)-ydat(i-1)) / (xdat(i)-xdat(i-1)))        &
      / (xdat(i+1)-xdat(i-1)) - fx * scr(i-1)) / fy
enddo
do j = n-1,2,-1
   yppdat(j) = yppdat(j) * yppdat(j+1) + scr(j)
enddo

deallocate(scr)

return
end

!===============================================================================

subroutine spline2(n,xdat,ydat,yppdat,x,y)

implicit none

integer, intent(in) :: n
real, intent(in) :: xdat(n),ydat(n),yppdat(n)
real, intent(in) :: x
real, intent(out) :: y

integer :: j,jhi,jlo
real  :: a,b,h

jlo = 1
jhi = n
1 if (jhi-jlo > 1) then
   j = (jhi+jlo) / 2
   if (xdat(j) > x) then
      jhi = j
   else
      jlo = j
   endif
   go to 1
endif

h = xdat(jhi) - xdat(jlo)
if (h == 0) stop 'bad xdat input in s2'
a = (xdat(jhi)-x) / h
b = (x-xdat(jlo)) / h
y = a * ydat(jlo) + b * ydat(jhi)  &
  + ((a**3-a) * yppdat(jlo) + (b**3-b) * yppdat(jhi)) * h**2 / 6.
return
end

