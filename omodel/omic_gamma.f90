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
real function gammp(a,x)

implicit none

real,intent(in) :: a,x

real :: gln,gammcf

if (x < a + 1.) then
   call gser(gammp,a,x,gln)
else
   call gcf(gammcf,a,x,gln)
   gammp = 1. - gammcf
endif

return
end function gammp

!===============================================================================

real function gammq(a,x)

implicit none

real, intent(in) :: a,x

real :: gamser,gln

if (x < a + 1.) then
   call gser(gamser,a,x,gln)
   gammq = 1. - gamser
else
   call gcf(gammq,a,x,gln)
endif

return
end function gammq

!===============================================================================

subroutine gcf(gammcf,a,x,gln)

implicit none

real, intent(out) :: gammcf,gln
real, intent(in) :: a,x

integer, parameter :: itmax=100
real, parameter :: eps=3.e-7
real,external :: gammln

real :: gold,a0,a1,b0,b1,fac,an,ana,anf,gaccel
integer :: n

gln = gammln(a)
gold = 0.
a0 = 1.
a1 = x
b0 = 0.
b1 = 1.
fac = 1.

do n = 1, itmax
    an = float(n)
    ana = an - a
    a0 = (a1 + a0 * ana) * fac
    b0 = (b1 + b0 * ana) * fac
    anf = an * fac
    a1 = x * a0 + anf * a1
    b1 = x * b0 + anf * b1
    if (a1 /= 0.) then
        fac = 1. / a1
        gaccel = b1 * fac
        if (abs((gaccel - gold) / gaccel) < eps) goto 20
        gold = gaccel
    endif
enddo

20 continue

gammcf = exp(-x + a * alog(x) - gln) * gaccel

if ((-x + a * log(x) - gln) > -38.) then
  gammcf = exp(-x + a * alog(x) - gln) * gaccel
else
  gammcf = 0.
endif

return
end subroutine gcf

!===============================================================================

subroutine gser(gamser,a,x,gln)

implicit none

real, intent(out) :: gamser,gln
real, intent(in) :: a,x

integer, parameter :: itmax=100
real, parameter :: eps=3.e-7
real,external ::gammln

real :: ap,sum,del
integer :: n

gln = gammln(a)

if (x <= 0.) then
    gamser = 0.
    return
endif

ap = a
sum = 1. / a
del = sum

do n = 1, itmax
    ap = ap + 1.
    del = del * x / ap
    sum =  sum  +  del
    if (abs(del) < abs(sum) * eps) goto 20
enddo

20 continue

if ((-x + a * log(x) - gln) > -38.) then
  gamser = sum * exp(-x + a * log(x) - gln)
else
  gamser = 0.
endif

return
end subroutine gser

!===============================================================================

real function gammln(xx)

implicit none

real, intent(in) :: xx

real(kind=8) :: cof(6),stp
data cof, stp/76.18009173d0, -86.50532033d0, 24.01409822d0,  &
     -1.231739516d0, .120858003d-2, -.536382d-5, 2.50662827465d0/
real(kind=8), parameter :: half=0.5d0, one=1.0d0, fpf=5.5d0

real :: x,tmp,ser
integer :: j

x = xx - one
tmp = x + fpf
tmp = (x + half) * log(tmp) - tmp
ser = one
do j = 1,6
    x = x + one
    ser = ser + cof(j) / x
enddo

gammln = tmp + log(stp * ser)

return
end function gammln
