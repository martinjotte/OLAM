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
subroutine tri_neighbors()

use mem_ijtabs, only: itab_md, itab_ud, itab_wd
use mem_grid,   only: nua, nwa
use misc_coms,  only: io6

implicit none

integer :: iu,iw
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9
integer :: iw1,iw2,iw3,iw4,iw5,iw6
integer :: iw1_iu1,iw1_iu2,iw1_iu3
integer :: iw2_iu1,iw2_iu2,iw2_iu3
integer :: j,iunow,iu0,npoly,im

! Loop over W points

do iw = 2,nwa
   itab_wd(iw)%npoly = 3

   iu1 = itab_wd(iw)%iu(1)
   iu2 = itab_wd(iw)%iu(2)
   iu3 = itab_wd(iw)%iu(3)

! Fill M points for current W point

   if (iw == itab_ud(iu2)%iw(1)) then          ! New sequence???
      itab_wd(iw)%im(1)   = itab_ud(iu2)%im(1) ! [1 = 2]
      itab_wd(iw)%im(3)   = itab_ud(iu2)%im(2) ! [2 = 1]
   else
      itab_wd(iw)%im(1)   = itab_ud(iu2)%im(2) ! [1 = 1]
      itab_wd(iw)%im(3)   = itab_ud(iu2)%im(1) ! [2 = 2]
   endif

   if (iw == itab_ud(iu1)%iw(1)) then
      itab_wd(iw)%im(2)   = itab_ud(iu1)%im(2) ! [3 = 2]
   else
      itab_wd(iw)%im(2)   = itab_ud(iu1)%im(1) ! [3 = 1]
   endif

! Fill inner W points and U directions for current W point

   if (iw == itab_ud(iu1)%iw(1)) then
      itab_wd(iw)%iw(1)   = itab_ud(iu1)%iw(2)
      itab_wd(iw)%diru(1) = -1.
   else
      itab_wd(iw)%iw(1)   = itab_ud(iu1)%iw(1)
      itab_wd(iw)%diru(1) = 1.
   endif

   if (iw == itab_ud(iu2)%iw(1)) then
      itab_wd(iw)%iw(2)   = itab_ud(iu2)%iw(2)
      itab_wd(iw)%diru(2) = -1.
   else
      itab_wd(iw)%iw(2)   = itab_ud(iu2)%iw(1)
      itab_wd(iw)%diru(2) = 1.
   endif

   if (iw == itab_ud(iu3)%iw(1)) then
      itab_wd(iw)%iw(3)   = itab_ud(iu3)%iw(2)
      itab_wd(iw)%diru(3) = -1.
   else
      itab_wd(iw)%iw(3)   = itab_ud(iu3)%iw(1)
      itab_wd(iw)%diru(3) = 1.
   endif
enddo

! Fill outer U and W points for current W point

do iw = 2,nwa
   iw1 = itab_wd(iw)%iw(1)
   iw2 = itab_wd(iw)%iw(2)
   iw3 = itab_wd(iw)%iw(3)

! This should work ok for iw1 = 1, but may bypass if desired

   if (iw == itab_wd(iw1)%iw(1)) then

      itab_wd(iw)%iu(4) = itab_wd(iw1)%iu(2)
      itab_wd(iw)%iu(5) = itab_wd(iw1)%iu(3)

      itab_wd(iw)%iw(4) = itab_wd(iw1)%iw(2)
      itab_wd(iw)%iw(5) = itab_wd(iw1)%iw(3)

   elseif (iw == itab_wd(iw1)%iw(2)) then

      itab_wd(iw)%iu(4) = itab_wd(iw1)%iu(3)
      itab_wd(iw)%iu(5) = itab_wd(iw1)%iu(1)

      itab_wd(iw)%iw(4) = itab_wd(iw1)%iw(3)
      itab_wd(iw)%iw(5) = itab_wd(iw1)%iw(1)

   else

      itab_wd(iw)%iu(4) = itab_wd(iw1)%iu(1)
      itab_wd(iw)%iu(5) = itab_wd(iw1)%iu(2)

      itab_wd(iw)%iw(4) = itab_wd(iw1)%iw(1)
      itab_wd(iw)%iw(5) = itab_wd(iw1)%iw(2)

   endif

! This should work ok for iw2 = 1, but may bypass if desired

   if (iw == itab_wd(iw2)%iw(1)) then

      itab_wd(iw)%iu(6) = itab_wd(iw2)%iu(2)
      itab_wd(iw)%iu(7) = itab_wd(iw2)%iu(3)

      itab_wd(iw)%iw(6) = itab_wd(iw2)%iw(2)
      itab_wd(iw)%iw(7) = itab_wd(iw2)%iw(3)

   elseif (iw == itab_wd(iw2)%iw(2)) then

      itab_wd(iw)%iu(6) = itab_wd(iw2)%iu(3)
      itab_wd(iw)%iu(7) = itab_wd(iw2)%iu(1)

      itab_wd(iw)%iw(6) = itab_wd(iw2)%iw(3)
      itab_wd(iw)%iw(7) = itab_wd(iw2)%iw(1)

   else

      itab_wd(iw)%iu(6) = itab_wd(iw2)%iu(1)
      itab_wd(iw)%iu(7) = itab_wd(iw2)%iu(2)

      itab_wd(iw)%iw(6) = itab_wd(iw2)%iw(1)
      itab_wd(iw)%iw(7) = itab_wd(iw2)%iw(2)

   endif

! This should work ok for iw3 = 1, but may bypass if desired

   if (iw == itab_wd(iw3)%iw(1)) then

      itab_wd(iw)%iu(8) = itab_wd(iw3)%iu(2)
      itab_wd(iw)%iu(9) = itab_wd(iw3)%iu(3)

      itab_wd(iw)%iw(8) = itab_wd(iw3)%iw(2)
      itab_wd(iw)%iw(9) = itab_wd(iw3)%iw(3)

   elseif (iw == itab_wd(iw3)%iw(2)) then

      itab_wd(iw)%iu(8) = itab_wd(iw3)%iu(3)
      itab_wd(iw)%iu(9) = itab_wd(iw3)%iu(1)

      itab_wd(iw)%iw(8) = itab_wd(iw3)%iw(3)
      itab_wd(iw)%iw(9) = itab_wd(iw3)%iw(1)

   else

      itab_wd(iw)%iu(8) = itab_wd(iw3)%iu(1)
      itab_wd(iw)%iu(9) = itab_wd(iw3)%iu(2)

      itab_wd(iw)%iw(8) = itab_wd(iw3)%iw(1)
      itab_wd(iw)%iw(9) = itab_wd(iw3)%iw(2)

   endif

enddo

! Loop over U points

do iu = 2,nua

   iw1 = itab_ud(iu)%iw(1)
   iw2 = itab_ud(iu)%iw(2)

   itab_ud(iu)%mrlu = max(itab_wd(iw1)%mrlw,itab_wd(iw2)%mrlw)

! This should work ok for iw1 = 1, but may bypass if desired

   iw1_iu1 = itab_wd(iw1)%iu(1)
   iw1_iu2 = itab_wd(iw1)%iu(2)
   iw1_iu3 = itab_wd(iw1)%iu(3)

! Fill IU1 and IU2 for current U point

   if (iw1_iu1 == iu) then
      itab_ud(iu)%iu(1) = iw1_iu2
      itab_ud(iu)%iu(2) = iw1_iu3
   elseif (iw1_iu2 == iu) then
      itab_ud(iu)%iu(1) = iw1_iu3
      itab_ud(iu)%iu(2) = iw1_iu1
   else
      itab_ud(iu)%iu(1) = iw1_iu1
      itab_ud(iu)%iu(2) = iw1_iu2
   endif

! This should work ok for iw2 = 1, but may bypass if desired

   iw2_iu1 = itab_wd(iw2)%iu(1)
   iw2_iu2 = itab_wd(iw2)%iu(2)
   iw2_iu3 = itab_wd(iw2)%iu(3)

! Fill IU3 and IU4 for current U point

   if (iw2_iu1 == iu) then
      itab_ud(iu)%iu(3) = iw2_iu3
      itab_ud(iu)%iu(4) = iw2_iu2
   elseif (iw2_iu2 == iu) then
      itab_ud(iu)%iu(3) = iw2_iu1
      itab_ud(iu)%iu(4) = iw2_iu3
   else
      itab_ud(iu)%iu(3) = iw2_iu2
      itab_ud(iu)%iu(4) = iw2_iu1
   endif

   iu1 = itab_ud(iu)%iu(1)
   iu2 = itab_ud(iu)%iu(2)
   iu3 = itab_ud(iu)%iu(3)
   iu4 = itab_ud(iu)%iu(4)

! Fill IW3 and DIRU1 for current U point
! This should work ok for iu1 = 1, but may bypass if desired

   if (itab_ud(iu1)%iw(1) == iw1) then
      itab_ud(iu)%iw(3) = itab_ud(iu1)%iw(2)
      itab_ud(iu)%diru(1) = -1.
   else
      itab_ud(iu)%iw(3) = itab_ud(iu1)%iw(1)
      itab_ud(iu)%diru(1) = 1.
   endif

! Fill IW4 and DIRU2 for current U point
! This should work ok for iu2 = 1, but may bypass if desired

   if (itab_ud(iu2)%iw(1) == iw1) then
      itab_ud(iu)%iw(4) = itab_ud(iu2)%iw(2)
      itab_ud(iu)%diru(2) = -1.
   else
      itab_ud(iu)%iw(4) = itab_ud(iu2)%iw(1)
      itab_ud(iu)%diru(2) = 1.
   endif

! Fill IW5 and DIRU3 for current U point
! This should work ok for iu3 = 1, but may bypass if desired

   if (itab_ud(iu3)%iw(1) == iw2) then
      itab_ud(iu)%iw(5) = itab_ud(iu3)%iw(2)
      itab_ud(iu)%diru(3) = -1.
   else
      itab_ud(iu)%iw(5) = itab_ud(iu3)%iw(1)
      itab_ud(iu)%diru(3) = 1.
   endif

! Fill IW6 and DIRU4 for current U point
! This should work ok for iu4 = 1, but may bypass if desired

   if (itab_ud(iu4)%iw(1) == iw2) then
      itab_ud(iu)%iw(6) = itab_ud(iu4)%iw(2)
      itab_ud(iu)%diru(4) = -1.
   else
      itab_ud(iu)%iw(6) = itab_ud(iu4)%iw(1)
      itab_ud(iu)%diru(4) = 1.
   endif

   iw3 = itab_ud(iu)%iw(3)
   iw4 = itab_ud(iu)%iw(4)
   iw5 = itab_ud(iu)%iw(5)
   iw6 = itab_ud(iu)%iw(6)

! Fill IU5 and IU6 for current U point
! This should work ok for iw3 = 1, but may bypass if desired

   if (iu1 == itab_wd(iw3)%iu(1)) then
      itab_ud(iu)%iu(5) = itab_wd(iw3)%iu(2)
      itab_ud(iu)%iu(6) = itab_wd(iw3)%iu(3)
   elseif (iu1 == itab_wd(iw3)%iu(2)) then
      itab_ud(iu)%iu(5) = itab_wd(iw3)%iu(3)
      itab_ud(iu)%iu(6) = itab_wd(iw3)%iu(1)
   else
      itab_ud(iu)%iu(5) = itab_wd(iw3)%iu(1)
      itab_ud(iu)%iu(6) = itab_wd(iw3)%iu(2)
   endif

! Fill IU7 and IU8 for current U point
! This should work ok for iw4 = 1, but may bypass if desired

   if (iu2 == itab_wd(iw4)%iu(1)) then
      itab_ud(iu)%iu(7) = itab_wd(iw4)%iu(2)
      itab_ud(iu)%iu(8) = itab_wd(iw4)%iu(3)
   elseif (iu2 == itab_wd(iw4)%iu(2)) then
      itab_ud(iu)%iu(7) = itab_wd(iw4)%iu(3)
      itab_ud(iu)%iu(8) = itab_wd(iw4)%iu(1)
   else
      itab_ud(iu)%iu(7) = itab_wd(iw4)%iu(1)
      itab_ud(iu)%iu(8) = itab_wd(iw4)%iu(2)
   endif

! Fill IU9 and IU10 for current U point
! This should work ok for iw5 = 1, but may bypass if desired

   if (iu3 == itab_wd(iw5)%iu(1)) then
      itab_ud(iu)%iu(9)  = itab_wd(iw5)%iu(3)
      itab_ud(iu)%iu(10) = itab_wd(iw5)%iu(2)
   elseif (iu3 == itab_wd(iw5)%iu(2)) then
      itab_ud(iu)%iu(9)  = itab_wd(iw5)%iu(1)
      itab_ud(iu)%iu(10) = itab_wd(iw5)%iu(3)
   else
      itab_ud(iu)%iu(9)  = itab_wd(iw5)%iu(2)
      itab_ud(iu)%iu(10) = itab_wd(iw5)%iu(1)
   endif

! Fill IU11 and IU12 for current U point
! This should work ok for iw6 = 1, but may bypass if desired

   if (iu4 == itab_wd(iw6)%iu(1)) then
      itab_ud(iu)%iu(11) = itab_wd(iw6)%iu(3)
      itab_ud(iu)%iu(12) = itab_wd(iw6)%iu(2)
   elseif (iu4 == itab_wd(iw6)%iu(2)) then
      itab_ud(iu)%iu(11) = itab_wd(iw6)%iu(1)
      itab_ud(iu)%iu(12) = itab_wd(iw6)%iu(3)
   else
      itab_ud(iu)%iu(11) = itab_wd(iw6)%iu(2)
      itab_ud(iu)%iu(12) = itab_wd(iw6)%iu(1)
   endif

enddo  ! end loop over U points

! Fill U and W points for M points (do this as loop over U points)

do iu = 2,nua
   do j = 1,2

      if (j == 1) im = itab_ud(iu)%im(1)
      if (j == 2) im = itab_ud(iu)%im(2)

      iw1 = itab_ud(iu)%iw(1)
      iw2 = itab_ud(iu)%iw(2)

      if (itab_md(im)%npoly == 0   .or.  &
         (iw1 == 1 .and. j == 1) .or.  &  ! This and next line added for walls
         (iw2 == 1 .and. j == 2)) then

         iunow = iu
         iu0 = 0
         npoly = 0

         do while (iunow /= iu0)

            iu0 = iu
!!                  npoly = npoly + 1

            if (itab_ud(iunow)%im(1) == im) then

               if (itab_ud(iunow)%iw(2) > 1) then

                  npoly = npoly + 1

                  itab_md(im)%iu(npoly) = iunow  ! NEW

                  itab_md(im)%iw(npoly) = itab_ud(iunow)%iw(2)
                  iunow = itab_ud(iunow)%iu(3)

               else

                  iunow = iu0  ! this section added for walls
               endif
            else

               if (itab_ud(iunow)%iw(1) > 1) then

                  npoly = npoly + 1

                  itab_md(im)%iu(npoly) = iunow  ! NEW

                  itab_md(im)%iw(npoly) = itab_ud(iunow)%iw(1)
                  iunow = itab_ud(iunow)%iu(2)

               else

                  iunow = iu0  ! this section added for walls
               endif
            endif

            itab_md(im)%npoly = npoly

            if (npoly > 10) stop 'stop tri_neighbors npoly'

         enddo

! Define extra (wrap-around) iu value for itab_md            

! next line probably no longer needed, so commented out
!         if (npoly > 0) itab_md(im)%iu(npoly+1) = iunow

      endif

   enddo
enddo

return
end subroutine tri_neighbors

!===============================================================================

subroutine matrix_3x3(a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3,x1,x2,x3)

implicit none

real, intent(in) :: a11,a21,a31,a12,a22,a32,a13,a23,a33,b1,b2,b3
real, intent(out) :: x1,x2,x3

integer :: i
real, dimension(4,3) :: abr

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
      call rchange(abr(i,1),abr(i,2))
   enddo
endif

if (abs(abr(1,3)) > abs(abr(1,1))) then
   do i = 1,4
      call rchange(abr(i,1),abr(i,3))
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
      call rchange(abr(i,2),abr(i,3))
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
end subroutine matrix_3x3

!===============================================================================

subroutine matrix_2x2(a11,a21,a12,a22,b1,b2,x1,x2)

implicit none

real, intent(in) :: a11,a21,a12,a22,b1,b2
real, intent(out) :: x1,x2

integer :: i
real, dimension(3,2) :: abr

abr(1,1) = a11
abr(2,1) = a21
abr(3,1) = b1

abr(1,2) = a12
abr(2,2) = a22
abr(3,2) = b2

! Interchange rows if necessary so that first row has 
! largest (magnitude) element of first column

if (abs(abr(1,2)) > abs(abr(1,1))) then
   do i = 1,3
      call rchange(abr(i,1),abr(i,2))
   enddo
endif

! Add -abr(1,2)/abr(1,1) times first row to second row

do i = 2,3
   abr(i,2) = abr(i,2) - abr(1,2)/abr(1,1)*abr(i,1)
enddo

! Back substitution

x2 = abr(3,2) / abr(2,2)
x1 = (abr(3,1) - abr(2,1) * x2) / abr(1,1)

return
end subroutine matrix_2x2

!===============================================================================

subroutine rchange(r1,r2)

implicit none

real, intent(inout) :: r1,r2
real :: c

c = r1
r1 = r2
r2 = c

return
end subroutine rchange

!===============================================================================

subroutine matinv3x3(a11,a21,a31,a12,a22,a32,a13,a23,a33  &
                    ,b11,b21,b31,b12,b22,b32,b13,b23,b33)

implicit none
 
real, intent(in)  :: a11,a21,a31,a12,a22,a32,a13,a23,a33
real, intent(out) :: b11,b21,b31,b12,b22,b32,b13,b23,b33

real :: det

det = a11*(a22*a33-a23*a32)  &
    + a21*(a32*a13-a12*a33)  &
    + a31*(a12*a23-a22*a13)
    
if (abs(det) < 1.e-12) stop 'stop singular matrix'    
    
b11 = (a22*a33-a32*a23)/det
b21 = (a31*a23-a21*a33)/det
b31 = (a21*a32-a31*a22)/det
b12 = (a32*a13-a12*a33)/det
b22 = (a11*a33-a31*a13)/det
b32 = (a31*a12-a11*a32)/det
b13 = (a12*a23-a22*a13)/det
b23 = (a21*a13-a11*a23)/det
b33 = (a11*a22-a21*a12)/det

return
end subroutine matinv3x3

!===============================================================================

subroutine matinv2x2(a11,a21,a12,a22,b11,b21,b12,b22)

implicit none

real, intent(in)  :: a11,a21,a12,a22
real, intent(out) :: b11,b21,b12,b22

real :: det

det = a11*a22 - a21*a12

if (abs(det) < 1.e-12) stop 'stop singular matrix'    

b11 =  a22/det
b21 = -a21/det
b12 = -a12/det
b22 =  a11/det

return
end subroutine matinv2x2

!===============================================================================

subroutine unit_normal(px,py,pz,qx,qy,qz,rx,ry,rz,vx,vy,vz)

implicit none

real, intent(in)  :: px,py,pz,qx,qy,qz,rx,ry,rz
real, intent(out) :: vx,vy,vz

real :: v

! Find components (vx,vy,vz) of unit normal vector to plane
! that passes through points (p,q,r)

! Un-normalized V = (PQ) X (PR):

vx = (qy - py) * (rz - pz) - (qz - pz) * (ry - py)
vy = (qz - pz) * (rx - px) - (qx - px) * (rz - pz)
vz = (qx - px) * (ry - py) - (qy - py) * (rx - px)

! Magnitude of V:

v = sqrt(vx * vx + vy * vy + vz * vz)

! Normalized components:

vx = vx / v
vy = vy / v
vz = vz / v

return
end subroutine unit_normal




