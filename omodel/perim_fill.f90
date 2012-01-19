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
subroutine perim_fill(ngr,mrloo,kma,kua,kwa,jm,ju,npts,ngrp,igsize)

use mem_grid,   only: nrows,mrows
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa
integer, intent(in) :: npts
integer, intent(inout) :: jm(nrows+1,npts)
integer, intent(inout) :: ju(nrows,npts)
integer, intent(in) :: ngrp
integer, intent(inout) :: igsize(npts)

integer :: jmo(nrows+1),juo(nrows)  ! temp storage for output row

integer :: irow, ig, iu, iw  ! last 2 special

! Determine width of transition zone at each point on perimeter

do irow = 1,mrows
   do ig = 1,ngrp
   
      if (igsize(ig) == 1) then
      
         call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) )
      
      elseif (igsize(ig) == 2) then
      
         call perim_fill_cent2(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,ju(2,ig),jm(3,ig)          &
                              ,jmo(1),juo(1),jmo(2)       )
         
      elseif (igsize(ig) == 3) then
      
         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,jmo(1),juo(1),jmo(2)       )

         call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                              ,jm(2,ig),ju(2,ig),jm(3,ig) )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(3,ig),ju(3,ig),jm(4,ig) &
                              ,jmo(2),juo(2),jmo(3)       )
      
      elseif (igsize(ig) == 4) then
      
         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,jmo(1),juo(1),jmo(2)       )

         call perim_fill_cent2(ngr,mrloo,kma,kua,kwa      &
                              ,jm(2,ig),ju(2,ig),jm(3,ig) &
                              ,ju(3,ig),jm(4,ig)          &
                              ,jmo(2),juo(2),jmo(3)       )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(4,ig),ju(4,ig),jm(5,ig) &
                              ,jmo(3),juo(3),jmo(4)       )
      
      elseif (igsize(ig) == 5) then
      
         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,jmo(1),juo(1),jmo(2)       )

         call perim_fill_right(ngr,mrloo,kma,kua,kwa      &
                              ,jm(2,ig),ju(2,ig),jm(3,ig) &
                              ,jmo(2),juo(2),jmo(3)       )

         call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                              ,jm(3,ig),ju(3,ig),jm(4,ig) )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(4,ig),ju(4,ig),jm(5,ig) &
                              ,jmo(3),juo(3),jmo(4)       )

         call perim_fill_left (ngr,mrloo,kma,kua,kwa      &
                              ,jm(5,ig),ju(5,ig),jm(6,ig) &
                              ,jmo(4),juo(4),jmo(5)       )

      endif
      
      igsize(ig) = igsize(ig) - 1
      
      jm(1:igsize(ig)+1,ig) = jmo(1:igsize(ig)+1)
      ju(1:igsize(ig)  ,ig) = juo(1:igsize(ig)  )
      
   enddo
      
enddo

return
end subroutine perim_fill

!===============================================================================

subroutine perim_fill2(ngr,mrloo,kma,kua,kwa,jm,ju,npts,ngrp,igsize,nwdivg)

use mem_grid,   only: nrows,mrows
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa
integer, intent(in) :: npts
integer, intent(inout) :: jm(nrows+1,npts)
integer, intent(inout) :: ju(nrows,npts)
integer, intent(in) :: ngrp
integer, intent(in) :: igsize(npts)
integer, intent(in) :: nwdivg(npts)

integer :: jmo(nrows+1),juo(nrows)  ! temp storage for output row

integer :: irow, ig, iu, iw  ! last 2 special

! Determine width of transition zone at each point on perimeter

do ig = 1,ngrp
   
   if (igsize(ig) == 1) then
      
      call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                           ,jm(1,ig),ju(1,ig),jm(2,ig) )
      
   elseif (nwdivg(ig) == 3) then

      call perim_fill_cent2(ngr,mrloo,kma,kua,kwa      &
                           ,jm(1,ig),ju(1,ig),jm(2,ig) &
                           ,ju(2,ig),jm(3,ig)          &
                           ,jmo(1),juo(1),jmo(2)       )
         
      call perim_fill_cent1(ngr,mrloo,kma,kua,kwa      &
                           ,jmo(1),juo(1),jmo(2) )
      
   elseif (nwdivg(ig) == 2) then

      call perim_fill_convex2(ngr,mrloo,kma,kua,kwa      &
                             ,jm(1,ig),ju(1,ig),jm(2,ig) &
                             ,ju(2,ig),jm(3,ig))

   else ! nwdivg(ig) = 4

      call perim_fill_concave2(ngr,mrloo,kma,kua,kwa      &
                              ,jm(1,ig),ju(1,ig),jm(2,ig) &
                              ,ju(2,ig),jm(3,ig))
   endif

enddo

return
end subroutine perim_fill2

!===========================================================================

subroutine perim_fill3(nper,imper,iuper)

use mem_ijtabs, only: itab_ud, itab_wd, ltab_ud, nest_ud, nest_wd
use mem_grid,   only: xem, yem, zem
use misc_coms,  only: io6

implicit none

integer, intent(in) :: nper, imper(nper), iuper(nper) 

integer :: jm1, jm2, jm3, ju1, ju2, ju3, ku1, ku2, ku3

integer :: im5,im12,im13,im17,im18,im19,im20,im24

integer :: iu15,iu16,iu25,iu26,iu33,iu34,iu35,iu41,iu42,iu43
integer :: iu44,iu45,iu48,iu49,iu50,iu51

integer :: iw6,iw7,iw8,iw9,iw19,iw20,iw21,iw27,iw29,iw31

integer :: iper

do iper = 1,nper,3

! Inventory of CM mesh points

   jm1 = imper(iper)
   jm2 = imper(iper+1)
   jm3 = imper(iper+2)

   ju1 = iuper(iper)
   ju2 = iuper(iper+1)
   ju3 = iuper(iper+2)

   if (jm1 == ltab_ud(ju1)%im(1)) then
      iu41 = ju1
      iu42 = nest_ud(ju1)%iu
      iw27 = ltab_ud(ju1)%iw(1)
   else
      iu41 = nest_ud(ju1)%iu
      iu42 = ju1
      iw27 = ltab_ud(ju1)%iw(2)
   endif

   if (jm2 == ltab_ud(ju2)%im(1)) then

      iu49 = ltab_ud(ju2)%iu(1)
      iu50 = ltab_ud(ju2)%iu(2)
      iu34 = ltab_ud(ju2)%iu(3)
      iu35 = ltab_ud(ju2)%iu(4)
      iu48 = ltab_ud(ju2)%iu(5)
      iu51 = ltab_ud(ju2)%iu(8)

      iw6 = ltab_ud(ju2)%iw(5)
      iw9 = ltab_ud(ju2)%iw(6)
      iw29 = ltab_ud(ju2)%iw(1)
      iw20 = ltab_ud(ju2)%iw(2)

   else

      iu49 = ltab_ud(ju2)%iu(4)
      iu50 = ltab_ud(ju2)%iu(3)
      iu34 = ltab_ud(ju2)%iu(2)
      iu35 = ltab_ud(ju2)%iu(1)
      iu48 = ltab_ud(ju2)%iu(12)
      iu51 = ltab_ud(ju2)%iu(9)

      iw6 = ltab_ud(ju2)%iw(4)
      iw9 = ltab_ud(ju2)%iw(3)
      iw29 = ltab_ud(ju2)%iw(2)
      iw20 = ltab_ud(ju2)%iw(1)

   endif

   if (jm3 == ltab_ud(ju3)%im(1)) then
      iu44 = ju3
      iu45 = nest_ud(ju3)%iu
      iw31 = ltab_ud(ju3)%iw(1)
   else
      iu44 = nest_ud(ju3)%iu
      iu45 = ju3
      iw31 = ltab_ud(ju3)%iw(2)
   endif

   im17 = nest_ud(ju1)%im
   im18 = jm2
   im19 = jm3
   im20 = nest_ud(ju3)%im
   iu43 = ju2

   ku1 = nest_wd(iw6)%iu(1)
   ku2 = nest_wd(iw6)%iu(2)
   ku3 = nest_wd(iw6)%iu(3)

   if (itab_ud(ku1)%im(1) > 1 .and. itab_ud(ku1)%im(2) > 1) then
      iu25 = ku2
      iu15 = ku3
   elseif (itab_ud(ku2)%im(1) > 1 .and. itab_ud(ku2)%im(2) > 1) then
      iu25 = ku3
      iu15 = ku1
   elseif (itab_ud(ku3)%im(1) > 1 .and. itab_ud(ku3)%im(2) > 1) then
      iu25 = ku1
      iu15 = ku2
   endif

   if (itab_ud(iu15)%iw(1) == iw6) then
      iw7 = itab_ud(iu15)%iw(2)
   else
      iw7 = itab_ud(iu15)%iw(1)
   endif

   if (itab_ud(iu25)%iw(1) == iw6) then
      iw19 = itab_ud(iu25)%iw(2)
      im12 = itab_ud(iu25)%im(2) ! to reposition xyzem(12)
   else
      iw19 = itab_ud(iu25)%iw(1)
      im12 = itab_ud(iu25)%im(1) ! to reposition xyzem(12)
   endif

   ku1 = nest_wd(iw9)%iu(1)
   ku2 = nest_wd(iw9)%iu(2)
   ku3 = nest_wd(iw9)%iu(3)

   if (itab_ud(ku1)%im(1) > 1 .and. itab_ud(ku1)%im(2) > 1) then
      iu16 = ku2
      iu26 = ku3
   elseif (itab_ud(ku2)%im(1) > 1 .and. itab_ud(ku2)%im(2) > 1) then
      iu16 = ku3
      iu26 = ku1
   elseif (itab_ud(ku3)%im(1) > 1 .and. itab_ud(ku3)%im(2) > 1) then
      iu16 = ku1
      iu26 = ku2
   endif

   if (itab_ud(iu16)%iw(1) == iw9) then
      iw8 = itab_ud(iu16)%iw(2)
   else
      iw8 = itab_ud(iu16)%iw(1)
   endif

   if (itab_ud(iu26)%iw(1) == iw9) then
      iw21 = itab_ud(iu26)%iw(2)
      im13 = itab_ud(iu26)%im(1) ! to reposition xyzem(13)
   else
      iw21 = itab_ud(iu26)%iw(1)
      im13 = itab_ud(iu26)%im(2) ! to reposition xyzem(13)
   endif

   if (im18 == itab_ud(iu49)%im(1)) then
      im24 = itab_ud(iu49)%im(2)
   else
      im24 = itab_ud(iu49)%im(1)
   endif 

! Fill neighbor indices:

   if (itab_ud(iu15)%im(1) == 1) then
      itab_ud(iu15)%im(1) = im18
   else
      itab_ud(iu15)%im(2) = im18
   endif

   if (itab_ud(iu16)%im(1) == 1) then
      itab_ud(iu16)%im(1) = im18
   else
      itab_ud(iu16)%im(2) = im18
   endif

   if (itab_ud(iu25)%im(1) == 1) then
      itab_ud(iu25)%im(1) = im18
   else
      itab_ud(iu25)%im(2) = im18
   endif

   if (itab_ud(iu26)%im(1) == 1) then
      itab_ud(iu26)%im(1) = im18
   else
      itab_ud(iu26)%im(2) = im18
   endif

   if (itab_ud(iu34)%im(1) == im18) then 
      itab_ud(iu34)%iw(1) = iw8
      itab_ud(iu34)%iw(2) = iw7
      im5 = itab_ud(iu34)%im(2)
   else
      itab_ud(iu34)%iw(1) = iw7
      itab_ud(iu34)%iw(2) = iw8
      im5 = itab_ud(iu34)%im(1)
   endif

   if (itab_ud(iu35)%im(1) == im19) then
      itab_ud(iu35)%iw(2) = iw19
      itab_ud(iu35)%iw(1) = iw21
      itab_ud(iu35)%im(2) = im18
   else
      itab_ud(iu35)%iw(1) = iw19
      itab_ud(iu35)%iw(2) = iw21
      itab_ud(iu35)%im(1) = im18
   endif

   if (itab_ud(iu41)%im(2) == im17) then
      itab_ud(iu41)%iw(1) = iw27
   else
      itab_ud(iu41)%iw(2) = iw27
   endif

   if (itab_ud(iu42)%im(1) == im17) then
      itab_ud(iu42)%im(2) = im19
      itab_ud(iu42)%iw(1) = iw20
   else
      itab_ud(iu42)%im(1) = im19
      itab_ud(iu42)%iw(2) = iw20
   endif

   if (itab_ud(iu43)%im(2) == im19) then
      itab_ud(iu43)%im(1) = im24
   else
      itab_ud(iu43)%im(2) = im24
   endif

   if (itab_ud(iu44)%im(1) == im19) then
      itab_ud(iu44)%iw(1) = iw29
   else
      itab_ud(iu44)%iw(2) = iw29
   endif

   if (itab_ud(iu45)%im(1) == im20) then
      itab_ud(iu45)%iw(1) = iw31
   else
      itab_ud(iu45)%iw(2) = iw31
   endif

   if (itab_ud(iu48)%iw(2) == iw27) then
      itab_ud(iu48)%im(2) = im17
   else
      itab_ud(iu48)%im(1) = im17
   endif

   if (itab_ud(iu49)%im(2) == im24) then
      itab_ud(iu49)%im(1) = im17
      itab_ud(iu49)%iw(2) = iw20
   else
      itab_ud(iu49)%im(2) = im17
      itab_ud(iu49)%iw(1) = iw20
   endif

   if (itab_ud(iu50)%im(1) == im24) then
      itab_ud(iu50)%im(2) = im20
   else
      itab_ud(iu50)%im(1) = im20
   endif

   if (itab_ud(iu51)%iw(2) == iw31) then
      itab_ud(iu51)%im(1) = im20
   else
      itab_ud(iu51)%im(2) = im20
   endif

   if (itab_wd(iw8)%iu(1) == iu16) then
      itab_wd(iw8)%iu(3) = iu34
   elseif (itab_wd(iw8)%iu(2) == iu16) then
      itab_wd(iw8)%iu(1) = iu34
   else
      itab_wd(iw8)%iu(2) = iu34
   endif

   if (itab_wd(iw19)%iu(1) == iu25) then
      iu33 = itab_wd(iw19)%iu(2)
      itab_wd(iw19)%iu(3) = iu35
   elseif (itab_wd(iw19)%iu(2) == iu25) then
      iu33 = itab_wd(iw19)%iu(3)
      itab_wd(iw19)%iu(1) = iu35
   else
      iu33 = itab_wd(iw19)%iu(1)
      itab_wd(iw19)%iu(2) = iu35
   endif

! special case

   if (itab_ud(iu33)%iw(2) == iw19) then
      itab_ud(iu33)%im(2) = im19
   else
      itab_ud(iu33)%im(1) = im19
   endif

   if (itab_wd(iw20)%iu(1) == iu43) then
      itab_wd(iw20)%iu(2) = iu42
      itab_wd(iw20)%iu(3) = iu49
   elseif (itab_wd(iw20)%iu(2) == iu43) then
      itab_wd(iw20)%iu(3) = iu42
      itab_wd(iw20)%iu(1) = iu49
   else
      itab_wd(iw20)%iu(1) = iu42
      itab_wd(iw20)%iu(2) = iu49
   endif

   if (itab_wd(iw27)%iu(1) == iu48) then
      itab_wd(iw27)%iu(2) = iu41
   elseif (itab_wd(iw27)%iu(2) == iu48) then
      itab_wd(iw27)%iu(3) = iu41
   else
      itab_wd(iw27)%iu(1) = iu41
   endif

   if (itab_wd(iw29)%iu(1) == iu50) then
      itab_wd(iw29)%iu(2) = iu44
      itab_wd(iw29)%iu(3) = iu43
   elseif (itab_wd(iw29)%iu(2) == iu50) then
      itab_wd(iw29)%iu(3) = iu44
      itab_wd(iw29)%iu(1) = iu43
   else
      itab_wd(iw29)%iu(1) = iu44
      itab_wd(iw29)%iu(2) = iu43
   endif

   if (itab_wd(iw31)%iu(1) == iu51) then
      itab_wd(iw31)%iu(3) = iu45
   elseif (itab_wd(iw31)%iu(2) == iu51) then
      itab_wd(iw31)%iu(1) = iu45
   else
      itab_wd(iw31)%iu(2) = iu45
   endif

! NEW M locations

   xem(im19) = .5 * (xem(im24) + xem(im5))
   yem(im19) = .5 * (yem(im24) + yem(im5))
   zem(im19) = .5 * (zem(im24) + zem(im5))

   xem(im18) = .5 * (xem(im19) + xem(im5))
   yem(im18) = .5 * (yem(im19) + yem(im5))
   zem(im18) = .5 * (zem(im19) + zem(im5))

   xem(im17) = .75 * xem(im17) + .25 * xem(im19)
   yem(im17) = .75 * yem(im17) + .25 * yem(im19)
   zem(im17) = .75 * zem(im17) + .25 * zem(im19)

   xem(im20) = .75 * xem(im20) + .25 * xem(im19)
   yem(im20) = .75 * yem(im20) + .25 * yem(im19)
   zem(im20) = .75 * zem(im20) + .25 * zem(im19)

   xem(im12) = .833 * xem(im12) + .167 * xem(im18)
   yem(im12) = .833 * yem(im12) + .167 * yem(im18)
   zem(im12) = .833 * zem(im12) + .167 * zem(im18)

   xem(im13) = .833 * xem(im13) + .167 * xem(im18)
   yem(im13) = .833 * yem(im13) + .167 * yem(im18)
   zem(im13) = .833 * zem(im13) + .167 * zem(im18)

enddo ! iper

end subroutine perim_fill3

!===========================================================================

subroutine perim_fill_cent1(ngr,mrloo,kma,kua,kwa,jm1,ju1,jm2)

use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

integer, intent(in) :: jm1,jm2  ! Original border M indices
integer, intent(in) :: ju1      ! Original border U index

integer :: ju4,ju7
integer :: iw1,iw2
integer :: iu1,iu2,iu3,iu4,iu5
integer :: im1,im2,im3,im4

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2

iu4 = kua + 1  ! Newly added point
iw2 = kwa + 1  ! Newly added point

! Increment indices for newly added points

kua = kua + 1
kwa = kwa + 1

! Determine orientation of positive ju1 direction

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
   iw1 = ltab_ud(ju1)%iw(1)
   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

   itab_ud(iu1)%iw(1) = iw1
   itab_ud(iu2)%iw(1) = iw2
else                              ! Positive ju1 points outward
   iw1 = ltab_ud(ju1)%iw(2)
   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = iw1
   itab_ud(iu2)%iw(2) = iw2
endif

if (ju1 == ltab_wd(iw1)%iu(1)) then
   iu3 = ltab_wd(iw1)%iu(2)
   iu5 = ltab_wd(iw1)%iu(3)
elseif (ju1 == ltab_wd(iw1)%iu(2)) then
   iu3 = ltab_wd(iw1)%iu(3)
   iu5 = ltab_wd(iw1)%iu(1)
else
   iu3 = ltab_wd(iw1)%iu(1)
   iu5 = ltab_wd(iw1)%iu(2)
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu3

if (iw1 == itab_ud(iu5)%iw(1)) then
   im4 = itab_ud(iu5)%im(2)
   itab_ud(iu5)%iw(1) = iw2
else
   im4 = itab_ud(iu5)%im(1)
   itab_ud(iu5)%iw(2) = iw2
endif

itab_ud(iu4)%im(1) = im2
itab_ud(iu4)%im(2) = im4
itab_ud(iu4)%iw(1) = iw1
itab_ud(iu4)%iw(2) = iw2

itab_wd(iw1)%iu(1) = iu1
itab_wd(iw1)%iu(2) = iu3
itab_wd(iw1)%iu(3) = iu4

itab_wd(iw2)%iu(1) = iu2
itab_wd(iw2)%iu(2) = iu4
itab_wd(iw2)%iu(3) = iu5

! Fill mrl for new W points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5

itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig
         
! Fill loop indices for new U points

itab_ud(iu4)%iup = iu4
call udloops('f',iu4, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu4,21,22,23, 0, 0, 0, 0, 0, 0, 0)

return

5 continue

write(io6,*) ''
write(io6, '(1x,A)')    'Error in subroutine perim_fill_cent1.'
write(io6, '(1x,A,I0)') 'Current nested grid ', ngr
write(io6, '(1x,A)')    'crosses pre-existing grid boundary.'
write(io6, '(1x,A,I0)') 'iw1 = ',iw1
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_cent1

!===========================================================================

subroutine perim_fill_cent2(ngr,mrloo,kma,kua,kwa,  &
                           jm1,ju1,jm2,ju2,jm3,jmo1,juo1,jmo2)

use mem_grid, only:  xem, yem, zem
use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa   ! Index values for latest added points

integer, intent(in) :: jm1,jm2,jm3   ! Original M indices

integer, intent(in) :: ju1,ju2       ! Original U indices

integer, intent(out) :: jmo1,juo1,jmo2  ! Original M,U indices for next row 

integer :: ju4,ju7

integer :: iw1,iw2,iw3,iw4,iw5,iw6          ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8  ! Temporary new U indices
integer :: iu9,iu10,iu11,iu12,iu13          ! Temporary new U indices
integer :: im1,im2,im3,im4,im5,im6,im7,im8  ! Temporary new M indices

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2
im4 = nest_ud(ju2)%im
im5 = jm3

im7 = kma + 1  ! Newly added point

iu6  = kua + 1  ! Newly added point
iu8  = kua + 2  ! Newly added point
iu10 = kua + 3  ! Newly added point

iw2 = kwa + 1  ! Newly added point
iw4 = kwa + 2  ! Newly added point
iw6 = kwa + 3  ! Newly added point

! Determine orientation of positive ju1 direction to get iw1 (=jw1)

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
   iw1 = ltab_ud(ju1)%iw(1)
   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

!   itab_ud(iu1)%iw(1) = iw1  ! not needed?
else                              ! Positive ju1 points outward
   iw1 = ltab_ud(ju1)%iw(2)
   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = iw1
endif

! Determine neighbors of iw1 (=jw1)

if (ju1 == ltab_wd(iw1)%iu(1)) then
   iu5 = ltab_wd(iw1)%iu(2)
   ju4 = ltab_wd(iw1)%iu(3)
elseif (ju1 == ltab_wd(iw1)%iu(2)) then
   iu5 = ltab_wd(iw1)%iu(3)
   ju4 = ltab_wd(iw1)%iu(1)
else
   iu5 = ltab_wd(iw1)%iu(1)
   ju4 = ltab_wd(iw1)%iu(2)
endif

iu7 = ju4

! Determine iw3 (=jw2)

if (jm2 == ltab_ud(ju4)%im(1)) then  ! Positive ju4 points to jw2
   im6 = ltab_ud(ju4)%im(2)
   iw3 = ltab_ud(ju4)%iw(2)
else                              ! Positive ju4 points to jw1
   im6 = ltab_ud(ju4)%im(1)
   iw3 = ltab_ud(ju4)%iw(1)
endif

! Revisit orientation of positive ju1 direction now that iw3 is known

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
   itab_ud(iu2)%iw(1) = iw3
else                              ! Positive ju1 points outward
   itab_ud(iu2)%iw(2) = iw3
endif

! Determine orientation of positive ju2 direction to get iw5 (=jw3)

if (jm2 == ltab_ud(ju2)%im(1)) then  ! Positive ju2 points inward
   iw5 = ltab_ud(ju2)%iw(1)

   iu3 = ju2
   iu4 = nest_ud(ju2)%iu

   itab_ud(iu3)%iw(1) = iw4
   itab_ud(iu4)%iw(1) = iw6
else                              ! Positive ju2 points outward
   iw5 = ltab_ud(ju2)%iw(2)

   iu3 = nest_ud(ju2)%iu
   iu4 = ju2

   itab_ud(iu3)%iw(2) = iw4
   itab_ud(iu4)%iw(2) = iw6
endif

! Determine neighbors of iw5 (=jw3)

if (ju2 == ltab_wd(iw5)%iu(1)) then
   iu9 = ltab_wd(iw5)%iu(2)
   iu11 = ltab_wd(iw5)%iu(3)
elseif (ju2 == ltab_wd(iw5)%iu(2)) then
   iu9 = ltab_wd(iw5)%iu(3)
   iu11 = ltab_wd(iw5)%iu(1)
else
   iu9 = ltab_wd(iw5)%iu(1)
   iu11 = ltab_wd(iw5)%iu(2)
endif

! Determine ju7 as neighbor of jw2 (=iw3)

if (ju4 == ltab_wd(iw3)%iu(1)) then
   ju7 = ltab_wd(iw3)%iu(2)
elseif (ju4 == ltab_wd(iw3)%iu(2)) then
   ju7 = ltab_wd(iw3)%iu(3)
else
   ju7 = ltab_wd(iw3)%iu(1)
endif

! Determine orientation of positive ju7 direction

if (im6 == ltab_ud(ju7)%im(1)) then  ! Positive ju7 points inward
   im8 = ltab_ud(ju7)%im(2)
   iu12 = ju7
   iu13 = kua + 4  ! Newly added point

   itab_ud(iu12)%im(1) = im6 ! not needed?
   itab_ud(iu12)%im(2) = im7
   itab_ud(iu12)%iw(2) = iw2

   itab_ud(iu13)%im(1) = im7
   itab_ud(iu13)%im(2) = im8
   itab_ud(iu13)%iw(2) = iw5
else
   im8 = ltab_ud(ju7)%im(1)
   iu12 = kua + 4  ! Newly added point
   iu13 = ju7

   itab_ud(iu12)%im(1) = im7
   itab_ud(iu12)%im(2) = im6
   itab_ud(iu12)%iw(1) = iw2

   itab_ud(iu13)%im(1) = im8
   itab_ud(iu13)%im(2) = im7
   itab_ud(iu13)%iw(1) = iw5
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu5

itab_ud(iu6)%im(1) = im2
itab_ud(iu6)%im(2) = im6
itab_ud(iu6)%iw(1) = iw1
itab_ud(iu6)%iw(2) = iw2

itab_ud(iu7)%im(1) = im2
itab_ud(iu7)%im(2) = im7
itab_ud(iu7)%iw(1) = iw2
itab_ud(iu7)%iw(2) = iw3

itab_ud(iu8)%im(1) = im3
itab_ud(iu8)%im(2) = im7
itab_ud(iu8)%iw(1) = iw3
itab_ud(iu8)%iw(2) = iw4

itab_ud(iu9)%im(1) = im4
itab_ud(iu9)%im(2) = im7
itab_ud(iu9)%iw(1) = iw4
itab_ud(iu9)%iw(2) = iw5

itab_ud(iu10)%im(1) = im4
itab_ud(iu10)%im(2) = im8
itab_ud(iu10)%iw(1) = iw5
itab_ud(iu10)%iw(2) = iw6

if (iw5 == ltab_ud(iu11)%iw(1)) then
   itab_ud(iu11)%iw(1) = iw6
else
   itab_ud(iu11)%iw(2) = iw6
endif

itab_wd(iw1)%iu(1) = iu1
itab_wd(iw1)%iu(2) = iu5
itab_wd(iw1)%iu(3) = iu6

itab_wd(iw2)%iu(1) = iu6
itab_wd(iw2)%iu(2) = iu12
itab_wd(iw2)%iu(3) = iu7

itab_wd(iw3)%iu(1) = iu2
itab_wd(iw3)%iu(2) = iu7
itab_wd(iw3)%iu(3) = iu8

itab_wd(iw4)%iu(1) = iu3
itab_wd(iw4)%iu(2) = iu8
itab_wd(iw4)%iu(3) = iu9

itab_wd(iw5)%iu(1) = iu9
itab_wd(iw5)%iu(2) = iu13
itab_wd(iw5)%iu(3) = iu10

itab_wd(iw6)%iu(1) = iu4
itab_wd(iw6)%iu(2) = iu10
itab_wd(iw6)%iu(3) = iu11

! Fill earth coordinates for new M point

xem(im7) = .5 * (xem(im6) + xem(im8))
yem(im7) = .5 * (yem(im6) + yem(im8))
zem(im7) = .5 * (zem(im6) + zem(im8))

! Fill mrl for new W points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5
if (itab_wd(iw3)%mrlw /= mrloo) go to 5
if (itab_wd(iw5)%mrlw /= mrloo) go to 5

itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig
         
itab_wd(iw4)%mrlw      = itab_wd(iw3)%mrlw
itab_wd(iw4)%mrlw_orig = itab_wd(iw3)%mrlw_orig
         
itab_wd(iw6)%mrlw      = itab_wd(iw5)%mrlw
itab_wd(iw6)%mrlw_orig = itab_wd(iw5)%mrlw_orig
         
! Fill loop indices for new U points

itab_ud(iu6)%iup = iu6
call udloops('f',iu6, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu6,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu8)%iup = iu8
call udloops('f',iu8, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu8,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu10)%iup = iu10
call udloops('f',iu10, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu10,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(kua+4)%iup = kua+4
call udloops('f',kua+4, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',kua+4,21,22,23, 0, 0, 0, 0, 0, 0, 0)

! Set JM and JU group indices in preparation for next row

jmo1 = im6
jmo2 = im8
juo1 = ju7

nest_ud(ju7)%im = im7
nest_ud(ju7)%iu = kua + 4

! Increment indices for newly added points

kma = kma + 1
kua = kua + 4
kwa = kwa + 3

return

5 continue

write(io6,*) 'In subroutine perim_fill_cent2, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_cent2

!===========================================================================

subroutine perim_fill_right(ngr,mrloo,kma,kua,kwa,jm1,ju1,jm2,jmo1,juo1,jmo2)

use mem_grid,   only:  xem, yem, zem
use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

integer, intent(in) :: jm1,jm2  ! Original M indices

integer, intent(in) :: ju1  ! Original U indices

integer, intent(out) :: jmo1,juo1,jmo2  ! Original M,U indices for next row 

integer :: ju3,ju5

integer :: iw1,iw2,iw3,iw4                      ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9  ! Temporary new U indices
integer :: im1,im2,im3,im4,im5,im6              ! Temporary new M indices

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2

im5 = kma + 1  ! Newly added point

iu4  = kua + 1  ! Newly added point
iu6  = kua + 2  ! Newly added point

iw2 = kwa + 1  ! Newly added point
iw4 = kwa + 2  ! Newly added point

! Determine orientation of positive ju1 direction to get iw1 (=jw1)

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
   iw1 = ltab_ud(ju1)%iw(1)
   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

   itab_ud(iu1)%iw(1) = iw1 ! not needed?
else                              ! Positive ju1 points outward
   iw1 = ltab_ud(ju1)%iw(2)
   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = iw1
endif

! Determine neighbors of iw1 (=jw1)

if (ju1 == ltab_wd(iw1)%iu(1)) then
   iu3 = ltab_wd(iw1)%iu(2)
   ju3 = ltab_wd(iw1)%iu(3)
elseif (ju1 == ltab_wd(iw1)%iu(2)) then
   iu3 = ltab_wd(iw1)%iu(3)
   ju3 = ltab_wd(iw1)%iu(1)
else
   iu3 = ltab_wd(iw1)%iu(1)
   ju3 = ltab_wd(iw1)%iu(2)
endif

iu5 = ju3

! Determine iw3 (=jw2)

if (jm2 == ltab_ud(ju3)%im(1)) then  ! Positive ju3 points to jw2
   im4 = ltab_ud(ju3)%im(2)
   iw3 = ltab_ud(ju3)%iw(2)
else                              ! Positive ju3 points to jw1
   im4 = ltab_ud(ju3)%im(1)
   iw3 = ltab_ud(ju3)%iw(1)
endif

! Revisit orientation of positive ju1 direction now that iw3 is known

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
   itab_ud(iu2)%iw(1) = iw3
else                              ! Positive ju1 points outward
   itab_ud(iu2)%iw(2) = iw3
endif

! Determine neighbors of iw3 (=jw2)

if (ju3 == ltab_wd(iw3)%iu(1)) then
   ju5 = ltab_wd(iw3)%iu(2)
   iu7 = ltab_wd(iw3)%iu(3)
elseif (ju3 == ltab_wd(iw3)%iu(2)) then
   ju5 = ltab_wd(iw3)%iu(3)
   iu7 = ltab_wd(iw3)%iu(1)
else
   ju5 = ltab_wd(iw3)%iu(1)
   iu7 = ltab_wd(iw3)%iu(2)
endif

! Determine orientation of positive ju5 direction

if (im4 == ltab_ud(ju5)%im(1)) then  ! Positive ju5 points inward
   im6 = ltab_ud(ju5)%im(2)
   iu8 = ju5
   iu9 = kua + 3  ! Newly added point

   itab_ud(iu8)%im(1) = im4
   itab_ud(iu8)%im(2) = im5
   itab_ud(iu8)%iw(2) = iw2

   itab_ud(iu9)%im(1) = im5
   itab_ud(iu9)%im(2) = im6
   itab_ud(iu9)%iw(2) = iw4
else
   im6 = ltab_ud(ju5)%im(1)
   iu8 = kua + 3  ! Newly added point
   iu9 = ju5

   itab_ud(iu8)%im(1) = im5
   itab_ud(iu8)%im(2) = im4
   itab_ud(iu8)%iw(1) = iw2

   itab_ud(iu9)%im(1) = im6
   itab_ud(iu9)%im(2) = im5
   itab_ud(iu9)%iw(1) = iw4
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu3

itab_ud(iu4)%im(1) = im2
itab_ud(iu4)%im(2) = im4
itab_ud(iu4)%iw(1) = iw1
itab_ud(iu4)%iw(2) = iw2

itab_ud(iu5)%im(1) = im2
itab_ud(iu5)%im(2) = im5
itab_ud(iu5)%iw(1) = iw2
itab_ud(iu5)%iw(2) = iw3

itab_ud(iu6)%im(1) = im3
itab_ud(iu6)%im(2) = im5
itab_ud(iu6)%iw(1) = iw3
itab_ud(iu6)%iw(2) = iw4

if (iw3 == ltab_ud(iu7)%iw(1)) then
   itab_ud(iu7)%iw(1) = iw4
else
   itab_ud(iu7)%iw(2) = iw4
endif

itab_wd(iw1)%iu(1) = iu1
itab_wd(iw1)%iu(2) = iu3
itab_wd(iw1)%iu(3) = iu4

itab_wd(iw2)%iu(1) = iu4
itab_wd(iw2)%iu(2) = iu8
itab_wd(iw2)%iu(3) = iu5

itab_wd(iw3)%iu(1) = iu2
itab_wd(iw3)%iu(2) = iu5
itab_wd(iw3)%iu(3) = iu6

itab_wd(iw4)%iu(1) = iu6
itab_wd(iw4)%iu(2) = iu9
itab_wd(iw4)%iu(3) = iu7

! Fill earth coordinates for new M point

xem(im5) = .5 * (xem(im4) + xem(im6))
yem(im5) = .5 * (yem(im4) + yem(im6))
zem(im5) = .5 * (zem(im4) + zem(im6))

! Fill mrlfor new W points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5
if (itab_wd(iw3)%mrlw /= mrloo) go to 5

itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig
         
itab_wd(iw4)%mrlw      = itab_wd(iw3)%mrlw
itab_wd(iw4)%mrlw_orig = itab_wd(iw3)%mrlw_orig
         
! Fill loop indices for new U points

itab_ud(iu4)%iup = iu4
call udloops('f',iu4, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu4,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu6)%iup = iu6
call udloops('f',iu6, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu6,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(kua+3)%iup = kua+3
call udloops('f',kua+3, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',kua+3,21,22,23, 0, 0, 0, 0, 0, 0, 0)

! Set JM and JU group indices in preparation for next row

jmo1 = im4
jmo2 = im6
juo1 = ju5

nest_ud(ju5)%im = im5
nest_ud(ju5)%iu = kua + 3

! Increment indices for newly added points

kma = kma + 1
kua = kua + 3
kwa = kwa + 2

return

5 continue

write(io6,*) 'In subroutine perim_fill_right, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_right

!===========================================================================

subroutine perim_fill_left(ngr,mrloo,kma,kua,kwa,jm1,ju1,jm2,jmo1,juo1,jmo2)

use mem_grid,   only:  xem, yem, zem
use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa  ! Index values for latest added points

integer, intent(in) :: jm1,jm2  ! Original M indices
integer, intent(in) :: ju1      ! Original U indices

integer, intent(out) :: jmo1,juo1,jmo2  ! Original M,U indices for next row 


integer :: ju3,ju5

integer :: iw1,iw2,iw3,iw4                      ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9  ! Temporary new U indices
integer :: im1,im2,im3,im4,im5,im6              ! Temporary new M indices

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2

im5 = kma + 1  ! Newly added point

iu4  = kua + 1  ! Newly added point
iu6  = kua + 2  ! Newly added point

iw2 = kwa + 1  ! Newly added point
iw4 = kwa + 2  ! Newly added point

! Determine orientation of positive ju1 direction to get iw3 (=jw2)

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points inward
   iw3 = ltab_ud(ju1)%iw(1)
   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

   itab_ud(iu1)%iw(1) = iw2
   itab_ud(iu2)%iw(1) = iw4
else                              ! Positive ju1 points outward
   iw3 = ltab_ud(ju1)%iw(2)
   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = iw2
   itab_ud(iu2)%iw(2) = iw4
endif

! Determine neighbors of iw3 (=jw2)

if (ju1 == ltab_wd(iw3)%iu(1)) then
   ju3 = ltab_wd(iw3)%iu(2)
   iu7 = ltab_wd(iw3)%iu(3)
elseif (ju1 == ltab_wd(iw3)%iu(2)) then
   ju3 = ltab_wd(iw3)%iu(3)
   iu7 = ltab_wd(iw3)%iu(1)
else
   ju3 = ltab_wd(iw3)%iu(1)
   iu7 = ltab_wd(iw3)%iu(2)
endif

iu5 = ju3

! Determine iw1 (=jw1)

if (jm1 == ltab_ud(ju3)%im(1)) then  ! Positive ju3 points to jw2
   im6 = ltab_ud(ju3)%im(2)
   iw1 = ltab_ud(ju3)%iw(1)
else                              ! Positive ju3 points to jw1
   im6 = ltab_ud(ju3)%im(1)
   iw1 = ltab_ud(ju3)%iw(2)
endif

! Determine neighbors of iw1 (=jw1)

if (ju3 == ltab_wd(iw1)%iu(1)) then
   iu3 = ltab_wd(iw1)%iu(2)
   ju5 = ltab_wd(iw1)%iu(3)
elseif (ju3 == ltab_wd(iw1)%iu(2)) then
   iu3 = ltab_wd(iw1)%iu(3)
   ju5 = ltab_wd(iw1)%iu(1)
else
   iu3 = ltab_wd(iw1)%iu(1)
   ju5 = ltab_wd(iw1)%iu(2)
endif

! Determine orientation of positive ju5 direction

if (im6 == ltab_ud(ju5)%im(2)) then  ! Positive ju5 points inward
   im4 = ltab_ud(ju5)%im(1)
   iu8 = ju5
   iu9 = kua + 3  ! Newly added point

   itab_ud(iu8)%im(1) = im4
   itab_ud(iu8)%im(2) = im5
   itab_ud(iu8)%iw(2) = iw1

   itab_ud(iu9)%im(1) = im5
   itab_ud(iu9)%im(2) = im6
   itab_ud(iu9)%iw(2) = iw3
else
   im4 = ltab_ud(ju5)%im(2)
   iu8 = kua + 3  ! Newly added point
   iu9 = ju5

   itab_ud(iu8)%im(1) = im5
   itab_ud(iu8)%im(2) = im4
   itab_ud(iu8)%iw(1) = iw1

   itab_ud(iu9)%im(1) = im6
   itab_ud(iu9)%im(2) = im5
   itab_ud(iu9)%iw(1) = iw3
endif

! Fill remaining neighbor indices for U and W points

! nothing needed for iu3

itab_ud(iu4)%im(1) = im1
itab_ud(iu4)%im(2) = im5
itab_ud(iu4)%iw(1) = iw1
itab_ud(iu4)%iw(2) = iw2

itab_ud(iu5)%im(1) = im2
itab_ud(iu5)%im(2) = im5
itab_ud(iu5)%iw(1) = iw2
itab_ud(iu5)%iw(2) = iw3

itab_ud(iu6)%im(1) = im2
itab_ud(iu6)%im(2) = im6
itab_ud(iu6)%iw(1) = iw3
itab_ud(iu6)%iw(2) = iw4

if (iw3 == ltab_ud(iu7)%iw(1)) then
   itab_ud(iu7)%iw(1) = iw4
else
   itab_ud(iu7)%iw(2) = iw4
endif

itab_wd(iw1)%iu(1) = iu3
itab_wd(iw1)%iu(2) = iu8
itab_wd(iw1)%iu(3) = iu4

itab_wd(iw2)%iu(1) = iu1
itab_wd(iw2)%iu(2) = iu4
itab_wd(iw2)%iu(3) = iu5

itab_wd(iw3)%iu(1) = iu5
itab_wd(iw3)%iu(2) = iu9
itab_wd(iw3)%iu(3) = iu6

itab_wd(iw4)%iu(1) = iu2
itab_wd(iw4)%iu(2) = iu6
itab_wd(iw4)%iu(3) = iu7

! Fill earth coordinates for new M point

xem(im5) = .5 * (xem(im4) + xem(im6))
yem(im5) = .5 * (yem(im4) + yem(im6))
zem(im5) = .5 * (zem(im4) + zem(im6))

! Fill mrlfor new W points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5
if (itab_wd(iw3)%mrlw /= mrloo) go to 5

itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig
         
itab_wd(iw4)%mrlw      = itab_wd(iw3)%mrlw
itab_wd(iw4)%mrlw_orig = itab_wd(iw3)%mrlw_orig
         
! Fill loop indices for new U points

itab_ud(iu4)%iup = iu4
call udloops('f',iu4, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu4,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu6)%iup = iu6
call udloops('f',iu6, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu6,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(kua+3)%iup = kua+3
call udloops('f',kua+3, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',kua+3,21,22,23, 0, 0, 0, 0, 0, 0, 0)

! Set JM and JU group indices in preparation for next row

jmo1 = im4
jmo2 = im6
juo1 = ju5

nest_ud(ju5)%im = im5
nest_ud(ju5)%iu = kua + 3

! Increment indices for newly added points

kma = kma + 1
kua = kua + 3
kwa = kwa + 2

return

5 continue

write(io6,*) 'In subroutine perim_fill_left, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_left

!===========================================================================

subroutine perim_fill_concave2(ngr,mrloo,kma,kua,kwa,  &
                               jm1,ju1,jm2,ju2,jm3)

use mem_grid, only:  xem, yem, zem
use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa   ! Index values for latest added points

integer, intent(in) :: jm1,jm2,jm3   ! Original M indices

integer, intent(in) :: ju1,ju2       ! Original U indices

integer :: jm4
integer :: ju3,ju4,ju5
integer :: jw1,jw2

integer :: im1,im2,im3,im4,im5,im6             ! Temporary new M indices
integer :: iw1,iw2,iw3,iw4                     ! Temporary new W indices
integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9 ! Temporary new U indices

! Newly added edges

iu6 = kua + 1
iu8 = kua + 2

! Newly added areas

iw2 = kwa + 1
iw3 = kwa + 2

! Increment indices for newly added points

kua = kua + 2
kwa = kwa + 2

! Determine orientation of positive ju1 direction to get iw1 (=jw1)

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points into fine mesh
   jw1 = ltab_ud(ju1)%iw(1)
   jw2 = ltab_ud(ju1)%iw(4)

   ju3 = ltab_ud(ju1)%iu(1)
   ju4 = ltab_ud(ju1)%iu(2)
   ju5 = ltab_ud(ju1)%iu(7)

   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

   itab_ud(iu1)%iw(1) = jw1  ! not needed?
   itab_ud(iu2)%iw(1) = iw3
else                              ! Positive ju1 points out from fine mesh
   jw1 = ltab_ud(ju1)%iw(2)
   jw2 = ltab_ud(ju1)%iw(5)

   ju3 = ltab_ud(ju1)%iu(4)
   ju4 = ltab_ud(ju1)%iu(3)
   ju5 = ltab_ud(ju1)%iu(10)

   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = jw1
   itab_ud(iu2)%iw(2) = iw3
endif

!write(io6,'(a,20i7)') 'pf11 ',jm1,jm2,jm3,ju1,ju2, &
!   ltab_ud(ju1)%im(1),ltab_ud(ju1)%im(2),ltab_ud(ju2)%im(1),ltab_ud(ju2)%im(2), &
!   ltab_ud(ju1)%iw(1),ltab_ud(ju1)%iw(2),ltab_ud(ju2)%iw(1),ltab_ud(ju2)%iw(2), &
!   nest_ud(ju1)%iu

! Determine orientation of positive ju2 direction to get iw4 (=jw2)

if (jm2 == ltab_ud(ju2)%im(1)) then  ! Positive ju2 points inward
   iu3 = ju2
   iu4 = nest_ud(ju2)%iu

   itab_ud(iu3)%iw(1) = iw3
   itab_ud(iu4)%iw(1) = jw2
else                              ! Positive ju2 points outward
   iu3 = nest_ud(ju2)%iu
   iu4 = ju2

   itab_ud(iu3)%iw(2) = iw3  ! not needed?
   itab_ud(iu4)%iw(2) = jw2
endif

! Determine neighbors of ju4

if (jm2 == ltab_ud(ju4)%im(1)) then
   jm4 = ltab_ud(ju4)%im(2)
else
   jm4 = ltab_ud(ju4)%im(1)
endif

! Assign remaining I indices that already existed as J indices

! (iu1:iu4 already assigned above)

! Assign new I indices for all M, W points, and for U points that need no checks

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2
im4 = nest_ud(ju2)%im
im5 = jm3
im6 = jm4

iu5 = ju3
iu7 = ju4
iu9 = ju5

iw1 = jw1
iw4 = jw2

!print*, ' '
!write(io6,'(a,20i7)') 'c2a.1 ',im1,im2,im3,im4,im5,im6
!write(io6,'(a,20i7)') 'c2a.2 ',iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9
!write(io6,'(a,20i7)') 'c2a.3 ',iw1,iw2,iw3,iw4
!write(io6,'(a,20i7)') 'c2a.4  ',itab_ud(iu1)%im(1),itab_ud(iu1)%im(2), &
!                                itab_ud(iu1)%iw(1),itab_ud(iu1)%iw(2)
!write(io6,'(a,20i7)') 'c2a.5  ',itab_ud(iu2)%im(1),itab_ud(iu2)%im(2), &
!                                itab_ud(iu2)%iw(1),itab_ud(iu2)%iw(2)
!write(io6,'(a,20i7)') 'c2a.6  ',itab_ud(iu3)%im(1),itab_ud(iu3)%im(2), &
!                                itab_ud(iu3)%iw(1),itab_ud(iu3)%iw(2)
!write(io6,'(a,20i7)') 'c2a.7  ',itab_ud(iu4)%im(1),itab_ud(iu4)%im(2), &
!                                itab_ud(iu4)%iw(1),itab_ud(iu4)%iw(2)
!write(io6,'(a,20i7)') 'c2a.8  ',itab_ud(iu5)%im(1),itab_ud(iu5)%im(2), &
!                                itab_ud(iu5)%iw(1),itab_ud(iu5)%iw(2)
!write(io6,'(a,20i7)') 'c2a.9  ',itab_ud(iu6)%im(1),itab_ud(iu6)%im(2), &
!                                itab_ud(iu6)%iw(1),itab_ud(iu6)%iw(2)
!write(io6,'(a,20i7)') 'c2a.10 ',itab_ud(iu7)%im(1),itab_ud(iu7)%im(2), &
!                                itab_ud(iu7)%iw(1),itab_ud(iu7)%iw(2)
!write(io6,'(a,20i7)') 'c2a.11 ',itab_ud(iu8)%im(1),itab_ud(iu8)%im(2), &
!                                itab_ud(iu8)%iw(1),itab_ud(iu8)%iw(2)
!write(io6,'(a,20i7)') 'c2a.12 ',itab_ud(iu9)%im(1),itab_ud(iu9)%im(2), &
!                                itab_ud(iu9)%iw(1),itab_ud(iu9)%iw(2)
!write(io6,'(a,20i7)') 'c2a.13 ',itab_wd(iw1)%iu(1),itab_wd(iw1)%iu(2), &
!                                itab_wd(iw1)%iu(3)
!write(io6,'(a,20i7)') 'c2a.14 ',itab_wd(iw2)%iu(1),itab_wd(iw2)%iu(2), &
!                                itab_wd(iw2)%iu(3)
!write(io6,'(a,20i7)') 'c2a.15 ',itab_wd(iw3)%iu(1),itab_wd(iw3)%iu(2), &
!                                itab_wd(iw3)%iu(3)
!write(io6,'(a,20i7)') 'c2a.16 ',itab_wd(iw4)%iu(1),itab_wd(iw4)%iu(2), &
!                                itab_wd(iw4)%iu(3)
!print*, ' '



! Fill remaining neighbor indices for U and W points

! nothing needed for iu5

itab_ud(iu6)%im(1) = im2
itab_ud(iu6)%im(2) = im6
itab_ud(iu6)%iw(1) = iw1
itab_ud(iu6)%iw(2) = iw2

itab_ud(iu7)%im(1) = im2
itab_ud(iu7)%im(2) = im4
itab_ud(iu7)%iw(1) = iw2
itab_ud(iu7)%iw(2) = iw3

itab_ud(iu8)%im(1) = im4
itab_ud(iu8)%im(2) = im6
itab_ud(iu8)%iw(1) = iw2
itab_ud(iu8)%iw(2) = iw4

! nothing needed for iu9

itab_wd(iw1)%iu(1) = iu1
itab_wd(iw1)%iu(2) = iu5
itab_wd(iw1)%iu(3) = iu6

itab_wd(iw2)%iu(1) = iu6
itab_wd(iw2)%iu(2) = iu8
itab_wd(iw2)%iu(3) = iu7

itab_wd(iw3)%iu(1) = iu2
itab_wd(iw3)%iu(2) = iu7
itab_wd(iw3)%iu(3) = iu3

itab_wd(iw4)%iu(1) = iu4
itab_wd(iw4)%iu(2) = iu8
itab_wd(iw4)%iu(3) = iu9

! Fill mrl for new W points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5
if (itab_wd(iw4)%mrlw /= mrloo) go to 5

itab_wd(iw2)%mrlw      = itab_wd(iw1)%mrlw
itab_wd(iw2)%mrlw_orig = itab_wd(iw1)%mrlw_orig
         
itab_wd(iw3)%mrlw      = itab_wd(iw1)%mrlw
itab_wd(iw3)%mrlw_orig = itab_wd(iw1)%mrlw_orig
         
! Fill loop indices for new U points

itab_ud(iu6)%iup = iu6
call udloops('f',iu6, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu6,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu8)%iup = iu8
call udloops('f',iu8, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu8,21,22,23, 0, 0, 0, 0, 0, 0, 0)

return

5 continue

write(io6,*) 'In subroutine perim_fill_concave2, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_concave2

!===========================================================================

subroutine perim_fill_concave4(ngr,mrloo,kma,kua,kwa,  &
                               jm1,ju1,jm2,ju2,jm3,ju3,jm4,ju4,jm5)

use mem_grid, only:  xem, yem, zem
use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa   ! Index values for latest added points

integer, intent(in) :: jm1,jm2,jm3,jm4,jm5   ! Original M indices

integer, intent(in) :: ju1,ju2,ju3,ju4       ! Original U indices

! Original J indices to be determined

integer :: jm6,jm7,jm8,jm9

integer :: ju5,ju6,ju7,ju8,ju9,ju10
integer :: ju11,ju12,ju13,ju14,ju15,ju16

integer :: jw1,jw2,jw3,jw4,jw5,jw6,jw7,jw8

! I indices

integer :: iw1,iw2,iw3,iw4,iw5,iw6,iw7,iw8,iw9,iw10
integer :: iw11,iw12,iw13,iw14,iw15,iw16

integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10
integer :: iu11,iu12,iu13,iu14,iu15,iu16,iu17,iu18,iu19,iu20
integer :: iu21,iu22,iu23,iu24,iu25,iu26,iu27,iu28,iu29,iu30

integer :: im1,im2,im3,im4,im5,im6,im7,im8,im9,im10
integer :: im11,im12,im13,im14,im15

! Newly added vertices

im11 = kma + 1
im13 = kma + 2

! Newly added edges

iu10 = kua + 1
iu12 = kua + 2
iu14 = kua + 3
iu16 = kua + 4
iu18 = kua + 5
iu20 = kua + 6
iu23 = kua + 7
iu24 = kua + 8
iu28 = kua + 9
iu29 = kua + 10

! Newly added areas

iw3  = kwa + 1
iw4  = kwa + 2
iw6  = kwa + 3
iw7  = kwa + 4
iw9  = kwa + 5
iw10 = kwa + 6
iw13 = kwa + 7
iw15 = kwa + 8

! Increment indices for newly added points

kma = kma + 2
kua = kua + 10
kwa = kwa + 8

! Determine orientation of positive ju1 direction to get jw1

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points into fine mesh
   jw1 = ltab_ud(ju1)%iw(1)
   jw2 = ltab_ud(ju1)%iw(4)
   ju5 = ltab_ud(ju1)%iu(1)
   ju6 = ltab_ud(ju1)%iu(2)
   ju7 = ltab_ud(ju1)%iu(8)
   ju12 = ltab_ud(ju1)%iu(7)

   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

   itab_ud(iu1)%iw(1) = jw1  ! not needed?
   itab_ud(iu2)%iw(1) = iw3
else
   jw1 = ltab_ud(ju1)%iw(2)
   jw2 = ltab_ud(ju1)%iw(5)
   ju5 = ltab_ud(ju1)%iu(4)
   ju6 = ltab_ud(ju1)%iu(3)
   ju7 = ltab_ud(ju1)%iu(9)
   ju12 = ltab_ud(ju1)%iu(10)

   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = jw1
   itab_ud(iu2)%iw(2) = iw3
endif

! Determine orientation of positive ju2 direction to get jw3

if (jm2 == ltab_ud(ju2)%im(1)) then  ! Positive ju2 points inward
   jw3 = ltab_ud(ju2)%iw(1)
   ju8 = ltab_ud(ju2)%iu(2)

   iu3 = ju2
   iu4 = nest_ud(ju2)%iu

   itab_ud(iu3)%iw(1) = jw3
   itab_ud(iu4)%iw(1) = iw6
else                              ! Positive ju2 points outward
   jw3 = ltab_ud(ju2)%iw(2)
   ju8 = ltab_ud(ju2)%iu(3)

   iu3 = nest_ud(ju2)%iu
   iu4 = ju2

   itab_ud(iu3)%iw(2) = jw3  ! not needed?
   itab_ud(iu4)%iw(2) = iw6
endif

! Determine orientation of positive ju3 direction to get jw4

if (jm3 == ltab_ud(ju3)%im(1)) then  ! Positive ju3 points inward
   jw4 = ltab_ud(ju3)%iw(1)

   iu5 = ju3
   iu6 = nest_ud(ju3)%iu

   itab_ud(iu5)%iw(1) = iw6
   itab_ud(iu6)%iw(1) = jw4
else                              ! Positive ju3 points outward
   jw4 = ltab_ud(ju3)%iw(2)

   iu5 = nest_ud(ju3)%iu
   iu6 = ju3

   itab_ud(iu5)%iw(2) = iw6  ! not needed?
   itab_ud(iu6)%iw(2) = jw4
endif

! Determine orientation of positive ju4 direction to get jw6

if (jm4 == ltab_ud(ju4)%im(1)) then  ! Positive ju4 points inward
   jw6 = ltab_ud(ju4)%iw(1)
   jw5 = ltab_ud(ju4)%iw(3)
   ju10 = ltab_ud(ju4)%iu(1)
   ju11 = ltab_ud(ju4)%iu(2)
   ju9 = ltab_ud(ju4)%iu(5)
   ju13 = ltab_ud(ju4)%iu(6)

   iu7 = ju4
   iu8 = nest_ud(ju4)%iu

   itab_ud(iu7)%iw(1) = iw10
   itab_ud(iu8)%iw(1) = jw6
else                              ! Positive ju4 points outward
   jw6 = ltab_ud(ju4)%iw(2)
   jw5 = ltab_ud(ju4)%iw(6)
   ju10 = ltab_ud(ju4)%iu(4)
   ju11 = ltab_ud(ju4)%iu(3)
   ju9 = ltab_ud(ju4)%iu(12)
   ju13 = ltab_ud(ju4)%iu(11)

   iu7 = nest_ud(ju4)%iu
   iu8 = ju4

   itab_ud(iu7)%iw(2) = iw10  ! not needed?
   itab_ud(iu8)%iw(2) = jw6
endif

! Determine neighbors of ju12

if (jw2 == ltab_ud(ju12)%iw(2)) then
   jm6 = ltab_ud(ju12)%im(1)
   jm7 = ltab_ud(ju12)%im(2)
   jw7 = ltab_ud(ju12)%iw(1)
   jw8 = ltab_ud(ju12)%iw(4)

   ju14 = ltab_ud(ju12)%iu(1)
   ju15 = ltab_ud(ju12)%iu(2)
   ju16 = ltab_ud(ju12)%iu(7)
else
   jm6 = ltab_ud(ju12)%im(2)
   jm7 = ltab_ud(ju12)%im(1)
   jw7 = ltab_ud(ju12)%iw(2)
   jw8 = ltab_ud(ju12)%iw(5)

   ju14 = ltab_ud(ju12)%iu(4)
   ju15 = ltab_ud(ju12)%iu(3)
   ju16 = ltab_ud(ju12)%iu(10)
endif

! Determine neighbors of ju16

if (jw8 == ltab_ud(ju16)%iw(1)) then
   jm8 = ltab_ud(ju16)%im(1)
   jm9 = ltab_ud(ju16)%im(2)
else
   jm8 = ltab_ud(ju16)%im(2)
   jm9 = ltab_ud(ju16)%iw(1)
endif

! Assign remaining I indices that already existed as J indices

! (iu1:iu8 already assigned above)

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2
im4 = nest_ud(ju2)%im
im5 = jm3
im6 = nest_ud(ju3)%im
im7 = jm4
im8 = nest_ud(ju4)%im
im9 = jm5

im10 = jm6
im12 = jm7
im14 = jm8
im15 = jm9

iu9 = ju5
iu11 = ju6
iu13 = ju7
iu15 = ju8
iu17 = ju9
iu19 = ju10
iu21 = ju11
iu22 = ju12
iu25 = ju13
iu26 = ju15
iu27 = ju14
iu30 = ju16

iw1 = jw1
iw2 = jw2
iw5 = jw3
iw8 = jw4
iw11 = jw5
iw12 = jw6
iw14 = jw7
iw16 = jw8

! Fill remaining itab_u and itab_w members

! Nothing needed for iu9

itab_ud(iu10)%im(1) = im2
itab_ud(iu10)%im(2) = im10
itab_ud(iu10)%iw(1) = iw1
itab_ud(iu10)%iw(2) = iw2

itab_ud(iu11)%im(1) = im2
itab_ud(iu11)%im(2) = im11
itab_ud(iu11)%iw(1) = iw2
itab_ud(iu11)%iw(2) = iw3

itab_ud(iu12)%im(1) = im3
itab_ud(iu12)%im(2) = im11
itab_ud(iu12)%iw(1) = iw3
itab_ud(iu12)%iw(2) = iw4

itab_ud(iu13)%im(1) = im3
itab_ud(iu13)%im(2) = im12
itab_ud(iu13)%iw(1) = iw4
itab_ud(iu13)%iw(2) = iw5

itab_ud(iu14)%im(1) = im4
itab_ud(iu14)%im(2) = im12
itab_ud(iu14)%iw(1) = iw5
itab_ud(iu14)%iw(2) = iw7

itab_ud(iu15)%im(1) = im4
itab_ud(iu15)%im(2) = im6
itab_ud(iu15)%iw(1) = iw7
itab_ud(iu15)%iw(2) = iw6

itab_ud(iu16)%im(1) = im6
itab_ud(iu16)%im(2) = im12
itab_ud(iu16)%iw(1) = iw7
itab_ud(iu16)%iw(2) = iw8

itab_ud(iu17)%im(1) = im7
itab_ud(iu17)%im(2) = im12
itab_ud(iu17)%iw(1) = iw8
itab_ud(iu17)%iw(2) = iw9

itab_ud(iu18)%im(1) = im7
itab_ud(iu18)%im(2) = im13
itab_ud(iu18)%iw(1) = iw9
itab_ud(iu18)%iw(2) = iw10

itab_ud(iu19)%im(1) = im8
itab_ud(iu19)%im(2) = im13
itab_ud(iu19)%iw(1) = iw10
itab_ud(iu19)%iw(2) = iw11

itab_ud(iu20)%im(1) = im8
itab_ud(iu20)%im(2) = im14
itab_ud(iu20)%iw(1) = iw11
itab_ud(iu20)%iw(2) = iw12

! Nothing needed for iu21

itab_ud(iu22)%im(1) = im10
itab_ud(iu22)%im(2) = im11
itab_ud(iu22)%iw(1) = iw14
itab_ud(iu22)%iw(2) = iw2

itab_ud(iu23)%im(1) = im11
itab_ud(iu23)%im(2) = im12
itab_ud(iu23)%iw(1) = iw13
itab_ud(iu23)%iw(2) = iw4

itab_ud(iu24)%im(1) = im12
itab_ud(iu24)%im(2) = im13
itab_ud(iu24)%iw(1) = iw13
itab_ud(iu24)%iw(2) = iw9

itab_ud(iu25)%im(1) = im13
itab_ud(iu25)%im(2) = im14
itab_ud(iu25)%iw(1) = iw16
itab_ud(iu25)%iw(2) = iw11

itab_ud(iu26)%im(1) = im11
itab_ud(iu26)%im(2) = im13
itab_ud(iu26)%iw(1) = iw15
itab_ud(iu26)%iw(2) = iw13

! Nothing needed for iu27

itab_ud(iu28)%im(1) = im11
itab_ud(iu28)%im(2) = im15
itab_ud(iu28)%iw(1) = iw14
itab_ud(iu28)%iw(2) = iw15

itab_ud(iu29)%im(1) = im13
itab_ud(iu29)%im(2) = im15
itab_ud(iu29)%iw(1) = iw15
itab_ud(iu29)%iw(2) = iw16

! Nothing needed for iu30

itab_wd(iw1)%iu(1) = iu1
itab_wd(iw1)%iu(2) = iu9
itab_wd(iw1)%iu(3) = iu10

itab_wd(iw2)%iu(1) = iu10
itab_wd(iw2)%iu(2) = iu22
itab_wd(iw2)%iu(3) = iu11

itab_wd(iw3)%iu(1) = iu2
itab_wd(iw3)%iu(2) = iu11
itab_wd(iw3)%iu(3) = iu12

itab_wd(iw4)%iu(1) = iu12
itab_wd(iw4)%iu(2) = iu23
itab_wd(iw4)%iu(3) = iu13

itab_wd(iw5)%iu(1) = iu3
itab_wd(iw5)%iu(2) = iu13
itab_wd(iw5)%iu(3) = iu14

itab_wd(iw6)%iu(1) = iu4
itab_wd(iw6)%iu(2) = iu15
itab_wd(iw6)%iu(3) = iu5

itab_wd(iw7)%iu(1) = iu15
itab_wd(iw7)%iu(2) = iu14
itab_wd(iw7)%iu(3) = iu16

itab_wd(iw8)%iu(1) = iu6
itab_wd(iw8)%iu(2) = iu16
itab_wd(iw8)%iu(3) = iu17

itab_wd(iw9)%iu(1) = iu17
itab_wd(iw9)%iu(2) = iu24
itab_wd(iw9)%iu(3) = iu18

itab_wd(iw10)%iu(1) = iu7
itab_wd(iw10)%iu(2) = iu18
itab_wd(iw10)%iu(3) = iu19

itab_wd(iw11)%iu(1) = iu19
itab_wd(iw11)%iu(2) = iu25
itab_wd(iw11)%iu(3) = iu20

itab_wd(iw12)%iu(1) = iu8
itab_wd(iw12)%iu(2) = iu20
itab_wd(iw12)%iu(3) = iu21

itab_wd(iw13)%iu(1) = iu23
itab_wd(iw13)%iu(2) = iu26
itab_wd(iw13)%iu(3) = iu24

itab_wd(iw14)%iu(1) = iu22
itab_wd(iw14)%iu(2) = iu27
itab_wd(iw14)%iu(3) = iu28

itab_wd(iw15)%iu(1) = iu26
itab_wd(iw15)%iu(2) = iu28
itab_wd(iw15)%iu(3) = iu29

itab_wd(iw16)%iu(1) = iu29
itab_wd(iw16)%iu(2) = iu30
itab_wd(iw16)%iu(3) = iu25

! Fill earth coordinates for new M points

xem(im11) = .33 * xem(im6) + .67 * xem(im10)
yem(im11) = .33 * yem(im6) + .67 * yem(im10)
zem(im11) = .33 * zem(im6) + .67 * zem(im10)

xem(im13) = .33 * xem(im4) + .67 * xem(im14)
yem(im13) = .33 * yem(im4) + .67 * yem(im14)
zem(im13) = .33 * zem(im4) + .67 * zem(im14)

! Shift earth coordinates for IM7 point

xem(im12) = .33 * xem(im5) + .67 * xem(im12)
yem(im12) = .33 * yem(im5) + .67 * yem(im12)
zem(im12) = .33 * zem(im5) + .67 * zem(im12)

! Check mrlw for old points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5
if (itab_wd(iw2)%mrlw /= mrloo) go to 5
if (itab_wd(iw5)%mrlw /= mrloo) go to 5
if (itab_wd(iw8)%mrlw /= mrloo) go to 5
if (itab_wd(iw11)%mrlw /= mrloo) go to 5
if (itab_wd(iw12)%mrlw /= mrloo) go to 5
if (itab_wd(iw14)%mrlw /= mrloo) go to 5
if (itab_wd(iw16)%mrlw /= mrloo) go to 5

! Fill mrlw for new W points

itab_wd(iw3)%mrlw      = itab_wd(jw1)%mrlw
itab_wd(iw3)%mrlw_orig = itab_wd(jw1)%mrlw_orig
         
itab_wd(iw4)%mrlw      = itab_wd(jw2)%mrlw
itab_wd(iw4)%mrlw_orig = itab_wd(jw2)%mrlw_orig
         
itab_wd(iw6)%mrlw      = itab_wd(jw3)%mrlw
itab_wd(iw6)%mrlw_orig = itab_wd(jw3)%mrlw_orig
         
itab_wd(iw7)%mrlw      = itab_wd(jw3)%mrlw
itab_wd(iw7)%mrlw_orig = itab_wd(jw3)%mrlw_orig
         
itab_wd(iw9)%mrlw      = itab_wd(jw5)%mrlw
itab_wd(iw9)%mrlw_orig = itab_wd(jw5)%mrlw_orig
         
itab_wd(iw10)%mrlw      = itab_wd(jw6)%mrlw
itab_wd(iw10)%mrlw_orig = itab_wd(jw6)%mrlw_orig

itab_wd(iw13)%mrlw      = itab_wd(jw7)%mrlw
itab_wd(iw13)%mrlw_orig = itab_wd(jw7)%mrlw_orig
         
itab_wd(iw15)%mrlw      = itab_wd(jw7)%mrlw
itab_wd(iw15)%mrlw_orig = itab_wd(jw7)%mrlw_orig
         
! Fill loop indices for new U points

itab_ud(iu10)%iup = iu10
call udloops('f',iu10, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu10,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu12)%iup = iu12
call udloops('f',iu12, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu12,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu14)%iup = iu14
call udloops('f',iu14, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu14,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu16)%iup = iu16
call udloops('f',iu16, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu16,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu18)%iup = iu18
call udloops('f',iu18, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu18,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu20)%iup = iu20
call udloops('f',iu20, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu20,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu23)%iup = iu23
call udloops('f',iu23, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu23,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu24)%iup = iu24
call udloops('f',iu24, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu24,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu28)%iup = iu28
call udloops('f',iu28, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu28,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu29)%iup = iu29
call udloops('f',iu29, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu29,21,22,23, 0, 0, 0, 0, 0, 0, 0)

return

5 continue

write(io6,*) 'In subroutine perim_fill_concave4, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_concave4

!===========================================================================

subroutine perim_fill_convex2(ngr,mrloo,kma,kua,kwa,  &
                              jm1,ju1,jm2,ju2,jm3)

use mem_grid, only:  xem, yem, zem
use mem_ijtabs, only: itab_wd, itab_ud, ltab_wd, ltab_ud, nest_wd, nest_ud
use misc_coms,  only: io6

implicit none

integer, intent(in) :: ngr,mrloo
integer, intent(inout) :: kma,kua,kwa   ! Index values for latest added points

integer, intent(in) :: jm1,jm2,jm3   ! Original M indices

integer, intent(in) :: ju1,ju2       ! Original U indices

! Original J indices to be determined

integer :: jm4,jm5,jm6,jm7,jm8,jm9

integer :: ju3,ju4,ju5,ju6,ju7,ju8,ju9,ju10
integer :: ju11,ju12,ju13,ju14,ju15,ju16

integer :: jw1,jw2,jw3,jw4,jw5,jw6,jw7,jw8

! I indices

integer :: iw1,iw2,iw3,iw4,iw5,iw6,iw7,iw8,iw9,iw10
integer :: iw11,iw12

integer :: iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8,iu9,iu10
integer :: iu11,iu12,iu13,iu14,iu15,iu16,iu17,iu18,iu19,iu20
integer :: iu21,iu22,iu23

integer :: im1,im2,im3,im4,im5,im6,im7,im8,im9,im10
integer :: im11,im12

! Newly added vertices

im8 = kma + 1

! Newly added edges

iu6  = kua + 1
iu9  = kua + 2
iu11 = kua + 3
iu14 = kua + 4
iu19 = kua + 5

! Newly added areas

iw3  = kwa + 1
iw4  = kwa + 2
iw5  = kwa + 3
iw10 = kwa + 4

! Increment indices for newly added points

kma = kma + 1
kua = kua + 5
kwa = kwa + 4

! Determine orientation of positive ju1 direction to get jw1

if (jm1 == ltab_ud(ju1)%im(1)) then  ! Positive ju1 points into fine mesh
   jw1 = ltab_ud(ju1)%iw(1)
   jw2 = ltab_ud(ju1)%iw(4)
   ju3 = ltab_ud(ju1)%iu(1)
   ju4 = ltab_ud(ju1)%iu(2)
   ju5 = ltab_ud(ju1)%iu(8)
   ju8 = ltab_ud(ju1)%iu(7)

   iu1 = ju1
   iu2 = nest_ud(ju1)%iu

   itab_ud(iu1)%iw(1) = jw1  ! not needed?
   itab_ud(iu2)%iw(1) = iw3
else
   jw1 = ltab_ud(ju1)%iw(2)
   jw2 = ltab_ud(ju1)%iw(5)
   ju3 = ltab_ud(ju1)%iu(4)
   ju4 = ltab_ud(ju1)%iu(3)
   ju5 = ltab_ud(ju1)%iu(9)
   ju8 = ltab_ud(ju1)%iu(10)

   iu1 = nest_ud(ju1)%iu
   iu2 = ju1

   itab_ud(iu1)%iw(2) = jw1
   itab_ud(iu2)%iw(2) = iw3
endif

! Determine orientation of positive ju2 direction to get jw4

if (jm2 == ltab_ud(ju2)%im(1)) then  ! Positive ju2 points inward
   jw4 = ltab_ud(ju2)%iw(1)
   jw3 = ltab_ud(ju2)%iw(3)
   ju6 = ltab_ud(ju2)%iu(1)
   ju7 = ltab_ud(ju2)%iu(2)
   ju9 = ltab_ud(ju2)%iu(6)

   iu3 = ju2
   iu4 = nest_ud(ju2)%iu

   itab_ud(iu3)%iw(1) = iw5
   itab_ud(iu4)%iw(1) = jw4
else                              ! Positive ju2 points outward
   jw4 = ltab_ud(ju2)%iw(2)
   jw3 = ltab_ud(ju2)%iw(6)
   ju6 = ltab_ud(ju2)%iu(4)
   ju7 = ltab_ud(ju2)%iu(3)
   ju9 = ltab_ud(ju2)%iu(11)

   iu3 = nest_ud(ju2)%iu
   iu4 = ju2

   itab_ud(iu3)%iw(2) = iw5  ! not needed?
   itab_ud(iu4)%iw(2) = jw4
endif

! Determine neighbors of ju8

if (jw2 == ltab_ud(ju8)%iw(2)) then
   jm4 = ltab_ud(ju8)%im(1)
   jm5 = ltab_ud(ju8)%im(2)

   jw5 = ltab_ud(ju8)%iw(1)
   jw6 = ltab_ud(ju8)%iw(4)

   ju10 = ltab_ud(ju8)%iu(1)
   ju11 = ltab_ud(ju8)%iu(2)
   ju12 = ltab_ud(ju8)%iu(8)
   ju15 = ltab_ud(ju8)%iu(7)
else
   jm4 = ltab_ud(ju8)%im(2)
   jm5 = ltab_ud(ju8)%im(1)

   jw5 = ltab_ud(ju8)%iw(2)
   jw6 = ltab_ud(ju8)%iw(5)

   ju10 = ltab_ud(ju8)%iu(4)
   ju11 = ltab_ud(ju8)%iu(3)
   ju12 = ltab_ud(ju8)%iu(9)
   ju15 = ltab_ud(ju8)%iu(10)
endif

! Determine neighbors of ju9

if (jw3 == ltab_ud(ju9)%iw(1)) then
   jm6 = ltab_ud(ju9)%im(1)

   jw8 = ltab_ud(ju9)%iw(2)
   jw7 = ltab_ud(ju9)%iw(6)

   ju13 = ltab_ud(ju9)%iu(4)
   ju14 = ltab_ud(ju9)%iu(3)
   ju16 = ltab_ud(ju9)%iu(11)
else
   jm6 = ltab_ud(ju9)%im(2)

   jw8 = ltab_ud(ju9)%iw(1)
   jw7 = ltab_ud(ju9)%iw(3)

   ju13 = ltab_ud(ju9)%iu(1)
   ju14 = ltab_ud(ju9)%iu(2)
   ju16 = ltab_ud(ju9)%iu(6)
endif

! Determine neighbors of ju15

if (jw6 == ltab_ud(ju15)%iw(1)) then
   jm7 = ltab_ud(ju15)%im(2)
   jm8 = ltab_ud(ju15)%im(1)
else
   jm7 = ltab_ud(ju15)%im(1)
   jm8 = ltab_ud(ju15)%im(2)
endif

! Determine neighbors of ju16

if (jw7 == ltab_ud(ju16)%iw(1)) then
   jm9 = ltab_ud(ju16)%im(1)
else
   jm9 = ltab_ud(ju16)%im(2)
endif

! Assign remaining I indices that already existed as J indices

! (iu1:iu4 already assigned above)

im1 = jm1
im2 = nest_ud(ju1)%im
im3 = jm2
im4 = nest_ud(ju2)%im
im5 = jm3
im6 = jm4
im7 = jm5
im9 = jm6

im10 = jm7
im11 = jm8
im12 = jm9

iu5 = ju3
iu7 = ju4
iu8 = ju5
iu10 = ju6
iu12 = ju7
iu13 = ju8
iu15 = ju9
iu16 = ju10
iu17 = ju11
iu18 = ju12
iu20 = ju13
iu21 = ju14
iu22 = ju15
iu23 = ju16

iw1 = jw1
iw2 = jw2
iw6 = jw3
iw7 = jw4
iw8 = jw5
iw9 = jw6
iw11 = jw7
iw12 = jw8

! Fill remaining itab_u and itab_w members

! Nothing needed for iu5

itab_ud(iu6)%im(1) = im2
itab_ud(iu6)%im(2) = im6
itab_ud(iu6)%iw(1) = iw1
itab_ud(iu6)%iw(2) = iw2

itab_ud(iu7)%im(1) = im2
itab_ud(iu7)%im(2) = im7
itab_ud(iu7)%iw(1) = iw2
itab_ud(iu7)%iw(2) = iw3

itab_ud(iu8)%im(1) = im3
itab_ud(iu8)%im(2) = im7
itab_ud(iu8)%iw(1) = iw3
itab_ud(iu8)%iw(2) = iw4

itab_ud(iu9)%im(1) = im3
itab_ud(iu9)%im(2) = im8
itab_ud(iu9)%iw(1) = iw4
itab_ud(iu9)%iw(2) = iw5

itab_ud(iu10)%im(1) = im4
itab_ud(iu10)%im(2) = im8
itab_ud(iu10)%iw(1) = iw5
itab_ud(iu10)%iw(2) = iw6

itab_ud(iu11)%im(1) = im4
itab_ud(iu11)%im(2) = im9
itab_ud(iu11)%iw(1) = iw6
itab_ud(iu11)%iw(2) = iw7

! Nothing needed for iu12

itab_ud(iu13)%im(1) = im6
itab_ud(iu13)%im(2) = im7
itab_ud(iu13)%iw(1) = iw8
itab_ud(iu13)%iw(2) = iw2

itab_ud(iu14)%im(1) = im7
itab_ud(iu14)%im(2) = im8
itab_ud(iu14)%iw(1) = iw10
itab_ud(iu14)%iw(2) = iw4

itab_ud(iu15)%im(1) = im8
itab_ud(iu15)%im(2) = im9
itab_ud(iu15)%iw(1) = iw12
itab_ud(iu15)%iw(2) = iw6

! Nothing needed for iu16

itab_ud(iu17)%im(1) = im7
itab_ud(iu17)%im(2) = im10
itab_ud(iu17)%iw(1) = iw8
itab_ud(iu17)%iw(2) = iw9

itab_ud(iu18)%im(1) = im7
itab_ud(iu18)%im(2) = im11
itab_ud(iu18)%iw(1) = iw9
itab_ud(iu18)%iw(2) = iw10

itab_ud(iu19)%im(1) = im8
itab_ud(iu19)%im(2) = im11
itab_ud(iu19)%iw(1) = iw10
itab_ud(iu19)%iw(2) = iw11

itab_ud(iu20)%im(1) = im8
itab_ud(iu20)%im(2) = im12
itab_ud(iu20)%iw(1) = iw11
itab_ud(iu20)%iw(2) = iw12

! Nothing needed for iu21
! Nothing needed for iu22
! Nothing needed for iu23

itab_wd(iw1)%iu(1) = iu1
itab_wd(iw1)%iu(2) = iu5
itab_wd(iw1)%iu(3) = iu6

itab_wd(iw2)%iu(1) = iu6
itab_wd(iw2)%iu(2) = iu13
itab_wd(iw2)%iu(3) = iu7

itab_wd(iw3)%iu(1) = iu2
itab_wd(iw3)%iu(2) = iu7
itab_wd(iw3)%iu(3) = iu8

itab_wd(iw4)%iu(1) = iu8
itab_wd(iw4)%iu(2) = iu14
itab_wd(iw4)%iu(3) = iu9

itab_wd(iw5)%iu(1) = iu3
itab_wd(iw5)%iu(2) = iu9
itab_wd(iw5)%iu(3) = iu10

itab_wd(iw6)%iu(1) = iu10
itab_wd(iw6)%iu(2) = iu15
itab_wd(iw6)%iu(3) = iu11

itab_wd(iw7)%iu(1) = iu4
itab_wd(iw7)%iu(2) = iu11
itab_wd(iw7)%iu(3) = iu12

itab_wd(iw8)%iu(1) = iu13
itab_wd(iw8)%iu(2) = iu16
itab_wd(iw8)%iu(3) = iu17

itab_wd(iw9)%iu(1) = iu17
itab_wd(iw9)%iu(2) = iu22
itab_wd(iw9)%iu(3) = iu18

itab_wd(iw10)%iu(1) = iu14
itab_wd(iw10)%iu(2) = iu18
itab_wd(iw10)%iu(3) = iu19

itab_wd(iw11)%iu(1) = iu19
itab_wd(iw11)%iu(2) = iu23
itab_wd(iw11)%iu(3) = iu20

itab_wd(iw12)%iu(1) = iu15
itab_wd(iw12)%iu(2) = iu20
itab_wd(iw12)%iu(3) = iu21

! Fill earth coordinates for new M point

xem(im8) = .33 * xem(im9) + .67 * xem(im7)
yem(im8) = .33 * yem(im9) + .67 * yem(im7)
zem(im8) = .33 * zem(im9) + .67 * zem(im7)

! Shift earth coordinates for IM7 point

xem(im7) = .33 * xem(im6) + .67 * xem(im7)
yem(im7) = .33 * yem(im6) + .67 * yem(im7)
zem(im7) = .33 * zem(im6) + .67 * zem(im7)

! Check mrlw for old points

if (itab_wd(iw1)%mrlw /= mrloo) go to 5
if (itab_wd(iw2)%mrlw /= mrloo) go to 5
if (itab_wd(iw6)%mrlw /= mrloo) go to 5
if (itab_wd(iw7)%mrlw /= mrloo) go to 5
if (itab_wd(iw8)%mrlw /= mrloo) go to 5
if (itab_wd(iw9)%mrlw /= mrloo) go to 5
if (itab_wd(iw11)%mrlw /= mrloo) go to 5
if (itab_wd(iw12)%mrlw /= mrloo) go to 5

! Fill mrlw for new W points

itab_wd(iw3)%mrlw      = itab_wd(jw1)%mrlw
itab_wd(iw3)%mrlw_orig = itab_wd(jw1)%mrlw_orig
         
itab_wd(iw4)%mrlw      = itab_wd(jw2)%mrlw
itab_wd(iw4)%mrlw_orig = itab_wd(jw2)%mrlw_orig
         
itab_wd(iw5)%mrlw      = itab_wd(jw4)%mrlw
itab_wd(iw5)%mrlw_orig = itab_wd(jw4)%mrlw_orig
         
itab_wd(iw10)%mrlw      = itab_wd(jw6)%mrlw
itab_wd(iw10)%mrlw_orig = itab_wd(jw6)%mrlw_orig

! Fill loop indices for new U points

itab_ud(iu6)%iup = iu6
call udloops('f',iu6, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu6,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu9)%iup = iu9
call udloops('f',iu9, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu9,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu11)%iup = iu11
call udloops('f',iu11, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu11,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu14)%iup = iu14
call udloops('f',iu14, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu14,21,22,23, 0, 0, 0, 0, 0, 0, 0)

itab_ud(iu19)%iup = iu19
call udloops('f',iu19, 1, 4, 7, 8,11,12,13,14,16,20)
call udloops('n',iu19,21,22,23, 0, 0, 0, 0, 0, 0, 0)

!print*, ' '

!write(io6,'(a,30i6)') 'jms ',jm1,jm2,jm3,jm4,jm5,jm6,jm7,jm8,jm9

!write(io6,'(a,30i6)') 'jus ',ju1,ju2,ju3,ju4,ju5,ju6,ju7,ju8, &
!                             ju9,ju10,ju11,ju12,ju13,ju14,ju15,ju16

!write(io6,'(a,30i6)') 'g1 ',im7,im10,im8,im11

!write(io6,'(a,30i6)') 'g2 ',itab_ud(iu17)%im(1), &
!                            itab_ud(iu17)%im(2), &
!                            itab_ud(iu19)%im(1), &
!                            itab_ud(iu19)%im(2)

!write(io6,'(a,12f10.0)') 'g3 ',xem(itab_ud(iu17)%im(1)), &
!                               xem(itab_ud(iu17)%im(2)), &
!                               xem(itab_ud(iu19)%im(1)), &
!                               xem(itab_ud(iu19)%im(2)), &
!                               yem(itab_ud(iu17)%im(1)), &
!                               yem(itab_ud(iu17)%im(2)), &
!                               yem(itab_ud(iu19)%im(1)), &
!                               yem(itab_ud(iu19)%im(2)), &
!                               zem(itab_ud(iu17)%im(1)), &
!                               zem(itab_ud(iu17)%im(2)), &
!                               zem(itab_ud(iu19)%im(1)), &
!                               zem(itab_ud(iu19)%im(2))
return

5 continue

write(io6,*) 'In subroutine perim_fill_convex2, current nested grid ',ngr
write(io6,*) 'crosses pre-existing grid boundary.'
write(io6,*) 'iw1 = ',iw1
write(io6,*) 'stopping model'
stop 'stop - Nested grid out of bounds'

end subroutine perim_fill_convex2

