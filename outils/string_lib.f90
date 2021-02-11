!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

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

!===============================================================================

module string_lib

contains

function to_upper(string) result(upper)
  
  implicit none

  character(len=*), intent(in) :: string
  character(len=len(string))   :: upper

  integer :: j

  do j = 1, len(string)
     if(string(j:j) >= "a" .and. string(j:j) <= "z") then
        upper(j:j) = achar(iachar(string(j:j)) - 32)
     else
        upper(j:j) = string(j:j)
     endif
  enddo

end function to_upper

!===============================================================================

function to_lower(string) result(lower)
  
  implicit none

  character(len=*), intent(in) :: string
  character(len=len(string))   :: lower

  integer :: j

  do j = 1, len(string)
     if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
        lower(j:j) = achar(iachar(string(j:j)) + 32)
     else
        lower(j:j) = string(j:j)
     endif
  enddo

end function to_lower

!===============================================================================

end module string_lib
