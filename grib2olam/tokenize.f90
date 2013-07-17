!
! Copyright (C) 1991-2003  ; All Rights Reserved ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================

subroutine tokenize_ws(str,tokens,ntok)

implicit none

integer :: ntok
character(len=*) :: str,tokens(*)
character(len=1) :: sep,tab

integer :: ntokmax=100,npt,nch,nc,ntbeg,ntend,n

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

! Modified 12/11/07: this version tokenizes by white space (spaces and tabs)
sep = char(32)
tab = char(9)

ntok=0
npt=1
nch=len_trim(str)
nc=1
do ntok=1,ntokmax
   do n=nc,nch
      if(str(n:n) /= sep .and. str(n:n) /= tab) then
         ntbeg=n
         goto 21
      endif
   enddo
   21 continue
   do n=ntbeg,nch
      if(str(n:n) == sep .or. str(n:n) == tab) then
         ntend=n-1
         goto 22
      endif
      if(n == nch) then
         ntend=n
         goto 22
      endif
   enddo
   22 continue
   tokens(ntok)=str(ntbeg:ntend)
   nc=ntend+1
   if(nc.ge.nch) goto 25
enddo

25 continue

return
end



subroutine tokenize(str,tokens,ntok,sep)
character*(*) str,tokens(*)
character sep*1
data ntokmax/100/

! this routine "parses" character string str into different pieces
! or tokens by looking for  possible token separators (toks
! str contains nch characters.  the number of tokens identified is nto
! the character string tokens are stored in tokens.

!sep=' '
ntok=0
npt=1
nch=len_trim(str)
nc=1
do ntok=1,ntokmax
   do n=nc,nch
      if(str(n:n).ne.sep) then
         ntbeg=n
         goto 21
      endif
   enddo
   21 continue
   do n=ntbeg,nch
      if(str(n:n).eq.sep) then
         ntend=n-1
         goto 22
      endif
      if(n.eq.nch) then
         ntend=n
         goto 22
      endif
   enddo
   22 continue
   tokens(ntok)=str(ntbeg:ntend)
   nc=ntend+1
   if(nc.ge.nch) goto 25
enddo

25 continue

return
end
