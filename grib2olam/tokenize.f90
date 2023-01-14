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

module token_mod

contains

subroutine tokenize_ws(str,tokens,ntok)

  implicit none

  integer          :: ntok
  character(len=*) :: str, tokens(:)
  character(len=1) :: sep,tab

  integer :: ntokmax,nch,nc,ntbeg,ntend,n,nn

! This routine "parses" character string str into different pieces
! or tokens by looking for possible token separators. The number of
! tokens identified is returned as ntok, and each token is returned
! in the array of strings tokens.

! Modified 12/11/07: this version tokenizes by white space (spaces and tabs)
  sep = char(32)
  tab = char(9)

  ntok=0
  nch=len_trim(str)
  nc=1

  ntokmax = size(tokens,dim=1)

  do nn=1,ntokmax
     ntok=ntok+1
     do n=nc,nch
        if (str(n:n) /= sep .and. str(n:n) /= tab) then
           ntbeg=n
           exit
        endif
     enddo

     do n=ntbeg,nch
        if (str(n:n) == sep .or. str(n:n) == tab) then
           ntend=n-1
           exit
        endif
        if (n == nch) then
           ntend=n
           exit
        endif
     enddo

     tokens(ntok)=str(ntbeg:ntend)
     nc=ntend+1
     if(nc.ge.nch) exit
  enddo

  tokens(ntok+1:ntokmax) = ' '

end subroutine tokenize_ws



subroutine tokenize(str,tokens,ntok,sep)

  implicit none

  integer          :: ntok
  character(len=*) :: str, tokens(:)
  character(len=1) :: sep

  integer :: ntokmax,nch,nc,ntbeg,ntend,n,nn

! This routine "parses" character string str into different pieces
! or tokens by looking for the token separator sep. The number of
! tokens identified is returned as ntok, and each token is returned
! in the array of strings tokens.

  ntok=0
  nch=len_trim(str)
  nc=1

  ntokmax = size(tokens,dim=1)

  do nn=1,ntokmax
     ntok=ntok+1
     do n=nc,nch
        if (str(n:n) /= sep) then
           ntbeg=n
           exit
        endif
     enddo

     do n=ntbeg,nch
        if (str(n:n) == sep) then
           ntend=n-1
           exit
        endif
        if (n == nch) then
           ntend=n
           exit
        endif
     enddo

     tokens(ntok)=str(ntbeg:ntend)
     nc=ntend+1
     if(nc.ge.nch) exit
  enddo

  tokens(ntok+1:ntokmax) = ' '

end subroutine tokenize

end module token_mod
