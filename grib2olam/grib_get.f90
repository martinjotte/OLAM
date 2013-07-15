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

module grib_get_mod

  integer :: nx,ny,nrec,iscan,longdate
  integer, allocatable :: levs(:),idates(:),ifhrs(:)
  character(len=20), allocatable :: fields(:),irecs(:)
  character(len=20) :: projection,cpole
  real :: alat1,alon1,alov,aorient,dx,dy,reflat1,reflat2

  integer :: lenout
  integer, parameter :: lenbuff=20000000,maxrecs=1000000
  character(len=lenbuff) :: buffer
  character(len=256) :: cmd,line
  character(len=20) :: tokens(100)
  character(len=1) :: funit

  character(len=128) :: lines(maxrecs), wgrib_exe
  integer :: grib_ver

  character(4) :: sdunit

contains

!!*****************************************************************

  subroutine grib_get_recf(filein,inrec,a,nxy)

    implicit none

    integer,      intent( in) :: nxy
    real,         intent(out) :: a(nxy)
    character(*), intent( in) :: filein,inrec

    integer :: n,istart,iend ! ,nb,ne

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -d " // trim(inrec) // &
             " -text -nh -o - " // trim(filein) // char(0)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) // " -d " // trim(inrec) // &
             " -text - -no_header -inv /dev/null -order we:sn "// trim(filein) // char(0)
       iscan = 1   ! switch output order to s-n
    endif

    ! print*,'recf cmd: ',trim(cmd)

    call ir_popen(lenbuff,buffer,cmd,lenout)
    ! print*,'returned: ',lenout!,buffer(1:lenout)

    ! For ifort, the \n can't be included in the substring for the float read,
    !     so subtract 2

    istart=1
    do n=1,nxy
       iend=istart+index(buffer(istart:),char(10))-2
       !if(n<10) print*,'aa:',istart,iend,index(buffer(istart:),char(10))
       lines(n) = buffer(istart:iend)
       !if(n<10)   print*,'nnn:',n,istart,iend,ichar(buffer(4:4)),buffer(istart:iend)
       istart=iend+2
       !print*,'nnn:',nxy,n,trim(lines(n))
    enddo
    
!    do n = 1, nxy
!       read(lines(n),*) a(n)
!    enddo

    read(lines(1:nxy),*) a(1:nxy)

    ! print*,'done    : ',trim(cmd)
   
return
end subroutine grib_get_recf

!*****************************************************************

  subroutine grib_get_proj(filein)

    implicit none

    character(*), intent(in) :: filein
    integer                  :: n,istart,iend,ntok,nb,ne,noff

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -V -d 1 " // trim(filein) // char(0)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) //" -scan -grid -d 1 "// trim(filein) // char(0)
    endif
    call ir_popen(lenbuff,buffer,cmd,lenout)

    if (grib_ver == 1) then

       ! Take first seven lines from the output. Assuming all records
       ! have the same projection as first record.
       
       istart=1
       do n=1,7
          iend=istart+index(buffer(istart:),char(10))-2
          lines(n) = buffer(istart:iend)
          istart=iend+2
       enddo
   
       call tokenize_ws(lines(3),tokens,ntok)
       read(tokens(10),*) nx
       read(tokens(12),*) ny
   
       call tokenize_ws(lines(5),tokens,ntok)
       read(tokens(1),*) projection

       if(projection(1:5) == 'polar') then
          call tokenize_ws(lines(5),tokens,ntok)
          read(tokens(4),*) alat1
          read(tokens(6),*) alon1
          read(tokens(8),*) aorient
          call tokenize_ws(lines(6),tokens,ntok)
          read(tokens(1),*) cpole
          read(tokens(7),*) dx
          read(tokens(9),*) dy
          read(tokens(11),*) iscan
       elseif(projection(1:5) == 'Lambe') then
          call tokenize_ws(lines(5),tokens,ntok)
          read(tokens(4),*) alat1
          read(tokens(6),*) alon1
          read(tokens(8),*) alov
          call tokenize_ws(lines(6),tokens,ntok)
          read(tokens(2),*) reflat1
          read(tokens(4),*) reflat2
          call tokenize_ws(lines(7),tokens,ntok)
          read(tokens(7),*) dx
          read(tokens(9),*) dy
          read(tokens(11),*) iscan
       elseif(projection(1:5) == 'latlo') then
          call tokenize_ws(lines(5),tokens,ntok)
          read(tokens(3),*) alat1
          read(tokens(7),*) dx
          call tokenize_ws(lines(6),tokens,ntok)
          read(tokens(2),*) alon1
          read(tokens(6),*) dy
          read(tokens(11),*) iscan
       endif

    elseif (grib_ver == 2) then

       ! Parse 5 lines
   
       istart=1
       do n=1,5
          iend=istart+index(buffer(istart:),char(10))-2
          lines(n) = buffer(istart:iend)
          istart=iend+2
       enddo
   
       tokens = ' '
       call tokenize_ws(lines(1),tokens,ntok)
       nb = len_trim(tokens(1))
       ne = len_trim(tokens(1))
       read(tokens(1)(nb:ne),*) iscan
   
       tokens = ' '
       call tokenize_ws(lines(2),tokens,ntok)
       read(tokens(1),*) projection
       noff=0
       if (projection(1:7) == 'Lambert') noff = 1
       nb = index(tokens(2+noff),'(')+1
       ne = len_trim(tokens(2+noff))
       read(tokens(2+noff)(nb:ne),*) nx
       nb = 1
       ne = index(tokens(4+noff),')')-1
       read(tokens(4+noff)(nb:ne) ,*) ny

       !print*,'ppp0:', iscan,nx,ny,projection
       !print*,'ppp0:', alat1,dx,alon1,dy
   
       if(projection(1:5) == 'polar') then
          tokens = ' '
          call tokenize_ws(lines(3),tokens,ntok)
          read(tokens(2),*) alat1
          read(tokens(4),*) alon1
          read(tokens(6),*) aorient
          tokens = ' '
          call tokenize_ws(lines(4),tokens,ntok)
          read(tokens(1),*) cpole
          read(tokens(7),*) dx
          read(tokens(9),*) dy
          read(tokens(11),*) iscan
       elseif(projection(1:7) == 'Lambert') then
          !print*,'3--'//trim(lines(3))//'--'
          tokens = ' '
          call tokenize_ws(lines(3),tokens,ntok)
          read(tokens(2),*) alat1
          read(tokens(4),*) alon1
          read(tokens(6),*) alov
          tokens = ' '
          !print*,'4--'//trim(lines(4))//'--'
          call tokenize_ws(lines(4),tokens,ntok)
          read(tokens(2),*) reflat1
          read(tokens(4),*) reflat2
          tokens = ' '
          call tokenize_ws(lines(5),tokens,ntok)
          !print*,ntok,'5--'//trim(lines(5))//'--'
          !print*,'5--'//trim(tokens(1)), len_trim(tokens(1)),ichar(tokens(1))
          read(tokens(7),*) dx
          read(tokens(9),*) dy
          read(tokens(11),*) iscan
       elseif(projection(1:7) == 'lat-lon') then
          projection = 'latlon'
          tokens = ' '
          call tokenize_ws(lines(3),tokens,ntok)
          read(tokens(2),*) alat1
          read(tokens(6),*) dx
          tokens = ' '
          call tokenize_ws(lines(4),tokens,ntok)
          read(tokens(2),*) alon1
          read(tokens(6),*) dy
       else
          print*,'unknown projection:', '+'//trim(projection)//'+'
          stop 'bad proj'
       endif
   
    endif

    !print*,'ppp:', iscan,nx,ny,projection
    !print*,'ppp:', alat1,dx,alon1,dy

  end subroutine grib_get_proj

!**************************************************************

  subroutine grib_get_date(filein)

    implicit none

    character(*), intent(in) :: filein
    integer                  :: n,istart,iend,ntok,nb,ne

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -V -d 1 " // trim(filein) // char(0)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) // " -t -d 1 " // trim(filein) // char(0)
    endif
    
    call ir_popen(lenbuff,buffer,cmd,lenout)

    if (grib_ver == 1) then

       ! Parse 1 line

       istart=1
       do n=1,1
          iend=istart+index(buffer(istart:),char(10))-2
          lines(n) = buffer(istart:iend)
          istart=iend+2
       enddo
   
       call tokenize_ws(lines(1),tokens,ntok)
       read(tokens(3),*) longdate

    elseif (grib_ver == 2) then

       ! Parse 1 line
   
       istart=1
       do n=1,1
          iend=istart+index(buffer(istart:),char(10))-2
          lines(n) = buffer(istart:iend)
          istart=iend+2
       enddo
       
       call tokenize_ws(lines(1),tokens,ntok)
       nb = index(tokens(1),'=')+1
       ne = len_trim(tokens(1))
       read(tokens(1)(nb:ne),*) longdate
 
    endif

  end subroutine grib_get_date

!************************************************************************

  subroutine grib_get_fields (filein)
  
    implicit none

    character(*), intent(in) :: filein

    integer :: istart, n, iend, nlines, ntok, nc,nr
    character(len=20) :: clevel,chour,stokens(20), sdepths(20)
    real :: ztop, zbot

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -s "// trim(filein) // char(0)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) // " -s "// trim(filein) // char(0)
    endif
    call ir_popen(lenbuff,buffer,cmd,lenout)

!print*,lenout,buffer(1:lenout)

    ! Sepatate output into lines

    istart=1
    nlines=0
    do while(.true.)
       iend=istart+index(buffer(istart:),char(10))-2
       !print*,'lines:',nlines,istart,iend,lenout
       if (iend > lenout .or. istart > iend) exit
       nlines=nlines+1
       if (nlines > maxrecs) then
          print*, 'Number of records in grib file exceeds MAXRECS.'
          print*, '  MAXRECS=',maxrecs
          print*, '  last record read=',lines(nlines)
          print*, '  Program stopped.'
       endif
       lines(nlines) = buffer(istart:iend)
       istart=iend+2
    enddo

    allocate(irecs(nlines),levs(nlines),idates(nlines),  &
             ifhrs(nlines),fields(nlines))

    nrec=0
    do n=1,nlines

       ! Limit extraction to fields that are the surface or pressure levels
       if ( index(lines(n),' mb:' ) > 0         .or. &
            index(lines(n),':sfc:') > 0         .or. &
            index(lines(n),':surface:') > 0     .or. &
            index(lines(n),'below ground:') > 0 .or. &
            index(lines(n),' down:') > 0        ) then
          
          nrec=nrec+1
          call tokenize(lines(n),tokens,ntok,':')

          irecs(nrec) = tokens(1)
          read(tokens(3)(3:),*) idates(nrec)
          fields(nrec) = tokens(4) 
          clevel = tokens(5) 
          chour = tokens(6) 

          ! Skip accumulation or average fields
          if ( (index(chour,'acc') > 0) .or. (index(chour,'ave') > 0) ) then
             nrec = nrec - 1
             cycle
          endif

          ! The CFSR SOILM variable messes up the soil layers
          if ( fields(nrec) == 'SOILM' ) then
             nrec = nrec - 1
             cycle
          endif

          if(index(chour,'anl') > 0) then
             funit='h'
             ifhrs(nrec)=0
          elseif(index(chour,'hr') > 0 .or. index(chour,'hour') > 0) then
             funit='h'
             nc=index(chour,'h')-1
             read(chour(1:nc),*) ifhrs(nrec)
          elseif(index(chour,'min') > 0) then
             funit='m'
             nc=index(chour,'m')-1
             read(chour(1:nc),*) ifhrs(nrec)
          elseif(index(chour,'sec') > 0) then
             funit='s'
             nc=index(chour,'s')-1
             read(chour(1:nc),*) ifhrs(nrec)
          elseif(index(chour,'day') > 0) then
             funit='d'
             nc=index(chour,'d')-1
             read(chour(1:nc),*) ifhrs(nrec)
          else
             write(*,*) "Unknown field type", fields(nrec)
             nrec = nrec - 1
             cycle
          endif
      
          if ( index(lines(n),':sfc:') > 0 .or. index(lines(n),':surface:') > 0 ) then

             levs(nrec)=0

          elseif ( index(lines(n),'below ground:') > 0 .or. index(lines(n),' down:') > 0 ) then

             call tokenize(clevel, stokens, ntok, " ")
             call tokenize(stokens(1), sdepths, ntok, "-")

             if (ntok >= 2) then

                read(sdepths(1),*) ztop
                read(sdepths(2),*) zbot

                if (stokens(2) == 'cm') then
                
                   levs(nrec) = - nint(0.5*(ztop + zbot))
                   !write(*,*) "Found cm", levs(nrec), trim(fields(nrec)), ztop, zbot

                elseif (stokens(2) == 'm') then
                
                   levs(nrec) = - nint(50.0*(ztop + zbot))
                   !write(*,*) "Found m", levs(nrec), trim(fields(nrec))

                endif

             else

                write(*,*) "Skipping badly formatted soil variable " // trim(fields(nrec))
                nrec = nrec - 1
                cycle

             endif

          else

             read(clevel,*) levs(nrec)

          endif
      
          ! Check if this is a duplicate field. We will use the first occurence...
          do nr=1,nrec-1
             if ( idates(nrec) == idates(nr) .and. fields(nrec) == fields(nr) .and. &
                  levs(nrec) == levs(nr) .and. ifhrs(nrec) == ifhrs(nr) ) then
                print*, 'duplicate:',idates(nrec),fields(nrec),levs(nrec),ifhrs(nrec)
                nrec = nrec-1
                exit
             endif
          enddo
      
       endif
       
    enddo
    
  end subroutine grib_get_fields

!**********************************************************************

end module grib_get_mod
