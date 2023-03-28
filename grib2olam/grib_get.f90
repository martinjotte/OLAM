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
  real :: alat1,alat2,alon1,alov,aorient,dx,dy,reflat1,reflat2

  integer :: lenout

! integer, parameter :: lenbuff=20000000,maxrecs=1000000

  integer, parameter :: lenbuff = 200000
  integer, parameter :: maxrecs =  10000

  character(len=lenbuff) :: buffer
  character(len=256) :: cmd,line
  character(len=1) :: funit

  character(len=128) :: lines(maxrecs), wgrib_exe, wgrib1_exe, wgrib2_exe
  integer :: grib_ver

  character(4) :: sdunit

  real, allocatable :: alats(:)
  real, allocatable :: alons(:)

  real, parameter :: amiss  = -2.e20
  real, parameter :: amiss0 = -1.e20


  character(len=20), allocatable :: units


  ! Interface to C function ir_popen

  interface
     subroutine ir_popenc(nbuff,buff,cmd,len) bind(C,NAME="ir_popen")
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char
       integer  (kind=c_int )               :: nbuff
       type(c_ptr),            value        :: buff
       character(kind=c_char), dimension(*) :: cmd
       integer  (kind=c_int )               :: len
     end subroutine ir_popenc
  end interface

contains

!*****************************************************************

  subroutine grib_get_recf(filein,inrec,a,nx,ny)

    implicit none

    integer,      intent( in) :: nx, ny
    real,         intent(out) :: a(nx,ny)
    character(*), intent( in) :: filein, inrec
    integer                   :: isize

    if (grib_ver == 1) then

       cmd = trim(wgrib_exe) // " -d " // trim(inrec) // &
             " -bin -nh -o - " // trim(filein) // char(0)

    elseif (grib_ver == 2) then

       cmd = trim(wgrib_exe) // " -d " // trim(inrec) // &
             " -bin - -no_header -inv /dev/null -order we:sn "// trim(filein)

       iscan = 1   ! wgrib2 always defaults to s-n
    endif

    isize = nx * ny
    call ir_popen_real(isize, a, cmd, lenout)

  end subroutine grib_get_recf

!*****************************************************************

  subroutine grib_get_proj(filein)

    use gaussian, only: gauaw
    use token_mod

    implicit none

    character(*), intent(in) :: filein
    integer                  :: n,istart,iend,ntok,nb,ne,noff
    character(len=20)        :: tokens(20)

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -V -d 1 " // trim(filein) // char(0)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) //" -scan -grid -d 1 "// trim(filein)
    endif
    call ir_popen_chars(lenbuff,buffer,cmd,lenout)

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

       ne = len_trim(projection)
       if (projection(ne:ne) == ':') projection = projection(1:ne-1)

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
          read(tokens(5),*) alat2
          read(tokens(7),*) dy

          ! Handle case of reversed lats (N-S)
          alat1 = min(alat1,alat2)

          allocate(alats(ny))
          alats(1) = alat1
          do n=2,ny
             alats(n) = alats(n-1) + dy
          enddo

          call tokenize_ws(lines(6),tokens,ntok)
          read(tokens(2),*) alon1
          read(tokens(6),*) dx
          read(tokens(11),*) iscan

          allocate(alons(nx))
          alons(1) = alon1
          do n = 2, nx
             alons(n) = alons(n-1) + dx
          enddo

       elseif(projection(1:5) == 'gauss') then

          allocate(alats(ny))
          call gauaw(alats, ny)
          alat1 = alats(1)
          dy = 180.0 / ny

          call tokenize_ws(lines(6),tokens,ntok)
          read(tokens(2),*) alon1
          read(tokens(6),*) dx
          read(tokens(11),*) iscan

          allocate(alons(nx))
          alons(1) = alon1
          do n = 2, nx
             alons(n) = alons(n-1) + dx
          enddo

       endif

    elseif (grib_ver == 2) then

       ! Parse 5 lines

       istart=1
       do n=1,5
          iend=istart+index(buffer(istart:),char(10))-2
          lines(n) = buffer(istart:iend)
          istart=iend+2
       enddo

       iscan = 1

!!       tokens = ' '
!!       call tokenize_ws(lines(1),tokens,ntok)
!!       nb = len_trim(tokens(1))
!!       ne = len_trim(tokens(1))
!!       read(tokens(1)(nb:ne),*) iscan

       tokens = ' '
       call tokenize_ws(lines(2),tokens,ntok)
       read(tokens(1),*) projection
       noff=0
       if (projection(1:7) == 'Lambert')  noff = 1
       if (projection(1:8) == 'Gaussian') noff = 1
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
!         read(tokens(11),*) iscan
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
!         read(tokens(11),*) iscan

       elseif(projection(1:7) == 'lat-lon') then

          projection = 'latlon'
          tokens = ' '
          call tokenize_ws(lines(3),tokens,ntok)
          read(tokens(2),*) alat1
          read(tokens(4),*) alat2
          read(tokens(6),*) dy

          ! wgrib2 always outputs data S to N, regardless of how it is stored
          ! in the grib file. So, the south lat can either be alat1 or alat2!
          alat1 = min(alat1,alat2)
          allocate(alats(ny))
          alats(1) = alat1
          do n=2,ny
             alats(n) = alats(n-1) + dy
          enddo

          tokens = ' '
          call tokenize_ws(lines(4),tokens,ntok)
          read(tokens(2),*) alon1
          read(tokens(6),*) dx
          allocate(alons(nx))
          alons(1) = alon1
          do n = 2, nx
             alons(n) = alons(n-1) + dx
          enddo

       elseif(projection(1:7) == 'Gaussian') then

          projection = 'gaussian'
          tokens = ' '
          allocate(alats(ny))
          call gauaw(alats, ny)
          alat1 = alats(1)
          dy = 180.0 / ny

          call tokenize_ws(lines(5),tokens,ntok)
          read(tokens(2),*) alon1
          read(tokens(6),*) dx
          allocate(alons(nx))
          alons(1) = alon1
          do n = 2, nx
             alons(n) = alons(n-1) + dx
          enddo

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

    use token_mod
    implicit none

    character(*), intent(in) :: filein
    integer                  :: n,istart,iend,ntok,nb,ne
    character(len=20)        :: tokens(20)

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -V -d 1 " // trim(filein)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) // " -t -d 1 " // trim(filein)
    endif

    call ir_popen_chars(lenbuff,buffer,cmd,lenout)

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

    use token_mod
    implicit none

    character(*), intent(in) :: filein

    integer :: istart, n, iend, nlines, ntok, ntok2, nc, nr, nvlev
    character(len=20) :: chour, tokens(20), sdepths(20), tok(10), vlevs(10)
    real :: ztop, zbot, pval

    if (grib_ver == 1) then
       cmd = trim(wgrib_exe) // " -s "// trim(filein)
    elseif (grib_ver == 2) then
       cmd = trim(wgrib_exe) // " -s "// trim(filein)
    endif
    call ir_popen_chars(lenbuff,buffer,cmd,lenout)

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
          chour = tokens(6)

          call tokenize(tokens(5), vlevs, nvlev, " ")

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

             call tokenize(vlevs(1), sdepths, ntok2, "-")

             ! NCEP gives a range for soil layers (like 0.5-1.0 m)
             ! Need to check how other reanalysis products handle this
             if (ntok2 > 1 .and. (vlevs(2) == 'cm' .or. vlevs(2) == 'm')) then

                read(sdepths(1),*) ztop
                read(sdepths(2),*) zbot

                if (vlevs(2) == 'cm') then

                   levs(nrec) = - nint(0.5*(ztop + zbot))

                elseif (vlevs(2) == 'm') then

                   ! convert soil depth to cm
                   levs(nrec) = - nint(50.0*(ztop + zbot))

                endif

             else

                write(*,*) "Skipping badly formatted soil variable " // trim(fields(nrec))
                nrec = nrec - 1
                cycle

             endif

          elseif (vlevs(2) == 'mb') then
             ! TODO: Are any reanalysis in height or not in mb?

             read(vlevs(1),*) pval

             ! convert to Pa
             pval = pval * 100.0

             if (pval < 0.99) then

                ! Skip pressure levels below 1 Pa, since we curently
                ! store pressure as an integer value
                nrec = nrec - 1

             else

                levs(nrec) = nint( pval )

             endif

          else

             ! Skip this record since we didn't match any fields
             nrec = nrec - 1

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

  subroutine grib_get_vers(filein)

    implicit none

    character(*), intent(in) :: filein

    write(*,*) "Checking grib file version..."

    cmd = trim(wgrib1_exe) // " 2>&1 "
    call ir_popen_chars(lenbuff,buffer,cmd,lenout)
    if (index(buffer, "wgrib ") == 0) then
       write(*,*)
       write(*,'(A)') "The wgrib command " // trim(wgrib1_exe) // " does not appear to be working."
       write(*,'(A)') "Compile the wgrib executable and specify the proper path in the namelist."
       stop
    endif

    grib_ver = 1
    cmd = trim(wgrib1_exe) // " -V -d 1 " // trim(filein) // " 2>&1 "

    call ir_popen_chars(lenbuff,buffer,cmd,lenout)

    if (index(buffer, "grib2 message ignored") /= 0) grib_ver = 0
    if (index(buffer, "missing GRIB record")   /= 0) grib_ver = 0

    if (grib_ver == 1) then

       write(*,*)
       write(*,*) "We have a grib1 file, using " // trim(wgrib1_exe) // " to read file."
       wgrib_exe = wgrib1_exe
       return
    else

       write(*,*)
       write(*,*) "We do not have a grib1 file that can be read with wgrib."
       write(*,*) "Checking if we have a grib2 file that can be read with wgrib2..."

       cmd = trim(wgrib2_exe) // " 2>&1 "
       call ir_popen_chars(lenbuff,buffer,cmd,lenout)
       if (index(buffer, "wgrib2 ") == 0) then
          write(*,*)
          write(*,'(A)') "The wgrib2 command " // trim(wgrib2_exe) // " does not appear to be working."
          write(*,'(A)') "Compile the included wgrib2 executable and specify the proper path in the namelist."
          stop
       endif

       cmd = trim(wgrib2_exe) // " -t -d 1 " // trim(filein) // " 2>&1 "
       call ir_popen_chars(lenbuff,buffer,cmd,lenout)

       if (index(buffer, "FATAL ERROR") == 0) then

          write(*,*)
          write(*,*) "We have a grib2 file, using " // trim(wgrib2_exe) // " to read file."
          grib_ver = 2
          wgrib_exe = wgrib2_exe
          return

       else

          write(*,*)
          write(*,*) "Cannot recognize file " // trim(filein)
          write(*,*) "Will not process grib file."
          stop

       endif

    endif

  end subroutine grib_get_vers

!************************************************************************

  function f2cstring(string)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    implicit none

    character(len=*), intent(in) :: string
    character(len=len_trim(string)+1, kind=c_char) :: f2cstring

    f2cstring = trim(string) // c_null_char
  end function f2cstring

!************************************************************************

  subroutine ir_popen_integer(nbuff, buff, cmd, len)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_loc, c_sizeof
    implicit none

    integer,      intent(in)            :: nbuff
    integer,      intent(inout), target :: buff(*)
    character(*), intent(in)            :: cmd
    integer,      intent(in)            :: len
    type(c_ptr)                         :: buffc
    integer                             :: ibuff
    character(len=len_trim(cmd)+1, kind=c_char) :: cstring

    buffc = c_loc(buff(1))
    ibuff = nbuff * c_sizeof(1)

    cstring = f2cstring(cmd)
    call ir_popenc(ibuff, buffc, cstring, len)
  end subroutine ir_popen_integer

!************************************************************************

  subroutine ir_popen_real(nbuff, buff, cmd, len)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_loc, c_sizeof
    implicit none

    integer,      intent(in)            :: nbuff
    real,         intent(inout), target :: buff(*)
    character(*), intent(in)            :: cmd
    integer,      intent(in)            :: len
    type(c_ptr)                         :: buffc
    integer                             :: ibuff
    character(len=len_trim(cmd)+1, kind=c_char) :: cstring

    buffc = c_loc(buff(1))
    ibuff = nbuff * c_sizeof(1.0)

    cstring = f2cstring(cmd)
    call ir_popenc(ibuff, buffc, cstring, len)
  end subroutine ir_popen_real

!************************************************************************

  subroutine ir_popen_chars(nbuff, buff, cmd, len)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_loc
    implicit none

    integer,      intent(in)                          :: nbuff
    character,    intent(inout), dimension(*), target :: buff
    character(*), intent(in)                          :: cmd
    integer,      intent(in)                          :: len
    type(c_ptr)                                       :: buffc
    character(len=len_trim(cmd)+1, kind=c_char)       :: cstring

    buffc = c_loc(buff(1))

    cstring = f2cstring(cmd)
    call ir_popenc(nbuff, buffc, cstring, len)
  end subroutine ir_popen_chars

end module grib_get_mod
