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

subroutine history_start(action)

  ! This routine initializes the model from the history file

  use misc_coms,  only: io6, hfilin, time8, time_istp8, runtype, &
                        iparallel, isubdomain
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_para,   only: myrank, mgroupsize
  use max_dims,   only: pathlen

  implicit none

  character(*), intent(in) :: action
  logical                  :: exans

  ! List of all allowable runtype, seq-para run, histfile read combinations:

  ! 'HISTORY'     , sequential run, sequential histfile  (isubdomain = 0)
  ! 'HISTORY'     , parallel   run, sequential histfile  (isubdomain = 1)
  ! 'PLOTONLY'    , sequential run, sequential histfile  (isubdomain = 0)
  ! 'PLOTONLY'    , parallel   run, sequential histfile  (isubdomain = 1)

  !future 'EXTRACT'     , sequential run, sequential histfile  (isubdomain = 1)
  !future 'PLOTEXTRACT' , sequential run, extracted  histfile  (isubdomain = 1)

  ! Check if history files exist
  
  inquire(file=hfilin, exist=exans)   ! global restart file

  if (exans) then

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening history file '//trim(hfilin)//' for '//trim(action)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(hfilin,'R')

     if (trim(action) == 'COMMIO') then

        ! Read the common variables
        call commio('READ')

     else

        ! Read the model fields
        call shdf5_irec(1, (/1/), 'time8', dvars=time8)

        if (isubdomain == 1) then
           call hist_read_subd()
        else
           call hist_read_glob()
        endif

        write(io6,*) 'back from hist_read'

     endif

     time_istp8 = time8  ! time_istp8 is used for plotting time
     call shdf5_close()

  else

     ! History files do not exist, stop model.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open history file file:'
     write(io6,*) '!!!   '//trim(hfilin)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'in history_start'

  endif

end subroutine history_start

!=========================================================================

subroutine hist_read_glob()

  ! Sequential run reading all points from a global history file

  use misc_coms,  only: io6, hfilin, runtype
  use var_tables, only: num_var, vtab_r, get_vtab_dims
  use hdf5_utils, only: shdf5_open, shdf5_info, shdf5_irec, shdf5_close

  implicit none

  integer           :: nv, nvcnt, ndims, ndims_m, idims(3), idims_m(3)
  character(len=32) :: varn

  nvcnt =  0

  ! Loop through all variables in the model vtables

  do nv = 1, num_var

     ndims = -1
     idims =  0

     ! Skip to next variable if we don't want the current one

     if (.not. vtab_r(nv)%ihist) cycle
     if (vtab_r(nv)%nread .and. runtype == 'HISTORY') cycle
   
     varn = trim(vtab_r(nv)%name)
     call get_vtab_dims(nv, ndims_m, idims_m)

     ! We want it...read it if it's in the history file

     call shdf5_info(varn, ndims, idims)

     ! Skip to next variable if the current one is not in the history file

     if (ndims < 0) then
        write(io6,*)
        write(io6,*) 'Variable '//trim(varn)//' is not in the history file, skipping'
        write(io6,*)
        cycle
     endif

     ! Report an error if the variable dimensions are different

     if (ndims /= ndims_m .or. any(idims(1:ndims) /= idims_m(1:ndims))) then
        write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(io6,*) '!!!   Trying to read variable '//trim(varn)
        write(io6,*) '!!!   from the history file '//trim(hfilin)
        write(io6,*) '!!!   but the array size is different. '
        write(io6,*) '!!!   ndims, ndims_m = ',ndims,ndims_m
        write(io6,*) '!!!   idims  (1:ndims) = ',idims  (1:ndims)
        write(io6,*) '!!!   idims_m(1:ndims) = ',idims_m(1:ndims)
        write(io6,*) '!!!   The run is ended'
        write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        stop    'stop in hist_read'
     endif

     if     (associated(vtab_r(nv)%ivar1_p)) then
        call shdf5_irec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar1_p)
     elseif (associated(vtab_r(nv)%ivar2_p)) then
        call shdf5_irec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar2_p)
     elseif (associated(vtab_r(nv)%ivar3_p)) then
        call shdf5_irec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar3_p)

     elseif (associated(vtab_r(nv)%rvar1_p)) then
        call shdf5_irec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar1_p)
     elseif (associated(vtab_r(nv)%rvar2_p)) then
        call shdf5_irec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar2_p)
     elseif (associated(vtab_r(nv)%rvar3_p)) then
        call shdf5_irec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar3_p)

     elseif (associated(vtab_r(nv)%dvar1_p)) then
        call shdf5_irec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar1_p)
     elseif (associated(vtab_r(nv)%dvar2_p)) then
        call shdf5_irec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar2_p)
     elseif (associated(vtab_r(nv)%dvar3_p)) then
        call shdf5_irec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar3_p)
     endif
   
     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)
  enddo

end subroutine hist_read_glob

!=========================================================================

subroutine hist_read_subd()

  ! Parallel node reading portion of a global history file,
  ! or a serial run reading a subdomain

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza
  use mem_ijtabs,  only: itabg_w, itab_w, itab_v, itab_m
  use misc_coms,   only: io6, runtype
  use var_tables,  only: num_var, vtab_r, get_vtab_dims
  use hdf5_utils,  only: shdf5_info, shdf5_irec
  use mem_leaf,    only: itab_wl
  use leaf_coms,   only: nwl, mwl, nzg
  use mem_sea,     only: itab_ws
  use sea_coms,    only: nws, mws
  use mem_para,    only: myrank
  use consts_coms, only: r8
  use mem_nudge,   only: nwnud, mwnud, itab_wnud

  implicit none

  integer          :: nv, nvcnt, ndims, idims(3)
  character(32)    :: varn
  character (2)    :: stagpt
  integer, pointer :: ilocal(:)

  nvcnt =  0

  ! Loop through all variables in the model vtables

  do nv = 1, num_var

     ndims = -1
     idims =  0

     ! Skip to next variable if we don't want the current one

     if (.not. vtab_r(nv)%ihist) cycle
     if (vtab_r(nv)%nread) cycle

     varn    = trim(vtab_r(nv)%name)
     stagpt  = vtab_r(nv)%stagpt

     ! We want it...read it if it's in the history file

     call shdf5_info(varn, ndims, idims)

     ! Skip to next variable if the current one is not in the history file

     if (ndims < 0) then
        write(io6,*)
        write(io6,*) 'Variable '//trim(varn)//' is not in the history file, skipping'
        write(io6,*)
        cycle
     endif

     ! Identify the points we want to read from the history file

     if     (stagpt == 'AW' .and. idims(ndims) == nwa) then
        ilocal => itab_w(:)%iwglobe
        idims(ndims) = mwa
     elseif (stagpt == 'AV' .and. idims(ndims) == nva) then
        ilocal => itab_v(:)%ivglobe
        idims(ndims) = mva
     elseif (stagpt == 'AM' .and. idims(ndims) == nma) then
        ilocal => itab_m(:)%imglobe
        idims(ndims) = mma
     elseif (stagpt == 'LW' .and. idims(ndims) == nwl) then
        ilocal => itab_wl(:)%iwglobe
        idims(ndims) = mwl
     elseif (stagpt == 'SW' .and. idims(ndims) == nws) then
        ilocal => itab_ws(:)%iwglobe
        idims(ndims) = mws
     elseif (stagpt == 'AN' .and. idims(ndims) == nwnud) then
        ilocal => itab_wnud(:)%iwnudglobe
        idims(ndims) = mwnud
     else

        ! TODO: Const values
        ! TODO: Land U and M; Sea U and M? (probably not!)

        stop "invalid array size in history_start_s"
     endif

     if     (associated(vtab_r(nv)%ivar1_p)) then
        call shdf5_irec(ndims, idims, varn, ivara=vtab_r(nv)%ivar1_p, points=ilocal)
     elseif (associated(vtab_r(nv)%ivar2_p)) then
        call shdf5_irec(ndims, idims, varn, ivara=vtab_r(nv)%ivar2_p, points=ilocal)
     elseif (associated(vtab_r(nv)%rvar1_p)) then
        call shdf5_irec(ndims, idims, varn, rvara=vtab_r(nv)%rvar1_p, points=ilocal)
     elseif (associated(vtab_r(nv)%rvar2_p)) then
        call shdf5_irec(ndims, idims, varn, rvara=vtab_r(nv)%rvar2_p, points=ilocal)
     elseif (associated(vtab_r(nv)%dvar1_p)) then
        call shdf5_irec(ndims, idims, varn, dvara=vtab_r(nv)%dvar1_p, points=ilocal)
     elseif (associated(vtab_r(nv)%dvar2_p)) then
        call shdf5_irec(ndims, idims, varn, dvara=vtab_r(nv)%dvar2_p, points=ilocal)
     endif

     nvcnt = nvcnt + 1
     write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
           'Read: ', nvcnt, trim(varn), idims(1:ndims)

  enddo

end subroutine hist_read_subd
