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
subroutine history_write(vtype)

  use var_tables, only: num_var, vtab_r, get_vtab_dims
  use misc_coms,  only: io6, ioutput, hfilepref, current_time, iclobber, iparallel
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use max_dims,   only: pathlen
  use mem_para,   only: myrank

  implicit none

! This routine writes the chosen variables on the history file.

  character(*), intent(in) :: vtype  ! not used yet - 'state'

  character(pathlen) :: hnamel
  character(32)      :: varn
  character(10)      :: post
  logical            :: exans
  integer            :: nv, nvcnt
  integer            :: ndims, idims(3)

  if (ioutput == 0) return

! Set filename post depending on whether run is parallel

  if (iparallel == 0) then
     post = '$'
  else
     write(post,'(i10)') myrank
     post = 'r'//trim(adjustl(post))
  endif

! Construct h5 file name and open the file

  call makefnam(hnamel, hfilepref, current_time, 'H', post, 'h5')

  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open file name :'
     write(io6,*) '!!!       '//trim(hnamel)
     write(io6,*) '!!!   but it already exists. run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'history_write'
  endif

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'history_write: opening file: ',trim(hnamel)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(hnamel,'W',iclobber)

! Write the common fields

  call commio('WRITE')

! Loop through the main variable table and write those variables
! with the correct flag set

  nvcnt = 0
  do nv = 1,num_var

     if (vtab_r(nv)%ihist) then

        varn = vtab_r(nv)%name
        call get_vtab_dims(nv, ndims, idims)

        write (io6, '(1x,a,2(I0,1x),a,3(1x,I0))')  &
             'Writing: ', nv, num_var, trim(varn), idims(1:ndims)

        if     (associated(vtab_r(nv)%ivar0_p)) then
           call shdf5_orec(ndims, idims, trim(varn), ivars=vtab_r(nv)%ivar0_p)
        elseif (associated(vtab_r(nv)%ivar1_p)) then
           call shdf5_orec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar1_p)
        elseif (associated(vtab_r(nv)%ivar2_p)) then
           call shdf5_orec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar2_p)
        elseif (associated(vtab_r(nv)%ivar3_p)) then
           call shdf5_orec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar3_p)

        elseif (associated(vtab_r(nv)%rvar0_p)) then
           call shdf5_orec(ndims, idims, trim(varn), rvars=vtab_r(nv)%rvar0_p)
        elseif (associated(vtab_r(nv)%rvar1_p)) then
           call shdf5_orec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar1_p)
        elseif (associated(vtab_r(nv)%rvar2_p)) then
           call shdf5_orec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar2_p)
        elseif (associated(vtab_r(nv)%rvar3_p)) then
           call shdf5_orec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar3_p)

        elseif (associated(vtab_r(nv)%dvar0_p)) then
           call shdf5_orec(ndims, idims, trim(varn), dvars=vtab_r(nv)%dvar0_p)
        elseif (associated(vtab_r(nv)%dvar1_p)) then
           call shdf5_orec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar1_p)
        elseif (associated(vtab_r(nv)%dvar2_p)) then
           call shdf5_orec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar2_p)
        elseif (associated(vtab_r(nv)%dvar3_p)) then
           call shdf5_orec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar3_p)
        endif

        nvcnt = nvcnt + 1

     endif

  enddo

  call shdf5_close()

end subroutine history_write
