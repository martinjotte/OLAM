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
Module hdf5_utils

Contains

subroutine shdf5_open(locfn, access, idelete)

  use hdf5_f2f
  implicit none

  character(*),      intent(in) :: locfn   ! file name
  character(*),      intent(in) :: access  ! File access ('R','W','RW')
  integer, optional, intent(in) :: idelete ! If W, overwrite file if exists?
                                           !    1=yes, 0=no

  integer                       :: hdferr  ! Error flag
  integer                       :: iaccess ! int access flag
  logical                       :: exists  ! File existence

! Check for existence of RAMS file.

  inquire(file=locfn, exist=exists)

! Create a new file or open an existing RAMS file.

  if (access(1:1) == 'R') then

     if (.not. exists) then
        print*, 'shdf5_open:'
        print*, '   Attempt to open a file for reading that does not exist.'
        print*, '   Filename: ', trim(locfn)
        stop    'shdf5_open: no file'
     else
        if (access == 'R ') iaccess = 1
        if (access == 'RW') iaccess = 2
        call fh5f_open(locfn, iaccess, hdferr)
        
        if (hdferr < 0) then
           print*, 'shdf5_open:'
           print*, '   Error opening hdf5 file - error -', hdferr
           print*, '   Filename: ',trim(locfn)
           stop    'shdf5_open: open error'      
        endif
     endif
     
  elseif (access(1:1) == 'W') then

     if (.not. exists) then
        iaccess = 2
        call fh5f_create(locfn, iaccess, hdferr)
     else
        if (.not. present(idelete) ) then
           print*, 'shdf5_open: idelete not specified when access=W'
           stop    'shdf5_open: no idelete'
        endif
      
        if (idelete == 0) then
           print*, 'In shdf5_open:'
           print*, '   Attempt to open an existing file for writing,'
           print*, '      but overwrite is disabled. idelete=', idelete
           print*, '   Filename: ', trim(locfn)
           stop    'shdf5_open'
        else
           iaccess = 1
           call fh5f_create(locfn, iaccess, hdferr)
        endif

     endif
   
     if(hdferr < 0) then
        print*, 'HDF5 file create failed:', hdferr
        print*, 'file name:', trim(locfn), ' ', trim(access), idelete
        stop    'shdf5_open: bad create'
     endif
  endif

end subroutine shdf5_open

!===============================================================================

subroutine shdf5_info(dsetname, ndims, dims)
  use hdf5_f2f

  implicit none

  character(*), intent(in)    :: dsetname ! Dataset name
  integer,      intent(inout) :: dims(:)
  integer,      intent(inout) :: ndims    ! Dataset rank (in file)
  integer                     :: hdferr   ! Error flag

! Open the dataset.

  call fh5d_open(dsetname, hdferr)

  if (hdferr < 0) then
     print*, 'In shdf5_info:'
     print*, 'Variable ', trim(dsetname), ' is not in the currently opened hdf5 file'
     ndims   = 0
     dims(1) = 0
     return
  endif

! Get dataset's dimensions

  call fh5s_get_ndims(ndims)
  call fh5s_get_dims(dims)

  call fh5d_close(hdferr)

end subroutine shdf5_info

!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,  &
                                          ivars,rvars,cvars,dvars,lvars,  &
                                          nglobe, lpoints, gpoints,       &
                                          units, long_name, positive,     &
                                          imissing, rmissing, dmissing,   &
                                          isdim, dimnames, standard_name, &
                                          cell_methods                    )

  use oname_coms,  only: nl
  use consts_coms, only: r8
  use hdf5_f2f
    
  implicit none

  character(*), intent(in) :: dsetname ! Variable label
  integer,      intent(in) :: ndims    ! Number of dimensions or rank
  integer,      intent(in) :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer,      intent(in), optional :: ivara(*), ivars
  real,         intent(in), optional :: rvara(*), rvars
  character,    intent(in), optional :: cvara(*), cvars
  real(r8),     intent(in), optional :: dvara(*), dvars
  logical,      intent(in), optional :: lvara(*), lvars

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional :: lpoints(:), gpoints(:), nglobe

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional :: dimnames(:)

! Local variables
  integer :: hdferr  ! Error flag

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec:', ndims, dims(1:ndims)
     stop    'shdf5_orec: bad dims'
  endif
     
! Prepare memory and options for the write

  call fh5_prepare_write(ndims, dims, hdferr, nl%icompress, &
       mcoords=lpoints, fcoords=gpoints, ifsize=nglobe)

  if (hdferr /= 0) then
     print*, "shdf5_orec: can't prepare requested field:", trim(dsetname)
     return
  endif

! Write the dataset.

      if (present(ivars)) then ; call fh5_write(ivars, dsetname, hdferr)
  elseif (present(rvars)) then ; call fh5_write(rvars, dsetname, hdferr)
  elseif (present(cvars)) then ; call fh5_write(cvars, dsetname, hdferr)
  elseif (present(dvars)) then ; call fh5_write(dvars, dsetname, hdferr)
  elseif (present(lvars)) then ; call fh5_write(lvars, dsetname, hdferr)
  elseif (present(ivara)) then ; call fh5_write(ivara, dsetname, hdferr)
  elseif (present(rvara)) then ; call fh5_write(rvara, dsetname, hdferr)
  elseif (present(cvara)) then ; call fh5_write(cvara, dsetname, hdferr)
  elseif (present(dvara)) then ; call fh5_write(dvara, dsetname, hdferr)
  elseif (present(lvara)) then ; call fh5_write(lvara, dsetname, hdferr)
  else
     print*, 'Incorrect or missing data field argument in shdf5_orec'
     stop    'shdf5_orec: bad data field'
  endif

  if (hdferr /= 0) then
     print*, 'In shdf5_orec: hdf5 write error =', hdferr
     stop    'shdf5_orec: hdf5 write error'
  endif

! Indicate if this variable is a dimension

  if (present(isdim)) then
     if (isdim) call fh5f_create_dim(dsetname, hdferr)
  endif

! Link each variable to its dimensions
  
  if ((present(dimnames)) .and. (.not. present(isdim))) then
     if (size(dimnames) >= ndims) call fh5f_attach_dims(ndims, dimnames, hdferr)
  endif

! Write any dataset attributes to match common NetCDF conventions

  if (present(units)) then
     if (len_trim(units) > 0) call fh5f_write_attribute("units", cvalue=units)
  endif

  if (present(long_name)) then
     if (len_trim(long_name) > 0) call fh5f_write_attribute("long_name", cvalue=long_name)
  endif

  if (present(standard_name)) then
     if (len_trim(standard_name) > 0) call fh5f_write_attribute("standard_name", cvalue=standard_name)
  endif

  if (present(cell_methods)) then
     if (len_trim(cell_methods) > 0) call fh5f_write_attribute("cell_methods", cvalue=cell_methods)
  endif

  if (present(positive)) then
     if (len_trim(positive) > 0) call fh5f_write_attribute("positive", cvalue=positive)
  endif

  if ((present(ivars) .or. present(ivara)) .and. present(imissing)) then
     call fh5f_write_attribute("missing_value", ivalue=imissing)
  endif

  if ((present(rvars) .or. present(rvara)) .and. present(rmissing)) then
     call fh5f_write_attribute("missing_value", rvalue=rmissing)
  endif

  if ((present(dvars) .or. present(dvara)) .and. present(dmissing)) then
     call fh5f_write_attribute("missing_value", dvalue=dmissing)
  endif

! Close the dataset, the dataspace for the dataset, and the dataspace properties.

  call fh5_close_write(hdferr)

end subroutine shdf5_orec

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                         ,ivars,rvars,cvars,dvars,lvars  &
                                         ,points,start,counts)
  use hdf5_f2f
  implicit none

  character(*), intent(IN) :: dsetname ! Dataset name
  integer,      intent(IN) :: ndims    ! Number of dimensions or rank
  integer,      intent(IN) :: dims(:)  ! Dataset dimensions

! Array and scalar arguments for different types. Only specify one in each call.
  integer,      intent(OUT), optional :: ivara(*), ivars
  real,         intent(OUT), optional :: rvara(*), rvars
  character,    intent(OUT), optional :: cvara(*), cvars
  real(kind=8), intent(OUT), optional :: dvara(*), dvars
  logical,      intent(OUT), optional :: lvara(*), lvars

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(IN),  optional :: points(:)
  integer,      intent(IN),  optional :: start (:)
  integer,      intent(IN),  optional :: counts(:)

! Local variables
  integer :: hdferr  ! Error flag

! Check dimensions

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_irec:', ndims, dims(1:ndims)
     stop    'shdf5_irec: bad dims'
  endif
    
! Prepare file and memory space for the read

  call fh5_prepare_read(dsetname, ndims, dims, hdferr, coords=points, &
                        start=start, counts=counts)
  if (hdferr < 0) then
     print*,'shdf5_irec: can''t prepare requested field:',trim(dsetname)
     return
  endif

! Read data from hyperslab in the file into the hyperslab in memory.

      if (present(ivars)) then ; call fh5d_read(ivars, hdferr)
  elseif (present(rvars)) then ; call fh5d_read(rvars, hdferr)
  elseif (present(cvars)) then ; call fh5d_read(cvars, hdferr)
  elseif (present(dvars)) then ; call fh5d_read(dvars, hdferr)
  elseif (present(lvars)) then ; call fh5d_read(lvars, hdferr)
  elseif (present(ivara)) then ; call fh5d_read(ivara, hdferr)
  elseif (present(rvara)) then ; call fh5d_read(rvara, hdferr)
  elseif (present(cvara)) then ; call fh5d_read(cvara, hdferr)
  elseif (present(dvara)) then ; call fh5d_read(dvara, hdferr)
  elseif (present(lvara)) then ; call fh5d_read(lvara, hdferr)
  else
     print*,'Incorrect or missing data field argument in shdf5_irec'
     print*, 'field = ', dsetname
     stop    'shdf5_irec: bad data field'
  endif
  
  if (hdferr /= 0) then
     print*, 'shdf5_irec: call fh5d_read: hdf5 error =', hdferr
     print*, 'Error reading ', trim(dsetname)
     print*, 'ndims = ', ndims
     print*, 'dims  = ', dims(1:ndims)
     stop
  endif

! Close the dataset, the dataspace for the dataset, and the memory space.

  call fh5_close_read(hdferr)

end subroutine shdf5_irec

!===============================================================================

subroutine shdf5_close()
  use hdf5_f2f
  implicit none

  integer :: hdferr  ! Error flags

! Close hdf file.

  call fh5f_close(hdferr)
end  subroutine shdf5_close

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara  &
                                              ,ivars,rvars,cvars,dvars,lvars)
  use hdf5_f2f
  implicit none

  character(*), intent(in)              :: dsetname, action
  integer,      intent(in)              :: ndims,    dims(:)
  integer,      intent(inout), optional :: ivara(*), ivars
  real,         intent(inout), optional :: rvara(*), rvars
  character,    intent(inout), optional :: cvara(*), cvars
  real(kind=8), intent(inout), optional :: dvara(*), dvars
  logical,      intent(inout), optional :: lvara(*), lvars
 
  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (action == 'READ') then
     
         if (present(ivars)) then ; call shdf5_irec(ndims, dims, dsetname, ivars=ivars)
     elseif (present(rvars)) then ; call shdf5_irec(ndims, dims, dsetname, rvars=rvars)
     elseif (present(cvars)) then ; call shdf5_irec(ndims, dims, dsetname, cvars=cvars)
     elseif (present(dvars)) then ; call shdf5_irec(ndims, dims, dsetname, dvars=dvars)
     elseif (present(lvars)) then ; call shdf5_irec(ndims, dims, dsetname, lvars=lvars)
     elseif (present(ivara)) then ; call shdf5_irec(ndims, dims, dsetname, ivara=ivara)
     elseif (present(rvara)) then ; call shdf5_irec(ndims, dims, dsetname, rvara=rvara)
     elseif (present(cvara)) then ; call shdf5_irec(ndims, dims, dsetname, cvara=cvara)
     elseif (present(dvara)) then ; call shdf5_irec(ndims, dims, dsetname, dvara=dvara)
     elseif (present(lvara)) then ; call shdf5_irec(ndims, dims, dsetname, lvara=lvara)
     else
        print*,'Incorrect or missing data field argument in shdf5_io'
        print*, 'field = ', dsetname
        stop    'shdf5_io: bad data field'
     endif
     
  elseif (action == 'WRITE') then
     
         if (present(ivars)) then ; call shdf5_orec(ndims, dims, dsetname, ivars=ivars)
     elseif (present(rvars)) then ; call shdf5_orec(ndims, dims, dsetname, rvars=rvars)
     elseif (present(cvars)) then ; call shdf5_orec(ndims, dims, dsetname, cvars=cvars)
     elseif (present(dvars)) then ; call shdf5_orec(ndims, dims, dsetname, dvars=dvars)
     elseif (present(lvars)) then ; call shdf5_orec(ndims, dims, dsetname, lvars=lvars)
     elseif (present(ivara)) then ; call shdf5_orec(ndims, dims, dsetname, ivara=ivara)
     elseif (present(rvara)) then ; call shdf5_orec(ndims, dims, dsetname, rvara=rvara)
     elseif (present(cvara)) then ; call shdf5_orec(ndims, dims, dsetname, cvara=cvara)
     elseif (present(dvara)) then ; call shdf5_orec(ndims, dims, dsetname, dvara=dvara)
     elseif (present(lvara)) then ; call shdf5_orec(ndims, dims, dsetname, lvara=lvara)
     else
        print*,'Incorrect or missing data field argument in shdf5_io'
        print*, 'field = ', dsetname
        stop    'shdf5_io: bad data field'
     endif

  else
     
     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif
  
end subroutine shdf5_io

!===============================================================================

subroutine shdf5_orec_ll(ndims,dims,dsetname,ivara,rvara,cvara,dvara,lvara,  &
                                             ivars,rvars,cvars,dvars,lvars,  &
                                             gpoints,                        &
                                             units, long_name, positive,     &
                                             imissing, rmissing, dmissing,   &
                                             isdim, dimnames, standard_name, &
                                             cell_methods                    )

  use oname_coms,  only: nl
  use consts_coms, only: r8
  use hdf5_f2f

  implicit none

  character(*), intent(in) :: dsetname ! Variable label
  integer,      intent(in) :: ndims    ! Number of dimensions or rank
  integer,      intent(in) :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer,      intent(in), optional :: ivara(*), ivars
  real,         intent(in), optional :: rvara(*), rvars
  character,    intent(in), optional :: cvara(*), cvars
  real(r8),     intent(in), optional :: dvara(*), dvars
  logical,      intent(in), optional :: lvara(*), lvars

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional :: gpoints(:)

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional :: dimnames(:)

! Local variables
  integer :: hdferr  ! Error flag

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec:', ndims, dims(1:ndims)
     stop    'shdf5_orec: bad dims'
  endif
     
! Prepare memory and options for the write

  call fh5_prepare_write_ll(ndims, dims, hdferr, nl%icompress, fcoords=gpoints)

  if (hdferr /= 0) then
     print*, "shdf5_orec: can't prepare requested field:", trim(dsetname)
     return
  endif

! Write the dataset.

      if (present(ivars)) then ; call fh5_write(ivars, dsetname, hdferr)
  elseif (present(rvars)) then ; call fh5_write(rvars, dsetname, hdferr)
  elseif (present(cvars)) then ; call fh5_write(cvars, dsetname, hdferr)
  elseif (present(dvars)) then ; call fh5_write(dvars, dsetname, hdferr)
  elseif (present(lvars)) then ; call fh5_write(lvars, dsetname, hdferr)
  elseif (present(ivara)) then ; call fh5_write(ivara, dsetname, hdferr)
  elseif (present(rvara)) then ; call fh5_write(rvara, dsetname, hdferr)
  elseif (present(cvara)) then ; call fh5_write(cvara, dsetname, hdferr)
  elseif (present(dvara)) then ; call fh5_write(dvara, dsetname, hdferr)
  elseif (present(lvara)) then ; call fh5_write(lvara, dsetname, hdferr)
  else
     print*, 'Incorrect or missing data field argument in shdf5_orec'
     stop    'shdf5_orec: bad data field'
  endif

  if (hdferr /= 0) then
     print*, 'In shdf5_orec: hdf5 write error =', hdferr
     stop    'shdf5_orec: hdf5 write error'
  endif

! Indicate if this variable is a dimension

  if (present(isdim)) then
     if (isdim) call fh5f_create_dim(dsetname, hdferr)
  endif

! Link each variable to its dimensions
  
  if ((present(dimnames)) .and. (.not. present(isdim))) then
     if (size(dimnames) >= ndims) call fh5f_attach_dims(ndims, dimnames, hdferr)
  endif

! Write any dataset attributes to match common NetCDF conventions

  if (present(units)) then
     if (len_trim(units) > 0) call fh5f_write_attribute("units", cvalue=units)
  endif

  if (present(long_name)) then
     if (len_trim(long_name) > 0) call fh5f_write_attribute("long_name", cvalue=long_name)
  endif

  if (present(standard_name)) then
     if (len_trim(standard_name) > 0) call fh5f_write_attribute("standard_name", cvalue=standard_name)
  endif

  if (present(cell_methods)) then
     if (len_trim(cell_methods) > 0) call fh5f_write_attribute("cell_methods", cvalue=cell_methods)
  endif

  if (present(positive)) then
     if (len_trim(positive) > 0) call fh5f_write_attribute("positive", cvalue=positive)
  endif

  if ((present(ivars) .or. present(ivara)) .and. present(imissing)) then
     call fh5f_write_attribute("missing_value", ivalue=imissing)
  endif

  if ((present(rvars) .or. present(rvara)) .and. present(rmissing)) then
     call fh5f_write_attribute("missing_value", rvalue=rmissing)
  endif

  if ((present(dvars) .or. present(dvara)) .and. present(dmissing)) then
     call fh5f_write_attribute("missing_value", dvalue=dmissing)
  endif

! Close the dataset, the dataspace for the dataset, and the dataspace properties.

  call fh5_close_write(hdferr)

end subroutine shdf5_orec_ll

end module
