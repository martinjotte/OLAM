Module hdf5_utils

  use hdf5_f2f, FORTRAN_REAL4_TYPE => H5T_IEEE_F32LE, &
                FORTRAN_REAL8_TYPE => H5T_IEEE_F64LE, &
                FORTRAN_INT1_TYPE  => H5T_STD_I8LE,   &
                FORTRAN_INT2_TYPE  => H5T_STD_I16LE,  &
                FORTRAN_INT4_TYPE  => H5T_STD_I32LE

  implicit none

  integer, parameter :: r8 = selected_real_kind(13,300)
  integer, parameter :: i1 = selected_int_kind(2)

  private
  public :: shdf5_open, shdf5_info, shdf5_orec, shdf5_irec, shdf5_close, &
            shdf5_io, shdf5_write_global_attribute, &
            FORTRAN_REAL4_TYPE, FORTRAN_REAL8_TYPE, &
            FORTRAN_INT1_TYPE, FORTRAN_INT2_TYPE, FORTRAN_INT4_TYPE

Contains

!===============================================================================

subroutine shdf5_open(locfn, access, idelete)

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
  implicit none

  character(*), intent(in)    :: dsetname ! Dataset name
  integer,      intent(inout) :: dims(:)
  integer,      intent(inout) :: ndims    ! Dataset rank (in file)
  integer                     :: hdferr   ! Error flag

! Open the dataset.

  call fh5d_open(dsetname, hdferr)

  if (hdferr < 0) then
     ndims   = -1
     dims(1) =  0
     return
  endif

! Get dataset's dimensions

  call fh5s_get_ndims(ndims)
  call fh5s_get_dims(dims)

  call fh5d_close(hdferr)

end subroutine shdf5_info

!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars,  &
                                          bvar1,ivar1,rvar1,dvar1,lvar1,  &
                                          bvar2,ivar2,rvar2,dvar2,lvar2,  &
                                          bvar3,ivar3,rvar3,dvar3,lvar3,  &
                                          bvar4,ivar4,rvar4,dvar4,        &
                                          nglobe, lpoints, gpoints,       &
                                          units, long_name, positive,     &
                                          imissing, rmissing, dmissing,   &
                                          isdim, dimnames, standard_name, &
                                          cell_methods, storage_type,     &
                                          icompress, dims_chunk           )
  implicit none

  character(*), intent(in) :: dsetname    ! Variable label
  integer,      intent(in) :: ndims       ! Number of dimensions or rank
  integer,      intent(in) :: dims(ndims) ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional :: bvars
  integer,      intent(in), optional :: ivars
  real,         intent(in), optional :: rvars
  real(r8),     intent(in), optional :: dvars
  logical,      intent(in), optional :: lvars

  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional, contiguous :: lpoints(:), gpoints(:)
  integer,      intent(in), optional             :: nglobe

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional :: dimnames(ndims)

! Compression/chunking options
  integer,      intent(in), optional :: icompress         ! Compression level (0-9)
  integer,      intent(in), optional :: dims_chunk(ndims) ! Compression dimensions

! Type of variable
  integer(HID_T), intent(in), optional :: storage_type

! Local variables
  integer        :: hdferr  ! Error flag
  integer(HID_T) :: stype   ! Data storage kind

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec:', ndims, dims(1:ndims)
     stop    'shdf5_orec: bad dims'
  endif

! Prepare memory and options for the write

  if (present(storage_type)) then
     stype = storage_type
  elseif (present(bvars) .or. present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
     stype = FORTRAN_INT1_TYPE
  elseif (present(ivars) .or. present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
     stype = FORTRAN_INT4_TYPE
  elseif (present(rvars) .or. present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
     stype = FORTRAN_REAL4_TYPE
  elseif (present(dvars) .or. present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
     stype = FORTRAN_REAL8_TYPE
  elseif (present(lvars) .or. present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
     stype = FORTRAN_INT1_TYPE
  endif

  call fh5_prepare_write(ndims, dims, dsetname, stype, hdferr, icompress=icompress, &
       mcoords=lpoints, fcoords=gpoints, ifsize=nglobe, dims_chunk=dims_chunk)

  if (hdferr /= 0) then
     print*, "shdf5_orec: can't prepare requested field:", trim(dsetname)
     return
  endif

! Write the dataset.

      if (present(ivars)) then ; call fh5_write(ivars, hdferr)
  elseif (present(rvars)) then ; call fh5_write(rvars, hdferr)
  elseif (present(dvars)) then ; call fh5_write(dvars, hdferr)
  elseif (present(lvars)) then ; call fh5_write(lvars, hdferr)
  elseif (present(bvars)) then ; call fh5_write(bvars, hdferr)

  elseif (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
  elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
  elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
  elseif (present(lvar1)) then ; call fh5_write(lvar1, hdferr)
  elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

  elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
  elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
  elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
  elseif (present(lvar2)) then ; call fh5_write(lvar2, hdferr)
  elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

  elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
  elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
  elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
  elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)

  elseif (present(ivar4)) then ; call fh5_write(ivar4, hdferr)
  elseif (present(rvar4)) then ; call fh5_write(rvar4, hdferr)
  elseif (present(dvar4)) then ; call fh5_write(dvar4, hdferr)
  elseif (present(bvar4)) then ; call fh5_write(bvar4, hdferr)

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

  if (present(imissing)) then
     call fh5f_write_attribute("missing_value", ivalue=imissing)
  endif

  if (present(rmissing)) then
     call fh5f_write_attribute("missing_value", rvalue=rmissing)
  endif

  if (present(dmissing)) then
     call fh5f_write_attribute("missing_value", dvalue=dmissing)
  endif

! Close the dataset, the dataspace for the dataset, and the dataspace properties.

  call fh5_close_write(hdferr)

end subroutine shdf5_orec

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars,  &
                                          bvar1,ivar1,rvar1,dvar1,lvar1,  &
                                          bvar2,ivar2,rvar2,dvar2,lvar2,  &
                                          bvar3,ivar3,rvar3,dvar3,lvar3,  &
                                          bvar4,ivar4,rvar4,      dvar4,  &
                                          points, start, counts)
  implicit none

  character(*), intent(IN) :: dsetname ! Dataset name
  integer,      intent(IN) :: ndims    ! Number of dimensions or rank
  integer,      intent(IN) :: dims(:)  ! Dataset dimensions

! Array and scalar arguments for different types. Only specify one in each call.
  integer(i1),  intent(inout), optional :: bvars, bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(inout), optional :: ivars, ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(inout), optional :: rvars, rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  real(r8),     intent(inout), optional :: dvars, dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(inout), optional :: lvars, lvar1(:), lvar2(:,:), lvar3(:,:,:)

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

      if (present(ivars)) then ; call fh5_read(ivars, hdferr)
  elseif (present(rvars)) then ; call fh5_read(rvars, hdferr)
  elseif (present(dvars)) then ; call fh5_read(dvars, hdferr)
  elseif (present(lvars)) then ; call fh5_read(lvars, hdferr)
  elseif (present(bvars)) then ; call fh5_read(bvars, hdferr)

  elseif (present(ivar1)) then ; call fh5_read(ivar1, hdferr)
  elseif (present(rvar1)) then ; call fh5_read(rvar1, hdferr)
  elseif (present(dvar1)) then ; call fh5_read(dvar1, hdferr)
  elseif (present(lvar1)) then ; call fh5_read(lvar1, hdferr)
  elseif (present(bvar1)) then ; call fh5_read(bvar1, hdferr)

  elseif (present(ivar2)) then ; call fh5_read(ivar2, hdferr)
  elseif (present(rvar2)) then ; call fh5_read(rvar2, hdferr)
  elseif (present(dvar2)) then ; call fh5_read(dvar2, hdferr)
  elseif (present(lvar2)) then ; call fh5_read(lvar2, hdferr)
  elseif (present(bvar2)) then ; call fh5_read(bvar2, hdferr)

  elseif (present(ivar3)) then ; call fh5_read(ivar3, hdferr)
  elseif (present(rvar3)) then ; call fh5_read(rvar3, hdferr)
  elseif (present(dvar3)) then ; call fh5_read(dvar3, hdferr)
  elseif (present(lvar3)) then ; call fh5_read(lvar3, hdferr)
  elseif (present(bvar3)) then ; call fh5_read(bvar3, hdferr)

  elseif (present(ivar4)) then ; call fh5_read(ivar4, hdferr)
  elseif (present(rvar4)) then ; call fh5_read(rvar4, hdferr)
  elseif (present(dvar4)) then ; call fh5_read(dvar4, hdferr)
  elseif (present(bvar4)) then ; call fh5_read(bvar4, hdferr)

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
  implicit none

  integer :: hdferr  ! Error flags

! Close hdf file.

  call fh5f_close(hdferr)
end subroutine shdf5_close

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars, &
                                               bvar1,ivar1,rvar1,dvar1,lvar1, &
                                               bvar2,ivar2,rvar2,dvar2,lvar2, &
                                               bvar3,ivar3,rvar3,dvar3,lvar3, &
                                               bvar4,ivar4,rvar4,dvar4        )
  implicit none

  character(*), intent(in)           :: dsetname, action
  integer,      intent(in)           :: ndims,    dims(:)
  integer(i1),  intent(inout), optional :: bvars, bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(inout), optional :: ivars, ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(inout), optional :: rvars, rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  real(r8),     intent(inout), optional :: dvars, dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(inout), optional :: lvars, lvar1(:), lvar2(:,:), lvar3(:,:,:)

  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (action == 'READ') then

     call shdf5_irec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars, &
                                         bvar1,ivar1,rvar1,dvar1,lvar1, &
                                         bvar2,ivar2,rvar2,dvar2,lvar2, &
                                         bvar3,ivar3,rvar3,dvar3,lvar3, &
                                         bvar4,ivar4,rvar4,dvar4        )
  elseif (action == 'WRITE') then

     call shdf5_orec(ndims,dims,dsetname,bvars,ivars,rvars,dvars,lvars, &
                                         bvar1,ivar1,rvar1,dvar1,lvar1, &
                                         bvar2,ivar2,rvar2,dvar2,lvar2, &
                                         bvar3,ivar3,rvar3,dvar3,lvar3, &
                                         bvar4,ivar4,rvar4,dvar4        )
  else

     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif

end subroutine shdf5_io

subroutine shdf5_write_global_attribute(name, ivalue, rvalue, dvalue, cvalue)
  implicit none

  character(*),           intent(in) :: name
  integer,      optional, intent(in) :: ivalue
  real,         optional, intent(in) :: rvalue
  real(r8),     optional, intent(in) :: dvalue
  character(*), optional, intent(in) :: cvalue

  if (present(ivalue)) then
     call fh5f_write_global_attribute(name, ivalue=ivalue)
  elseif (present(rvalue)) then
     call fh5f_write_global_attribute(name, rvalue=rvalue)
  elseif (present(dvalue)) then
     call fh5f_write_global_attribute(name, dvalue=dvalue)
  elseif (present(cvalue)) then
     call fh5f_write_global_attribute(name, cvalue=cvalue)
  endif

end subroutine shdf5_write_global_attribute

end module
