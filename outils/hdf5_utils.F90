Module hdf5_utils

  use mem_para,    only: myrank
  use consts_coms, only: r8, i1
  use max_dims,    only: pathlen
  use hdf5_f2f

  character(1) :: prevaccess = ' '

  private
  public :: shdf5_open, shdf5_info, shdf5_orec, shdf5_irec, shdf5_close, &
            shdf5_io, shdf5_orec_ll, shdf5_write_global_attribute, &
            FORTRAN_INT1_TYPE, FORTRAN_INT2_TYPE, FORTRAN_INTEGER_TYPE, &
            FORTRAN_REAL_TYPE, FORTRAN_REAL8_TYPE, FORTRAN_LOGICAL_TYPE, &
            FORTRAN_CHARACTER_TYPE

Contains

!===============================================================================

subroutine shdf5_open(locfn, access, idelete, trypario)

  implicit none

  character(*),      intent(in) :: locfn   ! file name
  character(*),      intent(in) :: access  ! File access ('R','W','RW')
  integer, optional, intent(in) :: idelete ! If W, overwrite file if exists?
                                           !    1=yes, 0=no
  logical, optional, intent(in) :: trypario! Use parallel I/O if MPI run and
                                           ! compiled with PHDF5

  integer                       :: hdferr  ! Error flag
  integer                       :: iaccess ! int access flag
  logical                       :: exists  ! File existence

  prevaccess = access(1:1)

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)
  if (prevaccess == 'W' .and. myrank /= 0) return
#endif

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
        call fh5f_open(locfn, iaccess, hdferr, pario=trypario)

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
        call fh5f_create(locfn, iaccess, hdferr, pario=trypario)
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
           call fh5f_create(locfn, iaccess, hdferr, pario=trypario)
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

subroutine shdf5_info(dsetname, ndims, dims, dimname, attached_dimnames, units)

  implicit none

  character(*),        intent(in)    :: dsetname ! Dataset name
  integer, contiguous, intent(inout) :: dims(:)
  integer,             intent(inout) :: ndims    ! Dataset rank (in file)

! Optional arrays to read common NetCDF convention attributes or dimension information
  character(*), intent(inout), optional :: units
  character(*), intent(inout), optional :: dimname
  character(*), intent(inout), optional :: attached_dimnames(:)

! Open the dataset.

  call fh5_get_info(dsetname, ndims, dims)

! Return if there was an error

  if (ndims < 0) return

! Is this variable a dimension?

  if (present(dimname)) then
     call fh5f_query_dimname(dimname)
  endif

! Does the variable have dimension scales?

  if (present(attached_dimnames)) then
     attached_dimnames(:) = ' '
     call fh5_get_attached_scales(attached_dimnames)
  endif

! Read any attributes

  if (present(units)) then
     units = ' '
     call fh5f_read_attribute("units", cvalue=units)
  endif

! Close the dataset

  call fh5_close_info()

end subroutine shdf5_info

!===============================================================================

subroutine shdf5_orec(ndims,dims,dsetname,bvars,ivars,rvars,cvars,dvars,lvars, &
                                          bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                          bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                          bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                          bvar4,ivar4,rvar4,      dvar4,       &
                                          nglobe, lpoints, gpoints, stagpt,    &
                                          units, long_name, positive,          &
                                          imissing, rmissing, dmissing,        &
                                          isdim, dimnames, standard_name,      &
                                          cell_methods, storage_type           )

  use oname_coms,  only: nl
  use misc_coms,   only: iparallel

  implicit none

  character(*), intent(in)             :: dsetname ! Variable label
  integer,      intent(in)             :: ndims    ! Number of dimensions or rank
  integer,      intent(in), contiguous :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional :: bvars
  integer,      intent(in), optional :: ivars
  real,         intent(in), optional :: rvars
  character,    intent(in), optional :: cvars
  real(r8),     intent(in), optional :: dvars
  logical,      intent(in), optional :: lvars

  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  character,    intent(in), optional, contiguous :: cvar1(:), cvar2(:,:), cvar3(:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional             :: nglobe
  integer,      intent(in), optional, contiguous :: lpoints(:), gpoints(:)

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional, contiguous :: dimnames(:)

! Type of variable
  character(2),   intent(in), optional :: stagpt
  integer(HID_T), intent(in), optional :: storage_type

! Local variables
  integer        :: hdferr
  integer(HID_T) :: stype

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec:', ndims, dims(1:ndims)
     stop    'shdf5_orec: bad dims'
  endif

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)
  if ( (.not. present(bvars)) .and. (.not. present(ivars)) .and. (.not. present(rvars)) .and. &
       (.not. present(cvars)) .and. (.not. present(dvars)) .and. (.not. present(lvars)) .and. &
       (iparallel == 1) .and. present(gpoints) .and. present(nglobe) ) then

     call shdf5_orec2(ndims,dims,dsetname,bvar1,ivar1,rvar1,cvar1,dvar1,lvar1,  &
                                          bvar2,ivar2,rvar2,cvar2,dvar2,lvar2,  &
                                          bvar3,ivar3,rvar3,cvar3,dvar3,lvar3,  &
                                          bvar4,ivar4,rvar4,      dvar4,        &
                                          nglobe, lpoints, gpoints,       &
                                          units, long_name, positive,     &
                                          imissing, rmissing, dmissing,   &
                                          isdim, dimnames, standard_name, &
                                          cell_methods, storage_type      )
     return
  else
     if (myrank /= 0) return
  endif
#endif

  ! Prepare memory and options for the write

  if (present(storage_type)) then
     stype = storage_type
  elseif (present(bvars) .or. present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
     stype = FORTRAN_INT1_TYPE
  elseif (present(ivars) .or. present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
     stype = FORTRAN_INTEGER_TYPE
  elseif (present(rvars) .or. present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
     stype = FORTRAN_REAL_TYPE
  elseif (present(dvars) .or. present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
     stype = FORTRAN_REAL8_TYPE
  elseif (present(cvars) .or. present(cvar1) .or. present(cvar2) .or. present(cvar3)) then
     stype = FORTRAN_CHARACTER_TYPE
  elseif (present(lvars) .or. present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
     stype = FORTRAN_LOGICAL_TYPE
  endif

  call fh5_prepare_write(ndims, dims, dsetname, stype, hdferr, &
                         icompress=nl%icompress, type=stagpt, &
                         mcoords=lpoints, fcoords=gpoints, ifsize=nglobe)

  if (hdferr < 0) then
     print*, "shdf5_orec: can't prepare requested field:", trim(dsetname)
     stop
  endif

! Write the dataset.

      if (present(ivars)) then ; call fh5_write(ivars, hdferr)
  elseif (present(rvars)) then ; call fh5_write(rvars, hdferr)
  elseif (present(cvars)) then ; call fh5_write(cvars, hdferr)
  elseif (present(dvars)) then ; call fh5_write(dvars, hdferr)
  elseif (present(lvars)) then ; call fh5_write(lvars, hdferr)
  elseif (present(bvars)) then ; call fh5_write(bvars, hdferr)

  elseif (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
  elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
  elseif (present(cvar1)) then ; call fh5_write(cvar1, hdferr)
  elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
  elseif (present(lvar1)) then ; call fh5_write(lvar1, hdferr)
  elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

  elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
  elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
  elseif (present(cvar2)) then ; call fh5_write(cvar2, hdferr)
  elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
  elseif (present(lvar2)) then ; call fh5_write(lvar2, hdferr)
  elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

  elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
  elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
  elseif (present(cvar3)) then ; call fh5_write(cvar3, hdferr)
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

!================================================================

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)

subroutine shdf5_orec2(ndims,dims,dsetname,bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                           bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                           bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                           bvar4,ivar4,rvar4,      dvar4,       &
                                           nglobe, lpoints, gpoints,       &
                                           units, long_name, positive,     &
                                           imissing, rmissing, dmissing,   &
                                           isdim, dimnames, standard_name, &
                                           cell_methods, storage_type      )

  use oname_coms,  only: nl
  use mem_para,    only: mgroupsize, myrank
  use mpi

  implicit none

  character(*), intent(in)             :: dsetname ! Variable label
  integer,      intent(in)             :: ndims    ! Number of dimensions or rank
  integer,      intent(in), contiguous :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  character,    intent(in), optional, contiguous :: cvar1(:), cvar2(:,:), cvar3(:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional             :: nglobe
  integer,      intent(in), optional, contiguous :: lpoints(:), gpoints(:)

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional, contiguous :: dimnames(:)

! Type of variable
  integer(HID_T), intent(in), optional :: storage_type

! Local variables
  integer        :: hdferr  ! Error flag
  integer(HID_T) :: stype

  integer,     allocatable :: ibuff(:), points(:)
  real,        allocatable :: rbuff(:)
  real(r8),    allocatable :: dbuff(:)
  logical,     allocatable :: lbuff(:)
  integer(i1), allocatable :: bbuff(:)

  integer :: dimsn(5)
  integer :: nu, ier, base, n, i, is
  integer :: maxbuff, locbuff
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec2:', ndims, dims(1:ndims)
     stop    'shdf5_orec2: bad dims'
  endif

  nu = size(gpoints)
  call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

  base = maxval(nus(2:mgroupsize))
  maxbuff = base
  locbuff = nu
  is      = 1
  do n = ndims-1, 1, -1
     maxbuff = maxbuff * dims(n)
     locbuff = locbuff * dims(n)
     is      = is      * dims(n)
  enddo

  if (myrank == 0) then

     allocate(points(base))
     if     (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
        allocate(ibuff(maxbuff))
        stype = FORTRAN_INTEGER_TYPE
     elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
        allocate(rbuff(maxbuff))
        stype = FORTRAN_REAL_TYPE
     elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
        allocate(dbuff(maxbuff))
        stype = FORTRAN_REAL8_TYPE
     elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
        allocate(bbuff(maxbuff))
        stype = FORTRAN_INT1_TYPE
     elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)                    ) then
        allocate(lbuff(maxbuff))
        stype = FORTRAN_LOGICAL_TYPE
     endif

     if (present(storage_type)) then
        stype = storage_type
     endif

   endif

  if (myrank > 0 .and. nu > 0) then

     call MPI_Send(gpoints, nu, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, ier)

     if (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then

        allocate(ibuff(locbuff))
        do n = 1, nu
           if (present(lpoints)) then
              i = lpoints(n)
           else
              i = n
           endif

           if (present(ivar1)) ibuff( (n-1)*is+1 : n*is) =    ivar1(i)
           if (present(ivar2)) ibuff( (n-1)*is+1 : n*is) =    ivar2(:,i)
           if (present(ivar3)) ibuff( (n-1)*is+1 : n*is) = (/ ivar3(:,:,i)   /)
           if (present(ivar4)) ibuff( (n-1)*is+1 : n*is) = (/ ivar4(:,:,:,i) /)

        enddo
        call MPI_Send(ibuff, locbuff, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, ier)

     elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then

        allocate(rbuff(locbuff))
        do n = 1, nu
           if (present(lpoints)) then
              i = lpoints(n)
           else
              i = n
           endif

           if (present(rvar1)) rbuff( (n-1)*is+1 : n*is) =    rvar1(i)
           if (present(rvar2)) rbuff( (n-1)*is+1 : n*is) =    rvar2(:,i)
           if (present(rvar3)) rbuff( (n-1)*is+1 : n*is) = (/ rvar3(:,:,i)   /)
           if (present(rvar4)) rbuff( (n-1)*is+1 : n*is) = (/ rvar4(:,:,:,i) /)

        enddo
        call MPI_Send(rbuff, locbuff, MPI_REAL, 0, itag, MPI_COMM_WORLD, ier)

     elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then

        allocate(dbuff(locbuff))
        do n = 1, nu
           if (present(lpoints)) then
              i = lpoints(n)
           else
              i = n
           endif

           if (present(dvar1)) dbuff( (n-1)*is+1 : n*is) =    dvar1(i)
           if (present(dvar2)) dbuff( (n-1)*is+1 : n*is) =    dvar2(:,i)
           if (present(dvar3)) dbuff( (n-1)*is+1 : n*is) = (/ dvar3(:,:,i)   /)
           if (present(dvar4)) dbuff( (n-1)*is+1 : n*is) = (/ dvar4(:,:,:,i) /)

        enddo
        call MPI_Send(dbuff, locbuff, MPI_REAL8, 0, itag, MPI_COMM_WORLD, ier)

     elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then

        allocate(bbuff(locbuff))
        do n = 1, nu
           if (present(lpoints)) then
              i = lpoints(n)
           else
              i = n
           endif

           if (present(bvar1)) bbuff( (n-1)*is+1 : n*is) =    bvar1(i)
           if (present(bvar2)) bbuff( (n-1)*is+1 : n*is) =    bvar2(:,i)
           if (present(bvar3)) bbuff( (n-1)*is+1 : n*is) = (/ bvar3(:,:,i)   /)
           if (present(bvar4)) bbuff( (n-1)*is+1 : n*is) = (/ bvar4(:,:,:,i) /)

        enddo
        call MPI_Send(bbuff, locbuff, MPI_INTEGER1, 0, itag, MPI_COMM_WORLD, ier)

     elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then

        allocate(lbuff(locbuff))
        do n = 1, nu
           if (present(lpoints)) then
              i = lpoints(n)
           else
              i = n
           endif

            if (present(lvar1)) lbuff( (n-1)*is+1 : n*is) =    lvar1(i)
            if (present(lvar2)) lbuff( (n-1)*is+1 : n*is) =    lvar2(:,i)
            if (present(lvar2)) lbuff( (n-1)*is+1 : n*is) = (/ lvar3(:,:,i) /)

        enddo
        call MPI_Send(lbuff, locbuff, MPI_LOGICAL, 0, itag, MPI_COMM_WORLD, ier)

     endif

  endif

  if (myrank /= 0) return

  do n = 1, mgroupsize

     if (nus(n) < 1) cycle

     ! Prepare memory and options for the write

     if (n == 1) then
        call fh5_prepare_write(ndims, dims, dsetname, stype, hdferr, &
                               mcoords=lpoints, fcoords=gpoints, ifsize=nglobe)
     else
        call MPI_Recv(points, base, MPI_INTEGER, n-1, itag, MPI_COMM_WORLD, &
                      MPI_STATUS_IGNORE, ier)
        dimsn(1:ndims) = (/ dims(1:ndims-1), nus(n) /)
        call fh5_prepare_write(ndims, dimsn, dsetname, stype, hdferr, &
                               fcoords=points(1:nus(n)), ifsize=nglobe)
     endif

     if (hdferr /= 0) then
        print*, "shdf5_orec2: can't prepare requested field:", trim(dsetname)
        return
     endif

     ! Write the dataset.

     if (present(ivar1) .or. present(ivar2) .or. present(ivar3) .or. present(ivar4)) then
        if (n == 1) then
               if (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
           elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
           elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
           elseif (present(ivar4)) then ; call fh5_write(ivar4, hdferr)
           endif
        else
           call MPI_Recv( ibuff, maxbuff, MPI_INTEGER, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(ibuff, hdferr)
        endif
     elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3) .or. present(rvar4)) then
        if (n == 1) then
               if (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
           elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
           elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
           elseif (present(rvar4)) then ; call fh5_write(rvar4, hdferr)
           endif
        else
           call MPI_Recv( rbuff, maxbuff, MPI_REAL,    n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(rbuff, hdferr)
        endif
     elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3) .or. present(dvar4)) then
        if (n == 1) then
               if (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
           elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
           elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
           elseif (present(dvar4)) then ; call fh5_write(dvar4, hdferr)
           endif
        else
           call MPI_Recv( dbuff, maxbuff, MPI_REAL8,   n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(dbuff, hdferr)
        endif
     elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3) .or. present(bvar4)) then
        if (n == 1) then
               if (present(bvar1)) then ; call fh5_write(bvar1, hdferr)
           elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)
           elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)
           elseif (present(bvar4)) then ; call fh5_write(bvar4, hdferr)
           endif
        else
           call MPI_Recv( bbuff, maxbuff, MPI_INTEGER1, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(bbuff, hdferr)
        endif
     elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
        if (n == 1) then
               if (present(lvar1)) then ; call fh5_write(lvar1, hdferr)
           elseif (present(lvar2)) then ; call fh5_write(lvar2, hdferr)
           elseif (present(lvar3)) then ; call fh5_write(lvar3, hdferr)
           endif
        else
           call MPI_Recv( lbuff, maxbuff, MPI_LOGICAL, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(lbuff, hdferr)
        endif
     else
        print*, 'Incorrect or missing data field argument in shdf5_orec2'
        stop    'shdf5_orec2: bad data field'
     endif

     if (hdferr /= 0) then
        print*, 'In shdf5_orec2: hdf5 write error =', hdferr
        stop    'shdf5_orec2: hdf5 write error'
     endif

     call fh5_close_write1(hdferr)

  enddo

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

  call fh5_close_write2(hdferr)

end subroutine shdf5_orec2

#endif

!===============================================================================

subroutine shdf5_irec(ndims,dims,dsetname,bvars,ivars,rvars,cvars,dvars,lvars,  &
                                          bvar1,ivar1,rvar1,cvar1,dvar1,lvar1,  &
                                          bvar2,ivar2,rvar2,cvar2,dvar2,lvar2,  &
                                          bvar3,ivar3,rvar3,cvar3,dvar3,lvar3,  &
                                          bvar4,ivar4,rvar4,      dvar4,        &
                                          points, start, counts, stagpt,        &
                                          imissing, rmissing, dmissing,         &
                                          dimname, attached_dimnames, units,    &
                                          standard_name, long_name,             &
                                          rscale, roffset, dscale, doffset,     &
                                          have_missing, have_scaling)

  implicit none

  character(*), intent(IN)             :: dsetname ! Dataset name
  integer,      intent(IN)             :: ndims    ! Number of dimensions or rank
  integer,      intent(IN), contiguous :: dims(:)  ! Dataset dimensions

! Array and scalar arguments for different types. Only specify one in each call.
  integer(i1), intent(inout), optional :: bvars
  integer,     intent(inout), optional :: ivars
  real,        intent(inout), optional :: rvars
  character,   intent(inout), optional :: cvars
  real(r8),    intent(inout), optional :: dvars
  logical,     intent(inout), optional :: lvars

  integer(i1), intent(inout), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,     intent(inout), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,        intent(inout), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  character,   intent(inout), optional, contiguous :: cvar1(:), cvar2(:,:), cvar3(:,:,:)
  real(r8),    intent(inout), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,     intent(inout), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(IN), optional, contiguous :: points(:)
  integer,      intent(IN), optional, contiguous :: start (:)
  integer,      intent(IN), optional, contiguous :: counts(:)
  character(2), intent(IN), optional             :: stagpt

! Optional arrays to read common NetCDF convention attributes or dimension information
  character(*), intent(inout), optional :: units
  character(*), intent(inout), optional :: dimname
  character(*), intent(inout), optional :: attached_dimnames(:)
  character(*), intent(inout), optional :: standard_name, long_name
  integer,      intent(inout), optional :: imissing
  real,         intent(inout), optional :: rmissing, rscale, roffset
  real(r8),     intent(inout), optional :: dmissing, dscale, doffset
  logical,      intent(inout), optional :: have_missing, have_scaling

! Local variables
  integer :: hdferr  ! Error flag
  logical :: exists

! Check dimensions

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_irec:', ndims, dims(1:ndims)
     stop    'shdf5_irec: bad dims'
  endif

! Prepare file and memory space for the read

  call fh5_prepare_read(dsetname, ndims, dims, hdferr, coords=points, &
                        type=stagpt, start=start, counts=counts)
  if (hdferr < 0) then
     print*,'shdf5_irec: can''t prepare requested field:',trim(dsetname)
     return
  endif

! Read data from hyperslab in the file into the hyperslab in memory.

      if (present(ivars)) then ; call fh5_read(ivars, hdferr)
  elseif (present(rvars)) then ; call fh5_read(rvars, hdferr)
  elseif (present(cvars)) then ; call fh5_read(cvars, hdferr)
  elseif (present(dvars)) then ; call fh5_read(dvars, hdferr)
  elseif (present(lvars)) then ; call fh5_read(lvars, hdferr)
  elseif (present(bvars)) then ; call fh5_read(bvars, hdferr)

  elseif (present(ivar1)) then ; call fh5_read(ivar1, hdferr)
  elseif (present(rvar1)) then ; call fh5_read(rvar1, hdferr)
  elseif (present(cvar1)) then ; call fh5_read(cvar1, hdferr)
  elseif (present(dvar1)) then ; call fh5_read(dvar1, hdferr)
  elseif (present(lvar1)) then ; call fh5_read(lvar1, hdferr)
  elseif (present(bvar1)) then ; call fh5_read(bvar1, hdferr)

  elseif (present(ivar2)) then ; call fh5_read(ivar2, hdferr)
  elseif (present(rvar2)) then ; call fh5_read(rvar2, hdferr)
  elseif (present(cvar2)) then ; call fh5_read(cvar2, hdferr)
  elseif (present(dvar2)) then ; call fh5_read(dvar2, hdferr)
  elseif (present(lvar2)) then ; call fh5_read(lvar2, hdferr)
  elseif (present(bvar2)) then ; call fh5_read(bvar2, hdferr)

  elseif (present(ivar3)) then ; call fh5_read(ivar3, hdferr)
  elseif (present(rvar3)) then ; call fh5_read(rvar3, hdferr)
  elseif (present(cvar3)) then ; call fh5_read(cvar3, hdferr)
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

! Is this variable a dimension?

  if (present(dimname)) then
     call fh5f_query_dimname(dimname)
  endif

! Does the variable have dimension scales?

  if (present(attached_dimnames)) then
     attached_dimnames(:) = ' '
     call fh5_get_attached_scales(attached_dimnames)
  endif

! Read any attributes

  if (present(have_missing)) have_missing = .false.
  if (present(have_scaling)) have_scaling = .false.

  if (present(units)) then
     units = ' '
     call fh5f_read_attribute("units", cvalue=units)
  endif

  if (present(imissing)) then
     imissing = -huge(imissing)
     call fh5f_read_attribute("missing_value", ivalue=imissing, exists=exists)
     if (present(have_missing)) have_missing = exists
  endif

  if (present(rmissing)) then
     rmissing = -huge(rmissing)
     call fh5f_read_attribute("missing_value", rvalue=rmissing, exists=exists)
     if (present(have_missing)) have_missing = exists
  endif

  if (present(dmissing)) then
     dmissing = -huge(dmissing)
     call fh5f_read_attribute("missing_value", dvalue=dmissing, exists=exists)
     if (present(have_missing)) have_missing = exists
  endif

  if (present(rscale)) then
     rscale = 1.0
     call fh5f_read_attribute("scale_factor", rvalue=rscale, exists=exists)
     if (present(have_scaling) .and. exists) have_scaling = .true.
  endif

  if (present(dscale)) then
     dscale = 1.0_r8
     call fh5f_read_attribute("scale_factor", dvalue=dscale, exists=exists)
     if (present(have_scaling) .and. exists) have_scaling = .true.
  endif

  if (present(roffset)) then
     roffset = 0.0
     call fh5f_read_attribute("add_offset", rvalue=roffset, exists=exists)
     if (present(have_scaling) .and. exists) have_scaling = .true.
  endif

  if (present(doffset)) then
     doffset = 0.0_r8
     call fh5f_read_attribute("add_offset", dvalue=doffset, exists=exists)
     if (present(have_scaling) .and. exists) have_scaling = .true.
  endif

  if (present(long_name)) then
     long_name = ' '
     call fh5f_read_attribute("long_name", cvalue=long_name)
  endif

  if (present(standard_name)) then
     standard_name = ' '
     call fh5f_read_attribute("standard_name", cvalue=standard_name)
  endif

! Close the dataset, the dataspace for the dataset, and the memory space.

  call fh5_close_read(hdferr)

end subroutine shdf5_irec

!===============================================================================

subroutine shdf5_close()

  use misc_coms, only: io6

  implicit none

  integer :: hdferr  ! Error flags

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)
  if (prevaccess == 'W' .and. myrank /= 0) return
#endif

  ! Close hdf file.
  call fh5f_close(hdferr)

end subroutine shdf5_close

!===============================================================================

subroutine shdf5_io(action,ndims,dims,dsetname,bvars,ivars,rvars,cvars,dvars,lvars, &
                                               bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                               bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                               bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                               bvar4,ivar4,rvar4,      dvar4        )
  implicit none

  character(*), intent(in)             :: dsetname, action
  integer,      intent(in)             :: ndims
  integer,      intent(in), contiguous :: dims(:)

  integer(i1),  intent(inout), optional :: bvars
  integer,      intent(inout), optional :: ivars
  real,         intent(inout), optional :: rvars
  character,    intent(inout), optional :: cvars
  real(r8),     intent(inout), optional :: dvars
  logical,      intent(inout), optional :: lvars

  integer(i1),  intent(inout), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:), bvar4(:,:,:,:)
  integer,      intent(inout), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:), ivar4(:,:,:,:)
  real,         intent(inout), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:), rvar4(:,:,:,:)
  character,    intent(inout), optional, contiguous :: cvar1(:), cvar2(:,:), cvar3(:,:,:)
  real(r8),     intent(inout), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:), dvar4(:,:,:,:)
  logical,      intent(inout), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

  ! THIS ROUTINE CALLS SHDF5_IREC OR SHDF5_OREC TO READ OR WRITE A VARIABLE
  ! DEPENDING ON WHETHER 'ACTION' EQUALS 'READ' OR 'WRITE'

  if (action == 'READ') then

     call shdf5_irec(ndims,dims,dsetname,bvars,ivars,rvars,cvars,dvars,lvars, &
                                         bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                         bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                         bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                         bvar4,ivar4,rvar4,      dvar4        )
  elseif (action == 'WRITE') then

     call shdf5_orec(ndims,dims,dsetname,bvars,ivars,rvars,cvars,dvars,lvars, &
                                         bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                         bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                         bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                         bvar4,ivar4,rvar4,      dvar4        )
  else

     print *, "Illegal action in shdf5_io."
     print *, "Action should be 'READ' or 'WRITE'"
     stop     "Ending model run"

  endif

end subroutine shdf5_io

!===============================================================================

subroutine shdf5_orec_ll(ndims,dims,dsetname,bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                             bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                             bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                             gpoints, stagpt,                &
                                             units, long_name, positive,     &
                                             imissing, rmissing, dmissing,   &
                                             isdim, dimnames, standard_name, &
                                             cell_methods, storage_type      )

  use oname_coms,  only: nl
  use misc_coms,   only: iparallel, io6

  implicit none

  character(*), intent(in)             :: dsetname ! Variable label
  integer,      intent(in)             :: ndims    ! Number of dimensions or rank
  integer,      intent(in), contiguous :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:)
  character,    intent(in), optional, contiguous :: cvar1(:), cvar2(:,:), cvar3(:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional, contiguous :: gpoints(:)

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Type of variable
  character(2),   intent(in), optional :: stagpt
  integer(HID_T), intent(in), optional :: storage_type

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional, contiguous :: dimnames(:)

! Local variables
  integer        :: hdferr
  integer(HID_T) :: stype

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec_ll:', ndims, dims(1:ndims)
     stop    'shdf5_orec_ll: bad dims'
  endif

  write(io6,'(A,8(1x,I0))') " Writing: "//trim(dsetname), dims(1:ndims)

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)
  if (present(gpoints) .and. iparallel == 1) then

     call shdf5_orec_ll2(ndims,dims,dsetname,bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                             bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                             bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                             gpoints,                        &
                                             units, long_name, positive,     &
                                             imissing, rmissing, dmissing,   &
                                             isdim, dimnames, standard_name, &
                                             cell_methods, storage_type      )
     return
  else
     if (myrank /= 0) return
  endif
#endif

! Prepare memory and options for the write

  if (present(storage_type)) then
     stype = storage_type
  elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3)) then
     stype = FORTRAN_INT1_TYPE
  elseif (present(ivar1) .or. present(ivar2) .or. present(ivar3)) then
     stype = FORTRAN_INTEGER_TYPE
  elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3)) then
     stype = FORTRAN_REAL_TYPE
  elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3)) then
     stype = FORTRAN_REAL8_TYPE
  elseif (present(cvar1) .or. present(cvar2) .or. present(cvar3)) then
     stype = FORTRAN_CHARACTER_TYPE
  elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
     stype = FORTRAN_LOGICAL_TYPE
  endif

  call fh5_prepare_write_ll( ndims, dims, dsetname, stype, hdferr, &
                             icompress=nl%icompress, type=stagpt, fcoords=gpoints )

  if (hdferr < 0) then
     print*, "shdf5_orec_ll: can't prepare requested field:", trim(dsetname)
     return
  endif

! Write the dataset

      if (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
  elseif (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
  elseif (present(cvar1)) then ; call fh5_write(cvar1, hdferr)
  elseif (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
  elseif (present(lvar1)) then ; call fh5_write(lvar1, hdferr)
  elseif (present(bvar1)) then ; call fh5_write(bvar1, hdferr)

  elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
  elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
  elseif (present(cvar2)) then ; call fh5_write(cvar2, hdferr)
  elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
  elseif (present(lvar2)) then ; call fh5_write(lvar2, hdferr)
  elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)

  elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
  elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
  elseif (present(cvar3)) then ; call fh5_write(cvar3, hdferr)
  elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
  elseif (present(lvar3)) then ; call fh5_write(lvar3, hdferr)
  elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)
  else
     print*, 'Incorrect or missing data field argument in shdf5_orec_ll'
     stop    'shdf5_orec_ll: bad data field'
  endif

  if (hdferr /= 0) then
     print*, 'In shdf5_orec_ll: hdf5 write error =', hdferr
     stop    'shdf5_orec_ll: hdf5 write error'
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

end subroutine shdf5_orec_ll

!===============================================================================

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)

subroutine shdf5_orec_ll2(ndims,dims,dsetname,bvar1,ivar1,rvar1,cvar1,dvar1,lvar1, &
                                              bvar2,ivar2,rvar2,cvar2,dvar2,lvar2, &
                                              bvar3,ivar3,rvar3,cvar3,dvar3,lvar3, &
                                              gpoints,                        &
                                              units, long_name, positive,     &
                                              imissing, rmissing, dmissing,   &
                                              isdim, dimnames, standard_name, &
                                              cell_methods, storage_type      )

  use mem_para, only: mgroupsize, myrank
  use mpi

  implicit none

  character(*), intent(in)             :: dsetname ! Variable label
  integer,      intent(in)             :: ndims    ! Number of dimensions or rank
  integer,      intent(in), contiguous :: dims(:)  ! Dataset dimensions.

! Array and scalar arguments for different types. Only specify one in each call
  integer(i1),  intent(in), optional, contiguous :: bvar1(:), bvar2(:,:), bvar3(:,:,:)
  integer,      intent(in), optional, contiguous :: ivar1(:), ivar2(:,:), ivar3(:,:,:)
  real,         intent(in), optional, contiguous :: rvar1(:), rvar2(:,:), rvar3(:,:,:)
  character,    intent(in), optional, contiguous :: cvar1(:), cvar2(:,:), cvar3(:,:,:)
  real(r8),     intent(in), optional, contiguous :: dvar1(:), dvar2(:,:), dvar3(:,:,:)
  logical,      intent(in), optional, contiguous :: lvar1(:), lvar2(:,:), lvar3(:,:,:)

! Optional arrays to determine cells for partial/parallel IO
  integer,      intent(in), optional :: gpoints(:)

! Optional arrays to write common NetCDF convention attributes
  character(*), intent(in), optional :: units, long_name, positive
  character(*), intent(in), optional :: standard_name, cell_methods
  integer,      intent(in), optional :: imissing
  real,         intent(in), optional :: rmissing
  real(r8),     intent(in), optional :: dmissing

! Type of variable
  integer(HID_T), intent(in), optional :: storage_type

! Indicates if this variable is a global dimension
  logical,      intent(in), optional :: isdim

! Indicate names of each dimension
  character(*), intent(in), optional, contiguous :: dimnames(:)

! Local variables
  integer        :: hdferr  ! Error flag
  integer(HID_T) :: stype

  integer,     allocatable :: ibuff(:), points(:)
  real,        allocatable :: rbuff(:)
  real(r8),    allocatable :: dbuff(:)
  logical,     allocatable :: lbuff(:)
  integer(i1), allocatable :: bbuff(:)

  integer :: nu, ier, base, n
  integer :: maxbuff, locbuff
  integer :: nus(mgroupsize)
  integer, parameter :: itag = 40

! Check dimensions and set compression chunk size

  if (ndims <=0 .or. minval(dims(1:ndims)) <=0) then
     print*, 'Dimension error in shdf5_orec_ll2:', ndims, dims(1:ndims)
     stop    'shdf5_orec_ll2: bad dims'
  endif

  nu = size(gpoints)
  call MPI_Gather(nu, 1, MPI_INTEGER, nus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

  base = maxval(nus(2:mgroupsize))
  maxbuff = base
  locbuff = nu
  do n = 3, ndims
     maxbuff = maxbuff * dims(n)
     locbuff = locbuff * dims(n)
  enddo

  if (myrank == 0) then

     allocate(points(base))
     if     (present(ivar1) .or. present(ivar2) .or. present(ivar3)) then
        allocate(ibuff(maxbuff))
        stype = FORTRAN_INTEGER_TYPE
     elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3)) then
        allocate(rbuff(maxbuff))
        stype = FORTRAN_REAL_TYPE
     elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3)) then
        allocate(dbuff(maxbuff))
        stype = FORTRAN_REAL8_TYPE
     elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
        allocate(lbuff(maxbuff))
        stype = FORTRAN_LOGICAL_TYPE
     elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3)) then
        allocate(bbuff(maxbuff))
        stype = FORTRAN_INT1_TYPE
     endif

     if (present(storage_type)) then
        stype = storage_type
     endif

  endif

  if (myrank > 0 .and. nu > 0) then
     call MPI_Send(gpoints, nu, MPI_INTEGER, 0, itag, MPI_COMM_WORLD, ier)

     if     (present(ivar1)) then ; call MPI_Send(ivar1, locbuff, MPI_INTEGER,  0, itag, MPI_COMM_WORLD, ier)
     elseif (present(rvar1)) then ; call MPI_Send(rvar1, locbuff, MPI_REAL,     0, itag, MPI_COMM_WORLD, ier)
     elseif (present(dvar1)) then ; call MPI_Send(dvar1, locbuff, MPI_REAL8,    0, itag, MPI_COMM_WORLD, ier)
     elseif (present(lvar1)) then ; call MPI_Send(lvar1, locbuff, MPI_LOGICAL,  0, itag, MPI_COMM_WORLD, ier)
     elseif (present(bvar1)) then ; call MPI_Send(bvar1, locbuff, MPI_INTEGER1, 0, itag, MPI_COMM_WORLD, ier)

     elseif (present(ivar2)) then ; call MPI_Send(ivar2, locbuff, MPI_INTEGER,  0, itag, MPI_COMM_WORLD, ier)
     elseif (present(rvar2)) then ; call MPI_Send(rvar2, locbuff, MPI_REAL,     0, itag, MPI_COMM_WORLD, ier)
     elseif (present(dvar2)) then ; call MPI_Send(dvar2, locbuff, MPI_REAL8,    0, itag, MPI_COMM_WORLD, ier)
     elseif (present(lvar2)) then ; call MPI_Send(lvar2, locbuff, MPI_LOGICAL,  0, itag, MPI_COMM_WORLD, ier)
     elseif (present(bvar2)) then ; call MPI_Send(bvar2, locbuff, MPI_INTEGER1, 0, itag, MPI_COMM_WORLD, ier)

     elseif (present(ivar3)) then ; call MPI_Send(ivar3, locbuff, MPI_INTEGER,  0, itag, MPI_COMM_WORLD, ier)
     elseif (present(rvar3)) then ; call MPI_Send(rvar3, locbuff, MPI_REAL,     0, itag, MPI_COMM_WORLD, ier)
     elseif (present(dvar3)) then ; call MPI_Send(dvar3, locbuff, MPI_REAL8,    0, itag, MPI_COMM_WORLD, ier)
     elseif (present(lvar3)) then ; call MPI_Send(lvar3, locbuff, MPI_LOGICAL,  0, itag, MPI_COMM_WORLD, ier)
     elseif (present(bvar3)) then ; call MPI_Send(bvar3, locbuff, MPI_INTEGER1, 0, itag, MPI_COMM_WORLD, ier)

     endif
  endif

  if (myrank /= 0) return

  do n = 1, mgroupsize

     if (nus(n) < 1) cycle

! Prepare memory and options for the write

     if (n == 1) then
        call fh5_prepare_write_ll(ndims, dims, dsetname, stype, hdferr, fcoords=gpoints)
     else
        call MPI_Recv(points, base, MPI_INTEGER, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        call fh5_prepare_write_ll(ndims, dims, dsetname, stype, hdferr, fcoords=points(1:nus(n)))
     endif

     if (hdferr /= 0) then
        print*, "shdf5_orec_ll2: can't prepare requested field:", trim(dsetname)
        return
     endif

     ! Write the dataset.

     if (present(ivar1) .or. present(ivar2) .or. present(ivar3)) then
        if (n == 1) then
               if (present(ivar1)) then ; call fh5_write(ivar1, hdferr)
           elseif (present(ivar2)) then ; call fh5_write(ivar2, hdferr)
           elseif (present(ivar3)) then ; call fh5_write(ivar3, hdferr)
           endif
        else
           call MPI_Recv( ibuff, maxbuff, MPI_INTEGER, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(ibuff, hdferr)
        endif
     elseif (present(rvar1) .or. present(rvar2) .or. present(rvar3)) then
        if (n == 1) then
               if (present(rvar1)) then ; call fh5_write(rvar1, hdferr)
           elseif (present(rvar2)) then ; call fh5_write(rvar2, hdferr)
           elseif (present(rvar3)) then ; call fh5_write(rvar3, hdferr)
           endif
        else
           call MPI_Recv( rbuff, maxbuff, MPI_REAL,    n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(rbuff, hdferr)
        endif
     elseif (present(dvar1) .or. present(dvar2) .or. present(dvar3)) then
        if (n == 1) then
               if (present(dvar1)) then ; call fh5_write(dvar1, hdferr)
           elseif (present(dvar2)) then ; call fh5_write(dvar2, hdferr)
           elseif (present(dvar3)) then ; call fh5_write(dvar3, hdferr)
           endif
        else
           call MPI_Recv( dbuff, maxbuff, MPI_REAL8,   n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(dbuff, hdferr)
        endif
     elseif (present(lvar1) .or. present(lvar2) .or. present(lvar3)) then
        if (n == 1) then
               if (present(lvar1)) then ; call fh5_write(lvar1, hdferr)
           elseif (present(lvar2)) then ; call fh5_write(lvar2, hdferr)
           elseif (present(lvar3)) then ; call fh5_write(lvar3, hdferr)
           endif
       else
           call MPI_Recv( lbuff, maxbuff, MPI_LOGICAL, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(lbuff, hdferr)
        endif
     elseif (present(bvar1) .or. present(bvar2) .or. present(bvar3)) then
        if (n == 1) then
               if (present(bvar1)) then ; call fh5_write(bvar1, hdferr)
           elseif (present(bvar2)) then ; call fh5_write(bvar2, hdferr)
           elseif (present(bvar3)) then ; call fh5_write(bvar3, hdferr)
           endif
       else
           call MPI_Recv( bbuff, maxbuff, MPI_INTEGER1, n-1, itag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
           call fh5_write(bbuff, hdferr)
        endif
     else
        print*, 'Incorrect or missing data field argument in shdf5_orec_ll2'
        stop    'shdf5_orec_ll2: bad data field'
     endif

     if (hdferr /= 0) then
        print*, 'In shdf5_orec_ll2: hdf5 write error =', hdferr, myrank
        stop    'shdf5_orec_ll2: hdf5 write error'
     endif

     call fh5_close_write1(hdferr)

  enddo

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

  call fh5_close_write2(hdferr)

end subroutine shdf5_orec_ll2

#endif

!===============================================================================

subroutine shdf5_write_global_attribute(name, ivalue, rvalue, dvalue, cvalue)

  implicit none

  character(*),           intent(in) :: name
  integer,      optional, intent(in) :: ivalue
  real,         optional, intent(in) :: rvalue
  real(r8),     optional, intent(in) :: dvalue
  character(*), optional, intent(in) :: cvalue

#if defined(OLAM_MPI) && !defined(OLAM_PARALLEL_HDF5)
  if (myrank /= 0) return
#endif

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

!===============================================================================

end module hdf5_utils
