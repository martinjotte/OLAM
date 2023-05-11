subroutine commio(action)

! THIS ROUTINE READS OR WRITES NON-HORIZONTALLY VARYING FIELDS
! THAT ARE COMMON TO ALL PROCESSES. CALL THIS ROUTINE WITH
! 'ACTION' EQUALS 'READ' OR 'WRITE' TO INPUT OR OUTPUT THE FIELDS

  use max_dims,   only: maxngrdll
  use misc_coms,  only: io6, itime1, idate1, imonth1, iyear1, nxp, ngrids, &
                        ngrdll, grdrad, grdlat, grdlon, nzp, &
                        mdomain, deltax, runtype, &
                        itopoflg, time8, ndz, hdz, dz, current_time
  use leaf_coms,  only: nzs, ivegflg, isfcl
  use mem_land,   only: nzg, landgrid_dztop, landgrid_depth
  use hdf5_utils, only: shdf5_orec, shdf5_irec, shdf5_io
  use hcane_rz,   only: hlat_hist, hlon_hist, icycle_hurrinit_hist

  implicit none
  integer          :: ndims, idims(2)
  character(len=*) :: action

  if (action /= "READ" .and. action /= "WRITE") then
     write(io6,*) 'Illegal action in routine commio'
     write(io6,*) 'action must be "READ" or "WRITE"'
     stop     'Stopping run'
  endif

! The following are namelist variables that are read from a history file
! anytime subroutine commio is called, which occurs for runtype = 'HISTORY',
! 'HISTREGRID', or 'PLOTONLY'.  This overrides values copied from the namelist,
! ensuring against inadvertent changes in their namelist values, which is not
! compatible with history restarts.

  ndims = 1
  idims(1) = 1
  call shdf5_io(action, ndims, idims, 'nl%itime1',   ivars=itime1)
  call shdf5_io(action, ndims, idims, 'nl%idate1' ,  ivars=idate1)
  call shdf5_io(action, ndims, idims, 'nl%imonth1',  ivars=imonth1)
  call shdf5_io(action, ndims, idims, 'nl%iyear1',   ivars=iyear1)
  call shdf5_io(action, ndims, idims, 'nl%nzp',      ivars=nzp)
  call shdf5_io(action, ndims, idims, 'nl%nzg',      ivars=nzg)
  call shdf5_io(action, ndims, idims, 'nl%nzs',      ivars=nzs)
  call shdf5_io(action, ndims, idims, 'nl%mdomain',  ivars=mdomain)
  call shdf5_io(action, ndims, idims, 'nl%isfcl',    ivars=isfcl)
  call shdf5_io(action, ndims, idims, 'nl%itopoflg', ivars=itopoflg)
  call shdf5_io(action, ndims, idims, 'nl%ivegflg',  ivars=ivegflg)
  call shdf5_io(action, ndims, idims, 'nl%ndz',      ivars=ndz)

  call shdf5_io(action, ndims, idims, 'nl%landgrid_dztop', rvars=landgrid_dztop)
  call shdf5_io(action, ndims, idims, 'nl%landgrid_depth', rvars=landgrid_depth)

  ndims = 1
  idims(1) = ndz
  call shdf5_io(action, ndims, idims, 'nl%hdz',    rvar1=hdz)
  call shdf5_io(action, ndims, idims, 'nl%dz',     rvar1=dz)

! The following are not namelist variables but are read from a history file
! anytime that subroutine commio is called.

  ndims=1
  idims(1) = 1
  call shdf5_io(action, ndims, idims, 'time8'               , dvars=time8)
  call shdf5_io(action, ndims, idims, 'cur%year'            , ivars=current_time%year)
  call shdf5_io(action, ndims, idims, 'cur%month'           , ivars=current_time%month)
  call shdf5_io(action, ndims, idims, 'cur%date'            , ivars=current_time%date)
  call shdf5_io(action, ndims, idims, 'cur%time'            , dvars=current_time%time)
  call shdf5_io(action, ndims, idims, 'hlat_hist'           , rvars=hlat_hist)
  call shdf5_io(action, ndims, idims, 'hlon_hist'           , rvars=hlon_hist)
  call shdf5_io(action, ndims, idims, 'icycle_hurrinit_hist', ivars=icycle_hurrinit_hist)

  ! The following are namelist variables that specify the horizontal structure
  ! of the ATM grid, and they are read from a history file (to maintain their
  ! value from previous model runs) unless the horizontal grid structure is
  ! being changed, in which case runtype = 'HISTREGRID'.

  if (action == 'WRITE' .or. &
      runtype == 'HISTORY' .or. runtype == 'PLOTONLY') then

     ndims=1
     idims(1) = 1
     call shdf5_io(action, ndims, idims, 'nl%nxp',      ivars=nxp)
     call shdf5_io(action, ndims, idims, 'nl%deltax',   rvars=deltax)
     call shdf5_io(action, ndims, idims, 'nl%ngrids',  ivars=ngrids)

     ndims = 1
     idims(1) = ngrids
     call shdf5_io(action, ndims, idims, 'nl%ngrdll', ivar1=ngrdll)

     ndims = 2
     idims(1) = ngrids
     idims(2) = maxngrdll

     call shdf5_io(action, ndims, idims, 'nl%grdrad', rvar2=grdrad(1:ngrids,:))
     call shdf5_io(action, ndims, idims, 'nl%grdlat', rvar2=grdlat(1:ngrids,:))
     call shdf5_io(action, ndims, idims, 'nl%grdlon', rvar2=grdlon(1:ngrids,:))

  endif

end subroutine commio
