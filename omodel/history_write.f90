subroutine history_write(vtype)

  use var_tables, only: num_var, vtab_r, get_vtab_dims
  use misc_coms,  only: io6, ioutput, hfilepref, current_time, iclobber
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close, mpi_does_parallel_io
  use max_dims,   only: pathlen
  use mem_grid,   only: nma, nva, nwa
  use mem_sfcg,   only: nwsfc, nvsfc, nmsfc
  use mem_land,   only: nland
  use mem_lake,   only: nlake
  use mem_sea,    only: nsea
  use mem_nudge,  only: nwnud
  use hcane_rz,   only: htc0
  use mem_para,   only: iva_globe_primary, iva_local_primary, &
                        iwa_globe_primary, iwa_local_primary, &
                        ima_globe_primary, ima_local_primary, &
                        iwsfc_globe_primary, iwsfc_local_primary, &
                        ivsfc_globe_primary, ivsfc_local_primary, &
                        imsfc_globe_primary, imsfc_local_primary, &
                        iland_globe_primary, iland_local_primary, &
                        ilake_globe_primary, ilake_local_primary, &
                        isea_globe_primary,  isea_local_primary, &
                        iwnud_globe_primary, iwnud_local_primary
  implicit none

! This routine writes the chosen variables on the history file.

  character(*), intent(in) :: vtype  ! not used yet - 'state'

  character(pathlen) :: hnamel
  character(32)      :: varn
  character(2)       :: stagpt
  logical            :: exans
  integer            :: nv, nvcnt
  integer            :: ndims, idims(3)
  integer            :: nglobe

  integer, pointer, contiguous :: ilpts(:), igpts(:)

  ! Always write output for hurricane initialization even if ioutput = 0

  if (vtype(1:3) /= 'HTC') then
     if (ioutput == 0) return
  endif

! Construct h5 file name and open the file

  if (trim(vtype) == 'HTC0') then
     call makefnam(hnamel, hfilepref, current_time, 'HTC0', '$', 'h5')
     htc0 = hnamel
  elseif (trim(vtype) == 'HTC1') then
     call makefnam(hnamel, hfilepref, current_time, 'HTC1', '$', 'h5')
  else
     call makefnam(hnamel, hfilepref, current_time, 'H', '$', 'h5')
  endif

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

  call shdf5_open(hnamel,'W',iclobber,trypario=.true.)

! Write the common fields

  call commio('WRITE')

! Loop through the main variable table and write those variables
! with the correct flag set

  nvcnt = 0
  do nv = 1, num_var

     if (vtab_r(nv)%ihist) then

        varn = vtab_r(nv)%name
        call get_vtab_dims(nv, ndims, idims)

        write (io6, '(1x,a,2(I0,1x),a,3(1x,I0))')  &
             'Writing: ', nv, num_var, trim(varn), idims(1:ndims)

        stagpt = vtab_r(nv)%stagpt

        ! Identify which points will be written to disk

        if     (stagpt == 'AW') then
           ilpts => iwa_local_primary
           igpts => iwa_globe_primary
           nglobe = nwa
        elseif (stagpt == 'AV') then
           ilpts => iva_local_primary
           igpts => iva_globe_primary
           nglobe = nva
        elseif (stagpt == 'AM') then
           ilpts => ima_local_primary
           igpts => ima_globe_primary
           nglobe = nma
        elseif (stagpt == 'CW') then ! Common surface cells
           ilpts => iwsfc_local_primary
           igpts => iwsfc_globe_primary
           nglobe = nwsfc
        elseif (stagpt == 'CV') then ! Common surface cells (V pts)
           ilpts => ivsfc_local_primary
           igpts => ivsfc_globe_primary
           nglobe = nvsfc
        elseif (stagpt == 'CM') then ! Common surface cells (M pts)
           ilpts => imsfc_local_primary
           igpts => imsfc_globe_primary
           nglobe = nmsfc
        elseif (stagpt == 'LW') then ! Lake cells
           ilpts => iland_local_primary
           igpts => iland_globe_primary
           nglobe = nland
        elseif (stagpt == 'RW') then ! Lake cells
           ilpts => ilake_local_primary
           igpts => ilake_globe_primary
           nglobe = nlake
        elseif (stagpt == 'SW') then ! Sea cells
           ilpts => isea_local_primary
           igpts => isea_globe_primary
           nglobe = nsea
        elseif (stagpt == 'AN') then
           ilpts => iwnud_local_primary
           igpts => iwnud_globe_primary
           nglobe = nwnud
        else
           ! TODO: Const values
           ! TODO: Land U and M; Sea U and M? (probably not!)
           ! See history_start too
           stop "invalid array size in history_write for parallel output"
        endif

        if (idims(1) == 1 .and. ndims > 1) then
           idims(1:ndims-1) = idims(2:ndims)
           ndims = ndims - 1
        endif

        ! Now do the writes

        if     (associated(vtab_r(nv)%ivar1_p)) then
           call shdf5_orec(ndims, idims, varn, ivar1=vtab_r(nv)%ivar1_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        elseif (associated(vtab_r(nv)%ivar2_p)) then
           call shdf5_orec(ndims, idims, varn, ivar2=vtab_r(nv)%ivar2_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        elseif (associated(vtab_r(nv)%ivar3_p)) then
           call shdf5_orec(ndims, idims, varn, ivar3=vtab_r(nv)%ivar3_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)

        elseif (associated(vtab_r(nv)%rvar1_p)) then
           call shdf5_orec(ndims, idims, varn, rvar1=vtab_r(nv)%rvar1_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        elseif (associated(vtab_r(nv)%rvar2_p)) then
           call shdf5_orec(ndims, idims, varn, rvar2=vtab_r(nv)%rvar2_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        elseif (associated(vtab_r(nv)%rvar3_p)) then
           call shdf5_orec(ndims, idims, varn, rvar3=vtab_r(nv)%rvar3_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)

        elseif (associated(vtab_r(nv)%dvar1_p)) then
           call shdf5_orec(ndims, idims, varn, dvar1=vtab_r(nv)%dvar1_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        elseif (associated(vtab_r(nv)%dvar2_p)) then
           call shdf5_orec(ndims, idims, varn, dvar2=vtab_r(nv)%dvar2_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        elseif (associated(vtab_r(nv)%dvar3_p)) then
           call shdf5_orec(ndims, idims, varn, dvar3=vtab_r(nv)%dvar3_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, stagpt=stagpt)
        endif

        nvcnt = nvcnt + 1

     endif

  enddo

  call shdf5_close()

end subroutine history_write
