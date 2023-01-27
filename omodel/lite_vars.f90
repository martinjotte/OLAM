module lite_vars

  use max_dims,   only: maxlite
  use var_tables, only: num_lite

  implicit none

  integer :: lite_map(maxlite)
  private :: maxlite, num_lite

contains

!=========================================================================

subroutine calc_lite_vars()

  ! This routine is a placeholder for computing extra derived variables
  ! to be output in the lite files

! use mem_grid,    only: mza, xew, yew, zew, lpw, mwa
! use mem_ijtabs,  only: jtab_w, jtw_prog
! use mem_basic,   only: vxe, vye, vze, rho
! use consts_coms, only: erad
! use misc_coms,   only: mdomain

  implicit none

end subroutine calc_lite_vars

!=========================================================================

subroutine prepare_lite()

  use var_tables, only: num_var, vtab_r
  use oname_coms, only: nl
  use max_dims,   only: maxlite
  use misc_coms,  only: io6

  implicit none

  integer :: n, j

  num_lite = 0

  if (nl%ioutput_lite /= 1) return

  ! allocate memory for any derived lite variables and add to variable tables

  call alloc_lite()

  ! loop through variable tables and flag variables for lite output

  do n = 1, maxlite

     if (nl%lite_vars(n) == '') exit

     do j = 1, num_var

        if (nl%lite_vars(n) == vtab_r(j)%name) then
           vtab_r(j)%ilite = .true.
           num_lite = num_lite + 1
           lite_map(num_lite) = j
           exit
        endif

        if (j == num_var) then
           write(io6,*) "No match for lite variable " // trim(nl%lite_vars(n)) // ", skipping."
        endif

     enddo
  enddo

end subroutine prepare_lite

!=========================================================================

subroutine alloc_lite()

  use oname_coms, only: nl
  use var_tables, only: increment_vtable
  use mem_grid,   only: mza, mwa

  implicit none

  integer :: n

  ! This routine allocates space and sets the variable tables for any extra
  ! derived quantities that we may want to output in lite files.

  do n = 1, maxlite

     if (nl%lite_vars(n) == '') exit

     select case ( nl%lite_vars(n) )

!     case('UE')
!
!        allocate(ue(mza,mwa))
!        call increment_vtable('UE', 'AW', noread=.true., hist=.false., &
!                              lite=.true., rvar2=ue)
!
!     case('VE')
!
!        allocate(ve(mza,mwa))
!        call increment_vtable('VE', 'AW', noread=.true., hist=.false., &
!                              lite=.true., rvar2=ve)
!
     end select

  enddo

end subroutine alloc_lite

!=========================================================================

subroutine lite_write()

  use oname_coms, only: nl
  use var_tables, only: num_var, vtab_r, get_vtab_dims
  use misc_coms,  only: io6, current_time, iclobber
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close, mpi_does_parallel_io
  use max_dims,   only: pathlen
  use mem_grid,   only: nma, nva, nwa
  use mem_sfcg,   only: nwsfc, nvsfc, nmsfc
  use mem_land,   only: nland
  use mem_lake,   only: nlake
  use mem_sea,    only: nsea
  use mem_nudge,  only: nwnud
  use mem_para,   only: iva_globe_primary,   iva_local_primary,   &
                        iwa_globe_primary,   iwa_local_primary,   &
                        ima_globe_primary,   ima_local_primary,   &
                        iwsfc_globe_primary, iwsfc_local_primary, &
                        ivsfc_globe_primary, ivsfc_local_primary, &
                        imsfc_globe_primary, imsfc_local_primary, &
                        iland_globe_primary, iland_local_primary, &
                        ilake_globe_primary, ilake_local_primary, &
                        isea_globe_primary,  isea_local_primary,  &
                        iwnud_globe_primary, iwnud_local_primary
  implicit none

  ! This routine writes the chosen variables to the lite files

  character(pathlen) :: hnamel
  character(32)      :: varn
  character(2)       :: stagpt
  logical            :: exans
  integer            :: nv, j
  integer            :: ndims, idims(3)
  integer            :: nglobe

  integer, pointer, contiguous :: ilpts(:), igpts(:)

  if (nl%ioutput_lite == 0 .or. num_lite < 1) return

  ! compute any derived lite variables

  call calc_lite_vars()

  ! Construct h5 file name and open the file

  call makefnam(hnamel,  nl%hfilepref, current_time, 'L', '$', 'h5')

  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open file name :'
     write(io6,*) '!!!       '//trim(hnamel)
     write(io6,*) '!!!   but it already exists. run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'lite_write'
  endif

  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(io6,*) 'lite_write: opening file: ',trim(hnamel)
  write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  call shdf5_open(hnamel,'W',iclobber,trypario=.true.)

  ! Write the common fields

  call commio('WRITE')

  ! Loop through the main variable table and write those variables
  ! with the correct flag set

  do j = 1, num_lite
     nv = lite_map(j)

     varn = vtab_r(nv)%name
     call get_vtab_dims(nv, ndims, idims)

     write (io6, '(1x,a,2(I0,1x),a,3(1x,I0))')  &
          'Writing: ', nv, num_var, trim(varn), idims(1:ndims)

     stagpt = vtab_r(nv)%stagpt

     ! Choose points to be written

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
     elseif (stagpt == 'CW') then
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
     elseif (stagpt == 'LW') then
        ilpts => iland_local_primary
        igpts => iland_globe_primary
        nglobe = nland
     elseif (stagpt == 'RW') then
        ilpts => ilake_local_primary
        igpts => ilake_globe_primary
        nglobe = nlake
     elseif (stagpt == 'SW') then
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
        stop "invalid array size in lite_write for parallel output"
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

     elseif (associated(vtab_r(nv)%rvar0_p)) then
        call shdf5_orec(ndims, idims, varn, rvars=vtab_r(nv)%rvar0_p, &
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

     elseif (associated(vtab_r(nv)%dvar0_p)) then
        call shdf5_orec(ndims, idims, varn, dvars=vtab_r(nv)%dvar0_p, &
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

  enddo

  call shdf5_close()

end subroutine lite_write

!=========================================================================

subroutine lite_read(litefile)

  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma
  use mem_ijtabs,  only: itab_w, itab_v, itab_m
  use misc_coms,   only: io6, time8, time_istp8
  use var_tables,  only: vtab_r, get_vtab_dims, num_lite
  use hdf5_utils,  only: shdf5_info, shdf5_irec, shdf5_open, shdf5_close
  use mem_sfcg,    only: itab_wsfc, nwsfc, mwsfc, itab_vsfc, nvsfc, mvsfc
  use mem_land,    only: itab_land, nland, mland
  use mem_lake,    only: itab_lake, nlake, mlake
  use mem_sea,     only: itab_sea, nsea, msea

  use mem_nudge,   only: nwnud, mwnud, itab_wnud

  implicit none

  character(*), intent(in) :: litefile

  logical       :: exans
  integer       :: j, nv, ns, ndims, idims(3)
  character(32) :: varn
  character (2) :: stagpt
  integer       :: ilocal(max(mwa,mva,mma,mwsfc,mvsfc,mland,mlake,msea,mwnud))

  inquire(file=trim(litefile),exist=exans)

  if (exans) then

! Lite file exists.  Open, read, and close file.

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening lite file ', trim(litefile)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(trim(litefile),'R',trypario=.true.)

     ! Read the common variables
     call commio('READ')

     ! Loop through all variables in the model vtables

     do j = 1, num_lite
        nv = lite_map(j)

        ndims = -1
        idims =  0

        varn    = trim(vtab_r(nv)%name)
        stagpt  = vtab_r(nv)%stagpt

        call shdf5_info(varn, ndims, idims)

        ! Skip to next variable if the current one is not in the lite file

        if (ndims < 0) then
           write(io6,*)
           write(io6,*) 'Variable '//trim(varn)//' is not in the input file, skipping'
           write(io6,*)
           cycle
        endif

        ! Identify the points we want to read from the history file

        if     (stagpt == 'AW' .and. idims(ndims) == nwa) then
           ilocal(1:mwa) = itab_w(1:mwa)%iwglobe
           idims(ndims) = mwa
        elseif (stagpt == 'AV' .and. idims(ndims) == nva) then
           ilocal(1:mva) = itab_v(1:mva)%ivglobe
           idims(ndims) = mva
        elseif (stagpt == 'AM' .and. idims(ndims) == nma) then
           ilocal(1:mma) = itab_m(1:mma)%imglobe
           idims(ndims) = mma
        elseif (stagpt == 'CW' .and. idims(ndims) == nwsfc) then
           ilocal(1:mwsfc) = itab_wsfc(1:mwsfc)%iwglobe
           idims(ndims) = mwsfc
        elseif (stagpt == 'CV' .and. idims(ndims) == nvsfc) then
           ilocal(1:mvsfc) = itab_vsfc(1:mvsfc)%ivglobe
           idims(ndims) = mvsfc
        elseif (stagpt == 'LW' .and. idims(ndims) == nland) then
           ilocal(1:mland) = itab_land(1:mland)%iwglobe
           idims(ndims) = mland
        elseif (stagpt == 'RW' .and. idims(ndims) == nlake) then
           ilocal(1:mlake) = itab_lake(1:mlake)%iwglobe
           idims(ndims) = mlake
        elseif (stagpt == 'SW' .and. idims(ndims) == nsea) then
           ilocal(1:msea) = itab_sea(1:msea)%iwglobe
           idims(ndims) = msea
        elseif (stagpt == 'AN' .and. idims(ndims) == nwnud) then
           ilocal(1:mwnud) = itab_wnud(1:mwnud)%iwnudglobe
           idims(ndims) = mwnud
        else

           ! TODO: Const values
           ! TODO: Land U and M; Sea U and M? (probably not!)

           stop "invalid array size in history_read"
        endif

        ns = idims(ndims)

        if     (associated(vtab_r(nv)%ivar1_p)) then
           call shdf5_irec(ndims, idims, varn, ivar1=vtab_r(nv)%ivar1_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        elseif (associated(vtab_r(nv)%ivar2_p)) then
           call shdf5_irec(ndims, idims, varn, ivar2=vtab_r(nv)%ivar2_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        elseif (associated(vtab_r(nv)%ivar3_p)) then
           call shdf5_irec(ndims, idims, varn, ivar3=vtab_r(nv)%ivar3_p, &
                           points=ilocal(1:ns), stagpt=stagpt)

        elseif (associated(vtab_r(nv)%rvar1_p)) then
           call shdf5_irec(ndims, idims, varn, rvar1=vtab_r(nv)%rvar1_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        elseif (associated(vtab_r(nv)%rvar2_p)) then
           call shdf5_irec(ndims, idims, varn, rvar2=vtab_r(nv)%rvar2_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        elseif (associated(vtab_r(nv)%rvar3_p)) then
           call shdf5_irec(ndims, idims, varn, rvar3=vtab_r(nv)%rvar3_p, &
                           points=ilocal(1:ns), stagpt=stagpt)

        elseif (associated(vtab_r(nv)%dvar1_p)) then
           call shdf5_irec(ndims, idims, varn, dvar1=vtab_r(nv)%dvar1_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        elseif (associated(vtab_r(nv)%dvar2_p)) then
           call shdf5_irec(ndims, idims, varn, dvar2=vtab_r(nv)%dvar2_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        elseif (associated(vtab_r(nv)%dvar3_p)) then
           call shdf5_irec(ndims, idims, varn, dvar3=vtab_r(nv)%dvar3_p, &
                           points=ilocal(1:ns), stagpt=stagpt)
        endif

        write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
             'Read: ', nv, trim(varn), idims(1:ndims)

     enddo

     time_istp8 = time8  ! time_istp8 is used for plotting time
     call shdf5_close()

  else

   ! Lite file does not exist, stop model.

     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open lite file:'
     write(io6,*) '!!!   '//trim(litefile)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'in lite_vars'

  endif

end subroutine lite_read

end module lite_vars
