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

module lite_vars
  implicit none

  real, allocatable :: ue(:,:)
  real, allocatable :: ve(:,:)

contains

!=========================================================================

subroutine calc_lite_vars()

  ! This routine is an example for computing extra derived variables
  ! that are to be output in lite files

  use mem_grid,    only: mza, xew, yew, zew, lpw, mwa
  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_basic,   only: vxe, vye, vze, rho
  use consts_coms, only: erad
  use misc_coms,   only: mdomain

  implicit none

  real    :: raxis
  integer :: j, iw, k, n

  ! Zonal and meridional wind speeds

  if (allocated(ue) .or. allocated(ve)) then

     if (mdomain < 2) then

        !$omp parallel do private(iw,k,raxis)
        do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

           if (allocated(ue)) ue(1:lpw(iw)-1,iw) = 0.
           if (allocated(ve)) ve(1:lpw(iw)-1,iw) = 0.

           raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

           if (raxis > 1.e3) then

              do k = lpw(iw), mza

                 if (allocated(ue)) then
                    ue(k,iw) = (vye(k,iw) * xew(iw) - vxe(k,iw) * yew(iw)) / raxis
                 endif

                 if (allocated(ve)) then
                    ve(k,iw) = vze(k,iw) * raxis / erad &
                             - (vxe(k,iw) * xew(iw) + vye(k,iw) * yew(iw)) &
                               * zew(iw) / (raxis * erad)
                 endif

              enddo

           else

              if (allocated(ue)) ue(lpw(iw):mza,iw) = 0.
              if (allocated(ve)) ve(lpw(iw):mza,iw) = 0.

           endif

        enddo
        !$omp end parallel do

     else
        if (allocated(ue)) ue = vxe
        if (allocated(ve)) ve = vye
     endif

  endif

end subroutine calc_lite_vars

!=========================================================================

subroutine prepare_lite()

  use var_tables, only: num_var, num_lite, vtab_r
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

  use max_dims,   only: maxlite
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

     case('UE')

        allocate(ue(mza,mwa))
        call increment_vtable('UE', 'AW', noread=.true., hist=.false., &
                              lite=.true., rvar2=ue)

     case('VE')

        allocate(ve(mza,mwa))
        call increment_vtable('VE', 'AW', noread=.true., hist=.false., &
                              lite=.true., rvar2=ve)

     end select

  enddo

end subroutine alloc_lite

!=========================================================================

subroutine lite_write()

  use oname_coms, only: nl
  use var_tables, only: num_var, vtab_r, get_vtab_dims, num_lite
  use misc_coms,  only: io6, ioutput, current_time, iclobber, &
                        iparallel
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close, mpi_does_parallel_io
  use max_dims,   only: pathlen
  use mem_grid,   only: nma, nua, nva, nwa, mza
  use mem_sfcg,   only: nwsfc
  use mem_land,   only: nland
  use mem_lake,   only: nlake
  use mem_sea,    only: nsea
  use mem_nudge,  only: nwnud
  use mem_para,   only: iva_globe_primary, iva_local_primary, mva_primary, &
                        iwa_globe_primary, iwa_local_primary, mwa_primary, &
                        ima_globe_primary, ima_local_primary, mma_primary, &
                        iwsfc_globe_primary, iwsfc_local_primary, mwsfc_primary, &
                        iland_globe_primary, iland_local_primary, mland_primary, &
                        ilake_globe_primary, ilake_local_primary, mlake_primary, &
                        isea_globe_primary, isea_local_primary, msea_primary, &
                        iwnud_globe_primary, iwnud_local_primary, mwnud_primary
  use hdf5_f2f,   only: fh5_close_cache

  implicit none

  ! This routine writes the chosen variables to the lite files

  character(pathlen) :: hnamel
  character(32)      :: varn
  character(2)       :: stagpt
  logical            :: exans
  integer            :: nv, nvcnt
  integer            :: ndims, idims(3)
  integer            :: nglobe, id, hdferr

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

  call shdf5_open(hnamel,'W',iclobber)

  ! Write the common fields

  call commio('WRITE')

  ! Store parallel I/O file mappings for repeated use

  if (mpi_does_parallel_io) then
     call cache_lite_writes()
  endif

  ! Loop through the main variable table and write those variables
  ! with the correct flag set

  nvcnt = 0
  do nv = 1, num_var

     if (vtab_r(nv)%ilite) then

        id = 1

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

!          if (ndims == 1) id = 2
           if (ndims == 2 .and. idims(1) == mza) id = 3

! Other sections of this IF block need "id" definitions??

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

        if     (associated(vtab_r(nv)%ivar0_p)) then
           call shdf5_orec(ndims, idims, varn, ivars=vtab_r(nv)%ivar0_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        elseif (associated(vtab_r(nv)%ivar1_p)) then
           call shdf5_orec(ndims, idims, varn, ivar1=vtab_r(nv)%ivar1_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        elseif (associated(vtab_r(nv)%ivar2_p)) then
           call shdf5_orec(ndims, idims, varn, ivar2=vtab_r(nv)%ivar2_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, cache_id=id)
        elseif (associated(vtab_r(nv)%ivar3_p)) then
           call shdf5_orec(ndims, idims, varn, ivar3=vtab_r(nv)%ivar3_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)

        elseif (associated(vtab_r(nv)%rvar0_p)) then
           call shdf5_orec(ndims, idims, varn, rvars=vtab_r(nv)%rvar0_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        elseif (associated(vtab_r(nv)%rvar1_p)) then
           call shdf5_orec(ndims, idims, varn, rvar1=vtab_r(nv)%rvar1_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        elseif (associated(vtab_r(nv)%rvar2_p)) then
           call shdf5_orec(ndims, idims, varn, rvar2=vtab_r(nv)%rvar2_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, cache_id=id)
        elseif (associated(vtab_r(nv)%rvar3_p)) then
           call shdf5_orec(ndims, idims, varn, rvar3=vtab_r(nv)%rvar3_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)

        elseif (associated(vtab_r(nv)%dvar0_p)) then
           call shdf5_orec(ndims, idims, varn, dvars=vtab_r(nv)%dvar0_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        elseif (associated(vtab_r(nv)%dvar1_p)) then
           call shdf5_orec(ndims, idims, varn, dvar1=vtab_r(nv)%dvar1_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        elseif (associated(vtab_r(nv)%dvar2_p)) then
           call shdf5_orec(ndims, idims, varn, dvar2=vtab_r(nv)%dvar2_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe, cache_id=id)
        elseif (associated(vtab_r(nv)%dvar3_p)) then
           call shdf5_orec(ndims, idims, varn, dvar3=vtab_r(nv)%dvar3_p, &
                           lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
        endif

        nvcnt = nvcnt + 1

     endif

  enddo

  if (mpi_does_parallel_io) then
!    call fh5_close_cache(2, hdferr)
     call fh5_close_cache(3, hdferr)
  endif

  call shdf5_close()

end subroutine lite_write



subroutine cache_lite_writes()

  use misc_coms,  only: ioutput
  use mem_grid,   only: mwa, nwa, mva, nva, mza
  use mem_para,   only: iwa_globe_primary, iwa_local_primary, mwa_primary
  use hdf5_f2f,   only: fh5_cache_write

  implicit none

  integer :: ndims, herr
  integer :: dims(3)

  dims = 0

!  ndims   = 1
!  dims(1) = mwa
!
!  call fh5_cache_write(ndims, dims, 2, herr, mcoords=iwa_local_primary, &
!                                 ifsize=nwa, fcoords=iwa_globe_primary)

  ndims   = 2
  dims(1) = mza
  dims(2) = mwa

  call fh5_cache_write(ndims, dims, 3, herr, mcoords=iwa_local_primary, &
                                 ifsize=nwa, fcoords=iwa_globe_primary)

! Need to add code for variables at other grid stagger points.

end subroutine cache_lite_writes

!=========================================================================

subroutine lite_read(litefile)

  use oname_coms,  only: nl
  use mem_grid,    only: nwa, nva, nma, mwa, mva, mma, mza
  use mem_ijtabs,  only: itab_w, itab_v, itab_m
  use misc_coms,   only: io6, runtype, time8, time_istp8
  use var_tables,  only: num_var, vtab_r, get_vtab_dims, num_lite
  use hdf5_utils,  only: shdf5_info, shdf5_irec, shdf5_open, shdf5_close
  use mem_sfcg,    only: itab_wsfc, nwsfc, mwsfc
  use mem_land,    only: itab_land, nland, mland, nzg
  use mem_lake,    only: itab_lake, nlake, mlake
  use mem_sea,     only: itab_sea, nsea, msea

  use mem_nudge,   only: nwnud, mwnud, itab_wnud

  implicit none

  character(*), intent(in) :: litefile

  logical       :: exans
  integer       :: nv, nvcnt, ns, ndims, idims(3)
  character(32) :: varn
  character (2) :: stagpt
  integer       :: ilocal(max(mwa,mva,mma,mwsfc,mland,mlake,msea,mwnud))

  inquire(file=trim(litefile),exist=exans)

  if (exans) then
  
! Lite file exists.  Open, read, and close file.

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening lite file ', trim(litefile)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(trim(litefile),'R')

     ! Read the common variables
     call commio('READ')

     nvcnt =  0

     ! Loop through all variables in the model vtables

     do nv = 1, num_var

  !      if (vtab_r(nv)%ilite) then
        if (vtab_r(nv)%name == 'WC') then

           ndims = -1
           idims =  0

           varn    = trim(vtab_r(nv)%name)
           stagpt  = vtab_r(nv)%stagpt

           call shdf5_info(varn, ndims, idims)

           ! Skip to next variable if the current one is not in the lite file

           if (ndims < 0) then
              write(io6,*)
              write(io6,*) 'Variable '//trim(varn)//' is not in the lite file, skipping'
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
              call shdf5_irec(ndims, idims, varn, ivar1=vtab_r(nv)%ivar1_p, points=ilocal(1:ns))
           elseif (associated(vtab_r(nv)%ivar2_p)) then
              call shdf5_irec(ndims, idims, varn, ivar2=vtab_r(nv)%ivar2_p, points=ilocal(1:ns))
           elseif (associated(vtab_r(nv)%ivar3_p)) then
              call shdf5_irec(ndims, idims, varn, ivar3=vtab_r(nv)%ivar3_p, points=ilocal(1:ns))

           elseif (associated(vtab_r(nv)%rvar1_p)) then
              call shdf5_irec(ndims, idims, varn, rvar1=vtab_r(nv)%rvar1_p, points=ilocal(1:ns))
           elseif (associated(vtab_r(nv)%rvar2_p)) then
              call shdf5_irec(ndims, idims, varn, rvar2=vtab_r(nv)%rvar2_p, points=ilocal(1:ns))
           elseif (associated(vtab_r(nv)%rvar3_p)) then
              call shdf5_irec(ndims, idims, varn, rvar3=vtab_r(nv)%rvar3_p, points=ilocal(1:ns))

           elseif (associated(vtab_r(nv)%dvar1_p)) then
              call shdf5_irec(ndims, idims, varn, dvar1=vtab_r(nv)%dvar1_p, points=ilocal(1:ns))
           elseif (associated(vtab_r(nv)%dvar2_p)) then
              call shdf5_irec(ndims, idims, varn, dvar2=vtab_r(nv)%dvar2_p, points=ilocal(1:ns))
           elseif (associated(vtab_r(nv)%dvar3_p)) then
              call shdf5_irec(ndims, idims, varn, dvar3=vtab_r(nv)%dvar3_p, points=ilocal(1:ns))

           endif

           nvcnt = nvcnt + 1
           write(io6, '(1x,a,I0,1x,a,3(1x,I0))')  &
                 'Read: ', nvcnt, trim(varn), idims(1:ndims)

        endif

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
