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

  real, allocatable, target :: ue(:,:)
  real, allocatable, target :: ve(:,:)

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
  integer :: j, iw, k, n, isf, ilf

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
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use max_dims,   only: pathlen
  use mem_grid,   only: nma, nua, nva, nwa
  use leaf_coms,  only: nwl
  use sea_coms,   only: nws
  use mem_nudge,  only: nwnud
  use mem_para,   only: iva_globe_primary, iva_local_primary, mva_primary, &
                        iwa_globe_primary, iwa_local_primary, mwa_primary, &
                        ima_globe_primary, ima_local_primary, mma_primary, &
                        iwl_globe_primary, iwl_local_primary, mwl_primary, &
                        iws_globe_primary, iws_local_primary, mws_primary, &
                        iwnud_globe_primary, iwnud_local_primary, mwnud_primary
  implicit none

  ! This routine writes the chosen variables to the lite files

  character(pathlen) :: hnamel
  character(32)      :: varn
  character(2)       :: stagpt
  logical            :: exans
  integer            :: nv, nvcnt
  integer            :: ndims, idims(3)
  integer, pointer   :: ilpts(:), igpts(:)
  integer            :: nglobe

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

  ! Loop through the main variable table and write those variables
  ! with the correct flag set

  nvcnt = 0
  do nv = 1, num_var

     if (vtab_r(nv)%ilite) then

        varn = vtab_r(nv)%name
        call get_vtab_dims(nv, ndims, idims)

        write (io6, '(1x,a,2(I0,1x),a,3(1x,I0))')  &
             'Writing: ', nv, num_var, trim(varn), idims(1:ndims)

        if (iparallel == 1) then

           stagpt = vtab_r(nv)%stagpt

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
           elseif (stagpt == 'LW') then
              ilpts => iwl_local_primary
              igpts => iwl_globe_primary
              nglobe = nwl
           elseif (stagpt == 'SW') then
              ilpts => iws_local_primary
              igpts => iws_globe_primary
              nglobe = nws
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

           if     (associated(vtab_r(nv)%ivar0_p)) then
              call shdf5_orec(ndims, idims, trim(varn), ivars=vtab_r(nv)%ivar0_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%ivar1_p)) then
              call shdf5_orec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar1_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%ivar2_p)) then
              call shdf5_orec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar2_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%ivar3_p)) then
              call shdf5_orec(ndims, idims, trim(varn), ivara=vtab_r(nv)%ivar3_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
              
           elseif (associated(vtab_r(nv)%rvar0_p)) then
              call shdf5_orec(ndims, idims, trim(varn), rvars=vtab_r(nv)%rvar0_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%rvar1_p)) then
              call shdf5_orec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar1_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%rvar2_p)) then
              call shdf5_orec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar2_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%rvar3_p)) then
              call shdf5_orec(ndims, idims, trim(varn), rvara=vtab_r(nv)%rvar3_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)

           elseif (associated(vtab_r(nv)%dvar0_p)) then
              call shdf5_orec(ndims, idims, trim(varn), dvars=vtab_r(nv)%dvar0_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%dvar1_p)) then
              call shdf5_orec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar1_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%dvar2_p)) then
              call shdf5_orec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar2_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           elseif (associated(vtab_r(nv)%dvar3_p)) then
              call shdf5_orec(ndims, idims, trim(varn), dvara=vtab_r(nv)%dvar3_p, &
                              lpoints=ilpts, gpoints=igpts, nglobe=nglobe)
           endif
           
        else

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

        endif

        nvcnt = nvcnt + 1

     endif

  enddo

  call shdf5_close()

end subroutine lite_write

end module lite_vars
