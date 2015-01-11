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
subroutine history_write(vtype)

  use var_tables, only: num_var, vtab_r, get_vtab_dims
  use misc_coms,  only: io6, ioutput, hfilepref, current_time, iclobber, &
                        iparallel
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use max_dims,   only: pathlen
  use mem_grid,   only: nma, nva, nwa
  use leaf_coms,  only: nwl
  use sea_coms,   only: nws
  use mem_nudge,  only: nwnud
  use mem_para,   only: iva_globe_primary, iva_local_primary, mva_primary, &
                        iwa_globe_primary, iwa_local_primary, mwa_primary, &
                        ima_globe_primary, ima_local_primary, mma_primary, &
                        iwl_globe_primary, iwl_local_primary, mwl_primary, &
                        iws_globe_primary, iws_local_primary, mws_primary, &
                        iwnud_globe_primary, iwnud_local_primary, mwnud_primary, &
                        myrank
  implicit none

! This routine writes the chosen variables on the history file.

  character(*), intent(in) :: vtype  ! not used yet - 'state'

  character(pathlen) :: hnamel
  character(32)      :: varn
  character(2)       :: stagpt
  logical            :: exans
  integer            :: nv, nvcnt
  integer            :: ndims, idims(3)
  integer, pointer   :: ilpts(:), igpts(:)
  integer            :: nglobe

  if (ioutput == 0) return

! Construct h5 file name and open the file

  call makefnam(hnamel, hfilepref, current_time, 'H', '$', 'h5')

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
  do nv = 1, num_var

     if (vtab_r(nv)%ihist) then

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
              stop "invalid array size in history_write for parallel output"
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

end subroutine history_write
