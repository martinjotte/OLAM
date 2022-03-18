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
subroutine isan_driver(iaction)

  use isan_coms,  only: innpr, ihour, idate, imonth, iyear, nfgfiles, ifgfile, &
                        ctotdate_fg, fnames_fg, s1900_fg, lzon_bot, npd, &
                        pcol_p, pcol_z, nprz, o_rho, o_press, o_theta, o_rrw, &
                        o_uzonal, o_umerid, o_ozone, pnpr, levpr, glat
  use misc_coms,  only: io6, runtype, s1900_init, s1900_sim, rinit, i_o3
  use mem_zonavg, only: zonavg_init, zonp_vect
  use mem_grid,   only: mza, mwa
  use mem_nudge,  only: nudflag, nudnxp, o3nudflag

  implicit none

  integer, intent(in) :: iaction

  integer :: nf
  real    :: pmin

  ! Check type of call to isan_driver

  if (iaction == 0) then

! Model run is being started, either for initialization or history restart

! Do inventory of isan file names and times

     call ISAN_file_inv()

! Loop over isan files and search for the one that corresponds to current
! or most recent model time.

     ifgfile = 0
     do nf = 1,nfgfiles
        if (s1900_fg(nf) <= s1900_sim) then
           ifgfile = nf
        endif
     enddo

     if (ifgfile < 1) then
        write(io6,*) ' '
        write(io6,*) 'Unable to find isan file for current model time '
        write(io6,*) 'Stopping model '
        stop 'stop: no current isan file'
     elseif (runtype == 'INITIAL' .and. &
             s1900_fg(ifgfile) < s1900_init - 1800.d0) then
        write(io6,*) ' '
        write(io6,*) 'Available isan file is more than 1800 seconds '
        write(io6,*) 'prior to simulation initialization time. '
        write(io6,*) 'Stopping model '
        stop 'stop: initial isan file'
     endif

  elseif (iaction == 1) then

! Processing next isan file (only called with iaction = 1 if nudflag = 1)

     ifgfile = ifgfile + 1

     if (ifgfile > nfgfiles) then
        write(io6,*) ' '
        write(io6,*) 'No future isan file is available for nudging '
        write(io6,*) 'Stopping model '
        stop 'stop: no future isan file'
     endif

  endif

! Process current isan file

  call date_unmake_big(iyear,imonth,idate,ihour,ctotdate_fg(ifgfile))
  ihour = ihour / 100

  write(io6,*)
  write(io6,*) 'Reading ISAN file ifgfile = ',ifgfile
  write(io6,*) fnames_fg(ifgfile)
  write(io6,'(a,4i6)') ctotdate_fg(ifgfile),iyear,imonth,idate,ihour

  innpr = fnames_fg(ifgfile)

! Fill zonavg arrays for current time

  call zonavg_init(idate,imonth,iyear)

! Read header information from gridded pressure-level files for this file time.
! This information includes input data array dimensions.

  call read_press_header()

! Determine if (and how many) additional levels above the reanalysis we need
! from the ZONAVG arrays

  pmin = 0.5 * zonp_vect(22) + zonp_vect(21)

  if ( pnpr(nprz) < pmin ) then

     lzon_bot = 23
     npd      = nprz + 2

  else

     ! Determine index of lowest ZONAVG pressure level that is at least 1/2
     ! ZONAVG pressure level higher than highest input pressure data level
     ! (i.e., maximum zonp_vect value that is less than 82.5% of levpr(nprz),
     ! which is in hPa)

     lzon_bot = min(23, nint(31. - 6. *  log10( pnpr(nprz) )) + 1)
     npd      = nprz + 25 - lzon_bot

  endif

! Allocate memory for ISAN processing

  allocate( pcol_p  (npd) )     ; pcol_p   = rinit
  allocate( pcol_z  (npd,mwa) ) ; pcol_z   = rinit

  allocate( o_rho   (mza,mwa) ) ; o_rho    = rinit
  allocate( o_press (mza,mwa) ) ; o_press  = rinit
  allocate( o_theta (mza,mwa) ) ; o_theta  = rinit
  allocate( o_rrw   (mza,mwa) ) ; o_rrw    = rinit
  allocate( o_uzonal(mza,mwa) ) ; o_uzonal = rinit
  allocate( o_umerid(mza,mwa) ) ; o_umerid = rinit

  if (i_o3 > 0) then
     allocate( o_ozone (mza,mwa) ) ; o_ozone  = rinit
  endif

! Read gridded pressure-level data and add any ZONAVG fields as necessary

  call pressure_stage()

! Interpolate data to OLAM grid

  call isnstage(iaction)

! If nudging, prepare observational nudging fields

  if (nudflag > 0) then
     if (nudnxp == 0) then
        call nudge_prep_obs (iaction, o_rho, o_theta, o_rrw, o_uzonal, o_umerid)
     else
        call nudge_prep_spec(iaction, o_rho, o_theta, o_rrw, o_uzonal, o_umerid)
     endif
  endif

  if (o3nudflag == 1) then
     call nudge_prep_o3(iaction, o_ozone)
  endif

! Deallocate ISAN arrays

  deallocate( pcol_p  )
  deallocate( pcol_z  )

  deallocate( o_rho   )
  deallocate( o_press )
  deallocate( o_theta )
  deallocate( o_rrw   )
  deallocate( o_uzonal)
  deallocate( o_umerid)

  if (allocated( o_ozone )) deallocate( o_ozone )

  ! These were allocated in read_press_header
  if (allocated( levpr )) deallocate( levpr )
  if (allocated( pnpr  )) deallocate( pnpr  )
  if (allocated( glat  )) deallocate( glat  )

end subroutine isan_driver
