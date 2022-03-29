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

Module mem_sfcnud

  use consts_coms, only: r8
  use max_dims,    only: pathlen

  implicit none

  integer :: nsfcnudfiles, isfcnudfile

  real(r8),          allocatable ::  s1900_sfcnud(:)
  character(len=80), allocatable :: fnames_sfcnud(:)

  character(pathlen) :: gw_spinup_sfcgfile, gw_spinup_histfile

  real, allocatable ::  sfcwat_nud(:)
  real, allocatable :: sfctemp_nud(:)
  real, allocatable :: fracliq_nud(:)
  real, allocatable :: wts_sfcnud(:,:)

  integer, allocatable :: iws_sfcnud(:,:)

  integer :: nzg_nl ! namelist value of nzg --- Used in groundwater spin-up run when
  integer :: nzg_sp ! spin-up value of nzg  --- nzg is reduced during initialization

  integer, allocatable :: kspm(:)      ! maps soil layer vertical indices from
                                       ! groundwater spin-up run to regular run

  ! Data structures for file-input sfcnud arrays.  If they have different size
  ! than sfcgrid in current model run, their values will be interpolated.

  Type sfcnudin_vars
     integer :: nvsfc
     integer :: nwsfc

     real, allocatable :: dnv(:)
     real, allocatable :: xew(:)
     real, allocatable :: yew(:)
     real, allocatable :: zew(:)
     real, allocatable :: glatw(:)
     real, allocatable :: glonw(:)
     real, allocatable ::  sfcwat(:)
     real, allocatable :: sfctemp(:)
     real, allocatable :: fracliq(:)
  End type sfcnudin_vars

  Type itab_wsfcnudin_vars
     integer :: npoly = 0
     integer :: ivn(7)
     integer :: iwn(7)
  End type itab_wsfcnudin_vars

  type       (sfcnudin_vars),              target ::       sfcnudin
  type (itab_wsfcnudin_vars), allocatable, target :: itab_wsfcnudin(:)

Contains

!===============================================================================

  subroutine alloc_sfcnud()

  use mem_sfcg, only: mwsfc

  implicit none

  allocate (  sfcwat_nud(mwsfc)) ;  sfcwat_nud(:) = 0.
  allocate ( sfctemp_nud(mwsfc)) ; sfctemp_nud(:) = 0.
  allocate ( fracliq_nud(mwsfc)) ; fracliq_nud(:) = 0.

  end subroutine alloc_sfcnud

!===============================================================================

  subroutine sfcnud_write()

  use misc_coms,  only: io6, current_time, simtime, hfilepref, iclobber
  use max_dims,   only: pathlen
  use mem_para,   only: myrank
  use mem_sfcg,   only: sfcg, nwsfc, nvsfc, itab_wsfc
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close

  use mem_plot, only: time8_prev0,         time8_prev1, &
               sfluxr_accum_prev0,  sfluxr_accum_prev1, & ! fast can nud
                  pcp_accum_prev0,     pcp_accum_prev1, & ! fast can nud
              sfctemp_accum_prev0, sfctemp_accum_prev1, & ! fast can nud
              fracliq_accum_prev0, fracliq_accum_prev1    ! fast can nud

  implicit none

  real :: pcp_dif2, sfluxr_dif2
  integer :: ndims, idims(2), iwsfc
  character(pathlen) :: hnamel
  logical, save :: firstime = .true.
  type(simtime) :: ctime

  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  ! Write surface nudging fields (sfcwat mass, temperature, and fracliq) to a
  ! file.  This subroutine is only executed in a PLOTONLY run, with iparallel = 0

  do iwsfc = 2,nwsfc
        pcp_dif2 =    pcp_accum_prev0(iwsfc) -    pcp_accum_prev1(iwsfc)
     sfluxr_dif2 = sfluxr_accum_prev0(iwsfc) - sfluxr_accum_prev1(iwsfc)

      sfcwat_nud(iwsfc) = pcp_dif2 - sfluxr_dif2
     sfctemp_nud(iwsfc) = sfctemp_accum_prev0(iwsfc) - sfctemp_accum_prev1(iwsfc)
     fracliq_nud(iwsfc) = fracliq_accum_prev0(iwsfc) - fracliq_accum_prev1(iwsfc)

     if (abs(time8_prev0 - time8_prev1) > .99) then
         sfcwat_nud(iwsfc) =  sfcwat_nud(iwsfc) / (time8_prev0 - time8_prev1)
        sfctemp_nud(iwsfc) = sfctemp_nud(iwsfc) / (time8_prev0 - time8_prev1)
        fracliq_nud(iwsfc) = fracliq_nud(iwsfc) / (time8_prev0 - time8_prev1)
     endif
  enddo

  ! Since nudging files will be used cyclicly over a period of 1 year, give
  ! each one a year designation of 0000

  ctime = current_time
  ctime%year = 0

  if (myrank == 0) write(io6,'(/,a)') "Writing surface nudge fields to disk..."

  call makefnam(hnamel, hfilepref, ctime, 'SN', '$', 'h5')
  call shdf5_open(hnamel,'W',iclobber) 

  ndims    = 1
  idims(1) = nwsfc

  call shdf5_orec(ndims, idims, 'SFCWAT_NUD'  , rvar1=sfcwat_nud)
  call shdf5_orec(ndims, idims, 'SFCTEMP_NUD' , rvar1=sfctemp_nud)
  call shdf5_orec(ndims, idims, 'FRACLIQ_NUD' , rvar1=fracliq_nud)

  call shdf5_close()

  ! On first call to this subroutine, write out sfcnud grid file

  if (firstime) then
     firstime = .false.

     call makefnam(hnamel, hfilepref, ctime, 'GN', '$', 'h5')
     call shdf5_open(hnamel,'W',iclobber) 

     ndims = 1
     idims(1) = 1

     call shdf5_orec(ndims, idims, 'nvsfc'  , ivars=nvsfc)
     call shdf5_orec(ndims, idims, 'nwsfc'  , ivars=nwsfc)

     idims(1) = nvsfc

     call shdf5_orec(ndims, idims, 'sfcg%dnv'   , rvar1=sfcg%dnv)

     idims(1) = nwsfc

     call shdf5_orec(ndims, idims, 'sfcg%xew'   , rvar1=sfcg%xew)
     call shdf5_orec(ndims, idims, 'sfcg%yew'   , rvar1=sfcg%yew)
     call shdf5_orec(ndims, idims, 'sfcg%zew'   , rvar1=sfcg%zew)
     call shdf5_orec(ndims, idims, 'sfcg%glatw' , rvar1=sfcg%glatw)
     call shdf5_orec(ndims, idims, 'sfcg%glonw' , rvar1=sfcg%glonw)

     allocate (iscr1(nwsfc))
     do iwsfc = 1,nwsfc
        iscr1(iwsfc) = itab_wsfc(iwsfc)%npoly
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%npoly' ,ivar1=iscr1)
     deallocate(iscr1)

     ndims = 2
     idims(1) = 7
     idims(2) = nwsfc

     allocate (iscr2(7,nwsfc))
     do iwsfc = 1,nwsfc
       iscr2(1:7,iwsfc) = itab_wsfc(iwsfc)%ivn(1:7)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%ivn',ivar2=iscr2)

     do iwsfc = 1,nwsfc
        iscr2(1:7,iwsfc) = itab_wsfc(iwsfc)%iwn(1:7)
     enddo
     call shdf5_orec(ndims,idims,'itab_wsfc%iwn',ivar2=iscr2)
     deallocate(iscr2)

     call shdf5_close()

  endif  ! firstime

  end subroutine sfcnud_write

!===============================================================================

  subroutine sfcnud_read_init()

  use misc_coms,  only: io6, current_time, s1900_sim
  use max_dims,   only: pathlen
  use mem_para,   only: myrank

  implicit none

  integer :: isfcnudy,isfcnudm,isfcnudd,isfcnudh
  integer :: iyears, imonths, idates, ihours
  integer :: nf

  ! SFCNUD input files are assumed to be cyclic with a 1-year period

  nsfcnudfiles = 12 ! For the case of one file per month

  allocate (  fnames_sfcnud(nsfcnudfiles))
  allocate (   s1900_sfcnud(nsfcnudfiles))

  ! Convert current model time from s1900 to years, months, dates, hours

  call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

  ! Fill fnames_sfcnud array with sfcnud file names, including rel or abs path

  fnames_sfcnud( 1) = './hist/clim4-SN-0000-02-01-000000.h5'
  fnames_sfcnud( 2) = './hist/clim4-SN-0000-03-01-000000.h5'
  fnames_sfcnud( 3) = './hist/clim4-SN-0000-04-01-000000.h5'
  fnames_sfcnud( 4) = './hist/clim4-SN-0000-05-01-000000.h5'
  fnames_sfcnud( 5) = './hist/clim4-SN-0000-06-01-000000.h5'
  fnames_sfcnud( 6) = './hist/clim4-SN-0000-07-01-000000.h5'
  fnames_sfcnud( 7) = './hist/clim4-SN-0000-08-01-000000.h5'
  fnames_sfcnud( 8) = './hist/clim4-SN-0000-09-01-000000.h5'
  fnames_sfcnud( 9) = './hist/clim4-SN-0000-10-01-000000.h5'
  fnames_sfcnud(10) = './hist/clim4-SN-0000-11-01-000000.h5'
  fnames_sfcnud(11) = './hist/clim4-SN-0000-12-01-000000.h5'
  fnames_sfcnud(12) = './hist/clim4-SN-0000-01-01-000000.h5'

  ! Fill s1900_sfcnud array with time for each sfcnud file

  call date_abs_secs2(iyears,  02,01,000000,s1900_sfcnud( 1))
  call date_abs_secs2(iyears,  03,01,000000,s1900_sfcnud( 2))
  call date_abs_secs2(iyears,  04,01,000000,s1900_sfcnud( 3))
  call date_abs_secs2(iyears,  05,01,000000,s1900_sfcnud( 4))
  call date_abs_secs2(iyears,  06,01,000000,s1900_sfcnud( 5))
  call date_abs_secs2(iyears,  07,01,000000,s1900_sfcnud( 6))
  call date_abs_secs2(iyears,  08,01,000000,s1900_sfcnud( 7))
  call date_abs_secs2(iyears,  09,01,000000,s1900_sfcnud( 8))
  call date_abs_secs2(iyears,  10,01,000000,s1900_sfcnud( 9))
  call date_abs_secs2(iyears,  11,01,000000,s1900_sfcnud(10))
  call date_abs_secs2(iyears,  12,01,000000,s1900_sfcnud(11))
  call date_abs_secs2(iyears+1,01,01,000000,s1900_sfcnud(12))

  ! Loop over number of sfcnud_DATABASE file times and search for the first one
  ! that corresponds to future time

  isfcnudfile = 0
  do nf = 1, nsfcnudfiles
     write(io6,*) 'nsfcnudf0 ',nf,s1900_sfcnud(nf),' ',s1900_sim

     if (s1900_sfcnud(nf) > s1900_sim) then
        isfcnudfile = nf
        exit
     endif
  enddo

  if (isfcnudfile < 1) then
     write(io6,*) ' '
     write(io6,*) 'Unable to find future sfcnud file for current'
     write(io6,*) 'model time.  Stopping model.'
     stop 'stop: no future sfcnud file'
  endif

  ! Read first file

  call sfcnud_read(0)

  end subroutine sfcnud_read_init

!===============================================================================

  subroutine sfcnud_read(inext)

  use misc_coms,  only: io6, current_time, simtime, hfilepref
  use max_dims,   only: pathlen
  use mem_para,   only: myrank
  use mem_sfcg,   only: itab_wsfc, nwsfc, mwsfc
  use hdf5_utils, only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info

  implicit none

  integer, intent(in) :: inext

  integer :: isfcnudy, isfcnudm, isfcnudd, isfcnudh
  integer :: nf, iw, iwsfc, iwglobe
  integer :: ndims, idims(2)
  character(pathlen) :: hnamel
  type(simtime) :: ctime

  logical, save :: firstime = .true.

  integer, allocatable :: iscr1(:)
  integer, allocatable :: iscr2(:,:)

  ! Processing next sfcnud file

  if (inext > 0) isfcnudfile = isfcnudfile + 1

  if (isfcnudfile > nsfcnudfiles) then

     ! Advance sfcnud files to next year

     isfcnudfile = 1

     do nf = 1, nsfcnudfiles
        call date_secs_ymdt(s1900_sfcnud(nf),isfcnudy,isfcnudm,isfcnudd,isfcnudh)
        call date_abs_secs2(isfcnudy+1,isfcnudm,isfcnudd,isfcnudh,s1900_sfcnud(nf))
     enddo

  endif

  ! On first call to this subroutine, read sfcnud grid file and allocate input arrays

  if (firstime) then
     firstime = .false.

     hnamel = './hist/clim4-GN-0000-01-01-000000.h5'

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening sfcnud file '//trim(hnamel)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(hnamel,'R')

     ndims = 1
     idims(1) = 1

     call shdf5_irec(ndims, idims, 'nvsfc'  , ivars=sfcnudin%nvsfc)
     call shdf5_irec(ndims, idims, 'nwsfc'  , ivars=sfcnudin%nwsfc)

     allocate (sfcnudin%sfcwat (sfcnudin%nwsfc))
     allocate (sfcnudin%sfctemp(sfcnudin%nwsfc))
     allocate (sfcnudin%fracliq(sfcnudin%nwsfc))

     ! If sfcnud data read from file did not use the same grid as the current
     ! model run, allocate additional grid coordinate arrays, read their values
     ! from the grid file, and compute interpolation coefficients in preparation
     ! for interpolating input values to grid cells in current model run.

     if (nwsfc /= sfcnudin%nwsfc) then ! Assume this check to be sufficient

        allocate (sfcnudin%dnv  (sfcnudin%nvsfc))
        allocate (sfcnudin%xew  (sfcnudin%nwsfc))
        allocate (sfcnudin%yew  (sfcnudin%nwsfc))
        allocate (sfcnudin%zew  (sfcnudin%nwsfc))
        allocate (sfcnudin%glatw(sfcnudin%nwsfc))
        allocate (sfcnudin%glonw(sfcnudin%nwsfc))

        allocate (itab_wsfcnudin(sfcnudin%nwsfc))

        idims(1) = sfcnudin%nvsfc

        call shdf5_irec(ndims, idims, 'sfcg%dnv'    , rvar1=sfcnudin%dnv)

        idims(1) = sfcnudin%nwsfc

        call shdf5_irec(ndims, idims, 'sfcg%xew'    , rvar1=sfcnudin%xew)
        call shdf5_irec(ndims, idims, 'sfcg%yew'    , rvar1=sfcnudin%yew)
        call shdf5_irec(ndims, idims, 'sfcg%zew'    , rvar1=sfcnudin%zew)
        call shdf5_irec(ndims, idims, 'sfcg%glatw'  , rvar1=sfcnudin%glatw)
        call shdf5_irec(ndims, idims, 'sfcg%glonw'  , rvar1=sfcnudin%glonw)

        allocate (iscr1(sfcnudin%nwsfc))
        call shdf5_irec(ndims,idims,'itab_wsfc%npoly' ,ivar1=iscr1)
        do iw = 1,sfcnudin%nwsfc
           itab_wsfcnudin(iw)%npoly = iscr1(iw)
        enddo
        deallocate(iscr1)

        ndims = 2
        idims(1) = 7
        idims(2) = sfcnudin%nwsfc

        allocate (iscr2(7,sfcnudin%nwsfc))
        call shdf5_irec(ndims,idims,'itab_wsfc%ivn',ivar2=iscr2)
        do iw = 1,sfcnudin%nwsfc
           itab_wsfcnudin(iw)%ivn(1:7) = iscr2(1:7,iw)
        enddo

        call shdf5_irec(ndims,idims,'itab_wsfc%iwn',ivar2=iscr2)
        do iw = 1,sfcnudin%nwsfc
           itab_wsfcnudin(iw)%iwn(1:7) = iscr2(1:7,iw)
        enddo
        deallocate(iscr2)

        allocate (iws_sfcnud(nwsfc,3)); iws_sfcnud = 0
        allocate (wts_sfcnud(nwsfc,3)); wts_sfcnud = 0.
 
        call find_3iws_sfcnud()

     endif

     call shdf5_close()

  endif  ! firstime

  ! Open and read sfcnud_database file

  write(io6,*) 'reading sfcnud file ', isfcnudfile, trim(fnames_sfcnud(isfcnudfile))

  call shdf5_open(fnames_sfcnud(isfcnudfile),'R')

  ndims    = 1
  idims(1) = sfcnudin%nwsfc

  call shdf5_irec(ndims, idims, 'SFCWAT_NUD' , rvar1=sfcnudin%sfcwat)
  call shdf5_irec(ndims, idims, 'SFCTEMP_NUD', rvar1=sfcnudin%sfctemp)
  call shdf5_irec(ndims, idims, 'FRACLIQ_NUD', rvar1=sfcnudin%fracliq)

  call shdf5_close()

  ! If sfcnud data read from file used the same grid as the current model run,
  ! copy file values to sfcnud arrays.  Otherwise, interpolate them.

  if (nwsfc == sfcnudin%nwsfc) then ! Assume this check to be sufficient

     do iwsfc = 1, mwsfc
        iwglobe = itab_wsfc(iwsfc)%iwglobe
        sfcwat_nud (iwsfc) = sfcnudin%sfcwat (iwglobe)
        sfctemp_nud(iwsfc) = sfcnudin%sfctemp(iwglobe)
        fracliq_nud(iwsfc) = sfcnudin%fracliq(iwglobe)
     enddo

  else

     do iwsfc = 1, mwsfc
        sfcwat_nud(iwsfc)  = wts_sfcnud(iwsfc,1) * sfcnudin%sfcwat(iws_sfcnud(iwsfc,1)) &
                           + wts_sfcnud(iwsfc,2) * sfcnudin%sfcwat(iws_sfcnud(iwsfc,2)) &
                           + wts_sfcnud(iwsfc,3) * sfcnudin%sfcwat(iws_sfcnud(iwsfc,3))

        sfctemp_nud(iwsfc) = wts_sfcnud(iwsfc,1) * sfcnudin%sfctemp(iws_sfcnud(iwsfc,1)) &
                           + wts_sfcnud(iwsfc,2) * sfcnudin%sfctemp(iws_sfcnud(iwsfc,2)) &
                           + wts_sfcnud(iwsfc,3) * sfcnudin%sfctemp(iws_sfcnud(iwsfc,3))

        fracliq_nud(iwsfc) = wts_sfcnud(iwsfc,1) * sfcnudin%fracliq(iws_sfcnud(iwsfc,1)) &
                           + wts_sfcnud(iwsfc,2) * sfcnudin%fracliq(iws_sfcnud(iwsfc,2)) &
                           + wts_sfcnud(iwsfc,3) * sfcnudin%fracliq(iws_sfcnud(iwsfc,3))

     enddo

  endif

  end subroutine sfcnud_read

!================================================================================

  subroutine find_3iws_sfcnud()

  use mem_sfcg,    only: sfcg, mwsfc
  use consts_coms, only: pio180

  implicit none

  integer :: iw, npoly, j, iwn, j1, j2, iwsfc

  real :: coswlon, sinwlon
  real :: coswlat, sinwlat
  real :: dxe, dye, dze

  real :: xwn(7), ywn(7), dot00(7), dot01(7), dot11(7), denomi(7), dn(7)
  real :: qx, qy, dnvmax, dist, distn
  real :: dot02,dot12,u,v

  real, parameter :: fuzz = 0.001

  ! Loop over all SFCNUD points in arrays read from file

  do iw = 2, sfcnudin%nwsfc

     sinwlat = sin(sfcnudin%glatw(iw) * pio180)
     coswlat = cos(sfcnudin%glatw(iw) * pio180)
     sinwlon = sin(sfcnudin%glonw(iw) * pio180)
     coswlon = cos(sfcnudin%glonw(iw) * pio180)

     ! Find max distance to neighbor W points

     npoly = itab_wsfcnudin(iw)%npoly
     dnvmax = maxval(sfcnudin%dnv(itab_wsfcnudin(iw)%ivn(1:npoly)))

     ! Loop over neighbor W points

     do j = 1,npoly
        iwn = itab_wsfcnudin(iw)%iwn(j)

        ! Transform neighbor W points to PS coordinates tangent at IW point

        dxe = sfcnudin%xew(iwn) - sfcnudin%xew(iw)
        dye = sfcnudin%yew(iwn) - sfcnudin%yew(iw)
        dze = sfcnudin%zew(iwn) - sfcnudin%zew(iw)

        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,xwn(j),ywn(j))
     enddo

     ! Loop over each pair of consecutive neighbor W points, and set up
     ! triangle-check coefficients that depend only on W points

     do j1 = 1,npoly
        j2 = j1 + 1
        if (j1 == npoly) j2 = 1

        dot00(j1) = xwn(j1) * xwn(j1) + ywn(j1) * ywn(j1)
        dot01(j1) = xwn(j1) * xwn(j2) + ywn(j1) * ywn(j2)
        dot11(j1) = xwn(j2) * xwn(j2) + ywn(j2) * ywn(j2)

        denomi(j1) = 1. / (dot00(j1) * dot11(j1) - dot01(j1) * dot01(j1))
        dn(j1)     = 1. / (xwn(j2) * ywn(j1) - xwn(j1) * ywn(j2))
     enddo
           
     ! Loop over all iwsfc points and determine which are closer to current
     ! input IW point than to any other input IW point on globe.  It is
     ! sufficient to show that a iwsfc point is closer to current input IW
     ! point than to any neighbor input IW point, AND that it is closer to
     ! current input IW point than most distant neighbor input IW point is.

     do iwsfc = 2, mwsfc
        if (abs(sfcg%zew(iwsfc) - sfcnudin%zew(iw)) > dnvmax) cycle

        ! Skip this iwsfc point if its nearest input IW point has already been found

        if (iws_sfcnud(iwsfc,1) > 1) cycle

        ! Skip this iwsfc point if any of its earth coordinates differs from
        ! input IW point by more than dnvmax

        if (abs(sfcg%xew(iwsfc) - sfcnudin%xew(iw)) > dnvmax) cycle
        if (abs(sfcg%yew(iwsfc) - sfcnudin%yew(iw)) > dnvmax) cycle

        ! Compute distance between IWSFC point and IW point

        dist = sqrt((sfcg%xew(iwsfc)-sfcnudin%xew(iw))**2 &
                  + (sfcg%yew(iwsfc)-sfcnudin%yew(iw))**2 &
                  + (sfcg%zew(iwsfc)-sfcnudin%zew(iw))**2)

        ! Skip this iwsfc point if it is farther from input IW point than most
        ! distant neighbor IWN point is

        if (dist > dnvmax) cycle


        ! Loop over neighbor W points

        do j = 1,npoly
           iwn = itab_wsfcnudin(iw)%iwn(j)

           ! Compute distance between lat/lon point and IWN point

           distn = sqrt((sfcg%xew(iwsfc)-sfcnudin%xew(iwn))**2 &
                      + (sfcg%yew(iwsfc)-sfcnudin%yew(iwn))**2 &
                      + (sfcg%zew(iwsfc)-sfcnudin%zew(iwn))**2)

           ! If iwsfc point is closer to IWN point than to input IW point, move
           ! on to next iwsfc point.  Bias is used to reduce chance of
           ! iwsfc point being rejected by all input IW points in domain; this 
           ! might lead to a few iwsfc values being interpolated on
           ! multiple MPI subdomains, but this is sorted out later.

           if (distn < 0.999999 * dist) goto 10
        enddo

        ! If this point was reached, current iwsfc point is inside input IW cell.
        ! Store input IW index for the iwsfc point and transform iwsfc point
        ! to PS coordinates tangent at input IW point.

        dxe = sfcg%xew(iwsfc) - sfcnudin%xew(iw)
        dye = sfcg%yew(iwsfc) - sfcnudin%yew(iw)
        dze = sfcg%zew(iwsfc) - sfcnudin%zew(iw)

        call de_ps(dxe,dye,dze,coswlat,sinwlat,coswlon,sinwlon,qx,qy)

        ! Loop over each pair of consecutive neighbor input IW points 

        do j1 = 1,npoly
           j2 = j1 + 1
           if (j1 == npoly) j2 = 1

           ! Set up triangle-check coefficients that depend on lat/lon points

           dot02 = xwn(j1) * qx + ywn(j1) * qy
           dot12 = xwn(j2) * qx + ywn(j2) * qy

           u = (dot11(j1) * dot02 - dot01(j1) * dot12) * denomi(j1)
           v = (dot00(j1) * dot12 - dot01(j1) * dot02) * denomi(j1)

           if (u > -fuzz .and. v > -fuzz .and. u + v < 1.0 + fuzz) then

              ! lat/lon point is inside current triangle; store indices of
              ! both IWN neighbors

              iws_sfcnud(iwsfc,1) = iw
              iws_sfcnud(iwsfc,2) = itab_wsfcnudin(iw)%iwn(j1)
              iws_sfcnud(iwsfc,3) = itab_wsfcnudin(iw)%iwn(j2)

              ! Compute and store 3 interpolation weights

              wts_sfcnud(iwsfc,2) = dn(j1) * (-ywn(j2) * qx + xwn(j2) * qy)
              wts_sfcnud(iwsfc,3) = dn(j1) * ( ywn(j1) * qx - xwn(j1) * qy)
              wts_sfcnud(iwsfc,1) = 1. - wts_sfcnud(iwsfc,2) - wts_sfcnud(iwsfc,3)

              exit

           endif

        enddo  ! j1 loop

        10 continue

     enddo ! iwsfc loop

  enddo  ! iw loop

  end subroutine find_3iws_sfcnud

!===============================================================================

  subroutine read_gw_spinup()

  ! Initialize soil water and energy, and lake energy, from results of spin-up
  ! simulation

  use misc_coms,   only: io6
  use mem_sfcg,    only: itab_wsfc, nwsfc, mwsfc
  use mem_land,    only: land, itab_land, nland, mland, omland, nzg, slzt
  use mem_lake,    only: lake, itab_lake, nlake, mlake
  use leaf4_soil,  only: soil_pot2wat
  use consts_coms, only: cice1000, cliq1000, alli1000
  use therm_lib,   only: qwtk
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec

  implicit none

  real,    parameter :: rhow = 1000. ! density of liquid water [kg/m^3]

  character(80) :: fname
  integer :: ndims, idims(2)

  real, allocatable :: soil_water_sp (:,:)
  real, allocatable :: soil_energy_sp(:,:)
  real, allocatable :: head_sp       (:,:)

  character(2) :: type

  integer :: iland, iwsfc, k, ksp

  real :: tempk, tempc, fracliq, psi

  ! Pointers to the global index of the local point

  integer :: lglake(mlake)
  integer :: lgland(mland)

  lglake = itab_lake(1:mlake)%iwglobe
  lgland = itab_land(1:mland)%iwglobe

  ! Open and read sfcgrid file from groundwater spinup simulation

  fname = trim(gw_spinup_sfcgfile)

  write(io6,*) 'reading gw_spinup sfcgfile ', fname

  call shdf5_open(fname,'R')

  ndims = 1
  idims(1) = 1

  call shdf5_irec(ndims, idims, 'nzg_nl' , ivars=nzg_nl)
  call shdf5_irec(ndims, idims, 'nzg_sp' , ivars=nzg_sp)

  if (nzg_nl /= nzg) then
     print*, 'nzg    = ',nzg
     print*, 'nzg_nl = ',nzg_nl
     print*, 'nzg_sp = ',nzg_sp
     print*, 'nzg of this gw-initialized simulation does not equal '
     print*, 'nzg_nl of the gw-spinup simulation '
     stop 'stop nzg_nl '
  endif

  allocate (kspm(nzg_nl))

  idims(1) = nzg_nl

  call shdf5_irec(ndims, idims, 'kspm'   , ivar1=kspm)

  call shdf5_close()

  ! Open and read history file from groundwater spinup simulation

  fname = trim(gw_spinup_histfile)

  write(io6,*) 'reading gw_spinup histfile ', trim(fname)

  call shdf5_open(fname,'R')

  ndims    = 1
  idims(1) = mlake
  type     = 'RW'

  call shdf5_irec(ndims, idims, 'LAKE%DEPTH',       rvar1=lake%depth,       points=lglake, stagpt=type)
  call shdf5_irec(ndims, idims, 'LAKE%LAKE_ENERGY', rvar1=lake%lake_energy, points=lglake, stagpt=type)

  ndims    = 2
  idims(1) = nzg_sp
  idims(2) = mland
  type     = 'LW'

  allocate (soil_water_sp (nzg_sp,mland))
  allocate (soil_energy_sp(nzg_sp,mland))
  allocate (head_sp       (nzg_sp,mland))

  call shdf5_irec(ndims, idims, 'LAND%SOIL_WATER' , rvar2=soil_water_sp,  points=lgland, stagpt=type)
  call shdf5_irec(ndims, idims, 'LAND%SOIL_ENERGY', rvar2=soil_energy_sp, points=lgland, stagpt=type)
  call shdf5_irec(ndims, idims, 'LAND%HEAD',        rvar2=head_sp,        points=lgland, stagpt=type)

  call shdf5_close()

  ! Horizontal loop over land points

  do iland = 2,mland
     iwsfc = iland + omland

     do k = 1,nzg
        ksp = kspm(k)

        if (ksp == k) then

           ! Copy soil water and energy directly from spin-up layer to current 
           ! layer if layers are identical

           land%soil_water (k,iland) = soil_water_sp (ksp,iland)
           land%soil_energy(k,iland) = soil_energy_sp(ksp,iland)

        else

           ! Here in the upper soil layers, k /= ksp, so spun-up soil layers
           ! are thicker than layers in the present simulation, and they may
           ! have different soil properties.

           ! Assign initial soil moisture in the present simulation such that
           ! its head value equals that of the spun-up simulation.

           psi = head_sp(ksp,iland) - slzt(k)

           call soil_pot2wat(psi, land%wresid_vg(k,iland), land%wsat_vg(k,iland), &
                             land%alpha_vg(k,iland), land%en_vg(k,iland), &
                             land%soil_water(k,iland))

           land%soil_water(k,iland) = max(land%wresid_vg(k,iland), &
                                      min(land%wsat_vg(k,iland),   &
                                      land%soil_water(k,iland)))

           ! Diagnose spun-up soil temperature and fractional liquid water phase.
           ! Since specifheat_drysoil is not available from the spin-up simulation
           ! (although it could be made available with some effort if necessary),
           ! assume that it is the same as that in level k of the current simulation.

           call qwtk(soil_energy_sp(ksp,iland), soil_water_sp(ksp,iland) * rhow, &
                     land%specifheat_drysoil(k,iland), tempk, fracliq)

           tempc = tempk - 273.15

           ! Diagnose corresponding soil energy

           if (tempc > 0.) then
              land%soil_energy(k,iland) =   tempc * land%specifheat_drysoil(k,iland)    &
                                        +   tempc * land%soil_water(k,iland) * cliq1000 &
                                        + fracliq * land%soil_water(k,iland) * alli1000
           else
              land%soil_energy(k,iland) =   tempc * land%specifheat_drysoil(k,iland)    &
                                        +   tempc * land%soil_water(k,iland) * cice1000 &
                                        + fracliq * land%soil_water(k,iland) * alli1000
           endif

        endif

     enddo

  enddo

  deallocate (soil_water_sp, soil_energy_sp, head_sp)

  end subroutine read_gw_spinup

End Module mem_sfcnud
