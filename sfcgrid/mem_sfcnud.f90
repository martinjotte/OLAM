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

implicit none

  integer :: nsfcnudfiles, isfcnudfile

  real(r8),          allocatable ::  s1900_sfcnud(:)
  character(len=80), allocatable :: fnames_sfcnud(:)

  real, allocatable ::  sfcwat_nud(:)
  real, allocatable :: sfctemp_nud(:)
  real, allocatable :: fracliq_nud(:)

Contains

!===============================================================================

  subroutine alloc_sfcnud()

  use mem_grid,  only: mza, mva, mwa
  use mem_basic, only: press, rho, theta, vc
  use mem_addsc, only: addsc
  use mem_sfcg,  only: mwsfc

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
  use mem_sfcg,   only: mwsfc
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close

  use mem_plot, only: time8_prev0,         time8_prev1, &
               sfluxr_accum_prev0,  sfluxr_accum_prev1, & ! fast can nud
                  pcp_accum_prev0,     pcp_accum_prev1, & ! fast can nud
              sfctemp_accum_prev0, sfctemp_accum_prev1, & ! fast can nud
              fracliq_accum_prev0, fracliq_accum_prev1    ! fast can nud

  implicit none

  real :: pcp_dif2, sfluxr_dif2
  integer :: ndims, idims(1), iwsfc
  character(pathlen) :: hnamel
  type(simtime) :: ctime

  ! Write surface nudging fields (sfcwat mass, temperature, and fracliq) to a
  ! file.  This subroutine is only executed in a PLOTONLY run, with iparallel = 0

  do iwsfc = 2,mwsfc
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
  idims(1) = mwsfc

  call shdf5_orec(ndims, idims, 'SFCWAT_NUD'  , rvar1=sfcwat_nud)
  call shdf5_orec(ndims, idims, 'SFCTEMP_NUD' , rvar1=sfctemp_nud)
  call shdf5_orec(ndims, idims, 'FRACLIQ_NUD' , rvar1=fracliq_nud)

  call shdf5_close()

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

  integer :: ndims, idims(1)

  ! SFCNUD input files are assumed to be cyclic with a 1-year period

  nsfcnudfiles = 12 ! For the case of one file per month

  allocate (  fnames_sfcnud(nsfcnudfiles))
  allocate (   s1900_sfcnud(nsfcnudfiles))

  ! Convert current model time from s1900 to years, months, dates, hours

  call date_secs_ymdt(s1900_sim,iyears,imonths,idates,ihours)

  ! Fill fnames_sfcnud array with sfcnud file names, including rel or abs path

  fnames_sfcnud( 1) = './hist/climstats1-SN-0000-02-01-000000.h5'
  fnames_sfcnud( 2) = './hist/climstats1-SN-0000-03-01-000000.h5'
  fnames_sfcnud( 3) = './hist/climstats1-SN-0000-04-01-000000.h5'
  fnames_sfcnud( 4) = './hist/climstats1-SN-0000-05-01-000000.h5'
  fnames_sfcnud( 5) = './hist/climstats1-SN-0000-06-01-000000.h5'
  fnames_sfcnud( 6) = './hist/climstats1-SN-0000-07-01-000000.h5'
  fnames_sfcnud( 7) = './hist/climstats1-SN-0000-08-01-000000.h5'
  fnames_sfcnud( 8) = './hist/climstats1-SN-0000-09-01-000000.h5'
  fnames_sfcnud( 9) = './hist/climstats1-SN-0000-10-01-000000.h5'
  fnames_sfcnud(10) = './hist/climstats1-SN-0000-11-01-000000.h5'
  fnames_sfcnud(11) = './hist/climstats1-SN-0000-12-01-000000.h5'
  fnames_sfcnud(12) = './hist/climstats1-SN-0000-01-01-000000.h5'

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

  use misc_coms,  only: io6, current_time
  use hdf5_utils, only: shdf5_open, shdf5_close, shdf5_irec, shdf5_info
  use max_dims,   only: pathlen
  use mem_para,   only: myrank
  use mem_sfcg,   only: itab_wsfc, nwsfc, mwsfc

  implicit none

  integer, intent(in) :: inext

  integer :: isfcnudy, isfcnudm, isfcnudd, isfcnudh
  integer :: nf
  integer :: ndims, idims(3)

  integer, allocatable :: lpoints(:)

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

  ! Open and read sfcnud_database file

  write(io6,*) 'reading sfcnud file ', isfcnudfile, trim(fnames_sfcnud(isfcnudfile))

  call shdf5_open(fnames_sfcnud(isfcnudfile),'R')

  ndims    = 1
  idims(1) = nwsfc

  allocate(lpoints(nwsfc))
  lpoints = itab_wsfc(:)%iwglobe

  call shdf5_irec(ndims, idims, 'SFCWAT_NUD' , rvar1=sfcwat_nud,  points=lpoints)
  call shdf5_irec(ndims, idims, 'SFCTEMP_NUD', rvar1=sfctemp_nud, points=lpoints)
  call shdf5_irec(ndims, idims, 'FRACLIQ_NUD', rvar1=fracliq_nud, points=lpoints)

  deallocate(lpoints)

  call shdf5_close()

  end subroutine sfcnud_read

!===============================================================================

  subroutine read_gw_spinup()

  ! Initialize soil water and energy, and lake energy, from results of spin-up
  ! simulation

  use misc_coms,   only: io6
  use mem_sfcg,    only: itab_wsfc, nwsfc, mwsfc
  use mem_land,    only: land, itab_land, nland, mland, omland, nzg
  use mem_lake,    only: lake, itab_lake, nlake, mlake
  use consts_coms, only: cice1000, cliq1000, alli1000
  use therm_lib,   only: qwtk
  use hdf5_utils,  only: shdf5_open, shdf5_close, shdf5_irec

  implicit none

  ! Set nzg_spinup to nzg value used in spin-up simulation

  integer, parameter :: nzg_sp = 20
  real, parameter :: rhow = 1000. ! density of liquid water [kg/m^3]

  character(80) :: fname
  integer :: ndims, idims(2)

  real, allocatable :: soil_water_sp (:,:)
  real, allocatable :: soil_energy_sp(:,:)

  integer, allocatable :: lpoints(:)

  ! Map 25 current soil layers into 20 spin-up soil layers

  integer :: kspm(25) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,18,19,19,20,20,20,20 /)

  integer :: iland, iwsfc, k, ksp

  real :: tempk, tempc, fracliq

  ! Open and read groundwater spinup file

  fname = './hist/nudgerun1-H-2100-01-01-000000.h5'

  write(io6,*) 'reading gw_spinup file ', trim(fname)

  call shdf5_open(fname,'R')

  ndims    = 1
  idims(1) = nlake

  allocate(lpoints(nlake))
  lpoints = itab_lake(:)%iwglobe

  call shdf5_irec(ndims, idims, 'LAKE%LAKE_ENERGY', rvar1=lake%lake_energy, points=lpoints)

  deallocate(lpoints)

  ndims    = 2
  idims(1) = nzg_sp
  idims(2) = nland

  allocate (soil_water_sp (nzg_sp,mland))
  allocate (soil_energy_sp(nzg_sp,mland))

  allocate(lpoints(nland))
  lpoints = itab_land(:)%iwglobe

  call shdf5_irec(ndims, idims, 'LAND%SOIL_WATER' , rvar2=soil_water_sp,  points=lpoints)
  call shdf5_irec(ndims, idims, 'LAND%SOIL_ENERGY', rvar2=soil_energy_sp, points=lpoints)

  deallocate(lpoints)

  call shdf5_close()

  ! Horizontal loop over land points

  do iland = 2,mland
     iwsfc = iland + omland

     do k = 1,nzg
        ksp = kspm(k)

        if (ksp < 18) then

           ! Copy soil water and energy directly from spin-up layer to current 
           ! layer if layers are identical (hardwired at 18 for specific case)

           land%soil_water (k,iland) = soil_water_sp (ksp,iland)
           land%soil_energy(k,iland) = soil_energy_sp(ksp,iland)

        else

           ! Diagnose spin-up soil temperature and fractional liquid water phase.
           ! Since specifheat_drysoil is not available from the spin-up simulation
           ! (although it could be made available with some effort if necessary),
           ! assume that it is the same as that in level k of the current simulation.

           call qwtk(soil_energy_sp(ksp,iland), soil_water_sp(ksp,iland) * rhow, &
                     land%specifheat_drysoil(k,iland), tempk, fracliq)

           tempc = tempk - 273.15

           ! If we're here, k /= ksp and soil_water_sp(ksp,iland) exceeds the
           ! porosity of level k.  Consequently, set soil_water(k,iland) to
           ! capacity.

           land%soil_water(k,iland) = min(land%wsat_vg(k,iland),soil_water_sp(ksp,iland))

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

  deallocate (soil_water_sp, soil_energy_sp)

  end subroutine read_gw_spinup

End Module mem_sfcnud
