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

subroutine read_nl(file)

use oname_coms, only: nl, cmdlne_runtype, cmdlne_fields, numcf
use misc_coms,  only: io6
use max_dims,   only: pathlen

implicit none

character(*), intent(in) :: file
character(pathlen)       :: fs
integer                  :: il, ios
logical                  :: fexists

namelist /OLAMIN/ nl

! OPEN THE NAMELIST FILE

inquire(file=file, exist=fexists)
if (.not. fexists) then
   write(io6,*) "The namelist file "//trim(file)//" is missing."
   stop "Stopping model run."
endif
open(10, status='OLD', file=file)
 
! READ GRID POINT, MODEL OPTIONS, AND PLOTTING INFORMATION FROM THE NAMELIST

read(10, nml=OLAMIN)
close(10)

! OVERWRITE NAMELIST WITH COMMAND LINE SWITCHES

fexists = .false.

do il = 1, numcf
   fs = "&OLAMIN NL%" // adjustl(trim(cmdlne_fields(il))) // " /"

   write(io6,*) trim(fs)

   read(fs, nml=OLAMIN, iostat=ios)

   if (ios /= 0) then
      fexists = .true.
      write(io6,*)
      write(io6,*) "Error setting name list values from command line:"
      write(io6,*) trim(cmdlne_fields(il))
      write(io6,*) ios
      write(io6,*) trim(fs)
   endif
enddo
if (fexists) stop "Stopping model run."

! OVERWRITE NAMELIST WITH COMMAND LINE RUNTYPE

if (len_trim(cmdlne_runtype) > 1) then
   nl%runtype = cmdlne_runtype
endif

end subroutine read_nl

!===============================================================================

subroutine copy_nl()

  use max_dims,    only: maxgrds, nzgmax, maxisdirs
  use consts_coms, only: r8
  use oname_coms,  only: nl
  use misc_coms,   only: io6, expnme, runtype, timeunit, timmax8, ndtrat, &
                         nacoust, idiffk, csz, csx, akmin, &
                         dtlong, initial, zonclim, topo_database, &
                         gridfile, hfilin, ioutput, hfilepref, iclobber, &
                         frqstate, naddsc, ilwrtyp, iswrtyp, radfrq, &
                         icfrac, cfracrh1, cfracrh2, cfraccup, nqparm, confrq, &
                         nsndg, ipsflg, itsflg, irtsflg, iusflg, &
                         hs, p_sfc, us, vs, ts, ps, rts, &
                         itime1, idate1, imonth1, iyear1, ngrids, ngrids_old, &
                         nzp, mdomain, itopoflg, nxp, &
                         ngrdll, grdrad, grdlat, grdlon, deltax, ndz, hdz, dz, &
                         current_time, debug_fp, init_nans, do_chem

  use micro_coms,  only: miclevel, icloud, idriz, irain, ipris, isnow, iaggr, &
                         igraup, ihail, iccn, igccn, iifn, &
                         rparm, sparm, aparm, gparm, hparm, &
                         ccnparm, gccnparm, ifnparm

  use mem_co2,     only: co2flag, co2_initppm

  use leaf_coms,   only: nvgcon, nslcon, isoilflg, ndviflg, &
                         isfcl, ivegflg, nzg, nzs, slz, &
                         veg_database, soil_database, &
                         ndvi_database, iupdndvi, landusefile, &
                         isoilstateinit, iwatertabflg, watertab_db

  use sea_coms,    only: isstflg, sst_database, seatmp, seafile, iupdsst, &
                         iseaiceflg, seaice_database, seaice, iupdseaice, &
                         iseagrid

  use oplot_coms,  only: op
  use isan_coms,   only: iapr
  use mem_nudge,   only: tnudcent, nudflag, nudnxp,  &
                         o3nudflag, tnudi_o3, o3nudpress
  use mem_rayf,    only: rayf_zmin,    rayf_distim,   rayf_expon,    &
                         rayfw_zmin,   rayfw_distim,   rayfw_expon,  &
                         rayfdiv_zmin, rayfdiv_distim, rayfdiv_expon
  use ed_misc_coms, only: ed2_active, ed2_namelist

  implicit none

  integer  :: i,j
  real(r8) :: tfact

! The variables in this section are always copied from the namelist for any
! RUNTYPE.  This allows some model options to be changed on a history start.

  expnme   = nl%expnme
  runtype  = nl%runtype
  timeunit = nl%timeunit
  timmax8  = nl%timmax

!----------------------------------------------------------
! Convert timmax8 units to seconds if necessary

  if (timeunit == 'd' .or. timeunit == 'D') tfact = 86400.0_r8
  if (timeunit == 'h' .or. timeunit == 'H') tfact =  3600.0_r8
  if (timeunit == 'm' .or. timeunit == 'M') tfact =    60.0_r8
  if (timeunit == 's' .or. timeunit == 'S') tfact =     1.0_r8

  timmax8 = timmax8 * tfact
!----------------------------------------------------------

  do i = 1,maxgrds
     ndtrat(i)  = nl%ndtrat(i)
     nacoust(i) = nl%nacoust(i)
     idiffk(i)  = nl%idiffk(i)
     csz(i)     = nl%csz(i)
     csx(i)     = nl%csx(i)
     akmin(i)   = nl%akmin(i)
     nqparm(i)  = nl%nqparm(i)
  enddo

  topo_database = nl%topo_database
  sst_database  = nl%sst_database
  seaice_database  = nl%seaice_database
  veg_database  = nl%veg_database
  soil_database = nl%soil_database
  ndvi_database = nl%ndvi_database
  watertab_db   = nl%watertab_db

  ed2_active = nl%ed2_active
  ed2_namelist = nl%ed2_namelist

  dtlong        = nl%dtlong
  initial       = nl%initial
  zonclim       = nl%zonclim
  tnudcent      = nl%tnudcent
  nudflag       = nl%nudflag
  hfilin        = nl%hfilin
  ioutput       = nl%ioutput
  hfilepref     = nl%hfilepref
  iclobber      = nl%iclobber
  frqstate      = nl%frqstate
  gridfile      = nl%gridfile
  landusefile   = nl%landusefile
  seafile       = nl%seafile
  iupdsst       = nl%iupdsst
  iupdseaice    = nl%iupdseaice
  iupdndvi      = nl%iupdndvi
  naddsc        = nl%naddsc
  debug_fp      = nl%debug_fp
  init_nans     = nl%init_nans
  rayf_zmin     = nl%rayf_zmin
  rayf_distim   = nl%rayf_distim
  rayf_expon    = nl%rayf_expon
  rayfw_zmin    = nl%rayfw_zmin
  rayfw_distim  = nl%rayfw_distim
  rayfw_expon   = nl%rayfw_expon
  rayfdiv_zmin  = nl%rayfdiv_zmin
  rayfdiv_distim= nl%rayfdiv_distim
  rayfdiv_expon = nl%rayfdiv_expon
  ilwrtyp       = nl%ilwrtyp
  iswrtyp       = nl%iswrtyp
  radfrq        = nl%radfrq
  icfrac        = nl%icfrac
  cfracrh1      = nl%cfracrh1
  cfracrh2      = nl%cfracrh2
  cfraccup      = nl%cfraccup
  confrq        = nl%confrq
  seatmp        = nl%seatmp
  seaice        = nl%seaice
  miclevel      = nl%miclevel
  icloud        = nl%icloud
  idriz         = nl%idriz
  irain         = nl%irain
  ipris         = nl%ipris
  isnow         = nl%isnow
  iaggr         = nl%iaggr
  igraup        = nl%igraup
  ihail         = nl%ihail
  iccn          = nl%iccn
  igccn         = nl%igccn
  iifn          = nl%iifn
  rparm         = nl%rparm
  sparm         = nl%sparm
  aparm         = nl%aparm
  gparm         = nl%gparm
  hparm         = nl%hparm
  ccnparm       = nl%ccnparm
  gccnparm      = nl%gccnparm
  ifnparm       = nl%ifnparm
  co2flag       = nl%co2flag
  co2_initppm   = nl%co2_ppmv_init
  nvgcon        = nl%nvgcon
  nslcon        = nl%nslcon
  nsndg         = nl%nsndg
  ipsflg        = nl%ipsflg
  itsflg        = nl%itsflg
  irtsflg       = nl%irtsflg
  iusflg        = nl%iusflg
  nudnxp        = nl%nudnxp
  isstflg       = nl%isstflg
  iseaiceflg    = nl%iseaiceflg
  isoilflg      = nl%isoilflg
  isoilstateinit= nl%isoilstateinit
  iwatertabflg  = nl%iwatertabflg
  ndviflg       = nl%ndviflg

  do_chem       = nl%do_chem
  o3nudflag     = nl%o3nudflag
  tnudi_o3      = 1.0 / max( nl%o3tnudcent, real(nl%dtlong) )
  o3nudpress    = nl%o3nudpress * 100.0  ! mb to Pa

  iapr(1:maxisdirs) = nl%iapr(1:maxisdirs)

  hs(1) = nl%hs
  p_sfc = nl%p_sfc
  do i = 1,nl%nsndg
     ps(i)  = nl%sounding(1,i)
     ts(i)  = nl%sounding(2,i)
     rts(i) = nl%sounding(3,i)
     us(i)  = nl%sounding(4,i)
     vs(i)  = nl%sounding(5,i)
  enddo

! Variables from $MODEL_PLOT namelist

  op%nplt       = nl%nplt
  op%nplt_files = nl%nplt_files
  op%frqplt     = nl%frqplt
  op%dtvec      = nl%dtvec
  op%headspeed  = nl%headspeed
  op%stemlength = nl%stemlength
  op%pltname    = nl%pltname
  op%plttype    = nl%plttype
  op%pltorient  = nl%pltorient
  op%vec_maxmrl = nl%vec_maxmrl
  op%prtval_size= nl%prtval_size

  do i = 1,op%nplt
     op%fldname(i)   = nl%plotspecs(i)%fldname
     op%projectn(i)  = nl%plotspecs(i)%projectn
     op%icolortab(i) = nl%plotspecs(i)%icolortab
     op%slabloc(i)   = nl%plotspecs(i)%slabloc

     if (index(nl%plotspecs(i)%pltspec2,'T') > 0) op%contrtyp(i) = 'T'
     if (index(nl%plotspecs(i)%pltspec2,'F') > 0) op%contrtyp(i) = 'F'
     if (index(nl%plotspecs(i)%pltspec2,'L') > 0) op%contrtyp(i) = 'L'
     if (index(nl%plotspecs(i)%pltspec2,'O') > 0) op%contrtyp(i) = 'O'

     if (index(nl%plotspecs(i)%pltspec2,'P') > 0) op%prtval(i) = 'P'

     if (index(nl%plotspecs(i)%pltspec2,'I') > 0) op%pltindx(i) = 'I'
     if (index(nl%plotspecs(i)%pltspec2,'J') > 0) op%pltindx(i) = 'J'

     if (index(nl%plotspecs(i)%pltspec2,'B') > 0) op%vectbarb(i) = 'B'
     if (index(nl%plotspecs(i)%pltspec2,'V') > 0) op%vectbarb(i) = 'V'
     if (index(nl%plotspecs(i)%pltspec2,'w') > 0) op%vectbarb(i) = 'w'

     if (index(nl%plotspecs(i)%pltspec2,'G') > 0) op%pltgrid(i) = 'G'
     if (index(nl%plotspecs(i)%pltspec2,'g') > 0) op%pltgrid_landsea(i) = 'g'

     if (index(nl%plotspecs(i)%pltspec2,'D') > 0) op%pltdualgrid(i) = 'D'

     if (index(nl%plotspecs(i)%pltspec2,'b') > 0) op%pltborder(i) = 'b'
     if (index(nl%plotspecs(i)%pltspec2,'t') > 0) op%pltborder(i) = 't'

     if (index(nl%plotspecs(i)%pltspec2,'n') > 0) op%labelbar(i) = 'n'
     if (index(nl%plotspecs(i)%pltspec2,'i') > 0) op%labelbar(i) = 'i'

     if (index(nl%plotspecs(i)%pltspec2,'c') > 0) op%colorbar(i) = 'c'

     if (index(nl%plotspecs(i)%pltspec2,'M') > 0) op%maptyp(i) = 'M'
     if (index(nl%plotspecs(i)%pltspec2,'m') > 0) op%maptyp(i) = 'm'

     if (index(nl%plotspecs(i)%pltspec2,'C') > 0) op%pltcone(i) = 'C'

     if (index(nl%plotspecs(i)%pltspec2,'p') > 0) op%pltlev(i) = 'p'
     if (index(nl%plotspecs(i)%pltspec2,'s') > 0) op%pltlev(i) = 's'

     if (index(nl%plotspecs(i)%pltspec2,'W') > 0) op%windowin(i) = 'W'

     if (index(nl%plotspecs(i)%pltspec2,'e') > 0) op%ext(i) = 'e'

     if (index(nl%plotspecs(i)%pltspec2,'f') > 0) op%frameoff(i) = 'f'
     if (index(nl%plotspecs(i)%pltspec2,'h') > 0) op%frameoff(i) = 'h'

     if (index(nl%plotspecs(i)%pltspec2,'1') > 0) op%panel(i) = '1'
     if (index(nl%plotspecs(i)%pltspec2,'2') > 0) op%panel(i) = '2'
     if (index(nl%plotspecs(i)%pltspec2,'3') > 0) op%panel(i) = '3'
     if (index(nl%plotspecs(i)%pltspec2,'4') > 0) op%panel(i) = '4'
     if (index(nl%plotspecs(i)%pltspec2,'5') > 0) op%panel(i) = '5'
     if (index(nl%plotspecs(i)%pltspec2,'6') > 0) op%panel(i) = '6'
     if (index(nl%plotspecs(i)%pltspec2,'7') > 0) op%panel(i) = '7'
     if (index(nl%plotspecs(i)%pltspec2,'8') > 0) op%panel(i) = '8'
     if (index(nl%plotspecs(i)%pltspec2,'9') > 0) op%panel(i) = '9'

     if (index(nl%plotspecs(i)%pltspec2,'u') > 0) op%noundrg(i) = 'u'
  enddo

  do i = 1,op%nplt_files
     op%plt_files(i) = nl%plt_files(i)
  enddo

! Variables in the following section specify the configuration of added
! grids (local mesh refinements) and must be copied from the namelist when
! a history file is not being read or when the namelist value is intended
! to replace the history file values.

  if (runtype == 'MAKEGRID'     .or. &
      runtype == 'INITIAL'      .or. &
      runtype == 'HISTADDGRID') then

     ngrids     = nl%ngrids
     ngrids_old = nl%ngrids_old

     do i = 1,maxgrds
        ngrdll(i) = nl%ngrdll(i)

        do j = 1,ngrdll(i)
           grdrad(i,j) = nl%grdrad(i,j)
           grdlat(i,j) = nl%grdlat(i,j)
           grdlon(i,j) = nl%grdlon(i,j)
        enddo
     enddo

  endif

! Variables in the following section either must not be changed on a history
! restart or changing them would be irrelevant.  Thus, they are only copied
! from the namelist if a history file is not being read.

  if (runtype == 'MAKEGRID' .or. runtype == 'INITIAL') then

     itime1    = nl%itime1
     idate1    = nl%idate1
     imonth1   = nl%imonth1
     iyear1    = nl%iyear1
     nzp       = nl%nzp
     nzg       = nl%nzg
     nzs       = nl%nzs
     nxp       = nl%nxp
     mdomain   = nl%mdomain
     isfcl     = nl%isfcl
     itopoflg  = nl%itopoflg
     ivegflg   = nl%ivegflg
     iseagrid  = nl%iseagrid
     deltax    = nl%deltax
     ndz       = nl%ndz
   
     hdz(1:ndz) = nl%hdz(1:ndz)
     dz (1:ndz) = nl%dz (1:ndz)

     slz(1:nzgmax) = nl%slz(1:nzgmax)

     ! set current time to initial time here.  If this is a history run,
     ! reset current time in subroutine history_start.
     current_time%year = iyear1
     current_time%month = imonth1
     current_time%date = idate1
     current_time%time = int(itime1 * 0.01) * 3600.0 &
          + (itime1 * 0.01 - int(itime1*0.01))*100.0*60.0

  endif

end subroutine copy_nl

!===============================================================================

subroutine namelist_print()

use oname_coms, only: nl
use misc_coms,  only: io6

implicit none

namelist /OLAMIN/ nl

write(io6, nml=OLAMIN)

end subroutine namelist_print
