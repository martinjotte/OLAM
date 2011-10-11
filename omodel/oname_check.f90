!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

subroutine oname_check()

! THIS ROUTINE CHECKS THE OPTION SPECIFICATIONS OF THE NAMELIST
! FILE OLAMIN FOR CONSISTENCY, AND OVERRIDES THE SETTINGS OF ICLOUD,
! IDRIZ, IRAIN, IPRIS, ISNOW, IAGGR, IGRAUP, IHAIL, AND ICCNLEV, SETTING THEM
! ALL TO ZERO IF LEVEL IS LESS THAN 3.

! EACH ERROR WILL BE ASSIGNED A SEVERITY.  FATAL ERRORS WILL CAUSE
! THE RUN TO STOP IMMEDIATELY AND WARNING ERRORS WILL BE LISTED.

use max_dims,    only: nzgmax, maxgrds, maxsndg, maxnplt, maxisdirs, &
                       maxpltfiles, maxngrdll
use oname_coms,  only: nl
use consts_coms, only: erad, pi1, pi2, p00
use misc_coms,   only: io6

implicit none

integer      :: nfatal, nwarn, ifm, ng, i, k, iplt, i_huge, izaux
real         :: r_huge, r_min, r_max, r_tiny, dzxmin, zb_min, zb_max
real(kind=8) :: d_huge

character(len=1), dimension(8) :: tunits = (/ 'd','D','h','H','m','M','s','S' /)

nfatal = 0
nwarn  = 0

r_tiny = spacing(1.0)       ! real value small compared to 1.0
r_huge = huge(1.0) - 1.0    ! largest allowable real value
i_huge = huge(1) - 1        ! largest allowable integer
d_huge = huge(nl%timmax) - 1.d0 ! largest allowable real*8 value
dzxmin = 0.001              ! minimum grid spacing (horiz. and vertical)

zb_min = -2.e3     ! min. vertical base (allow gridpoints below MSL)
if (nl%mdomain == 0) then
   zb_max = 0.0    ! global run - vertical base at MSL
else
   zb_max = 2.e3  ! max. vertical base (allow above MSL for a regional run)
endif

write(io6,*) ' '
write(io6,*) '  CHECKING NAMELIST VALUES'
write(io6,*) ' '

!--------------------------------------------------------------------------
! RUNTYPE
!--------------------------------------------------------------------------

if ( (nl%runtype /= 'MAKESFC'   ) .and. &
     (nl%runtype /= 'MAKEGRID'  ) .and. &
     (nl%runtype /= 'INITIAL'   ) .and. &
     (nl%runtype /= 'HISTORY'   ) .and. &
     (nl%runtype /= 'PLOTONLY'  ) .and. &
     (nl%runtype /= 'PARCOMBINE') ) then
   write(io6,*) " -- FATAL -- RUNTYPE = "//trim(nl%runtype)
   write(io6,*) "             RUNTYPE must be either 'MAKESFC', 'MAKEGRID',"
   write(io6,*) "             'INITIAL', 'HISTORY', 'PLOTONLY', or"
   write(io6,*) "             'PARCOMBINE'"
   nfatal = nfatal + 1
endif

!--------------------------------------------------------------------------
! SIMULATION ENDING TIME
!--------------------------------------------------------------------------

if (.not. any(nl%timeunit == tunits)) then
   write(io6,*) " -- FATAL -- TIMEUNIT = "//nl%timeunit
   write(io6,*) "             TIMEUNIT must be either 'd', 'h', 'm', or 's'"
   nfatal = nfatal + 1
endif

call dchk_bnds( nl%timmax, "TIMMAX", 0.d0,  d_huge, 0, nfatal, nwarn )

!--------------------------------------------------------------------------
! START OF SIMULATION OR ISAN PROCESSING
!--------------------------------------------------------------------------

call ichk_bnds( nl%itime1,   "ITIME1",       0,    2359, 0, nfatal, nwarn )
call ichk_bnds( nl%idate1,   "IDATE1",       1,      31, 0, nfatal, nwarn )
call ichk_bnds( nl%imonth1, "IMONTH1",       1,      12, 0, nfatal, nwarn )
call ichk_bnds( nl%iyear1,   "IYEAR1",       0,    9999, 0, nfatal, nwarn )

!--------------------------------------------------------------------------
! GRID SPECIFICATIONS
!--------------------------------------------------------------------------

call ichk_bnds( nl%mdomain, "MDOMAIN",       0,       4, 0, nfatal, nwarn )
call ichk_bnds( nl%meshtype,"MESHTYPE",      1,       2, 0, nfatal, nwarn )

call ichk_bnds( nl%ngrids,   "NGRIDS",       1, maxgrds, 0, nfatal, nwarn, &
     msgmax="Increase maxgrds in max_dims.f90 if more nests are needed." )

call ichk_bnds( nl%nzp,         "NZP",       3,   10000, 0, nfatal, nwarn, &
     msgmin="At least 3 vertical levels are needed for OLAM." )

call ichk_bnds( nl%nxp,         "NXP",       1,   10000, 0, nfatal, nwarn )
call rchk_bnds( nl%dtlong,   "DTLONG",  r_tiny,  r_huge, 0, nfatal, nwarn )
call rchk_bnds( nl%deltax,   "DELTAX",  dzxmin,  r_huge, 0, nfatal, nwarn )
call rchk_bnds( nl%deltaz,   "DELTAZ",      0.,  r_huge, 0, nfatal, nwarn )
call rchk_bnds( nl%zbase,    "ZBASE",   zb_min,  zb_max, 0, nfatal, nwarn )
call rchk_bnds( nl%dzbase,   "DZBASE",     0.0,  r_huge, 0, nfatal, nwarn )
call rchk_bnds( nl%ztop,     "ZTOP",       0.0,  r_huge, 0, nfatal, nwarn )
call rchk_bnds( nl%dztop,    "DZTOP",      0.0,  r_huge, 0, nfatal, nwarn )
call ichk_bnds( nl%nzaux,    "NZAUX",       -1,      10, 0, nfatal, nwarn )

do izaux=1, nl%nzaux
   call rchk_bnds( nl%zaux(izaux),   "ZAUX",   0.0,  r_huge, 0, nfatal, nwarn )
   call rchk_bnds( nl%dzaux(izaux),  "DZAUX",  0.0,  r_huge, 0, nfatal, nwarn )

   if (nl%zaux(izaux) <= nl%zbase) then
      write(io6,*) 'FATAL - zaux(',izaux,') must be larger than zbase.'
      nfatal = nfatal + 1
   endif

   if (nl%zaux(izaux) >= nl%ztop) then
      write(io6,*) 'FATAL - zaux(',izaux,') must be less than ztop.'
      nfatal = nfatal + 1
   endif

   if (izaux > 1 .and. nl%zaux(izaux) <= nl%zaux(izaux-1) ) then
      write(io6,*) 'FATAL - zaux(',izaux,') must be larger than zaux(',izaux-1,').'
      nfatal = nfatal + 1
   endif
enddo

if (nl%deltaz >= dzxmin) then

   call rchk_bnds( nl%dzrat,  "DRZAT",     1.0,    10.0, 0, nfatal, nwarn )
   call rchk_bnds( nl%dzrat,  "DRZAT",     0.0,     1.2, 1, nfatal, nwarn, &
        msgmax="Large DZRAT degrades 2nd-order accuracy in the vertical &
              & differencing")
   call rchk_bnds( nl%dzmax,  "DZMAX",  dzxmin,  r_huge, 0, nfatal, nwarn )

else

   call rchk_bnds( nl%zz(1),     "ZZ",  zb_min,  zb_max, 0, nfatal, nwarn )
   do k=2,nl%nzp
      r_min = nl%zz(k-1) + dzxmin
      call rchk_bnds( nl%zz(k),  "ZZ",   r_min,  r_huge, 0, nfatal, nwarn )
   enddo

endif

!--------------------------------------------------------------------------
! NESTED GRID DEFINITION
!--------------------------------------------------------------------------

do ng=2, nl%ngrids
   call ichk_bnds(nl%ngrdll(ng),  "NGRDLL",    1, maxngrdll, 0, nfatal, nwarn )
   call rchk_bnds(nl%grdrad(ng),  "GRDRAD", dzxmin, erad*2., 0, nfatal, nwarn )
enddo

if (nl%mdomain < 2) then

   do i = 1,nl%ngrdll(ng)
      do ng=2, nl%ngrids
         call rchk_bnds( nl%grdlat(ng,i), "GRDLAT",  -90.,  90., 0, nfatal, nwarn )
         call rchk_bnds( nl%grdlon(ng,i), "GRDLON", -180., 180., 0, nfatal, nwarn )
      enddo
   enddo

endif

!--------------------------------------------------------------------------
! TIMESTEP RATIOS
!--------------------------------------------------------------------------

do ng=1, nl%ngrids

 ! call ichk_bnds( nl%ndtrat(ng), "NDTRAT", 1,  2, 0, nfatal, nwarn ) ! future
   call ichk_bnds( nl%ndtrat(ng), "NDTRAT", 1,  1, 2, nfatal, nwarn,   &
        msgboth="NDTRAT must = 1 for all MRLs in this version of OLAM.")

   call ichk_bnds( nl%nacoust(ng), "NACOUST", 1, 20, 0, nfatal, nwarn )
enddo

!--------------------------------------------------------------------------
! VARIABLE INITIALIZATION INPUT
!--------------------------------------------------------------------------

call ichk_bnds( nl%initial,   "INITIAL", 1,  3, 0, nfatal, nwarn )

!--------------------------------------------------------------------------
! NUDGING PARAMETERS
!--------------------------------------------------------------------------

if (nl%initial == 2) then

   call ichk_bnds( nl%nudflag, "NUDFLAG", 0,  1, 2, nfatal, nwarn )

   if (nl%nudflag == 1) then

      call ichk_bnds( nl%nudnxp, "NUDNXP", 0, 10000, 2, nfatal, nwarn )
      call rchk_bnds( nl%tnudcent, "TNUDCENT", nl%dtlong, r_huge, 2, nfatal, &
                      nwarn, msgmin="Nudging time must be larger than dtlong")
   endif

endif

!--------------------------------------------------------------------------
! HISTORY FILE OUTPUT
!--------------------------------------------------------------------------

call ichk_bnds( nl%ioutput,   "IOUTPUT", 0,  1, 2, nfatal, nwarn )
call ichk_bnds( nl%iclobber, "ICLOBBER", 0,  1, 2, nfatal, nwarn )
call rchk_bnds( nl%frqstate, "FRQSTATE", nl%dtlong,  r_huge, 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! Topography
!--------------------------------------------------------------------------

call ichk_bnds( nl%itopoflg, "ITOPOFLG", 1,  2, 0, nfatal, nwarn )

!--------------------------------------------------------------------------
! MODEL OPTIONS / NUMERICAL SCHEMES
!--------------------------------------------------------------------------

call ichk_bnds( nl%naddsc,   "NADDSC", 0, 1000, 0, nfatal, nwarn )
call ichk_bnds( nl%icorflg, "ICORFLG", 0,    1, 0, nfatal, nwarn )

!--------------------------------------------------------------------------
! RAYLEIGH FRICTION PARAMETERS
!--------------------------------------------------------------------------

call rchk_bnds( nl%rayf_distim, "RAYF_DISTIM", 0.0, r_huge, 0, nfatal, nwarn )
if (nl%rayf_distim > r_tiny) &
     call rchk_bnds( nl%rayf_distim, "RAYF_DISTIM", nl%dtlong, r_huge, 2, &
                     nfatal, nwarn )
call rchk_bnds( nl%rayf_expon, "RAYF_EXPON", 0.0,   5.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayf_zmin,   "RAYF_ZMIN", 0.0, r_huge, 2, nfatal, nwarn )

call rchk_bnds( nl%rayfw_distim, "RAYFW_DISTIM", 0.0, r_huge, 2, nfatal, nwarn )
if (nl%rayfw_distim > r_tiny) &
     call rchk_bnds( nl%rayfw_distim, "RAYFW_DISTIM", nl%dtlong, r_huge, 2, &
     nfatal, nwarn )
call rchk_bnds( nl%rayfw_expon, "RAYFW_EXPON", 0.0,    5.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfw_zmin,   "RAYFW_ZMIN", 0.0, r_huge, 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! RADIATION PARAMETERIZATION PARAMETERS
!--------------------------------------------------------------------------

call ichk_bnds( nl%iswrtyp, "ISWRTYP", 0,  3, 0, nfatal, nwarn )
call ichk_bnds( nl%ilwrtyp, "ILWRTYP", 0,  3, 0, nfatal, nwarn )
call rchk_bnds( nl%radfrq,   "RADFRQ", nl%dtlong, r_huge, 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! CUMULUS PARAMETERIZATION PARAMETERS
!--------------------------------------------------------------------------

do ng=1, nl%ngrids
   call ichk_bnds( nl%nqparm(ng), "NQPARM", 0, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%nqparm_sh(ng), "NQPARM_SH", 0, 1, 0, nfatal, nwarn )
enddo
call rchk_bnds( nl%confrq, "CONFRQ", nl%dtlong, r_huge, 2, nfatal, nwarn )
call rchk_bnds( nl%wcldbs, "WCLDBS", -10., 10., 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! EDDY DIFFUSION PARAMETERS
!--------------------------------------------------------------------------

do ng=1, nl%ngrids
   call ichk_bnds( nl%idiffk(ng), "IDIFFK", 0,     7, 0, nfatal, nwarn )
   call rchk_bnds( nl%csx(ng),    "CSX",    0.,  10., 0, nfatal, nwarn )
   call rchk_bnds( nl%csz(ng),    "CSZ",    0.,  10., 0, nfatal, nwarn )
   call rchk_bnds( nl%zkhkm(ng),  "ZKHKM",  0., 100., 0, nfatal, nwarn )
   call rchk_bnds( nl%xkhkm(ng),  "XKHKM",  0., 100., 0, nfatal, nwarn )
   call rchk_bnds( nl%akmin(ng),  "AKMIN",  0.,  10., 0, nfatal, nwarn )
enddo

!--------------------------------------------------------------------------
! MICROPHYSICS PARAMETERS
!--------------------------------------------------------------------------

call ichk_bnds( nl%level, "LEVEL", 0, 3, 0, nfatal, nwarn )

if (nl%level <= 2) then
   nl%icloud  = 0
   nl%idriz   = 0
   nl%irain   = 0
   nl%ipris   = 0
   nl%isnow   = 0
   nl%iaggr   = 0
   nl%igraup  = 0
   nl%ihail   = 0
   nl%iccnlev = 0

elseif (nl%level == 3) then

   if (.not. any(nl%icloud == (/0, 1, 4, 5, 6, 7/) )) then
      write(io6,*) 'FATAL - icloud must be set to 0, 1, 4, 5, 6, or 7.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%idriz == (/0, 1, 4, 5, 6, 7/) )) then
      write(io6,*) 'FATAL - idriz must be set to 0, 1, 4, 5, 6, or 7.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%irain == (/0, 1, 2, 5/) )) then
      write(io6,*) 'FATAL - irain must be set to 0, 1, 2, or 5.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%ipris == (/0, 5, 6, 7/) )) then
      write(io6,*) 'FATAL - ipris must be set to 0, 5, 6, or 7.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%isnow == (/0, 1, 2, 5/) )) then
      write(io6,*) 'FATAL - isnow must be set to 0, 1, 2, or 5.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%iaggr == (/0, 1, 2, 5/) )) then
      write(io6,*) 'FATAL - iaggr must be set to 0, 1, 2, or 5.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%igraup == (/0, 1, 2, 5/) )) then
      write(io6,*) 'FATAL - igraup must be set to 0, 1, 2, or 5.'
      nfatal = nfatal + 1
   endif
   
   if (.not. any(nl%ihail == (/0, 1, 2, 5/) )) then
      write(io6,*) 'FATAL - ihail must be set to 0, 1, 2, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%iccnlev == (/0, 1/) )) then
      write(io6,*) 'FATAL - iccnlev must be set to 0 or 1.'
      nfatal = nfatal + 1
   endif

   if (any(nl%icloud == (/4, 5/) )) &   
      call rchk_bnds( nl%cparm, "CPARM", 0., 1.e10, 0, nfatal, nwarn )

   if (any(nl%idriz == (/4, 5/) )) &   
      call rchk_bnds( nl%dparm, "DPARM", 0.,  1.e3, 0, nfatal, nwarn )

   if (nl%irain == 2) &
      call rchk_bnds( nl%rparm, "RPARM", 0., 1.e-2, 0, nfatal, nwarn )

   if (nl%isnow == 2) &
      call rchk_bnds( nl%sparm, "SPARM", 0., 1.e-2, 0, nfatal, nwarn )

   if (nl%iaggr == 2) &
      call rchk_bnds( nl%aparm, "APARM", 0., 1.e-2, 0, nfatal, nwarn )

   if (nl%igraup == 2) &
      call rchk_bnds( nl%gparm, "GPARM", 0., 1.e-2, 0, nfatal, nwarn )

   if (nl%ihail == 2) &
      call rchk_bnds( nl%hparm, "HPARM", 0., 3.e-2, 0, nfatal, nwarn )

   call rchk_bnds( nl%cnparm, "CNPARM", 0., 1.e-7, 0, nfatal, nwarn )
   call rchk_bnds( nl%gnparm, "GNPARM", 0., 1.e-5, 0, nfatal, nwarn )

endif

!--------------------------------------------------------------------------
! SOUNDING SPECIFICATION
!--------------------------------------------------------------------------

if (nl%initial == 1) then

   call ichk_bnds( nl%nsndg, "NSNDG", 2, maxsndg, 0, nfatal, nwarn, &
        msgmin="Sounding must have at least 2 levels.", &
        msgmax="Increase maxsndg in max_dims.f90 if more levels are needed." )

   call ichk_bnds( nl%ipsflg,  "IPSFLG",  0, 1, 0, nfatal, nwarn )
   call ichk_bnds( nl%itsflg,  "ITSFLG",  0, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%irtsflg, "IRTSFLG", 0, 4, 0, nfatal, nwarn )
   call ichk_bnds( nl%iusflg,  "IUSFLG",  0, 1, 0, nfatal, nwarn )

   if (nl%ipsflg == 0) &
        call rchk_bnds( nl%hs, "HS", -1.e3, 1.e5, 0, nfatal, nwarn )

   if (nl%ipsflg == 1) then
      r_min = 0.2 * p00/100.0 
      r_max = 1.5 * p00/100.0 
      call rchk_bnds( nl%p_sfc, "P_SFC", r_min, r_max, 0, nfatal, nwarn )
   endif

endif

!--------------------------------------------------------------------------
! SOUNDING
!--------------------------------------------------------------------------

! if (nl%initial == 1) then
!    todo: add checks for input sounding
! endif

!--------------------------------------------------------------------------
!  LEAF VARIABLES
!--------------------------------------------------------------------------

call ichk_bnds( nl%isfcl, "ISFCL", 0, 1, 0, nfatal, nwarn )

if (nl%isfcl == 1) then

   call ichk_bnds( nl%nzg, "NZG", 2, nzgmax, 0, nfatal, nwarn, &
        msgmin="At least 2 soil levels are needed for soil model.", &
        msgmax="Increase nzgmax in max_dims.f90 if more soil layers are needed")

   call ichk_bnds( nl%nzs, "NZS", 0, 10, 0, nfatal, nwarn )

   call ichk_bnds( nl%ivegflg,     "IVEGFLG",     1, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%isoilflg,    "ISOILFLG",    1, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%ndviflg,     "NDVIFLG",     1, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%isstflg,     "ISSTFLG",     0, 2, 0, nfatal, nwarn )
   CALL ichk_bnds( nl%iseaiceflg,  "ISEAICEFLG",  0, 2, 0, nfatal, nwarn )

   call ichk_bnds( nl%isoilstateinit, "ISOILSTATEINIT", 0, 1, 0, nfatal,nwarn )
   call ichk_bnds( nl%isoildepthflg,  "ISOILDEPTHFLG",  0, 1, 0, nfatal,nwarn )

   call ichk_bnds( nl%iupdndvi,  "IUPDNDVI",    0, 1, 0, nfatal, nwarn )
   call ichk_bnds( nl%iupdsst,   "IUPDSST",     0, 1, 0, nfatal, nwarn )
   call ichk_bnds( nl%iupdseaice,"IUPDSEAICE",  0, 1, 0, nfatal, nwarn )

   if (nl%ivegflg == 2) &
        call ichk_bnds( nl%nvgcon, "NVGCON", 0, 20, 0, nfatal, nwarn )
   
   if (nl%isoilflg == 2) &
        call ichk_bnds( nl%nslcon, "NSLCON", 1, 12, 0, nfatal, nwarn )

   if (nl%isstflg == 2) &
        call rchk_bnds( nl%seatmp, "SEATMP", 100., 500., 0, nfatal, nwarn )

   call rchk_bnds( nl%slz(nl%nzg), "SLZ", -r_huge, -0.02, 0, nfatal, nwarn )
   do k=nl%nzg-1, 1, -1
      r_max =  nl%slz(k+1) - 0.02
      call rchk_bnds( nl%slz(k), "SLZ", -r_huge, -0.02, 0, nfatal, nwarn )
   enddo

   do k=1,nl%nzg
      call rchk_bnds( nl%slmstr(k), "SLMSTR", 0., 1., 0, nfatal, nwarn )
   enddo

   ! Set iupdndvi to 0 if not reading NDVI files
   if (nl%ndviflg == 2) then
      if (nl%iupdndvi == 1) then
         write(io6,*) "Turning off NDVI update for ndviflg == 2"
         nl%iupdndvi = 0
         nwarn = nwarn + 1
      endif
   endif

   ! Set iupdsst to 0 if not reading SST files
   if (nl%isstflg == 2) then
      if (nl%iupdsst == 1) then
         write(io6,*) "Turning off SST update for isstflg == 2"
         nl%iupdsst = 0
         nwarn = nwarn + 1
      endif
   endif

   ! Set iupdseaice to 0 if not reading SEAICE files
   if (nl%iseaiceflg == 2) then
      if (nl%iupdseaice == 1) then
         write(io6,*) "Turning off SEAICE update for iseaiceflg == 2"
         nl%iupdseaice = 0
         nwarn = nwarn + 1
      endif
   endif

endif
      
!--------------------------------------------------------------------------
! ISENTROPIC CONTROL 
!--------------------------------------------------------------------------

CALL ichk_bnds( nl%isdirs,   "ISANDIRS", -1, maxisdirs, 0, nfatal, nwarn, &
     msgmax="Increase maxisdirs in maxdims.f90 if you need to read from more &
          & directories" )

!--------------------------------------------------------------------------
! PLOTTING PARAMETERS
!--------------------------------------------------------------------------

call ichk_bnds( nl%nplt, "NPLT", 0, maxnplt, 0, nfatal, nwarn )

if (nl%runtype == 'PLOTONLY') &
     call ichk_bnds( nl%nplt_files, "NPLT_FILES", 0, maxpltfiles, 0, nfatal, &
     nwarn, msgmax="Increase maxpltfiles in maxdims.f90 if you need to plot &
     & from more files" )

if ((nl%runtype == 'INITIAL') .or. (nl%runtype == 'HISTORY')) &
     call rchk_bnds( nl%frqplt, "FRQPLT", nl%dtlong, r_huge, 2, nfatal, nwarn )

call rchk_bnds( nl%dtvec,     "DTVEC",     1.e-3, r_huge, 2, nfatal, nwarn )
call rchk_bnds( nl%headspeed, "HEADSPEED", 1.e-3, r_huge, 2, nfatal, nwarn )
call rchk_bnds( nl%stemlength, "STEMLENGTH",1.e-3, r_huge, 2, nfatal, nwarn )
call ichk_bnds( nl%plttype,   "PLTTYPE",   0, 2, 2, nfatal, nwarn )
call ichk_bnds( nl%pltorient, "PLTORIENT", 0, 1, 2, nfatal, nwarn )
      
!------------------------------------------------------------------------
! PLOTTING SPECIFICATIONS
!------------------------------------------------------------------------

! If there are less than NPLT specifications, reset NPLT
do iplt = 1, nl%nplt
   if (nl%plotspecs(iplt)%fldname == '') exit
enddo
iplt = iplt - 1
if (iplt < nl%nplt) then
   write(io6, '(A,I0,A)') " Only ", iplt, " plotting specifications found."
   write(io6, '(A,I0,A)') " Resetting NPLT to ", iplt, "."
   nwarn = nwarn + 1
   nl%nplt = iplt
endif

! Make sure that plot projections are compatible with MDOMAIN value
do iplt = 1, nl%nplt
   if ((nl%plotspecs(iplt)%projectn == 'V'  .or.   &
        nl%plotspecs(iplt)%projectn == 'Z') .and.  &
        nl%mdomain < 2) then
        
        write(io6,*) 'When MDOMAIN < 2, may not use V or Z plot projection'
        nfatal = nfatal + 1
   endif
   
   if ((nl%plotspecs(iplt)%projectn == 'L'  .or.   &
        nl%plotspecs(iplt)%projectn == 'P'  .or.   &
        nl%plotspecs(iplt)%projectn == 'O'  .or.   &
        nl%plotspecs(iplt)%projectn == 'C') .and.  &
        nl%mdomain > 1) then
        
        write(io6,*) 'When mdomain > 1, may not use L, P, O, or C plot projection'
        nfatal = nfatal + 1
   endif
enddo    

! TODO - add checks for the plotting variables

!------------------------------------------------------------------------
! OTHER CONSTANCY CHECKS
!------------------------------------------------------------------------

! TEMPORARY!! ONLY ALLOW MDOMAIN VALUES OF 0, 3, and 4

if (nl%mdomain /= 0 .and. nl%mdomain /= 3 .and. nl%mdomain /= 4) then
   write(io6,*) ' FATAL - MDOMAIN temporarily restricted to values of 0, 3, or 4.'
   nfatal = nfatal + 1
endif

! TEMPORARY!! NACOUST MUST BE SAME ON ALL MESHES FOR THIS VERSION

if (nl%ngrids > 1) then
   if (any(nl%nacoust(2:nl%ngrids) /= nl%nacoust(1))) then
      write(io6,*) 'FATAL - NACOUST must be at least the same for each mesh'
      write(io6,*) 'refinement level IN THIS VERSION OF OLAM.'
      nfatal = nfatal + 1
   endif
endif

! IF THIS IS A GLOBAL SIMULATION, PRINT MESSAGE THAT DELTAX WILL BE
! REDEFINED BY NXP.

if (nl%mdomain == 0 .and. nfatal ==0) then

   write(io6,*) ' '
   write(io6,*) 'Since MDOMAIN is set to 0, this is configured as a global&
          & simulation.'
   write(io6,*) 'Hence, DELTAX is ignored and horizontal grid spacing will be&
          & computed'
   write(io6,*) 'from nxp and the Earth (or other planetary) radius.'
   write(io6,*) ' '

endif

! RUN WITH LONGITUDINALLY HOMOGENEOUS INITIALIZATION MUST BE GLOBAL

if (nl%initial == 3 .and. nl%mdomain /= 0) then
   write(io6,*) ' FATAL  - mdomain must be 0 if INITIAL = 3.'
   nfatal = nfatal + 1
endif
  
! CONVECTIVE PARAMETERIZATION MUST HAVE AT LEAST WATER VAPOR

if (any(nl%nqparm(1:nl%ngrids) > 0) .and. nl%level == 0) then
   write(io6,*) ' FATAL - LEVEL must be at least 1 for the cumulus parameterization'
   nfatal = nfatal + 1
endif

! MAKE SURE THAT RADIATION SCHEMES ACTIVATED IF SURFACE MODEL USED

if (nl%isfcl > 0 .and. (nl%ilwrtyp == 0 .or. nl%iswrtyp == 0)) then
   write(io6,'(A)') ' FATAL - longwave and shortwave radiation schemes must be&
        & activated when using the surface model.'
   nfatal = nfatal + 1
endif

! MAKE SURE THAT SURFACE MODEL IS ACTIVATED IF RADIATION SCHEMES ARE USED

if (nl%isfcl == 0 .and. (nl%ilwrtyp /= 0 .or. nl%iswrtyp /= 0)) then
   write(io6,'(A)') ' FATAL - surface model must be activated when using&
        & longwave or shortwave radiation schemes.'
   nfatal = nfatal + 1
endif

! MUST HAVE WATER VAPOR IF USING VARIABLE INITIALIZATION

if (nl%initial == 2 .and. nl%level == 0) then
   write(io6,*) 'FATAL - Moisture complexity LEVEL must be 1 or larger if'
   write(io6,*) 'variable initialization is used (i.e., if INITIAL = 2).'
   nfatal = nfatal + 1
endif

! ALWAYS SET IOUTPUT=1 FOR A PARCOMBINE RUN

if (nl%runtype == 'PARCOMBINE' .and. nl%ioutput == 0) then
   write(io6,*) ' WARNING - setting IOUTPUT=1 for a PARCOMBINE run'
   nwarn = nwarn + 1
endif

! STOP THE RUN IF THERE ARE ANY FATAL ERRORS, AND LIST HOW MANY
! FATAL AND WARNING MESSAGES

write(io6,*) ''
write(io6,*) ' ----------oname_check--------------------------'
write(io6,*) ' FATAL     errors - ', nfatal
write(io6,*) ' WARNING   errors - ', nwarn
write(io6,*) ' -----------------------------------------------'
write(io6,*) ''
if (nfatal > 0) stop 'ONAME_CHECK'


contains
  
  
  subroutine ichk_bnds(ivar, name, iminv, imaxv, iflag, nfatal, nwarn,  &
                       msgmin, msgmax, msgboth)
    implicit none
    integer,       intent(inout) :: ivar, nfatal, nwarn
    integer,          intent(in) :: imaxv, iminv, iflag
    character(len=*), intent(in) :: name
    character(len=*),   optional :: msgmin, msgmax, msgboth

    ! CHECK THE BOUNDS OF AN INTEGER NAMELIST VARIABLE. IF THE BOUNDS ARE
    ! EXCEEDED, THEN IT RESPONDS BASED ON THE VALUE OF IFLAG.
    !
    ! IFLAG:
    !  0 - REPORT A FATAL ERROR IF BOUNDS EXCEEDED, RUN STOPS
    !  1 - REPORT A WARNING IF BOUNDS EXCEEDED, RUN CONTINUES
    !  2 - IF BOUNDS EXCEEDED, RESET THE VARIABLE TO WITHIN THE GIVEN BOUNDS
    !       AND REPORT A WARNING, RUN CONTINUES
    !
    ! ALSO PRINT THE OPTIONAL STRINGS MSGMIN, MSGMAX OR MSGBOTH IF THE
    ! MINIMUM AND/OR MAXIMUM LIMITS ARE EXCEEDED

    if (iflag .eq. 0) then

       if ((ivar < iminv) .or. (ivar > imaxv)) then
          nfatal = nfatal + 1
          write(io6, "(A,I0)") &
               '-- FATAL -- Input variable '//trim(name)//' = ', ivar
          write(io6, "(A,I0,A,I0)") &
               '            Allowable range is ', iminv, ' to ', imaxv
       endif

    elseif (iflag .eq. 1) then

       if ((ivar < iminv) .or. (ivar > imaxv)) then
          nwarn = nwarn + 1
          write(io6, "(A,I0)") &
               '-- WARNING -- Input variable '//trim(name)//' = ', ivar
          write(io6, "(A,I0,A,I0/)") &
               '              Acceptable range is ', iminv, ' to ', imaxv
       endif

    elseif (iflag .eq. 2) then

       if (ivar < iminv) then

          write(io6, "(A,I0)") &
               '-- WARNING -- Input variable '//trim(name)//' = ', ivar
          write(io6, "(A,I0)") &
               '              is less than the minimum acceptable value ', iminv
          write(io6, "(A,I0,A/)") &
               '   Will set '//trim(name)//' = ', iminv, ' for this run.'
          ivar = iminv
          nwarn = nwarn + 1

       elseif (ivar > imaxv) then

          write(io6, "(A,I0)") &
               '-- WARNING -- Input variable '//trim(name)//' = ', ivar
          write(io6, "(A,I0)") &
               '              is larger than the maximum acceptable value ', imaxv
          write(io6, "(A,I0,A/)") &
               '   Will set '//trim(name)//' = ', imaxv, ' for this run.'
          ivar = imaxv
          nwarn = nwarn + 1

       endif
    endif

    if (present(msgmin) .and. (ivar < iminv)) write(io6,*) msgmin
    if (present(msgmax) .and. (ivar > imaxv)) write(io6,*) msgmax
    if (present(msgboth) .and. ((ivar < iminv) .or. (ivar > imaxv))) &
         write(io6,*) msgboth
    if ((ivar < iminv) .or. (ivar > imaxv)) write(io6, *) ''

  end subroutine ichk_bnds



  subroutine rchk_bnds(rvar, name, rminv, rmaxv, iflag, nfatal, nwarn, &
                       msgmin, msgmax, msgboth)
    implicit none
    integer,       intent(inout) :: nfatal, nwarn
    integer,       intent(in)    :: iflag
    real,          intent(inout) :: rvar
    real,          intent(in)    :: rmaxv, rminv

    character(len=*), intent(in) :: name
    character(len=*),   optional :: msgmin, msgmax, msgboth
    character(len=6)             :: varformat, minformat, maxformat

    ! CHECK THE BOUNDS OF A REAL NAMELIST VARIABLE. IF THE BOUNDS ARE
    ! EXCEEDED, THEN IT RESPONDS BASED ON THE VALUE OF IFLAG.
    !
    ! IFLAG:
    !  0 - REPORT A FATAL ERROR IF BOUNDS EXCEEDED, RUN STOPS
    !  1 - REPORT A WARNING IF BOUNDS EXCEEDED, RUN CONTINUES
    !  2 - IF BOUNDS EXCEEDED, RESET THE VARIABLE TO WITHIN THE GIVEN BOUNDS
    !       AND REPORT A WARNING, RUN CONTINUES
    !
    ! ALSO PRINT THE OPTIONAL STRINGS MSGMIN OR MSGMAX IF THE MINIMUM OR
    ! MAXIMUM LIMITS ARE EXCEEDED

    call rsetformat(varformat, rvar)
    call rsetformat(minformat, rminv)
    call rsetformat(maxformat, rmaxv)

    if (iflag .eq. 0) then

       if ((rvar < rminv) .or. (rvar > rmaxv)) then
          nfatal = nfatal + 1
          write(io6, "(A,"//varformat//")") &
               '-- FATAL -- Input variable '//trim(name)//' = ', rvar
          write(io6, "(A,"//minformat//",A,"//maxformat//")") &
               '            Allowable range is ', rminv, ' to ', rmaxv
       endif

    elseif (iflag .eq. 1) then

       if ((rvar < rminv) .or. (rvar > rmaxv)) then
          nwarn = nwarn + 1
          write(io6, "(A,I0)") &
               '-- WARNING -- Input variable '//trim(name)//' = ', rvar
          write(io6, "(A,"//minformat//",A,"//maxformat//")") &
               '              Acceptable range is ', rminv, ' to ', rmaxv
       endif

    elseif (iflag .eq. 2) then

       if (rvar < rminv) then

          write(io6, "(A,"//varformat//")") &
               '-- WARNING -- Input variable '//trim(name)//' = ', rvar
          write(io6, "(A,"//minformat//")") &
               '              is less than the minimum acceptable value ', rminv
          write(io6, "(A,"//minformat//",A/)") &
               '   Will set '//trim(name)//' = ', rminv, ' for this run.'
          rvar = rminv
          nwarn = nwarn + 1

       elseif (rvar > rmaxv) then

          write(io6, "(A,"//varformat//")") &
               '-- WARNING -- Input variable '//trim(name)//' = ', rvar
          write(io6, "(A,"//maxformat//")") &
               '              is larger than the maximum acceptable value ', rmaxv
          write(io6, "(A,"//maxformat//",A/)") &
               '   Will set '//trim(name)//' = ', rmaxv, ' for this run.'
          rvar = rmaxv
          nwarn = nwarn + 1

       endif
    endif

    if (present(msgmin) .and. (rvar < rminv)) write(io6,*) msgmin
    if (present(msgmax) .and. (rvar > rmaxv)) write(io6,*) msgmax
    if (present(msgboth) .and. ((rvar < rminv) .or. (rvar > rmaxv))) &
         write(io6,*) msgboth
    if ((rvar < rminv) .or. (rvar > rmaxv)) write(io6,*) ''

  end subroutine rchk_bnds



  subroutine rsetformat(fstring, var)
    character(len=*), intent(out) :: fstring
    real,             intent(in)  :: var
    if ((abs(var) >= 1.0e7) .or. (abs(var) <= 1.0e-4)) then
       fstring = "ES12.5"
    else
       fstring = "F0.7"
    endif
    
  end subroutine rsetformat



  subroutine dchk_bnds(dvar, name, dminv, dmaxv, iflag, nfatal, nwarn, &
                       msgmin, msgmax, msgboth)
    implicit none
    integer,       intent(inout) :: nfatal, nwarn
    integer,       intent(in)    :: iflag
    real(kind=8),  intent(inout) :: dvar
    real(kind=8),  intent(in)    :: dmaxv, dminv

    character(len=*), intent(in) :: name
    character(len=*),   optional :: msgmin, msgmax, msgboth
    character(len=6)             :: varformat, minformat, maxformat

    ! CHECK THE BOUNDS OF A REAL NAMELIST VARIABLE. IF THE BOUNDS ARE
    ! EXCEEDED, THEN IT RESPONDS BASED ON THE VALUE OF IFLAG.
    !
    ! IFLAG:
    !  0 - REPORT A FATAL ERROR IF BOUNDS EXCEEDED, RUN STOPS
    !  1 - REPORT A WARNING IF BOUNDS EXCEEDED, RUN CONTINUES
    !  2 - IF BOUNDS EXCEEDED, RESET THE VARIABLE TO WITHIN THE GIVEN BOUNDS
    !       AND REPORT A WARNING, RUN CONTINUES
    !
    ! ALSO PRINT THE OPTIONAL STRINGS MSGMIN OR MSGMAX IF THE MINIMUM OR
    ! MAXIMUM LIMITS ARE EXCEEDED

    call dsetformat(varformat, dvar)
    call dsetformat(minformat, dminv)
    call dsetformat(maxformat, dmaxv)

    if (iflag .eq. 0) then

       if ((dvar < dminv) .or. (dvar > dmaxv)) then
          nfatal = nfatal + 1
          write(io6, "(A,"//varformat//")") &
               '-- FATAL -- Input variable '//trim(name)//' = ', dvar
          write(io6, "(A,"//minformat//",A,"//maxformat//")") &
               '            Allowable range is ', dminv, ' to ', dmaxv
       endif

    elseif (iflag .eq. 1) then

       if ((dvar < dminv) .or. (dvar > dmaxv)) then
          nwarn = nwarn + 1
          write(io6, "(A,I0)") &
               '-- WARNING -- Input variable '//trim(name)//' = ', dvar
          write(io6, "(A,"//minformat//",A,"//maxformat//")") &
               '              Acceptable range is ', dminv, ' to ', dmaxv
       endif

    elseif (iflag .eq. 2) then

       if (dvar < dminv) then

          write(io6, "(A,"//varformat//")") &
               '-- WARNING -- Input variable '//trim(name)//' = ', dvar
          write(io6, "(A,"//minformat//")") &
               '              is less than the minimum acceptable value ', dminv
          write(io6, "(A,"//minformat//",A/)") &
               '   Will set '//trim(name)//' = ', dminv, ' for this run.'
          dvar = dminv
          nwarn = nwarn + 1

       elseif (dvar > dmaxv) then

          write(io6, "(A,"//varformat//")") &
               '-- WARNING -- Input variable '//trim(name)//' = ', dvar
          write(io6, "(A,"//maxformat//")") &
               '              is larger than the maximum acceptable value ', dmaxv
          write(io6, "(A,"//maxformat//",A/)") &
               '   Will set '//trim(name)//' = ', dmaxv, ' for this run.'
          dvar = dmaxv
          nwarn = nwarn + 1

       endif
    endif

    if (present(msgmin) .and. (dvar < dminv)) write(io6,*) msgmin
    if (present(msgmax) .and. (dvar > dmaxv)) write(io6,*) msgmax
    if (present(msgboth) .and. ((dvar < dminv) .or. (dvar > dmaxv))) &
         write(io6,*) msgboth
    if ((dvar < dminv) .or. (dvar > dmaxv)) write(io6,*) ''
    
  end subroutine dchk_bnds
    


  subroutine dsetformat(fstring, var)
    character(len=*), intent(out) :: fstring
    real(kind=8),     intent(in)  :: var
    if ((abs(var) >= 1.0d7) .or. (abs(var) <= 1.0d-4)) then
       fstring = "ES12.5"
    else
       fstring = "F0.7"
    endif
      
  end subroutine dsetformat


end subroutine oname_check





