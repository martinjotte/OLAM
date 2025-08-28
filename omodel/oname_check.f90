subroutine oname_check()

! THIS ROUTINE CHECKS THE OPTION SPECIFICATIONS OF THE NAMELIST FILE OLAMIN
! FOR CONSISTENCY, AND OVERRIDES THE SETTINGS OF ICLOUD, IDRIZ, IRAIN, IPRIS,
! ISNOW, IAGGR, IGRAUP, and IHAIL, SETTING THEM ALL TO ZERO IF MICLEVEL < 3.

! EACH ERROR WILL BE ASSIGNED A SEVERITY.  FATAL ERRORS WILL CAUSE
! THE RUN TO STOP IMMEDIATELY AND WARNING ERRORS WILL BE LISTED.

use max_dims,    only: maxgrds, maxsndg, maxnplt, maxisdirs, &
                       maxpltfiles, maxngrdll
use oname_coms,  only: nl
use consts_coms, only: erad2, pi1, pi2, p00, r8
use misc_coms,   only: io6
use mem_para,    only: olam_stop

implicit none

integer  :: nfatal, nwarn, ng, i, k, iplt, i_huge, idz
real     :: r_huge, r_min, r_max, r_tiny, dzxmin, zb_min, zb_max, dtlong4
real(r8) :: d_huge, d_tiny

character(len=1), dimension(8) :: tunits = [ 'd','D','h','H','m','M','s','S' ]

nfatal = 0
nwarn  = 0

r_tiny = spacing(1.0)       ! real value small compared to 1.0
r_huge = huge(1.0) - 1.0    ! largest allowable real value
i_huge = huge(1) - 1        ! largest allowable integer
d_tiny = spacing(nl%timmax)     ! smallest real*8 value
d_huge = huge(nl%timmax) - 1.d0 ! largest allowable real*8 value
dzxmin = 0.001              ! minimum grid spacing (horiz. and vertical)

zb_min = -2.e3     ! min. vertical base (allow gridpoints below MSL)
if (nl%mdomain == 0) then
   zb_max = 0.0    ! global run - vertical base at MSL
else
   zb_max = 2.e3  ! max. vertical base (allow above MSL for a regional run)
endif

dtlong4 = real(nl%dtlong)

write(io6,*) ' '
write(io6,*) '  CHECKING NAMELIST VALUES'
write(io6,*) ' '

!--------------------------------------------------------------------------
! RUNTYPE
!--------------------------------------------------------------------------

if ( (nl%runtype /= 'MAKEGRID'     ) .and. &
     (nl%runtype /= 'INITIAL'      ) .and. &
     (nl%runtype /= 'HISTORY'      ) .and. &
     (nl%runtype /= 'MAKEREGRID'  ) .and. &
     (nl%runtype /= 'HISTREGRID'  ) .and. &
     (nl%runtype /= 'PLOTONLY'     ) .and. &
     (nl%runtype /= 'MAKEGRID_PLOT') ) then
   write(io6,*) " -- FATAL -- RUNTYPE = "//trim(nl%runtype)
   write(io6,*) "             RUNTYPE must be either 'MAKEGRID', 'INITIAL', "
   write(io6,*) "             'HISTORY', 'MAKEREGRID', 'HISTREGRID', "
   write(io6,*) "             'MAKEGRID_PLOT', or 'PLOTONLY'"
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

call ichk_bnds( nl%mdomain, "MDOMAIN",       0,       5, 0, nfatal, nwarn )

call ichk_bnds( nl%nzp,         "NZP",       2,   10000, 0, nfatal, nwarn, &
     msgmin="At least 3 vertical levels are needed for OLAM." )

call ichk_bnds( nl%nxp,         "NXP",       1,   10000, 0, nfatal, nwarn )
call dchk_bnds( nl%dtlong,   "DTLONG",  d_tiny,  d_huge, 0, nfatal, nwarn )
call rchk_bnds( nl%deltax,   "DELTAX",  dzxmin,  r_huge, 0, nfatal, nwarn )
call ichk_bnds( nl%ndz,         "NDZ",       1,      10, 0, nfatal, nwarn )

do idz=1, nl%ndz
   call rchk_bnds( nl%hdz(idz), "HDZ", -2000.0, r_huge, 0, nfatal, nwarn )
   call rchk_bnds( nl%dz (idz),  "DZ",     0.0, r_huge, 0, nfatal, nwarn )

   if (idz > 1) then
      if (nl%hdz(idz) <= nl%hdz(idz-1) ) then
         write(io6,*) 'FATAL - hdz(',idz,') must be larger than hdz(',idz-1,').'
         nfatal = nfatal + 1
      endif
   endif
enddo

if (nl%ndz == 1) then

   call rchk_bnds( nl%zz(1),     "ZZ",  zb_min,  zb_max, 0, nfatal, nwarn )
   do k=2,nl%nzp
      r_min = nl%zz(k-1) + dzxmin
      call rchk_bnds( nl%zz(k),  "ZZ",   r_min,  r_huge, 0, nfatal, nwarn )
   enddo

endif

!--------------------------------------------------------------------------
! NESTED GRID DEFINITION
!--------------------------------------------------------------------------

call ichk_bnds( nl%ngrids,   "NGRIDS",       1, maxgrds, 0, nfatal, nwarn, &
     msgmax="Increase maxgrds in max_dims.f90 if more nests are needed." )

call ichk_bnds( nl%gridplot_base, "GRIDPLOT_BASE", 2, maxgrds, 2, nfatal, nwarn )

! global mesh requires nxp to be divisible by 3

if (nl%mdomain < 2 .and. mod(nl%nxp,3) /= 0) then
   write(io6,*) 'FATAL - NXP must be divisible by 3 for a global run'
   nfatal = nfatal + 1
endif

do ng = 2, nl%ngrids
   call ichk_bnds(nl%ngrdll(ng), "NGRDLL", 1, maxngrdll, 0, nfatal, nwarn )
enddo

if (nl%mdomain < 2) then
   do ng = 2, nl%ngrids
      do i = 1, nl%ngrdll(ng)
         call rchk_bnds( nl%grdrad(ng,i), "GRDRAD", dzxmin, erad2, 0, nfatal, nwarn )
         call rchk_bnds( nl%grdlat(ng,i), "GRDLAT",   -90.,   90., 0, nfatal, nwarn )
         call rchk_bnds( nl%grdlon(ng,i), "GRDLON",  -180.,  180., 0, nfatal, nwarn )
      enddo
   enddo
endif

!--------------------------------------------------------------------------
! SURFACE NESTED GRID DEFINITION
!--------------------------------------------------------------------------

call ichk_bnds( nl%nsfcgrids,   "NSFCGRIDS", 0, maxgrds, 0, nfatal, nwarn, &
     msgmax="Increase maxgrds in max_dims.f90 if more sfc nests are needed." )

call ichk_bnds( nl%sfcgridplot_base, "SFCGRIDPLOT_BASE", 1, maxgrds, 2, nfatal, nwarn )

! If greater than 1, sfcgrid_res_factor must have prime factors of 2 and/or 3

ng = nl%sfcgrid_res_factor
if (ng > 1) then

   do while( mod(ng,2) == 0 )
      ng = ng / 2
   enddo

   do while( mod(ng,3) == 0 )
      ng = ng / 3
   enddo

   if (ng /= 1) then
      stop "Error: sfcgrid_res_factor must be 1 or have prime factors of only 2 and/or 3."
   endif
elseif (ng < 1) then
   stop "Error: sfcgrid_res_factor must be positive."
endif

do ng = 1, nl%nsfcgrids
   call ichk_bnds(nl%nsfcgrdll(ng), "NSFCGRDLL", 1, maxngrdll, 0, nfatal, nwarn )
enddo

if (nl%mdomain < 2) then
   do ng = 1, nl%nsfcgrids
      do i = 1, nl%nsfcgrdll(ng)
         call rchk_bnds( nl%sfcgrdrad(ng,i), "SFCGRDRAD", dzxmin, erad2, 0, nfatal, nwarn )
         call rchk_bnds( nl%sfcgrdlat(ng,i), "SFCGRDLAT",   -90.,   90., 0, nfatal, nwarn )
         call rchk_bnds( nl%sfcgrdlon(ng,i), "SFCGRDLON",  -180.,  180., 0, nfatal, nwarn )
      enddo
   enddo
endif

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
      call ichk_bnds( nl%max_nud_mrl, "MAX_NUD_MRL", 1, maxgrds, 2, nfatal, nwarn )
      call ichk_bnds( nl%nudnxp, "NUDNXP", 0, 10000, 2, nfatal, nwarn )
      call rchk_bnds( nl%tnudcent, "TNUDCENT", dtlong4, r_huge, 2, nfatal, &
                      nwarn, msgmin="Nudging time must be larger than dtlong")

      ! We don't need to do anything extra to preserve tracer mass and mixing ratio
      ! with NUD_PRESERVE_TOTAL_MASS active
      if (nl%nud_preserve_total_mass) nl%nud_preserve_mix_ratio = .false.
   endif

endif

!--------------------------------------------------------------------------
! HISTORY FILE OUTPUT
!--------------------------------------------------------------------------

call ichk_bnds( nl%ioutput       , "IOUTPUT"       , 0, 1, 2, nfatal, nwarn )
call ichk_bnds( nl%ioutput_mavg  , "IOUTPUT_MAVG"  , 0, 1, 2, nfatal, nwarn )
call ichk_bnds( nl%ioutput_davg  , "IOUTPUT_DAVG"  , 0, 1, 2, nfatal, nwarn )
call ichk_bnds( nl%ioutput_lite  , "IOUTPUT_LITE"  , 0, 1, 2, nfatal, nwarn )
call ichk_bnds( nl%ioutput_latlon, "IOUTPUT_LATLON", 0, 1, 2, nfatal, nwarn )

call ichk_bnds( nl%iclobber  , "ICLOBBER " , 0, 1, 2, nfatal, nwarn )
call ichk_bnds( nl%icompress , "ICOMPRESS" , 0, 9, 2, nfatal, nwarn )
call ichk_bnds( nl%latlonplot, "LATLONPLOT", 0, 1, 2, nfatal, nwarn )

call dchk_bnds( nl%frqstate , "FRQSTATE" , nl%dtlong, d_huge, 2, nfatal, nwarn )
call dchk_bnds( nl%frqlite  , "FRQLITE"  , nl%dtlong, d_huge, 2, nfatal, nwarn )
call dchk_bnds( nl%frqlatlon, "FRQLATLON", nl%dtlong, d_huge, 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! Topography
!--------------------------------------------------------------------------

call ichk_bnds( nl%itopoflg, "ITOPOFLG", 1,  2, 0, nfatal, nwarn )
call ichk_bnds( nl%ibathflg, "IBATHFLG", 1,  2, 0, nfatal, nwarn )

if (nl%nswmzons > 0) then
   if (nl%itopoflg == 2 .or. nl%ibathflg == 2) then
      write(*,*) "Shallow water tidal model requires both topography"
      write(*,*) "and bathymetry databases."
      stop "Error: topopgraphy and bathymetry must be present for nswmzons > 0"
   endif
endif

!--------------------------------------------------------------------------
! MODEL OPTIONS / NUMERICAL SCHEMES
!--------------------------------------------------------------------------

call ichk_bnds( nl%naddsc, "NADDSC", 0, i_huge, 2, nfatal, nwarn )

call ichk_bnds( nl%ithil_monot, "ITHIL_MONOT", 0, 1, 0, nfatal, nwarn )
call ichk_bnds( nl%iwind_monot, "IWIND_MONOT", 0, 1, 0, nfatal, nwarn )
call ichk_bnds( nl%iscal_monot, "ISCAL_MONOT", 0, 2, 0, nfatal, nwarn )

call ichk_bnds( nl%thil_horiz_adv_order, "THIL_HORIZ_ADV_ORDER", 1, 3, 2, nfatal, nwarn )
call ichk_bnds( nl%wind_horiz_adv_order, "WIND_HORIZ_ADV_ORDER", 1, 3, 2, nfatal, nwarn )
call ichk_bnds( nl%scal_horiz_adv_order, "SCAL_HORIZ_ADV_ORDER", 1, 3, 2, nfatal, nwarn )

call ichk_bnds( nl%acoust_timestep_level, "ACOUST_TIMESTEP_LEVEL", 2, 3, 2, nfatal, nwarn)
call ichk_bnds( nl%scalar_timestep_level, "SCALAR_TIMESTEP_LEVEL", 2, 3, 2, nfatal, nwarn)

call rchk_bnds( nl%divh_damp_fact,  "DIVH_DAMP_FACT",  0.0, 0.2, 2, nfatal, nwarn )
call ichk_bnds( nl%divh_damp_level, "DIVH_DAMP_LEVEL", 0,   3,   2, nfatal, nwarn )

call rchk_bnds( nl%vort_damp_fact,  "VORT_DAMP_FACT",  0.0, 0.1, 2, nfatal, nwarn )
call ichk_bnds( nl%vort_damp_level, "VORT_DAMP_LEVEL", 0,   3,   2, nfatal, nwarn )

!--------------------------------------------------------------------------
! RAYLEIGH FRICTION PARAMETERS
!--------------------------------------------------------------------------

call rchk_bnds( nl%rayf_fact,  "RAYF_FACT",  0.0,    1.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayf_expon, "RAYF_EXPON", 0.0,    5.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayf_zmin,  "RAYF_ZMIN" , 0.0, r_huge, 2, nfatal, nwarn )


call rchk_bnds( nl%rayfw_fact,  "RAYFW_FACT",  0.0, r_huge, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfw_expon, "RAYFW_EXPON", 0.0,    5.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfw_zmin,  "RAYFW_ZMIN" , 0.0, r_huge, 2, nfatal, nwarn )

call rchk_bnds( nl%rayfdiv_fact,  "RAYFDIV_FACT",  0.0,    0.3, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfdiv_expon, "RAYFDIV_EXPON", 0.0,    5.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfdiv_zmin,  "RAYFDIV_ZMIN" , 0.0, r_huge, 2, nfatal, nwarn )

call rchk_bnds( nl%rayfdiv_fact,  "RAYFMIX_FACT",  0.0,    2.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfdiv_expon, "RAYFMIX_EXPON", 0.0,    5.0, 2, nfatal, nwarn )
call rchk_bnds( nl%rayfdiv_zmin,  "RAYFMIX_ZMIN" , 0.0, r_huge, 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! RADIATION PARAMETERIZATION PARAMETERS
!--------------------------------------------------------------------------

! With only RRTMg, 0 turns off rad and positive turns on rad
call ichk_bnds( nl%iswrtyp, "ISWRTYP", 0, 3, 2, nfatal, nwarn )
call ichk_bnds( nl%ilwrtyp, "ILWRTYP", 0, 3, 2, nfatal, nwarn )

call dchk_bnds( nl%radfrq,    "RADFRQ", nl%dtlong, d_huge, 2, nfatal, nwarn )
call ichk_bnds( nl%icfrac,    "ICFRAC", 0,   6,  0, nfatal, nwarn )
call rchk_bnds( nl%cfracrh1,"CFRACRH1", 0.,  2., 0, nfatal, nwarn )
call rchk_bnds( nl%cfracrh2,"CFRACRH2", 0.,  2., 0, nfatal, nwarn )
call rchk_bnds( nl%cfraccup,"CFRACCUP", 0.1,  1., 0, nfatal, nwarn )

! Currently, radiation type 1 is invalid
!if ( nl%iswrtyp==1 .or. nl%ilwrtyp==1 ) then
!   write(io6,*) 'FATAL - Radiation type 1 is invalid'
!   nfatal = nfatal + 1
!endif

if (nl%icfrac == 2 .and. nl%cfracrh2 - nl%cfracrh1 < 1.e-6 ) then
   write(io6,*) 'FATAL - CFRACRH2 must be greater than CFRACRH1'
   nfatal = nfatal + 1
endif
!--------------------------------------------------------------------------
! CUMULUS PARAMETERIZATION PARAMETERS
!--------------------------------------------------------------------------

do ng=1, nl%ngrids
   call ichk_bnds( nl%nqparm(ng), "NQPARM", 0, 3, 0, nfatal, nwarn )
enddo

call dchk_bnds( nl%confrq, "CONFRQ", nl%dtlong, d_huge, 2, nfatal, nwarn )

call ichk_bnds( nl%conv_uv_mix,     "CONV_UV_MIX",     0, 1, 0, nfatal, nwarn )
call ichk_bnds( nl%conv_tracer_mix, "CONV_TRACER_MIX", 0, 1, 0, nfatal, nwarn )

!--------------------------------------------------------------------------
! EDDY DIFFUSION PARAMETERS
!--------------------------------------------------------------------------

do ng=1, nl%ngrids
   call ichk_bnds( nl%idiffk(ng), "IDIFFK", 0,     3, 0, nfatal, nwarn )
   call rchk_bnds( nl%csx(ng),    "CSX",    0.,  10., 0, nfatal, nwarn )
   call rchk_bnds( nl%csz(ng),    "CSZ",    0.,  10., 0, nfatal, nwarn )
   call rchk_bnds( nl%akmin(ng),  "AKMIN",  0.,  10., 0, nfatal, nwarn )
enddo

!--------------------------------------------------------------------------
! MICROPHYSICS PARAMETERS
!--------------------------------------------------------------------------

call ichk_bnds( nl%miclevel, "MICLEVEL", 0, 3, 0, nfatal, nwarn )

if (nl%miclevel <= 2) then
   nl%icloud  = 0
   nl%idriz   = 0
   nl%irain   = 0
   nl%ipris   = 0
   nl%isnow   = 0
   nl%iaggr   = 0
   nl%igraup  = 0
   nl%ihail   = 0

elseif (nl%miclevel == 3) then

   if (.not. any(nl%icloud == (/0, 4, 5/) )) then
      write(io6,*) 'FATAL - icloud must be set to 0, 4, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%idriz == (/0, 5/) )) then
      write(io6,*) 'FATAL - idriz must be set to 0 or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%irain == (/0, 2, 5/) )) then
      write(io6,*) 'FATAL - irain must be set to 0, 2, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%ipris == (/0, 5/) )) then
      write(io6,*) 'FATAL - ipris must be set to 0 or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%isnow == (/0, 2, 5/) )) then
      write(io6,*) 'FATAL - isnow must be set to 0, 2, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%iaggr == (/0, 2, 5/) )) then
      write(io6,*) 'FATAL - iaggr must be set to 0, 2, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%igraup == (/0, 2, 5/) )) then
      write(io6,*) 'FATAL - igraup must be set to 0, 2, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%ihail == (/0, 2, 5/) )) then
      write(io6,*) 'FATAL - ihail must be set to 0, 2, or 5.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%iccn == (/1, 2, 3, 4, 5, 6/) )) then
      write(io6,*) 'FATAL - iccn must be set to 1, 2, 3, 4, 5, or 6.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%igccn == (/1, 2/) )) then
      write(io6,*) 'FATAL - igccn must be set to 1 or 2.'
      nfatal = nfatal + 1
   endif

   if (.not. any(nl%iifn == (/1, 2/) )) then
      write(io6,*) 'FATAL - iifn must be set to 1 or 2.'
      nfatal = nfatal + 1
   endif

   if (nl%icloud == 5 .and. nl%iccn < 2) then
      write(io6,*) 'FATAL - if icloud = 5, iccn must be at least 2.'
      nfatal = nfatal + 1
   endif

   if (nl%icloud > 0 .and. nl%irain > 0 .and. nl%idriz /= 5) then
      write(io6,*) 'FATAL - idriz must set to 5 if both cloud and rain are turned on.'
      nfatal = nfatal + 1
   endif

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

   call rchk_bnds( nl%ccnparm,  "CCNPARM",  0., 2000.e6, 0, nfatal, nwarn )
   call rchk_bnds( nl%gccnparm, "GCCNPARM", 0., 100.e3, 0, nfatal, nwarn )
   call rchk_bnds( nl%ifnparm,  "IFNPARM",  0., 100.e3, 0, nfatal, nwarn )

   if (nl%icloud == 4) then

      if (nl%iccn /= 1) then
         write(io6,*) 'FATAL - ICCN must be set to 1 when ICLOUD = 4'
         nfatal = nfatal + 1
      endif

      if (nl%igccn /= 1) then
         write(io6,*) 'FATAL - IGCCN must be set to 1 when ICLOUD = 4'
         nfatal = nfatal + 1
      endif

   endif

endif

!--------------------------------------------------------------------------
! CO2 PARAMETERS
!--------------------------------------------------------------------------

 call ichk_bnds( nl%co2flag       , "CO2FLAG"       , 0, 1, 2, nfatal, nwarn )

!--------------------------------------------------------------------------
! HURRICANE DYNAMIC INITIALIZATION PARAMETERS
!--------------------------------------------------------------------------

   call ichk_bnds( nl%ncycle_hurrinit,  "NCYCLE_HURRINIT",  0, 20, 0, nfatal, nwarn )

if (  nl%ncycle_hurrinit > 0 .and. &
     (nl%test_case /= 0      .or.  &
      nl%initial   /= 2)   ) then
   write(io6,*) " -- FATAL -- Hurricane dynamic initialization requires that "
   write(io6,*) "             TEST_CASE must be 0 and INITIAL must be 2 "
   nfatal = nfatal + 1
endif

if (  nl%ncycle_hurrinit > 1 .and. &
      nl%runtype == 'HISTREGRID') then
   write(io6,*) " -- FATAL -- Hurricane relocation (for which ncycle_hurrinit > 1) "
   write(io6,*) "             cannot be done in a HISTREGRID run "
   nfatal = nfatal + 1
endif

!--------------------------------------------------------------------------
! SOUNDING SPECIFICATION
!--------------------------------------------------------------------------

if (nl%initial == 1) then

   call ichk_bnds( nl%nsndg, "NSNDG", 2, maxsndg, 0, nfatal, nwarn, &
        msgmin="Sounding must have at least 2 levels.", &
        msgmax="Increase maxsndg in max_dims.f90 if more levels are needed." )

   call ichk_bnds( nl%ipsflg,  "IPSFLG",  0, 1, 0, nfatal, nwarn )
   call ichk_bnds( nl%itsflg,  "ITSFLG",  0, 3, 0, nfatal, nwarn )
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

   call ichk_bnds( nl%nzg, "NZG", 2, 100, 0, nfatal, nwarn, &
        msgmin="At least 2 soil levels are needed for soil model.")
   call rchk_bnds( nl%landgrid_dztop, "LANDGRID_DZTOP", 0.02, 1.1, 0, nfatal, nwarn )
   call rchk_bnds( nl%landgrid_depth, "LANDGRID_DEPTH", 1.0, 1000.0, 0, nfatal, nwarn )

   call ichk_bnds( nl%nzs, "NZS", 0, 10, 0, nfatal, nwarn )

   call ichk_bnds( nl%isoilflg,    "ISOILFLG",    1, 3, 0, nfatal, nwarn )
   call ichk_bnds( nl%isoilptf,    "ISOILPTF",    1, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%ivegflg,     "IVEGFLG",     1, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%ndviflg,     "NDVIFLG",     1, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%isstflg,     "ISSTFLG",     0, 2, 0, nfatal, nwarn )
   call ichk_bnds( nl%iseaiceflg,  "ISEAICEFLG",  0, 2, 0, nfatal, nwarn )

   call ichk_bnds( nl%isoilstateinit, "ISOILSTATEINIT", 0, 2, 0, nfatal,nwarn )

   call ichk_bnds( nl%iupdndvi,  "IUPDNDVI",    0, 1, 0, nfatal, nwarn )
   call ichk_bnds( nl%iupdsst,   "IUPDSST",     0, 1, 0, nfatal, nwarn )
   call ichk_bnds( nl%iupdseaice,"IUPDSEAICE",  0, 1, 0, nfatal, nwarn )

   ! Accept either + or - values, but the code assumes negative values for depth
   nl%zbedrock = -abs(nl%zbedrock)

   call ichk_bnds( nl%ihoriz_gndwater_transport, "IHORIZ_GNDWATER_TRANSPORT", 0, 1, 0, nfatal, nwarn )

   call ichk_bnds( nl%iseasprayflg , "ISEASPRAYFLG", 0, 2, 0, nfatal, nwarn )
   call rchk_bnds( nl%seaspray_vmin, "SEASPRAY_VMIN", 15.1, r_huge, 2, nfatal, nwarn )

   if (nl%isoilflg == 3) &
        call ichk_bnds( nl%isoiltext, "ISOILTEXT", 1, 12, 0, nfatal, nwarn )

   if (nl%ivegflg == 2) &
        call ichk_bnds( nl%nvgcon, "NVGCON", 0, 20, 0, nfatal, nwarn )

   if (nl%isstflg == 0) &
        call rchk_bnds( nl%seatmp, "SEATMP", 100., 500., 0, nfatal, nwarn )

   if (nl%iseaiceflg == 0) &
        call rchk_bnds( nl%seaice, "SEAICE", 0., 1., 2, nfatal, nwarn )

   ! Set iupdndvi to 0 if not reading NDVI files
   if (nl%ndviflg == 2) then
      if (nl%iupdndvi == 1) then
         write(io6,*) "Turning off NDVI update for ndviflg == 2"
         nl%iupdndvi = 0
         nwarn = nwarn + 1
      endif
   endif

   ! Set iupdsst to 0 if not reading SST files
   if (nl%isstflg == 0) then
      if (nl%iupdsst == 1) then
         write(io6,*) "Turning off SST update for isstflg == 2"
         nl%iupdsst = 0
         nwarn = nwarn + 1
      endif
   endif

   ! Set iupdseaice to 0 if not reading SEAICE files
   if (nl%iseaiceflg == 0) then
      if (nl%iupdseaice == 1) then
         write(io6,*) "Turning off SEAICE update for iseaiceflg == 2"
         nl%iupdseaice = 0
         nwarn = nwarn + 1
      endif
   endif

   ! Set isoilstateinit to 0 if not doing 3D initialization, as it requires
   ! the 3D analysis files
   if (nl%isoilstateinit == 1) then
      if (nl%initial /= 2) then
         write(io6,*) "Turning off SOIL initialization when not not doing 3d initialization"
         nl%isoilstateinit = 0
         nwarn = nwarn + 1
      endif
   endif

endif

!--------------------------------------------------------------------------
! CMAQ CHEMISTRY PARAMETERS
!--------------------------------------------------------------------------

call ichk_bnds( nl%do_chem , "DO_CHEM ", 0,   1, 2, nfatal, nwarn )
call ichk_bnds( nl%ltng_nox, "LTNG_NOx", 0,   1, 2, nfatal, nwarn )
call dchk_bnds( nl%photfrq,  "PHOTFRQ", nl%dtlong, d_huge, 2, nfatal, nwarn )

if (nl%initial /= 2) then
   if (nl%o3nudflag /= 0) then
      write(io6,*) "Turning off ozone nudging when not not doing 3d initialization"
      nl%o3nudflag = 0
      nwarn = nwarn + 1
   endif
endif

if (nl%do_chem == 0) then
   if (nl%o3nudflag /= 0) then
      write(io6,*) "Turning off ozone nudging when not not chemistry"
      nl%o3nudflag = 0
      nwarn = nwarn + 1
   endif
endif

if (nl%initial == 2 .and. nl%do_chem == 1) then

   call ichk_bnds( nl%o3nudflag, "O3NUDFLAG", 0,  1, 2, nfatal, nwarn )

   if (nl%o3nudflag == 1) then

      call rchk_bnds( nl%o3tnudcent, "O3TNUDCENT", dtlong4, r_huge, 2, nfatal, &
                      nwarn, msgmin="Ozone nudging time must be larger than dtlong")

      call rchk_bnds( nl%o3nudpress, "O3NUDPRESS", 0.0, 1200.0, 2, nfatal, nwarn)

   endif

endif

!--------------------------------------------------------------------------
! ISENTROPIC CONTROL
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! PLOTTING PARAMETERS
!--------------------------------------------------------------------------

call ichk_bnds( nl%nplt, "NPLT", 0, maxnplt, 0, nfatal, nwarn )

if (nl%runtype == 'PLOTONLY') &
     call ichk_bnds( nl%nplt_files, "NPLT_FILES", 0, maxpltfiles, 0, nfatal, &
     nwarn, msgmax="Increase maxpltfiles in maxdims.f90 if you need to plot &
     & from more files" )

if ((nl%runtype == 'INITIAL') .or. &
    (nl%runtype == 'HISTORY') .or. &
    (nl%runtype == 'HISTREGRID')) &
     call dchk_bnds( nl%frqplt, "FRQPLT", nl%dtlong, d_huge, 2, nfatal, nwarn )

call rchk_bnds( nl%dtvec,     "DTVEC",     1.e-3, r_huge, 2, nfatal, nwarn )
call rchk_bnds( nl%headspeed, "HEADSPEED", 1.e-3, r_huge, 2, nfatal, nwarn )
call rchk_bnds( nl%stemlength,"STEMLENGTH",1.e-3, r_huge, 2, nfatal, nwarn )
call ichk_bnds( nl%plttype,   "PLTTYPE",   0, 2, 2, nfatal, nwarn )
call ichk_bnds( nl%pltorient, "PLTORIENT", 0, 1, 2, nfatal, nwarn )
call ichk_bnds( nl%mapcolor,  "MAPCOLOR",  0, 255, 2, nfatal, nwarn)
call ichk_bnds( nl%llcolor,   "LLCOLOR",   0, 255, 2, nfatal, nwarn)

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
        nl%plotspecs(iplt)%projectn == 'G'  .or.   &
        nl%plotspecs(iplt)%projectn == 'O') .and.  &
        nl%mdomain > 1) then

        write(io6,*) 'When mdomain > 1, may not use L, P, G, or O plot projection'
        nfatal = nfatal + 1
   endif

   if (nl%plotspecs(iplt)%projectn == 'C' .and.  &
       index(nl%plotspecs(iplt)%pltspec2,'T') > 0 .and. &
       nl%mdomain > 1) then

        write(io6,*) 'When mdomain > 1, may not use C plot projection for tileplot'
        nfatal = nfatal + 1

! Planning to relax this condition eventually.
! At present (12/20/2012), also may not use C plot projection for contouring at V
! point, but this is not checked for.

   endif
enddo

! TODO - add checks for the plotting variables

!------------------------------------------------------------------------
! OTHER CONSISTANCY CHECKS
!------------------------------------------------------------------------

! TEMPORARY!! ONLY ALLOW MDOMAIN VALUES OF 0, 3, 4, and 5

if (nl%mdomain /= 0 .and. nl%mdomain /= 3 .and. &
    nl%mdomain /= 4 .and. nl%mdomain /= 5) then
   write(io6,*) ' FATAL - MDOMAIN temporarily restricted to values of 0, 3, 4, or 5.'
   nfatal = nfatal + 1
endif

! IF THIS IS A GLOBAL SIMULATION, PRINT MESSAGE THAT DELTAX WILL BE
! REDEFINED BY NXP.

if (nl%mdomain == 0 .and. nfatal == 0) then

   write(io6,*) ' '
   write(io6,*) 'Since MDOMAIN is set to 0, this is configured as a global&
          & simulation.'
   write(io6,*) 'Hence, DELTAX is ignored and horizontal grid spacing will be&
          & computed'
   write(io6,*) 'from nxp and the Earth (or other planetary) radius.'
   write(io6,*) ' '

endif

! RUN WITH 3D-VARIABLE OR LONGITUDINALLY HOMOGENEOUS INITIALIZATION MUST BE GLOBAL

if ((nl%initial == 2 .or. nl%initial == 3) .and. nl%mdomain /= 0) then
   write(io6,*) ' FATAL  - mdomain must be 0 if INITIAL = 2 or 3.'
   nfatal = nfatal + 1
endif

! CONVECTIVE PARAMETERIZATION MUST HAVE AT LEAST WATER VAPOR

if (any(nl%nqparm(1:nl%ngrids) > 0) .and. nl%miclevel == 0) then
   write(io6,*) ' FATAL - MICLEVEL must be at least 1 for the cumulus parameterization'
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

if (nl%initial == 2 .and. nl%miclevel == 0) then
   write(io6,*) 'FATAL - Moisture complexity MICLEVEL must be 1 or larger if'
   write(io6,*) 'variable initialization is used (i.e., if INITIAL = 2).'
   nfatal = nfatal + 1
endif

! STOP THE RUN IF THERE ARE ANY FATAL ERRORS, AND LIST HOW MANY
! FATAL AND WARNING MESSAGES

write(io6,*) ''
write(io6,*) ' ----------oname_check--------------------------'
write(io6,*) ' FATAL     errors - ', nfatal
write(io6,*) ' WARNING   errors - ', nwarn
write(io6,*) ' -----------------------------------------------'
write(io6,*) ''
if (nfatal > 0) call olam_stop( 'ONAME_CHECK' )

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
    !  1 - REPORT A WARNING IF BOUNDS EXCEEDED, RUN CONTINUES UNCHANGED
    !  2 - IF BOUNDS EXCEEDED, RESET THE VARIABLE TO WITHIN THE GIVEN BOUNDS
    !      AND REPORT A WARNING, RUN CONTINUES
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
