subroutine read_nl(file)

  use oname_coms,   only: nl, cmdlne_runtype, cmdlne_fields, numcf
  use max_dims,     only: pathlen
  use mem_para,     only: olam_mpi_barrier, myrank

  implicit none

  character(*), intent(in) :: file
  character(pathlen)       :: fs
  integer                  :: il, ios
  logical                  :: fexists

  namelist /OLAMIN/ nl

! OPEN THE NAMELIST FILE

  inquire(file=file, exist=fexists)
  if (.not. fexists) then
     write(*,*) "The namelist file "//trim(file)//" is missing."
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

     ! unit io6 not setup yet
     if (myrank == 0) write(*,*) trim(fs)

     read(fs, nml=OLAMIN, iostat=ios)

     if (ios /= 0) then
        fexists = .true.

        ! unit io6 not setup yet
        if (myrank == 0) then
           write(*,*)
           write(*,*) "Error setting name list values from command line:"
           write(*,*) trim(cmdlne_fields(il))
           write(*,*) ios
           write(*,*) trim(fs)
           write(*,*) "Stopping model run."
        endif
        exit
     endif
  enddo

  call olam_mpi_barrier()
  if (fexists) stop

! OVERWRITE NAMELIST WITH COMMAND LINE RUNTYPE

  if (len_trim(cmdlne_runtype) > 1) then
     nl%runtype = cmdlne_runtype
  endif

  deallocate(cmdlne_fields)

end subroutine read_nl

!===============================================================================

subroutine copy_nl()

  use max_dims,    only: maxgrds, maxisdirs
  use consts_coms, only: r8
  use oname_coms,  only: nl
  use misc_coms,   only: expnme, runtype, timeunit, timmax8, &
                         nacoust, idiffk, csz, csx, akmin, &
                         dtlong, initial, zonclim, topo_database, bathym_database, &
                         gridfile, hfilin, ioutput, hfilepref, iclobber, &
                         frqstate, naddsc, ilwrtyp, iswrtyp, radfrq, &
                         icfrac, cfracrh1, cfracrh2, cfraccup, nqparm, confrq, &
                         nsndg, ipsflg, itsflg, irtsflg, iusflg, &
                         hs, p_sfc, us, vs, ts, ps, rts, &
                         itime1, idate1, imonth1, iyear1, ngrids, &
                         nzp, mdomain, itopoflg, ibathflg, nxp, &
                         ngrdll, grdrad, grdlat, grdlon, deltax, ndz, hdz, dz, &
                         current_time, initial_time, debug_fp, init_nans, do_chem, &
                         nrk_wrtv, nrk_scal, topodb_cutoff

  use micro_coms,  only: miclevel, icloud, idriz, irain, ipris, isnow, iaggr, &
                         igraup, ihail, iccn, igccn, iifn, &
                         rparm, sparm, aparm, gparm, hparm, &
                         ccnparm, gccnparm, ifnparm

  use mem_co2,     only: co2flag, co2_initppm

  use hcane_rz,    only: ncycle_hurrinit, timmax_hurrinit, hlat0, hlon0,     &
                         rad1_blend, rad2_blend, zcent_thpert, zhwid_thpert, &
                         rcent_thpert, rhwid_thpert, maxrate_thpert, &
                         vtan_targ, pmsl_targ

  use leaf_coms,   only: nvgcon, isoilflg, isoilptf, ndviflg, &
                         isfcl, ivegflg, &
                         veg_database, soil_database, soilgrids_database, &
                         glhymps_database, ndvi_database, iupdndvi, &
                         isoilstateinit, iwatertabflg, watertab_db

  use mem_land,    only: nzg, landgrid_dztop, landgrid_depth

  use sea_coms,    only: isstflg, sst_database, seatmp, iupdsst, &
                         iseaiceflg, seaice_database, seaice, iupdseaice, tide_database

  use oplot_coms,  only: op
  use isan_coms,   only: iapr
  use mem_nudge,   only: tnudcent, nudflag, nudnxp,  &
                         o3nudflag, tnudi_o3, o3nudpress
  use mem_rayf,    only: rayf_zmin,    rayf_fact,    rayf_expon,    &
                         rayfw_zmin,   rayfw_fact,   rayfw_expon,   &
                         rayfdiv_zmin, rayfdiv_fact, rayfdiv_expon, &
                         rayfmix_zmin, rayfmix_fact, rayfmix_expon

  use mem_sfcg,    only: nsfcgrids, sfcgrid_res_factor, nxp_sfc, &
                         nsfcgrdll, sfcgrdrad, sfcgrdlat, sfcgrdlon, sfcgfile, &
                         nswmzons, nswmzonll, swmzonrad, swmzonlat, swmzonlon

  use mem_sea,     only: npomzons, npomzonll, pomzonrad, pomzonlat, pomzonlon

  use sea_swm,     only: niter_swm

  use umwm_module, only: umwmflg

  use pom2k1d,     only: nzpom, pom_dztop, pom_depth

  use mem_sfcnud,  only: gw_spinup_sfcgfile, gw_spinup_histfile

  implicit none

  integer  :: i, ihrs, imns
  real(r8) :: tfact

! Copy namelist variables to their counterparts that are used within the model.
! If this model run is a 'HISTORY', 'HISTREGRID', or 'PLOTONLY' runtype, some
! of the copied namelist values will be overwritten by their values read from
! the history file in order to ensure that they are not changed from previous
! runs.

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

  itime1    = nl%itime1
  idate1    = nl%idate1
  imonth1   = nl%imonth1
  iyear1    = nl%iyear1

  ihrs = itime1 / 100
  imns = mod(itime1, 100)

  initial_time%year  = iyear1
  initial_time%month = imonth1
  initial_time%date  = idate1
  initial_time%time  = ihrs * 3600.0_r8 + imns * 60.0_r8

  ! If this is not a 'HISTORY', 'HISTREGRID', or 'PLOTONLY' runtype, set
  ! current time to initial time here.  Otherwise, current time will be
  ! set in subroutine history_start.

  if (runtype /= 'HISTORY' .and. runtype /= 'HISTREGRID' .and. &
      runtype /= 'PLOTONLY') then

     current_time%year  = iyear1
     current_time%month = imonth1
     current_time%date  = idate1
     current_time%time  = ihrs * 3600.0_r8 + imns * 60.0_r8
  endif

  nzp       = nl%nzp
  ndz       = nl%ndz

  hdz(1:ndz) = nl%hdz(1:ndz)
  dz (1:ndz) = nl%dz (1:ndz)

  mdomain   = nl%mdomain
  nxp       = nl%nxp
  deltax    = nl%deltax

  ngrids     = nl%ngrids

  ngrdll = nl%ngrdll
  grdrad = nl%grdrad
  grdlat = nl%grdlat
  grdlon = nl%grdlon

  sfcgrid_res_factor = nl%sfcgrid_res_factor

  nxp_sfc            = nl%sfcgrid_res_factor * nl%nxp

  nsfcgrids = nl%nsfcgrids
  nsfcgrdll = nl%nsfcgrdll

  sfcgrdrad = nl%sfcgrdrad
  sfcgrdlat = nl%sfcgrdlat
  sfcgrdlon = nl%sfcgrdlon

  nswmzons = nl%nswmzons

  if (nswmzons > 0) then
     nswmzonll = nl%nswmzonll

     swmzonrad = nl%swmzonrad
     swmzonlat = nl%swmzonlat
     swmzonlon = nl%swmzonlon
  endif

  npomzons = nl%npomzons

  if (npomzons > 0) then
     npomzonll = nl%npomzonll

     pomzonrad = nl%pomzonrad
     pomzonlat = nl%pomzonlat
     pomzonlon = nl%pomzonlon
  endif

  do i = 1,maxgrds
     idiffk(i)  = nl%idiffk(i)
     csz(i)     = nl%csz(i)
     csx(i)     = nl%csx(i)
     akmin(i)   = nl%akmin(i)
     nqparm(i)  = nl%nqparm(i)
  enddo

  dtlong        = nl%dtlong
  nacoust       = nl%nacoust
  initial       = nl%initial
  zonclim       = nl%zonclim
  nudflag       = nl%nudflag
  nudnxp        = nl%nudnxp
  tnudcent      = nl%tnudcent
  ioutput       = nl%ioutput
  hfilepref     = nl%hfilepref
  iclobber      = nl%iclobber
  frqstate      = nl%frqstate
  hfilin        = nl%hfilin

  naddsc        = nl%naddsc
  nrk_wrtv      = nl%acoust_timestep_level
  nrk_scal      = nl%scalar_timestep_level
  debug_fp      = nl%debug_fp
  init_nans     = nl%init_nans

  rayf_zmin     = nl%rayf_zmin
  rayf_fact     = nl%rayf_fact
  rayf_expon    = nl%rayf_expon
  rayfw_zmin    = nl%rayfw_zmin
  rayfw_fact    = nl%rayfw_fact
  rayfw_expon   = nl%rayfw_expon
  rayfdiv_zmin  = nl%rayfdiv_zmin
  rayfdiv_fact  = nl%rayfdiv_fact
  rayfdiv_expon = nl%rayfdiv_expon
  rayfmix_zmin  = nl%rayfmix_zmin
  rayfmix_fact  = nl%rayfmix_fact
  rayfmix_expon = nl%rayfmix_expon

  ilwrtyp       = nl%ilwrtyp
  iswrtyp       = nl%iswrtyp
  radfrq        = nl%radfrq
  icfrac        = nl%icfrac
  cfracrh1      = nl%cfracrh1
  cfracrh2      = nl%cfracrh2
  cfraccup      = nl%cfraccup
  confrq        = nl%confrq

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

  ncycle_hurrinit = nl%ncycle_hurrinit
  timmax_hurrinit = nl%timmax_hurrinit
  hlat0           = nl%hlat0
  hlon0           = nl%hlon0
  rad1_blend      = nl%rad1_blend
  rad2_blend      = nl%rad2_blend
  zcent_thpert    = nl%zcent_thpert
  zhwid_thpert    = nl%zhwid_thpert
  rcent_thpert    = nl%rcent_thpert
  rhwid_thpert    = nl%rhwid_thpert
  maxrate_thpert  = nl%maxrate_thpert
  vtan_targ       = nl%vtan_targ
  pmsl_targ       = nl%pmsl_targ * 1.e2  ! mb to Pa

  nsndg         = nl%nsndg
  ipsflg        = nl%ipsflg
  itsflg        = nl%itsflg
  irtsflg       = nl%irtsflg
  iusflg        = nl%iusflg

  hs(1) = nl%hs
  p_sfc = nl%p_sfc
  do i = 1,nl%nsndg
     ps(i)  = nl%sounding(1,i)
     ts(i)  = nl%sounding(2,i)
     rts(i) = nl%sounding(3,i)
     us(i)  = nl%sounding(4,i)
     vs(i)  = nl%sounding(5,i)
  enddo

  gridfile      = nl%gridfile
  sfcgfile      = nl%sfcgfile

  isfcl     = nl%isfcl
  gw_spinup_sfcgfile = nl%gw_spinup_sfcgfile
  gw_spinup_histfile = nl%gw_spinup_histfile
  nzg       = nl%nzg
  landgrid_dztop = nl%landgrid_dztop
  landgrid_depth = nl%landgrid_depth
  niter_swm = nl%niter_swm
  umwmflg   = nl%umwmflg
  nzpom     = nl%nzpom
  pom_dztop = nl%pom_dztop
  pom_depth = nl%pom_depth

  isoilflg      = nl%isoilflg
  isoilptf      = nl%isoilptf
  itopoflg  = nl%itopoflg
  ibathflg  = nl%ibathflg
  ivegflg   = nl%ivegflg
  ndviflg       = nl%ndviflg
  isstflg       = nl%isstflg
  iseaiceflg    = nl%iseaiceflg
  isoilstateinit= nl%isoilstateinit
  iwatertabflg  = nl%iwatertabflg
  topodb_cutoff      = nl%topodb_cutoff

  topo_database(:)   = nl%topo_database(:)
  bathym_database    = nl%bathym_database
  veg_database       = nl%veg_database
  soil_database      = nl%soil_database
  soilgrids_database = nl%soilgrids_database
  glhymps_database   = nl%glhymps_database
  ndvi_database      = nl%ndvi_database
  watertab_db        = nl%watertab_db
  sst_database       = nl%sst_database
  seaice_database    = nl%seaice_database
  tide_database      = nl%tide_database
  iupdndvi      = nl%iupdndvi
  iupdsst       = nl%iupdsst
  iupdseaice    = nl%iupdseaice
  seatmp        = nl%seatmp
  seaice        = nl%seaice
  nvgcon        = nl%nvgcon

  iapr(1:maxisdirs) = nl%iapr(1:maxisdirs)

  do_chem       = nl%do_chem
  o3nudflag     = nl%o3nudflag
  tnudi_o3      = 1.0 / max( nl%o3tnudcent, real(nl%dtlong) )
  o3nudpress    = nl%o3nudpress * 100.0  ! mb to Pa

! Variables from $MODEL_PLOT namelist

  op%nplt_files = nl%nplt_files

  do i = 1,op%nplt_files
     op%plt_files(i) = nl%plt_files(i)
  enddo

  op%nplt       = nl%nplt
  op%frqplt     = nl%frqplt
  op%dtvec      = nl%dtvec
  op%headspeed  = nl%headspeed
  op%stemlength = nl%stemlength
  op%vec_maxmrl = nl%vec_maxmrl
  op%prtval_size= nl%prtval_size
  op%plttype    = nl%plttype
  op%pltname    = nl%pltname
  op%pltorient  = nl%pltorient
  op%mapcolor   = nl%mapcolor
  op%llcolor    = nl%llcolor

  op%nx_grid    = nl%nx_grid
  op%nx_vect    = nl%nx_vect

  if (nl%ifont > 0) then
     op%ncarg_font = nl%ifont
  endif

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

     if (index(nl%plotspecs(i)%pltspec2,'I') > 0) op%pltindx1(i) = 'I'
     if (index(nl%plotspecs(i)%pltspec2,'J') > 0) op%pltindx1(i) = 'J'

     if (index(nl%plotspecs(i)%pltspec2,'B') > 0) op%vectbarb(i) = 'B'
     if (index(nl%plotspecs(i)%pltspec2,'V') > 0) op%vectbarb(i) = 'V'
     if (index(nl%plotspecs(i)%pltspec2,'w') > 0) op%vectbarb(i) = 'w'

     if (index(nl%plotspecs(i)%pltspec2,'Y') > 0) op%vectsea(i)  = 'Y'
     if (index(nl%plotspecs(i)%pltspec2,'y') > 0) op%vectsea(i)  = 'y'
     if (index(nl%plotspecs(i)%pltspec2,'z') > 0) op%vectsea(i)  = 'z'

     if (index(nl%plotspecs(i)%pltspec2,'G') > 0) op%pltgrid(i) = 'G'

     if (index(nl%plotspecs(i)%pltspec2,'g') > 0) op%pltgrid_sfc(i) = 'g'

     if (index(nl%plotspecs(i)%pltspec2,'D') > 0) op%pltdualgrid(i) = 'D'

     if (index(nl%plotspecs(i)%pltspec2,'i') > 0) op%pltindx2(i) = 'i'
     if (index(nl%plotspecs(i)%pltspec2,'j') > 0) op%pltindx2(i) = 'j'

     if (index(nl%plotspecs(i)%pltspec2,'b') > 0) op%pltborder(i) = 'b'
     if (index(nl%plotspecs(i)%pltspec2,'t') > 0) op%pltborder(i) = 't'
     if (index(nl%plotspecs(i)%pltspec2,'a') > 0) op%pltborder(i) = 'a'

     if (index(nl%plotspecs(i)%pltspec2,'n') > 0) op%labelbar(i) = 'n'
     if (index(nl%plotspecs(i)%pltspec2,'o') > 0) op%labelbar(i) = 'o'

     if (index(nl%plotspecs(i)%pltspec2,'c') > 0) op%colorbar(i) = 'c'
     if (index(nl%plotspecs(i)%pltspec2,'r') > 0) op%colorbar(i) = 'r'

     if (index(nl%plotspecs(i)%pltspec2,'m') > 0) op%maptyp(i) = 'm'

     if (index(nl%plotspecs(i)%pltspec2,'l') > 0) op%pltll(i) = 'l'

     if (index(nl%plotspecs(i)%pltspec2,'H') > 0) op%centplthurr(i) = 'H'

     if (index(nl%plotspecs(i)%pltspec2,'C') > 0) op%pltcone(i) = 'C'

     if (index(nl%plotspecs(i)%pltspec2,'p') > 0) op%pltlev(i) = 'p'
     if (index(nl%plotspecs(i)%pltspec2,'s') > 0) op%pltlev(i) = 's'
     if (index(nl%plotspecs(i)%pltspec2,'S') > 0) op%pltlev(i) = 'S'

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

     if (index(nl%plotspecs(i)%pltspec2,'R') > 0) op%gridded(i) = 'R'
  enddo

end subroutine copy_nl

!===============================================================================

subroutine namelist_print()

  use oname_coms, only: nl
  use misc_coms,  only: io6

  implicit none

  namelist /OLAMIN/ nl

  write(io6, nml=OLAMIN)

end subroutine namelist_print
