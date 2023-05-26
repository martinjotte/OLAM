Module oname_coms

  use max_dims,    only: maxsndg, maxgrds, maxisdirs, maxnplt,  &
                         maxpltfiles, maxngrdll, pathlen, maxlite, maxlatlon
  use consts_coms, only: r8

  implicit none

  ! Do not re-export symbols from other modules

  private :: maxsndg, maxgrds, maxisdirs, maxnplt, &
             maxpltfiles, maxngrdll, pathlen, maxlite, maxlatlon, r8

  ! Derived types to hold the components of the plot specification fields

  Type oname_colortab
     integer       :: icolortab   = 0
     character(30) :: palette     = ' '
     character(3)  :: cscale      = ' '
     real          :: cmin        = 0.0
     real          :: cmax        = 0.0
     real          :: cinc        = 0.0
     real          :: zerohalfwid = 0.0
  End Type oname_colortab

  Type oname_plot
     character(30) :: fldname    = ' '
     character(1)  :: projectn   = ' '
     integer       :: icolortab  = 0
     character(20) :: pltspec2   = 'N'
     real          :: plotcoord1 = 0.0
     real          :: plotcoord2 = 0.0
     real          :: slabloc    = 0.0
     real          :: plotwid    = 0.0
     real          :: viewazim   = 0.0
  End Type oname_plot

  ! Derived type to hold the namelist variables

  Type oname_vars

!!    RUNTYPE/NAME

     character(64) :: expnme  = 'OLAM test'
     character(16) :: runtype = ' '

!!    SIMULATION ENDING TIME

     character(1) :: timeunit = ' '
     real(r8)     :: timmax   = 0.0_r8

!!    START OF SIMULATION OR ISAN PROCESSING

     integer :: itime1  = 0
     integer :: idate1  = 0
     integer :: imonth1 = 0
     integer :: iyear1  = 0

!!    GRID SPECIFICATIONS

     integer :: mdomain  = 0
     integer :: nzp      = 0
     integer :: nxp      = 0

     real(r8) :: dtlong = 0.0_r8
     integer :: ndz = 0

     real :: hdz(10) = 0.
     real :: dz (10) = 0.

     real :: zz(maxsndg) = 0.0

!!    LIMITED-AREA RUN SPECIFICATION

     real :: deltax = 1000.0

     real :: u_geostrophic = 0.0
     real :: v_geostrophic = 0.0

     real :: rlat0 = 45.0
     real :: rlon0 =  0.0

     real(r8) :: les_diag_freq = 0.0_r8

!!    NESTED GRID DEFINITION

     integer :: ngrids     = 1
     integer :: gridplot_base = 2

     integer :: ngrdll(maxgrds) = 0

     real    :: grdrad(maxgrds,maxngrdll) = 0.0
     real    :: grdlat(maxgrds,maxngrdll) = 0.0
     real    :: grdlon(maxgrds,maxngrdll) = 0.0

     integer :: nsfcgrids          =  0
     integer :: sfcgrid_res_factor =  1
     integer :: sfcgridplot_base   =  1

     integer :: nsfcgrdll(maxgrds) = 0

     real    :: sfcgrdrad(maxgrds,maxngrdll) = 0.0
     real    :: sfcgrdlat(maxgrds,maxngrdll) = 0.0
     real    :: sfcgrdlon(maxgrds,maxngrdll) = 0.0

     integer :: nswmzons = 0

     integer :: nswmzonll(maxgrds) = 0

     real    :: swmzonrad(maxgrds,maxngrdll) = 0.0
     real    :: swmzonlat(maxgrds,maxngrdll) = 0.0
     real    :: swmzonlon(maxgrds,maxngrdll) = 0.0

     integer :: npomzons = 0

     integer :: npomzonll(maxgrds) = 0

     real    :: pomzonrad(maxgrds,maxngrdll) = 0.0
     real    :: pomzonlat(maxgrds,maxngrdll) = 0.0
     real    :: pomzonlon(maxgrds,maxngrdll) = 0.0

!!    TIMESTEP RATIOS

     integer :: nacoust = 1

!!    VARIABLE INITIALIZATION INPUT

     integer            :: initial  = 0
     character(pathlen) :: zonclim  = '../etc/ZONAVG_CLIMATE'

!!    NUDGING PARAMETERS

     integer :: nudflag     = 0
     integer :: max_nud_mrl = maxgrds
     integer :: nudnxp      = 0
     real    :: tnudcent    = 86400.0

     logical :: nud_preserve_mix_ratio = .true.

!!    HISTORY/OUTPUT FILES

     logical :: save_node_logs = .true.

     integer :: ioutput      = 1
     integer :: ioutput_mavg = 1
     integer :: ioutput_davg = 1
     integer :: ioutput_lite = 0
     integer :: ioutput_latlon = 0

     integer :: iblocksize = -1
     integer :: iclobber  = 0
     integer :: icompress = 0
     integer :: latlonplot = 0
     integer :: nlatlonpd = 1

     real(r8):: frqstate  = 3600.0_r8
     real(r8):: frqlite   = 3600.0_r8
     real(r8):: frqlatlon = 3600.0_r8

     real :: beglat =  -90.
     real :: endlat =   90.
     real :: beglon = -180.
     real :: endlon =  180.

     character(pathlen) :: hfilin    = ' '
     character(pathlen) :: hfilepref = 'hist/'
     character(pathlen) :: lfilepref = 'hist/l'
     character(32)      :: lite_vars(maxlite) = ' '
     character(32)      :: latlon_vars(maxlatlon) = ' '

!!    MODEL/NUMERICAL OPTIONS

     integer :: naddsc      = 0

     integer :: ithil_monot = 0
     integer :: iscal_monot = 0

     integer :: horiz_adv_order = 2

     integer :: acoust_timestep_level = 3
     integer :: scalar_timestep_level = 2

     real    :: akmin_vort     = 0.4
     real    :: divh_damp_fact = 0.1

     logical :: zero_neg_scalars = .true.

     logical :: debug_fp    = .false.
     logical :: init_nans   = .false.

!!    RAYLEIGH FRICTION PARAMETERS

     real :: rayf_zmin  = 30000.0
     real :: rayf_fact  = 0.0
     real :: rayf_expon = 1.0

     real :: rayfw_zmin  = 30000.0
     real :: rayfw_fact  = 0.0
     real :: rayfw_expon = 1.0

     real :: rayfdiv_zmin  = 30000.0
     real :: rayfdiv_fact  = 0.0
     real :: rayfdiv_expon = 1.0

     real :: rayfmix_zmin  = 30000.0
     real :: rayfmix_fact  = 0.0
     real :: rayfmix_expon = 1.0

!!    RADIATION PARAMETERIZATION PARAMETERS

     integer :: iswrtyp = 0
     integer :: ilwrtyp = 0
     real(r8):: radfrq  = 1800.0_r8

     character(pathlen) :: rrtmg_datadir = '../etc'

     integer :: icfrac   = 0
     real    :: cfracrh1 = 0.8
     real    :: cfracrh2 = 1.1
     real    :: cfraccup = 0.7

!!    CUMULUS PARAMETERIZATION PARAMETERS

     integer :: nqparm(maxgrds) = 0
     integer :: conv_uv_mix     = 0
     integer :: conv_tracer_mix = 0

     real(r8) :: confrq          = 1800.0_r8

!!    EDDY DIFFUSION PARAMETERS

     integer :: idiffk(maxgrds) = 2
     real    :: csx   (maxgrds) = 0.2
     real    :: csz   (maxgrds) = 0.2
     real    :: akmin (maxgrds) = 0.0
     integer :: moist_buoy      = 1

!!    MICROPHYSICS PARAMETERS

     integer :: miclevel  = 3
     integer :: icloud = 4
     integer :: idriz  = 5
     integer :: ipris  = 5
     integer :: irain  = 2
     integer :: isnow  = 2
     integer :: iaggr  = 2
     integer :: igraup = 2
     integer :: ihail  = 2

     integer :: iccn   = 1
     integer :: igccn  = 1
     integer :: iifn   = 1

     real :: rparm  = .001   ! [m]
     real :: sparm  = .003   ! [m]
     real :: aparm  = .003   ! [m]
     real :: gparm  = .003   ! [m]
     real :: hparm  = .01    ! [m]

     real :: ccnparm  = 300.e6 ! [#/kg]
     real :: gccnparm =   3.e0 ! [#/kg]
     real :: ifnparm  =   0.e0 ! [#/kg]

!!    CO2 PARAMETERS

     integer :: co2flag       = 0
     real    :: co2_ppmv_init = 400.0

!!    HURRICANE DYNAMIC INITIALIZATION PARAMETERS

     integer  :: ncycle_hurrinit
     real(r8) :: timmax_hurrinit

     real :: hlat0          ! Obs hurricane latitude (deg)
     real :: hlon0          ! Obs hurricane longitude (deg)

     real :: rad1_blend     ! Inner radius for relocation blending weights (m)
     real :: rad2_blend     ! Outer radius for relocation blending weights (m)

     real :: zcent_thpert   ! Center height of toroidal heating region (m)
     real :: zhwid_thpert   ! Vertical half-width of toroidal heating region (m)

     real :: rcent_thpert   ! Center radius of toroidal heating region (m)
     real :: rhwid_thpert   ! Radial half-width of toroidal heating region (m)

     real :: maxrate_thpert ! Maximum heating rate (K/s) in toroidal heating region (K/s)

     real :: vtan_targ      ! Target maximum tangential wind speed (to modulate heating rate)
     real :: pmsl_targ      ! Target minimum sea level pressure (to modulate heating rate)

!!    SOUNDING SPECIFICATION

     integer :: nsndg   = 0
     integer :: ipsflg  = 1
     integer :: itsflg  = 0
     integer :: irtsflg = 3
     integer :: iusflg  = 0

     real :: hs    = 0.0
     real :: p_sfc = 1000.0

     real :: sounding(5,maxsndg) = 0.0

!!    ATM AND SFC GRID FILES PATH/NAMES

     character(pathlen) :: gridfile   = 'sfcfiles/gridfile_0'
     character(pathlen) :: sfcgfile   = 'sfcfiles/sfcgfile_0'

!!    TOPOGRAPHY, LEAF, SEA VARIABLES

     integer :: isfcl       =  1
     integer :: igw_spinup  =  0
     integer :: nzs         =  1
     integer :: nzg         = 21
     integer :: nzpom       = 40
     integer :: niter_swm   =  1

     real :: landgrid_dztop = 0.05
     real :: landgrid_depth = 5.00

     real :: pom_dztop = 5.0
     real :: pom_depth = 4500.

     integer :: isoilflg   = 2
     integer :: isoilptf   = 2
     integer :: isoiltext  = 4
     integer :: itopoflg   = 2
     integer :: ibathflg   = 2
     integer :: ivegflg    = 2
     integer :: ndviflg    = 2
     integer :: isstflg    = 0
     integer :: iseaiceflg = 0

     real :: zbedrock  = -10.0
     real :: gnd_ksat  = 2.e-7
     real :: gnd_poros = 0.12

     integer :: isoilstateinit = 0
     integer :: iwatertabflg   = 0
     integer :: iorogslopeflg  = 0

     integer :: iupdndvi   = 0
     integer :: iupdsst    = 0
     integer :: iupdseaice = 0
     integer :: nvgcon     = 8

     integer :: ihoriz_gndwater_transport = 0

     real    :: seatmp = 280.0
     real    :: seaice =   0.0

     real    :: topodb_cutoff = 200.

     character(pathlen) :: gw_spinup_sfcgfile = ' '
     character(pathlen) :: gw_spinup_histfile = ' '
     character(pathlen) :: topo_database(2)   = ' '
     character(pathlen) :: bathym_database    = ' '
     character(pathlen) :: veg_database       = ' '
     character(pathlen) :: soil_database      = ' '
     character(pathlen) :: soilgrids_database = ' '
     character(pathlen) :: glhymps_database   = ' '
     character(pathlen) :: ndvi_database      = ' '
     character(pathlen) :: sst_database       = ' '
     character(pathlen) :: seaice_database    = ' '
     character(pathlen) :: watertab_db        = ' '
     character(pathlen) :: orog_slope_db      = ' '

!!    ISENTROPIC CONTROL

     character(pathlen) :: iapr(maxisdirs) = ' '

!!    CMAQ Chemistry

     integer :: do_chem               =   0
     integer :: ltng_nox              =   0
     real(r8):: photfrq               = 600.0_r8
     integer :: aerosol_optics_photol = 0
     integer :: do_emis               = 1
     integer :: do_drydep             = 1
     integer :: do_aesedi             = 1

     character(pathlen) :: emis_dir  = '../../olamdatah5/edgar42'
     character(pathlen) :: geia_emis_file = '../etc/geia_emis.nc4'

     logical            :: megan_use_pfts_file = .false.
     character(pathlen) :: megan_pfts_file = '../../olamdatah5/megan/pfts.nc4'

     integer :: o3nudflag  = 0
     real    :: o3nudpress = 150.0
     real    :: o3tnudcent = 86400.0

!!    MODEL_PLOT VARIABLES

     integer :: nplt        = 0
     integer :: nplt_files  = 0
     integer :: nmavg_files = 0
     integer :: ndavg_files = 0
     integer :: nlite_files = 0
     real(r8):: frqplt     = 3600.0_r8
     real    :: dtvec      = 1200.0
     real    :: headspeed  = 3.0
     real    :: stemlength = 3.0e3
     integer :: plttype    = 0
     integer :: pltorient  = 0
     integer :: vec_maxmrl = maxgrds
     real    :: zplot_min  = -1.0
     real    :: zplot_max  = -1.0
     integer :: mapcolor   = 13
     integer :: llcolor    = 13
     integer :: ncolortabs = 0

     character(pathlen) :: pltname     = 'gmeta'
     character(10)      :: prtval_size = 'medium'

!!    THE LIST OF FILES TO PLOT FROM

     character(pathlen) ::  plt_files(maxpltfiles) = ' '
     character(pathlen) :: mavg_files(maxpltfiles) = ' '
     character(pathlen) :: davg_files(maxpltfiles) = ' '
     character(pathlen) :: lite_files(maxpltfiles) = ' '

!!    THE ARRAYS OF FIELDS TO PLOT

     type(oname_colortab) :: colortabs(maxnplt)
     type(oname_plot)     :: plotspecs(maxnplt)

!!    NCAR DCMIP 2012 PARAMETERS

     integer :: test_case

  End Type oname_vars

  type (oname_vars) :: nl
  integer           :: numcf =  0

  character(40)              :: cmdlne_runtype = ' '
  character(99), allocatable :: cmdlne_fields(:)

End Module oname_coms
