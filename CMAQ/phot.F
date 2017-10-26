
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header: /project/yoj/arc/CCTM/src/phot/phot_inline/phot.F,v 1.7 2011/10/21 16:11:28 yoj Exp $

! what(1) key, module and SID; SCCS file; date and time of last delta:
! %W% %P% %G% %U%

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE PHOT ( )

      use phot_mod
      USE UTILIO_DEFN
      USE CGRID_SPCS
      USE AERO_PHOTDATA
      USE RXNS_DATA,   only: nphotab
      use hrdata,      only: rj
      use mem_radiate, only: nadd_rad
      use cmaq_cloud
      use mem_ijtabs,  only: jtab_w, jtw_prog, istp
      use mem_grid,    only: mza, mwa, glatw, lpw
      use misc_coms,   only: current_time, time8p, dtlong, mstp
      use consts_coms, only: r8
      use mem_mclat,   only: mclat_spline
      use oname_coms,  only: nl

      implicit none

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV

      integer :: jdate      ! current Julian date (YYYYDDD)
      integer :: jday       ! current Julian day (DDD)
      real    :: jyear      ! year
      integer :: iyear      ! year
      REAL    :: JDSTRT     ! current Julian day (DDD)
      integer :: jtime      ! current time (HHMMSS)
      real    :: jyfreq
      REAL    :: CURRHR     ! current GMT hour
      REAL    :: JULIAN_DAY ! time of year [days]

      CHARACTER( 16 ) :: VARNM
      CHARACTER( 60 ) :: XMSG
      CHARACTER( 16 ) :: PNAME = 'PHOT'

      integer :: j, iw, iwl, nlays_max
      integer :: koff, nlays, nlays_model

      real, parameter :: twopi = 2.0 * pi
      integer, external :: julday, timetohhmmss

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
         LOGDEV   = INIT3()

         allocate( rj( mza, nphotab, mwa ) ) ; rj = 0.0

         CALL LOAD_CSQY_DATA( )

         CALL LOAD_OPTICS_DATA( )

         nlays_max = mza + nadd_rad
         CALL INIT_AERO_DATA( nlays_max )

         call cmaq_cld_optics_init( )

         ! set pointers to species O3 and NO2 in CGRID

         VARNM = 'O3'
         LGC_O3 = INDEX1( VARNM, N_GC_SPC, GC_SPC )
         IF ( LGC_O3 .LE. 0 ) THEN
            XMSG = 'Could not find ' // VARNM // 'in species table'
            CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT3 )
         END IF

         VARNM = 'NO2'
         LGC_NO2 = INDEX1( VARNM, N_GC_SPC, GC_SPC )
         IF ( LGC_NO2 .LE. 0 ) THEN
            XMSG = 'Could not find ' // VARNM // 'in species table'
            CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT3 )
         END IF

         allocate( actinic_flux_sfc( nwl, mwa ) )
         actinic_flux_sfc = 0.0

         allocate( cloud_frac_2d( mwa ) )
         cloud_frac_2d = 0.0

         allocate( sinlats( mwa ) )
         allocate( coslats( mwa ) )

         sinlats(1) = 0.0
         coslats(1) = 0.0
         do iw = 2, mwa
            sinlats(iw) =  SIN( PI180 * glatw(iw) )
            coslats(iw) =  cos( PI180 * glatw(iw) )
	enddo

         allocate( conv_wm2( nwl ) )
         allocate( srayl   ( nwl ) )

         ! FSB Set up conversion factor for
         ! [photons / ( cm**2 s) ] to [Watts / m**2 ]
         ! THE 1.0E13 FACTO IS 1.0E9 * 1.0 E4
         ! The 1.0e9 is for the wavelength [ nm ] -> [ m ]
         ! The 1.0e4 is for the area [ cm **2 ] -> [ m**2 ]
         DO IWL = 1, NWL
            CONV_WM2( IWL ) = 1.0E13 * ( PLANCK_C * LIGHT_SPEED ) / WAVELENGTH( IWL )
         END DO

         ! get molecular scattering cross sections
         CALL GETSRAY ( NWL, WAVELENGTH, SRAYL )

      END IF  ! firstime

      ! Check whether it is time to update photolysis rates

      if ( (istp == 1 .and. mod(time8p, nl%photfrq) < dtlong) .or.
     &     (istp == 1 .and. mstp == 0) ) then

         ! Calculate current hour in GMT and julian day

         currhr     = current_time%time / 3600.0_r8
         jday       = julday( current_time%month, current_time%date, current_time%year )
         julian_day = real( jday )
         jyear      = real( current_time%year )
         jdstrt     = julian_day
         jdate      = current_time%year * 1000 + jday
         jtime      = timetohhmmss( current_time%time )
         iyear      = current_time%year

         if ( (mod(iyear,400) == 0) .or. (mod(iyear,100) /= 0 .and. mod(iyear,4) == 0) ) then
            jyfreq = twopi / 366.0
         else
            jyfreq = twopi / 365.0
         endif
 
         ! Interpolate Mclatchy soundings between summer and winter values,
         ! and prepare spline coefficients for interpolation by latitude.

         call mclat_spline(jday)

         ! Calculate the solar information.

         CALL SOLEFM3 ( JYEAR, JDSTRT, EQT, SINDEC, COSDEC, SOLDIST )

         RSQDI = 1.0 / (SOLDIST * SOLDIST)

         ! MAIN loop over all rows and columns

         !$omp parallel do private(iw,koff,nlays,nlays_model)
         do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

            ! K index offset for radiation column arrays
            koff = lpw(iw) - 1

            nlays       = mza - koff + nadd_rad
            nlays_model = mza - koff

            call phot1( iw, currhr, julian_day, jyfreq,
     &                  nlays, nlays_model, koff)

         enddo
         !$omp end parallel do

      endif

      end subroutine phot

!-----------------------------------------------------------------------

      subroutine phot1( iw, currhr, julian_day, jyfreq, 
     &                  nlays, nlays_model, koff)

! Function:  Calculates the photolysis rate constant to be used by the
!     chemical solver.  It calculates these rates at each gridcell using
!     codes adapted from JPROC.  Cloud correction now called within the
!     loops over MY-ROW & MY_COLS
!
! Preconditions: HGRD_INIT() called from PAR_INIT, which is called from
!     DRIVER
!
! Subroutines/Functions called: INIT3, M3EXIT, SUBHFILE, CGRID_MAP,
!     OPPHOT, LOAD_CSQY_DATA, LOAD_OPTICS_DATA, INITIALIZE_ALBEDO, 
!     GET_PHOT_MET, UPDATE_SUN, GET_ALBEDO, GET_DROPLET_OPTICS, 
!     GET_ICE_OPTICS, GET_AGGREGATE_OPTICS, CLEAR_HYDROMETEOR_OPTICS, 
!     GET_AERO_DATA, O3TOTCOL, and NEW_OPTICS
!
! Revision History.
!     Started 10/08/2004 with existing PHOT and JPROC coded by
!         Dr. Francis S. Binkowski
!         Carolina Environmental Program
!         University of North Carolina at Chapel Hill
!         email: frank_binkowski@unc.edu
!     August 2005, Sarav Arunachalam, CEP, UNC-CH
!       - Minor revisions while integrating with CMAQ
!       - Error check for NPHOTS added (this version works only for SAPRC-99)
!       - Added creation of new file CTM_RJ_1 to write out RJ values
!         for O3 and NO2 (both clear sky and cloud effects), and
!         ETOT_SFC, TAU_AERO, TAU_TOT and TAUO3_TOP values for 7 wavelengths
!     June 2007, David Wong
!       -- inline with CMAQ
!       - declare RJ as assumed shape array to match with the caller routine
!       - allow PE 0 only to open the output file
!       - output species: NO2_CLOUD and O3_CLOUD with AMISS value when all cells
!         are dark and JTIME_CHK = 0
!       - output species: NO2_CLOUD and O3_CLOUD with AMISS value when CLDATT is
!         0 and JTIME_CHK = 0
!     December 2007, Francis Binkowski
!         code has been modified to call the new on-line version that
!         has the cloud effects built in.  new photolysis routine to
!         replace PHOT in CMAQ
!     January 2008, Shawn Roselle
!       - reformatted for inclusion in CMAQ
!       - added additional 3-d photolysis rate diagnostic file
!       - moved code for opening the diagnostic files to a separate subroutine
!       - moved aerosol pointer evaluation to a FORTRAN module
!       - simplified code for writing the diagnostic file
!       - changed code to call NEW_OPTICS twice, once for clear sky and
!         another time for the cloudy fraction of the grid cell.  RJ's are
!         computed based on the cloud fraction weighting.
!      March 2011, Bill Hutzell
!       - enable wavelength dependent arrays to have an allocatable number
!         of wavelength bins
!       - added data structure and algorithm to compute a surface albedo that
!         depends on time and landuse catagory based on work by John Striecher
!         (AMAD/USEPA)
!       - revised writing to RJ1 file to include surface albedo
!       - moved photolysis and opacity data from CSQY module to an ASCII input
!         file
!       - added routine called LOAD_REF_DATA (inside the PHOT_MOD module) that i
!         reads this input file
!       - added call to a routine called AERO_PHOTDATA that returns opacity data
!         on the aerosol distribution
!       - revised NEW_OPTICS' arguments based on  aerosol redesign in CMAQ
!         version 5.0
!     March 29, 2011 S.Roselle
!       - Replaced I/O API include files with UTILIO_DEFN
!     07 Jul 14 B.Hutzell: replaced mechanism include file(s) with fortran module
!     26 Sep 14 B.Hutzell: 1) moved calculation of surface albedo to its own
!                             fortran module
!                          2) changed loading procedure for loading optical data;
!                             two files now used
!                          3) reading and calculation of met and geo data
!                             now acomplished by a fortran module
!                          4) changed description and accounting of cloud effects
!                             from 2D liquid water clouds to 3D resolved and subgrid 
!                             clouds with multi-phases of water
!                          5) inserted calculation of aerosol optical properites via
!                             fortran module to improve efficiency in radiative
!                             transfer solution
!                          6) moved the O3TOTCOL routine from the PHOT_MOD to simplify
!                             the NEW_OPTICS routine
!                          7) Several miscellaneous changes attempting to improve efficiency
!     June 10 15 J.Young: Modified diagnostic output timestamp to fix for other than one
!                         hour time steps.
!     Aug 12, 15 D. Wong: Replaced MYPE with IO_PE_INCLUSIVE for parallel I/O implementation

!----------------------------------------------------------------------

C...modules

      USE RXNS_DATA            ! chemistry varaibles and data
      USE CGRID_SPCS           ! CGRID species number and offsets
!     USE PCGRID_DEFN          ! get cgrid
      USE UTILIO_DEFN
      USE AERO_DATA            ! describes aerosol distribution
      USE PHOT_MOD             ! photolysis in-line module - inherits CSQY_DATA module
      USE AERO_PHOTDATA        ! arrays and routines for aerosol dimensions and refractive indices
!     USE PHOTOLYSIS_ALBEDO    ! surface albedo data and routines
!     USE PHOT_MET_DATA        ! Met and Grid data
 !    USE CLOUD_OPTICS         ! data and routines for optics of cloud hydrometeors

      use hrdata,      only: rj
      use mem_grid,    only: mza, nsw_max
      use mem_radiate, only: nadd_rad
      use mem_grid,    only: mza, glatw, glonw, zt, zm, lpw, dzt
      use cgrid_defn
      use mem_basic,   only: tair, press, rho, sh_v
      use consts_coms, only: r8, t00, alvl, eps_vap, cp
      use mem_mclat,   only: rad_mclat
      use cmaq_cloud
      use micro_coms,  only: ncat, rxmin, emb0, reffcof, pwmasi,
     &                       dmncof, jhabtab, emb2, zfactor_ccn
      use mem_micro,   only: cldnum
      use mem_cuparm,  only: kcutop, kcubot, qwcon, conprr, iactcu
      use therm_lib,   only: rhovsl

      IMPLICIT NONE

!...arguments

      integer, intent(in) :: iw
      REAL,    intent(in) :: CURRHR          ! current GMT hour
      REAL,    intent(in) :: JULIAN_DAY      ! time of year [days]
      real,    intent(in) :: jyfreq
      integer, intent(in) :: koff, nlays, nlays_model

!     REAL RJ( NCOLS,NROWS,NLAYS, NPHOTAB )
!     REAL,    INTENT( OUT ) :: RJ( :,:,:,: )  ! gridded J-values  (1/min units)

!     REAL CGRID( NCOLS,NROWS,NLAYS, * )  ! Conc array
!     REAL, SAVE, POINTER :: CGRID( :,:,:,: ) ! species concentrations

!...parameters

      LOGICAL, PARAMETER :: CLDATT = .TRUE.   ! include cloud attenuation

      REAL, PARAMETER :: DENS_CONV = ( 1.0E+03 * AVO / MWAIR ) * 1.0E-06  ! convert from kg/m**3 to #/cc
      REAL, PARAMETER :: PPM_MCM3  = 1.0E-06  ! convert from ppm to molecules / cc mol_Spec/mol_Air = ppm * 1E-06
      REAL, PARAMETER :: PRES_CONV = 1.0 / STDATMPA ! conversion factor Pa to atm
      REAL, PARAMETER :: ZTOA      = 50.0E3   ! height of top of atmosphere [ m ] (=50km)
                                              ! based a 2005 WRF model Documentation
                                                 
      REAL, PARAMETER :: EPSLON  = 1.0E-30  ! Small number

      real, parameter :: mwo3    = 48.0    ! molecular weight of ozone
      real, parameter :: o3_conv = mwair / mwo3 * 1.e6 ! ozone mixing ratio to ppmV
      real, parameter :: cosz_min= cos( 88.0 * PI180 )
      REAL, PARAMETER :: TWOPI   = 2.0 * PI

!...external functions: none

!...local variables

!     CHARACTER(   3 ), ALLOCATABLE, SAVE :: WLTXT( : )
!     CHARACTER(  16 )                    :: VARNM
!     CHARACTER(  16 ), SAVE              :: PNAME = 'PHOT'
!     CHARACTER(  16 ), SAVE              :: CTM_PHOTDIAG = 'CTM_PHOTDIAG'

!     CHARACTER(  80 ) :: VARDESC  ! environment variable description
!     CHARACTER( 240 ) :: XMSG = ' '

!     INTEGER,    SAVE :: TSTEP    ! output timestep in sec

!     INTEGER ESTAT                ! status from environment var check
      INTEGER IPHOT                ! photolysis rate loop index
!     INTEGER ROW
!     INTEGER COL
!     INTEGER LEV
!     INTEGER SPC
!     INTEGER IWL
      INTEGER L
!     INTEGER V, N, MODE
!     LOGICAL JTIME_CHK            ! To check for JTIME to write RJ values
!     INTEGER ODATE                ! output date
!     INTEGER OTIME                ! output time 
!     INTEGER ALLOCSTAT

      integer :: k, nlays_chem

!      integer :: jdate      ! current Julian date (YYYYDDD)
!      integer :: jday       ! current Julian day (DDD)
!      real    :: jyear      ! year
!      integer :: iyear      ! year
!      REAL    :: JDSTRT     ! current Julian day (DDD)
!      integer :: jtime      ! current time (HHMMSS)
!      real    :: jyfreq

      real    :: dens

      REAL CURRHR_LST      ! local standard time at each grid cell
!      REAL CTOP            ! cloud top in single dimension
!      REAL CBASE           ! cloud base in single dimension
!      REAL ZLEV            ! height in single dimension
!      REAL ZEN             ! cosine of zenith angle
!      REAL SINLAT          ! sine of latitude
!      REAL COSLAT          ! cosine of latitude
      REAL ZSFC            ! surface height (msl) [ m ]
      REAL COSZEN          ! working cosine of the solar zenith angle
      REAL SINZEN          ! working sine of the solar zenith angle
!      REAL LATCR           ! local latitude
!      REAL LONCR           ! local longitude
!      REAL OWATER_FRAC     ! Open water fraction
!      REAL SNOW_FRAC       ! Snow fractional coverage
!      REAL SEAICE_FRAC     ! Sea Ice fraction
!      REAL RES_SKY_REFLECT ! reflection coefficient based on resolved sky
!      REAL RES_SKY_TRANS   ! diffuse transmission coefficient based on resolved sky
!      REAL RES_SKY_TRANSD  ! direct transmission coefficient based on resolved sky
      REAL MSCALE          ! combined factor to scale ppm to Molecules / cm**3
                           ! and correct for ambient temperaure and pressure

      real :: sctop( nwl )
      real :: extop( nwl )
      real :: astop( nwl )

!     REAL, ALLOCATABLE, SAVE :: TAUO3_TOP( : )    ! optical depth of ozone above model domain
!     REAL, ALLOCATABLE, SAVE :: TAU_RAY  ( : )    ! Rayleigh optical depth above model domain

!     REAL :: ETOT_SFC       ( nwl )  ! total downward irradiance at sfc [ Watts / m**2  ]
!     REAL :: TAUC_AERO( nlays,nwl )  ! aerosol optical depth at layer bottom
!     REAL :: TAU_TOT  ( nlays,nwl )  ! total optical depth at layer bottom
!     REAL :: TAU_CLOUD( nlays,nwl )  ! cloud optical depth at layer bottom
!     REAL :: SSA            ( nwl )  ! aerosol single scattering albedo, column average

      real :: alb        ( nsw_max,nwl )  ! surface albedo
      real :: alb_dif    ( nsw_max,nwl ) 

! FSB new arrays for new on-line cloud version

!     REAL :: LWC    ( nlays )      ! cloud liquid water content [ g/m**3 ]
!     REAL :: RWC    ( nlays )      ! rain water content [ g/m**3 ]
!     REAL :: IWC    ( nlays )      ! ice liquid water content [ g/m**3 ]
!     REAL :: SWC    ( nlays )      ! snow content [ g/m**3 ]
!     REAL :: GWC    ( nlays )      ! graupel content [ g/m**3 ]
      REAL :: CLDFRAC( nlays )      ! fractional cloud cover
      REAL :: BLKPRS ( nlays )      ! Air pressure in [ Pa ]
      REAL :: BLKTA  ( nlays )      ! Air temperature [ K ]
      REAL :: BLKDENS( nlays )      ! Air density  [ molecules / m**3 ]
      REAL :: BLKZH  ( nlays )      ! layer half-height [ m ]
      REAL :: BLKO3  ( nlays )      ! O3 concentration [ molecules/cm**3 ]
      REAL :: BLKNO2 ( nlays )      ! NO2 concentration [ molecules/cm**3 ]
      REAL :: BLKZF  ( nlays+1 )    ! layer full-height [ m ]
      REAL :: BLKDZ  ( nlays )
      REAL :: BLKSHV ( nlays )
      REAL :: BLKRJ  ( nlays,nphotab ) ! photolysis rates

!     LOGICAL :: NEW_PROFILE    ! Has atmospheric temperature and density profile changed?
      
!...Variables for diagnostic outputs

!      REAL, ALLOCATABLE, SAVE :: TOC( :,: )              ! total ozone column
!      REAL, ALLOCATABLE, SAVE :: TRANSMIS_DIFFUSE( :,: ) ! diffuse transmission coefficient at surface
!      REAL, ALLOCATABLE, SAVE :: TRANSMIS_DIRECT( :,: )  ! direct transmission coefficient at surface
!      REAL, ALLOCATABLE, SAVE :: REFLECT_COEFF( :,: )    ! reflection coefficient at top of atmosphere
!      REAL, ALLOCATABLE, SAVE :: ETOT_SFC_WL ( :,:,: )   ! total downward irradiance at sfc [ Watts / m**2  ]
!      REAL, ALLOCATABLE, SAVE :: TAU_AERO_WL ( :,:,: )   ! total aerosol optical depth
!      REAL, ALLOCATABLE, SAVE :: TAU_CLOUD_WL( :,:,: )   ! total cloud optical depth
!      REAL, ALLOCATABLE, SAVE :: CLR_TRANSMISSION( :,: ) ! diffuse transmission coefficient of clouds
!      REAL, ALLOCATABLE, SAVE :: CLR_REFLECTION  ( :,: ) ! reflection coefficient of cloud
!      REAL, ALLOCATABLE, SAVE :: CLR_TRANS_DIRECT( :,: ) ! direct transmission coefficient of clouds
!#ifdef phot_debug      
!      REAL, ALLOCATABLE, SAVE :: ASY_CLOUD_WL( :,:,: ) ! columm average of cloud asymmetry factor
!      REAL, ALLOCATABLE, SAVE :: SSA_CLOUD_WL( :,:,: ) ! columm average of cloud single scattering albedo
!#endif      
!      REAL, ALLOCATABLE, SAVE :: TAU_TOT_WL  ( :,:,: ) ! total optical depth
!      REAL, ALLOCATABLE, SAVE :: TAUO3_TOP_WL( :,:,: ) ! optical depth of ozone above model domain
!
!      REAL, ALLOCATABLE, SAVE :: AERO_SSA  ( :,:,:,: ) ! aerosol single scattering albedo
!      REAL, ALLOCATABLE, SAVE :: AERO_ASYM ( :,:,:,: ) ! aerosol asymmetry factor
!      REAL, ALLOCATABLE, SAVE :: TAU       ( :,:,:,: ) ! optical depth
!      REAL, ALLOCATABLE, SAVE :: TAU_AERO  ( :,:,:,: ) ! aerosol optical depth
!      REAL, ALLOCATABLE, SAVE :: ACTINIC_FX( :,:,:,: ) ! net actinic flux [watts/m**2]

!     INTERFACE
!        SUBROUTINE O3TOTCOL ( LATITUDE, LONGITUDE, JDATE, OZONE )
!           INTEGER, INTENT( IN )    :: JDATE      ! Julian day of the year (yyyyddd)
!           REAL,    INTENT( IN )    :: LATITUDE   ! latitude of point on earth's surface
!           REAL,    INTENT( IN )    :: LONGITUDE  ! longitude of point on earth's surface
!           REAL,    INTENT( INOUT ) :: OZONE      ! total column ozone [DU]
!        END SUBROUTINE O3TOTCOL
!     END INTERFACE

      integer, external :: julday, timetohhmmss

      real, parameter :: rliqland  = 14.0
      real, parameter :: rliqocean =  8.0

      integer :: jhcat(mza,ncat) ! hydrom category table with ice habits

      real :: rhov (mza)        ! vapor density [kg_vap/m^3]
      real :: rx   (mza,ncat)   ! hydrom bulk spec dens [kg_hyd/kg_air]
      real :: cx   (mza,ncat)   ! hydrom bulk number [num_hyd/kg_air]
      real :: emb  (mza,ncat)   ! hydrom mean particle mass [kg/particle]

      real :: r_ef, dmean, watp, twc, prate, tc, rh

      real    :: frac(mza)
      logical :: dosnow

      real :: tauclds(nwl, nlays)
      real :: ssaclds(nwl, nlays)
      real :: asmclds(nwl, nlays)

      integer :: mc, mcat, ih, krad, ns, nt, ib

! Array kradcat maps RAMS/OLAM microphysics hydrometeor categories to those
! represented in the cloud optics tables according to the following numbering:

!     Lookup table category             OLAM Microphysics
! ----------------------------------------------------------------
!  1:   cloud drops                 1.  cloud drops
!  2:   drizzle                     2.  rain
!  3:   rain                        3.  pristine ice columns
!  4:   hail                        4.  snow columns
!  5:   aggregates                  5.  aggregates
!  6:   hollow columns              6.  graupel
!  7:   solid columns               7.  hail
!  8:   hexagonal plates            8.  drizzle
!  9:   rosettes                    9.  pristine ice hexagonal plates
!                                  10.  pristine ice dendrites
!                                  11.  pristine ice needles
!                                  12.  pristine ice rosettes
!                                  13.  snow hexagonal plates
!                                  14.  snow dendrites
!                                  15.  snow needles
!                                  16.  snow rosettes

                                          ! 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
      integer, parameter :: kradcat(16) = (/1,3,6,6,5,4,4,2,8, 8, 7, 9, 8, 8, 7, 9/)

! ----------------------------------------------------------------------

               RJ( :,:,iw ) = 0.0

               actinic_flux_sfc( :,iw ) = 0.0
               cloud_frac_2d( iw ) = 0.0

               CURRHR_LST = CURRHR + glonw( iw )  / 15.0

               CALL GETZEN2 ( CURRHR_LST , SINLATS(iw), COSLATS(iw),
     &                     SINDEC, COSDEC, EQT, COSZEN )

               ! If column is dark, skip to next column
               if (coszen < cosz_min) return

               do k = mza-1, lpw(iw), -1
                  if (press(k,iw) > 100.d2) exit
               enddo
               nlays_chem = (k+1) - koff

!               IF ( COSZEN .LE. 0.0 ) THEN
!                  RJ( COL,ROW, :,: ) = 0.0
!   
!                  IF ( JTIME_CHK ) THEN
!                     TOC( COL,ROW )             = 0.0                            
!                     ETOT_SFC_WL ( COL,ROW, : ) = 0.0
!                     TAUO3_TOP_WL( COL,ROW, : ) = 0.0
!                     TAU_AERO_WL ( COL,ROW, : ) = 0.0
!                     TAU_CLOUD_WL( COL,ROW, : ) = 0.0
!#ifdef phot_debug                   
!                     SSA_CLOUD_WL( COL,ROW, : ) = 0.0
!                     ASY_CLOUD_WL( COL,ROW, : ) = 0.0
!#endif                   
!                     TAU_TOT_WL( COL,ROW, : )   = 0.0
!                     TAU       ( COL,ROW, :,: ) = 0.0
!                     TAU_AERO  ( COL,ROW, :,: ) = 0.0
!                     AERO_SSA  ( COL,ROW, :,: ) = 0.0
!                     AERO_ASYM ( COL,ROW, :,: ) = 0.0
!                     ACTINIC_FX( COL,ROW, :,: ) = 0.0
!                   
!                     TRANSMIS_DIFFUSE( COL,ROW ) = 0.0
!                     TRANSMIS_DIRECT ( COL,ROW ) = 0.0
!                     REFLECT_COEFF   ( COL,ROW ) = 0.0
!                     CLR_TRANSMISSION( COL,ROW ) = 0.0
!                     CLR_TRANS_DIRECT( COL,ROW ) = 0.0
!                     CLR_REFLECTION  ( COL,ROW ) = 0.0 
!                  END IF
!
!                  CYCLE LOOP_COLS
!
!               END IF

!...initialize BLKRJ using F90 array operations.

               BLKRJ = 0.0
!              BLKRJ_ACM = 0.0

!...Set height of lowest level to zero

               BLKZF( 1 ) = 0.0

               ZSFC   = zm( lpw(iw)-1 ) ! surface height [m]
               SINZEN = SQRT( 1.0 - COSZEN * COSZEN ) ! sine of zenith angle

!...local latitude and longitude

!               LATCR = LAT( COL,ROW )
!               LONCR = LON( COL,ROW )

!...get total ozone column based on OMI observations
!               CALL O3TOTCOL ( LAT( COL,ROW ), LON( COL,ROW ), JDATE, TOTAL_O3_COLUMN )
            
!               IF ( USE_ACM_CLOUD .OR. CLDATT ) THEN
!                  OWATER_FRAC = MAX( ( 1.0 - SEAICE( COL,ROW ) ), 0.0 )
!     &                        * WATER_FRACTION( COL,ROW )
!                  SEAICE_FRAC = SEAICE( COL,ROW ) * WATER_FRACTION( COL,ROW )
!                  SNOW_FRAC   = SNOCOV( COL,ROW )
!                  COL_CLOUD   = PHOT_COL 
!                  ROW_CLOUD   = PHOT_ROW
!               END IF

               cldfrac = 0.0
               tauclds = 0.0
               ssaclds = 0.0
               asmclds = 0.0

!...loop over vertical layers ambient air conditions and gas concentration
               do k = lpw(iw), mza
                  l = k - koff

!...Fetch the grid cell ambient data at each layer.

                  BLKTA  ( L )   = tair ( k, iw ) ! temperature [K]
!                 BLKPRS ( L )   = press( k, iw ) * PRES_CONV ! [atmospheres]
                  BLKDENS( L )   = rho  ( k, iw ) * DENS_CONV ! [molecules / cm**3]
                  BLKZH  ( L )   = zt   ( k ) ! mid layer height [m]
                  BLKZF  ( L+1 ) = zm   ( k ) ! full layer height [m]

!...set scale factor for [ppm] -> [molecule / cm**3]
!...  To go from ppm to molecule/cc:
!...  molecule/cc = ppm *  1.0E-06 * DENS (given in molecule/cc)

                  MSCALE = BLKDENS( L ) * PPM_MCM3

!...fetch ozone and no2 and convert to [ molecules / cm **3 ]
!...  and adjust the volume for ambient temperature and pressure.

                  BLKO3 ( L ) = max(CGRID( k, iw, LGC_O3  ), 0.0) * MSCALE
                  BLKNO2( L ) = max(CGRID( k, iw, LGC_NO2 ), 0.0) * MSCALE

                  rhov   ( k ) = max(0.,sh_v(k,iw)) * rho(k,iw)
                  
               END DO           ! loop on layers ambient conditions and gases

!..calculate needed aerosol properties in column

!              IF ( CORE_SHELL ) THEN
                  CALL GET_AERO_DATA ( iw, nlays_model, koff )
!              ELSE
!                 CALL AERO_OPTICS_INTERNAL( COL,ROW, NLAYS, CGRID )
!              END IF

!...Fill ozone column and any extra radiation layers at model top

               blkprs(mza-koff) = press(mza,iw)
               call rad_mclat( iw, nlays, koff, glatw(iw), blkdens, blkprs, blkshv, blkta,
     &                         blko3, blkzf(2:nlays+1), blkzh, blkdz, o3col=.false. )

               sctop(1:nwl) = aero_scat_coef(mza-koff,1:nwl) / blkdens(mza-koff)
               extop(1:nwl) = aero_exti_coef(mza-koff,1:nwl) / blkdens(mza-koff)
               astop(1:nwl) = aero_asym_fac (mza-koff,1:nwl)
               
               do l = 2, nlays+1
                  blkzf(l) = blkzf(l) - zsfc
               enddo

               do l = mza-koff+1, nlays
                  dens         = blkdens(l)
                  blkdens(l)   = blkdens(l)   * dens_conv

                  mscale       = blkdens(l)   * ppm_mcm3
                  blko3  (l)   = blko3  (l)   * o3_conv * mscale / dens
                  blkno2 (l)   = 0.5e-3       * mscale

                  aero_scat_coef(l,1:nwl) = sctop(1:nwl) * blkdens(l)
                  aero_exti_coef(l,1:nwl) = extop(1:nwl) * blkdens(l)
                  aero_asym_fac (l,1:nwl) = astop(1:nwl)
               enddo

! Fill arrays rx, cx, and emb with hydrometeor properties

               call cloudprep_rad(iw,lpw(iw),mcat,jhcat,rhov,rx,cx,emb)

! Compute fractional cloudiness for the resolved microphysics moisture fields. 
! The cloud fraction estimated here will only be applied later in this routine
! if there are any resolved hydrometeors or subgrid cumulus

               call get_cloud_frac(iw, lpw(iw), frac)

! Get optical properties of resolved clouds

               do mc = 1, mcat
                  do k = lpw(iw), mza
                     krad = k - koff

                     if (rx(k,mc) >= rxmin(mc) .and. emb(k,mc) >= emb0(mc)) then

                        ! Set lower bound on frac(k) because there is condensate
                        cldfrac(krad) = max(frac(k), 0.1)

                        ! lookup table category
                        ih = jhcat(k,mc)
                        l  = kradcat(ih)

                        ! effective radius in microns
                        r_ef = 1.e6 * reffcof(ih) * emb(k,mc) ** pwmasi(ih)

                        ! ice or liquid water path in g/m^2
                        watp = rx(k,mc) * real(rho(k,iw)) * 1000. * dzt(k)
                        watp = watp / max( frac(k), 0.2 )

                        call lookup_cmaq_cld_optics( l, r_ef, watp, krad )

                     endif
                  enddo
               enddo

! Include the optical properties of subgrid convective clouds

               if (iactcu(iw) == 1) then

                  do k = kcubot(iw), kcutop(iw)
                     krad = k - koff

                     ! Set lower bound on frac(k) because there is condensate
                     cldfrac(krad) = max(frac(k), 0.1)

                     ! water path in g/m^2
                     watp = qwcon(k,iw) * real(rho(k,iw)) * 1000. * dzt(k)

                     tc = tair(k,iw) - t00

                     if (tc > -10.0) then

                        ! Add convective cloud water to cloud drops if warmer then -10C
                        ih = 1
                        l  = kradcat(ih)

                        ! effective radius in microns
                        r_ef = 1.e6 * reffcof(ih) 
     &                       * ( qwcon(k,iw) / (cldnum(iw) * zfactor_ccn(k)) ) ** pwmasi(ih)

                     else

                        ! Add convective cloud water to pristine ice. Diagnose habit
                        ! from temperature and humidity
                        rh = min( 1., rhov(k) / rhovsl(tc) )
                        ns = max( 1, nint(100. * rh) )
                        nt = max( 1, min(31,-nint(tc)) )
                        ih = jhabtab(nt,ns,1)
                        l  = kradcat(ih)

                        ! Mean maximum dimension of ice crystals as a function of T
                        ! (see Kristjansson et al., 2000, JGR)
                        dmean = 1030.7 * exp(0.05522*(tair(k,iw)-279.5))

                        ! Convert mean diameter to an effective radius using the 
                        ! microphysics power laws
                        r_ef = reffcof(ih) / dmncof(ih) * dmean

                     endif

                     call lookup_cmaq_cld_optics( l, r_ef, watp, krad )

                  enddo

! Also include the optical properties of convective rain/snow

                  if (conprr(iw) > 1.e-10) then

                     ! Estimate rain/snow from parameterization used in EPA's CMAQ model

                     prate = conprr(iw) * 3600.0 ! precip rate: mm / hr
                     twc   = 0.06 * prate**0.846 ! tot wat: g / m3

                     ! if entire cloud is below freezing, map to snow
                     dosnow = all(tair(kcubot(iw):kcutop(iw),iw) < 273.)

                     do k = lpw(iw), kcutop(iw)
                        krad = k - koff

                        tc = tair(k,iw) - t00

                        ! Set lower bound on frac(k) because there is condensate
                        cldfrac(krad) = max(frac(k), 0.1)

                        ! Water/ice path: g / m^2
                        watp  = twc * dzt(k)

                        if (dosnow) then

                           if (tc > 0.0) then
                              ! rain
                              mc = 2
                              ih = 2
                           else
                              ! snow
                              mc = 4
                              rh = min( 1., rhov(k) / rhovsl(tc) )
                              ns = max( 1, nint(100. * rh) )
                              nt = max( 1, min(31,-nint(tc)) )
                              ih = jhabtab(nt,ns,2)
                           endif

                        else

                           if (tc > -10.) then
                              ! rain
                              mc = 2
                              ih = 2
                           else
                              ! hail
                              mc = 7
                              ih = 7
                           endif

                        endif

                        ! cloud optics category
                        l  = kradcat(ih)

                        ! effective radius of convctive rain/snow
                        r_ef = 1.e6 * reffcof(ih) * emb2(mc) ** pwmasi(ih)

                        call lookup_cmaq_cld_optics( l, r_ef, watp, krad )

                     enddo

                  endif         ! convective rain/snow

               endif            ! convective clouds

! Combine cloud optical properties

               do k = lpw(iw), mza
                  krad = k - koff

                  if ( all(tauclds(1:nwl,krad) < 2.e-5) ) then
                     tauclds(1:nwl,krad) = 0.0
                     ssaclds(1:nwl,krad) = 1.0
                     asmclds(1:nwl,krad) = 0.0
                     cldfrac      (krad) = 0.0
                  else
                     do ib = 1, nwl
                        if (tauclds(ib,krad) > 1.e-9 .and. ssaclds(ib,krad) > 1.e-9) then
                           asmclds(ib,krad) = asmclds(ib,krad) / ssaclds(ib,krad)
                           ssaclds(ib,krad) = ssaclds(ib,krad) / tauclds(ib,krad)
                        endif
                     enddo
                  endif

               enddo

! Get optical properties of resolved clouds


!               IF ( CLDATT .AND. CFRAC_2D( COL,ROW ) .GT. 0.0 ) THEN
!                  DO L = 1, NLAYS
!                 
!                     IF ( CFRAC_3D( COL,ROW,L ) .GT. 0.0 ) THEN
!                        CLOUDS ( L )        = .TRUE.
!                        CLOUD_LAYERING( L ) = .TRUE.
!                        CLDFRAC( L )        = CFRAC_3D( COL,ROW,L )
!!... set hydrometeor concentrations for resolved cloud
!                        MSCALE   = 1.0E+3 * DENS ( COL,ROW,L )
!                        IWC( L ) = MSCALE * QI( COL,ROW,L )
!                        GWC( L ) = MSCALE * QG( COL,ROW,L )
!                        SWC( L ) = MSCALE * QS( COL,ROW,L )    
!                        LWC( L ) = MSCALE * QC( COL,ROW,L )
!                        RWC( L ) = MSCALE * QR( COL,ROW,L )
!                     ELSE
!                        CLOUDS ( L )        = .FALSE.
!                        CLOUD_LAYERING( L ) = .FALSE.
!                        CLDFRAC( L )        = 0.0
!                        IWC( L ) = 0.0
!                        GWC( L ) = 0.0
!                        SWC( L ) = 0.0    
!                        LWC( L ) = 0.0
!                        RWC( L ) = 0.0
!                     END IF
!                  END DO ! loop on layers clouds
!! get optical properties of resolved cloud hydrometeors
!                  CALL GET_DROPLET_OPTICS( NLAYS, BLKTA, OWATER_FRAC, SEAICE_FRAC, SNOW_FRAC, LWC )     
!                  CALL GET_ICE_OPTICS( NLAYS, BLKTA, IWC )
!                  CALL GET_AGGREGATE_OPTICS( NLAYS, RWC, SWC, GWC )     
!               ELSE
!                  CLOUDS         = .FALSE.
!!                 CLOUD_LAYERING = .FALSE.
!                  CLDFRAC        = 0.0
!  hydrometeor concentrations
!                  LWC = 0.0
!                  IWC = 0.0
!                  RWC = 0.0
!                  SWC = 0.0
!                  RWC = 0.0
!                 CALL CLEAR_HYDROMETEOR_OPTICS()
!               END IF
            
! set surface albedo

               call phot_albedo(iw, coszen, currhr_lst, julian_day, jyfreq, alb, alb_dif)

!               FORALL ( IWL = 1:NWL )
!                  ALB( IWL ) = SURFACE_ALBEDO( IWL, COL,ROW )
!               END FORALL
!               alb(:) = 0.2

!...calculate resolved-sky photolysis rates at all layers:

!              NEW_PROFILE    = .TRUE. 
!              ONLY_SOLVE_RAD = .FALSE. 

               CALL NEW_OPTICS ( iw, NLAYS, nlays_model, nlays_chem,
     &                           BLKTA, BLKDENS, BLKZF,
     &                           BLKO3, BLKNO2, alb, 
     &                           ZSFC, COSZEN, SINZEN,
     &                           CLDFRAC, BLKRJ,
     &                           tauclds, ssaclds, asmclds )

!...load diagnostic file arrays

!               IF ( JTIME_CHK  ) THEN
!                  TOC( COL,ROW ) = REAL( TOTAL_O3_COLUMN )
!                  TRANSMIS_DIFFUSE( COL,ROW ) = TRANSMISSION
!                  TRANSMIS_DIRECT( COL,ROW )  = TRANS_DIRECT
!                  REFLECT_COEFF( COL,ROW )    = REFLECTION
!                  FORALL( IWL = 1:NWL )
!                     ETOT_SFC_WL ( COL,ROW,IWL ) = IRRADIANCE( 1,IWL )
!                     TAUO3_TOP_WL( COL,ROW,IWL ) = TAUO3_TOP( IWL )
!                     TAU_AERO_WL ( COL,ROW,IWL ) = TAUC_AERO( 1,IWL )
!                     TAU_TOT_WL  ( COL,ROW,IWL ) = TAU_TOT  ( 1,IWL )
!                     TAU_CLOUD_WL( COL,ROW,IWL ) = TAU_CLOUD( 1,IWL )
!#ifdef phot_debug                      
!                     SSA_CLOUD_WL( COL,ROW,IWL ) = AVE_SSA_CLD  ( IWL )
!                     ASY_CLOUD_WL( COL,ROW,IWL ) = AVE_ASYMM_CLD( IWL )
!#endif                      
!                  END FORALL
!                  FORALL ( LEV = 1:NLAYS, IWL = 1:NWL  )
!                     ACTINIC_FX( COL,ROW,LEV,IWL ) = ACTINIC_FLUX( LEV,IWL )
!                  END FORALL
!                  
!                  DO L = 1, N_DIAG_WVL
!                     IWL = DIAG_WVL( L )
!                     FORALL ( LEV = 1:NLAYS )
!                        TAU     ( COL,ROW,LEV,L ) = TAU_TOT  ( LEV,IWL )
!                        TAU_AERO( COL,ROW,LEV,L ) = TAUC_AERO( LEV,IWL )
!                     END FORALL
!                     FORALL ( LEV = 1:NLAYS, AERO_EXTI_COEF( LEV,IWL ) .GT. EPSLON )
!                        AERO_SSA ( COL,ROW,LEV,L ) = AERO_SCAT_COEF( LEV,IWL )
!     &                                             / AERO_EXTI_COEF( LEV,IWL )
!                        AERO_ASYM( COL,ROW,LEV,L ) = AERO_ASYM_FAC( LEV,IWL )
!                     END FORALL
!                     FORALL ( LEV = 1:NLAYS, AERO_EXTI_COEF( LEV,IWL ) .LE. EPSLON )
!                        AERO_SSA ( COL,ROW,LEV,L ) = 1.0
!                        AERO_ASYM( COL,ROW,LEV,L ) = 0.0
!                     END FORALL
!                  END DO
!               END IF

!Set Photolysis rates to resolved sky values
!               FORALL ( L = 1:NLAYS, IPHOT = 1:NPHOTAB )
!                  RJ( COL,ROW, L,IPHOT ) = 60.0 * BLKRJ_RES( L,IPHOT ) 
!               END FORALL  ! Loop on layers and NPHOTAB

               do l = 1, nlays_chem
                  k = l + koff

                  DO IPHOT = 1, NPHOTAB
                     RJ( k, IPHOT, iw ) = 60.d0 * BLKRJ( L,IPHOT )
                  END DO        ! Loop on PHOT

               END DO           ! Loop on layers

!               IF ( USE_ACM_CLOUD ) THEN
!                  IF ( ACM_CLOUDS( COL,ROW ) .GT. 0.0 ) THEN
!! save resolved sky reflection and transmission coefficients for possible latter use
!                     RES_SKY_REFLECT = REFLECTION
!                     RES_SKY_TRANS   = TRANSMISSION
!                     RES_SKY_TRANSD  = TRANS_DIRECT
!!...find the highest layer of the sub-grid (convective) cloud
!                     DO LEV = NLAYS, 1, -1
!                        IF ( ACM_CFRAC( LEV, COL,ROW ) .GT. 0.0 ) EXIT
!                     END DO
!!...replace the lower layers with sub-grid cloud properties
!                     DO L = 1, LEV
!                        SWC( L ) = 0.0
!                        IF ( ACM_CFRAC( L,COL,ROW ) .GT. 0.0 ) THEN
!                           CLOUDS ( L ) = .TRUE.
!                           CLDFRAC( L ) = 1.0
!                           MSCALE = 1.0E+3 * DENS ( COL,ROW, L )
!                           LWC( L ) = MSCALE * ACM_QC( L,COL,ROW )
!                           IWC( L ) = MSCALE * ACM_QI( L,COL,ROW )
!                           RWC( L ) = MSCALE * ACM_QR( L,COL,ROW )
!                           GWC( L ) = MSCALE * ACM_QG( L,COL,ROW )
!                        ELSE
!                           CLOUDS( L )  = .FALSE.
!                           CLDFRAC( L ) = 0.0
!                           LWC( L )     = 0.0
!                           IWC( L )     = 0.0
!                           RWC( L )     = 0.0
!                           GWC( L )     = 0.0
!                        END IF
!                        CLOUD_LAYERING( L ) = .FALSE.
!                     END DO
!!                    write(logdev,*)'ACM cloud present fraction, cloud lwc(lev),iwc(lev),rwc(1),gwc(1) = ',
!!     &              ACM_CLOUDS( COL,ROW ),lwc(lev),iwc(lev),rwc(1),gwc(1)
!                 
!! get optical properties of of subgrid cloud hydrometeors
!                     CALL GET_DROPLET_OPTICS( LEV, BLKTA, OWATER_FRAC, SEAICE_FRAC, SNOW_FRAC, LWC )     
!                     CALL GET_ICE_OPTICS( LEV, BLKTA, IWC )     
!                     CALL GET_AGGREGATE_OPTICS( LEV, RWC, SWC, GWC )     
!
!!...calculate the acm-cloud photolysis rates for all layers:
!                     NEW_PROFILE = .FALSE.   
!                     CALL NEW_OPTICS ( JDATE, JTIME, NLAYS,
!     &                                 BLKTA, BLKPRS, BLKDENS, BLKZH, BLKZF,
!     &                                 BLKO3, BLKNO2,
!     &                                 ZSFC, COSZEN, SINZEN, RSQD,
!     &                                 NEW_PROFILE, CLOUDS, CLDFRAC,
!     &                                 BLKRJ_ACM, TAUC_AERO, TAU_TOT, TAUO3_TOP,
!     &                                 TAU_RAY, SSA, TAU_CLOUD, TOTAL_O3_COLUMN )
!
!!...load diagnostic file arrays
!!...compute a cloud-fraction weighted average of ETOT_SFC and TAU_TOT
!!...  note that both TAUC_AERO and TAUO3_TOP are the same for clear and
!!...  cloudy regions
!                     MSCALE = MAX( 1.0 - ACM_CLOUDS( COL,ROW ), 0.0 )
!
!                     IF ( JTIME_CHK ) THEN
!
!                        TRANSMIS_DIRECT( COL,ROW )  = MSCALE * TRANSMIS_DIRECT( COL,ROW )
!     &                                              + ACM_CLOUDS( COL,ROW ) * TRANS_DIRECT
!                        TRANSMIS_DIFFUSE( COL,ROW ) = MSCALE * TRANSMIS_DIFFUSE( COL,ROW )
!     &                                              + ACM_CLOUDS( COL,ROW ) * TRANSMISSION
!                        REFLECT_COEFF( COL,ROW )    = MSCALE * REFLECT_COEFF( COL,ROW )
!     &                                              + ACM_CLOUDS( COL,ROW ) * REFLECTION
!                        FORALL ( IWL = 1:NWL )
!                           ETOT_SFC_WL ( COL,ROW,IWL ) = MSCALE * ETOT_SFC_WL( COL,ROW,IWL )
!     &                                                 + ACM_CLOUDS( COL,ROW ) * IRRADIANCE( 1,IWL )
!                           TAU_TOT_WL  ( COL,ROW,IWL ) = MSCALE * TAU_TOT_WL( COL,ROW,IWL )
!     &                                                 + ACM_CLOUDS( COL,ROW ) * TAU_TOT( 1,IWL )
!                           TAU_CLOUD_WL( COL,ROW,IWL ) = MSCALE * TAU_CLOUD_WL( COL,ROW,IWL )
!     &                                                 + ACM_CLOUDS( COL,ROW ) * TAU_CLOUD( 1,IWL )
!#ifdef phot_debug
!                           SSA_CLOUD_WL( COL,ROW,IWL ) = MSCALE * SSA_CLOUD_WL( COL,ROW,IWL )
!     &                                                 + ACM_CLOUDS( COL,ROW ) * AVE_SSA_CLD ( IWL )
!                           ASY_CLOUD_WL( COL,ROW,IWL ) = MSCALE * ASY_CLOUD_WL( COL,ROW,IWL )
!     &                                                 + ACM_CLOUDS( COL,ROW ) * AVE_ASYMM_CLD( IWL )
!#endif
!                        END FORALL   ! iwl
!                        FORALL ( LEV = 1:NLAYS, IWL = 1:NWL )
!                           ACTINIC_FX( COL,ROW,LEV,IWL ) = MSCALE * ACTINIC_FX( COL,ROW,LEV,IWL )
!     &                                                   + ACM_CLOUDS( COL,ROW ) * ACTINIC_FLUX( LEV,IWL )
!                        END FORALL ! lev and iwl
!
!                        DO L = 1, N_DIAG_WVL
!                           IWL = DIAG_WVL( L )
!                           FORALL ( LEV = 1:NLAYS)
!                              TAU( COL,ROW,LEV,L ) = MSCALE * TAU( COL,ROW,LEV,L )
!     &                                             + ACM_CLOUDS( COL,ROW ) * TAU_TOT( LEV,IWL )
!                           END FORALL
!                        END DO
!                     END IF  ! photdiag
!!Photolysis rates become a weighted average of the values from resolved and ACM skies
!                     FORALL ( L = 1:NLAYS, IPHOT = 1:NPHOTAB )
!                         RJ( COL,ROW, L, IPHOT ) = 60.0 * ACM_CLOUDS( COL,ROW ) * BLKRJ_ACM( L,IPHOT )
!     &                                           + MSCALE * RJ( COL,ROW,L,IPHOT )
!                     END FORALL  ! Loop on layers and PHOT
!                  END IF
!               END IF ! not USE_ACM_CLOUD and  ACM_CLOUDS > 0

!               IF ( JTIME_CHK ) THEN ! compute clear sky reflection and transmission coefficients
!                  IF ( ANY( CLOUDS ) ) THEN
!                     IF ( CFRAC_2D( COL,ROW ) .GT. 0.0 ) THEN ! resolved and subgrid clouds exist
!                        CLOUDS         = .FALSE.
!                        NEW_PROFILE    = .FALSE.
!                        ONLY_SOLVE_RAD = .TRUE.
!                        CALL NEW_OPTICS ( JDATE, JTIME, NLAYS, 
!     &                                    BLKTA, BLKPRS, BLKDENS, BLKZH, BLKZF,
!     &                                    BLKO3, BLKNO2,
!     &                                    ZSFC, COSZEN, SINZEN, RSQD,
!     &                                    NEW_PROFILE, CLOUDS, CLDFRAC,
!     &                                    BLKRJ_RES, TAUC_AERO, TAU_TOT, TAUO3_TOP,
!     &                                    TAU_RAY, SSA, TAU_CLOUD, TOTAL_O3_COLUMN)
!                        CLR_REFLECTION  ( COL,ROW ) = REFLECTION
!                        CLR_TRANSMISSION( COL,ROW ) = TRANSMISSION
!                        CLR_TRANS_DIRECT( COL,ROW ) = TRANS_DIRECT
!                     ELSE ! only subgrid clouds exist
!                        CLR_REFLECTION  ( COL,ROW ) = RES_SKY_REFLECT
!                        CLR_TRANSMISSION( COL,ROW ) = RES_SKY_TRANS
!                        CLR_TRANS_DIRECT( COL,ROW ) = RES_SKY_TRANSD
!                     END IF
!                  ELSE ! no cloud in vertical column
!                     CLR_REFLECTION  ( COL,ROW ) = REFLECTION
!                     CLR_TRANSMISSION( COL,ROW ) = TRANSMISSION
!                     CLR_TRANS_DIRECT( COL,ROW ) = TRANS_DIRECT
!                  END IF
!               END IF

!            END DO
!         END DO

!      END IF

!...write diagnostic data to output file at the end of every output tstep
!
!      IF ( JTIME_CHK ) THEN
!
!         VARNM = 'COSZENS'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME,
!     &                     COSINE_ZENITH ) ) THEN
!             XMSG = 'Error writing variable ' // VARNM
!             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'OZONE_COLUMN'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, TOC ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'TRANS_DIFFUSE'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, TRANSMIS_DIFFUSE ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'TRANS_DIRECT'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, TRANSMIS_DIRECT ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'REFLECTION'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, REFLECT_COEFF ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'CLR_TRANS_DIF'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, CLR_TRANSMISSION ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'CLR_TRANS_DIR'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, CLR_TRANS_DIRECT ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'CLR_REFLECTION'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, CLR_REFLECTION ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'JNO2'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, RJ( :,:,1, LNO2 ) ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'JO3O1D'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, RJ( :,:,1,LO3O1D ) ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'RESOLVED_CFRAC'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, CFRAC_2D ) ) THEN
!            XMSG = 'Error writing variable ' // VARNM
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         VARNM = 'RESOLVED_WBAR'
!         IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, AVE_HYDROMETEORS ) ) THEN
!             XMSG = 'Error writing variable ' // VARNM
!             CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!
!         IF ( USE_ACM_CLOUD ) THEN
!             VARNM = 'SUBGRID_CFRAC'
!             IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, ACM_CLOUDS ) ) THEN
!                XMSG = 'Error writing variable ' // VARNM
!                CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!             END IF
!             VARNM = 'SUBGRID_WBAR'
!             IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME, ACM_AVE_H2O ) ) THEN
!                XMSG = 'Error writing variable ' // VARNM
!                CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!             END IF
!         END IF
!
!         DO IWL = 1, NWL
!
!            VARNM = 'ETOT_SFC_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, ETOT_SFC_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, ODATE, OTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'TAU_AERO_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, TAU_AERO_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'TAU_CLOUD_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, TAU_CLOUD_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!#ifdef phot_debug
!            VARNM = 'SSA_CLOUD_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, SSA_CLOUD_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'ASY_CLOUD_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, ASY_CLOUD_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!#endif
!
!            VARNM = 'TAU_TOT_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, TAU_TOT_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'TAUO3_TOP_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE,
!     &                         OTIME, TAUO3_TOP_WL( :,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'ALBEDO_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_1, VARNM, ODATE, OTIME,
!     &                         SURFACE_ALBEDO( IWL,:,: ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!         END DO  ! iwl
!
!         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
!     &          'RJ Values written to', CTM_RJ_1,
!     &          'for date and time', ODATE, OTIME
!
!         DO IPHOT = 1, NPHOTAB
!            IF ( .NOT. WRITE3( CTM_RJ_2, PHOTAB( IPHOT ), ODATE,
!     &                         OTIME, RJ( :,:,:,IPHOT ) ) ) THEN
!               XMSG = 'Could not write ' // CTM_RJ_2 // ' file'
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!         END DO
!
!         VARNM = 'CFRAC_3D'
!         IF ( .NOT. WRITE3( CTM_RJ_2, VARNM, ODATE, OTIME, CFRAC_3D ) ) THEN
!              XMSG = 'Could not write ' // TRIM( VARNM ) // ' to ' // CTM_RJ_2 // ' file'
!              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
!         
!         DO IWL = 1, NWL
!            VARNM = 'ACTINIC_FX_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_2, VARNM, ODATE, OTIME, ACTINIC_FX( :,:,:,IWL ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!         END DO
!
!         DO L = 1, N_DIAG_WVL
!            IWL = DIAG_WVL( L )
!
!            VARNM = 'AERO_SSA_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_2, VARNM, ODATE, OTIME, AERO_SSA( :,:,:,L ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'AERO_ASYM_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_2, VARNM, ODATE, OTIME, AERO_ASYM( :,:,:,L ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'TAU_AERO_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_2, VARNM, ODATE, OTIME, TAU_AERO( :,:,:,L ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!
!            VARNM = 'TAU_W' // WLTXT( IWL )
!            IF ( .NOT. WRITE3( CTM_RJ_2, VARNM, ODATE, OTIME, TAU( :,:,:,L ) ) ) THEN
!               XMSG = 'Error writing variable ' // VARNM
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!         END DO
!
!         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
!     &          'RJ and Optical Data written to', CTM_RJ_2,
!     &          'for date and time', ODATE, OTIME
!
!      END IF   ! if JTIME_CHK 

       contains

       subroutine lookup_cmaq_cld_optics( l, r_ef, watp, krad )

       integer, intent(in) :: l     ! lookup table category
       real,    intent(in) :: r_ef  ! effective radius (um)
       real,    intent(in) :: watp  ! ice or liquid water path (g/m^2)
       integer, intent(in) :: krad

       integer :: num, index
       real    :: rstart, rend, rscale, reff
       real    :: fint0, fint1
       real    :: tau, ssa, asm

       num    = cloud_props(l)%num
       rstart = cloud_props(l)%start
       rend   = cloud_props(l)%end

       reff   = max(rstart, min(rend, r_ef))

       rscale = (reff - rstart) / cloud_props(l)%delr
       index  = max(1, min(num-1, int(rscale) + 1))
       fint1  = rscale - real(index-1)
       fint0  = 1.0 - fint1

       do ib = 1, nwl
          tau = ( fint0 * cloud_props(l)%extsw(ib, index  )
     &        + fint1 * cloud_props(l)%extsw(ib, index+1) ) * watp

          ssa = fint0 * cloud_props(l)%ssasw(ib, index  )
     &        + fint1 * cloud_props(l)%ssasw(ib, index+1)

          asm = fint0 * cloud_props(l)%asysw(ib, index  )
     &        + fint1 * cloud_props(l)%asysw(ib, index+1)

          tauclds(ib,krad) = tauclds(ib,krad) + tau
          ssaclds(ib,krad) = ssaclds(ib,krad) + tau * ssa
          asmclds(ib,krad) = asmclds(ib,krad) + tau * ssa * asm
       enddo

       end subroutine lookup_cmaq_cld_optics

      END SUBROUTINE PHOT1
