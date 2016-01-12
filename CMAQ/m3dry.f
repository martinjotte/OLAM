
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

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/yoj/arc/CCTM/src/depv/m3dry/m3dry.F,v 1.12 2012/01/19 14:19:43 yoj Exp $

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE m3dry ( mrl, abflux, sfc_hono )

C-------------------------------------------------------------------------------
C Name:     Models-3 Dry Deposition
C Purpose:  Computes dry deposition velocities using Rst and Ra, and
C           elements of ADOM DD model.
C Revised:  21 Jan 1998  Original version.  (J. Pleim and A. Bourgeois)
C           18 Sep 2001  Made general for USGS 24-category system.
C                        (T. Otte, J. Pleim, and W. Hutzell)
C           14 Jan 2002  Added temperature dependence to Henry's Law
C                        constants.  Added temperature and pressure
C                        dependence to diffusivity.  Added new dry
C                        deposition species, methanol.  (Y. Wu and T. Otte)
C           18 Jan 2002  Changed the reference wet cuticle resistance.
C                        (J. Pleim)
C           09 Jun 2003  Added logic for modeling snow covered surfaces.
C                        Changed the reactivities for SO2, HNO3 and NH3.
C                        Changed pH values to have an east-west variation.
C                        Using the Henry's law constant function from CMAQ in
C                        place of local code.  Also changed the code for
C                        deposition to water to use a pH of 8.1 and the
C                        temperature of water in calculating the Henry's law
C                        constant.  Adjusted values of RSNOW0 = 1000 and
C                        A(NH3) = 20.  Added new dry deposition species: N2O5,
C                        NO3, Generic_aldehyde.  Corrected diffusivities of
C                        chemicals and water and viscosity of air to all be at
C                        the same temperature (273.15K).  Temperature and
C                        pressure adjustments to the values are not needed
C                        because the diffusivities and viscosity are always used
C                        as ratios, so the temperature-pressure dependence was
C                        removed.  Removed dry deposition species, ATRA and
C                        ATRAP, from output.  (D. Schwede, J. Pleim, and
C                        T. Otte)
C           28 Feb 2005  Added optional dry deposition species for chlorine
C                        and mercury.  (G. Sarwar, R. Bullock, and T. Otte)
C           02 Feb 2006  Added mesophyll resistance to dry deposition velocity
C                        calculation, and defined non-zero value for mercury.
C                        (D. Schwede, J. Pleim, and R. Bullock)
C           01 Aug 2007  Added a non-zero mesophyll resistance for NO, NO2, and
C                        CO.  Restored wet cuticle resistance for O3 based on
C                        field study measurements.  Added wet ground resistance.
C                        Changed ground resistance to include partitioning of
C                        wet and dry ground.  Updated pH of rain water for
C                        eastern United States and outside of North America.
C                        Changed reactivity for PAN.  Removed dry deposition
C                        velocity calculations for obsolete chlorine species
C                        ICL1 and ICL2.  Corrected error in the calculation of
C                        surface resistance over water where (Sc/Pr)**(2/3) had
C                        been inadvertently omitted from the numerator.
C                        Surface resistance over water is now a function of
C                        species.  Surface resistance over water now uses wet
C                        bulb temperature rather than ground (water) temperature
C                        in the calculation of the effective Henry's law
C                        constant, and the algorithm has been updated.  Changed
C                        (Sc/Pr)**(2/3) over water to a species-dependent,
C                        meteorologically dependent variable.  Effective Henry's
C                        law constant over land now uses 2-m temperature rather
C                        than layer 1 temperature. Changed ES
C                        into ES_AIR and ES_GRND, and changed QSS into QSS_AIR
C                        and QSS_GRND to clarify usage.  (J. Pleim, E. Cooter,
C                        J. Bash, T. Otte, and G. Sarwar)
C           07 Dec 2007  Add into CMAQ for in-line deposition velocities.
C                        (W. Hutzell, J. Young and T. Otte)
C           07 Jan 2008  Changed the value of d3, the scaling parameter used to
C                        estimate the friction velocity in surface waters from
C                        the atmospheric friction velocity to a value following
C                        Slinn et al. (1978) and Fairall et al. (2007).
C                        (J. Bash)
C           01 Feb 2008  Added bidirectional NH3 flux calculations. (J. Pleim and
C                        J. Young)
C           20 Mar 2008  Added a trap for undefined dry deposition velocities
C                        (e.g., NaN's).  (T. Otte)
C           21 Mar 2008  Added heterogeneous reaction for HONO. It affects HONO, NO2
C                        and HNO3 (G. Sarwar)
C           30 Apr 2008  Added five air toxic species to output.  (W. Hutzell
C                        and T. Otte)
C           August 2008  Applied a minimum value to ustg (0.001 m/s) to prevent 
C                        negative rbg values in the bidi calculation (J. Pleim) 
C           05 Oct 2009  Added condition that vegetation fraction must be
C                        greater than zero to be considered a land point.  This
C                        works around intermittent inconsistencies in surface
C                        fields in some WRF data sets.  (T. Otte)
C           Dec 2009     Revised bidirectional NH3 flux calculations to use soil 
C                        Gamma values read from gridded file.  Bidi flux calcs
C                        based on comparisons with Lillington, NC corn data 2007 (J. Pleim)
C           30 Mar 2010  Modified to output the NH3 bidi stomtal, cuticle and soil component 
C                        fluxes and chaged NH3 bidi variables used in estimating the compensation
C                        point double. (J. Bash)
C           16 Feb 2011 S.Roselle: replaced I/O API include files with UTILIO_DEFN
C           20 May 2011  D.Schwede: add MOSAIC processing
C           14 Jul 2011  Replaced dw25 calculation with Hayduk and Laudie method.
C                        LeBas molar volumes are from the Schroeder additive method
C                        with the exception of HGIIGAS (modeled as HgCl2) which was 
C                        obtained using the Tyn and Calus method. Also, ICL1 and ICL2
C                        were removed. (D. Schwede)
C           27 Jul 2011 J.Bash: Parmaterized the mesophyll resistance as a function
C                        of solubility following Wesely 1989 Atmos Environ.
C           15 Aug 2011  Modified HONO calculation so that deposition velocity for NO2
C                        that is output in DEPV file does not include the loss due to
C                        the heterogeneous reaction. This additional loss is now
C                        accounted for in vdiff.F (D. Schwede and G. Sarwar)
C           29 Aug 2011 Added NH3 bidirectional flux variables and modules and integrated
C                       NH3 bidi algorithms with MOSAIC algorithms. NH3 bidi routines now
C                       read in foratted EPIC output and maintain a soil NH4 budget and 
C                       fluxes are calculated for individual land cover types. 
C                       (J. Bash and D. Schwede)
C           22 Sep 2011 -- incorporated twoway model implementation
C                       -- removed non-use dluse array
C                          (David Wong)
C           26 Sep 2011 -- made the number of actual and dummy arguments the same in
C                          calling subroutine Init_ABFlux
C                          (David Wong)
C-------------------------------------------------------------------------------

      use depvvars
      use depv_defn,   only: depvel_gas_land, depvel_gas_sea, ie_hono, ic_no2
      use utilio_defn
!     Use ABFlux_Mod
!     Use Bidi_Mod,    only: lufrac
      use cgrid_defn,  only: cgrid, vdemis_gc
      use misc_coms,   only: io6, isubdomain
      use mem_ijtabs,  only: itabg_w
      use leaf_coms,   only: mwl
      use mem_leaf,    only: land, itab_wl
      use sea_coms,    only: mws
      use mem_sea,     only: sea,  itab_ws
      use mem_grid,    only: dzt, volti
      use mem_basic,   only: sh_v, theta, rho
      use mem_turb,    only: pblh
      use leaf_coms,   only: nzg, slcpd
      use hlconst_mod, only: hlconst, hlconst_spcs_init

      IMPLICIT NONE

C Includes:

      INCLUDE 'CONST.EXT'     ! constants
!     INCLUDE 'FILES_CTM.EXT' ! file name parameters

C Arguments:

      INTEGER,     INTENT( IN )  :: mrl
      LOGICAL,     INTENT( IN )  :: abflux
      LOGICAL,     INTENT( IN )  :: sfc_hono
      
C Local Variables:

      REAL,            PARAMETER :: a0         = 8.0     ! [dim'less]
      INTEGER                    :: c
      REAL,            PARAMETER :: d3         = 1.38564e-2 ! [dim'less]
                                               ! k*sqrt(rhoair/rhowater) from Slinn 78
      real                       :: delta
      real                       :: d_ice      ! dep. vel. over ice
      REAL                       :: dw
      REAL                       :: dw25       ! diffusivity of water at 298.15 k
      REAL,            PARAMETER :: dwat       = 0.2178  ! [cm^2/s] at 273.15K
      LOGICAL                    :: effective  ! true=compute effective Henry's Law const

      LOGICAL, SAVE              :: first_call = .TRUE.
      real                       :: glat
      real                       :: glon
      REAL                       :: heff                 ! effective Henry's Law constant
      REAL                       :: heff_ap              ! Henry's Law constant for leaf apoplast M/atm
!     REAL,            EXTERNAL  :: hlconst              ! [M / atm]
      REAL                       :: hplus
      REAL,            PARAMETER :: hplus_ap   = 1.0e-6  ! pH=6.0 leaf apoplast solution Ph (Massad et al 2008)
      REAL,            PARAMETER :: hplus_def  = 1.0e-5  ! pH=5.0
      REAL,            PARAMETER :: hplus_east = 1.0e-5  ! pH=5.0
      REAL,            PARAMETER :: hplus_h2o  = 7.94328e-9 ! 10.0**(-8.1)
      REAL,            PARAMETER :: hplus_west = 3.16228e-6 ! 10.0**(-5.5)
      real                       :: hveg
      LOGICAL                    :: ifurban
      LOGICAL                    :: ifsnow
      LOGICAL                    :: isocean
      integer                    :: iw
      integer                    :: iwl
      integer                    :: iws
      integer                    :: j
      REAL,            PARAMETER :: kvis       = 0.132   ! [cm^2 / s] at 273.15K
      REAL                       :: kvisw      ! kinematic viscosity of water [cm^2/s]
      INTEGER                    :: kw
      INTEGER                    :: l
      REAL                       :: lai                ! cell leaf area index
      INTEGER                    :: n
      INTEGER                    :: nlev_water
      CHARACTER( 16 ), PARAMETER :: pname      = 'M3DRY'
      REAL,            PARAMETER :: pr         = 0.709   ! [dim'less]
      real                       :: pvd
      INTEGER                    :: r
      real                       :: ra                   ! cell aerodynamic resistance
      real                       :: raice                ! aerodynamic resistance over ice
      real                       :: raw                  ! aerodynamic resistance over water
      REAL                       :: rac
      REAL                       :: rbc
      REAL                       :: rbsulf
      REAL                       :: rbw
      REAL                       :: rci
      REAL                       :: rcut
      REAL,            PARAMETER :: rcut0      = 3000.0  ! [s/m]
      REAL,            PARAMETER :: rcw0       = 125000.0 ! acc'd'g to Padro and
                                                          ! adapted from Slinn 78
      REAL,            PARAMETER :: rg0        = 1000.0  ! [s/m]
      REAL                       :: rgnd
      REAL                       :: rgndc
      REAL                       :: rgw                  ! resist for water-covered sfc
      REAL,            PARAMETER :: rgwet0     = 25000.0 ! [s/m]
      REAL                       :: rh_air               ! rel humidity (air)
      REAL                       :: rh_grnd              ! rel humidity (ground)
      REAL,            external  :: rhovsil
      REAL                       :: rinc
      REAL                       :: rs                   ! stomatal resistance [s/m]
      REAL,            PARAMETER :: rsndiff    = 10.0    ! snow diffusivity fac
      REAL                       :: rsnow
      REAL,            PARAMETER :: rsnow0     = 1000.0
      REAL                       :: rstomi
      REAL                       :: rsurf
      REAL                       :: rwet                 ! wet sfc resist (cuticle or grnd)
      REAL                       :: rwetsfc
      REAL                       :: scw_pr_23            ! (scw/pr)**2/3
      real                       :: snowfac
      real                       :: sst
      REAL,            PARAMETER :: svp2       = 17.67   ! from MM5 and WRF
      REAL,            PARAMETER :: svp3       = 29.65   ! from MM5 and WRF
      REAL,            PARAMETER :: rt25inK    = 1.0/(stdtemp + 25.0) ! 298.15K = 25C
      REAL                       :: tempcr               ! cell canopy temp
      REAL                       :: tempgr               ! cell ground temp
      REAL                       :: tice
      REAL                       :: tw
      REAL,            PARAMETER :: twothirds  = 2.0 / 3.0
      REAL                       :: ustar                ! cell friction velocity
      REAL                       :: ustari               ! friction velocity over ice
      REAL                       :: ustarw               ! friction velocity over water
      REAL                       :: vegfr                ! cell veg coverage fraction
      REAL                       :: vegfrac              ! cell veg coverage fraction including snow
      CHARACTER( 16 )            :: vname
      real                       :: wrcr
      REAL                       :: wrmax
      REAL                       :: wstar
      REAL                       :: wtv0
      REAL                       :: xm                   ! liquid water mass frac
      REAL                       :: xt                   ! liquid water mass frac
      CHARACTER( 96 )            :: xmsg = ' '
      CHARACTER( 16 ), SAVE      :: vname_ra, vname_rc, vname_rn, vname_rs

      INTEGER, SAVE              :: n_spc_m3dry = ltotg       ! from DEPVVARS module

      REAL                       :: ar       ( ltotg )        ! reactivity relative to HNO3
      REAL                       :: dif0     ( ltotg )        ! molecular diffusivity [cm2/s]
      REAL                       :: lebas    ( ltotg )        ! Le Bas molar volume [cm3/mol ]
      REAL                       :: meso     ( ltotg )        ! Exception for species that 
                                                              ! react with cell walls. fo in 
                                                              ! Wesely 1989 eq 6.
      REAL, SAVE                 :: scc_pr_23( ltotg )        ! (SCC/PR)**2/3, fn of DIF0

      CHARACTER( 16 )            :: subname  ( ltotg )        ! for subroutine HLCONST
      integer, save              :: hlspc    ( ltotg )        ! species index in hlconst
      INTEGER                    :: SPC

C-------------------------------------------------------------------------------
C For ozone exchange over the ocean
      REAL                       :: pchang      ! p value used in equation (12) in Chang et al( 2004)
      REAL                       :: kwchang     ! kw value used in Chang et al (2004)
      REAL                       :: ciodide     ! iodide concentration (umol/L)
      REAL                       :: qiodide     ! q in Chang et al (2004)

C-------------------------------------------------------------------------------
C For heterogenous hono on surfaces
      REAL                       :: kno2                      ! first order rate constant for the heterogenous reaction [1/s]  
      REAL                       :: surf_bldg                 ! Surface area of bldgs to volume of air [-]
      REAL                       :: surf_leaf                 ! Surface area of leaves to volume of air [-]
      REAL                       :: zf                        ! layer height [m]      
      
C-------------------------------------------------------------------------------
C For Ammonia bi-directional flux                      
!      REAL                       :: sltyp     ( ncols,nrows ) ! soil type
!      REAL                       :: wg        ( ncols,nrows ) ! soil moisture in top 1 cm (vol frc)
!      REAL                       :: w2        ( ncols,nrows ) ! soil moisture in top 1 m (vol frc)
!      REAL                       :: cnh3         ! Layer 1 NH3 concentration from CGRID [ppm]
!      REAL                       :: tpvd         ! temp pvd variable
!      REAL                       :: f_stom       ! stomatal flux component
!      REAL                       :: f_cut        ! cuticular flux component
!      REAL                       :: f_soil       ! soil flux component
!      REAL                       :: f_emis       ! Emission flux component
!      REAL                       :: f_dep        ! Deposition flux component
!      REAL                       :: f_ag         ! flux from agriculture
!      REAL                       :: f_nat        ! flux from natural systems
!      REAL                       :: f_wat        ! direct flux to surface water
!      REAL                       :: lnh3         ! loss part of ammonia flux
                                                  ! Note that lnh3 is stored in DEPVEL_GAS(nh3) for use in vdiff

C-------------------------------------------------------------------------------
C-- Chemical-Dependent Parameters (Original Source: Modified ADOM - Padro)
C
C                                                        at 298.15 K
C     Species   Dif0(cm2/s) Alphastar   Reactivity     -DKHOR [K]  KH0 [M/atm]
C     _______   ___________ _________   __________     _________  __________
C  1   SO2      0.1089 ~    1000.       8.00            3100.0 ~     1.2e00 ~
C  2   SULF      --       --           --                --           --
C  3   NO2      0.1361 ~    1.00      8.00 2.00*        2500.0 ~     1.2e-2 ~
C  4   NO       0.1802 ~    1         2.0+              1500.0 ~     1.9e-3 ~
C  5   O3       0.1444 ~    10.00       15.0  8@        2700.0 ~     1.2e-2 ~
C  6   HNO3     0.1628 ~    1E9       18.0 800* 8000**  8700.0 ~     2.6e+6 ~
C  7   H2O2     0.2402 ~    1.00      12.0 30*          6600.0 ~     9.7e+4 ~
C  8   ACET ALD 0.1525 ~    1         10.+              5700.0 ~     1.3e+1 ~
C  9   HCHO     0.1877 ~    1.00      10.+              6400.0 ~     7.0e+3 ~
C 10   OP       0.1525 ~    1         10.0+             5200.0 ~     3.1e+2 ~
C 11   PAA      0.1220 ~    1         20+               5300.0 ~     8.4e+2 ~
C 12   ORA      0.1525 ~    1         20+               5700.0 ~     3.7e+3 ~
C 13   NH3      0.1978 ~    1E5       10.0              4100.0 ~     5.8e+1 ~
C 14   PAN      0.0938 ~    1         16.0~~            5900.0 ~     2.9e00 ~
C 15   HONO     0.1525 ~    1         20+               4800.0 ~     4.9e+1 ~
C 16   CO       0.1807 ~    1         5.0               1600.0 ~     9.5e-4 ~
C 17   METHANOL 0.1363 ~    1.0 ~       2.0 ~           5200.0 ~     2.2e+2 ~
C --   CO2      0.1381 ~                                2400.0 ~     3.4e-2 ~
C 18   N2O5     0.0808 ^^             5000.0**
C 19   NO3      0.1153 ^^             5000.0**
C 20   GEN ALD  0.1525 ##   1.0       10.0##
C 21   CL2      0.1080 %              10.0 %
C 22   HOCL     0.1300 %              10.0 %
C 23   HCL      0.1510 %              8000.0 %
C 24   FMCL     0.1094 %              10.0 %
C 27   HG       0.1194 $              0.1 $
C 28   HGIIGAS  0.0976 $              8000.0 $
C 50   HEXAMETHYLENE DIISOCYANATE       10.0 <>
C 51   HYDRAZINE                        20.0 <>
C 52   MALEIC ANHYDRIDE                 10.0 <>
C 53   TOLUENE DIISOCYANATE             10.0 <>
C 54   TRIETHYLAMINE                    10.0 <>
C
C---------Notes
C  * Updates based on literature review 7/96 JEP
C  # Diff and H based on Wesely (1988) same as RADM
C  + Estimated by JEP 2/97
C  @ Updated by JEP 9/01
C  ~ Added by YW 1/02.  Dif0 based on Massman (1998).  Henry's Law constant
C    is defined here as: h=cg/ca, where cg is the concentration of a species
C    in gas-phase, and ca is its aqueous-phase concentration.  The smaller h,
C    the larger solubility.  Henry's Law constant in another definition (KH):
C    KH = ca/pg [M/atm], KH = KH0 * exp(-DKH/R(1/T-1/T0)), where KH0 and -DKH
C    values are from Rolf Sander (1999).  h=1/(KH*R*T).
C ** Update by DBS based on estimates by JEP 1/03
C ^^ From Bill Massman, personal communication 4/03
C ## Diffusivity calculated by SPARC, reactivity = other aldehydes
C ++ Dif0 in Massman is diffusivity at temperature 0C and 1 atm (101.325kPa), so
C    chemicals that were not in Massman's paper need to be adjusted.  We assume
C    JEP's original values were for 25C and 1 atm.
C  % Added by G. Sarwar (10/04)
C  $ Added by R. Bullock (02/05) HG diffusivity is from Massman (1999).
C    HGIIGAS diffusivity calculated from the HG value and a mol. wt. scaling
C    factor of MW**(-2/3) from EPA/600/3-87/015. ORD, Athens, GA.  HGIIGAS
C    mol.wt. used is that of HgCl2.  Reactivity of HG is 1/20th of NO and NO2
C    values based on general atmospheric lifetimes of each species.  Reactivity
C    of HGIIGAS is based on HNO3 surrogate.
C @@ Mesophyll resistances for NO, NO2, and CO added by J. Pleim (07/07) based
C    on values in Pleim, Venkatram, and Yamartino, 1984:  ADOM/TADAP Model
C    Development Program, Volume 4, The Dry Deposition Module.  ERT, Inc.,
C    Concord, MA (peer reviewed).
C ~~ Reactivity for PAN changed from 4.0 to 16.0 by J. Pleim (07/07) based on
C    comparisons with Turnipseed et al., JGR, 2006.
C %% Species ICL1 and ICL2 are removed, not used in CB05.  G. Sarwar (07/07)
C <> Hazardous Air Pollutants that are believed to undergo significant dry
C    deposition. Hydrazine and triethylamine reactivities are based on analogies
C    to NH3. Maleic anhydride reactivity is assumed similar to aldehydes.
C    Toluene diisocyanate and hexamethylene diisocyanate reactivities are
C    assumed to be similar to SO2. Diffusivities are calculated with standard
C    formulas.  W. Hutzell (04/08)
C-------------------------------------------------------------------------------

      DATA subname( 1), dif0( 1), ar( 1), meso( 1), lebas( 1) / 'SO2             ', 0.1089,   10.0,      0.0,  35.0/
      DATA subname( 2), dif0( 2), ar( 2), meso( 2), lebas( 2) / 'SULFATE         ', 0.0001,    0.0,      0.0,  49.0/
      DATA subname( 3), dif0( 3), ar( 3), meso( 3), lebas( 3) / 'NO2             ', 0.1361,    2.0,      0.1,  21.0/
      DATA subname( 4), dif0( 4), ar( 4), meso( 4), lebas( 4) / 'NO              ', 0.1802,    2.0,      0.0,  14.0/
      DATA subname( 5), dif0( 5), ar( 5), meso( 5), lebas( 5) / 'O3              ', 0.1444,    8.0,      1.0,  21.0/
      DATA subname( 6), dif0( 6), ar( 6), meso( 6), lebas( 6) / 'HNO3            ', 0.1067, 8000.0,      0.0,  35.0/
      DATA subname( 7), dif0( 7), ar( 7), meso( 7), lebas( 7) / 'H2O2            ', 0.1300,   30.0,      1.0,  28.0/
      DATA subname( 8), dif0( 8), ar( 8), meso( 8), lebas( 8) / 'ACETALDEHYDE    ', 0.1111,   10.0,      0.0,  56.0/
      DATA subname( 9), dif0( 9), ar( 9), meso( 9), lebas( 9) / 'FORMALDEHYDE    ', 0.1554,   10.0,      0.0,  35.0/
      DATA subname(10), dif0(10), ar(10), meso(10), lebas(10) / 'METHYLHYDROPEROX', 0.1179,   10.0,      0.1,  49.0/
      DATA subname(11), dif0(11), ar(11), meso(11), lebas(11) / 'PEROXYACETIC_ACI', 0.0868,   20.0,      0.1,  70.0/
      DATA subname(12), dif0(12), ar(12), meso(12), lebas(12) / 'ACETIC_ACID     ', 0.0944,   20.0,      0.0,  63.0/
      DATA subname(13), dif0(13), ar(13), meso(13), lebas(13) / 'NH3             ', 0.1978,   20.0,      0.0,  28.0/
      DATA subname(14), dif0(14), ar(14), meso(14), lebas(14) / 'PAN             ', 0.0687,   16.0,      0.1,  91.0/
      DATA subname(15), dif0(15), ar(15), meso(15), lebas(15) / 'HNO2            ', 0.1349,   20.0,      0.1,  28.0/
      DATA subname(16), dif0(16), ar(16), meso(16), lebas(16) / 'CO              ', 0.1807,    5.0,      0.0,  14.0/
      DATA subname(17), dif0(17), ar(17), meso(17), lebas(17) / 'METHANOL        ', 0.1329,    2.0,      0.0,  42.0/
      DATA subname(18), dif0(18), ar(18), meso(18), lebas(18) / 'N2O5            ', 0.0808, 5000.0,      0.0,  49.0/
      DATA subname(19), dif0(19), ar(19), meso(19), lebas(19) / 'NO3             ', 0.1153, 5000.0,      0.0,  28.0/
      DATA subname(20), dif0(20), ar(20), meso(20), lebas(20) / 'GENERIC_ALDEHYDE', 0.0916,   10.0,      0.0,  56.0/
      DATA subname(21), dif0(21), ar(21), meso(21), lebas(21) / 'CL2             ', 0.1080,   10.0,      0.0,  49.0/
      DATA subname(22), dif0(22), ar(22), meso(22), lebas(22) / 'HOCL            ', 0.1300,   10.0,      0.0,  38.5/
      DATA subname(23), dif0(23), ar(23), meso(23), lebas(23) / 'HCL             ', 0.1510, 8000.0,      0.0,  31.5/
      DATA subname(24), dif0(24), ar(24), meso(24), lebas(24) / 'FMCL            ', 0.1094,   10.0,      0.0,  45.5/
      DATA subname(25), dif0(25), ar(25), meso(25), lebas(25) / 'HG              ', 0.1194,    0.1,      0.0,  14.8/   ! lebas not used
      DATA subname(26), dif0(26), ar(26), meso(26), lebas(26) / 'HGIIGAS         ', 0.0976, 8000.0,      0.0,  95.0/   ! estimation from back calculating to get dw25 = 1.04e-5 (Garland et al, 1965)
      DATA subname(27), dif0(27), ar(27), meso(27), lebas(27) / 'TECDD_2378      ', 0.0525,    2.0,      0.0, 217.0/
      DATA subname(28), dif0(28), ar(28), meso(28), lebas(28) / 'PECDD_12378     ', 0.0508,    2.0,      0.0, 234.5/
      DATA subname(29), dif0(29), ar(29), meso(29), lebas(29) / 'HXCDD_123478    ', 0.0494,    2.0,      0.0, 252.0/
      DATA subname(30), dif0(30), ar(30), meso(30), lebas(30) / 'HXCDD_123678    ', 0.0494,    2.0,      0.0, 252.0/
      DATA subname(31), dif0(31), ar(31), meso(31), lebas(31) / 'HXCDD_123478    ', 0.0494,    2.0,      0.0, 252.0/
      DATA subname(32), dif0(32), ar(32), meso(32), lebas(32) / 'HPCDD_1234678   ', 0.0480,    2.0,      0.0, 269.5/
      DATA subname(33), dif0(33), ar(33), meso(33), lebas(33) / 'OTCDD           ', 0.0474,    2.0,      0.0, 287.0/
      DATA subname(34), dif0(34), ar(34), meso(34), lebas(34) / 'TECDF_2378      ', 0.0534,    2.0,      0.0, 210.0/
      DATA subname(35), dif0(35), ar(35), meso(35), lebas(35) / 'PECDF_12378     ', 0.0517,    2.0,      0.0, 227.5/
      DATA subname(36), dif0(36), ar(36), meso(36), lebas(36) / 'PECDF_23478     ', 0.0517,    2.0,      0.0, 227.5/
      DATA subname(37), dif0(37), ar(37), meso(37), lebas(37) / 'HXCDF_123478    ', 0.0512,    2.0,      0.0, 245.0/
      DATA subname(38), dif0(38), ar(38), meso(38), lebas(38) / 'HXCDF_123678    ', 0.0512,    2.0,      0.0, 245.0/
      DATA subname(39), dif0(39), ar(39), meso(39), lebas(39) / 'HXCDF_234678    ', 0.0512,    2.0,      0.0, 245.0/
      DATA subname(40), dif0(40), ar(40), meso(40), lebas(40) / 'HXCDF_123789    ', 0.0512,    2.0,      0.0, 245.0/
      DATA subname(41), dif0(41), ar(41), meso(41), lebas(41) / 'HPCDF_1234678   ', 0.0487,    2.0,      0.0, 262.5/
      DATA subname(42), dif0(42), ar(42), meso(42), lebas(42) / 'HPCDF_1234789   ', 0.0487,    2.0,      0.0, 262.5/
      DATA subname(43), dif0(43), ar(43), meso(43), lebas(43) / 'OTCDF           ', 0.0474,    2.0,      0.0, 280.0/
      DATA subname(44), dif0(44), ar(44), meso(44), lebas(44) / 'NAPHTHALENE     ', 0.0778,    4.0,      0.0, 119.0/
      DATA subname(45), dif0(45), ar(45), meso(45), lebas(45) / '1NITRONAPHTHALEN', 0.0692,    4.0,      0.0, 133.0/
      DATA subname(46), dif0(46), ar(46), meso(46), lebas(46) / '2NITRONAPHTHALEN', 0.0692,    4.0,      0.0, 133.0/
      DATA subname(47), dif0(47), ar(47), meso(47), lebas(47) / '14NAPHTHOQUINONE', 0.0780,    4.0,      0.0, 119.0/
      DATA subname(48), dif0(48), ar(48), meso(48), lebas(48) / 'HEXAMETHYLE_DIIS', 0.0380,   10.0,      0.0, 196.0/
      DATA subname(49), dif0(49), ar(49), meso(49), lebas(49) / 'HYDRAZINE       ', 0.4164,   20.0,      0.0,  42.0/
      DATA subname(50), dif0(50), ar(50), meso(50), lebas(50) / 'MALEIC_ANHYDRIDE', 0.0950,   10.0,      0.0,  70.0/
      DATA subname(51), dif0(51), ar(51), meso(51), lebas(51) / '24-TOLUENE_DIIS ', 0.0610,   10.0,      0.0, 154.0/
      DATA subname(52), dif0(52), ar(52), meso(52), lebas(52) / 'TRIETHYLAMINE   ', 0.0881,   20.0,      0.0, 154.0/
      DATA subname(53), dif0(53), ar(53), meso(53), lebas(53) / 'ORG_NTR         ', 0.0607,   16.0,      0.1, 160.0/  ! assumes 58.2% C5H11O4N and 41.8% C5H11O3N
      DATA subname(54), dif0(54), ar(54), meso(54), lebas(54) / 'HYDROXY_NITRATES', 0.0609,   16.0,      0.1, 156.1/
      DATA subname(55), dif0(55), ar(55), meso(55), lebas(55) / 'MPAN            ', 0.0580,   16.0,      0.1, 133.0/
      DATA subname(56), dif0(56), ar(56), meso(56), lebas(56) / 'PPN             ', 0.0631,   16.0,      0.1, 118.2/
      DATA subname(57), dif0(57), ar(57), meso(57), lebas(57) / 'MVK             ', 0.0810,    8.0,      1.0,  88.8/
      DATA subname(58), dif0(58), ar(58), meso(58), lebas(58) / 'DINTR           ', 0.0810,    8.0,      0.0,  88.8/
      DATA subname(59), dif0(59), ar(59), meso(59), lebas(59) / 'NTR_ALK         ', 0.0810,    8.0,      0.0,  88.8/
      DATA subname(60), dif0(60), ar(60), meso(60), lebas(60) / 'NTR_OH          ', 0.0810,    8.0,      0.0,  88.8/
      DATA subname(61), dif0(61), ar(61), meso(61), lebas(61) / 'NTR_PX          ', 0.0810,    8.0,      0.0,  88.8/
      DATA subname(62), dif0(62), ar(62), meso(62), lebas(62) / 'PROPNN          ', 0.0810,    8.0,      0.0,  88.8/
      DATA subname(63), dif0(63), ar(63), meso(63), lebas(63) / 'NITRYL_CHLORIDE ', 0.0888,    8.0,      0.0,  45.5/   ! dif0 estimated following Erickson III et al., JGR, 104, D7, 8347-8372, 1999
      DATA subname(64), dif0(64), ar(64), meso(64), lebas(64) / 'ISOPNN          ',0.0457,    8.0,      0.0,  206.8/  
      DATA subname(65), dif0(65), ar(65), meso(65), lebas(65) / 'MTNO3           ',0.0453,    8.0,      0.0,  251.2/  

      IF ( first_call ) THEN
         first_call = .FALSE.

         DO l = 1, n_spc_m3dry
            IF ( dif0( l ) > 0.0 ) THEN
               scc_pr_23( l ) = ( ( kvis / dif0( l ) ) / pr ) ** twothirds
            ELSE
               scc_pr_23( l ) = 0.0
            END IF
         END DO

         ! Set up Henry's law species indices

         hlspc( : ) = 0

         DO l = 1, n_spc_m3dry
            IF ( .NOT. use_depspc( l ) ) CYCLE
            call hlconst_spcs_init ( subname( l ), hlspc( l ) )
         enddo

      END IF   ! first_call
      
C-------------------------------------------------------------------------------
C Loop over land cells and calculate dry deposition.
C-------------------------------------------------------------------------------

      effective = .TRUE.

      do iwl = 2, mwl

         iw = itab_wl(iwl)%iw
         if (isubdomain == 1) then
            iw = itabg_w(iw)%iw_myrank
         endif
         kw = itab_wl(iwl)%kw

         ! If run is parallel, get local rank indices

         glat  = land%glatw(iwl)
         glon  = land%glonw(iwl)
         ustar = land%ustar(iwl)
         zf    = dzt(kw)

         ra    = 1.0 / land%ggaer(iwl)

         wtv0  = land%sfluxt(iwl) * (1. + .61 * sh_v(kw,iw))
     &         + land%sfluxr(iwl) * .61 * theta(kw,iw)

         lai   = land%veg_lai(iwl)
         vegfr = land%veg_fracarea(iwl)
         wrcr  = land%veg_water(iwl) * 1000.0  ! kg/m^2 -> m
         hveg  = land%veg_height(iwl)
         rs    = land%stom_resist(iwl)

         tempcr = land%cantemp(iwl)
         rh_air = 100.0 * land%canshv(iwl) * rho(kw,iw) / rhovsil( land%cantemp(iwl) - 273.15 )
         rh_air = min( 100.0, max( rh_air, 0.0 ) )

         nlev_water = land%nlev_sfcwater(iwl)
         if (nlev_water > 0) then
            snowfac = land%snowfac(iwl)
         else
            snowfac = 0.0
         endif

         ifurban   = ( any(land%leaf_class(iwl) == (/ 19, 21 /)) )

         if (wtv0 > 0.0) then
            wstar  = (grav * pblh(iw) * wtv0 / theta(kw,iw)) ** 0.33333333
         else
            wstar  = 0.0
         endif

         if ( nlev_water > 0 ) then

           ! Determine sfcwater temperature and liquid water mass fraction (0.0 to 0.5)
            call qtk( land%sfcwater_energy(nlev_water,iwl), tempgr, xm )

            if (xm < 0.9) then

               ifsnow = .true.
               xm = MIN (xm, 0.5)
               xm = MAX (xm, 0.0)

            else
               
               ifsnow = .false.
               xm     = 1.0

            endif

         else

            ! Determine the bare ground temperature
            call qwtk( land%soil_energy(nzg,iwl),land%soil_water(nzg,iwl)*1.e3,
     &                 slcpd(land%ntext_soil(nzg,iwl)), tempgr, xm )

            ifsnow = .false.
            xm     = 0.0

         endif

         ! Canopy Wetness

         wrmax = 0.2e-3 * vegfr * lai ! [m]
         IF ( wrcr .LE. 0.0 ) THEN
            delta  = 0.0
         ELSE
            delta = wrcr / wrmax  ! refer to SiB model
            delta = MIN( delta, 1.0 )
         END IF

         ! Assign a pH for rain water based on longitude if US simulation.
         ! Otherwise use default pH.  Use pH value in HPLUS calculation.

         IF ( ( glat .GE.   30.0 ) .AND. ( glat .LE.  45.0 ) .AND.
     &        ( glon .GE. -120.0 ) .AND. ( glon .LE. -70.0 ) ) THEN

            IF ( glon .GT. -100.0 ) THEN
               hplus = hplus_east
            ELSE
               hplus = hplus_west
            ENDIF
         ELSE
            hplus = hplus_def
         ENDIF

         ! Loop over species to calculate dry deposition velocities.

         depvel_gas_land( iwl,: ) = 0.0  ! initialize for this time period

         n = 0
         dloopl: DO l = 1, n_spc_m3dry
         
            IF ( .NOT. use_depspc( l ) ) CYCLE dloopl

            n = n + 1

!           IF ( depspc( l ) .EQ. 'SULF' ) THEN  ! Sulfate (SULF)
            IF ( l .EQ. l_sulf ) THEN  ! Sulfate (SULF)

         ! Sulfate calculation follows Wesely (1985), Eqn. 11.

               rbsulf = 500. * ustar / (ustar**2 + 0.24 * wstar**2)
               depvel_gas_land( iwl,n ) =  1.0 / ( ra + rbsulf )
               
            else

         ! Use CMAQ function for calculating the effective Henry's Law
         ! constant.  Note that original M3DRY wants inverse,
         ! non-dimensional Henry's Law (caq/cg).

               heff  = hlconst( tempcr, effective, hplus, hlspc( l ) )

         ! Make Henry's Law constant non-dimensional.
               
               heff  = heff * 0.08205 * tempcr

         ! Wet surface resistance.  (Note DELTA = CWC in ADOM lingo.)
         ! This now applies to cuticle and ground.

!              IF ( depspc( l ) .NE. 'O3' ) THEN
               IF ( l .NE. l_o3 ) THEN
                  rwet = rcw0 / heff ! wet cuticle
               ELSE
                  ! Set RCW/LAI = 200 s/m on basis of Keysburg exp for O3
                  rwet = 1250.0 ! s/m
               END IF

               rgw  = rgwet0 / heff ! wet ground

         ! Dry snow resistance.

               rsnow = rsnow0 * a0 / ar( l )

         ! If the surface is cold and wet, use dry snow.

               IF ( tempgr .LT. stdtemp ) THEN
                  rwetsfc = rsnow
               ELSE
                  rwetsfc = rwet
               END IF

         ! Dry cuticle resistance.

!              IF ( depspc( l ) .NE. 'NH3' ) THEN
               IF ( l .NE. l_nh3 ) THEN
                  rcut = rcut0 * a0 / ar( l )
               ELSE
                  rcut = 4000.0 * EXP( -0.054 * rh_air )
               END IF

         ! Dry ground resistance.  (revised according to Erisman)

               rinc  = 14.0 * lai * hveg / ustar
               rgnd  = rg0 * a0 / ar( l )
               rgndc = 1.0 / ( ( 1.0 - delta ) / rgnd + delta / rgw )
     &               + rinc          ! Add in-canopy part

         ! Bulk stomatal resistance; include mesophyll resistance.

               heff_ap = hlconst( tempcr, effective, hplus_ap, hlspc( l ) )

!              rstom = rs * dwat / dif0( l )
!    &               + 1.0 / ( heff_ap / 3000.0 + 100.0 * meso( l ) ) / lai

               rstomi = lai / ( lai * rs * dwat / dif0( l ) 
     &                          + 1.0 / ( heff_ap / 3000.0 + 100.0 * meso( l ) ) )

         ! Bulk surface resistance.

!                  rci = vegfr
!     &                * ( 1.0/rstom + (1.0-delta( c,r ) ) * lai / rcut
!     &                +   ( delta( c,r ) * lai / rwetsfc ) + 1.0 / rgndc )
!     &                + real( 1-ifsnow ) * ( (1.0 - vegfr) * ( (1.0-delta( c,r ) ) /
!     &                                  rgnd + delta( c,r ) / rgw ) )
!     &                + real( ifsnow ) * ( (1.0 - xm) / rsnow + xm / (rsndiff + rgw) )

               vegfrac = vegfr * (1.0 - snowfac)

               rci = vegfrac * ( rstomi + (1.0-delta ) * lai / rcut +
     &                    delta * lai / rwetsfc + 1.0 / rgndc )

               if ( ifsnow ) then
                  rci = rci + (1.0 - vegfrac) * ( (1.0-xm) / rsnow + xm / (rsndiff + rgw) )
               else
                  rci = rci + (1.0 - vegfrac) * ( (1.0-delta) / rgnd + delta / rgw )
               endif

               rsurf = 1.0 / rci

         ! Compute dry deposition velocity.

               rbc = 5.0 / ustar * scc_pr_23( l )
               rac = ra + rbc

               depvel_gas_land( iwl,n ) = 1.0 / ( rsurf + rac )

! TODO: SURFACE SOURCE OF AMMONIA
!
!              IF ( abflux ) THEN   ! Ammonia Bidirectional Flux
!
!                  IF ( depspc( l ) .EQ. 'NH3' ) THEN
!                    cnh3  = cgridl1( n,c,r )
!
!                    CALL Get_Flux( tempgcr,rh_air,cnh3,rwetsfc,rgw,wg(c,r),w2(c,r),
!     &                                 sltyp(c,r),dif0(l),r,c,l,tpvd,lnh3,
!     &                                 f_stom,f_cut,f_soil,f_emis,f_dep,f_ag,f_nat, f_wat,
!     &                                 dt(2) )
!                    pvd( n,c,r ) = tpvd
!                    depvel_gas( n,c,r ) = lnh3
!                    cmp(1,c,r) = f_emis
!                    cmp(2,c,r) = f_dep
!                    cmp(3,c,r) = f_stom
!                    cmp(4,c,r) = f_cut
!                    cmp(5,c,r) = f_soil
!                    cmp(6,c,r) = f_ag
!                    cmp(7,c,r) = f_nat
!                    cmp(8,c,r) = f_wat
!                 END IF   ! 'NH3'
!              END IF   ! abflux

               IF ( sfc_hono ) THEN

! HONO production via heterogeneous reaction on ground surfaces,
! 2NO2 = HONO + HNO3
! Rate constant for the reaction = (3.0E-3/60)* (A/V),
! where A/V is surface area/volume ratio
! HONO is produced and released into the atmosphere
! NO2 is lost via chemical reaction
! HNO3 is sticky and stays on the surfaces

! Calculate A/V for leaves.
! LAI was multiplied by 2 to account for the fact that surface area
! is provided by both sides of the leaves.
! Matthews Jones, Ammonia deposition to semi-natural vegetation,
! PhD dissertation, University of Dundee, Scotland, 2006

                  surf_leaf = 2.0 * lai / zf

! Calculate A/V for buildings and other structures.
! Buildings and other structures can provide additional surfaces in
! urban areas for the heterogeneous reaction to occur. However, such
! information is not readily available; in the absence of such information,
! it is scaled to purb(c,r). Svensson et al., (1987) suggests a typical value
! of 0.2 for A/V for buildings in urban environments. A maximum value of 0.2
! for A/V for buildings is assigned to the grid cell containing the highest
! purb(c,r) i.e., 100.0. A/V for buildings for other grid-cell is calculated
! as purb(c,r)*(0.2/100.0); Cai et al. (2006) used a value of 1.0 for their
! study at New York (total A/V)

                  if ( ifurban ) then
                     surf_bldg = 0.2
                  else
                     surf_bldg = 0.0
                  endif

! Calculate rate constant for the reaction (psudeo-first order reaction,
! unit per second). Calculate pseudo-first order rate constant using Eq 1
! of Vogel et al. (2003).  Unit of KNO2 is in 1/min in the paper; divide it
! by 60 to convert it into 1/sec.

                  kno2 = MAX( 0.0, 5.0E-5 * (surf_leaf + surf_bldg) )

! Determine NO2 concentration needed for HONO production term.

!                 IF ( depspc( l ) .EQ. 'NO2' ) THEN
                  IF ( l .EQ. l_no2 ) THEN

! Loss of NO2 via the heterogeneous reaction is accounted as additional
! depositional loss. Add the loss of NO2 via the heterogeneous reaction
! to the regular deposition velocity (increased dep. vel.).  This will
! reduce the NO2 conc. in the atmosphere. Dep vel is adjusted back to the
! original value in vdiffacm2 after NO2 conc is reduced but before calculating
! depositional loss.

                     depvel_gas_land( iwl,n ) = depvel_gas_land( iwl,n ) + 2.0 * kno2 * zf

                  END IF

!                 IF ( depspc( l ) .EQ. 'HONO' ) then
                  IF ( l .EQ. l_hono ) then

C Calculate production (pvd) for HONO; unit = ppm * m/s
                     
                     pvd = kno2 * cgrid(kw,iw,ic_no2) * zf
                     
! Add this to the HONO emissions for now when we apply it in vertical diffusion

                     vdemis_gc(kw,iw,ie_hono) = vdemis_gc(kw,iw,ie_hono) + pvd * land%area(iwl) * volti(kw,iw)

                  ENDIF

               END IF   ! sfc_hono
               
            END IF   ! special condition for sulfate (SULF)

            ! Check for negative values or NaN's

            IF ( depvel_gas_land( iwl,n ) .LT. 0.0 .OR.
     &           depvel_gas_land( iwl,n ) .NE. depvel_gas_land( iwl,n ) ) then

               xmsg = 'NEGATIVE or UNDEFINED Dry Deposition Velocity for ' // trim(depspc(l))
               call m3exit( pname, 0, 0, xmsg, xstat1 )

            END IF
            
         END DO dloopl   ! (l = 1, n_spc_m3dry)

      enddo   ! land

C-------------------------------------------------------------------------------
C Loop over sea cells and calculate dry deposition.
C-------------------------------------------------------------------------------

      do iws = 2, mws

         iw = itab_ws(iws)%iw
         if (isubdomain == 1) then
            iw = itabg_w(iw)%iw_myrank
         endif
         kw = itab_ws(iws)%kw

         isocean = (sea%leaf_class(iws) == 0)

         ustar  = sea%ustar(iws)
         ustarw = sea%sea_ustar(iws)
         raw    = 1.0 / sea%sea_ggaer(iws)

         sst   = sea%seatc(iws)
         tw    = sea%sea_cantemp(iws) ! water surface film temperature
         
         wtv0  = sea%sfluxt(iws) * (1. + .61 * sh_v(kw,iw))
     &         + sea%sfluxr(iws) * .61 * theta(kw,iw)

         if (wtv0 > 0.0) then
            wstar = (grav * pblh(iw) * wtv0 / theta(kw,iw)) ** 0.33333333
         else
            wstar = 0.0
         endif

         if  (sea%nlev_seaice(iws) > 0 ) then
            ustari = sea%ice_ustar(iws)
            raice  = 1.0 / sea%ice_ggaer(iws)
            tice   = sea%ice_cantemp(iws) ! ice surface film temperature

            ! Determine the seaice liquid water mass fraction (0.0 to 0.5).

            call qtk_sea( sea%seaice_energy(sea%nlev_seaice(iws),iws), xt, xm)
            xm = MIN (xm, 0.5)
            xm = MAX (xm, 0.0)
         endif

         ! Loop over species to calculate dry deposition velocities.

         depvel_gas_sea( iws,: ) = 0.0  ! initialize for this time period

         n = 0

         dloops: DO l = 1, n_spc_m3dry
         
            IF ( .NOT. use_depspc( l ) ) CYCLE dloops

            n = n + 1

!           IF ( depspc( l ) .EQ. 'SULF' ) THEN  ! Sulfate (SULF)
            IF ( l .EQ. l_sulf ) THEN  ! Sulfate (SULF)

         ! Sulfate calculation follows Wesely (1985), Eqn. 11.

               rbsulf = 500. * ustar / (ustar**2 + 0.24 * wstar**2)
               depvel_gas_sea( iws, n ) =  1.0 / ( raw + rbsulf )
               
            else

         ! Use CMAQ function for calculating the effective Henry's Law
         ! constant.  Note that original M3DRY wants inverse, non-
         ! dimensional Henry's Law (caq/cg).   Water pH is different
         ! than rain, and we need to use the water temperature.

               heff  = hlconst( tw, effective, hplus_h2o, hlspc( l ) )

         ! Make Henry's Law constant non-dimensional.

               heff  = heff * 0.08205 * tw

         ! from Hayduk and Laudie

               dw25      = 13.26e-5 / ( 0.8904**1.14 * lebas( l )**0.589 )
               kvisw     = 0.017 * EXP( -0.025 * ( tw - stdtemp ) )
               dw        = dw25 * ( tw * rt25inK ) * ( 0.009025 / kvisw )
               scw_pr_23 = ( ( kvisw / dw ) / pr ) ** twothirds

               if (l == l_o3) then !implement Chang et al(2004)
                  
c        pChang is a/H or alpha/H which would be 1/H in current model
c        note that in Chang et al (2004) and Garland et al (1980) their H is Cair/Cwater with is
c        the inverse of heff

                  pChang = 1.75
                  kwChang = d3 * ustarw / scw_pr_23

c        If a file of chlorophyll concentrations is provided, Iodide concentration are estimated from
c        a fit to the Rebello et al 1990 data. The slope and correlation are given in the paper
c        but not the intercept, so the data in Tables 3 & 4 were fit to get the relationship below.
c        The regression gives the concentration in umol/L and is converted to mol/L for use in Chang et al eq.
c        The slope and correlation are a slightly different than in Table 5.
c        If chlorophyll concs are not available, a constant value for [I-] of 100e-9 mol/l is used
c        Use ocean file variables to determine if the water cell is ocean or lake; method is only for ocean cells

                  if (isocean) then
c                    Iodide in sea-water based on SST  (mol /dm-3)
                     ciodide = 1.46E6 * EXP( -9134.0 / sst )
                     qiodide = sqrt ( 2.0e9 * ciodide * dw * 1.e-4 ) * heff
                  else
                     qiodide = 0.0
                  endif

                  rsurf   = 1.0 / ( pChang * kwchang + qiodide )

               else  ! ozone

                  rsurf = scw_pr_23 / ( heff * d3 * ustarw )

               endif

         ! Compute dry deposition velocity over sea

               rbc = 5.0 / ustarw * scc_pr_23( l )

               depvel_gas_sea( iws,n ) = 1.0 / ( rsurf + raw + rbc )

         ! Include fractional contribution of seaice if present
               
               if  (sea%nlev_seaice(iws) > 0 ) then
                  
         ! Use CMAQ function for calculating the effective Henry's Law
         ! constant.  Note that original M3DRY wants inverse,
         ! non-dimensional Henry's Law (caq/cg).

                  heff  = hlconst( tice, effective, hplus_h2o, hlspc( l ) )

         ! Make Henry's Law constant non-dimensional

                  heff  = heff * 0.08205 * tw

         ! Wet surface resistance

                  rgw  = rgwet0 / heff
            
         ! Dry snow resistance

                  rsnow = rsnow0 * a0 / ar( l )

         ! Assume seaice behaves similarly to snowcover
                  
                  rsurf = 1.0 / ( (1.0 - xm) / rsnow + xm / (rsndiff + rgw) )

         ! Compute dry deposition velocity contribution from seaice

                  rbc   = 5.0 / ustari * scc_pr_23( l )
                  d_ice = 1.0 / ( rsurf + raice + rbc )

                  depvel_gas_sea( iws,n ) = (1.0 - sea%seaicec(iws)) * depvel_gas_sea( iws,n )
     &                                    +        sea%seaicec(iws)  * d_ice

               endif  ! seaice
               
            endif  ! special condition for sulfate (SULF)

         ! Check for negative values or NaN's

            IF ( depvel_gas_sea( iws,n ) .LT. 0.0 .OR.
     &           depvel_gas_sea( iws,n ) .NE. depvel_gas_sea( iws,n ) ) then

               xmsg = 'NEGATIVE or UNDEFINED Dry Deposition Velocity for ' // trim(depspc(l))
               call m3exit( pname, 0, 0, xmsg, xstat1 )

            END IF
            
         END DO dloops   ! (l = 1, n_spc_m3dry)

      enddo   ! sea

      END SUBROUTINE m3dry
