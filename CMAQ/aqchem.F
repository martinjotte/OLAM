
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
C $Header: /project/yoj/arc/CCTM/src/cloud/cloud_acm_ae6_mp/aqchem.F,v 1.4 2012/03/28 15:39:52 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AQCHEM ( kw, iw, TEMP, PRES_PA, TAUCLD, WCAVG, GAS, AEROSOL,
     &                    GASWDEP, AIRM, HP, HPWDEP, EWASH, DARK )
C-----------------------------------------------------------------------
C  Description:
C    Compute concentration changes in cloud due to aqueous chemistry,
C    scavenging and wet deposition amounts.
C
C  Revision History:
C      No   Date   Who  What
C      -- -------- ---  -----------------------------------------
C      0  / /86    CW   BEGIN PROGRAM - Walceks's Original Code
C      1  / /86    RB   INCORPORATE INTO RADM
C      2  03/23/87 DH   REFORMAT
C      3  04/11/88 SJR  STREAMLINED CODE - ADDED COMMENTS
C      4  08/27/88 SJR  COMMENTS, MODIFIED FOR RPM
C      4a 03/15/96 FSB  Scanned hard copy to develop Models3
C                       Version.
C      5  04/24/96 FSB  Made into Models3 Format
C      6  02/18/97 SJR  Revisions to link with Models3
C      7  08/12/97 SJR  Revised for new concentration units (moles/mole)
C                       and new treatment of nitrate and nitric acid
C      8  01/15/98 sjr  revised to add new aitken mode scavenging
C                       and aerosol number scavenging
C      9  12/15/98 David Wong at LM:
C             -- change division of XL, TEMP to multiplication of XL, TEMP
C                reciprocal, respectively
C             -- change / TOTOX / TSIV to / ( TOTOX * TSIV )
C     10  03/18/99 David Wong at LM:
C             -- removed "* 1.0" redundant calculation at TEMP1 calculation
C     11  04/27/00 sjr  Added aerosol surface area as modeled species
C     12  12/02    sjr  changed calls to HLCONST and updated the dissociation
C                       constants
C     13  06/26/03 sjr  revised calculations of DTW based on CMAS website
C                       discussions
C     14  08/05/03 sjr  revision made to the coarse aerosol number washout
C     15  04/20/05  us  revisions to add sea salt species in the fine and
C                       coarse aerosol modes, and HCl dissolution/dissociation
C     16  10/13/05 sjr  fixed bug in the integration time step calculation
C                       (reported by Bonyoung Koo)
C     17  03/01/06 sjr  added elemental carbon aerosol; organic aerosols
C                       replaced with primary, secondary biogenic, and
C                       secondary anthropogenic; fixed 3rd moment calc to
C                       include EC and primary organics (not secondary);
C                       re-arranged logic for setting Cl & Na ending conc;
C                       added pointers/indirect addressing for arrays WETDEP
C                       and LIQUID
C     16  03/30/07 sjr  Limit integration timestep by cloud washout time
C     17  04/10/07 sjr  increased loop limits as follows: I20C <10000,
C                       I7777C <10000, I30C <10000, ICNTAQ <60000
C     18  01/10/07 agc  added organic chemistry for GLY and MGLY oxidation
C     19  09/10/07 sln  updated SOA species list for AE5
C     20  01/29/08 agc  updated DOHDT calculation
C     21  04/14/08 jtk  added coding for coarse NH4 and scavenging of
c                       coarse surface area
C     22  05/20/08 agc  for CB05, use the Henry's Law constant for glyoxal
C                       as a surrogate for methyl glyoxal
C     23  04/15/09 sjr& Several changes made to improve mass conservation in the
C                  agc  solver.  (1) OH concentration is now considered to be
C                       steady state; (2) only allow sulfur oxidation to affect
C                       time step; (3) implemented mass conservation checks -
C                       limit oxidation rates by the available mass for the
C                       specified timestep.
C   10 Oct 10 J.Young:  update to use aero_reeng by Steve Howard, Prakash Bhave,
C                       Jeff Young, Sergey Napelenok, and Shawn Roselle
C   01 Mar 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
C    9 Mar 11 S.Napelenok: update for AE6 - pH calculation now expanded to
C                       include Ca Mg K SOIL CORS SEAS
C   23 May 11 G.Sarwar: update S(VI) production rate via H2O2, O3, MHP, PAA
C                       pathways (Jacobson 1997)
C   23 May 11 G.Sarwar: update S(VI) production rate via O2 pathway (metal
C                       catalysis) (Martin and Goodman, 1991)
C   01 Jul 11 G.Sarwar: Incorporate day and night dependent Fe III oxidation
C                       state (Alexander et al.,  2009)
C   12 Aug 11 G.Sarwar: Revise Fe and Mn solubility based on
C                       Alexander et al., 2009
C
C    8 Mar 12 J.Bash:   FE_OX and MN_OX were calculated from FE and MN before
C                       a floor value of 0.0 was established for these
C                       concentrations sometimes resulting in negative
C                       concentrations and model crashes. The code used to
C                       estimate FE_OX and MN_OX was moved to be after a floor
C                       value for FE and MN was set. Also the washout rate was
C                       removed from the calculation of the estimate for doubling
C                       the time step based on sulfur oxidized < 5%.
C   28 Nov 12 G.Sarwar: Sulfate inhibition effect is implemented in the metal catalysis pathway
C   07 Jul 14 B.Hutzell: replaced mechanism include file(s) with fortran module
C   12 Feb 15 B.Hutzell: reduced number of exp(...) calculations for scavenging aitken
C                        aerosols to improve efficiency
C   15 Jun 15 J.Young:  Fixed bug found by Martin Otte in calculations for scavenging
C                       aitken aerosols
C   15 Apr 16 J.Young:  Use aerosol factors from AERO_DATA module named constants
C
C  References:
C     Walcek & Taylor, 1986, A theoretical Method for computing
C        vertical distributions of acidity and sulfate within cumulus
C        clouds, J. Atmos Sci.,  Vol. 43, no. 4 pp 339 - 355
C     Carlton, A.G., B.J. Turpin, K.E. Altieri, S.P. Seitzinger, R. Mathur,
C        S.J. Roselle, and R.J. Weber, CMAQ Model Performance Enhanced When
C        In-Cloud Secondary Organic Aerosol is Included:  Comparison of Organic
C        Carbon Predictions with Measurements, Environ. Sci. Technol., 42(23),
C        8798-8802, 2008.
C     Jacobson, M., Development and application of a new air pollution modeling
C        system II. Aerosol module structure and design, Atmospheric
C        Environment, 31, 131-144, 1997
C     Martin, R.L. and T.W. Good, catalyzed oxidation of sulfur dioxide in
C        solution: the iron-manganese synercism, Atmospheric Environment, 25A,
C        2395-2399, 1991
C     Alexander, B., R.J. Park, D.J. jacob, S. Gong, Transition metal-catalyzed
C        oxidation of atmospheric sulfur: global implications for the sulfur
C        budget, GRL, 114, D02309, 2009

C  Called by:  AQMAP

C  Calls the following subroutines:  none

C  Calls the following functions:  HLCONST

C  Arguments     Type      I/O       Description
C  ---------     ----  ------------  --------------------------------
C  GAS(ngas)     real  input&output  Concentration for species i=1,15
C  GASWDEP(ngas) real     output     wet deposition for species

C  AEROSOL(naer,nmodes) real input&output   Concentration for species i=1,51
C  AERWDEP(naer,nmodes) real     output     wet deposition for species
C-----------------------------------------------------------------------

      USE AQ_DATA         ! doesn't inherit; gets only n_aerospc, conmin from AERO_DATA
      USE AERO_DATA
      USE UTILIO_DEFN
      use hlconst_mod
      use const_data
      use mem_ijtabs, only: itab_w

      IMPLICIT NONE

      CHARACTER( 120 ) :: XMSG   ! Exit status message

C...........Parameters:

      INTEGER, PARAMETER :: NUMOX   = 5         ! number of oxidation reactions
      REAL(8), PARAMETER :: H2ODENS = 1000.0D0  ! water density at 20 C and 1 ATM (kg/m3)

      REAL, PARAMETER :: ONETHIRD = 1.0 / 3.0
      REAL, PARAMETER ::  LN10    = LOG(10.0)
      REAL, PARAMETER :: MLN10    = -LN10

C...........Arguments:

      integer, intent(in) :: kw, iw

      REAL,      INTENT( OUT ) :: HP        ! h+ concentration in cloud water (mol/liter)
      REAL,      INTENT( OUT ) :: HPWDEP    ! hydrogen wet deposition (mm mol/liter)
!     REAL,      INTENT( OUT ) :: EBETASO4T
      REAL,      INTENT( IN )  :: EWASH
      REAL,      INTENT( IN )  :: AIRM
      REAL,      INTENT( IN )  :: PRES_PA   ! pressure (Pa)
      REAL,      INTENT( IN )  :: TAUCLD    ! timestep for cloud (s)
      REAL,      INTENT( IN )  :: TEMP      ! temperature (K)
      REAL,      INTENT( IN )  :: WCAVG     ! liquid water content (kg/m3)
      REAL( 8 ), INTENT( INOUT ) :: GAS    ( NGAS )   ! gas phase concentrations (mol/molV)
      REAL     , INTENT( INOUT ) :: GASWDEP( NGAS )   ! gas phase wet deposition array (mm mol/liter)
      REAL( 8 ), INTENT( INOUT ) :: AEROSOL( NAER,2:NMODES ) ! aerosol concentrations (mol/molV)
!     REAL( 8 ), INTENT( INOUT ) :: AERWDEP( NAER,NMODES ) ! aerosol wet deposition array (mm mol/liter)
      LOGICAL,   INTENT( IN )    :: DARK           ! DARK = TRUE is night,  DARK = FALSE is day
!     REAL,      INTENT( IN )    :: CTHK1     ! cloud thickness (m)

C...........Local Variables (scalars):

!     LOGICAL, SAVE :: FIRSTIME = .TRUE. ! flag for first pass thru

      CHARACTER( 16 ), SAVE :: PNAME = 'AQCHEM'             ! driver program name

      integer :: I
      INTEGER :: ISPC            ! loop counter for species
      INTEGER :: I30C            ! loop counter for do loop 30
!     INTEGER :: LIQ             ! loop counter for liquid species
!     INTEGER :: IGAS            ! loop counter for gas species
!     INTEGER :: IOX             ! index over oxidation reactions
      INTEGER :: IMODE

!     REAL :: DEPSUM
!     REAL :: so4avail

      REAL( 8 ) :: A               ! iron's anion concentration
      REAL( 8 ) :: HPLUS           ! H+ concentration in cloudwater (mol/liter)
      REAL      :: ACT1            ! activity correction factor, single ions
      REAL      :: ACT2            ! activity factor correction, double ions
      REAL( 8 ) :: B               ! manganese's anion concentration
      REAL( 8 ) :: PRES_ATM        ! pressure (Atm)
      REAL      :: BB, bb1, bb2    ! cloudwater pH guesses
      REAL( 8 ) :: CA              ! Calcium conc in cloudwater (mol/liter)
      REAL( 8 ) :: CL              ! total Cl-  conc in cloudwater (mol/liter)
      REAL( 8 ) :: CLACC           ! fine Cl- in cloudwater (mol/liter)
      REAL( 8 ) :: CLCOR           ! coarse Cl-  conc in cloudwater (mol/liter)
      REAL( 8 ) :: CO2H            ! Henry's Law constant for CO2
!     REAL( 8 ) :: CO21            ! First dissociation constant for CO2
!     REAL( 8 ) :: CO22            ! Second dissociation constant for CO2
      REAL( 8 ) :: CO212           ! CO21*CO22
      REAL( 8 ) :: CO212H          ! CO2H*CO21*CO22
      REAL( 8 ) :: CO21H           ! CO2H*CO21
      REAL( 8 ) :: CO212A
      REAL( 8 ) :: CO212B
      REAL( 8 ) :: CO2L            ! CO2 conc in cloudwater (mol/liter)
      REAL( 8 ) :: CO3             ! CO3= conc in cloudwater (mol/liter)
      REAL( 8 ) :: DSIV_SCALE      ! mass conservation scale factor for S(IV)
!     REAL( 8 ) :: DTRMV           !
!     REAL( 8 ) :: DTS6            !
!     REAL( 8 ) :: DGLYDT          ! change in GLY (mol/liter/sec)
!     REAL( 8 ) :: DMGLYDT         ! change in MGLY (mol/liter/sec)
!     REAL( 8 ) :: DOHDT           ! change in OH
!     REAL( 8 ) :: DGLY1           ! change due to Rxn. in GLY for DTW(0) time step
!     REAL( 8 ) :: DMGLY1          ! change due to Rxn. in MGLY for DTW(0) time step
!     REAL( 8 ) :: DOH1            ! change in OH for DTW(0) time step
!     REAL( 8 ) :: DORGC           ! change in ORGC for DTW(0) time step (mol/liter)
!     REAL( 8 ) :: EBETASO4T       ! EXP( -BETASO4 * TAUCLD )
!     REAL( 8 ) :: EC              ! elemental carbon acc mode aerosol in cloudwater (mol/liter)
      REAL      :: FA              ! functional value ??
      REAL      :: FB              ! functional value ??
      REAL( 8 ) :: FCLCOR          ! frac weight of coarse CL to (acc+coarse) CL
      REAL      :: FN              ! functional value ??
      REAL( 8 ) :: FE              ! Fe+++ conc in cloudwater (mol/liter)
      REAL( 8 ) :: FNH3            ! frac weight of NH3 to total ammonia
      REAL( 8 ) :: FNH4ACC         ! frac weight of NH4 acc to total ammonia
      REAL( 8 ) :: FNH4COR         ! frac weight of coarse NH4 to (acc+coarse) NH4
      REAL( 8 ) :: FHNO3           ! frac weight of HNO3 to total NO3
      REAL( 8 ) :: FNO3ACC         ! frac weight of NO3 acc to total NO3
      REAL( 8 ) :: FNO3COR         ! frac weight of coarse NO3 to (acc+coarse) NO3
!     REAL( 8 ) :: FOA1            ! First dissociation constant for FOA (Formic Acid)
      REAL( 8 ) :: FOA1A
      REAL( 8 ) :: FOA1B
      REAL( 8 ) :: FOAH            ! Henry's Law constant for FOA
      REAL( 8 ) :: FOA1H           ! FOAH*FOA1
      REAL( 8 ) :: FOAL            ! FOA conc in cloudwater (mol/liter)
!     REAL( 8 ) :: FTST            !
      REAL( 8 ) :: GLYH            ! Henry's Law constant for glyoxal
      REAL( 8 ) :: GLYL            ! glyoxal conc in cloud water (mol/liter)
!     REAL      :: GM              !
      REAL      :: GM1             !
      REAL      :: GM1LOG          !
      REAL      :: GM2             ! activity correction factor
!     REAL( 8 ) :: GM2LOG          !
      REAL      :: HA              !
      REAL      :: HB              !
!     REAL( 8 ) :: H2OW            !
      REAL( 8 ) :: H2O2H           ! Henry's Law Constant for H2O2
      REAL( 8 ) :: H2O2L           ! H2O2 conc in cloudwater (mol/liter)
      REAL( 8 ) :: HCLH            ! Henry's Law Constant for HCL
!     REAL( 8 ) :: HCL1            ! First dissociation constant for HCL
      REAL( 8 ) :: HCL1H           ! HCL1*HCLH
      REAL( 8 ) :: HCLL            ! HCl  conc in  cloudwater (mol/liter)
      REAL( 8 ) :: HCLA
      REAL( 8 ) :: HCLB
      REAL( 8 ) :: HCO2            ! HCO2 conc in cloudwater (mol/liter)
      REAL( 8 ) :: HCO3            ! HCO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: HNO3H           ! Henry's Law Constant for HNO3
!     REAL( 8 ) :: HNO31           ! First dissociation constant for HNO3
      REAL( 8 ) :: HNO3A
      REAL( 8 ) :: HNO3B
      REAL( 8 ) :: HNO31H          !
      REAL( 8 ) :: HNO3L           ! HNO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: HOH             ! Henry's Law Constant for HO
      REAL( 8 ) :: HSO3            ! HSO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: HSO4            ! HSO4 concn in cloudwater (mol/liter)
!     REAL( 8 ) :: HSO4ACC         ! accumulation mode HSO4 concn in cloudwater (mol/liter)
!     REAL( 8 ) :: HSO4COR         ! coarse HSO4 concn in cloudwater (mol/liter)
!     REAL( 8 ) :: HTST            !
      REAL( 8 ) :: K               ! K conc in cloudwater (mol/liter)
!     REAL( 8 ) :: LGTEMP          ! log of TEMP
      REAL( 8 ) :: MG              !
      REAL( 8 ) :: MGLYH           ! Henry's Law Constant for methylglyoxal
      REAL( 8 ) :: MGLYL           ! MGLY conc in cloud water (mol/liter)
      REAL( 8 ) :: MHPH            ! Henry's Law Constant for MHP
      REAL( 8 ) :: MHPL            ! MHP conc in cloudwater (mol/liter)
      REAL( 8 ) :: MN              ! Mn++ conc in cloudwater (mol/liter)
      REAL( 8 ) :: NA              ! Na conc in cloudwater (mol/liter)
!     REAL( 8 ) :: NAACC           ! Na in cloudwater (mol/liter)
      REAL( 8 ) :: NACOR           ! coarse Na in cloudwater (mol/liter)
!     REAL( 8 ) :: NH31            ! First dissociation constant for NH3
      REAL( 8 ) :: NH3H            ! Henry's Law Constant for NH3
      REAL( 8 ) :: NH3DH20         !
      REAL( 8 ) :: NH31HDH         !
      REAL( 8 ) :: NH31A           !
      REAL( 8 ) :: NH31B           !
      REAL( 8 ) :: NH3L            ! NH3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: NH4             ! NH4+ conc in cloudwater (mol/liter)
      REAL( 8 ) :: NH4ACC          ! NH4 acc conc in cloudwater (mol/liter)
      REAL( 8 ) :: NH4COR          ! NH4 coarse conc in cloudwater (mol/liter)
!     REAL( 8 ) :: NITAER          ! total aerosol nitrate
      REAL( 8 ) :: NO3             ! NO3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: NO3ACC          ! NO3 acc conc in cloudwater (mol/liter)
      REAL( 8 ) :: NO3COR          ! NO3 coarse conc in cloudwater (mol/liter)
!     REAL( 8 ) :: NUMCOR          ! coarse aerosol number in cloudwater (mol/liter)
      REAL( 8 ) :: O3H             ! Henry's Law Constant for O3
      REAL( 8 ) :: O3L             ! O3 conc in cloudwater (mol/liter)
      REAL( 8 ) :: OH              ! OH conc in cloudwater (mol/liter)
      REAL( 8 ) :: OHL             ! OH radical conc in cloudwater (mol/liter)
!     REAL( 8 ) :: SOA             ! secondary organic aerosol in cloudwater (mol/liter)
      REAL( 8 ) :: ORGC            ! cloud-produced SOA in cloudwater (treated as primary)
!     REAL( 8 ) :: POA             ! primary organic aerosol in cloudwater (mol/liter)
      REAL( 8 ) :: PAAH            ! Henry's Law Constant for PAA
      REAL( 8 ) :: PAAL            ! PAA conc in cloudwater (mol/liter)
      REAL( 8 ) :: PCO2F           ! gas only CO2 partial pressure (atm)
      REAL( 8 ) :: PFOAF           ! gas only ORGANIC ACID partial press (atm)
      REAL( 8 ) :: PGLYF           ! gas only GLY partial pressure (atm)
      REAL( 8 ) :: PH2O2F          ! gas only H2O2 partial pressure (atm)
      REAL( 8 ) :: PHCLF           ! gas only HCL partial pressure (atm)
      REAL( 8 ) :: PHNO3F          ! gas only HNO3 partial pressure (atm)
      REAL( 8 ) :: PHOF            ! gas only HO partial pressure (atm)
      REAL( 8 ) :: PMGLYF          ! gas only MGLY parital pressure (atm)
      REAL( 8 ) :: PMHPF           ! gas only MHP partial pressure (atm)
      REAL( 8 ) :: PNH3F           ! gas only NH3 partial pressure (atm)
      REAL( 8 ) :: PO3F            ! gas only O3 partial pressure (atm)
      REAL( 8 ) :: PPAAF           ! gas only PAA partial pressure (atm)
!     REAL( 8 ) :: PRIM            ! PRIMARY acc aerosol in cloudwater (mol/liter)
!     REAL( 8 ) :: PRIMCOR         ! PRIMARY coarse aerosol in cloudwater (mol/liter)

      REAL( 8 ) :: PSO2F           ! gas only SO2 partial pressure (atm)
!     REAL( 8 ) :: RATE            !
!     REAL( 8 ) :: RECIPA1         !
!     REAL( 8 ) :: RECIPA2         !
      REAL( 8 ) :: RECIPAP1        ! one over pressure (/atm)
!     REAL( 8 ) :: RH2O2           !
!     REAL( 8 ) :: RMHP            !
!     REAL( 8 ) :: RPAA            !
      REAL( 8 ) :: RT              ! gas const * temperature (liter atm/mol)
      REAL( 8 ) :: SIV             ! dissolved so2 in cloudwater (mol/liter)
!     REAL( 8 ) :: SK6             !
      REAL( 8 ) :: SK6TS6          !
!     REAL( 8 ) :: SO21            ! First dissociation constant for SO2
!     REAL( 8 ) :: SO22            ! Second dissociation constant for SO2
      REAL( 8 ) :: SO2H            ! Henry's Law Constant for SO2
      REAL( 8 ) :: SO212           ! SO21*SO22
      REAL( 8 ) :: SO212H          ! SO21*SO22*SO2H
      REAL( 8 ) :: SO21H           ! SO21*SO2H
      REAL( 8 ) :: SO21A
      REAL( 8 ) :: SO21B
      REAL( 8 ) :: SO21C
      REAL( 8 ) :: SO2L            ! SO2 conc in cloudwater (mol/liter)
      REAL( 8 ) :: SO3             ! SO3= conc in cloudwater (mol/liter)
      REAL( 8 ) :: SO4             ! SO4= conc in cloudwater (mol/liter)
!     REAL( 8 ) :: SO4ACC          ! accumulation mode SO4= conc in cloudwater (mol/liter)
!     REAL( 8 ) :: SO4COR          ! coarse SO4= conc in cloudwater (mol/liter)
      REAL      :: STION           ! ionic strength
!     REAL( 8 ) :: TAC             !
      REAL( 8 ) :: TCLa            ! sum of accumulation and coarse mode chloride
      REAL      :: TEMP1           ! (1/T) - (1/298) (1/K)
!     REAL( 8 ) :: TIMEW           ! cloud chemistry clock (sec)
!     REAL( 8 ) :: THO             ! total hydroxyl radical available for oxidation
!     REAL( 8 ) :: TGLY            ! total glyoxal available for oxidation
!     REAL( 8 ) :: TMGLY           ! total methylglyoxal available for oxidation
!     REAL( 8 ) :: TOTOX           !
!     REAL( 8 ) :: TH2O2
!     REAL( 8 ) :: TO3
!     REAL( 8 ) :: TMHP
      REAL( 8 ) :: TNH4a           ! sum of accumulation and coarse mode ammonium
      REAL( 8 ) :: TNO3a           ! sum of accumulation and coarse mode nitrate
!     REAL( 8 ) :: TPAA
      REAL( 8 ) :: TOTAMM          ! total ammonium
      REAL( 8 ) :: TOTNIT          ! total nitrate (excluding coarse mode)
      REAL( 8 ) :: TS6             ! SO4 conc in cloudwater (mol/liter)
      REAL( 8 ) :: TS6ACC          ! SO4 acc conc in cloudwater (mol/liter)
      REAL( 8 ) :: TS6COR          ! coarse SO4 conc in cloudwater   (mol/liter)
      REAL( 8 ) :: TSIV            ! total S(iv) available for oxidation
      REAL      :: EFACT
!     REAL( 8 ) :: WETFAC          ! converts mol/l to mm-mol/l based on precip
      REAL      :: XC              ! (/mm)
!     REAL( 8 ) :: XC1             ! (/mm)
!     REAL( 8 ) :: XC2             ! (liter-atm/mol/mm)
      REAL( 8 ) :: XL              ! conversion factor (liter-atm/mol)
      REAL( 8 ) :: ONE_OVER_XL     ! 1.0 / XL
      REAL( 8 ) :: PRES_ATM_OVER_XL ! PRES_ATM / XL
!     REAL( 8 ) :: SCAVENGED       ! aitken scavenging factor by cloud water
      REAL( 8 ) :: XLCO2           !
      REAL( 8 ) :: XLH2O2          !
      REAL( 8 ) :: XLHCL           ! const in calc of HCL final partial pres
      REAL( 8 ) :: XLHNO3          !
      REAL( 8 ) :: XLMHP           !
      REAL( 8 ) :: XLNH3           !
      REAL( 8 ) :: XLO3            !
      REAL( 8 ) :: XLPAA           !
      REAL( 8 ) :: XLSO2           !
!     REAL( 8 ) :: CAACC           ! accumulation mode Calcium (AE6) SLN 16March2011
!     REAL( 8 ) :: MGACC           ! accumulation mode Magnesium (AE6) SLN 16March2011
!     REAL( 8 ) :: KACC            ! accumulation mode Potassium (AE6) SLN 16March2011
      REAL( 8 ) :: CACOR           ! coarse mode Calcium (AE6) SLN 16March2011
      REAL( 8 ) :: MGCOR           ! coarse mode Magnesium (AE6) SLN 16March2011
      REAL( 8 ) :: KCOR            ! coarse mode Potassium (AE6) SLN 16March2011
      REAL( 8 ) :: SOILCOR         ! coarse mode SOIL (AE6) SLN 16March2011
      REAL( 8 ) :: ANTHCOR         ! coarse mode CORS (AE6) SLN 16March2011
      REAL( 8 ) :: SEASCOR         ! coarse mode SEAS (AE6) SLN 16March2011
!     REAL( 8 ) :: FEACC           ! accumulation mode Fe (AE6) SLN 22March2011
!     REAL( 8 ) :: MNACC           ! accumulation mode Mn (AE6) SLN 22March2011
      REAL( 8 ) :: FECOR           ! coarse mode Fe (AE6) SLN 22March2011
      REAL( 8 ) :: MNCOR           ! coarse mode Mn (AE6) SLN 22March2011
      REAL( 8 ) :: FE_OX           ! Fe(III) available for sulfate oxidation
      REAL( 8 ) :: MN_OX           ! Mn(II) available for sulfate oxidation
      REAL( 8 ) :: FE_III          ! Fractional Fe(III) partitioning, GS - July 1, 2011
      REAL( 8 ) :: MN_II           ! Fractional Mn(II) partitioning, GS - July 1, 2011
      REAL( 8 ) :: BI0
      REAL( 8 ) :: ST0
      REAL( 8 ) :: FE_SOL          ! Fractional Fe solubility, GS - July 1, 2011
      REAL( 8 ) :: MN_SOL          ! Fractional Mn solubility, GS - July 1, 2011
      real( 8 ) :: f2o2f, fo3f, fmhpf, fpaaf, fglyf, fmglyf, fhof
      real      :: ts6_prev, so2_prev

!     real(8) :: c1, c2, c3, dhoh
!     real(8) :: s1, s2, s3
!     real(8) :: dso3, dso4, dhso3, dhso4, ds
!     real(8) :: dco3, dhco3, dnh4, dhco2, doh

!For Varaible used by TXHG Version
!     REAL( 8 ) :: TRACER           ! TRACER acc aerosol in cloudwater (mol/liter)
!     REAL( 8 ) :: TRACERCOR        ! TRACER coarse aerosol in cloudwater (mol/liter)
!     REAL( 8 ) :: HGFINE           ! mercury PM acc aerosol in cloudwater (mol/liter)
!     REAL( 8 ) :: HGCOR            ! mercury PM coarse aerosol in cloudwater (mol/liter)

C...........Local Variables (arrays):

      REAL( 8 ) :: LOADING( NAER, 2:NMODES ) ! aerosol loading (mol/liter)
      REAL( 8 ) :: INITGAS( NGAS ) ! initial gas partial pressure (atm)
!     REAL( 8 ) :: LIQUID( NLIQS ) ! wet deposition array (mm mol/liter)
!     REAL( 8 ) :: WETDEP( NLIQS ) ! wet deposition array (mm mol/liter)
      REAL      :: DSIVDT( 0:NUMOX+2 ) ! rate of so2 oxid incloud (mol/liter/sec)
      REAL( 8 ) :: DS4   ( 0:NUMOX+2 ) ! S(IV) oxidized over timestep DTW(0)
      REAL( 8 ) :: SSIV  ( 0:NUMOX+2 ) ! S(IV)

      real( 8 ) :: a1, a2, phclfoa1, phno3foa1, pfoafoa1
      real( 8 ) :: PSO2Foa1a2, PCO2Foa1a2, a1a2
      real( 8 ) :: tfoa, tco2

      real :: rates(15)

      real :: rr0(15) = (/ 1.80E-04, 1.02E-02, 1.30E-02, 6.60E-08, 4.30E-07, 4.68E-11,
     &                     1.00E-14, 1.70E-05, 1.74E+06, 1.54E+01, 7.45E+07, 1.90E+07,
     &                     3.60E+07, 3.70E+05, 1.50E+09 /)

      real :: rr1(15) = (/-2.00E+01, 2.72E+03, 1.96E+03, 1.50E+03,-1.00E+03,-1.76E+03,
     &                    -6.71E+03,-4.50E+02, 6.90E+03, 8.70E+03,-4.76E+03,-3.80E+03,
     &                    -4.00E+03,-5.53E+03,-5.28E+03 /)

      ! From Carlton et al. (2007)
      REAL, PARAMETER :: RGLY3  = 3.0E10   ! liter/(mol sec)
      REAL, PARAMETER :: RMGLY3 = 3.0E10   ! liter/(mol sec)

      REAL, PARAMETER :: T0  = 298.0
      REAL, PARAMETER :: T0M = 1.0 / T0

      REAL, PARAMETER :: ONEM = NEAREST(1.0, -1.0)
      REAL      :: EWASHM, XL_EWASH, XL_EWASHM
      LOGICAL   :: DO_WETDEP
      real      :: faa, fbb, haa, hba
!     real :: hcl1

C...initialize counters and compute several conversion factors

      RT         = ( MOLVOL / STDTEMP ) * TEMP         ! R * T (liter atm / mol)
      PRES_ATM   = PRES_PA /  STDATMPA                 ! pressure (atm)
      XL         = WCAVG * RT / H2ODENS                ! conversion factor (l-atm/mol)

      ONE_OVER_XL      = 1.0D0 / XL
      PRES_ATM_OVER_XL = PRES_ATM * ONE_OVER_XL

      ACT1      = 1.0
      ACT2      = 1.0
      GM2       = 1.0
      HPWDEP    = 0.0

      RECIPAP1  = 1.0D0 / PRES_ATM
!     XC        = WCAVG * CTHK1
!     XC1       = 1.0 / XC
!     XC2       = RT / ( 1000.0D0 * CTHK1 )

      EWASHM    = 1.0 - EWASH
      DO_WETDEP = (EWASH < ONEM)

      XL_EWASH  = XL * EWASH
      XL_EWASHM = XL * EWASHM
      XC        = XL_EWASHM * RECIPAP1 * AIRM

C...set equilibrium constants as a function of temperature
C...   Henry`s law constants

      SO2H  = HLCONST( TEMP, .FALSE., 0.0, hlso2  )
      CO2H  = HLCONST( TEMP, .FALSE., 0.0, hlco2  )
      NH3H  = HLCONST( TEMP, .FALSE., 0.0, hlnh3  )
      H2O2H = HLCONST( TEMP, .FALSE., 0.0, hlh2o2 )
      O3H   = HLCONST( TEMP, .FALSE., 0.0, hlo3   )
      HCLH  = HLCONST( TEMP, .FALSE., 0.0, hlhcl  )
      HNO3H = HLCONST( TEMP, .FALSE., 0.0, hlhno3 )
      MHPH  = HLCONST( TEMP, .FALSE., 0.0, hlmeth )
      PAAH  = HLCONST( TEMP, .FALSE., 0.0, hlperx )
      FOAH  = HLCONST( TEMP, .FALSE., 0.0, hlfora )
      GLYH  = HLCONST( TEMP, .FALSE., 0.0, hlglyo )
      MGLYH = HLCONST( TEMP, .FALSE., 0.0, hlmgly )
      HOH   = HLCONST( TEMP, .FALSE., 0.0, hloh   )

      TEMP1 = 1.0 / TEMP - T0M

      DO I = 1, 15
         RATES(I) = RR0(I) * EXP( RR1(I) * TEMP1 )
      ENDDO

C...make initializations

      LOADING = 0.0D0
!     INITGAS = 0.0D0

C...compute fractional weights for several species

      TOTNIT = GAS( LHNO3 ) + AEROSOL( LNO3, ACC )
      IF ( TOTNIT .GT. 0.0D0 ) THEN
        FHNO3   = GAS( LHNO3 ) / TOTNIT
        FNO3ACC = AEROSOL( LNO3, ACC ) / TOTNIT
      ELSE
        FHNO3   = 1.0D0
        FNO3ACC = 0.0D0
      END IF

      TOTAMM = GAS( LNH3 ) + AEROSOL( LNH4, ACC )
      IF ( TOTAMM .GT. 0.0D0 ) THEN
        FNH3    = GAS( LNH3 ) / TOTAMM
        FNH4ACC = AEROSOL( LNH4, ACC ) / TOTAMM
      ELSE
        FNH3    = 1.0D0
        FNH4ACC = 0.0D0
      END IF

      TNO3a = AEROSOL( LNO3, ACC ) + AEROSOL( LNO3, COR )
      IF ( TNO3a .GT. 0.0D0) THEN
         FNO3COR = AEROSOL( LNO3, COR ) / TNO3a
      ELSE
         FNO3COR = 0.0D0
      END IF

      TNH4a = AEROSOL( LNH4, ACC ) + AEROSOL( LNH4, COR )
      IF ( TNH4a .GT. 0.0D0) THEN
         FNH4COR = AEROSOL( LNH4, COR ) / TNH4a
      ELSE
         FNH4COR = 0.0D0
      END IF

      TCLa  = AEROSOL( LCL, ACC )  + AEROSOL( LCL, COR )
      IF ( TCLa .GT. 0.0D0) THEN
         FCLCOR = AEROSOL( LCL, COR ) / TCLa
      ELSE
         FCLCOR = 0.0D0
      END IF

C...Assign fraction partitioning of FE(III) and MN(II)

      IF ( DARK ) THEN
        FE_III = 0.9D0  ! Night time, GS 01July2011
      ELSE
        FE_III = 0.1D0  ! Day time, GS 01July2011
      END IF

      MN_II = 1.0D0                     ! Same for day and night, GS  01July2011

C...Assign solubility of Fe and Mn

      FE_SOL = 0.1D0                   ! GS 01July2011
      MN_SOL = 0.5D0                   ! GS 28July2011

C...initial concentration from accumulation-mode aerosol loading (mol/liter)
C...  an assumption is made that all of the accumulation and coarse mode
C...  aerosol masses are incorporated into the cloud droplets

      DO IMODE = 2, 3
         DO ISPC = 1, NAER
            LOADING( ISPC, IMODE ) = AEROSOL( ISPC, IMODE ) * PRES_ATM_OVER_XL
         END DO
      END DO

!     LOADING( LSO4, ACC ) = ( AEROSOL( LSO4, ACC ) + GAS( LH2SO4 ) ) * PRES_ATM_OVER_XL

C...set constant factors that will be used in later multiplications (moles/atm)

      XLH2O2  = H2O2H * XL
      XLO3    = O3H   * XL
      XLMHP   = MHPH  * XL
      XLPAA   = PAAH  * XL
      XLSO2   = SO2H  * XL
      XLNH3   = NH3H  * XL
      XLHCL   = HCLH  * XL
      XLHNO3  = HNO3H * XL
      XLCO2   = CO2H  * XL

      SO212   = rates(3)  * rates(4)
      SO21H   = rates(3)  * SO2H
      SO212H  = SO212     * SO2H
      SO21A   = 1.D0 + XLSO2
      SO21B   = XLSO2 * rates(3)
      SO21C   = XLSO2 * SO212

!     CO212   = rates(5)  * rates(6) * xlco2
      CO21H   = rates(5)  * CO2H
      CO212H  = rates(6)  * CO21H  !CO212     * CO2H
      CO212   = CO212H    * XL
      CO212A  = RATES(5)  * XLCO2
      CO212B  = 1.D0 + XLCO2

      NH3DH20 = rates(8)  / rates(7)
      NH31HDH = NH3H      * NH3DH20
      NH31A   = NH3DH20 * XLNH3
      NH31B   = 1.D0 + XLNH3

      FOA1H   = rates(1)  * FOAH
      FOA1A   = 1.d0 + XL * FOAH
      FOA1B   = FOA1H * XL

      HNO31H  = rates(10) * HNO3H
      HNO3A   = 1.D0 + XLHNO3
      HNO3B   = rates(10) * XLHNO3

      HCL1H   = rates(9)  * HCLH
      HCLA    = 1.D0 + XLHCL
      HCLB    = rates(9) * XLHCL

      f2o2f  = 1.d0 / ( 1.0D0 + XLH2O2 )
      fo3f   = 1.d0 / ( 1.0D0 + XLO3 )
      fmhpf  = 1.d0 / ( 1.0D0 + XLMHP )
      fpaaf  = 1.d0 / ( 1.0D0 + XLPAA )
      fglyf  = 1.d0 / ( 1.0D0 + GLYH * XL )
      fmglyf = 1.d0 / ( 1.0D0 + MGLYH * XL )
      fhof   = 1.d0 / ( 1.0D0 + HOH * XL)

C...loop If kinetic calculations are made, return to this point

      DO ISPC = 1, NGAS
         INITGAS( ISPC )  = GAS( ISPC ) * PRES_ATM
      ENDDO

      INITGAS( LNH3 )  = INITGAS( LNH3 )
     &                 + ( LOADING( LNH4, ACC ) + LOADING( LNH4, COR ) ) * XL

      INITGAS( LHNO3 ) = INITGAS( LHNO3 ) + 2.0 * INITGAS( LN2O5 )
     &                 + ( LOADING( LNO3, ACC ) + LOADING( LNO3, COR ) ) * XL

      INITGAS( LHCL )  = INITGAS( LHCL )                 ! new for sea salt
     &                 + ( LOADING(  LCL, ACC ) + LOADING(  LCL, COR ) ) * XL

      TS6COR  = LOADING( LSO4,  COR )
      NO3COR  = LOADING( LNO3,  COR )
      CLCOR   = LOADING( LCL,   COR )
      NH4COR  = LOADING( LNH4,  COR )
      SOILCOR = LOADING( LSOILC,COR )
      ANTHCOR = LOADING( LANTHC,COR )
      SEASCOR = LOADING( LSEASC,COR )

      FECOR   =                         SOIL_FE_FAC * SOILCOR + CORS_FE_FAC * ANTHCOR
      MNCOR   =                         SOIL_MN_FAC * SOILCOR + CORS_MN_FAC * ANTHCOR
      NACOR   = SEAS_NA_FAC * SEASCOR + SOIL_NA_FAC * SOILCOR + CORS_NA_FAC * ANTHCOR
      MGCOR   = SEAS_MG_FAC * SEASCOR + SOIL_MG_FAC * SOILCOR + CORS_MG_FAC * ANTHCOR
      CACOR   = SEAS_CA_FAC * SEASCOR + SOIL_CA_FAC * SOILCOR + CORS_CA_FAC * ANTHCOR
      KCOR    = SEAS_K_FAC  * SEASCOR + SOIL_K_FAC  * SOILCOR + CORS_K_FAC  * ANTHCOR

      TS6     = LOADING( LSO4,   ACC ) + TS6COR
      NA      = LOADING( LNA,    ACC ) + NACOR
      CA      = LOADING( LCAACC, ACC ) + CACOR
      MG      = LOADING( LMGACC, ACC ) + MGCOR
      K       = LOADING( LKACC,  ACC ) + KCOR
      FE      = LOADING( LFEACC, ACC ) + FECOR
      MN      = LOADING( LMNACC, ACC ) + MNCOR
!     SOA     = LOADING( LSOA,   ACC )
      ORGC    = LOADING( LORGC,  ACC )
!     POA     = LOADING( LPOA,   ACC )
!     EC      = LOADING( LEC,    ACC )
!     PRIM    = LOADING( LPRI,   ACC )
!     NUMCOR  = LOADING( LNUM,   COR )

      BI0   = NA + K + 2.0D0 * (CA + MG)

      FE_OX = FE_III * FE_SOL * FE  ! GS 01July2011
      MN_OX = MN_II  * MN_SOL * MN  ! GS 01July2011

      A = 3.0D0 * FE_OX
      B = 2.0D0 * MN_OX

      ST0 = 4.0d0 * (CA + MG + MN_OX) + 9.0D0 * FE_OX + NA + K + A + B

      SK6TS6 = rates(2) * TS6

C...Aerosol specific to TXHG Versions

!     if (have_tox) then
!        TRACER     = LOADING( LTRACER_ACC, ACC )
!        TRACERCOR  = LOADING( LTRACER_COR, COR )
!        HGFINE     = LOADING( LPHG_ACC, ACC )
!        HGCOR      = LOADING( LPHG_COR, COR )
!     endif

C...  Find solution of the equation using a method of reiterative
C...  bisections Make initial guesses for pH:   between  1 and  10

      HA    = 2.0
      HPLUS = 1.E-2
      FA    = BISECT()

      HB    =  9.2
      HPLUS = 0.6309573445e-9
      FB    = BISECT()

      IF (FA * FB > 0.0) THEN

         write(*,*) "PH NOT BETWEEN 2 AND 9.2 AT ", itab_w(iw)%iwglobe, kw
         write(*,*) "RESETTING BOUNDS TO 1 AND 10"

         HA    = 1.0
         HPLUS = 1.E-1
         FA    = BISECT()

         HB    = 10.0
         HPLUS = 1.E-10
         FB    = BISECT()

         IF (FA * FB > 0.0) THEN

            XMSG = 'PH NOT BETWEEN 1 AND 10'
            CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT2 )

         ENDIF

      ENDIF

      haa = ha
      faa = fa

      hba = hb
      fbb = fb

C...  Start iteration and bisection ****************

      DO I30C = 1, 15

         IF (I30C == 1) THEN
            BB    = 6.5
            HPLUS = .31622777E-6
         ELSE
            BB1 = (FB * HA - FA * HB) / (FB - FA)
            BB2 = 0.5 * ( HA + HB )
            BB  = 0.85 * BB1 + 0.15 * BB2
            HPLUS = EXP( MLN10 * BB )
         ENDIF

         FN = bisect()

C... Check convergence of solutions

!        if (abs(fn) < 1.e-20) exit
         if (I30C>1 .and. abs(bb-hb) <= .01) exit

         if ( fn * fb < 0.0 ) then

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            ha = hb
            fa = fb
            hb = bb
            fb = fn

         else

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            hb = bb
            fb = fn
            fa = 0.5 * fa

         endIf

C... Compute Ionic strength and activity coefficient by the Davies equation

         STION = 0.5D0 * ( ST0 + HPLUS + NH4 + OH + HCO3 + HSO3
     &         + 4.0D0 * ( SO4 + CO3 + SO3 ) + NO3 + HSO4 + CL + HCO2 ) ! KMF 08September2011

         GM1LOG = STION * ( 0.1018 - 0.509 / (STION + SQRT(STION)) )

         GM1  = EXP( LN10 * GM1LOG )
         ACT1 = MAX( GM1  * GM1,  1.0E-30 )
         GM2  = MAX( ACT1 * ACT1, 1.0E-30 )
         ACT2 = MAX( ACT1 * GM2,  1.0E-30 )


      END DO  ! I30C

C...  Make those concentration calculations which can be made outside
C...  of the function.

!     phclf  = phclfoa1   * a1
!     phno3f = phno3foa1  * a1
!     pfoaf  = pfoafoa1   * a1
      pso2f  = pso2foa1a2 * a1a2
!     pco2f  = pco2foa1a2 * a1a2

      SO2L = SO2H * PSO2F
      SIV  = SO3 + HSO3 + SO2L

C...Calculate final gas phase concentrations of oxidants (atm)

      PH2O2F = INITGAS( LH2O2 ) * f2o2f
      PO3F   = INITGAS( LO3 )   * fo3f
      PMHPF  = INITGAS( LMHP )  * fmhpf
      PPAAF  = INITGAS( LPAA )  * fpaaf
      PGLYF  = INITGAS( LGLY )  * fglyf
      PMGLYF = INITGAS( LMGLY ) * fmglyf
      PHOF   = INITGAS( LHO )   * fhof

C...Calculate liquid phase concentrations of oxidants (moles/liter)

      H2O2L = PH2O2F * H2O2H
      O3L   = PO3F   * O3H
      MHPL  = PMHPF  * MHPH
      PAAL  = PPAAF  * PAAH
!     FOAL  = PFOAF  * FOAH
!     NH3L  = PNH3F  * NH3H
!     CO2L  = PCO2F  * CO2H
!     HCLL  = PHCLF  * HCLH
!     HNO3L = PHNO3F * HNO3H
      GLYL  = PGLYF  * GLYH
      MGLYL = PMGLYF * MGLYH
      OHL   = PHOF   * HOH

C...compute modal concentrations

!     SO4COR  = rates(2) * TS6COR / ( HPLUS * GM2 + rates(2) )
!     HSO4COR = MAX( TS6COR - SO4COR, 0.0D0 )

!     TS6ACC  = MAX( TS6  - TS6COR,   0.0D0 )
!     SO4ACC  = MAX( SO4  - SO4COR,   0.0D0 )
!     HSO4ACC = MAX( HSO4 - HSO4COR,  0.0D0 )
!     NO3ACC  = MAX( NO3  - NO3COR,   0.0D0 )
!     NAACC   = MAX( NA   - NACOR,    0.0D0 )
!     CLACC   = MAX( CL   - CLCOR,    0.0D0 )
!     NH4ACC  = MAX( NH4  - NH4COR,   0.0D0 )
!     CAACC   = MAX( CA   - CACOR,    0.0D0 ) ! AE6
!     MGACC   = MAX( MG   - MGCOR,    0.0D0 ) ! AE6
!     KACC    = MAX( K    - KCOR,     0.0D0 ) ! AE6
!     FEACC   = MAX( FE   - FECOR,    0.0D0 ) ! AE6
!     MNACC   = MAX( MN   - MNCOR,    0.0D0 ) ! AE6

C...Define the total S(iv) available for oxidation

      TSIV = INITGAS( LSO2 ) * ONE_OVER_XL

C...Calculate sulfur iv oxidation rate due to H2O2 (Jacobson, 1997)

      DSIVDT( 1 ) = -rates(11) * H2O2L * HSO3 * HPLUS / ( 1.0D0 + 13.0D0 * HPLUS )
      SSIV  ( 1 ) = MIN( INITGAS( LH2O2 ) * ONE_OVER_XL, TSIV )

C...Calculate sulfur iv oxidation rate due to O3 (Jacobson, 1997)

      DSIVDT( 2 ) = -O3L * ( 2.40D04   * SO2L
     &                     + rates(14) * HSO3
     &                     + rates(15) *  SO3 )
      SSIV  ( 2 ) = MIN( INITGAS( LO3 ) * ONE_OVER_XL, TSIV )

C...Calculate sulfur iv oxidation rate due to 02 catalyzed by Mn++ and Fe+++
C...(Martin and Goodman, 1991)
C...Implement sulfate inhibition based on Martin and Good, 1991

      DSIVDT( 3 ) = - SIV * (  750.0D0  * MN_OX + ! GS 4May2011
     &                        2600.0D0  * FE_OX +                 ! GS 4May2011
     &                           1.0D10 * MN_OX * FE_OX )         ! GS 4May2011
     &              / ( 1.0 + 75.0 * (real(TS6)**onethird)**2 )  ! GS 28Nov2012
      SSIV  ( 3 ) = TSIV

C...Calculate sulfur oxidation rate due to MHP (Jacobson,  1997)

      DSIVDT( 4 ) = -rates(12) * HPLUS * MHPL * HSO3
      SSIV  ( 4 ) = MIN( INITGAS( LMHP ) * ONE_OVER_XL, TSIV )

C...Calculate sulfur oxidation due to PAA (Jacobson,  1997)

      DSIVDT( 5 ) = -( rates(13) * HPLUS + 7.00D2 ) * HSO3 * PAAL
      SSIV  ( 5 ) = MIN( INITGAS( LPAA ) * ONE_OVER_XL, TSIV )

C...CALCULATE IN-CLOUD SOA PRODUCTION
C...  Reference:  Carlton, A.G., B.J. Turpin, K.E. Altieri, A. Reff,
C...  S. Seitzinger, H.J. Lim, and B. Ervens (2007), Atmospheric Oxalic
C...  Acid and SOA Production from Glyoxal: Results of Aqueous
C...  Photooxidation Experiments, Atmos. Environ., 41(35), 7588-7602.

C...Calculate GLY oxidation due to OH

      DSIVDT( 6 ) = -RGLY3 * GLYL * OHL
      SSIV  ( 6 ) = INITGAS( LGLY ) * ONE_OVER_XL

C...Calculate MGLY oxidation due to OH

      DSIVDT( 7 ) = -RMGLY3 * MGLYL * OHL
      SSIV  ( 7 ) = INITGAS( LMGLY ) * ONE_OVER_XL

C...Calculate oxidation rates

      DO I = 1, NUMOX + 2
         EFACT = EXP( TAUCLD * DSIVDT(I) / REAL( MAX(SSIV(I), 1.D-30) ) )
         DS4(I) = SSIV(I) * ( REAL(EFACT,8) - 1.D0 )
      ENDDO

C...Calculate total sulfur iv oxidation rate

      DS4(0) = SUM( DS4( 1:NUMOX ) )

C...    If the total S(iv) oxidation rate over the timestep exceeds the amount
C...    of S(iv) available then scale the rates to conserve mass

      IF ( -DS4( 0 ) .GT. TSIV ) THEN
         DSIV_SCALE = -TSIV / DS4(0)
         DO I = 0, NUMOX
            DS4(I) = DS4(I) * DSIV_SCALE
         ENDDO
      ENDIF

      so2_prev = INITGAS( LSO2 )
      ts6_prev = ts6

      INITGAS( LSO2 )  = MAX( INITGAS( LSO2 )  + XL * DS4( 0 ), 0.D0 )
      INITGAS( LH2O2 ) = MAX( INITGAS( LH2O2 ) + XL * DS4( 1 ), 0.D0 )
      INITGAS( LO3 )   = MAX( INITGAS( LO3 )   + XL * DS4( 2 ), 0.D0 )
      INITGAS( LMHP )  = MAX( INITGAS( LMHP )  + XL * DS4( 4 ), 0.D0 )
      INITGAS( LPAA )  = MAX( INITGAS( LPAA )  + XL * DS4( 5 ), 0.D0 )
      INITGAS( LGLY )  = MAX( INITGAS( LGLY )  + XL * DS4( 6 ), 0.D0 )
      INITGAS( LMGLY ) = MAX( INITGAS( LMGLY ) + XL * DS4( 7 ), 0.D0 )

      TS6 = TS6 - DS4( 0 )

      SK6TS6 = rates(2) * TS6

C...Calculate AORGC Production:  4% SOAcld (ORGC) yield from glyoxal
C...  and methylglyoxal is assumed

      ORGC = ORGC - 0.04D0 * ( DS4(6) + DS4(7) )

C...  Redo iteration and bisection ****************

      HA    = haa
      FA    = faa

      HB    = hba
      FB    = fbb

      DO I30C = 1, 15

         IF (I30C > 1) THEN
            BB1 = (FB * HA - FA * HB) / (FB - FA)
            BB2 = 0.5 * ( HA + HB )
            BB  = 0.85 * BB1 + 0.15 * BB2
            HPLUS = EXP( MLN10 * BB )
         ENDIF

         FN = bisect()

C... Check convergence of solutions

!        if (abs(fn) < 1.e-20) exit
         if (I30C>1 .and. abs(bb-hb) <= .01) exit

         if ( fn * fb < 0.0 ) then

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            ha = hb
            fa = fb
            hb = bb
            fb = fn

         else

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            hb = bb
            fb = fn
            fa = 0.5 * fa

         endIf

C... Compute Ionic strength and activity coefficient by the Davies equation

         STION = 0.5D0 * ( ST0 + HPLUS + NH4 + OH + HCO3 + HSO3
     &         + 4.0D0 * ( SO4 + CO3 + SO3 ) + NO3 + HSO4 + CL + HCO2 ) ! KMF 08September2011

         GM1LOG = STION * ( 0.1018 - 0.509 / (STION + SQRT(STION)) )

         GM1  = EXP( LN10 * GM1LOG )
         ACT1 = MAX( GM1  * GM1,  1.0E-30 )
         GM2  = MAX( ACT1 * ACT1, 1.0E-30 )
         ACT2 = MAX( ACT1 * GM2,  1.0E-30 )

      END DO  ! I30C

C...  Make those concentration calculations which can be made outside
C...  of the function.

      phclf  = phclfoa1   * a1
      phno3f = phno3foa1  * a1
      pfoaf  = pfoafoa1   * a1
      pso2f  = pso2foa1a2 * a1a2
      pco2f  = pco2foa1a2 * a1a2

      SO2L = SO2H * PSO2F
      SIV  = SO3 + HSO3 + SO2L

C...Calculate final gas phase concentrations of oxidants (atm)

      PH2O2F = INITGAS( LH2O2 ) * f2o2f
      PO3F   = INITGAS( LO3 )   * fo3f
      PMHPF  = INITGAS( LMHP )  * fmhpf
      PPAAF  = INITGAS( LPAA )  * fpaaf
      PGLYF  = INITGAS( LGLY )  * fglyf
      PMGLYF = INITGAS( LMGLY ) * fmglyf
!     PHOF   = INITGAS( LHO )   * fhof

C...Calculate liquid phase concentrations of oxidants (moles/liter)

      H2O2L = PH2O2F * H2O2H
      O3L   = PO3F   * O3H
      MHPL  = PMHPF  * MHPH
      PAAL  = PPAAF  * PAAH
      FOAL  = PFOAF  * FOAH
      NH3L  = PNH3F  * NH3H
      CO2L  = PCO2F  * CO2H
      HCLL  = PHCLF  * HCLH
      HNO3L = PHNO3F * HNO3H
      GLYL  = PGLYF  * GLYH
      MGLYL = PMGLYF * MGLYH
!     OHL   = PHOF   * HOH

C...compute modal concentrations

!     SO4COR  = rates(2) * TS6COR / ( HPLUS * GM2 + rates(2) )
!     HSO4COR = MAX( TS6COR - SO4COR, 0.0D0 )

      TS6ACC  = MAX( TS6  - TS6COR,   0.0D0 )
!     SO4ACC  = MAX( SO4  - SO4COR,   0.0D0 )
!     HSO4ACC = MAX( HSO4 - HSO4COR,  0.0D0 )

!     NAACC   = MAX( NA   - NACOR,    0.0D0 )
!     CAACC   = MAX( CA   - CACOR,    0.0D0 ) ! AE6
!     MGACC   = MAX( MG   - MGCOR,    0.0D0 ) ! AE6
!     KACC    = MAX( K    - KCOR,     0.0D0 ) ! AE6
!     FEACC   = MAX( FE   - FECOR,    0.0D0 ) ! AE6
!     MNACC   = MAX( MN   - MNCOR,    0.0D0 ) ! AE6

C...Avoid adding mass when the coarse mode concentration is greater
C...  than the total amount left in the aqueous phase after redistribution
C...  of a species between the gas/aqueous phases

      IF ( NO3COR .GT. NO3 ) then
         NO3ACC = (1.0D0 - FNO3COR) * NO3
         NO3COR = FNO3COR * NO3
      ELSE
         NO3ACC  = MAX( NO3 - NO3COR,   0.0D0 )
      END IF

      IF ( CLCOR .GT. CL ) then
         CLACC = (1.0D0 - FCLCOR) * CL
         CLCOR = FCLCOR * CL
      ELSE
         CLACC   = MAX( CL - CLCOR,    0.0D0 )
      END IF

      IF ( NH4COR .GT. NH4 ) THEN
         NH4ACC = (1.0D0 - FNH4COR) * NH4
         NH4COR = FNH4COR * NH4
      ELSE
         NH4ACC  = MAX( NH4 - NH4COR,   0.0D0 )
      END IF

C...  compute the scavenging coefficient for SO4 which will be used for
C...  scavenging aerosol number in the accumulation mode

!     EBETASO4T = 1.0
!     EBETASO4T = EWASH

!     IF (DO_WETDEP) THEN
!        SO4AVAIL = LOADING( LSO4, ACC ) - DS4( 0 )
!        IF (SO4AVAIL > 1.E-30) THEN
!           DEPSUM = ( WETDEP(LSO4ACCL) + WETDEP(LHSO4ACCL) ) * XC1
!           EBETASO4T = EXP( -DEPSUM / SO4AVAIL )
!        ENDIF
!     ENDIF

C...Compute the output concentrations and wet deposition amounts

      TOTAMM = ( PNH3F  + ( NH4ACC + NH3L  ) * XL ) * RECIPAP1
      TOTNIT = ( PHNO3F + ( NO3ACC + HNO3L ) * XL ) * RECIPAP1

C...gas concentrations (mol/molV)

      tfoa = foal + hco2
      tco2 = co2l + hco3 + co3

      GAS( LSO2   ) = ( PSO2F  + XL_EWASH * SIV  ) * RECIPAP1
      GAS( LH2O2  ) = ( PH2O2F + XL_EWASH * H2O2L) * RECIPAP1
      GAS( LO3    ) = ( PO3F   + XL_EWASH * O3L  ) * RECIPAP1
      GAS( LCO2   ) = ( PCO2F  + XL_EWASH * TCO2 ) * RECIPAP1
      GAS( LFOA   ) = ( PFOAF  + XL_EWASH * TFOA ) * RECIPAP1
      GAS( LMHP   ) = ( PMHPF  + XL_EWASH * MHPL ) * RECIPAP1
      GAS( LPAA   ) = ( PPAAF  + XL_EWASH * PAAL ) * RECIPAP1
      GAS( LHCL   ) = ( PHCLF  + XL_EWASH * HCLL ) * RECIPAP1
      GAS( LGLY   ) = ( PGLYF  + XL_EWASH * GLYL ) * RECIPAP1
      GAS( LMGLY  ) = ( PMGLYF + XL_EWASH * MGLYL) * RECIPAP1
!     GAS( LHO    ) = ( PHOF   + XL_EWASH * OHL  ) * RECIPAP1

      GAS( LNH3   ) = FNH3  * TOTAMM -  NH3L * XL_EWASHM * RECIPAP1
      GAS( LHNO3  ) = FHNO3 * TOTNIT - HNO3L * XL_EWASHM * RECIPAP1
      GAS( LN2O5  ) = 0.0D0 ! assume all into aerosol
!     GAS( LH2SO4 ) = 0.0D0 ! assume all into aerosol

      HP = HPLUS

C...  Gas-phase species wet deposition (mm mol/lit)

      IF (DO_WETDEP) THEN
         gaswdep( lso2  ) = xc * siv
         gaswdep( lnh3  ) = xc * nh3l
         gaswdep( lh2o2 ) = xc * h2o2l
         gaswdep( lo3   ) = xc * o3l
         gaswdep( lco2  ) = xc * tco2
         gaswdep( lfoa  ) = xc * tfoa
         gaswdep( lmhp  ) = xc * mhpl
         gaswdep( lpaa  ) = xc * paal
         gaswdep( lhcl  ) = xc * hcll
         gaswdep( lhno3 ) = xc * hno3l
         gaswdep( lgly  ) = xc * glyl
         gaswdep( lmgly ) = xc * mglyl
!        gaswdep( ln2o5 ) = 0.0d0
!        gaswdep( lh2so4) = 0.0d0
!        gaswdep( lho   ) = xc * lohl

C...store the amount of hydrogen deposition
         HPWDEP = HPLUS * XC
      ENDIF

!         GASWDEP( LSO2   ) = WETDEP( LSO3L  ) + WETDEP( LHSO3L )
!     &                     + WETDEP( LSO2L  )
!         GASWDEP( LNH3   ) = WETDEP( LNH3L  )
!         GASWDEP( LH2O2  ) = WETDEP( LH2O2L )
!         GASWDEP( LO3    ) = WETDEP( LO3L   )
!         GASWDEP( LCO2   ) = WETDEP( LCO3L  ) + WETDEP( LHCO3L )
!     &                     + WETDEP( LCO2L  )
!         GASWDEP( LFOA   ) = WETDEP( LFOAL  ) + WETDEP( LHCO2L )
!         GASWDEP( LMHP   ) = WETDEP( LMHPL  )
!         GASWDEP( LPAA   ) = WETDEP( LPAAL  )
!         GASWDEP( LHCL   ) = WETDEP( LHCLL  )
!         GASWDEP( LHNO3  ) = WETDEP( LHNO3L )
!!        GASWDEP( LN2O5  ) = 0.0D0
!!        GASWDEP( LH2SO4 ) = 0.0D0
!         GASWDEP( LGLY   ) = WETDEP( LGLYL  )
!         GASWDEP( LMGLY  ) = WETDEP( LMGLYL )
!!        GASWDEP( LHO    ) = WETDEP( LOHL   )
!
!C...  Aerosol species wet deposition (mm mol/lit)
!C...  there is no wet deposition of aitken particles, they attached
!C...  to the accumulation mode particles
!
!         AERWDEP( LSO4, ACC ) = WETDEP( LSO4ACCL ) + WETDEP( LHSO4ACCL )
!         AERWDEP( LNH4, ACC ) = WETDEP( LNH4ACCL )
!         AERWDEP( LNO3, ACC ) = WETDEP( LNO3ACCL )
!         AERWDEP( LEC,  ACC ) = WETDEP( LECL     )
!         AERWDEP( LPRI, ACC ) = WETDEP( LPRIML   )
!
!         AERWDEP( LSOA, ACC ) = WETDEP( LSOAL )
!         AERWDEP( LORGC,ACC ) = WETDEP( LORGCL)
!         AERWDEP( LPOA, ACC ) = WETDEP( LPOAL )
!
!         AERWDEP( LSO4, COR ) = WETDEP( LTS6CORL  )
!         AERWDEP( LNO3, COR ) = WETDEP( LNO3CORL  )
!         AERWDEP( LNH4, COR ) = WETDEP( LNH4CORL  )
!!        AERWDEP( LPRICOR, COR ) = WETDEP( LPRIMCORL )
!
!         AERWDEP( LNA, ACC  ) = WETDEP( LNAACCL )
!         AERWDEP( LCL, ACC  ) = WETDEP( LCLACCL )
!!        AERWDEP( LNA, COR  ) = WETDEP( LNACORL )
!         AERWDEP( LCL, COR  ) = WETDEP( LCLCORL )
!
!!        AERWDEP( LK,     COR ) = WETDEP( LKCORL  )
!!        AERWDEP( LA3FE,  COR ) = WETDEP( LFECORL )
!!        AERWDEP( LB2MN,  COR ) = WETDEP( LMNCORL )
!!        AERWDEP( LCACO3, COR ) = WETDEP( LCACORL )
!!        AERWDEP( LMGCO3, COR ) = WETDEP( LMGCORL )
!
!         AERWDEP( LCAACC, ACC ) = WETDEP( LCAACCL ) ! AE6 - SLN 16March2011
!         AERWDEP( LMGACC, ACC ) = WETDEP( LMGACCL ) ! AE6 - SLN 16March2011
!         AERWDEP( LKACC,  ACC ) = WETDEP( LKACCL  ) ! AE6 - SLN 16March2011
!         AERWDEP( LSOILC, COR ) = WETDEP( LSOILCL ) ! AE6 - SLN 16March2011
!         AERWDEP( LANTHC, COR ) = WETDEP( LANTHCL ) ! AE6 - SLN 16March2011
!         AERWDEP( LSEASC, COR ) = WETDEP( LSEASCL ) ! AE6 - SLN 16March2011
!         AERWDEP( LFEACC, ACC ) = WETDEP( LFEACCL ) ! AE6 - SLN 22March2011
!         AERWDEP( LMNACC, ACC ) = WETDEP( LMNACCL ) ! AE6 - SLN 22March2011

!        AERWDEP( LNUM, ACC ) = 0.0D0
!        AERWDEP( LNUM, COR ) = 0.0D0
!        AERWDEP( LSRF, ACC ) = 0.0D0
!        AERWDEP( LSRF, COR ) = 0.0D0

C...Compute for Aerosol species in  TXHG Version

!         if (have_tox) then
!            AERWDEP( LTRACER_ACC, ACC ) = WETDEP( LTRACERL )
!            AERWDEP( LTRACER_COR, COR ) = WETDEP( LTRACERCORL )
!            AERWDEP( LPHG_ACC,    ACC ) = WETDEP( LPHGFINEL )
!            AERWDEP( LPHG_COR,    COR ) = WETDEP( LPHGCORL )
!         endif

!        AERWDEP( LH2O, ACC) = AEROSOL( LH2O, ACC ) * PRES_ATM_OVER_XL * XC * (1.0 - EBETASO4T)
!        AERWDEP( LH2O, COR) = AEROSOL( LH2O, COR ) * PRES_ATM_OVER_XL * XC * EWASHM

C...store the amount of hydrogen deposition

!         HPWDEP = HPLUS

!      ENDIF

C...aerosol concentrations (mol/molV)

      AEROSOL( LSO4,   ACC ) = TS6ACC  * XL * RECIPAP1
      AEROSOL( LORGC,  ACC ) = ORGC    * XL * RECIPAP1
      AEROSOL( LNH4,   ACC ) = FNH4ACC * TOTAMM
      AEROSOL( LNO3,   ACC ) = FNO3ACC * TOTNIT
      AEROSOL( LCL,    ACC ) = CLACC   * XL * RECIPAP1
!     AEROSOL( LNA,    ACC ) = NAACC   * XL * RECIPAP1

!     AEROSOL( LCAACC, ACC ) = CAACC   * XL * RECIPAP1 ! AE6 - SLN 16March2011
!     AEROSOL( LMGACC, ACC ) = MGACC   * XL * RECIPAP1 ! AE6 - SLN 16March2011
!     AEROSOL( LKACC,  ACC ) = KACC    * XL * RECIPAP1 ! AE6 - SLN 16March2011
!     AEROSOL( LFEACC, ACC ) = FEACC   * XL * RECIPAP1 ! AE6 - SLN 22March2011
!     AEROSOL( LMNACC, ACC ) = MNACC   * XL * RECIPAP1 ! AE6 - SLN 22March2011

!     AEROSOL( LEC,    ACC ) = EC      * XL * RECIPAP1  ! entirely unused!!!
!     AEROSOL( LPRI,   ACC ) = PRIM    * XL * RECIPAP1  ! entirely unused!!!
!     AEROSOL( LSOA,   ACC ) = SOA     * XL * RECIPAP1  ! entirely unused!!!
!     AEROSOL( LPOA,   ACC ) = POA     * XL * RECIPAP1  ! entirely unused!!!

!      AEROSOL( LH2O,   ACC ) = AEROSOL( LH2O, ACC ) * EBETASO4T
!      AEROSOL( LNUM,   ACC ) = AEROSOL( LNUM, ACC ) * EBETASO4T

      AEROSOL( LNO3,    COR ) = NO3COR  * XL * RECIPAP1
      AEROSOL( LCL,     COR ) = CLCOR   * XL * RECIPAP1
      AEROSOL( LNH4,    COR ) = NH4COR  * XL * RECIPAP1

      ! these rmeain unchanged
!     AEROSOL( LSO4,    COR ) = TS6COR  * XL * RECIPAP1
!     AEROSOL( LSOILC,  COR ) = SOILCOR * XL * RECIPAP1 ! AE6 - SLN 16March2011
!     AEROSOL( LANTHC,  COR ) = ANTHCOR * XL * RECIPAP1 ! AE6 - SLN 16March2011
!     AEROSOL( LSEASC,  COR ) = SEASCOR * XL * RECIPAP1 ! AE6 - SLN 16March2011

!     AEROSOL( LPRICOR, COR ) = PRIMCOR * XL * RECIPAP1
!     AEROSOL( LK,      COR ) = KCOR    * XL * RECIPAP1
!     AEROSOL( LA3FE,   COR ) = FECOR   * XL * RECIPAP1
!     AEROSOL( LB2MN,   COR ) = MNCOR   * XL * RECIPAP1
!     AEROSOL( LCACO3,  COR ) = CACOR   * XL * RECIPAP1
!     AEROSOL( LMGCO3,  COR ) = MGCOR   * XL * RECIPAP1
!     AEROSOL( LNA,     COR ) = NACOR   * XL * RECIPAP1

!     AEROSOL( LH2O,    COR ) = AEROSOL( LH2O, COR ) * EWASH
!     AEROSOL( LNUM,    COR ) = NUMCOR  * XL * RECIPAP1

C...Compute for Aerosol species in  TXHG Version

!     if (have_tox) then
!        AEROSOL( LTRACER_ACC, ACC ) = TRACER    * XL * RECIPAP1
!        AEROSOL( LTRACER_COR, COR ) = TRACERCOR * XL * RECIPAP1
!        AEROSOL( LPHG_ACC,    ACC ) = HGFINE * XL * RECIPAP1
!        AEROSOL( LPHG_COR,    COR ) = HGCOR * XL * RECIPAP1
!     endif

      CONTAINS


      REAL FUNCTION BISECT()
      IMPLICIT NONE

!     HPLUS = EXP( MLN10 * PH )

      A1   = HPLUS * ACT1
      A2   = HPLUS * HPLUS * ACT2
      a1a2 = a1 * a2

C...calculate final gas phase partial pressure of SO2, NH3, HCL, HNO3
C...  HCOOH, and CO2 (atm)

!      PSO2Foa1a2 = INITGAS( LSO2 ) / ( a1a2 + XLSO2
!     &     * ( a1a2 + rates(3) * a2 + SO212 * a1 ) )
      PSO2Foa1a2 = INITGAS( LSO2 ) / ( a1a2 * SO21A + SO21B * a2 + SO21C * a1 )

!     PNH3F = INITGAS( LNH3 ) / ( 1.0D0 + XLNH3 * ( 1.0D0 + NH3DH20 * HPLUS ) )
      PNH3F = INITGAS( LNH3 ) / ( NH31A * HPLUS + NH31B )

!     PHCLFoA1 = INITGAS( LHCL ) / ( A1 + XLHCL * ( A1 + rates(9) ) )
      PHCLFoA1 = INITGAS( LHCL ) / ( A1 * HCLA + HCLB )

!     PHNO3FoA1 = INITGAS( LHNO3 ) / ( A1 + XLHNO3 * ( A1 + rates(10) ) )
      PHNO3FoA1 = INITGAS( LHNO3 ) / ( A1 * HNO3A + HNO3B )

!     PFOAFoA1 = INITGAS( LFOA ) / ( A1 + XL * ( FOAH * A1 + FOA1H ) )
      PFOAFoA1 = INITGAS( LFOA ) / ( A1 * FOA1A  + FOA1B )

      PCO2FOA1A2 = INITGAS( LCO2 ) / ( A1A2 * CO212B + CO212A * A2 + CO212 * A1 )

C...calculate liquid phase concentrations (moles/liter)

      SO4  = SK6TS6 / ( HPLUS * GM2 + rates(2) )
      HSO4 = TS6 - SO4
      SO3  = SO212H  * PSO2FoA1A2 * A1
      HSO3 = SO21H   * PSO2FoA1A2 * A2
      CO3  = CO212H  * PCO2FoA1A2 * A1
      HCO3 = CO21H   * PCO2FoA1A2 * A2
      OH   = rates(7) / A1
      NH4  = NH31HDH * PNH3F  * HPLUS
      HCO2 = FOA1H   * PFOAFoA1
      NO3  = HNO31H  * PHNO3FoA1
      CL   = HCL1H   * PHCLFoA1 ! new for sea salt

C...compute functional value
!           FB = HPLUS + NH4 + NA + 2.0D0 * ( CA + MG - CO3 - SO3 - SO4 )
!    &         - OH - HCO3 - HSO3 - NO3 - HSO4 - HCO2 - CL

!     BISECT = HPLUS + NH4 + NA + K + 2.0D0 * ( CA + MG - CO3 - SO3 - SO4 ) ! SLN 16March2011
!    &       - OH - HCO3 - HSO3 - NO3 - HSO4 - HCO2 - CL

      BISECT = HPLUS + NH4 + BI0 - 2.0D0 * ( CO3 + SO3 + SO4 ) ! SLN 16March2011
     &       - OH - HCO3 - HSO3 - NO3 - HSO4 - HCO2 - CL

      END FUNCTION BISECT

      END SUBROUTINE AQCHEM
