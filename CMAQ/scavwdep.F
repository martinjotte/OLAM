
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
C $Header: /project/yoj/arc/CCTM/src/cloud/cloud_acm_ae5/scavwdep.F,v 1.7 2011/10/21 16:10:28 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE SCAVWDEP ( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                      CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                      REMOV, REMOVAC, ALFAC, RHOAIR )
C-----------------------------------------------------------------------
C  Description: Compute simplistic incloud scavenging and wet removal
 
C  Revision History:
C      No   Date   Who  What
C      -- -------- ---  -----------------------------------------
C       0 01/15/98 sjr  created program
C       1 07/27/99 sjr  added scavenging coefficient scaling factors
C       3 Dec 00   Jeff move CGRID_MAP into f90 module
C       4 Dec 02   sjr  revised calls to HLCONST
C       5 Jun 10   J.Young: convert for Namelist redesign
C       6 Mar 11   S.Roselle: replaced I/O API include files with UTILIO_DEFN;
C                       removed deprecated TRIMLEN
C       6 Mar 11   B.Hutzell: added TRACER string for AE_SPCs to omit while
C                       computing Aitken mode mass
 
C  Called by:  RADMCLD and RESCLD
 
C  Calls the following subroutines:  GETALPHA
 
C  Calls the following functions:  HLCONST
 
C  Arguments    Type      I/O       Description
C  ---------   -------  ------  --------------------------------
C    JDATE     integer   input  current model julian date (yyyyddd)
C    JTIME     integer   input  current model time (hhmmss)
C    WTBAR      real     input  avg total water content (kg/m3)
C    WCBAR      real     input  avg liquid water content (kg/m3)
C    TBARC      real     input  avg cloud temperature (K)
C    PBARC      real     input  avg cloud pressure (Pa)
C    CTHK1      real     input  cloud thickness (m)
C    AIRM       real     input  total air mass (moles/m2) in cloudy air
C    PRATE1     real     input  precipitation rate (mm/hr)
C    TAUCLD     real     input  cloud lifetime (s)
C    POLC       real     input  ave vert conc incloud (moles sp/ mole air)
C    CEND       real    output  ending incloud conc (moles/mole)
C    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
C    REMOVAC    real    output  variable storing H+ deposition
C    ALFA0      real    output  scav coef for aitken aerosol number
C    ALFA2      real    output  scav coef for aitken mode sfc area
C    ALFA3      real    output  scav coef for aitken aerosol mass
 
C-----------------------------------------------------------------------

      USE CGRID_SPCS                    ! CGRID mechanism species
      USE UTILIO_DEFN
      use hlconst_mod
      use misc_coms, only: io6

      IMPLICIT NONE

C...........Includes:

      INCLUDE 'CONST.EXT'               ! constants

      CHARACTER( 120 ) :: XMSG = ' '    ! Exit status message

C...........Parameters:

      REAL, PARAMETER :: H2ODENS = 1000.0  ! density of water at 20 C
                                           ! and 1 ATM (kg/m3)

C Number of species in CGRID
      INTEGER, SAVE :: MXSPCS

C Number of species scavenged
      INTEGER, SAVE :: N_CGRID_SCAV

      REAL, PARAMETER :: TWOTHIRDS = 2.0 / 3.0

      REAL, PARAMETER :: KGPG = 1.0E-03   ! kilograms per gram

C...........Arguments:

      INTEGER, INTENT(  IN ) :: JDATE   ! current model date, coded YYYYDDD
      INTEGER, INTENT(  IN ) :: JTIME   ! current model time, coded HHMMSS
      REAL,    INTENT(  IN ) :: WTBAR   ! total wat cont (kg/m2) int. thru cloud depth
      REAL,    INTENT(  IN ) :: WCBAR   ! liq water content of cloud (kg/m3)
      REAL,    INTENT(  IN ) :: TBARC   ! mean cloud temp (K)
      REAL,    INTENT(  IN ) :: PBARC   ! mean cloud pressure (Pa)
      REAL,    INTENT(  IN ) :: CTHK1   ! aq chem calc cloud thickness
      REAL,    INTENT(  IN ) :: AIRM    ! total air mass (moles/m2) in cloudy air
      REAL,    INTENT(  IN ) :: PRATE1  ! storm rainfall rate (mm/hr)
      REAL,    INTENT(  IN ) :: TAUCLD  ! cloud lifetime
      REAL,    INTENT( OUT ) :: ALFAC   ! scav coef for aerosol
      REAL,    INTENT( OUT ) :: REMOVAC ! variable storing H+ deposition
      REAL,    INTENT(  IN ) :: POLC ( : )  ! avg vert conc incloud (moles/mole)
      REAL,    INTENT( OUT ) :: CEND ( : )  ! ending incloud conc (moles/mole)
      REAL,    INTENT( OUT ) :: REMOV( : )  ! moles/m2 or mm*mol/lit scavenged
      REAL,    INTENT(  IN ) :: RHOAIR      ! air density in kg/m3

C...........Local Variables:

      LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! flag for first pass thru
      CHARACTER( 16 ), SAVE :: PNAME = 'SCAVWDEP'  ! program name

      CHARACTER( 16 ), ALLOCATABLE, SAVE :: CGRID_SCAV( : )  ! CGRID species scavenged
      INTEGER, ALLOCATABLE, SAVE :: CGRID_SCAV_MAP( : ) ! CGRID map to scavenged spc
      REAL, ALLOCATABLE,    SAVE :: CGRID_SCAV_FAC( : )  ! CGRID scav coef factors
      INTEGER, ALLOCATABLE, SAVE :: L_NUMAKN( : ) ! pointers to aitken aerosol #
      INTEGER, ALLOCATABLE, SAVE :: L_MASAKN( : ) ! pointers to aitken aerosols
      INTEGER, ALLOCATABLE, SAVE :: L_SRFAKN( : ) ! pntrs to aitken aerosol surface area
      INTEGER ASTAT
      integer, allocatable, save :: hlspc_map( : )

      INTEGER       I
      INTEGER, SAVE :: N_NUMAKN            ! # aitken aerosol number species
      INTEGER, SAVE :: N_MASAKN            ! # aitken aerosol mass species
      INTEGER, SAVE :: N_SRFAKN            ! # aitken aerosol sfc area species
      INTEGER       PNTR                ! relative pointer variable
      INTEGER       SPC                 ! liquid species loop counter
      INTEGER       VAR                 ! variable loop counter

      REAL          ALFA                ! scavenging coefficient (1/s)
      REAL          KH                  ! Henry's law constant (mol/l/atm)
      REAL          NUMAKN              ! Aitken mode aerosol # (#/m3)
      REAL          MASAKN              ! Total Aitken mode mass (ug/m3)
      REAL          SRFAKN
      REAL          ONE_OVER_TWASH      ! 1 / TWASH
      REAL          RTCH                ! chemical gas const times temp
      REAL          TWASH               ! washout time for clouds (sec) with low liq wat content
      REAL          TWF                 ! washout scaling factor (mol/l/atm)

      REAL, PARAMETER :: HPLUS = 1.0E-4 ! typical value hydrogen ion concentration [mol/l]

      REAL :: EALFAC
      REAL :: EALFAZ
      REAL :: EALFA
      logical :: do_rainout

      logical, allocatable, save :: cgrid_scav_isnum( : )
      logical, allocatable, save :: cgrid_scav_issrf( : )
      logical, allocatable, save :: cgrid_scav_ish2o( : )
      logical, allocatable, save :: cgrid_scav_istrc( : )
      logical, allocatable, save :: cgrid_scav_isaik( : )

C...........External Functions:

      INTEGER, EXTERNAL :: INDEXN

C-----------------------------------------------------------------------

C...INITIALIZATION SCAVWDEP module:
C...  event-statistics variables.  Open output files.

      !$omp critical
      IF ( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

        N_CGRID_SCAV = N_GC_SCAV + N_AE_SCAV + N_NR_SCAV

C...first check to make sure that some species in CGRID were specified
C...  for scavenging, otherwise notify the user and return

        IF ( N_CGRID_SCAV .LE. 0 ) THEN
          XMSG = 'No species were specified for scavenging by cloud ' //
     &           'or rain water...SCAVENGING WILL NOT BE PERFORMED!'
          CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
!         RETURN
        END IF

        ALLOCATE ( CGRID_SCAV    ( N_CGRID_SCAV ),
     &             CGRID_SCAV_MAP( N_CGRID_SCAV ),
     &             CGRID_SCAV_FAC( N_CGRID_SCAV ), STAT = ASTAT )
        IF ( ASTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating CGRID_SCAV, CGRID_SCAV_MAP or CGRID_SCAV_FAC'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        MXSPCS = N_GC_SPCD + N_AE_SPC + N_NR_SPC + N_TR_SPC

        ALLOCATE ( L_NUMAKN( MXSPCS ),
     &             L_MASAKN( MXSPCS ),
     &             L_SRFAKN( MXSPCS ), STAT = ASTAT )
        IF ( ASTAT .NE. 0 ) THEN
          XMSG = 'Failure allocating L_NUMAKN, L_MASAKN or L_SRFAKN'
          CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
        END IF

        allocate( hlspc_map( N_CGRID_SCAV ) )
        hlspc_map(:) = -1

C...prepare indices for scavenged species

C...  load the CGRID to scavenged species pointers for the gases

        SPC = 0
        DO VAR=1, N_GC_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = GC_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = GC_SCAV_MAP( VAR ) + GC_STRT - 1
          CGRID_SCAV_FAC( SPC ) = GC_SCAV_FAC( VAR )
          call hlconst_spcs_init( cgrid_scav( SPC ), hlspc_map( SPC ) )
        END DO

C...  load the CGRID to scavenged species pointers for the aerosols

        DO VAR=1, N_AE_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = AE_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = AE_SCAV_MAP( VAR ) + AE_STRT - 1
          CGRID_SCAV_FAC( SPC ) = AE_SCAV_FAC( VAR )
        END DO

C...  load the CGRID to scavenged species pointers for the non-reactives

        DO VAR=1, N_NR_SCAV
          SPC = SPC + 1
          CGRID_SCAV( SPC ) = NR_SCAV( VAR )
          CGRID_SCAV_MAP( SPC ) = NR_SCAV_MAP( VAR ) + NR_STRT - 1
          CGRID_SCAV_FAC( SPC ) = NR_SCAV_FAC( VAR )
          call hlconst_spcs_init( cgrid_scav( SPC ), hlspc_map( SPC ) )
        END DO

        do spc = 1, n_cgrid_scav
           if ( abs( CGRID_SCAV_FAC(SPC) - 1.0 ) > 1.e-7 ) then
              XMSG = 'Code assumes cgrid_scav_fac = 1'
              CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
           END IF
        enddo

C...create the pointers from CGRID to the species needed by AQCHEM

        N_NUMAKN = INDEXN( 'NUM_AITKEN      ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_NUMAKN )
        N_MASAKN = INDEXN( 'AITKEN', N_CGRID_SCAV, CGRID_SCAV, L_MASAKN )
        N_SRFAKN = INDEXN( 'SRF_AITKEN      ', N_CGRID_SCAV, CGRID_SCAV,
     &                     L_SRFAKN )

!...  store type of species

        allocate( cgrid_scav_isnum( N_CGRID_SCAV ) )
        allocate( cgrid_scav_issrf( N_CGRID_SCAV ) )
        allocate( cgrid_scav_ish2o( N_CGRID_SCAV ) )
        allocate( cgrid_scav_istrc( N_CGRID_SCAV ) )
        allocate( cgrid_scav_isaik( N_CGRID_SCAV ) )

        do var = 1,  N_CGRID_SCAV
           cgrid_scav_isnum( var ) = ( INDEX( CGRID_SCAV( var ), 'NUM'    ) .GT. 0 )
           cgrid_scav_issrf( var ) = ( INDEX( CGRID_SCAV( var ), 'SRF'    ) .GT. 0 )
           cgrid_scav_ish2o( var ) = ( INDEX( CGRID_SCAV( var ), 'H2O'    ) .GT. 0 )
           cgrid_scav_istrc( var ) = ( INDEX( CGRID_SCAV( var ), 'TRACER' ) .GT. 0 )
           cgrid_scav_isaik( var ) = ( INDEX( CGRID_SCAV( var ), 'AITKEN' ) .GT. 0 )
        enddo

      END IF  ! FIRSTIME
      !$omp end critical

C...for subsequent calls, check to make sure some species were
C...  specified, otherwise there is no need to perform scavenging

      IF ( N_CGRID_SCAV .LE. 0 ) THEN
        RETURN
      END IF

      RTCH = ( MOLVOL / STDTEMP ) * TBARC
      ONE_OVER_TWASH = PRATE1 / (WTBAR * CTHK1 * 3600.0)
      ONE_OVER_TWASH = MIN( ONE_OVER_TWASH, 1.0 / TAUCLD )
      TWF = H2ODENS / ( WTBAR * RTCH )
      REMOVAC = 0.0

      do_rainout = (ONE_OVER_TWASH * taucld > 1.e-7)

C...compute total Aitken mode number (#/m3)

      NUMAKN = 0.0
      DO I = 1, N_NUMAKN
        PNTR = CGRID_SCAV_MAP( L_NUMAKN( I ) )
        NUMAKN = NUMAKN + POLC( PNTR )
      END DO
      NUMAKN = NUMAKN * AIRM / CTHK1

C...compute total Aitken mode mass (ug/m3)

      MASAKN = 0.0
      DO I = 1, N_MASAKN
        PNTR = CGRID_SCAV_MAP( L_MASAKN( I ) )
        IF ( ( .not. cgrid_scav_isnum( L_MASAKN( I ) ) ) .AND.
     &       ( .not. cgrid_scav_issrf( L_MASAKN( I ) ) ) .AND.
     &       ( .not. cgrid_scav_ish2o( L_MASAKN( I ) ) ) .AND.
     &       ( .not. cgrid_scav_istrc( L_MASAKN( I ) ) ) ) THEN
          MASAKN = MASAKN + POLC( PNTR ) * AE_MOLWT( PNTR - AE_STRT + 1 )
        END IF
      END DO
      masakn = masakn * airm / cthk1 * 1.e6

C...compute total Aitken mode surface area (m2/m3)

      SRFAKN = 0.0
      DO I = 1, N_SRFAKN
        PNTR = CGRID_SCAV_MAP( L_SRFAKN( I ) )
        SRFAKN = SRFAKN + POLC( PNTR )
      END DO
      srfakn = srfakn * airm / cthk1

C...compute the scavenging coefficients for aitken mode aerosol mass
C...  and number
C...  NOTE:  for now, scavenging coefficients are computed for only
C...         the liquid water content, not on the total water content
C...         therefore, no ice phase scavenging is considered at this
C...         time, but it should be added in the future!

      CALL GETALPHA ( NUMAKN, MASAKN, SRFAKN, WCBAR, TBARC, PBARC,
     &                RHOAIR, ALFAC )

C...gas scavenging and wet deposition

      IF (DO_RAINOUT) THEN
         DO VAR = 1, N_GC_SCAV
            SPC = VAR
            PNTR = CGRID_SCAV_MAP( SPC )
            KH = HLCONST( TBARC, .TRUE., HPLUS, hlspc_map( spc ) )
            ALFA = ONE_OVER_TWASH * KH / ( KH + TWF )
            CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
            REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
         END DO
      ENDIF

C...aerosol scavenging and wet deposition

      EALFAC = EXP( -ALFAC * TAUCLD )
      if (do_rainout) then
         EALFAZ = EXP( -ONE_OVER_TWASH * TAUCLD )
      else
         EALFAZ = 1.0
      endif

      DO VAR = 1, N_AE_SCAV
        SPC = N_GC_SCAV + VAR

        PNTR = CGRID_SCAV_MAP( SPC )

        IF ( cgrid_scav_isaik( SPC ) ) THEN
           EALFA = EALFAC
        ELSE
           EALFA = EALFAZ
        END IF

        CEND ( PNTR ) = POLC( PNTR ) * EALFA
        REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM

      END DO

C...non-reactive scavenging and wet deposition

      IF (DO_RAINOUT) THEN
         DO VAR = 1, N_NR_SCAV
            SPC = N_GC_SCAV + N_AE_SCAV + VAR
            PNTR = CGRID_SCAV_MAP( SPC )
            KH = HLCONST( TBARC, .TRUE., HPLUS, hlspc_map( spc ) )
            ALFA = ONE_OVER_TWASH * KH / ( KH + TWF )
            CEND ( PNTR ) = POLC( PNTR ) * EXP( -ALFA * TAUCLD )
            REMOV( PNTR ) = ( POLC( PNTR ) - CEND( PNTR ) ) * AIRM
         END DO
      ENDIF

      END
