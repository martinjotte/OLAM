
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
C $Header: /project/yoj/arc/CCTM/src/cloud/cloud_acm_ae6/aq_map.F,v 1.7 2011/10/21 16:10:23 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE AQ_MAP ( JDATE, JTIME, WTBAR, WCBAR, TBARC, PBARC,
     &                    CTHK1, AIRM, PRATE1, TAUCLD, POLC, CEND,
     &                    REMOV, REMOVAC, ALFAC, DARK )

C-----------------------------------------------------------------------
C  DESCRIPTION:  This is an interface processor between the cloud dynamics
C     module(s) and the aqueous chemistry module.  It uses indices to
C     map species from their order in the array POLC, CEND, REMOV to
C     the species that are expected for AQCHEM, ie GAS and AEROSOL.
 
C  Revision History:
C      No   Date   Who  What
C      -- -------- ---  -----------------------------------------
C       0 01/15/98 sjr  created program
C       1 02/13/98 sjr  modified/generalized the mapping procedure
C         Dec 00   Jeff move CGRID_MAP into f90 module
C       3 06/07/05 sjr  added logic for coarse sulfate
C       4 04/11/08 jtk  added logic for coarse ammonium
C       5 10/10/10 yoj  update to use aero_reeng by Steve Howard, Prakash Bhave,
C                       Jeff Young, Sergey Napelenok, and Shawn Roselle
C       6 03/01/11 sjr  replaced I/O API include files with UTILIO_DEFN;
C       7 07/01/01 gs   calculate zenith angle to determine daytime and nightime 
C                       needed for sulfur oxidation via metal catalysis
C  Called by:  RADMCLD and RESCLD
 
C  Calls the following subroutines:  AQCHEM
 
C  ARGUMENTS    TYPE      I/O       DESCRIPTION
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
C    ALFAC      real     input  scav coef for aerosol
C    DARK       logical  input DARK = TRUE is night,  DARK = FALSE is day
C-----------------------------------------------------------------------

      USE CGRID_SPCS     ! CGRID mechanism species
      USE AQ_DATA        ! aqueous chemistry shared parameters
      USE AERO_DATA      ! aerosol shared parameters
      USE UTILIO_DEFN

      IMPLICIT NONE

!     CHARACTER( 120 ) :: XMSG = ' '    ! Exit status message

C Parameters:

!     INTEGER, SAVE :: MXSPCS   ! Number of species in CGRID

      REAL, PARAMETER :: ONETHIRD  = 1.0 / 3.0
      REAL, PARAMETER :: TWOTHIRDS = 2.0 / 3.0

C Arguments:

      INTEGER, INTENT( IN )    :: JDATE     ! current model date, coded YYYYDDD
      INTEGER, INTENT( IN )    :: JTIME     ! current model time, coded HHMMSS

      REAL,    INTENT( IN )    :: WTBAR     ! total wat cont (kg/m2) int. thru cld
      REAL,    INTENT( IN )    :: WCBAR     ! liq water content of cloud (kg/m3)
      REAL,    INTENT( IN )    :: TBARC     ! mean cloud temp (K)
      REAL,    INTENT( IN )    :: PBARC     ! mean cloud pressure (Pa)
      REAL,    INTENT( IN )    :: CTHK1     ! aq chem calc cloud thickness
      REAL,    INTENT( IN )    :: AIRM      ! total air mass (moles/m2) in cloudy air
      REAL,    INTENT( IN )    :: PRATE1    ! storm rainfall rate (mm/hr)
      REAL,    INTENT( IN )    :: TAUCLD    ! cloud lifetime
      REAL,    INTENT( IN )    :: POLC ( NSPCSD )  ! avg vert conc incloud (moles/mole)
      REAL,    INTENT( INOUT ) :: CEND ( NSPCSD )  ! ending incloud conc (moles/mole)
      REAL,    INTENT( INOUT ) :: REMOV( NSPCSD )  ! moles/m2 or mm*mol/lit scavenged
      REAL,    INTENT( INOUT ) :: REMOVAC   ! variable storing H+ deposition
      REAL,    INTENT( IN )    :: ALFAC     ! scav coef for aerosols
      LOGICAL, INTENT( IN )    :: DARK      ! DARK = TRUE is night,  DARK = FALSE is day

C Local Variables (scalars):

!     LOGICAL, SAVE :: FIRSTIME = .TRUE.      ! flag for first pass thru

!     CHARACTER(16), SAVE :: PNAME = 'AQ_MAP' ! program name

      INTEGER       IAER                ! aerosol loop counter
      INTEGER       IMODE               ! aerosol mode loop counter
      INTEGER       IGAS                ! gas loop counter
      INTEGER       ISRG                ! surrogate loop counter
      INTEGER       PNTR                ! relative pointer variable
!     INTEGER       SPC                 ! liquid species loop counter

      REAL       :: EBETASO4T
      REAL       :: EALFA2T             ! EXP( -ALFAC * TAUCLD )
      REAL       :: M3OLD( NMODES )     ! modal mass at time 0
      REAL       :: M3NEW( NMODES )     ! modal mass at time t
!     REAL( 8 )  :: SUM

      real       :: airm1
      REAL       :: HPWDEP                        ! hydrogen wet dep (mm mol/liter)

      REAL( 8 )  :: GAS    ( NGAS )               ! gas phase conc (mol/mol)
      REAL( 8 )  :: GASWDEP( NGAS )               ! gas phase wet dep array (mm mol/liter)
      REAL( 8 )  :: WSRGGAS( NGAS,MXSRG )         ! weights for surrogate

      REAL( 8 )  :: AEROSOL( MAX_NAER,NMODES )        ! aerosol conc (mol/mol)
      REAL( 8 )  :: AERWDEP( MAX_NAER,NMODES )        ! aerosol wet dep array (mm mol/liter)
      REAL( 8 )  :: WSRGAER( MAX_NAER,NMODES,MXSRG )  ! weights for surrogate
      real       :: ratio1, ratio2, ratio3

C-----------------------------------------------------------------------

C...  for subsequent calls, check to make sure some surrogates were
C...  specified, otherwise there is no need to perform aqueous chemistry

      IF (N_CGRID2AQ <= 0) RETURN

C...load gas-phase concentrations

      gaswdep = 0.0d0
      gas     = 0.0d0

      DO IGAS = 1, NGAS

         DO ISRG = 1, NSRGGAS( IGAS )
            PNTR = CGRID2AQ_MAP( LSRGGAS( IGAS, ISRG ) )
            GAS( IGAS ) = GAS( IGAS ) + POLC( PNTR )
         END DO

        IF      ( NSRGGAS( IGAS ) == 0 ) THEN
           GAS( IGAS ) = SRGGAS( IGAS )%BACKGND
        ELSE IF ( NSRGGAS( IGAS ) == 1 ) THEN
           WSRGGAS( IGAS,1 ) = 1.0D0
        ELSE
           DO ISRG = 1, NSRGGAS( IGAS )
              PNTR = CGRID2AQ_MAP( LSRGGAS( IGAS, ISRG ) )
              WSRGGAS( IGAS,ISRG ) = POLC( PNTR ) / GAS( IGAS )
           END DO
        END IF

      END DO

C...load aerosol concentrations

      airm1   = 1.0 / airm
      aerwdep = 0.0d0
      aerosol = 0.0d0

      DO IAER = 1, NAER
         DO IMODE = 1, NMODES

            IF ( srgaerm( iaer,imode ) ) THEN

               DO ISRG = 1, NSRGAER( IAER, IMODE )
                  PNTR = CGRID2AQ_MAP( LSRGAER( IAER, IMODE, ISRG ) )
                  AEROSOL( IAER, IMODE ) = AEROSOL( IAER, IMODE ) + POLC( PNTR )
               END DO

               IF      ( NSRGAER(IAER,IMODE) == 0 ) THEN
                  AEROSOL( IAER, IMODE ) = SRGAER( IAER )%BACKGND * CTHK1 * AIRM1
               ELSE IF ( NSRGAER(IAER,IMODE) == 1 ) THEN
                  WSRGAER( IAER, IMODE, 1 ) = 1.0D0
               ELSE
                  DO ISRG = 1, NSRGAER( IAER, IMODE )
                     PNTR = CGRID2AQ_MAP( LSRGAER( IAER, IMODE, ISRG ) )
                     WSRGAER( IAER, IMODE, ISRG ) = POLC( PNTR ) / AEROSOL( IAER, IMODE )
                  END DO
               END IF

            END IF
         END DO
      END DO

C *** extract grid cell concentrations of aero species from CGRID
C     into aerospc_conc in aero_data module

      CALL EXTRACT_AERO ( POLC, .FALSE., dry=.true. )

      m3old(2:3) = moment3_conc(2:3)

C...perform aqueous-phase chemistry calculations

      CALL AQCHEM ( JDATE, JTIME, TBARC, PBARC, TAUCLD, PRATE1,
     &              WCBAR, WTBAR, AIRM, ALFAC, GAS, AEROSOL, GASWDEP,
     &              AERWDEP, HPWDEP, EALFA2T, EBETASO4T, DARK, cthk1 )

C...store the amount of hydrogen deposition

      REMOVAC = HPWDEP

C...Now, re-apportion mass back into cend/remov (cgrid-type) array

      DO IGAS = 1, NGAS
        DO ISRG = 1, NSRGGAS( IGAS )
          PNTR = CGRID2AQ_MAP( LSRGGAS( IGAS, ISRG ) )
          CEND ( PNTR ) = GAS    ( IGAS ) * WSRGGAS( IGAS, ISRG )
          REMOV( PNTR ) = GASWDEP( IGAS ) * WSRGGAS( IGAS, ISRG )
        END DO
      END DO

      DO IAER = 1, NAER
         DO IMODE = 1, NMODES
            DO ISRG = 1, NSRGAER( IAER, IMODE )
               PNTR = CGRID2AQ_MAP( LSRGAER( IAER, IMODE, ISRG ) )
               CEND ( PNTR ) = AEROSOL( IAER, IMODE ) * WSRGAER( IAER, IMODE, ISRG )
               REMOV( PNTR ) = AERWDEP( IAER, IMODE ) * WSRGAER( IAER, IMODE, ISRG )
            END DO
         END DO
      END DO

C *** extract grid cell concentrations of aero species from CGRID
C     into aerospc_conc in aero_data module

      CALL EXTRACT_AERO ( CEND, .FALSE., dry=.true. )

      m3new(2:3) = moment3_conc(2:3)

C...Update surface area

      ratio1 = CEND( AERONUM_MAP(3) ) / POLC( AERONUM_MAP(3) )
      ratio2 = m3new(2) / m3old(2)
      ratio3 = m3new(3) / m3old(3)

      CEND( AEROSRF_MAP( 1 ) ) = POLC( AEROSRF_MAP( 1 ) ) * EALFA2T

      CEND( AEROSRF_MAP( 2 ) ) = POLC( AEROSRF_MAP( 2 ) )
     &                         * ( ebetaso4t * ratio2**2 ) ** ONETHIRD

      CEND( AEROSRF_MAP( 3 ) ) = POLC( AEROSRF_MAP( 3 ) )
     &                         * ( ratio1 * ratio3**2 ) ** ONETHIRD

      END SUBROUTINE AQ_MAP
