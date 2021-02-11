
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
      SUBROUTINE AQ_MAP ( k, iw, WCBAR, TBARC, PBARC, TAUCLD, CONC, AIRM,
     &                    ACTIV, REMOV, HPLUS, REMOVAC, EWASH, DARK )

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
C    WTBAR      real     input  avg total water content (kg/m3)
C    WCBAR      real     input  avg liquid water content (kg/m3)
C    TBARC      real     input  avg cloud temperature (K)
C    PBARC      real     input  avg cloud pressure (Pa)
C    TAUCLD     real     input  cloud lifetime (s)
C    CONC       real    output  incloud conc (moles/mole)
C    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
C    REMOVAC    real    output  variable storing H+ deposition
C    DARK       logical  input DARK = TRUE is night,  DARK = FALSE is day
C-----------------------------------------------------------------------

      USE CGRID_SPCS     ! CGRID mechanism species
      USE AQ_DATA        ! aqueous chemistry shared parameters
      USE AERO_DATA      ! aerosol shared parameters
      USE UTILIO_DEFN

      ! temp
      use cgrid_defn

      IMPLICIT NONE

!     CHARACTER( 120 ) :: XMSG = ' '    ! Exit status message

C Parameters:

!     INTEGER, SAVE :: MXSPCS   ! Number of species in CGRID

      REAL, PARAMETER :: ONETHIRD  = 1.0 / 3.0
      REAL, PARAMETER :: TWOTHIRDS = 2.0 / 3.0

C Arguments:

      integer, intent( in )    :: k, iw
      REAL,    INTENT( IN )    :: WCBAR     ! liq water content of cloud (kg/m3)
      REAL,    INTENT( IN )    :: TBARC     ! mean cloud temp (K)
      REAL,    INTENT( IN )    :: PBARC     ! mean cloud pressure (Pa)
      REAL,    INTENT( IN )    :: TAUCLD    ! cloud lifetime
      REAL,    INTENT( IN )    :: AIRM
      REAL,    INTENT( IN )    :: ACTIV     ! fraction of accum and coarse activated
      REAL,    INTENT( INOUT ) :: CONC ( NSPCSD+1 )  ! incloud conc (moles/mole)
      REAL,    INTENT( INOUT ) :: REMOV( NSPCSD+1 )  ! moles/m2 or mm*mol/lit scavenged
      REAL,    INTENT( INOUT ) :: REMOVAC   ! variable storing H+ deposition
      REAL,    INTENT( INOUT ) :: HPLUS     ! H+ concentration
      REAL,    INTENT( IN )    :: EWASH     ! scav coef for rainout
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

!     REAL       :: EBETASO4T
!     REAL       :: EALFA2T             ! EXP( -ALFAC * TAUCLD )
      REAL       :: M3OLD( 2:NMODES )        ! modal mass at time 0
      REAL       :: M3NEW( 2:NMODES )        ! modal mass at time t

      REAL       :: CONC2
      REAL       :: RATIO

      REAL( 8 )  :: GAS    ( NGAS )               ! gas phase conc (mol/mol)
      REAL       :: GASWDEP( NGAS )               ! gas phase wet dep array (mm mol/liter)
!     REAL( 8 )  :: WSRGGAS( NGAS,MXSRG )         ! weights for surrogate

      REAL( 8 )  :: AEROSOL(       NAER,2:NMODES )  ! aerosol conc (mol/mol)
      REAL( 8 )  :: WSRGAER( NXSRG,NAER,2:NMODES )  ! weights for surrogate

C-----------------------------------------------------------------------

C...  for subsequent calls, check to make sure some surrogates were
C...  specified, otherwise there is no need to perform aqueous chemistry

      IF (N_CGRID2AQ <= 0) RETURN

C...load gas-phase concentrations

      gaswdep = 0.0

      DO IGAS = 1, NGAS
         GAS( IGAS ) = CONC( LSRGGAS( IGAS ) )
      ENDDO

C...load aerosol concentrations

      aerosol = 0.0d0

      DO IMODE = 2, NMODES
         DO IAER = 1, NAER

            DO ISRG = 1, NSRGAER( IAER, IMODE )
               PNTR = LSRGAER( ISRG, IAER, IMODE )
               AEROSOL( IAER, IMODE ) = AEROSOL( IAER, IMODE ) + max(CONC( PNTR ), 1.e-36)
            END DO

            IF (AQ_UPDATE(IAER, IMODE) .AND. NSRGAER(IAER,IMODE) > 1) THEN
               DO ISRG = 1, NSRGAER( IAER, IMODE )
                  PNTR = LSRGAER( ISRG, IAER, IMODE )
                  WSRGAER( ISRG, IAER, IMODE ) = CONC( PNTR ) / AEROSOL( IAER, IMODE )
               END DO
            ENDIF

         END DO
      END DO

C...save original 3rd moments

      m3old = 0.0

      do iaer = 1, naccum
         isrg = iaccum_map(iaer)
         m3old(2) = m3old(2) + CONC(isrg) * accum_m3c(iaer)
      enddo

      do iaer = 1, ncoars
         isrg = icoars_map(iaer)
         m3old(3) = m3old(3) + CONC(isrg) * coars_m3c(iaer)
      enddo

      M3NEW = M3OLD

C...limit aerosol loading based on how much aerosol is activated
C   (for mixed-phase clouds)

      aerosol = aerosol * activ

C...perform aqueous-phase chemistry calculations

      CALL AQCHEM ( k, iw, TBARC, PBARC, TAUCLD, WCBAR, GAS, AEROSOL,
     &              GASWDEP, AIRM, HPLUS, REMOVAC, EWASH, DARK )

C...Now, re-apportion mass back into cend/remov (cgrid-type) array

      DO IGAS = 1, NGAS
         IF (LSRGGAS(IGAS) > 0) THEN
            PNTR = LSRGGAS( IGAS )
            CONC ( PNTR ) = MAX(GAS( IGAS ), 1.D-36)
            REMOV( PNTR ) = GASWDEP( IGAS )
         ENDIF
      END DO

      DO IMODE = 2, NMODES
         DO IAER = 1, NAER

            IF (AQ_UPDATE(IAER, IMODE)) THEN
               IF (NSRGAER(IAER,IMODE) > 1) THEN

                  DO ISRG = 1, NSRGAER( IAER, IMODE )
                     PNTR  = LSRGAER( ISRG, IAER, IMODE )

                     CONC2 = MAX(AEROSOL( IAER, IMODE ) * WSRGAER( ISRG, IAER, IMODE ), 1.D-36)
     &                     + (1. - activ) * CONC (PNTR)

                     M3NEW(IMODE) = M3NEW(IMODE) + (CONC2 - CONC(PNTR)) * ACCUM_M3C(IAER)
                     CONC (PNTR)  = CONC2
                  END DO

               ELSE

                  PNTR  = LSRGAER( 1, IAER, IMODE )
                  CONC2 = MAX(AEROSOL( IAER, IMODE ), 1.D-36) + (1. - ACTIV) * CONC (PNTR)

                  M3NEW(IMODE) = M3NEW(IMODE) + (CONC2 - CONC(PNTR)) * ACCUM_M3C(IAER)
                  CONC (PNTR)  = CONC2

               ENDIF
            ENDIF

         END DO
      END DO

C...Update surface area

      DO IMODE = 2, 3
         PNTR  = SRF_STR + IMODE - 1
         RATIO = ( (M3NEW(IMODE) / M3OLD(IMODE)) ** ONETHIRD ) ** 2

         CONC( PNTR ) = CONC( PNTR ) * RATIO
      ENDDO

      END SUBROUTINE AQ_MAP
