
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

      SUBROUTINE GAS_INCLOUD_RAINOUT( EWASH, WCBAR, TBARC, HPLUS, AIRM, CONC, REMOV )

!-----------------------------------------------------------------------
!  Description: Compute simplistic incloud scavenging and wet removal

!  Revision History:
!      No   Date   Who  What
!      -- -------- ---  -----------------------------------------
!       0 01/15/98 sjr  created program
!       1 07/27/99 sjr  added scavenging coefficient scaling factors
!       3 Dec 00   Jeff move CGRID_MAP into f90 module
!       4 Dec 02   sjr  revised calls to HLCONST
!       5 Jun 10   J.Young: convert for Namelist redesign
!       6 Mar 11   S.Roselle: replaced I/O API include files with UTILIO_DEFN;
!                       removed deprecated TRIMLEN
!       6 Mar 11   B.Hutzell: added TRACER string for AE_SPCs to omit while
!                       computing Aitken mode mass

!  Called by:  RADMCLD and RESCLD

!  Calls the following subroutines:  GETALPHA

!  Calls the following functions:  HLCONST

!  Arguments    Type      I/O       Description
!  ---------   -------  ------  --------------------------------
!    WCBAR      real     input  avg liquid water content (kg/m3)
!    TBARC      real     input  avg cloud temperature (K)
!    PBARC      real     input  avg cloud pressure (Pa)
!    CTHK1      real     input  cloud thickness (m)
!    AIRM       real     input  total air mass (moles/m2) in cloudy air
!    PRATE1     real     input  precipitation rate (mm/hr)
!    CONC       real     input  ave vert conc incloud (moles sp/ mole air)
!    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
!    REMOVAC    real    output  variable storing H+ deposition
!    ALFA0      real    output  scav coef for aitken aerosol number
!    ALFA2      real    output  scav coef for aitken mode sfc area
!    ALFA3      real    output  scav coef for aitken aerosol mass
!-----------------------------------------------------------------------

      USE CGRID_SPCS,  ONLY: NSPCSD
      USE CONST_DATA,  ONLY: MOLVOL, STDTEMP
      USE HLCONST_MOD, ONLY: HLCONST
      USE AQ_DATA,     ONLY: HLSPC_MAP, HLEFF_MAP, N_GAS_SCAV, GAS_SCAV_MAP

      IMPLICIT NONE

!...........Parameters:

      REAL, PARAMETER :: H2ODENS   = 1000.0   ! density of water at 20 C
                                              ! and 1 ATM (kg/m3)
      REAL, PARAMETER :: H2ODENSI  = 1.0 / H2ODENS

!...........Arguments:

      REAL,    INTENT( IN    ) :: EWASH   ! fraction liquid remaining after washout
      REAL,    INTENT( IN    ) :: WCBAR   ! liq water content of cloud (kg/m3)
      REAL,    INTENT( IN    ) :: TBARC   ! mean cloud temp (K)
      REAL,    INTENT( IN    ) :: HPLUS   ! hydrogen ion concentration [mol/l]
      REAL,    INTENT( IN    ) :: AIRM    ! total air mass (moles/m2) in cloudy air
      REAL,    INTENT( INOUT ) :: CONC ( NSPCSD+1 ) ! avg vert conc incloud (moles/mole)
      REAL,    INTENT(   OUT ) :: REMOV( NSPCSD+1 ) ! moles/m2 or mm*mol/lit scavenged

!...........Local Variables:

      INTEGER :: PNTR              ! relative pointer variable
      INTEGER :: SPC               ! liquid species loop counter

      REAL :: RTCH                 ! chemical gas const times temp
      REAL :: TWFI                 ! washout scaling factor (mol/l/atm)
      REAL :: KTWF
      REAL :: EWASHM               ! fraction liquid removed by washout
      REAL :: CONC0                ! starting concentration
      REAL :: LIQI(N_GAS_SCAV)     ! fraction of gas left after rainout
      REAL :: KH  (N_GAS_SCAV)     ! Henry's law constant (mol/l/atm)

      REAL, PARAMETER :: ONEM  = NEAREST(1.0, -1.0)
      REAL, PARAMETER :: RGAS  = MOLVOL / STDTEMP

!-----------------------------------------------------------------------

!...  for subsequent calls, check to make sure some species were
!...  specified, otherwise there is no need to perform scavenging

      IF ( N_GAS_SCAV .LE. 0) RETURN

!...  gas deposition only if there is liquid water and precipitation

      IF (WCBAR < 1.E-8 .OR. EWASH >= ONEM) RETURN

      EWASHM = 1.0 - EWASH

!...  get Henry's law constants for gases

      DO SPC = 1, N_GAS_SCAV
         KH(SPC) = HLCONST( TBARC, HLEFF_MAP( SPC ), HPLUS, HLSPC_MAP( SPC ) )
      ENDDO

!...compute loss term for gases

      RTCH = RGAS * TBARC
      TWFI = WCBAR * RTCH * H2ODENSI

      DO SPC = 1, N_GAS_SCAV
         KTWF = KH(SPC) * TWFI
         LIQI(SPC) = 1.0 - EWASHM * KTWF / ( 1.0 + KTWF )
      ENDDO

!...gas/non reactive wet deposition

      DO SPC = 1, N_GAS_SCAV
         PNTR = GAS_SCAV_MAP( SPC )

         CONC0       = CONC(PNTR)
         CONC (PNTR) = CONC(PNTR) * LIQI(SPC)
         REMOV(PNTR) = (CONC0 - CONC(PNTR)) * AIRM
      END DO

      END SUBROUTINE GAS_INCLOUD_RAINOUT



      SUBROUTINE GAS_INCLOUD_RAINOUT_SKIPAQ( EWASH, WCBAR, TBARC, HPLUS, AIRM, CONC, REMOV )

!-----------------------------------------------------------------------
!  Description: Compute simplistic incloud scavenging and wet removal

!  Revision History:
!      No   Date   Who  What
!      -- -------- ---  -----------------------------------------
!       0 01/15/98 sjr  created program
!       1 07/27/99 sjr  added scavenging coefficient scaling factors
!       3 Dec 00   Jeff move CGRID_MAP into f90 module
!       4 Dec 02   sjr  revised calls to HLCONST
!       5 Jun 10   J.Young: convert for Namelist redesign
!       6 Mar 11   S.Roselle: replaced I/O API include files with UTILIO_DEFN;
!                       removed deprecated TRIMLEN
!       6 Mar 11   B.Hutzell: added TRACER string for AE_SPCs to omit while
!                       computing Aitken mode mass

!  Called by:  RADMCLD and RESCLD

!  Calls the following subroutines:  GETALPHA

!  Calls the following functions:  HLCONST

!  Arguments    Type      I/O       Description
!  ---------   -------  ------  --------------------------------
!    WCBAR      real     input  avg liquid water content (kg/m3)
!    TBARC      real     input  avg cloud temperature (K)
!    PBARC      real     input  avg cloud pressure (Pa)
!    CTHK1      real     input  cloud thickness (m)
!    AIRM       real     input  total air mass (moles/m2) in cloudy air
!    PRATE1     real     input  precipitation rate (mm/hr)
!    CONC       real     input  ave vert conc incloud (moles sp/ mole air)
!    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
!    REMOVAC    real    output  variable storing H+ deposition
!    ALFA0      real    output  scav coef for aitken aerosol number
!    ALFA2      real    output  scav coef for aitken mode sfc area
!    ALFA3      real    output  scav coef for aitken aerosol mass
!-----------------------------------------------------------------------

      USE CGRID_SPCS,  ONLY: NSPCSD
      USE CONST_DATA,  ONLY: MOLVOL, STDTEMP
      USE HLCONST_MOD, ONLY: HLCONST
      USE AQ_DATA,     ONLY: HLSPC_NOAQ_MAP, HLEFF_NOAQ_MAP,
     &                       N_GAS_NOAQ_SCAV, GAS_NOAQ_SCAV_MAP

      IMPLICIT NONE

!...........Parameters:

      REAL, PARAMETER :: H2ODENS   = 1000.0   ! density of water at 20 C
                                              ! and 1 ATM (kg/m3)
      REAL, PARAMETER :: H2ODENSI  = 1.0 / H2ODENS

!...........Arguments:

      REAL,    INTENT( IN    ) :: EWASH   ! fraction liquid remaining after washout
      REAL,    INTENT( IN    ) :: WCBAR   ! liq water content of cloud (kg/m3)
      REAL,    INTENT( IN    ) :: TBARC   ! mean cloud temp (K)
      REAL,    INTENT( IN    ) :: HPLUS   ! hydrogen ion concentration [mol/l]
      REAL,    INTENT( IN    ) :: AIRM    ! total air mass (moles/m2) in cloudy air
      REAL,    INTENT( INOUT ) :: CONC ( NSPCSD+1 ) ! avg vert conc incloud (moles/mole)
      REAL,    INTENT(   OUT ) :: REMOV( NSPCSD+1 ) ! moles/m2 or mm*mol/lit scavenged

!...........Local Variables:

      INTEGER :: PNTR               ! relative pointer variable
      INTEGER :: SPC                ! liquid species loop counter

      REAL :: RTCH                  ! chemical gas const times temp
      REAL :: TWFI                  ! washout scaling factor (mol/l/atm)
      REAL :: KTWF
      REAL :: EWASHM                ! fraction liquid removed by washout
      REAL :: CONC0                 ! starting concentration
      REAL :: LIQI(N_GAS_NOAQ_SCAV) ! fraction of gas left after rainout
      REAL :: KH  (N_GAS_NOAQ_SCAV) ! Henry's law constant (mol/l/atm)

      REAL, PARAMETER :: ONEM  = NEAREST(1.0, -1.0)
      REAL, PARAMETER :: RGAS  = MOLVOL / STDTEMP

!-----------------------------------------------------------------------

!...  for subsequent calls, check to make sure some species were
!...  specified, otherwise there is no need to perform scavenging

      IF ( N_GAS_NOAQ_SCAV .LE. 0) RETURN

!...  gas deposition only if there is liquid water and precipitation

      IF (WCBAR < 1.E-8 .OR. EWASH >= ONEM) RETURN

      EWASHM = 1.0 - EWASH

!...  get Henry's law constants for gases

      DO SPC = 1, N_GAS_NOAQ_SCAV
         KH(SPC) = HLCONST( TBARC, HLEFF_NOAQ_MAP( SPC ), HPLUS, HLSPC_NOAQ_MAP( SPC ) )
      ENDDO

!...compute loss term for gases

      RTCH = RGAS * TBARC
      TWFI = WCBAR * RTCH * H2ODENSI

      DO SPC = 1, N_GAS_NOAQ_SCAV
         KTWF = KH(SPC) * TWFI
         LIQI(SPC) = 1.0 - EWASHM * KTWF / ( 1.0 + KTWF )
      ENDDO

!...gas/non reactive wet deposition

      DO SPC = 1, N_GAS_NOAQ_SCAV
         PNTR = GAS_NOAQ_SCAV_MAP( SPC )

         CONC0       = CONC(PNTR)
         CONC (PNTR) = CONC(PNTR) * LIQI(SPC)
         REMOV(PNTR) = (CONC0 - CONC(PNTR)) * AIRM
      END DO

      END SUBROUTINE GAS_INCLOUD_RAINOUT_SKIPAQ



      SUBROUTINE AERO_INCLOUD_RAINOUT( EWASH, AIRM, FR_ACT, CONC, REMOV )

!-----------------------------------------------------------------------
!  Description: Compute incloud aerosol removal by rainout

!  Revision History:
!      No   Date   Who  What
!      -- -------- ---  -----------------------------------------
!       0 01/15/98 sjr  created program
!       1 07/27/99 sjr  added scavenging coefficient scaling factors
!       3 Dec 00   Jeff move CGRID_MAP into f90 module
!       4 Dec 02   sjr  revised calls to HLCONST
!       5 Jun 10   J.Young: convert for Namelist redesign
!       6 Mar 11   S.Roselle: replaced I/O API include files with UTILIO_DEFN;
!                       removed deprecated TRIMLEN
!       6 Mar 11   B.Hutzell: added TRACER string for AE_SPCs to omit while
!                       computing Aitken mode mass

!  Calls the following functions:  HLCONST

!  Arguments    Type      I/O       Description
!  ---------   -------  ------  --------------------------------
!    EWASH      real     input  fraction of cloud remaining after rainout
!    AIRM       real     input  total air mass (moles/m2) in cloudy air
!    CONC       real    in/out  ave vert conc incloud (moles sp/ mole air)
!    REMOV      real    output  moles/m2 or mm*mol/lit scavenged
!-----------------------------------------------------------------------

      USE CGRID_SPCS,  ONLY: N_AE_SCAV, N_GC_SCAV, NSPCSD
      USE AQ_DATA,     ONLY: CGRID_SCAV_ISAIK, CGRID_SCAV_MAP

      IMPLICIT NONE

!...........Arguments:

      REAL,    INTENT( IN )    :: EWASH   ! fraction liquid remaining after washout
      REAL,    INTENT( IN )    :: AIRM    ! total air mass (moles/m2) in cloudy air
      REAL,    INTENT( IN )    :: FR_ACT  ! fraction of aerosol activated into clouds
      REAL,    INTENT( INOUT ) :: CONC ( NSPCSD+1 )  ! avg vert conc incloud (moles/mole)
      REAL,    INTENT(   OUT ) :: REMOV( NSPCSD+1 )  ! moles/m2 or mm*mol/lit scavenged

!...........Local Variables:

      INTEGER         :: PNTR   ! relative pointer variable
      INTEGER         :: SPC    ! liquid species loop counter
      INTEGER         :: VAR    ! variable loop counter
      REAL            :: CONC0  ! starting concentration
      REAL            :: REMAIN ! fraction aerosol remaining

!...........Parameters:

      REAL, PARAMETER :: ONEM  = NEAREST(1.0, -1.0)
      REAL, PARAMETER :: SMALL = EPSILON(1.0)

!-----------------------------------------------------------------------

!...  for subsequent calls, check to make sure some species were
!...  specified, otherwise there is no need to perform scavenging

      IF ( N_AE_SCAV .LE. 0 ) RETURN
      IF ( EWASH >= ONEM .or. FR_ACT <= SMALL ) RETURN

!...aerosol rainout/snowout
!...Aitken mode is always interstitial and not rained out

      REMAIN = 1.0 - (1.0 - EWASH) * FR_ACT

      IF (REMAIN >= ONEM) RETURN

      DO VAR = 1, N_AE_SCAV
         SPC = N_GC_SCAV + VAR
         IF ( .not. CGRID_SCAV_ISAIK( SPC ) ) THEN
            PNTR = CGRID_SCAV_MAP( SPC )

            CONC0       = CONC(PNTR)
            CONC (PNTR) = CONC(PNTR) * REMAIN
            REMOV(PNTR) = (CONC0 - CONC(PNTR)) * AIRM
         ENDIF
      END DO

      END SUBROUTINE AERO_INCLOUD_RAINOUT
