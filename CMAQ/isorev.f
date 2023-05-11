
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

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE ISRP1R (WI, RHI, TEMPI)
      use isrpia

      implicit none

      REAL(8), INTENT(INOUT) :: WI(NCOMP)
      REAL,    INTENT(IN)    :: RHI, TEMPI
C
C *** INITIALIZE COMMON BLOCK VARIABLES *********************************
C
      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      IF (RH.GE.DRNH42S4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
         SULRATW = GETASR(WAER(2), RHI)     ! AEROSOL SULFATE RATIO
      ELSE
         SULRATW = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
      ENDIF
      SULRAT  = WAER(3)/WAER(2)         ! SULFATE RATIO
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR
C
      IF (SULRAT.GE.SULRATW) THEN

         SCASE = ISS
         CALL CALCS2                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (SULRAT.GE.1.D0 .AND. SULRAT.LT.SULRATW) THEN
         W(2) = WAER(2)
         W(3) = WAER(3)

         SCASE = IBS
         CALL CALCB4                 ! Only liquid (metastable)
         CALL CALCNH3P               ! Compute NH3(g)
C
C *** SULFATE RICH (FREE ACID)
C
      ELSE
         W(2) = WAER(2)
         W(3) = WAER(3)

         SCASE = ICS
         CALL CALCC2                 ! Only liquid (metastable)
         CALL CALCNH3P

      ENDIF

      END SUBROUTINE ISRP1R

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP2R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE ISRP2R (WI, RHI, TEMPI)
      use isrpia

      implicit none

      REAL(8), INTENT(INOUT) :: WI(NCOMP)
      REAL,    INTENT(IN)    :: RHI, TEMPI
      LOGICAL                :: TRYLIQ
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      TRYLIQ = .TRUE.             ! Assume liquid phase, sulfate poor limit

10    CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN ! *** WET AEROSOL
         SULRATW = GETASR(WAER(2), RHI)     ! LIMITING SULFATE RATIO
      ELSE
         SULRATW = 2.0D0                    ! *** DRY AEROSOL
      ENDIF
      SULRAT = WAER(3)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR
C
      IF (SULRAT.GE.SULRATW) THEN

         SCASE = INS
         CALL CALCN3                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
C     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
C     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE
C     AEROSOL EQUILIBRIUM.
C
      ELSEIF (SULRAT.GE.1.D0 .AND. SULRAT.LT.SULRATW) THEN
         W(2) = WAER(2)
         W(3) = WAER(3)
         W(4) = WAER(4)

         SCASE = IBS
         CALL CALCB4                 ! Only liquid (metastable)
C
C *** Add the NO3 to the solution now and calculate partitioning.
C
         MOLAL(7) = WAER(4)            ! There is always water, so NO3(aer) is NO3-
         MOLAL(1) = MOLAL(1) + WAER(4) ! Add H+ to balance out
         CALL CALCNAP                  ! HNO3, NH3 dissolved
         CALL CALCNH3P
C
C *** SULFATE RICH (FREE ACID)
C
C     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
C     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE
C     AEROSOL EQUILIBRIUM.
C
      ELSE
         W(2) = WAER(2)
         W(3) = WAER(3)
         W(4) = WAER(4)

         SCASE = ICS
         CALL CALCC2                 ! Only liquid (metastable)
C
C *** Add the NO3 to the solution now and calculate partitioning.
C
         MOLAL(7) = WAER(4)            ! There is always water, so NO3(aer) is NO3-
         MOLAL(1) = MOLAL(1) + WAER(4) ! Add H+ to balance out

         CALL CALCNAP                   ! HNO3, NH3 dissolved
         CALL CALCNH3P
      ENDIF
C
C *** IF SULRATW < SULRAT < 2.0 and WATER = 0 => SULFATE RICH CASE.
C
      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.D0
     &                      .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF
C
      END SUBROUTINE ISRP2R


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP3R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE ISRP3R (WI, RHI, TEMPI)

      USE ISRPIA
      IMPLICIT NONE

      REAL(8), INTENT(INOUT) :: WI(NCOMP)
      REAL,    INTENT(IN)    :: RHI, TEMPI

      REAL(8)                :: FRSO4, SRI
      LOGICAL                :: TRYLIQ
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
c     WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
c     WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** INITIALIZE ALL VARIABLES ******************************************
C
      TRYLIQ = .TRUE.             ! Use liquid phase sulfate poor limit

10    CALL ISOINIT3 (WI, RHI, TEMPI) ! COMMON block variables
C
C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
C
c     REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
c     IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
c        WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
c     ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN  ! ** WET AEROSOL
         FRSO4   = WAER(2) - WAER(1)/2.0D0     ! SULFATE UNBOUND BY SODIUM
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR(FRSO4, RHI)          ! SULFATE RATIO FOR NH4+
         SULRATW = (WAER(1)+FRSO4*SRI)/WAER(2) ! LIMITING SULFATE RATIO
         SULRATW = MIN (SULRATW, 2.0D0)
      ELSE
         SULRATW = 2.0D0                     ! ** DRY AEROSOL
      ENDIF
      SULRAT = (WAER(1)+WAER(3))/WAER(2)
      SODRAT = WAER(1)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM POOR
C
      IF (SULRAT.GE.SULRATW .AND. SODRAT.LT.2.D0) THEN

         SCASE = IQS
         CALL CALCQ5                 ! Only liquid (metastable)
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.D0) THEN

         SCASE = IRS
         CALL CALCR6                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (SULRAT.GE.1.D0 .AND. SULRAT.LT.SULRATW) THEN

         W(1:NCOMP) = WAER(1:NCOMP)

         SCASE = IIS
         CALL CALCI6                 ! Only liquid (metastable)
         CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
         CALL CALCNH3P
C
C *** SULFATE RICH (FREE ACID)
C
      ELSE

         W(1:NCOMP) = WAER(1:NCOMP)

         SCASE = IJS
         CALL CALCJ3                 ! Only liquid (metastable)
         CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
         CALL CALCNH3P

      ENDIF
C
C *** IF AFTER CALCULATIONS, SULRATW < SULRAT < 2.0
C                            and WATER = 0          => SULFATE RICH CASE.
C
      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.D0
     &                      .AND. WATER.LE.TINY4) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      END SUBROUTINE ISRP3R



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE ISRP4R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTTASIUM-MAGNESIUM AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE ISRP4R (WI, RHI, TEMPI)
      USE ISRPIA

      IMPLICIT NONE

      REAL(8), INTENT(INOUT) :: WI(NCOMP)
      REAL, INTENT(IN)       :: RHI, TEMPI

      REAL(8)                :: FRSO4, SRI
      LOGICAL                :: TRYLIQ
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
c     WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
c     WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** INITIALIZE ALL VARIABLES ******************************************
C
      TRYLIQ  = .TRUE.             ! Use liquid phase sulfate poor limit
10    CALL INIT4 (WI, RHI, TEMPI) ! COMMON block variables
C
C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
C
c     REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
c     IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
c        WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
c     ENDIF
C
C *** CALCULATE SULFATE, CRUSTAL & SODIUM RATIOS ***********************
C
      IF (TRYLIQ) THEN                               ! ** WET AEROSOL
         ! SULFATE UNBOUND BY SODIUM,CALCIUM,POTTASIUM,MAGNESIUM
         FRSO4   = WAER(2) - WAER(1)/2.0D0
     &           - WAER(6) - WAER(7)/2.0D0 - WAER(8)
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR(FRSO4, RHI)                ! SULFATE RATIO FOR NH4+
         SULRATW = (WAER(1)+FRSO4*SRI+WAER(6)
     &              +WAER(7)+WAER(8))/WAER(2)       ! LIMITING SULFATE RATIO
         SULRATW = MIN (SULRATW, 2.0D0)
      ELSE
         SULRATW = 2.0D0                     ! ** DRY AEROSOL
      ENDIF
      SO4RAT  = (WAER(1)+WAER(3)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
      CRNARAT = (WAER(1)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
      CRRAT   = (WAER(6)+WAER(7)+WAER(8))/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM+CRUSTALS POOR
C
      IF (SO4RAT.GE.SULRATW .AND. CRNARAT.LT.2.D0) THEN

         SCASE = IVS
         CALL CALCV7                 ! Only liquid (metastable)
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
      ELSEIF (SO4RAT.GE.SULRATW .AND. CRNARAT.GE.2.D0) THEN

       IF (CRRAT.LE.2.D0) THEN

         SCASE = IUS
         CALL CALCU8                 ! Only liquid (metastable)
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
       ELSE

         SCASE = IWS
         CALL CALCW13                 ! Only liquid (metastable)

      ENDIF
C
C *** SULFATE RICH (NO ACID): 1<Rso4<2;
C
      ELSEIF (SO4RAT.GE.1.D0 .AND. SO4RAT.LT.SULRATW) THEN

         W(1:NCOMP) = WAER(1:NCOMP)

         SCASE = ILS
         CALL CALCL9                 ! Only liquid (metastable)
         CALL CALCNHP                ! MINOR SPECIES: HNO3, HCl
         CALL CALCNH3P               !                NH3
C
C *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
C
      ELSE

         W(1:NCOMP) = WAER(1:NCOMP)

         SCASE = IKS
         CALL CALCK4                 ! Only liquid (metastable)
         CALL CALCNHP                  ! MINOR SPECIES: HNO3, HCl
         CALL CALCNH3P                 !                NH3
C
      ENDIF
C
C *** IF AFTER CALCULATIONS, SULRATW < SO4RAT < 2.0
C                            and WATER = 0          => SULFATE RICH CASE.
C
      IF (SULRATW.LE.SO4RAT .AND. SO4RAT.LT.2.D0
     &                      .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF

      END SUBROUTINE ISRP4R


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCS2
C *** CASE S2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCS2
      use isrpia

      implicit none

      real   :: nh4i, so4i, hso4i, hi, ohi, nh3aq, del
      real   :: a2, akw, hi0, ohi0, molalr(4:4)
      integer :: i
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU   =.TRUE.     ! Outer loop activity calculation flag
      FRST     =.TRUE.
      CALAIN   =.TRUE.
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLALR(4)= MIN(WAER(2), 0.5d0*WAER(3))
      WATER    = MOLALR(4)*M0I(4)  ! ZSR correlation
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
         NH4I = WAER(3)
         SO4I = WAER(2)
         HSO4I= 0.0

         IF (I == 1) THEN
            CALL CALCPH (2.*SO4I - NH4I, HI0, OHI0) ! GET PH
            if (HI0 <= OHI0) AKW = SQRT(XKW*RH)*WATER
         ENDIF

         HI  = HI0
         OHI = OHI0

         NH3AQ = 0.0                              ! AMMONIA EQUILIBRIUM
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ (NH4I, OHI, DEL)
            NH4I  = MAX (NH4I-DEL, 0.0)
            OHI   = MAX (OHI -DEL, TINY4)
            NH3AQ = DEL
            HI    = (AKW/OHI) * AKW
         ENDIF

         CALL CALCHS4 (HI, SO4I, 0.0, DEL)         ! SULFATE EQUILIBRIUM
         SO4I  = SO4I - DEL
         HI    = HI   - DEL
         HSO4I = DEL
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL(1) = HI
         MOLAL(3) = NH4I
         MOLAL(5) = SO4I
         MOLAL(6) = HSO4I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            EXIT
         ENDIF

      ENDDO

      A2 = XK2*R*TEMP*GAMA(8)**2 / (XKW*RH*GAMA(9)**2)

      COH      = OHI
      GASAQ(1) = NH3AQ
      GNH3     = NH4I / (HI * A2)

      END SUBROUTINE CALCS2


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCN3
C *** CASE N3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCN3
      use isrpia

      implicit none

      real    :: NH4I, NO3I, NH3AQ, NO3AQ, a2, a3i, akw, so4i, hso4i
      real    :: hi0, ohi0, hi, ohi, del, gg, molalr(4:5), rtw
      real(8) :: aml5
      integer :: i
C
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU =.TRUE.              ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** AEROSOL WATER CONTENT
C
      MOLALR(4) = MIN(WAER(2),0.5d0*WAER(3))       ! (NH4)2SO4
      AML5      = MAX(WAER(3)-2.D0*MOLALR(4),ZERO) ! "free" NH4
      MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO)     ! NH4NO3=MIN("free",NO3)

      WATER     = MOLALR(4)*M0I(4) + MOLALR(5)*M0I(5)
      WATER     = MAX(WATER, TINY4)
      AKW       = SQRT(XKW*RH)*WATER
      RTW       = R*TEMP*WATER
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
C ION CONCENTRATIONS
C
         NH4I  = WAER(3)
         NO3I  = WAER(4)
         SO4I  = WAER(2)
         HSO4I = ZERO4

         IF (I == 1) THEN
            GG = 2.*SO4I + NO3I - NH4I
            CALL CALCPH (GG, HI0, OHI0)
         ENDIF

         HI  = HI0
         OHI = OHI0
C
C AMMONIA ASSOCIATION EQUILIBRIUM
C
         NH3AQ = ZERO4
         NO3AQ = ZERO4
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = (AKW/OHI) * AKW
         ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) ! HNO3
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
            CALL CALCHS4 (HI, SO4I, 0.0, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = (AKW/HI) * AKW
         ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
         MOLAL (1) = HI
         MOLAL (3) = NH4I
         MOLAL (5) = SO4I
         MOLAL (6) = HSO4I
         MOLAL (7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            EXIT
         ENDIF

      ENDDO

      A2  = XK2*R*TEMP*GAMA(8)**2 / (XKW*RH*GAMA(9)**2)
      A3I = GAMA(10) / (RTW*XK4*WATER)

      GHNO3     = HI*NO3I*A3I
      GNH3      = NH4I/(HI*A2)  !   NH3AQ/A21

      GASAQ(1)  = NH3AQ
      GASAQ(3)  = NO3AQ

      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      COH       = OHI

      END SUBROUTINE CALCN3



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCQ5
C *** CASE Q5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
C     2. LIQUID AND SOLID PHASES ARE POSSIBLE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCQ5
      use isrpia

      implicit none

      real    :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, so4i, cli, gg, ohi
      real    :: hi, sdd, hso4i, ggno3, ggcl, del, akw, a2, rtw
      real    :: hi0, ohi0, a3i, a4i, molalr(2:6)
      integer :: i
C
C *** SETUP PARAMETERS ************************************************
C
      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
      CALL CALCQ1A

      MOLALR(2) = CNA2SO4                                  ! NA2SO4
      MOLALR(4) = CNH42S4                                  ! (NH4)2SO4
      MOLALR(5) = CNH4NO3                                  ! NH4NO3
      MOLALR(6) = CNH4CL                                   ! NH4CL

      WATER  = SUM(MOLALR(4:6)*M0I(4:6)) + MOLALR(2)*M0I(2)
      RTW    = R*TEMP*WATER
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      AKW = sqrt(XKW*RH) * WATER

      DO I = 1, NSWEEP
C
C ION CONCENTRATIONS
C
      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
C
C SOLUTION ACIDIC OR BASIC?
C
      IF (I == 1) THEN
         GG  = 2.*SO4I + NO3I + CLI - NAI - NH4I
         SDD = SQRT(GG*GG + 4.*AKW*AKW)
         SDD = MAX(SDD, ABS(GG), 2.*AKW)

         IF (GG > 0.0) THEN           ! H+ in excess
            HI0 = 0.5 * (SDD + GG)
            OHI0= (AKW / HI0) * AKW
         ELSE                         ! OH- in excess
            OHI0= 0.5 * (SDD - GG)
            HI0 = (AKW / OHI0) * AKW
         ENDIF
      ENDIF

      HI  = HI0
      OHI = OHI0
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      NH3AQ  = ZERO4
      NO3AQ  = ZERO4
      CLAQ   = ZERO4

      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = (AKW/OHI) * AKW
         HSO4I = 0.0
      ELSE
         GGNO3 = MAX(GG-CLI, 0.0)
         GGCL  = MAX(GG-GGNO3, 0.0)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = 0.0
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, 0.0, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = (AKW/HI)*AKW
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         EXIT
      ENDIF

      ENDDO
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2      = XK2*R*TEMP*GAMA(10)**2/(XKW*GAMA(5)**2) ! NH3  <==> NH4+
      A3I     = GAMA(10)**2 / (RTW*XK4*WATER)    ! HNO3 <==> NO3-
      A4I     = GAMA(11)**2 / (RTW*XK3*WATER)    ! HCL  <==> CL-
C
      GNH3    = NH4I/(HI*A2)
      GHNO3   = HI*NO3I*A3I
      GHCL    = HI*CLI *A4I
C
      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ
C
      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = ZERO
C
      END SUBROUTINE CALCQ5



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCQ1A
C *** CASE Q1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCQ1A
      use isrpia

      implicit none

      REAL(8) :: FRSO4, FRNH3
C
C *** CALCULATE SOLIDS **************************************************
C
      CNA2SO4 = 0.5d0*WAER(1)
      FRSO4   = MAX (WAER(2)-CNA2SO4, ZERO)
C
      CNH42S4 = MAX (MIN(FRSO4,0.5d0*WAER(3)), TINY)
      FRNH3   = MAX (WAER(3)-2.D0*CNH42S4, ZERO)
C
      CNH4NO3 = MIN (FRNH3, WAER(4))
      FRNH3   = MAX (FRNH3-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRNH3, WAER(5))
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      END SUBROUTINE CALCQ1A



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCR6
C *** CASE R6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCR6
      use isrpia

      implicit none

      real    :: NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, rtw
      real    :: gg, sdd, hi, ohi, so4i, cli, hso4i, ggno3, ggcl, del
      real    :: a2, a3i, a4i, akw, hi0, ohi0, molalr(6)
      integer :: i
C
C *** SETUP PARAMETERS ************************************************
C
      CALL CALCR1A

      FRST   = .TRUE.
      CALAIN = .TRUE.
      CALAOU = .TRUE.
C
C *** CALCULATE WATER **************************************************
C
      MOLALR(1) = CNACL         ! NACL
      MOLALR(2) = CNA2SO4       ! NA2SO4
      MOLALR(3) = CNANO3        ! NANO3
      MOLALR(4) = ZERO          ! (NH4)2SO4
      MOLALR(5) = CNH4NO3       ! NH4NO3
      MOLALR(6) = CNH4CL        ! NH4CL

      WATER = SUM( MOLALR(1:6) * M0I(1:6) )
      RTW   = R*TEMP*WATER
C
C *** SETUP LIQUID CONCENTRATIONS **************************************
C
      AKW = SQRT(XKW*RH)*WATER                           ! H2O    <==> H+
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
C
C SOLUTION ACIDIC OR BASIC?
C
      if (I == 1) then
         GG  = 2.*SO4I + NO3I + CLI - NAI - NH4I
         SDD = SQRT(GG*GG + 4.*AKW*AKW)
         SDD = MAX(SDD, ABS(GG), 2.*AKW)

         IF (GG > 0.0) THEN           ! H+ in excess
            HI0 = 0.5 * (SDD + GG)
            OHI0= (AKW/HI0) * AKW
         ELSE                         ! OH- in excess
            OHI0= 0.5 * (SDD - GG)
            HI0 = (AKW/OHI0) * AKW
         ENDIF
      ENDIF

      HI  = HI0
      OHI = OHI0
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      HSO4I  = ZERO4
      NH3AQ  = ZERO4
      NO3AQ  = ZERO4
      CLAQ   = ZERO4

      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI = (AKW/OHI) * AKW
      ELSE
         GGNO3 = MAX(GG-CLI, 0.0)
         GGCL  = MAX(GG-GGNO3, 0.0)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = 0.0
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, 0.0, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = (AKW/HI) * AKW
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         EXIT
      ENDIF

      ENDDO
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2       = XK2*R*TEMP*GAMA(10)**2/(XKW*GAMA(5)**2)  ! NH3  <==> NH4+
      A3I      = GAMA(10)**2 / (RTW*XK4*WATER)      ! HNO3 <==> NO3-
      A4I      = GAMA(11)**2 / (RTW*XK3*WATER)      ! HCL  <==> CL-

      GNH3     = NH4I/(HI*A2)
      GHNO3    = HI*NO3I*A3I
      GHCL     = HI*CLI *A4I

      GASAQ(1) = NH3AQ
      GASAQ(2) = CLAQ
      GASAQ(3) = NO3AQ

      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = ZERO
C
      END SUBROUTINE CALCR6


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCR1A
C *** CASE R1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4, NANO3, NACL
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCR1A
      use isrpia

      implicit none

      REAL(8) :: FRNO3, FRNA, FRCL, FRNH3
C
C *** CALCULATE SOLIDS **************************************************
C
      CNA2SO4 = WAER(2)
      FRNA    = MAX (WAER(1)-2*CNA2SO4, ZERO)
C
      CNH42S4 = ZERO
C
      CNANO3  = MIN (FRNA, WAER(4))
      FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
      FRNA    = MAX (FRNA-CNANO3, ZERO)
C
      CNACL   = MIN (FRNA, WAER(5))
      FRCL    = MAX (WAER(5)-CNACL, ZERO)
      FRNA    = MAX (FRNA-CNACL, ZERO)
C
      CNH4NO3 = MIN (FRNO3, WAER(3))
      FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
      FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRCL, FRNH3)
      FRCL    = MAX (FRCL-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      END SUBROUTINE CALCR1A


C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCV7
C *** CASE V7
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCV7
      use isrpia

      implicit none

      real    :: NH4I, NAI, NO3I, CAI, KI, MGI
      real    :: gg, sdd, hi, ohi, so4i, cli, hso4i, akw, del, rtw
      real    :: a2, a3i, a4i, hi0, ohi0, molalr(21)
      integer :: i
C
C *** SETUP PARAMETERS ************************************************
C
      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
      CALL CALCV1A

      MOLALR(2) = CNA2SO4       ! NA2SO4
      MOLALR(4) = CNH42S4       ! (NH4)2SO4
      MOLALR(5) = CNH4NO3       ! NH4NO3
      MOLALR(6) = CNH4CL        ! NH4CL
      MOLALR(17)= CK2SO4        ! K2SO4
      MOLALR(21)= CMGSO4        ! MGSO4

      WATER  = SUM(MOLALR(4:6)*M0I(4:6)) + MOLALR(2)*M0I(2)
     &       + MOLALR(17)*M0I(17) + MOLALR(21)*M0I(21)

      RTW = R*TEMP*WATER
      AKW = SQRT(XKW*RH) * WATER
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
C ION CONCENTRATIONS
C
      NAI    = WAER(1)
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = 0.0
      KI     = WAER(7)
      MGI    = WAER(8)
C
C SOLUTION ACIDIC OR BASIC?
C
      IF (I == 1) THEN
         GG  = 2.*SO4I + NO3I + CLI - NAI - NH4I - KI - 2.*MGI
         SDD = SQRT(GG*GG + 4.*AKW*AKW)
         SDD = MAX(SDD, ABS(GG), 2.*AKW)

         IF (GG > 0.0) THEN               ! H+ in excess
            HI0 = 0.5 * (SDD + GG)
            OHI0= (AKW/HI0) * AKW
         ELSE                             ! OH- in excess
            OHI0= 0.5 * (SDD - GG)
            HI0 = (AKW/OHI0) * AKW
         ENDIF
      ENDIF

      HI  = HI0
      OHI = OHI0
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.GT.OHI) THEN
         ! CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
         CALL CALCHS4 (HI, SO4I, 0.0, DEL)
      ELSE
         DEL= 0.0
      ENDIF

      SO4I  = SO4I  - DEL
      HSO4I = DEL
      HI    = HI    - DEL
      HI    = MAX(HI, AKW)
      OHI   = (AKW/HI) * AKW
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         EXIT
      ENDIF

      ENDDO
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2      = XK2*R*TEMP*GAMA(10)**2/(XKW*GAMA(5)**2) ! NH3  <==> NH4+
      A3I     = GAMA(10)**2/(RTW*XK4*WATER)       ! HNO3 <==> NO3-
      A4I     = GAMA(11)**2/(RTW*XK3*WATER)      ! HCL  <==> CL-

      GNH3    = NH4I/(HI*A2)
      GHNO3   = HI*NO3I*A3I
      GHCL    = HI*CLI *A4I

      GASAQ(1)= ZERO
      GASAQ(2)= ZERO
      GASAQ(3)= ZERO

      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CMGSO4  = ZERO
      CK2SO4  = ZERO
      CCASO4  = MIN (WAER(6), WAER(2))

      END SUBROUTINE CALCV7



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCV1A
C *** CASE V1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCV1A

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: SO4FR, CAFR, FRK, NAFR, FRMG, FRNH3
C
C *** CALCULATE SOLIDS **************************************************
C
      CCASO4  = MIN (WAER(6), WAER(2))                     ! CCASO4
      SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
      CAFR    = MAX (WAER(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)                 ! CK2SO4
      FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*WAER(1), SO4FR)                 ! CNA2SO4
      NAFR    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CNA2SO4, ZERO)
      CMGSO4  = MIN (WAER(8), SO4FR)                       ! CMGSO4
      FRMG    = MAX (WAER(8) - CMGSO4, ZERO)
      SO4FR   = MAX (SO4FR - CMGSO4, ZERO)
      CNH42S4 = MAX (MIN (SO4FR , 0.5d0*WAER(3)) , TINY)
      FRNH3   = MAX (WAER(3) - 2.D0*CNH42S4, ZERO)
      CNH4NO3 = MIN (FRNH3, WAER(4))
      FRNH3   = MAX (FRNH3 - CNH4NO3, ZERO)
      CNH4CL  = MIN (FRNH3, WAER(5))
      FRNH3   = MAX (FRNH3 - CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      END SUBROUTINE CALCV1A




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCU8
C *** CASE U8
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCU8
      use isrpia
      implicit none
C
      real    :: NH4I, NAI, NO3I, CAI, KI, MGI
      real    :: gg, sdd, hi, ohi, so4i, cli, hso4i, akw, del, rtw
      real    :: hi0, ohi0, a2, a3i, a4i, molalr(21)
      integer :: i
C
C
C *** SETUP PARAMETERS ************************************************
C
      CALL CALCU1A

      FRST   = .TRUE.
      CALAIN = .TRUE.
      CALAOU = .TRUE.
C
C *** CALCULATE WATER **************************************************
C
      MOLALR(1) = CNACL         ! NACL
      MOLALR(2) = CNA2SO4       ! NA2SO4
      MOLALR(3) = CNANO3        ! NANO3
      MOLALR(5) = CNH4NO3       ! NH4NO3
      MOLALR(6) = CNH4CL        ! NH4CL
      MOLALR(17)= CK2SO4        ! K2SO4
      MOLALR(21)= CMGSO4        ! MGSO4

      WATER  = SUM(MOLALR(1:3)*M0I(1:3)) + SUM(MOLALR(5:6)*M0I(5:6))
     &       + MOLALR(17)*M0I(17) + MOLALR(21)*M0I(21)

      RTW = R*TEMP*WATER
C
C *** SETUP LIQUID CONCENTRATIONS **************************************
C
      AKW = SQRT(XKW*RH)*WATER                       ! H2O    <==> H+
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP

      NAI    = WAER(1)
      SO4I   = MAX(WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = 0.0
      KI     = WAER(7)
      MGI    = WAER(8)
C
C SOLUTION ACIDIC OR BASIC?
C
      IF (I == 1) THEN
         GG   = 2.*SO4I + NO3I + CLI - NAI - NH4I - KI - 2.*MGI
         SDD = SQRT(GG*GG + 4.*AKW*AKW)
         SDD = MAX(SDD, ABS(GG), 2.*AKW)

         IF (GG > 0.0) THEN                        ! H+ in excess
            HI0 = 0.5 * (SDD + GG)
            OHI0= (AKW/HI0) * AKW
         ELSE                                      ! OH- in excess
            OHI0= 0.5 * (SDD - GG)
            HI0 = (AKW/OHI0) * AKW
         ENDIF
      ENDIF

      HI  = HI0
      OHI = OHI0
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.GT.OHI) THEN
C         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
C         HI    = AKW/OHI
C         HSO4I = ZERO
C      ELSE
C         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
C     &           - KI - 2.D0*MGI, ZERO)
C         GGCL  = MAX(GG-GGNO3, ZERO)
C         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
C         IF (GGNO3.GT.TINY) THEN
C            IF (GGCL.LE.TINY) HI = ZERO
C            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
C         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
         CALL CALCHS4 (HI, SO4I, 0.0, DEL)
      ELSE
        DEL= ZERO
      ENDIF

      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL
      HI    = MAX(HI,AKW)
      OHI   = (AKW/HI) * AKW
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         EXIT
      ENDIF

      ENDDO
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2       = XK2*R*TEMP*GAMA(10)**2/(XKW*GAMA(5)**2) ! NH3  <==> NH4+
      A3I      = GAMA(10) / (RTW*XK4*WATER)        ! HNO3 <==> NO3-
      A4I      = GAMA(11) / (RTW*XK3*WATER)        ! HCL  <==> CL-

      GNH3     = NH4I/(HI*A2)
      GHNO3    = HI*NO3I*A3I
      GHCL     = HI*CLI *A4I

      GASAQ(1) = ZERO
      GASAQ(2) = ZERO
      GASAQ(3) = ZERO

      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = ZERO
      CMGSO4   = ZERO
      CK2SO4   = ZERO
      CCASO4   = MIN (WAER(6), WAER(2))

      END SUBROUTINE CALCU8



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCU1A
C *** CASE U1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
C     2. THERE IS ONLY A SOLID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCU1A

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: SO4FR, CAFR, FRK, FRMG, FRNA, FRNO3, FRCL, FRNH3
C
C *** CALCULATE SOLIDS *************************************************
C
      CCASO4  = MIN (WAER(6), WAER(2))                 ! CCASO4
      SO4FR   = MAX(WAER(2) - CCASO4, ZERO)
      CAFR    = MAX(WAER(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)             ! CK2SO4
      FRK     = MAX(WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
      CMGSO4  = MIN (WAER(8), SO4FR)                   ! CMGSO4
      FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
      CNA2SO4 = MAX (SO4FR, ZERO)                      ! CNA2SO4
      FRNA    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)

      CNH42S4 = ZERO

      CNANO3  = MIN (FRNA, WAER(4))
      FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
      FRNA    = MAX (FRNA-CNANO3, ZERO)

      CNACL   = MIN (FRNA, WAER(5))
      FRCL    = MAX (WAER(5)-CNACL, ZERO)
      FRNA    = MAX (FRNA-CNACL, ZERO)

      CNH4NO3 = MIN (FRNO3, WAER(3))
      FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
      FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)

      CNH4CL  = MIN (FRCL, FRNH3)
      FRCL    = MAX (FRCL-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      END SUBROUTINE CALCU1A
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCW13
C *** CASE W13
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
C                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCW13

      USE ISRPIA
      IMPLICIT NONE

      real    :: NH4I, NAI, NO3I, CAI, KI, MGI
      real    :: gg, sdd, hi, ohi, so4i, cli, hso4i, akw, del, rtw
      real    :: hi0, ohi0, a2, a3i, a4i, molalr(23)
      integer :: i
C
C *** SETUP PARAMETERS ************************************************
C
      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
      CALL CALCW1A

      MOLALR(1) = CNACL         ! NACL
      MOLALR(3) = CNANO3        ! NANO3
      MOLALR(5) = CNH4NO3       ! NH4NO3
      MOLALR(6) = CNH4CL        ! NH4CL
      MOLALR(15)= CCANO32       ! CANO32
      MOLALR(16)= CCACL2        ! CACL2
      MOLALR(17)= CK2SO4        ! K2SO4
      MOLALR(19)= CKNO3         ! KNO3
      MOLALR(20)= CKCL          ! KCL
      MOLALR(21)= CMGSO4        ! MGSO4
      MOLALR(22)= CMGNO32       ! MGNO32
      MOLALR(23)= CMGCL2        ! MGCL2

      WATER  = SUM(MOLALR(05:06)*M0I(05:06)) + SUM(MOLALR(15:17)*M0I(15:17))
     &       + SUM(MOLALR(19:23)*M0I(19:23))
     &       + MOLALR(1)*M0I(1) + MOLALR(3)*M0I(3)

      RTW = R*TEMP*WATER
      AKW = SQRT(XKW*RH)*WATER                    ! H2O       <==> H+
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
C ION CONCENTRATIONS
C
      NAI    = WAER(1)
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = 0.0
      KI     = WAER(7)
      MGI    = WAER(8)
C
C SOLUTION ACIDIC OR BASIC?
C
      IF (I == 1) THEN
         GG   = 2.*SO4I + NO3I + CLI - NAI - NH4I - KI - 2.*MGI
         SDD = SQRT(GG*GG + 4.*AKW*AKW)
         SDD = MAX(SDD, ABS(GG), 2.*AKW)

         IF (GG > 0.0) THEN                        ! H+ in excess
            HI0 = 0.5 * (SDD + GG)
            OHI0= (AKW/HI0) * AKW
         ELSE                                      ! OH- in excess
            OHI0= 0.5 * (SDD - GG)
            HI0 = (AKW/OHI0) * AKW
         ENDIF
      ENDIF

      HI  = HI0
      OHI = OHI0
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.GT.OHI) THEN
C         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
C         HI    = AKW/OHI
C         HSO4I = ZERO
C      ELSE
C         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
C     &           - KI - 2.D0*MGI, ZERO)
C         GGCL  = MAX(GG-GGNO3, ZERO)
C         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
C         IF (GGNO3.GT.TINY) THEN
C            IF (GGCL.LE.TINY) HI = ZERO
C            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
C         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, 0.0, DEL)
      ELSE
        DEL= ZERO
      ENDIF

      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL
      HI    = MAX(HI,AKW)
      OHI   = (AKW/HI) * AKW
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         EXIT
      ENDIF

      ENDDO
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      A2       = XK2*R*TEMP*GAMA(10)**2/(XKW*GAMA(5)**2) ! NH3  <==> NH4+
      A3I      = GAMA(10) / (RTW*XK4*WATER)        ! HNO3 <==> NO3-
      A4I      = GAMA(11) / (RTW*XK3*WATER)        ! HCL  <==> CL-

      GNH3    = NH4I/(HI*A2)
      GHNO3   = HI*NO3I*A3I
      GHCL    = HI*CLI *A4I

      GASAQ(1)= ZERO
      GASAQ(2)= ZERO
      GASAQ(3)= ZERO

      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CMGSO4  = ZERO
      CK2SO4  = ZERO
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = ZERO
      CKCL    = ZERO
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO

      END SUBROUTINE CALCW13



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCW1A
C *** CASE W1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
C                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCW1A

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: CAFR, SO4FR, FRK, FRMG, FRNA, CLFR, FRNO3
C
C *** CALCULATE SOLIDS **************************************************
C
      CCASO4  = MIN (WAER(2), WAER(6))              !SOLID CASO4
      CAFR    = MAX (WAER(6) - CCASO4, ZERO)
      SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
      CK2SO4  = MIN (SO4FR, 0.5D0*WAER(7))          !SOLID K2SO4
      FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CMGSO4  = SO4FR                               !SOLID MGSO4
      FRMG    = MAX (WAER(8) - CMGSO4, ZERO)
      CNACL   = MIN (WAER(1), WAER(5))              !SOLID NACL
      FRNA    = MAX (WAER(1) - CNACL, ZERO)
      CLFR    = MAX (WAER(5) - CNACL, ZERO)
      CCACL2  = MIN (CAFR, 0.5D0*CLFR)              !SOLID CACL2
      CAFR    = MAX (CAFR - CCACL2, ZERO)
      CLFR    = MAX (WAER(5) - 2.D0*CCACL2, ZERO)
      CCANO32 = MIN (CAFR, 0.5D0*WAER(4))           !SOLID CA(NO3)2
      CAFR    = MAX (CAFR - CCANO32, ZERO)
      FRNO3   = MAX (WAER(4) - 2.D0*CCANO32, ZERO)
      CMGCL2  = MIN (FRMG, 0.5D0*CLFR)              !SOLID MGCL2
      FRMG    = MAX (FRMG - CMGCL2, ZERO)
      CLFR    = MAX (CLFR - 2.D0*CMGCL2, ZERO)
      CMGNO32 = MIN (FRMG, 0.5D0*FRNO3)             !SOLID MG(NO3)2
      FRMG    = MAX (FRMG - CMGNO32, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CMGNO32, ZERO)
      CNANO3  = MIN (FRNA, FRNO3)                   !SOLID NANO3
      FRNA    = MAX (FRNA - CNANO3, ZERO)
      FRNO3   = MAX (FRNO3 - CNANO3, ZERO)
      CKCL    = MIN (FRK, CLFR)                     !SOLID KCL
      FRK     = MAX (FRK - CKCL, ZERO)
      CLFR    = MAX (CLFR - CKCL, ZERO)
      CKNO3   = MIN (FRK, FRNO3)                    !SOLID KNO3
      FRK     = MAX (FRK - CKNO3, ZERO)
      FRNO3   = MAX (FRNO3 - CKNO3, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO

      END SUBROUTINE CALCW1A
