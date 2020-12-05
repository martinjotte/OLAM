
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

C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS

C=======================================================================

      SUBROUTINE ISRP1F (WI, RHI, TEMPI)
      USE ISRPIA

      IMPLICIT NONE

      REAL(8) :: WI(NCOMP)
      REAL    :: RHI, TEMPI
      REAL(8) :: DC

C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR
C
      IF (SULRAT.GE.2.D0) THEN
         DC   = W(3) - 2.001D0*W(2) ! For numerical stability
         W(3) = W(3) + MAX(-DC, ZERO)

         SCASE = IAS
         CALL CALCA2            ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (SULRAT.GE.1.D0 .AND. SULRAT.LT.2.D0) THEN

         SCASE = IBS
         CALL CALCB4                 ! Only liquid (metastable)
         CALL CALCNH3
C
C *** SULFATE RICH (FREE ACID)
C
      ELSE

         SCASE = ICS
         CALL CALCC2                 ! Only liquid (metastable)
         CALL CALCNH3

      ENDIF
C
      END SUBROUTINE ISRP1F


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP2F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE ISRP2F (WI, RHI, TEMPI)
      USE ISRPIA

      IMPLICIT NONE
      REAL(8) :: WI(NCOMP)
      REAL    :: RHI, TEMPI
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR
C
      IF (SULRAT.GT.2.D0) THEN
C
         SCASE = IDS
         CALL CALCD3                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
C     SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (SULRAT.GT.1.D0 .AND. SULRAT.LT.2.D0) THEN
C
         SCASE = IBS
         CALL CALCB4                 ! Only liquid (metastable)
         CALL CALCNA                 ! HNO3(g) DISSOLUTION
C
C *** SULFATE RICH (FREE ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES,
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
C     SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSE
C
         SCASE = ICS
         CALL CALCC2                 ! Only liquid (metastable)
         CALL CALCNA                 ! HNO3(g) DISSOLUTION

      ENDIF
C
      END SUBROUTINE ISRP2F


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP3F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE ISRP3F (WI, RHI, TEMPI)
      USE ISRPIA

      IMPLICIT NONE

      REAL(8) :: WI(NCOMP)
      REAL    :: RHI, TEMPI
      REAL(8) :: REST
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
C
      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
      ENDIF
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL ISOINIT3 (WI, RHI, TEMPI)
C
C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
C
      REST = 2.D0*W(2) + W(4) + W(5)
      IF (W(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
         W(1) = (ONE-1.D-6)*REST         ! Adjust Na amount
      ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

C *** SULFATE POOR ; SODIUM POOR
C
      IF (SULRAT.GE.2.D0 .AND. SODRAT.LT.2.D0) THEN
C
         SCASE = IGS
         CALL CALCG5                 ! Only liquid (metastable)
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.2.D0 .AND. SODRAT.GE.2.D0) THEN
C
         SCASE = IHS
         CALL CALCH6                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (SULRAT.GE.1.D0 .AND. SULRAT.LT.2.D0) THEN
C
         SCASE = IIS
         CALL CALCI6                 ! Only liquid (metastable)
         CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
         CALL CALCNH3                !                NH3
C
C *** SULFATE RICH (FREE ACID)
C
      ELSE
C
         SCASE = IJS
         CALL CALCJ3                 ! Only liquid (metastable)
         CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
         CALL CALCNH3                !                NH3
      ENDIF
C
C *** END OF SUBROUTINE ISRP3F *****************************************
C
      END SUBROUTINE ISRP3F


C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE ISRP4F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTASSIUM-MAGNESIUM
C     AEROSOL SYSTEM.
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
      SUBROUTINE ISRP4F (WI, RHI, TEMPI)
      USE ISRPIA

      IMPLICIT NONE

      REAL(8) :: WI(NCOMP)
      REAL    :: RHI, TEMPI
      REAL(8) :: NAFRI, NO3FRI, REST, REST1, REST2, REST3, CCASO4I
      REAL(8) :: CASO4I, FRSO4I, CAFRI, CCANO32I, CCACL2I, CLFRI, CMGCL2I
      REAL(8) :: CNA2SO4I, CNACLI, CNANO3I, NO3FR, CMGSO4I, FRMGI, CMGNO32I

C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
C      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
C      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
C
C      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
C         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
C         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
C      ENDIF
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT4 (WI, RHI, TEMPI)
C
C *** CHECK IF TOO MUCH SODIUM+CRUSTALS ; ADJUST AND ISSUE ERROR MESSAGE
C
      REST = 2.D0*W(2) + W(4) + W(5)
C
      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
C
      CCASO4I  = MIN (W(2),W(6))
      FRSO4I   = MAX (W(2) - CCASO4I, ZERO)
      CAFRI    = MAX (W(6) - CCASO4I, ZERO)
      CCANO32I = MIN (CAFRI, 0.5D0*W(4))
      CAFRI    = MAX (CAFRI - CCANO32I, ZERO)
      NO3FRI   = MAX (W(4) - 2.D0*CCANO32I, ZERO)
      CCACL2I  = MIN (CAFRI, 0.5D0*W(5))
      CLFRI    = MAX (W(5) - 2.D0*CCACL2I, ZERO)
      REST1    = 2.D0*FRSO4I + NO3FRI + CLFRI
C
      CNA2SO4I = MIN (FRSO4I, 0.5D0*W(1))
      FRSO4I   = MAX (FRSO4I - CNA2SO4I, ZERO)
      NAFRI    = MAX (W(1) - 2.D0*CNA2SO4I, ZERO)
      CNACLI   = MIN (NAFRI, CLFRI)
      NAFRI    = MAX (NAFRI - CNACLI, ZERO)
      CLFRI    = MAX (CLFRI - CNACLI, ZERO)
      CNANO3I  = MIN (NAFRI, NO3FRI)
      NO3FR    = MAX (NO3FRI - CNANO3I, ZERO)
      REST2    = 2.D0*FRSO4I + NO3FRI + CLFRI
C
      CMGSO4I  = MIN (FRSO4I, W(8))
      FRMGI    = MAX (W(8) - CMGSO4I, ZERO)
      FRSO4I   = MAX (FRSO4I - CMGSO4I, ZERO)
      CMGNO32I = MIN (FRMGI, 0.5D0*NO3FRI)
      FRMGI    = MAX (FRMGI - CMGNO32I, ZERO)
      NO3FRI   = MAX (NO3FRI - 2.D0*CMGNO32I, ZERO)
      CMGCL2I  = MIN (FRMGI, 0.5D0*CLFRI)
      CLFRI    = MAX (CLFRI - 2.D0*CMGCL2I, ZERO)
      REST3    = 2.D0*FRSO4I + NO3FRI + CLFRI
C
         IF (W(6).GT.REST) THEN                       ! Ca > 2*SO4+CL+NO3 ?
             W(6) = (ONE-1D-6)*REST              ! Adjust Ca amount
             W(1)= ZERO                          ! Adjust Na amount
             W(7)= ZERO                          ! Adjust K amount
             W(8)= ZERO                          ! Adjust Mg amount
!            CALL PUSHERR (0051, 'ISRP4F')       ! Warning error: Ca, Na, K, Mg in excess
C
         ELSE IF (W(1).GT.REST1) THEN                 ! Na > 2*FRSO4+FRCL+FRNO3 ?
             W(1) = (ONE-1D-6)*REST1             ! Adjust Na amount
             W(7)= ZERO                          ! Adjust K amount
             W(8)= ZERO                          ! Adjust Mg amount
!            CALL PUSHERR (0052, 'ISRP4F')       ! Warning error: Na, K, Mg in excess
C
         ELSE IF (W(8).GT.REST2) THEN                 ! Mg > 2*FRSO4+FRCL+FRNO3 ?
             W(8) = (ONE-1D-6)*REST2             ! Adjust Mg amount
             W(7)= ZERO                          ! Adjust K amount
!            CALL PUSHERR (0053, 'ISRP4F')       ! Warning error: K, Mg in excess
C
         ELSE IF (W(7).GT.REST3) THEN                 ! K > 2*FRSO4+FRCL+FRNO3 ?
             W(7) = (ONE-1D-6)*REST3             ! Adjust K amount
!            CALL PUSHERR (0054, 'ISRP4F')       ! Warning error: K in excess
         ENDIF
      ENDIF
C
C *** CALCULATE RATIOS *************************************************
C
      SO4RAT  = (W(1)+W(3)+W(6)+W(7)+W(8))/W(2)
      CRNARAT = (W(1)+W(6)+W(7)+W(8))/W(2)
      CRRAT   = (W(6)+W(7)+W(8))/W(2)
C
C *** FIND CALCULATION REGIME FROM (SO4RAT, CRNARAT, CRRAT, RRH) ********
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) POOR: R(Cr+Na)<2
C
      IF (SO4RAT.GE.2.D0 .AND. CRNARAT.LT.2.D0) THEN
C
         SCASE = IOS
         CALL CALCO7                 ! Only liquid (metastable)
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
      ELSEIF (SO4RAT.GE.2.D0 .AND. CRNARAT.GE.2.D0) THEN
C
       IF (CRRAT.LE.2.D0) THEN
C
          SCASE = IMS
          CALL CALCM8                 ! Only liquid (metastable)
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
       ELSE
C
          SCASE = IPS
          CALL CALCP13                 ! Only liquid (metastable)

       ENDIF
C
C *** SULFATE RICH (NO ACID): 1<Rso4<2;
C
      ELSEIF (SO4RAT.GE.1.D0 .AND. SO4RAT.LT.2.D0) THEN
C
         SCASE = ILS
         CALL CALCL9                 ! Only liquid (metastable)
         CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
         CALL CALCNH3                !                NH3
C
C *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
C
      ELSE
C
         SCASE = IKS
         CALL CALCK4                 ! Only liquid (metastable)
         CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
         CALL CALCNH3                !                NH3
C
      ENDIF
C
      END SUBROUTINE ISRP4F
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCA2
C *** CASE A2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT >= 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS X, THE
C     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
C     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
C     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
C     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCA2
      USE ISRPIA

      IMPLICIT NONE
      REAL(8) :: OMELO, OMEHI, X1, Y1, ABS, X2, Y2, X3, Y3, A20
      REAL    :: DYM
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU    =.TRUE.       ! Outer loop activity calculation flag
      OMELO     = TINY        ! Low  limit: SOLUTION IS VERY BASIC
      OMEHI     = 2.0D0*W(2)  ! High limit: FROM NH4+ -> NH3(g) + H+(aq)

      A20 = XK2 * R * TEMP / XKW
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLAL(5) = W(2)
      WATER    = W(2) * M0I(4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = OMELO
      CALL RSTGAMP
      Y1 = FUNCA2 (X1)
      IF (ABS(Y1).LE.1.e-30) RETURN

      X2 = OMEHI
      CALL RSTGAMP
      Y2 = FUNCA2 (X2)
      IF (ABS(Y2).LE.1.e-30) RETURN

      IF (Y1 * Y2 > ZERO) THEN
         IF (ABS(Y1) < ABS(Y2)) THEN
            CALL RSTGAMP
            Y1 = FUNCA2 (X1)
         ENDIF
         RETURN
      ENDIF

      DO I = 1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3*X3 + FX12*(X1+X2)

         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCA2 (X3)

         if (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*X3) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF
      ENDDO


      CONTAINS


      REAL(8) FUNCTION FUNCA2 (OMEGI)

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: OMEGI
      INTEGER             :: I
      REAL(8)             :: LAMDA, NH4I, SO4I, A1, A2, A3, ZETA, PSI
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI    = W(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
         A1    = XK1 * WATER * GAMA(8)**2 / GAMA(7)**3
         A2    = A20 * (GAMA(8) / GAMA(9))**2
         A3    = XKW * RH * REAL(WATER,8)**2
C
         LAMDA = PSI * OMEGI / (A1 + OMEGI)
         ZETA  = A3 / OMEGI

         SO4I = MAX(PSI-LAMDA,TINY)
         NH4I = MAX( W(3)*A2*OMEGI/(ONE+A2*OMEGI), 2.D0*SO4I )
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = OMEGI                                        ! HI
         MOLAL (3) = NH4I                                         ! NH4I
         MOLAL (5) = SO4I                                         ! SO4I
         MOLAL (6) = LAMDA                                        ! HSO4I
         GNH3      = MAX (W(3)-NH4I, TINY)                        ! NH3GI
         COH       = ZETA                                         ! OHI
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
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
      FUNCA2 = OMEGI + NH4I - 2.D0*SO4I - LAMDA
C
      END FUNCTION FUNCA2

      END SUBROUTINE CALCA2


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB4
C *** CASE B4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
C     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
C     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCB4

      USE ISRPIA
      IMPLICIT NONE

      INTEGER :: I
      REAL(8) :: AK1, BET, GAM, BB, CC, DD, SO4I, HSO4I
      REAL    :: MOLALR(13)
C
C *** SOLVE EQUATIONS **************************************************
C
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.
C
C *** CALCULATE WATER CONTENT ******************************************
C
      CALL CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.

      MOLALR(13) = CLC
      MOLALR(9)  = CNH4HS4
      MOLALR(4)  = CNH42S4
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)*M0I(13) + MOLALR(9)*M0I(9) + MOLALR(4)*M0I(4)

      MOLAL(3)   = W(3)   ! NH4I

      DO I = 1, NSWEEP
         AK1   = XK1 * WATER * GAMA(8)**2 / GAMA(7)**3
         BET   = W(2)
         GAM   = MOLAL(3)

         BB    = GAM - BET - AK1
         CC    = AK1*BET
         DD    = BB*BB + 4.D0*CC
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (5) = MAX(TINY,MIN(0.5D0*(BB + SQRT(DD)), W(2))) ! SO4I
         MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))          ! HSO4I
         MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2)))  ! HI

         SO4I  = MOLAL(5) - MOLAL(1)     ! CORRECT FOR HSO4 DISSOCIATION
         HSO4I = MOLAL(6) + MOLAL(1)

         IF (SO4I.LT.HSO4I) THEN
            MOLALR(13) = SO4I                   ! [LC] = [SO4]
            MOLALR(9)  = MAX(HSO4I-SO4I, ZERO)  ! NH4HSO4
            MOLALR(4)  = 0.0
         ELSE
            MOLALR(13) = HSO4I                  ! [LC] = [HSO4]
            MOLALR(4)  = MAX(SO4I-HSO4I, ZERO)  ! (NH4)2SO4
            MOLALR(9)  = 0.0
         ENDIF

         WATER = MOLALR(13)*M0I(13) + MOLALR(9)*M0I(9) + MOLALR(4)*M0I(4)
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) EXIT
         CALL CALCACT
      ENDDO

      END SUBROUTINE CALCB4


C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB1A
C *** CASE B1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH
C     2. THERE IS NO LIQUID PHASE
C     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
C                         BUT NOT BOTH)
C
C     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
C     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
C     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT
C     IS MIXED WITH THE LC.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCB1A

      USE ISRPIA
      IMPLICIT NONE
      REAL(8) :: X, Y
C
C *** SETUP PARAMETERS ************************************************
C
      X = 2.D0*W(2)-W(3)    ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4
C
C *** CALCULATE COMPOSITION *******************************************
C
      IF (X.LE.Y) THEN      ! LC is the MIN (x,y)
         CLC     = X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = Y-X
      ELSE
         CLC     = Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = X-Y
         CNH42S4 = ZERO
      ENDIF
C
      END SUBROUTINE CALCB1A



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCC2
C *** CASE C2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================

      SUBROUTINE CALCC2

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: LAMDA, KAPA, PSI, PARM, BB, CC
      REAL    :: MOLALR(7:9)
      INTEGER :: I

      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** SOLVE EQUATIONS **************************************************
C
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION

      MOLAL (3) = LAMDA       ! NH4I
      MOLALR(9) = LAMDA       ! NH4HSO4
      MOLALR(7) = PSI         ! H2SO4
      WATER = MOLALR(9)*M0I(9) + MOLALR(7)*M0I(7)

      DO I = 1, NSWEEP
         PARM  = WATER * XK1 * GAMA(8)**2 / GAMA(7)**3
         BB    = PSI + PARM
         CC    = PARM * (LAMDA+PSI)
         KAPA  = 0.5D0 * (SQRT(BB*BB + 4.D0*CC) - BB)
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = PSI+KAPA                               ! HI
         MOLAL (5) = KAPA                                   ! SO4I
         MOLAL (6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) EXIT
         CALL CALCACT
      ENDDO

      END SUBROUTINE CALCC2



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD3
C *** CASE D3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS OLNY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCD3

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: X1, X2, X3, Y1, Y2, Y3, CHI1, CHI2, CHI3, CHI4
      REAL(8) :: PSI4LO, PSI4HI, P4, YY, A40, A340
      REAL    :: DELTA, DYM, MOLALR(4:5)
      INTEGER :: I, NCHECK
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL CALCD1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3

      MOLAL(5) = CHI2              ! Include initial amount in water calc
      MOLAL(3) = CHI1
      MOLAL(7) = CHI1

      MOLALR(4) = CHI2                        ! (NH4)2SO4
      MOLALR(5) = MAX(CHI1 - 2.D0*CHI2, ZERO) ! NH4NO3 = MIN("free", NO3)

      WATER = MOLALR(4)*M0I(4) + MOLALR(5)*M0I(5)

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit

      CNH4NO3   = ZERO             ! Solid NH4NO3
      CNH42S4   = ZERO             ! Solid (NH4)2SO4

      NCHECK = 0

      A40                   = XK2 * R * TEMP / XKW
      IF (CHI3 > TINY) A340 = XK4 * R * TEMP * A40
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
 60   X1 = PSI4LO
      CALL RSTGAMP
      Y1 = FUNCD3 (X1)
      IF (ABS(Y1) < 1.D-30) GOTO 50

      X2 = PSI4HI
      CALL RSTGAMP
      Y2 = FUNCD3 (X2)
      IF (ABS(Y2) < 1.D-30) GOTO 50

      NCHECK = NCHECK + 1

C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
C Physically I dont know when this might happen, but I have put this
C branch in for completeness. I assume there is no solution; all NO3 goes to the
C gas phase.

      IF (Y1.LT.ZERO .AND. Y2.LT.ZERO) THEN

         P4 = TINY ! PSI4LO ! CHI4
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCD3(P4)
         GOTO 50

C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
C This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
C and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
C and proceed again with root tracking.
C
      ELSE IF (Y1.GT.ZERO .AND. Y2.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(CHI1+CHI2) ! No solution; some NH3 evaporates
         IF (PSI4LO.LT.-(CHI1+CHI2) .or. ncheck > 2) THEN
!            CALL PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
            RETURN
         ELSE
            MOLAL(5) = CHI2              ! Include sulfate in initial water calculation
            MOLAL(3) = CHI1
            MOLAL(7) = CHI1

            MOLALR(4) = CHI2
            MOLALR(5) = MAX(CHI1 - 2.D0*CHI2, ZERO)
            WATER     = MOLALR(4)*M0I(4) + MOLALR(5)*M0I(5)
            GOTO 60                      ! Redo root tracking
         ENDIF
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO I = 1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3*X3 + FX12*(X1+X2)

         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCD3 (X3)

         if (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*ABS(X3)) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF
      ENDDO

 50   IF (MOLAL(1).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), 0.0, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF


      CONTAINS


      REAL(8) FUNCTION FUNCD3 (P4)
      IMPLICIT NONE

      REAL(8), INTENT(IN) :: P4
      INTEGER :: I
      REAL(8) :: PSI4, A34, A4, A7, PSI3, BB, AHI, NH4I
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4

      GNH3   = CHI4 - PSI4                  ! Gas NH3
      NH4I   = CHI1 + PSI4 + 2.D0*CHI2
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP

         IF (CHI3 > TINY) THEN
            A34  = A340 * GNH3 * (WATER / GAMA(5))**2
            PSI3 = (A34*CHI3 - CHI1*NH4I) / (A34 + NH4I)
            PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
         ELSE
            PSI3 = CHI3
         ENDIF

         A7 = XKW * RH * REAL(WATER,8)**2
         BB = PSI4 - PSI3

         IF (BB < 0.0) THEN
            AHI = 0.5D0   * (SQRT(BB*BB + 4.D0*A7) - BB)
         ELSE
            AHI = 2.D0*A7 / (SQRT(BB*BB + 4.D0*A7) + BB)
         ENDIF
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = NH4I                            ! NH4I
         MOLAL (7) = PSI3 + CHI1                     ! NO3I
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         MOLALR(5) = MIN(CHI1 + PSI4, PSI3 + CHI1)
         WATER     = MOLALR(4)*M0I(4) + MOLALR(5)*M0I(5)
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
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
      A4     = A40 * (GAMA(10) / GAMA(5))**2
      FUNCD3 = NH4I - AHI * GNH3 * A4
C
      END FUNCTION FUNCD3

      END SUBROUTINE CALCD3



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD1A
C *** CASE D1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
C     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
C     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
C     THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCD1A

      USE ISRPIA
      IMPLICIT NONE
      REAL(8) :: PARM, X, PS, OM, OMPS, DIAK, ZE
C
C *** SETUP PARAMETERS ************************************************
C
      PARM    = XK10 / (R*TEMP)**2
C
C *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
C
      CNH42S4 = W(2)
      X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  ! MAX NH4NO3
      PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)
C
      OMPS    = OM + PS
      DIAK    = SQRT(OMPS*OMPS + 4.D0*PARM)             ! DIAKRINOUSA
      ZE      = MIN(X, 0.5D0*(DIAK - OMPS))             ! THETIKI RIZA
C
C *** SPECIATION *******************************************************
C
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3
C
      END SUBROUTINE CALCD1A



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCG5
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCG5

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: X1, X2, X3, Y1, Y2, Y3, PSI6LO, PSI6HI
      REAL(8) :: CHI1, CHI2, CHI4, CHI5, CHI6
      REAL(8) :: A650, A4M0, A60
      REAL    :: DELTA, DYM, WATER0, MOLALR(6)
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.
      CHI1   = 0.5D0*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)

      MOLAL (2) = 2.0D0*CHI1                     ! NAI
      MOLAL (5) = CHI2 + CHI1                    ! SO4I

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      MOLALR(2) = CHI1          ! NA2SO4
      MOLALR(4) = CHI2          ! (NH4)2SO4
      WATER     = MOLALR(2)*M0I(2) + MOLALR(4)*M0I(4)
      WATER0    = WATER

      IF (CHI5 > TINY) A650 = XK3 / XK4
      IF (w(2) > TINY) A4M0 = XKW / (XK2 * R * TEMP)
      A60 = XK3*R*TEMP
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCG5A (X1)
      IF (CHI6 <= TINY .OR. ABS(Y1) <= 1.D-30) GOTO 50

      X2 = PSI6HI
      CALL RSTGAMP
      Y2 = FUNCG5A (X2)
      IF (ABS(Y2) <= 1.E-30) GOTO 50

      IF (Y2*Y1 > ZERO) THEN
         IF (ABS(Y1) < ABS(Y2)) THEN
            CALL RSTGAM
            Y1 = FUNCG5A (X1)
         ENDIF
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      CALL RSTGAM

      DO I = 1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3*X3 + FX12*(X1+X2)

         Y3 = FUNCG5A (X3)

         if (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*X3) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF
      ENDDO
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
 50   IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN ! If quadrat.called
         CALL CALCHS4 (MOLAL(1), MOLAL(5), 0.0, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
         MOLAL(6) = DELTA                               ! HSO4 EFFECT
      ENDIF


      CONTAINS


      REAL(8) FUNCTION FUNCG5A (PSI6)
      IMPLICIT NONE

      REAL(8), INTENT(IN) :: PSI6
      REAL(8)             :: PSI4, PSI5, BB, CC, DD, A6_A5, A6, A4M, FRNH4
      REAL                :: SMIN, HI, OHI, WATER0
      INTEGER             :: I
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.

      GHCL     = MAX(CHI6 - PSI6, TINY)                 ! Gas HCl
      MOLAL(4) = PSI6                                 ! CLI
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5>TINY) THEN
         A6_A5= A650 * (GAMA(10) / GAMA(11))**2
         PSI5 = PSI6*CHI5 / (A6_A5*GHCL + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
      IF(W(2).GT.TINY) THEN                     ! Accounts for NH3 evaporation
         A4M  = A4M0 * (GAMA(5) / GAMA(10))**2
         BB   = (CHI4 + PSI6 + PSI5 + A4M)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*CHI2*A4M
         DD   = MAX(BB*BB-4.d0*CC,ZERO)         ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 = 0.5d0*(BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (3) = 2.0D0*CHI2 + PSI4                   ! NH4I
      MOLAL (7) = PSI5                                ! NO3I

      SMIN = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3

      MOLALR(5) = MIN(PSI4, PSI5)                     ! NH4NO3
      FRNH4     = MAX(PSI4 - PSI5, ZERO)
      MOLALR(6) = MIN(PSI6, FRNH4)                    ! NH4CL

      WATER = WATER0 + MOLALR(5)*M0I(5) + MOLALR(6)*M0I(6)
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      A6      = A60 * (WATER/GAMA(11))**2
      FUNCG5A = HI * PSI6 - GHCL * A6

      END FUNCTION FUNCG5A

      END SUBROUTINE CALCG5



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCH6
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCH6

      USE ISRPIA
      IMPLICIT NONE

      INTEGER :: I
      REAL    :: DELTA, DYM, WATER0, MOLALR(6)
      REAL(8) :: FRNA, PSI6LO, PSI6HI, X1, Y1, X2, Y2, X3, Y3
      REAL(8) :: CHI1, CHI4, CHI5, CHI6, CHI7, CHI8, A640, A650, A4M0
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.
      CHI1   = W(2)                                ! CNA2SO4
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      MOLAL (2) = CHI8 + CHI7 + 2.D0*CHI1          ! NAI
      MOLAL (5) = CHI1                             ! SO4I

      MOLALR(1) = CHI7          ! NACL
      MOLALR(2) = CHI1          ! NA2SO4
      MOLALR(3) = CHI8          ! NANO3

      WATER     = SUM( MOLALR(1:3)*M0I(1:3) )
      WATER0    = WATER

      IF (CHI5 > TINY) A650 = XK3 / XK4
      IF (CHI4 > TINY) A4M0 = XKW / (XK2 * R * TEMP)
      A640 = XK3*XK2*(R*TEMP/XKW)**2
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCH6A (X1)
      IF (CHI6 <= TINY .OR. ABS(Y1) <= 1.D-30) GOTO 50

      X2 = PSI6HI
      CALL RSTGAMP
      Y2 = FUNCH6A (X2)
      IF (ABS(Y2) <= 1.e-30) GOTO 50

      IF (Y2*Y1 > ZERO) THEN
         IF (ABS(Y1) < ABS(Y2)) THEN
            CALL RSTGAM
            Y1 = FUNCH6A (X1)
         ENDIF
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      CALL RSTGAM

      DO I = 1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3*X3 + FX12*(X1+X2)

         Y3 = FUNCH6A (X3)

         if (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*X3) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF
      ENDDO
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
 50   IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), 0.0, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF



      CONTAINS



      REAL(8) FUNCTION FUNCH6A (PSI6)

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: PSI6
      REAL(8)             :: A4, A6_A5, PSI5, A4M, BB, CC, DD, PSI4, FRNH4, A64
      REAL                :: SMIN, HI, OHI, WATER0
      INTEGER             :: I
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      GHCL   = MAX(CHI6 - PSI6, TINY)
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5 > TINY) THEN
         A6_A5 = A650 * GHCL * (GAMA(10) / GAMA(11))**2
         PSI5 = CHI5*(PSI6+CHI7) - A6_A5*CHI8
         PSI5 = PSI5 / (A6_A5 + PSI6 + CHI7)
         PSI5 = MIN(MAX(PSI5, TINY),CHI5)
      ELSE
         PSI5 = CHI5
      ENDIF
C
      IF (CHI4.GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         A4M  = A4M0 * (GAMA(5) / GAMA(10))**2
         BB   = CHI4 + PSI6 + PSI5 + A4M
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 = 0.5d0*(BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (3) = PSI4                                  ! NH4I
      MOLAL (4) = PSI6 + CHI7                           ! CLI
      MOLAL (7) = PSI5 + CHI8                           ! NO3I

      SMIN      = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)

      MOLALR(5) = MIN(PSI4, PSI5)                  ! NH4NO3
      FRNH4     = MAX(PSI4 - PSI5, ZERO)           ! "FREE" NH3
      MOLALR(6) = MIN(PSI6, FRNH4)                 ! NH4CL

      WATER = WATER0 + SUM( MOLALR(5:6)*M0I(5:6) )
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      A64     = A640 * (WATER * GAMA(10) / GAMA(11) * GAMA(5))**2
      FUNCH6A = MOLAL(3)*MOLAL(4) - GHCL*GNH3*A64
      FUNCH6A = MAX(FUNCH6A, -CHI6 * CHI4)
C
      END FUNCTION FUNCH6A

      END SUBROUTINE CALCH6



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI6
C *** CASE I6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCI6

      USE ISRPIA

      IMPLICIT NONE
      REAL(8) :: A6, BB, CC, DD, PSI6
      REAL(8) :: CHI1, CHI2, CHI3, CHI4, CHI5
      REAL    :: MOLALR(13)
      INTEGER :: I
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL CALCI1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.

      MOLAL (2) = 2.D0*CHI4 + CHI3              ! NAI
      MOLAL (3) = 3.D0*CHI2 + 2.D0*CHI5 + CHI1  ! NH4I

      MOLALR(02) = CHI4         ! NA2SO4
      MOLALR(04) = CHI5         ! (NH4)2SO4
      MOLALR(09) = CHI1         ! NH4HSO4
      MOLALR(12) = CHI3         ! NAHSO4
      MOLALR(13) = CHI2         ! LC

      WATER = MOLALR(02)*M0I(02) + MOLALR(04)*M0I(04) + MOLALR(09)*M0I(09)
     &      + MOLALR(12)*M0I(12) + MOLALR(13)*M0I(13)

      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = ZERO
      CNH42S4   = ZERO
      CNH4HS4   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I=1,NSWEEP
         A6 = XK1 * WATER * GAMA(8)**2 / GAMA(7)**3
C
C  CALCULATE DISSOCIATION QUANTITIES
C
         BB   = CHI2 + CHI4 + CHI5 + A6
         CC   = A6*(CHI2 + CHI3 + CHI1)
         DD   = BB*BB + 4.D0*CC
         PSI6 = 0.5D0*(SQRT(DD) - BB)
C
C *** CALCULATE SPECIATION ********************************************
C
         MOLAL (1) = PSI6                                 ! HI
         MOLAL (5) = CHI2 + CHI4 + CHI5 + PSI6            ! SO4I
         MOLAL (6) = CHI2 + CHI3 + CHI1 - PSI6            ! HSO4I
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
      END SUBROUTINE CALCI6



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI1A
C *** CASE I1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCI1A

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: FRSO4, FRNH4
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CNA2SO4 = 0.5D0*W(1)
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      FRSO4   = MAX(W(2)-CNA2SO4, ZERO)
C
      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
C
      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
      ENDIF
C
C *** CALCULATE GAS SPECIES *********************************************
C
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO
C
C *** END OF SUBROUTINE CALCI1A *****************************************
C
      END SUBROUTINE CALCI1A



C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCJ3
C *** CASE J3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCJ3

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: LAMDA, KAPA, A3, BB, CC, DD, CHI1, CHI2
      REAL    :: MOLALR(7:12)
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      CHI1   = W(1)                           ! NA TOTAL as NaHSO4
      CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4

      MOLAL (2) = CHI1                        ! NAI
      MOLAL (3) = CHI2                        ! NH4I

      MOLALR(07) = LAMDA                      ! H2SO4
      MOLALR(09) = CHI2                       ! NH4HSO4
      MOLALR(12) = CHI1                       ! NAHSO4

      WATER = MOLALR(07)*M0I(07) + MOLALR(09)*M0I(09) + MOLALR(12)*M0I(12)
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
         A3 = XK1 * WATER * GAMA(8)**2 / GAMA(7)**3
C
C  CALCULATE DISSOCIATION QUANTITIES
C
         BB   = A3 + LAMDA
         CC   = A3 * (LAMDA + CHI1 + CHI2)
         DD   = BB*BB + 4.D0*CC
         KAPA = 0.5D0*(SQRT(DD) - BB)
C
C *** CALCULATE SPECIATION ********************************************
C
         MOLAL (1) = LAMDA + KAPA               ! HI
         MOLAL (5) = KAPA                       ! SO4I
         MOLAL (6) = LAMDA + CHI1 + CHI2 - KAPA ! HSO4I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         EXIT
      ENDIF

      ENDDO

      END SUBROUTINE CALCJ3



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCO7
C *** CASE O7
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCO7

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: SO4FR, CAFR, FRK, NAFR, FRMG, PSI6LO, PSI6HI
      REAL(8) :: X1, Y1, X2, Y2, X3, Y3, A650, A4M0, A60
      REAL(8) :: CHI1, CHI2, CHI4, CHI5, CHI6, CHI7, CHI8, CHI9
      REAL    :: DELTA, DYM, WATER0, MOLALR(21)
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     ! CCASO4
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              ! CK2SO4
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              ! NA2SO4
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    ! CMGSO4
      FRMG   = MAX(W(8) - CHI8, ZERO)
      SO4FR  = MAX(SO4FR - CHI8, ZERO)
      CCASO4 = CHI9

      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)

      PSI6LO = TINY
      PSI6HI = CHI6-TINY

      MOLAL (2) = 2.0D0*CHI1                       ! Na+
      MOLAL (5) = CHI1 + CHI2 + CHI7 + CHI8        ! SO4I
      MOLAL (9) = 2.0D0*CHI7                       ! KI
      MOLAL (10)= CHI8                             ! Mg

      MOLALR(2) = CHI1          ! NA2SO4
      MOLALR(4) = CHI2          ! (NH4)2SO4
      MOLALR(17)= CHI7          ! K2SO4
      MOLALR(21)= CHI8          ! MGSO4

      WATER  = MOLALR(02)*M0I(02) + MOLALR(04)*M0I(04) + MOLALR(17)*M0I(17)
     &       + MOLALR(21)*M0I(21)
      WATER0 = WATER

      IF (CHI5 > TINY) A650 = XK3 / XK4
      IF (CHI4 > TINY .AND. W(2)>TINY) A4M0 = XKW / (XK2 * R * TEMP)
      A60 = XK3*R*TEMP
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCO7 (X1)
      IF (CHI6 <= TINY .OR. abs(y1) <= 1.d-30) GOTO 50

      X2 = PSI6HI
      CALL RSTGAMP
      Y2 = FUNCO7 (X2)
      IF (ABS(Y2) <= 1.D-30) GOTO 50

      IF (Y1 * Y2 > ZERO) THEN
         IF (ABS(Y1) <= ABS(Y2)) THEN
            CALL RSTGAM
            Y1 = FUNCO7 (X1)
         ENDIF
         GOTO 50
      ENDIF

C
C *** PERFORM BISECTION ***********************************************
C
      CALL RSTGAM

      DO I=1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3*X3 + FX12*(X1+X2)

         Y3 = FUNCO7 (X3)

         IF (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*X3) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF

      ENDDO
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
 50   IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4( MOLAL(1), MOLAL(5), 0.0, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
         MOLAL(6) = DELTA                               ! HSO4 EFFECT
      ENDIF


      CONTAINS


      REAL(8) FUNCTION FUNCO7 (PSI6)

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: PSI6
      REAL(8)             :: A6_A5, PSI5, A4M, BB, CC, DD, PSI4, FRNH4, A6
      INTEGER             :: I
      REAL                :: HI, OHI, SMIN
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.

      GHCL     = MAX(CHI6 - PSI6, TINY)
      MOLAL(4) = PSI6                       ! CLI
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I=1,NSWEEP
C
      IF (CHI5.GT.TINY) THEN
         A6_A5= A650 * (GAMA(10) / GAMA(11))**2
         PSI5 = PSI6*CHI5 / (A6_A5*GHCL + PSI6)
         PSI5 = MIN (PSI5,CHI5)
      ELSE
         PSI5 = TINY
      ENDIF
C
      IF(W(2)>TINY .AND. CHI4>TINY) THEN  ! Accounts for NH3 evaporation
         A4M  = A4M0 * (GAMA(5) / GAMA(10))**2
         BB   = CHI4 + PSI6 + PSI5 + A4M
         CC   = CHI4*(PSI5+PSI6) - 2.d0*CHI2*A4M
         DD   = MAX(BB*BB-4.d0*CC, zero)   ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 = 0.5d0 * (BB - SQRT(DD))
         PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL (3) = 2.0D0*CHI2 + PSI4                ! NH4I
      MOLAL (7) = PSI5                             ! NO3I
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
      SMIN      = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
C
C *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
C
      MOLALR(5) = MIN(PSI4, PSI5)
      FRNH4     = MAX(PSI4 - PSI5, ZERO)
      MOLALR(6) = MIN(PSI6, FRNH4)                  ! NH4CL

      WATER = WATER0 + SUM( MOLALR(5:6)*M0I(5:6) )
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      A6     = A60 * (WATER/GAMA(11))**2
      FUNCO7 = HI * psi6 - GHCL * A6

      END FUNCTION FUNCO7

      END SUBROUTINE CALCO7



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCM8
C *** CASE M8
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCM8

      USE ISRPIA
      IMPLICIT NONE

      REAL    :: DELTA, DYM, WATER0, MOLALR(21)
      REAL(8) :: SO4FR, CAFR, FRK, FRMG, FRNA, PSI6LO, PSI6HI
      REAL(8) :: CHI1, CHI4, CHI5, CHI6, CHI7, CHI8, CHI9
      REAL(8) :: CHI10, CHI11, X1, Y1, X2, Y2, X3, Y3, A650, A4M0, A60
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    ! CCASO4
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             ! CK2S04
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  ! CMGSO4
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    ! CNA2SO4
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)

      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)

      CCASO4 = CHI11

      MOLAL (2) = CHI8 + CHI7 + 2.D0*CHI1          ! NAI
      MOLAL (5) = CHI1 + CHI9 + CHI10              ! SO4I
      MOLAL (9) = 2.D0*CHI9                        ! KI
      MOLAL (10)= CHI10                            ! MGI

      MOLALR(1) = CHI7                             ! NACL
      MOLALR(2) = CHI1                             ! NA2SO4
      MOLALR(3) = CHI8                             ! NANO3
      MOLALR(17)= CHI9                             ! K2SO4
      MOLALR(21)= CHI10                            ! MGSO4

      WATER  = SUM(MOLALR(1:3)*M0I(1:3)) + MOLALR(17)*M0I(17)
     &       + MOLALR(21)*M0I(21)
      WATER0 = WATER

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      IF (CHI5 > TINY) A650 = XK3 / XK4
      IF (CHI4 > TINY) A4M0 = XKW / (XK2 * R * TEMP)
      A60 = XK3 * R * TEMP
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCM8 (X1)
      IF (CHI6 <= TINY .OR. ABS(Y1) <= 1.D-30) GOTO 50

      X2 = PSI6HI
      CALL RSTGAMP
      Y2 = FUNCM8 (X2)
      IF (ABS(Y2) <= 1.D-30) GOTO 50

      IF (Y1 * Y2 > ZERO) THEN
         IF (ABS(Y1) <= ABS(Y2)) THEN
            CALL RSTGAM
            Y1 = FUNCM8 (X1)
         ENDIF
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      CALL RSTGAM

      DO I=1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3*X3 + FX12*(X1+X2)

         Y3 = FUNCM8 (X3)

         if (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*X3) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF

      ENDDO
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
 50   IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), 0.0, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF


      CONTAINS


      REAL(8) FUNCTION FUNCM8 (PSI6)
      IMPLICIT none

      REAL(8), INTENT(IN) :: PSI6
      REAL                :: SMIN, HI, OHI
      REAL(8)             :: CLI, A6_A5, PSI4, PSI5, A4M, BB, CC, DD, FRNH4, A6
      INTEGER             :: I
C
C *** SETUP PARAMETERS ************************************************
C
      FRST     = .TRUE.
      CALAIN   = .TRUE.
      CLI      = CHI7 + PSI6
      GHCL     = MAX(CHI6 - PSI6, TINY)
      MOLAL(4) = CLI
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5 > TINY) THEN
         A6_A5 = A650 * GHCL * (GAMA(10) / GAMA(11))**2
         PSI5 = CHI5*CLI - A6_A5*CHI8
         PSI5 = PSI5 / (A6_A5 + CLI)
         PSI5 = MIN(MAX(PSI5, TINY),CHI5)
      ELSE
         PSI5 = CHI5
      ENDIF

      IF (CHI4.GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         A4M  = A4M0 * (GAMA(5) / GAMA(10))**2
         BB   = CHI4 + PSI6 + PSI5 + A4M
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5D0*(BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      SMIN       = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)

      MOLAL (1) = HI
      MOLAL (3) = PSI4                                  ! NH4I
      MOLAL (7) = PSI5 + CHI8                           ! NO3I

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)

      MOLALR(5) = MIN(PSI4, PSI5)                       ! NH4NO3
      FRNH4     = MAX(PSI4 - PSI5, ZERO)                ! "FREE" NH3
      MOLALR(6) = MIN(PSI6, FRNH4)                         ! NH4CL

      WATER = WATER0 + SUM( MOLALR(5:6)*M0I(5:6) )
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      A6     = A60 * (WATER / GAMA(11))**2
      FUNCM8 = MOLAL(1)*MOLAL(4) - GHCL*A6

      END FUNCTION FUNCM8

      END SUBROUTINE CALCM8



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCP13
C *** CASE P13
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
      SUBROUTINE CALCP13

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: FRCA, FRSO4, FRK, FRMG, FRNA, FRCL, FRNO3
      REAL(8) :: CHI4, CHI5, CHI6, CHI7, CHI8, CHI9, CHI10, CHI11
      REAL(8) :: CHI12, CHI13, CHI14, CHI15, CHI16, CHI17
      REAL(8) :: X1, Y1, X2, Y2, X3, Y3, PSI6LO, PSI6HI
      REAL(8) :: A650, A4M0, A60
      REAL    :: DELTA, DYM, WATER0, MOLALR(23)
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    ! CCASO4
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             ! CK2SO4
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               ! CMGSO4
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    ! CNACL
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              ! CCANO32
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              ! CCACL2
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             ! CMGNO32
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              ! CMGCL2
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   ! CNANO3
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     ! CKCL
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    ! CKNO3
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)

      CHI5    = FRNO3                               ! HNO3(g)
      CHI6    = FRCL                                ! HCL(g)
      CHI4    = W(3)                                ! NH3(g)

      CCASO4    = CHI11
      MOLAL (2) = CHI8 + CHI7                                     ! NAI
      MOLAL (5) = CHI9 + CHI10                                    ! SO4I
      MOLAL (8) = CHI12 + CHI17                                   ! CAI
      MOLAL (9) = 2.D0*CHI9 + CHI13 + CHI14                       ! KI
      MOLAL (10)= CHI10 + CHI15 + CHI16                           ! MGI

      MOLALR(1) = CHI7          ! NACL
      MOLALR(3) = CHI8          ! NANO3
      MOLALR(15)= CHI12         ! CANO32
      MOLALR(16)= CHI17         ! CACL2
      MOLALR(17)= CHI9          ! K2SO4
      MOLALR(19)= CHI13         ! KNO3
      MOLALR(20)= CHI14         ! KCL
      MOLALR(21)= CHI10         ! MGSO4
      MOLALR(22)= CHI15         ! MGNO32
      MOLALR(23)= CHI16         ! MGCL2

      WATER  = MOLALR(1)*M0I(1) + MOLALR(3)*M0I(3) + SUM(MOLALR(15:23)*M0I(15:23))
      WATER0 = WATER

      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)

      IF (CHI5 > TINY) A650 = XK3 / XK4
      IF (CHI4 > TINY) A4M0 = XKW / (XK2 * R * TEMP)
      A60 = XK3*R*TEMP
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCP13 (X1)
      IF (CHI6 <= TINY .OR. ABS(Y1) <= 1.D-30) GOTO 50

      X2 = PSI6HI
      CALL RSTGAMP
      Y2 = FUNCP13 (X2)
      IF (ABS(Y2) <= 1.D-30) GOTO 50

      IF (Y1 * Y2 > ZERO) THEN
         IF (ABS(Y1) <= ABS(Y2)) THEN
            CALL RSTGAM
            Y1 = FUNCP13 (X1)
         ENDIF
         GOTO 50
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
      CALL RSTGAM

      DO I=1, MAXIT

         DYM = 1.0 / REAL(Y2 - Y1)
         X3 = (Y2 * X1 - Y1 * X2) * DYM
         X3 = FX3 *X3 + FX12*(X1+X2)

         Y3 = FUNCP13 (X3)

         IF (ABS(Y3) <= 1.D-30) EXIT
         IF (I>1 .AND. ABS(X3-X2) <= EPS*X3) EXIT

         IF (Y3 * Y2 < 0.D0) THEN

            ! ROOT WAS TRAPPED, SO USE REGULA FALSI
            X1 = X2
            Y1 = Y2
            X2 = X3
            Y2 = Y3

         ELSE

            ! ROOT WAS NOT TRAPPED, SO USE ILLINOIS MODIFICATION
            X2 = X3
            Y2 = Y3
            Y1 = 0.5 * Y1

         ENDIF

      ENDDO
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
 50   IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), 0.0, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF


      CONTAINS


      REAL(8) FUNCTION FUNCP13 (PSI6)

      IMPLICIT NONE
      REAL(8), INTENT(IN) :: PSI6
      REAL(8)             :: PSI4, PSI5, A4M, BB, CC, DD, FRNH4
      REAL(8)             :: A6, A6_A5
      REAL                :: SMIN, HI, OHI
      INTEGER             :: I
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      GHCL   = MAX(CHI6 - PSI6, TINY)

      MOLAL(4) = PSI6 + CHI7 + CHI14 + 2.D0*CHI16 + 2.D0*CHI17  ! CLI
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I=1,NSWEEP
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5 > TINY) THEN
         A6_A5= A650 * GHCL * (GAMA(10) / GAMA(11))**2
         PSI5 = CHI5*(PSI6 + CHI7 + CHI14 + 2.D0*CHI16 + 2.D0*CHI17) -
     &        A6_A5*(CHI8+2.D0*CHI12+CHI13+2.D0*CHI15)
         PSI5 = PSI5 / (A6_A5 + PSI6 + CHI7 + CHI14 + 2.D0*CHI16 + 2.D0*CHI17)
         PSI5 = MIN(MAX(PSI5, TINY),CHI5)
      ELSE
         PSI5 = CHI5
      ENDIF

      IF (CHI4.GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         A4M  = A4M0 * (GAMA(5) / GAMA(10))**2
         BB   = CHI4 + PSI6 + PSI5 + A4M
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 = 0.5d0*(BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION *********************************************
C
      MOLAL (3) = PSI4                                            ! NH4I
      MOLAL (7) = PSI5 + CHI8 + 2.D0*CHI12 + CHI13 + 2.D0*CHI15   ! NO3I

      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
C
C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
C
      SMIN = PSI5 + PSI6 - PSI4
      CALL CALCPH (SMIN, HI, OHI)

      MOLAL (1) = HI
      MOLALR(5) = MIN(PSI4, PSI5)         ! NH4NO3
      FRNH4     = MAX(PSI4 - PSI5, ZERO)  ! "FREE" NH3
      MOLALR(6) = MIN(PSI6, FRNH4)        ! NH4CL

      WATER = WATER0 + SUM( MOLALR(5:6)*M0I(5:6) )
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
      A6      = A60 * (WATER/GAMA(11))**2
      FUNCP13 = MOLAL(1)*MOLAL(4) - GHCL * A6
C
      END FUNCTION FUNCP13

      END SUBROUTINE CALCP13



C======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCL9
C *** CASE L9
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : CASO4
C     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4, NA2SO4, K2SO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCL9

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: A9, BB, CC, DD, LAMDA, SO4A, HSO4A
      REAL(8) :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8
      REAL    :: MOLALR(21)
      INTEGER :: I

C *** FIND DRY COMPOSITION **********************************************
C
      CALL CALCL1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCL1 run
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4

      CLC     = ZERO
      CNAHSO4 = ZERO
      CNA2SO4 = ZERO
      CNH42S4 = ZERO
      CNH4HS4 = ZERO
      CK2SO4  = ZERO
      CMGSO4  = ZERO
      CKHSO4  = ZERO

      MOLAL(2) = 2.D0*CHI4 + CHI3               ! NAI
      MOLAL(3) = 3.D0*CHI2 + 2.D0*CHI5 + CHI1   ! NH4I
      MOLAL(9) = CHI8 + 2.0D0*CHI6              ! KI
      MOLAL(10)= CHI7                           ! MGI

      MOLALR(02) = CHI4         ! NA2SO4
      MOLALR(04) = CHI5         ! (NH4)2SO4
      MOLALR(09) = CHI1         ! NH4HSO4
      MOLALR(12) = CHI3         ! NAHSO4
      MOLALR(13) = CHI2         ! LC
      MOLALR(17) = CHI6         ! K2SO4
      MOLALR(18) = CHI8         ! KHSO4
      MOLALR(21) = CHI7         ! MGSO4

      WATER = MOLALR(02)*M0I(02) + MOLALR(04)*M0I(04) + MOLALR(09)*M0I(09)
     &      + MOLALR(12)*M0I(12) + MOLALR(13)*M0I(13) + MOLALR(17)*M0I(17)
     &      + MOLALR(18)*M0I(18) + MOLALR(21)*M0I(21)

      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      SO4A  = CHI2 + CHI4 + CHI5 + CHI6 + CHI7
      HSO4A = CHI1 + CHI2 + CHI3 + CHI8

      DO I = 1, NSWEEP

         A9 = XK1 * WATER * GAMA(8)**2 / GAMA(7)**3
C
C  CALCULATE DISSOCIATION QUANTITIES
C
         BB   = SO4A + A9
         DD   = BB*BB + 4.D0* A9 * HSO4A
         LAMDA= 0.5D0*(SQRT(DD) - BB)
         LAMDA= MIN( MAX( LAMDA, TINY), HSO4A)
C
C *** CALCULATE SPECIATION ********************************************
C
         MOLAL(1) = LAMDA                 ! HI
         MOLAL(5) = SO4A  + LAMDA         ! SO4I
         MOLAL(6) = HSO4A - LAMDA         !HSO4I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            EXIT
         ENDIF

      ENDDO

      END SUBROUTINE CALCL9


C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCL1A
C *** CASE L1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCL1A

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: FRSO4, FRK, FRMG, FRNH4, CAFR, FRNA
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CCASO4  = MIN (W(6), W(2))                    ! CCASO4
      FRSO4   = MAX(W(2) - CCASO4, ZERO)
      CAFR    = MAX(W(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*W(7), FRSO4)             ! CK2SO4
      FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
      FRSO4   = MAX(FRSO4 - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*W(1), FRSO4)             ! CNA2SO4
      FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
      FRSO4   = MAX(FRSO4 - CNA2SO4, ZERO)
      CMGSO4  = MIN (W(8), FRSO4)                   ! CMGSO4
      FRMG    = MAX(W(8) - CMGSO4, ZERO)
      FRSO4   = MAX(FRSO4 - CMGSO4, ZERO)
C
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      CKHSO4  = ZERO
C
      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
C
      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
         IF (CK2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CKHSO4 = 2.D0*FRSO4
            CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
       ENDIF
      ENDIF
C
C *** CALCULATE GAS SPECIES ********************************************
C
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO

      END SUBROUTINE CALCL1A



C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCK4
C *** CASE K4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CASO4
C
C *** COPYRIGHT 1996-2012, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C *** UPDATE|ADJOINT BY SHANNON CAPPS
C
C=======================================================================
C
      SUBROUTINE CALCK4

      USE ISRPIA
      IMPLICIT NONE

      REAL(8) :: LAMDA, KAPA, A4, BB, CC, DD
      REAL(8) :: CHI1, CHI2, CHI3, CHI4
      REAL    :: MOLALR(7:21)
      INTEGER :: I
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
      CHI2   = W(1)                !  Total NA initially as NaHSO4
      CHI3   = W(7)                !  Total K initially as KHSO4
      CHI4   = W(8)                !  Total Mg initially as MgSO4
C
      CCASO4 = W(6)
      LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4

      MOLAL (2) = CHI2                                            ! NAI
      MOLAL (3) = CHI1                                            ! NH4I
      MOLAL (9) = CHI3                                            ! KI
      MOLAL (10)= CHI4                                            ! MGI

      MOLALR(07) = LAMDA
      MOLALR(09) = CHI1                        ! NH4HSO4
      MOLALR(12) = CHI2                        ! NAHSO4
      MOLALR(18) = CHI3                        ! KHSO4
      MOLALR(21) = CHI4                        ! MGSO4

      WATER = MOLALR(07)*M0I(07) + MOLALR(09)*M0I(09) + MOLALR(12)*M0I(12)
     &      + MOLALR(18)*M0I(18) + MOLALR(21)*M0I(21)
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO I = 1, NSWEEP

         A4 = XK1 * WATER * GAMA(8)**2 / GAMA(7)**3

         BB   = A4 + LAMDA + CHI4
         CC   = A4*(LAMDA + CHI3 + CHI2 + CHI1) - LAMDA*CHI4
         DD   = MAX(BB*BB + 4.D0*CC, ZERO)
         KAPA = 0.5D0*(SQRT(DD) - BB)
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
         MOLAL (1) = MAX(LAMDA + KAPA, TINY)                         ! HI
         MOLAL (5) = MAX(KAPA + CHI4, ZERO)                          ! SO4I
         MOLAL (6) = MAX(LAMDA + CHI1 + CHI2 + CHI3 - KAPA, ZERO)    ! HSO4I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT
         ELSE
            EXIT
         ENDIF

      ENDDO

      END SUBROUTINE CALCK4
