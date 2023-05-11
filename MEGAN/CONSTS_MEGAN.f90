module CONSTS_MEGAN

  implicit none

  INTEGER, PARAMETER :: N_MGN_SPC = 20
  INTEGER, PARAMETER :: Layers = 5

  REAL, PARAMETER ::                      &
       ConvertWm2toUmolm2s = 4.766,       & ! Convert radiation from [W/m2] to [umol/m2/s1]
       SolarConstant       = 1367,        & ! Solar constant [W/m2]
       WaterAirRatio       = 18.016/28.97   ! Ratio between water and air molecules

  REAL, PARAMETER :: pstd_Sun = 200.0, pstd_Shade = 50.0
  REAL, PARAMETER :: Cce =0.56

!=======================================================================
!  CANOPY.EXT
!  This include file contains MEGAN species
!
!  Who                   When       What
!  ---------------------------------------------------------------------
!  Xuemei Wang          06/16/2009 - Creates this file
!=======================================================================

      CHARACTER(16) :: MGN_SPC(N_MGN_SPC)
      REAL          :: CLeo ( N_MGN_SPC)
      REAL          :: Ctm1 ( N_MGN_SPC)

      DATA     MGN_SPC(  1) , CLeo(1) ,   Ctm1(1)  &
               / 'ISOP ' , 2.0, 95.0   /
      DATA     MGN_SPC(  2), CLeo(2) ,   Ctm1(2)   &
               / 'MYRC ', 1.83, 80.0   /
      DATA     MGN_SPC(  3), CLeo(3) ,   Ctm1(3)   &
               / 'SABI ', 1.83, 80.0   /
      DATA     MGN_SPC(  4), CLeo(4) ,   Ctm1(4)   &
              / 'LIMO  ',  1.83, 80.0  /
      DATA     MGN_SPC(  5), CLeo(5) ,   Ctm1(5)   &
              / 'A_3CAR ',  1.83, 80.0    /
      DATA     MGN_SPC(  6), CLeo(6) ,   Ctm1(6)   &
              / 'OCIM    ', 1.83, 80.0    /
      DATA     MGN_SPC(  7), CLeo(7) ,   Ctm1(7)   &
              / 'BPIN    ',  1.83, 80.0    /
      DATA     MGN_SPC(  8), CLeo(8) ,   Ctm1(8)   &
              / 'APIN   ',  1.83, 80.0    /
      DATA     MGN_SPC(  9), CLeo(9) ,   Ctm1(9)   &
              / 'OMTP   ',   1.83, 80.0   /
      DATA     MGN_SPC( 10), CLeo(10) ,   Ctm1(10) &
              / 'FARN   ',  2.37,130.0    /
      DATA     MGN_SPC( 11), CLeo(11),   Ctm1(11)  &
              / 'BCAR   ',  2.37,130.0    /
      DATA     MGN_SPC( 12), CLeo(12) ,   Ctm1(12) &
              / 'OSQT   ',  2.37,130.0    /
      DATA     MGN_SPC( 13), CLeo(13) ,   Ctm1(13) &
              / 'MBO    ',  2.0, 95.0    /
      DATA     MGN_SPC( 14), CLeo(14) ,   Ctm1(14) &
              / 'MEOH   ',  1.6, 60.0    /
      DATA     MGN_SPC( 15), CLeo(15) ,   Ctm1(15) &
              / 'ACTO   ',  1.83, 80.0    /
      DATA     MGN_SPC( 16), CLeo(16) ,   Ctm1(16) &
              / 'CO     ',  1.6, 60.0    /
      DATA     MGN_SPC( 17), CLeo(17) ,   Ctm1(17) &
              / 'NO     ',  1.83, 80.0    /
      DATA     MGN_SPC( 18), CLeo(18) ,   Ctm1(18) &
              / 'BIDER  ',  2.0, 95.0    /
      DATA     MGN_SPC( 19), CLeo(19) ,   Ctm1(19) &
              / 'STRESS ', 1.83, 80.0    /
      DATA     MGN_SPC( 20), CLeo(20) ,   Ctm1(20) &
              / 'OTHER  ', 1.83, 80.0     /


        INTEGER, PARAMETER :: NMON = 12
        INTEGER, PARAMETER :: NrTyp = 16, NrCha = 16

    ! 1  = canopy depth
    ! 2  = leaf width
    ! 3  = leaf length
    ! 4  = canopy height
    ! 5  = scattering coefficient for PPFD
    ! 6  = scattering coefficient for near IR
    ! 7  = reflection coefficient for diffuse PPFD
    ! 8  = reflection coefficient for diffuse near IR
    ! 9  = clustering coefficient (accounts for leaf clumping influence on mean
    !    projected leaf area in the direction of the suns beam)
    !    use 0.85 for default, corn=0.4-0.9; Pine=0.6-1.0; oak=0.53-0.67;
    !    tropical rainforest=1.1
    ! 10 = leaf IR emissivity
    ! 11 = leaf stomata and cuticle factor: 1=hypostomatous, 2=amphistomatous,
    !     1.25=hypostomatous but with some transpiration through cuticle
    ! 12 = daytime temperature lapse rate (K m-1)
    ! 13 = nighttime temperature lapse rate (K m-1)
    ! 14 = warm (>283K) canopy total humidity change (Pa)
    ! 15 = cool (>= 283K) canopy total humidity change (Pa)
    ! 16 = normalized canopy depth where wind is negligible

    !     NT NT NT TF BT TF BT BT SB SB SB HB HB HB CR CR
        REAL,DIMENSION(NrCha,NrTyp) :: Canopychar = RESHAPE( &
          (/ 16.   , 16.   , 16.  ,    16.,   &
             16.   , 16.   , 16.  ,    16. ,  &
              1.   , 1.    ,  1.  ,    0.756, &
              0.756,  0.756 ,  1.  ,    1. ,  &
              0.05,  0.05 ,   0.05,    0.05,  &
              0.05,  0.05 ,   0.05,    0.05,  &
              0.015,   0.015, 0.015,   0.015, &
              0.015,   0.015, 0.02 ,   0.02 , &
              0.1  ,  0.1  ,  0.1  ,   0.1  , &
              0.1  ,  0.1  ,  0.1  ,   0.1  , &
              0.1  ,  0.1  ,  0.1  ,   0.15 , &
              0.15 ,  0.15 ,  0.15 ,   0.15 , &
              24.  , 24.   ,  24.,     24.  , &
              24.  , 24.   ,  24.,     24.  , &
              2.   , 2.    ,  2. ,     0.75 , &
              0.75 , 0.75  ,  1. ,     1.   , &
              0.2  ,  0.2  ,  0.2  ,   0.2  , &
              0.2  ,  0.2  ,  0.2  ,   0.2  , &
              0.2  ,  0.2  ,  0.2  ,   0.2  , &
              0.2  ,  0.2  ,  0.2  ,   0.2  , &
              0.8  ,  0.8  ,  0.8  ,   0.8  , &
              0.8  ,  0.8  ,  0.8  ,   0.8  , &
              0.8  ,  0.8  ,  0.8  ,   0.8  , &
              0.8  ,  0.8  ,  0.8  ,   0.8  , &
              0.057,  0.057,  0.057,   0.057, &
              0.057,  0.057,  0.057,   0.057, &
              0.057,  0.057,  0.057,   0.057, &
              0.057,  0.057,  0.057,   0.057, &
              0.389,  0.389,  0.389,   0.389, &
              0.389,  0.389,  0.389,   0.389, &
              0.389,  0.389,  0.389,   0.389, &
              0.389,  0.389,  0.389,   0.389, &
              0.85 ,  0.85 ,  0.85 ,   1.1  , &
              0.95 ,  1.1  ,  0.95 ,   0.95 , &
              0.85 ,  0.85 ,  0.85 ,   0.76 , &
              0.76 ,  0.76 ,  0.65 ,   0.65 , &
              0.95 ,  0.95 ,  0.95 ,   0.95 , &
              0.95 ,  0.95 ,  0.95 ,   0.95 , &
              0.95 ,  0.95 ,  0.95 ,   0.95 , &
              0.95 ,  0.95 ,  0.95 ,   0.95 , &
              1.25 ,  1.25 ,  1.25 ,   1.25 , &
              1.25 ,  1.25 ,  1.25 ,   1.25 , &
              1.00 ,  1.00 ,  1.00 ,   1.25 , &
              1.25 ,  1.25 ,  1.25 ,   1.25 , &
              0.06 ,  0.06 ,  0.06 ,   0.06 , &
              0.06 ,  0.06 ,  0.06 ,   0.06 , &
              0.06 ,  0.06 ,  0.06 ,   0.06 , &
              0.06 ,  0.06 ,  0.06 ,   0.06 , &
             -0.06 , -0.06 , -0.06 ,  -0.06 , &
             -0.06 , -0.06 , -0.06 ,  -0.06 , &
             -0.06 , -0.06 , -0.06 ,  -0.06 , &
             -0.06 , -0.06 , -0.06 ,  -0.06 , &
              700. ,  700. ,  700. ,   700. , &
              700. ,  700. ,  700. ,   700. , &
              700. ,  700. ,  700. ,   700. , &
              700. ,  700. ,  700. ,   700. , &
              150. ,  150. ,  150. ,   150. , &
              150. ,  150. ,  150. ,   150. , &
              150. ,  150. ,  150. ,   150. , &
              150. ,  150. ,  150. ,   150. , &
              0.7  ,  0.7  ,  0.7  ,   0.7  , &
              0.7  ,  0.7  ,  0.7  ,   0.7  , &
              0.7  ,  0.7  ,  0.7  ,   0.7  , &
              0.7  ,  0.7  ,  0.7  ,   0.7    /) &
            ,SHAPE=(/NrCha,NrTyp/)               &
            ,ORDER=(/2,1/)                      )

!=======================================================================
!  LD_FCT.EXT
!  This include file contains "light dependent" factors.
!
!
!
!  MEGAN v2.02
!  INPUT version 210
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Guenther A.  08/11/07 - Move from MEGAN v2.0 to MEGAN v2.02 with updates.
!                          See the update document.
!  Jiang X.     05/07/12 - Update LDFs
!=======================================================================

      INTEGER, PARAMETER :: N_LDF_SPC = 20
      CHARACTER(16)      :: LDF_SPC(N_LDF_SPC)
      REAL               :: LDF_FCT(N_LDF_SPC)
      INTEGER            :: LDF_MAP(N_LDF_SPC)

      DATA     LDF_SPC(  1)      , LDF_FCT(  1), LDF_MAP(  1) &
             / 'ISOP            ', 0.999      , 1            /

      DATA     LDF_SPC(  2)      , LDF_FCT(  2), LDF_MAP(  2) &
             / 'MYRC            ', 0.6        , 2            /

      DATA     LDF_SPC(  3)      , LDF_FCT(  3), LDF_MAP(  3) &
             / 'SABI            ', 0.6         , 3            /

      DATA     LDF_SPC(  4)      , LDF_FCT(  4), LDF_MAP(  4) &
             / 'LIMO            ', 0.4        , 4            /

      DATA     LDF_SPC(  5)      , LDF_FCT(  5), LDF_MAP(  5) &
             / 'A_3CAR          ', 0.4        , 5            /

      DATA     LDF_SPC(  6)      , LDF_FCT(  6), LDF_MAP(  6) &
             / 'OCIM            ', 0.4         , 6            /

      DATA     LDF_SPC(  7)      , LDF_FCT(  7), LDF_MAP(  7) &
             / 'BPIN            ', 0.4         , 7            /

      DATA     LDF_SPC(  8)      , LDF_FCT(  8), LDF_MAP(  8) &
             / 'APIN            ', 0.6         , 8            /

      DATA     LDF_SPC(  9)      , LDF_FCT(  9), LDF_MAP(  9) &
             / 'OMTP            ', 0.4         , 9            /

      DATA     LDF_SPC( 10)      , LDF_FCT( 10), LDF_MAP( 10) &
             / 'FARN            ', 0.5         , 10           /

      DATA     LDF_SPC( 11)      , LDF_FCT( 11), LDF_MAP( 11) &
             / 'BCAR            ', 0.5         , 11           /

      DATA     LDF_SPC( 12)      , LDF_FCT( 12), LDF_MAP( 12) &
             / 'OSQT            ', 0.5         , 12           /

      DATA     LDF_SPC( 13)      , LDF_FCT( 13), LDF_MAP( 13) &
             / 'MBO             ', 0.999      , 13           /

      DATA     LDF_SPC( 14)      , LDF_FCT( 14), LDF_MAP( 14) &
             / 'MEOH            ', 0.8        , 14           /

      DATA     LDF_SPC( 15)      , LDF_FCT( 15), LDF_MAP( 15) &
             / 'ACTO            ', 0.2        , 15           /

      DATA     LDF_SPC( 16)      , LDF_FCT( 16), LDF_MAP( 16) &
             / 'CO              ', 0.999      , 16           /

      DATA     LDF_SPC( 17)      , LDF_FCT( 17), LDF_MAP( 17) &
             / 'NO              ', 0.0         , 17           /

      DATA     LDF_SPC( 18)      , LDF_FCT( 18), LDF_MAP( 18) &
             / 'BIDER           ', 0.8         , 18           /

      DATA     LDF_SPC( 19)      , LDF_FCT( 19), LDF_MAP( 19) &
             / 'STRESS          ', 0.8         , 19           /

      DATA     LDF_SPC( 20)      , LDF_FCT( 20), LDF_MAP( 20) &
             / 'OTHER           ', 0.2         , 20           /

!=======================================================================
!  TEMPD_PRM.EXT
!  This include file contains "temperature dependent" parameter for
!  light-independent emissions.
!
!
!  MEGAN v2.02
!  INPUT version 210
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Guenther A.  08/11/07 - Move from MEGAN v2.0 to MEGAN v2.02 with updates.
!                          See the update document.
!=======================================================================

      INTEGER, PARAMETER :: N_TDF_SPC = 20
      CHARACTER(16)      :: TDF_SPC(N_TDF_SPC)
      REAL               :: TDF_PRM(N_TDF_SPC)
      INTEGER            :: TDF_MAP(N_TDF_SPC)

      DATA     TDF_SPC(  1)      , TDF_PRM(  1), TDF_MAP(  1) &
             / 'ISOP            ', 0.13        , 1           /

      DATA     TDF_SPC(  2)      , TDF_PRM(  2), TDF_MAP(  2) &
             / 'MYRC            ', 0.1         , 2           /

      DATA     TDF_SPC(  3)      , TDF_PRM(  3), TDF_MAP(  3) &
             / 'SABI            ', 0.1         , 3           /

      DATA     TDF_SPC(  4)      , TDF_PRM(  4), TDF_MAP(  4) &
             / 'LIMO            ', 0.1         , 4           /

      DATA     TDF_SPC(  5)      , TDF_PRM(  5), TDF_MAP(  5) &
             / 'A_3CAR          ', 0.1         , 5           /

      DATA     TDF_SPC(  6)      , TDF_PRM(  6), TDF_MAP(  6) &
             / 'OCIM            ', 0.1         , 6           /

      DATA     TDF_SPC(  7)      , TDF_PRM(  7), TDF_MAP(  7) &
             / 'BPIN            ', 0.1         , 7           /
      DATA     TDF_SPC(  8)      , TDF_PRM(  8), TDF_MAP(  8) &
             / 'APIN            ', 0.1         , 8           /

      DATA     TDF_SPC(  9)      , TDF_PRM(  9), TDF_MAP(  9) &
             / 'OMTP            ', 0.1         , 9           /

      DATA     TDF_SPC( 10)      , TDF_PRM( 10), TDF_MAP( 10) &
             / 'FARN            ', 0.17        , 10          /

      DATA     TDF_SPC( 11)      , TDF_PRM( 11), TDF_MAP( 11) &
             / 'BCAR            ', 0.17        , 11          /

      DATA     TDF_SPC( 12)      , TDF_PRM( 12), TDF_MAP( 12) &
             / 'OSQT            ', 0.17        , 12          /

      DATA     TDF_SPC( 13)      , TDF_PRM( 13), TDF_MAP( 13) &
             / 'MBO             ', 0.13        , 13          /

      DATA     TDF_SPC( 14)      , TDF_PRM( 14), TDF_MAP( 14) &
             / 'MEOH            ', 0.08        , 14          /
      DATA     TDF_SPC( 15)      , TDF_PRM( 15), TDF_MAP( 15) &
             / 'ACTO            ', 0.10        , 15          /

      DATA     TDF_SPC( 16)      , TDF_PRM( 16), TDF_MAP( 16) &
             / 'CO              ', 0.08        , 16          /

      DATA     TDF_SPC( 17)      , TDF_PRM( 17), TDF_MAP( 17) &
             / 'NO              ', 0.10        , 17          /

      DATA     TDF_SPC( 18)      , TDF_PRM( 18), TDF_MAP( 18) &
             / 'BIDER           ', 0.13        , 18          /

      DATA     TDF_SPC( 19)      , TDF_PRM( 19), TDF_MAP( 19) &
             / 'STRESS          ', 0.1        , 19          /

      DATA     TDF_SPC( 20)      , TDF_PRM( 20), TDF_MAP( 20) &
             / 'OTHER           ', 0.1        , 20          /
!------------------------------------------------------------------------
!  relative emission activity index
!------------------------------------------------------------------------
      INTEGER, PARAMETER :: N_REA_SPC = 20
      CHARACTER(16)      :: REA_SPC  (N_REA_SPC)
      INTEGER            :: REA_INDEX(N_REA_SPC)

      DATA     REA_SPC(  1)      , REA_INDEX(  1) &
             / 'ISOP            ', 5                  /

      DATA     REA_SPC(  2)      , REA_INDEX(  2) &
             / 'MYRC            ', 2                  /

      DATA     REA_SPC(  3)      , REA_INDEX(  3) &
             / 'SABI            ', 2                  /

      DATA     REA_SPC(  4)      , REA_INDEX(  4) &
             / 'LIMO            ', 2                  /

      DATA     REA_SPC(  5)      , REA_INDEX(  5) &
             / 'A_3CAR          ', 2                  /

      DATA     REA_SPC(  6)      , REA_INDEX(  6) &
             / 'OCIM            ', 2                   /

      DATA     REA_SPC(  7)      , REA_INDEX(  7) &
             / 'BPIN            ', 2                   /

      DATA     REA_SPC(  8)      , REA_INDEX(  8)  &
             / 'APIN            ', 2                   /

      DATA     REA_SPC(  9)      , REA_INDEX(  9) &
             / 'OMTP            ', 2                   /

      DATA     REA_SPC( 10)      , REA_INDEX( 10) &
             / 'FARN            ', 3                   /

      DATA     REA_SPC( 11)      , REA_INDEX( 11) &
             / 'BCAR            ', 3                   /

      DATA     REA_SPC( 12)      , REA_INDEX( 12)         &
             / 'OSQT            ', 3                    /

      DATA     REA_SPC( 13)      , REA_INDEX( 13) &
             / 'MBO             ', 5                    /

      DATA     REA_SPC( 14)      , REA_INDEX( 14) &
             / 'MEOH            ', 4                     /

      DATA     REA_SPC( 15)      , REA_INDEX( 15) &
             / 'ACTO            ', 1                     /

      DATA     REA_SPC( 16)      , REA_INDEX( 16) &
             / 'CO             ', 1                     /

      DATA     REA_SPC( 17)      , REA_INDEX( 17) &
             / 'NO              ', 1                     /

      DATA     REA_SPC( 18)      , REA_INDEX( 18) &
             / 'BIDER           ', 1                     /

      DATA     REA_SPC( 19)      , REA_INDEX( 19) &
             / 'STRESS          ', 1                     /

      DATA     REA_SPC( 20)      , REA_INDEX( 20) &
             / 'OTHER           ', 1                     /

!=======================================================================
!  REL_EM_ACT.EXT
!  This include file contains "produciton and loss within canopy"
!  factors.
!
!
!  MEGAN v2.02
!  INPUT version 200
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
!=======================================================================

      INTEGER, PARAMETER :: N_CAT = 5
      REAL               :: Anew(N_CAT)
      REAL               :: Agro(N_CAT)
      REAL               :: Amat(N_CAT)
      REAL               :: Aold(N_CAT)

      DATA    Anew(  1),  Agro(  1),  Amat(  1),  Aold(  1) &
           /  1.0      ,  1.0      ,  1.0      ,  1.0       /

      DATA    Anew(  2),  Agro(  2),  Amat(  2),  Aold(  2) &
           /  2.0      ,  1.8      ,  1.0     ,  1.05       /

      DATA    Anew(  3),  Agro(  3),  Amat(  3),  Aold(  3) &
           /  0.4      ,  0.6      ,  1.0    ,  0.95       /

      DATA    Anew(  4),  Agro(  4),  Amat(  4),  Aold(  4) &
           /  3.5      ,  3.0      ,  1.0     ,  1.2       /

      DATA    Anew(  5),  Agro(  5),  Amat(  5),  Aold(  5) &
           /  0.05     ,  0.6      ,  1.0    ,  0.9       /

!=======================================================================
!  PDT_LOS_CP.EXT
!  This include file contains "produciton and loss within canopy"
!  factors.
!
!
!  MEGAN v2.02
!  INPUT version 200
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
!=======================================================================

      INTEGER, PARAMETER :: N_RHO_SPC = 20
      CHARACTER(16)      :: RHO_SPC(N_RHO_SPC)
      REAL               :: RHO_FCT(N_RHO_SPC)
      INTEGER            :: RHO_MAP(N_RHO_SPC)

      DATA     RHO_SPC(  1)      , RHO_FCT(  1), RHO_MAP(  1) &
             / 'ISOP            ',   1.0       , 1           /

      DATA     RHO_SPC(  2)      , RHO_FCT(  2), RHO_MAP(  2) &
             / 'MYRC            ',   1.0       , 2           /

      DATA     RHO_SPC(  3)      , RHO_FCT(  3), RHO_MAP(  3) &
             / 'SABI            ',   1.0       , 3           /
      DATA     RHO_SPC(  4)      , RHO_FCT(  4), RHO_MAP(  4) &
             / 'LIMO            ',   1.0       , 4           /

      DATA     RHO_SPC(  5)      , RHO_FCT(  5), RHO_MAP(  5) &
             / 'A_3CAR          ',   1.0       , 5           /

      DATA     RHO_SPC(  6)      , RHO_FCT(  6), RHO_MAP(  6) &
             / 'OCIM            ',   1.0       , 6           /

      DATA     RHO_SPC(  7)      , RHO_FCT(  7), RHO_MAP(  7) &
             / 'BPIN            ',   1.0       , 7           /

      DATA     RHO_SPC(  8)      , RHO_FCT(  8), RHO_MAP(  8) &
             / 'APIN            ',   1.0       , 8           /

      DATA     RHO_SPC(  9)      , RHO_FCT(  9), RHO_MAP(  9) &
             / 'OMTP            ',   1.0       , 9           /

      DATA     RHO_SPC( 10)      , RHO_FCT( 10), RHO_MAP( 10) &
             / 'FARN            ',   1.0       , 10          /

      DATA     RHO_SPC( 11)      , RHO_FCT( 11), RHO_MAP( 11) &
             / 'BCAR            ',   1.0       , 11          /

      DATA     RHO_SPC( 12)      , RHO_FCT( 12), RHO_MAP( 12) &
             / 'OSQT            ',   1.0       , 12          /

      DATA     RHO_SPC( 13)      , RHO_FCT( 13), RHO_MAP( 13) &
             / 'MBO             ',   1.0       , 13          /
      DATA     RHO_SPC( 14)      , RHO_FCT( 14), RHO_MAP( 14) &
             / 'MEOH            ',   1.0       , 14          /

      DATA     RHO_SPC( 15)      , RHO_FCT( 15), RHO_MAP( 15) &
             / 'ACTO            ',   1.0       , 15          /

      DATA     RHO_SPC( 16)      , RHO_FCT( 16), RHO_MAP( 16) &
             / 'CO             ',   1.0       , 16          /

      DATA     RHO_SPC( 17)      , RHO_FCT( 17), RHO_MAP( 17) &
             / 'NO              ',   1.0       , 17          /

      DATA     RHO_SPC( 18)      , RHO_FCT( 18), RHO_MAP( 18) &
             / 'BIDER           ',   1.0       , 18          /

      DATA     RHO_SPC( 19)      , RHO_FCT( 19), RHO_MAP( 19) &
             / 'STRESS          ',   1.0       , 19          /

      DATA     RHO_SPC( 20)      , RHO_FCT( 20), RHO_MAP( 20) &
             / 'OTHER           ',   1.0       , 20          /

!  TEMPD_PRM
!  This include file contains MEGAN species
!
!
!
!  MEGAN v2.1
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Xuemei      06/15/2009 - Creates this file
!=======================================================================

      INTEGER, PARAMETER :: N_PRM_SPC = 20
      REAL               :: CT1 (N_PRM_SPC)
      REAL               :: Cceo(N_PRM_SPC)

      DATA    CT1(1),Cceo (1)    /  95.0, 2.0  /
      DATA    CT1(2),Cceo (2)    /  80.0, 1.83  /
      DATA    CT1(3),Cceo (3)    /  80.0, 1.83  /
      DATA    CT1(4),Cceo (4)    /  80.0, 1.83  /
      DATA    CT1(5),Cceo (5)    /  80.0, 1.83  /
      DATA    CT1(6),Cceo (6)    /  80.0, 1.83  /
      DATA    CT1(7),Cceo (7)    /  80.0, 1.83  /
      DATA    CT1(8),Cceo (8)    /  80.0, 1.83  /
      DATA    CT1(9),Cceo (9)    /  80.0, 1.83  /
      DATA    CT1(10),Cceo (10)    / 130.0, 2.37  /
      DATA    CT1(11),Cceo (11)    / 130.0, 2.37  /
      DATA    CT1(12),Cceo (12)    / 130.0, 2.37  /
      DATA    CT1(13),Cceo (13)    /  95.0, 2.0  /
      DATA    CT1(14),Cceo (14)    /  60.0, 1.6  /
      DATA    CT1(15),Cceo (15)    /  80.0, 1.83  /
      DATA    CT1(16),Cceo (16)    /  60.0, 1.6  /
      DATA    CT1(17),Cceo (17)    /  80.0, 1.83  /
      DATA    CT1(18),Cceo (18)    /  95.0, 2.0  /
      DATA    CT1(19),Cceo (19)    /  80.0, 1.83  /
      DATA    CT1(20),Cceo (20)    /  80.0, 1.83  /

!=======================================================================
!  PFT_MGN.EXT
!  This include file contains MEGAN species
!
!
!
!  MEGAN v2.10
!  INPUT version XXX
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          07/07/11 - Creates this file for MEGANv2.10
!=======================================================================

      INTEGER, PARAMETER :: N_MGN_PFT = 36
      CHARACTER(16)      :: MGN_PFT(N_MGN_PFT)
      CHARACTER(36)      :: MGN_NAM(N_MGN_PFT)
      REAL               :: MGN_DUM(N_MGN_PFT)

      DATA     MGN_PFT(  1), MGN_DUM(  1) / 'NT_EG_TEMP      ', 1.0    /
      DATA     MGN_PFT(  2), MGN_DUM(  2) / 'NT_DC_BORL      ', 1.0    /
      DATA     MGN_PFT(  3), MGN_DUM(  3) / 'NT_EG_BORL      ', 1.0    /
      DATA     MGN_PFT(  4), MGN_DUM(  4) / 'BT_EG_TROP      ', 1.0    /
      DATA     MGN_PFT(  5), MGN_DUM(  5) / 'BT_EG_TEMP      ', 1.0    /
      DATA     MGN_PFT(  6), MGN_DUM(  6) / 'BT_DC_TROP      ', 1.0    /
      DATA     MGN_PFT(  7), MGN_DUM(  7) / 'BT_DC_TEMP      ', 1.0    /
      DATA     MGN_PFT(  8), MGN_DUM(  8) / 'BT_DC_BORL      ', 1.0    /
      DATA     MGN_PFT(  9), MGN_DUM(  9) / 'SB_EG_TEMP      ', 1.0    /
      DATA     MGN_PFT( 10), MGN_DUM( 10) / 'SB_DC_TEMP      ', 1.0    /
      DATA     MGN_PFT( 11), MGN_DUM( 11) / 'SB_DC_BORL      ', 1.0    /
      DATA     MGN_PFT( 12), MGN_DUM( 12) / 'GS_C3_COLD      ', 1.0    /
      DATA     MGN_PFT( 13), MGN_DUM( 13) / 'GS_C3_COOL      ', 1.0    /
      DATA     MGN_PFT( 14), MGN_DUM( 14) / 'GS_C3_WARM      ', 1.0    /
      DATA     MGN_PFT( 15), MGN_DUM( 15) / 'CORN            ', 1.0    /
      DATA     MGN_PFT( 16), MGN_DUM( 16) / 'CROP            ', 1.0    /

      DATA     MGN_NAM(  1) /'Needleaf evergreen temperate tree   '/
      DATA     MGN_NAM(  2) /'Needleaf deciduous boreal tree      '/
      DATA     MGN_NAM(  3) /'Needleaf evergreen boreal tree      '/
      DATA     MGN_NAM(  4) /'Broadleaf evergreen tropical tree   '/
      DATA     MGN_NAM(  5) /'Broadleaf evergreen temperate tree  '/
      DATA     MGN_NAM(  6) /'Broadleaf deciduous tropical tree   '/
      DATA     MGN_NAM(  7) /'Broadleaf deciduous temperate tree  '/
      DATA     MGN_NAM(  8) /'Broadleaf deciduous boreal tree     '/
      DATA     MGN_NAM(  9) /'Broadleaf evergreen temperate shrub '/
      DATA     MGN_NAM( 10) /'Broadleaf deciduous temperate shrub '/
      DATA     MGN_NAM( 11) /'Broadleaf deciduous boreal shrub    '/
      DATA     MGN_NAM( 12) /'Cold C3 grass                       '/
      DATA     MGN_NAM( 13) /'Cool C3 grass                       '/
      DATA     MGN_NAM( 14) /'Warm C3 grass                       '/
      DATA     MGN_NAM( 15) /'Corn                                '/
      DATA     MGN_NAM( 16) /'Other crops                         '/


!=======================================================================
!  EF_MGN20.EXT
!  This include file contains EF for 20 MEGAN species.  The values in
!  this file must be in the same order as in SPC_MGN.EXT
!
!  MEGAN v2.1.0
!  INPUT version 210
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/2006 - Creates this file
!  Guenther A.  08/11/2007 - Creates this file again with updates and move
!                          from v2.0 to v2.02
!  Xuemei Wang and Alex 26/07/2011-Extend EFs to 16 PFTs
!  Jiang X.      05/07/2012 - Updates EFs for 20 SPC with new values from Guenther
!=======================================================================

      INTEGER, PARAMETER :: N_EF_SPC = 20     ! Number of categories

!     EF_NT_EG_TEMP EF_NT_DC_BORL EF_NT_EG_BORL EF_BT_EG_TROP
!     EF_BT_EG_TEMP EF_BT_DC_TROP EF_BT_DC_TEMP EF_BT_DC_BORL
!     EF_SB_EG_TEMP EF_SB_DC_TEMP EF_SB_DC_BORL EF_GS_C3_COLD
!     EF_GS_C3_COOL EF_GS_C3_WARM EF_CROP EF_CORN

      REAL,DIMENSION(16,N_EF_SPC) :: EF_ALL = RESHAPE( (/  &
          600.,       1.,    3000.,    7000.,  &
        10000.,    7000.,   10000.,   11000.,  &
         2000.,    4000.,    4000.,    1600.,  &
          800.,     200.,      50.,       1.,  &

           70.,      60.,      70.,      80.,  &
           30.,      80.,      30.,      30.,  &
           30.,      50.,      30.,     0.3,  &
          0.3,     0.3,     0.3,     0.3,  &

           70.,      40.,      70.,      80.,  &
           50.,      80.,      50.,      50.,  &
           50.,      70.,      50.,     0.7,  &
          0.7,     0.7,     0.7,     0.7,  &

          100.,     130.,     100.,      80.,  &
           80.,      80.,      80.,      80.,  &
           60.,     100.,      60.,     0.7,  &
          0.7,     0.7,     0.7,     0.7,  &

          160.,      80.,     160.,      40.,  &
           30.,      40.,      30.,      30.,  &
           30.,     100.,      30.,     0.3,  &
          0.3,     0.3,     0.3,     0.3,  &

           70.,      60.,      70.,     150.,  &
          120.,     150.,     120.,     120.,  &
           90.,     150.,      90.,       2.,  &
            2.,       2.,       2.,       2.,  &

          300.,     200.,     300.,     120.,  &
          130.,     120.,     130.,     130.,  &
          100.,     150.,     100.,     1.5,  &
          1.5,     1.5,     1.5,     1.5,  &

          500.,     510.,     500.,     600.,  &
          400.,     600.,     400.,     400.,  &
          200.,     300.,     200.,       2.,  &
            2.,       2.,       2.,       2.,  &

          180.,     170.,     180.,     150.,  &
          150.,     150.,     150.,     150.,  &
          110.,     200.,     110.,       5.,  &
            5.,       5.,       5.,       5.,  &

           40.,      40.,      40.,      60.,  &
           40.,      60.,      40.,      40.,  &
           40.,      40.,      40.,       3.,  &
            3.,       3.,       4.,       4.,  &

           80.,      80.,      80.,      60.,  &
           40.,      60.,      40.,      40.,  &
           50.,      50.,      50.,       1.,  &
            1.,       1.,       2.,       4.,  &

          120.,     120.,     120.,     120.,  &
          100.,     120.,     100.,     100.,  &
          100.,     100.,     100.,       2.,  &
            2.,       2.,       2.,       2.,  &

          200.,    0.01,      10.,    0.01,  &
         0.01,    0.01,    0.01,    0.01,  &
            2.,    0.01,    0.01,    0.01,  &
         0.01,    0.01,    0.01,    0.01,  &

          900.,     900.,     900.,     500.,  &
          900.,     500.,     900.,     900.,  &
          900.,     900.,     900.,     500.,  &
          500.,     500.,     900.,     900.,  &

          240.,     240.,     240.,     240.,  &
          240.,     240.,     240.,     240.,  &
          240.,     240.,     240.,      80.,  &
           80.,      80.,      80.,      80.,  &

          600.,     600.,     600.,     600.,  &
          600.,     600.,     600.,     600.,  &
          600.,     600.,     600.,     600.,  &
          600.,     600.,     600.,     600.,  &

            2.,       2.,       2.,       2.,  &
            2.,       2.,       2.,       2.,  &
            2.,       2.,       2.,      27.,  &
           27.,      27.,      40.,      68.,  &

          500.,     500.,     500.,     500.,  &
          500.,     500.,     500.,     500.,  &
          500.,     500.,     500.,      80.,  &
           80.,      80.,      80.,      80.,  &

          300.,     300.,     300.,     300.,  &
          300.,     300.,     300.,     300.,  &
          300.,     300.,     300.,     300.,  &
          300.,     300.,     300.,     300.,  &

          140.,     140.,     140.,     140.,  &
          140.,     140.,     140.,     140.,  &
          140.,     140.,     140.,     140.,  &
          140.,     140.,     140.,     140.  /),  &
             SHAPE=(/16,N_EF_SPC/),  &
             ORDER=(/1,2/)                                )

!=======================================================================
!  EFFS_MGN20T150.EXT
!  This include file contains EF fractions for speciation from 20 MEGAN
!  categories to 150 species.  The values in this file must be in the
!  same order as in MAP_MGN20T150.EXT
!
!  MEGAN v2.1.0
!  INPUT version 210
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Guenther A.  08/11/07 - Move from MEGAN v2.0 to MEGAN v2.02 with update on
!                          Nitrogen gas.
!  Guenther A.  26/07/2011-Extend EFs for 16 PFTs
!  Jiang X.     05/07/12 - Update EF factions with new values from Guenther
!=======================================================================

      INTEGER, PARAMETER :: N_EFFS_SPC = 150  ! Number of chemical species

!     EFFS_NT_EG_TEMP EFFS_NT_DC_BORL EFFS_NT_EG_BORL EFFS_BT_EG_TROP EFFS_BT_EG_TEMP
!     EFFS_BT_DC_TROP EFFS_BT_DC_TEMP EFFS_BT_DC_BORL EFSF_SB_EG_TEMP EFFS_SB_DC_TEMP
!     EFFS_SB_DC_BORL EFFS_GS_C3_COLD EFFS_GS_C3_COOL EFFS_GS_C3_WARM EFFS_CROP EFFS_CORN


      REAL,DIMENSION(16,N_EFFS_SPC) :: EFFS_ALL = RESHAPE( (/  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       0.00600 , 0.00600 , 0.00600 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &
       0.00900 , 0.00900 , 0.00900 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.05500 , 0.05500 , 0.05500 , 0.05700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 ,  &
       0.04600 , 0.04600 , 0.04600 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.01700 , 0.01700 , 0.01700 , 0.03400 , 0.03400 , 0.03400 , 0.03400 , 0.03400 ,  &
       0.02800 , 0.02800 , 0.02800 , 0.03100 , 0.03100 , 0.03100 , 0.03100 , 0.03100 ,  &

       0.05500 , 0.05500 , 0.05500 , 0.04600 , 0.04600 , 0.04600 , 0.04600 , 0.04600 ,  &
       0.04600 , 0.04600 , 0.04600 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.03300 , 0.03300 , 0.03300 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &
       0.03700 , 0.03700 , 0.03700 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.05500 , 0.05500 , 0.05500 , 0.05700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 ,  &
       0.04600 , 0.04600 , 0.04600 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.05500 , 0.05500 , 0.05500 , 0.05700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 ,  &
       0.04600 , 0.04600 , 0.04600 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.06700 , 0.06700 , 0.06700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 ,  &
       0.05500 , 0.05500 , 0.05500 , 0.06200 , 0.06200 , 0.06200 , 0.06200 , 0.06200 ,  &

       0.16000 , 0.16000 , 0.16000 , 0.05700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 ,  &
       0.09200 , 0.09200 , 0.09200 , 0.10400 , 0.10400 , 0.10400 , 0.10400 , 0.10400 ,  &

       0.24800 , 0.24800 , 0.24800 , 0.18000 , 0.18000 , 0.18000 , 0.18000 , 0.18000 ,  &
       0.18600 , 0.18600 , 0.18600 , 0.14600 , 0.14600 , 0.14600 , 0.14600 , 0.14600 ,  &

       0.00500 , 0.00500 , 0.00500 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.00800 , 0.00800 , 0.00800 , 0.00800 , 0.00800 , 0.00800 , 0.00800 , 0.00800 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &
       0.00300 , 0.00300 , 0.00300 , 0.00400 , 0.00400 , 0.00400 , 0.00400 , 0.00400 ,  &

       0.00600 , 0.00600 , 0.00600 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &
       0.00900 , 0.00900 , 0.00900 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.02200 , 0.02200 , 0.02200 , 0.04600 , 0.04600 , 0.04600 , 0.04600 , 0.04600 ,  &
       0.03700 , 0.03700 , 0.03700 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.00600 , 0.00600 , 0.00600 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &
       0.00900 , 0.00900 , 0.00900 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &
       0.00300 , 0.00300 , 0.00300 , 0.00400 , 0.00400 , 0.00400 , 0.00400 , 0.00400 ,  &

       0.03300 , 0.03300 , 0.03300 , 0.03400 , 0.03400 , 0.03400 , 0.03400 , 0.03400 ,  &
       0.04600 , 0.04600 , 0.04600 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.00600 , 0.00600 , 0.00600 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &
       0.00900 , 0.00900 , 0.00900 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &
       0.00300 , 0.00300 , 0.00300 , 0.00400 , 0.00400 , 0.00400 , 0.00400 , 0.00400 ,  &

       0.02800 , 0.02800 , 0.02800 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &
       0.04600 , 0.04600 , 0.04600 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.00600 , 0.00600 , 0.00600 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00900 , 0.00900 , 0.00900 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.01100 , 0.01100 , 0.01100 , 0.05700 , 0.05700 , 0.05700 , 0.05700 , 0.05700 ,  &
       0.03700 , 0.03700 , 0.03700 , 0.04200 , 0.04200 , 0.04200 , 0.04200 , 0.04200 ,  &

       0.00400 , 0.00400 , 0.00400 , 0.00800 , 0.00800 , 0.00800 , 0.00800 , 0.00800 ,  &
       0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.13800 , 0.13800 , 0.13800 , 0.13800 , 0.13800 ,  &
       0.11100 , 0.11100 , 0.11100 , 0.12500 , 0.12500 , 0.12500 , 0.12500 , 0.12500 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00700 , 0.00700 , 0.00700 , 0.00700 , 0.00700 ,  &
       0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 ,  &

       0.01700 , 0.01700 , 0.01700 , 0.03400 , 0.03400 , 0.03400 , 0.03400 , 0.03400 ,  &
       0.02800 , 0.02800 , 0.02800 , 0.03100 , 0.03100 , 0.03100 , 0.03100 , 0.03100 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00700 , 0.00700 , 0.00700 , 0.00700 , 0.00700 ,  &
       0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 ,  &

       0.01700 , 0.01700 , 0.01700 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &
       0.02400 , 0.02400 , 0.02400 , 0.02700 , 0.02700 , 0.02700 , 0.02700 , 0.02700 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00300 , 0.00300 , 0.00300 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       0.01600 , 0.01600 , 0.01600 , 0.01900 , 0.01900 , 0.01900 , 0.01900 , 0.01900 ,  &
       0.01900 , 0.01900 , 0.01900 , 0.02200 , 0.02200 , 0.02200 , 0.02200 , 0.02200 ,  &

       0.00600 , 0.00600 , 0.00600 , 0.00700 , 0.00700 , 0.00700 , 0.00700 , 0.00700 ,  &
       0.00800 , 0.00800 , 0.00800 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &

       0.14400 , 0.14400 , 0.14400 , 0.08400 , 0.08400 , 0.08400 , 0.08400 , 0.08400 ,  &
       0.09600 , 0.09600 , 0.09600 , 0.09800 , 0.09800 , 0.09800 , 0.09800 , 0.09800 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.12000 , 0.12000 , 0.12000 , 0.05600 , 0.05600 , 0.05600 , 0.05600 , 0.05600 ,  &
       0.06700 , 0.06700 , 0.06700 , 0.07700 , 0.07700 , 0.07700 , 0.07700 , 0.07700 ,  &

       0.02400 , 0.02400 , 0.02400 , 0.02800 , 0.02800 , 0.02800 , 0.02800 , 0.02800 ,  &
       0.02900 , 0.02900 , 0.02900 , 0.03300 , 0.03300 , 0.03300 , 0.03300 , 0.03300 ,  &

       0.01200 , 0.01200 , 0.01200 , 0.01400 , 0.01400 , 0.01400 , 0.01400 , 0.01400 ,  &
       0.01400 , 0.01400 , 0.01400 , 0.01600 , 0.01600 , 0.01600 , 0.01600 , 0.01600 ,  &

       0.00800 , 0.00800 , 0.00800 , 0.00900 , 0.00900 , 0.00900 , 0.00900 , 0.00900 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &

       0.00500 , 0.00500 , 0.00500 , 0.00600 , 0.00600 , 0.00600 , 0.00600 , 0.00600 ,  &
       0.00600 , 0.00600 , 0.00600 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00800 , 0.00800 , 0.00800 , 0.00900 , 0.00900 , 0.00900 , 0.00900 , 0.00900 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &

       0.01200 , 0.01200 , 0.01200 , 0.01400 , 0.01400 , 0.01400 , 0.01400 , 0.01400 ,  &
       0.01400 , 0.01400 , 0.01400 , 0.01600 , 0.01600 , 0.01600 , 0.01600 , 0.01600 ,  &

       0.00800 , 0.00800 , 0.00800 , 0.00900 , 0.00900 , 0.00900 , 0.00900 , 0.00900 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &

       0.01600 , 0.01600 , 0.01600 , 0.01900 , 0.01900 , 0.01900 , 0.01900 , 0.01900 ,  &
       0.01900 , 0.01900 , 0.01900 , 0.02200 , 0.02200 , 0.02200 , 0.02200 , 0.02200 ,  &

       0.23400 , 0.23400 , 0.23400 , 0.27600 , 0.27600 , 0.27600 , 0.27600 , 0.27600 ,  &
       0.28500 , 0.28500 , 0.28500 , 0.22400 , 0.22400 , 0.22400 , 0.22400 , 0.22400 ,  &

       0.00800 , 0.00800 , 0.00800 , 0.00900 , 0.00900 , 0.00900 , 0.00900 , 0.00900 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01100 , 0.01100 , 0.01100 , 0.01100 , 0.01100 ,  &

       0.02400 , 0.02400 , 0.02400 , 0.02800 , 0.02800 , 0.02800 , 0.02800 , 0.02800 ,  &
       0.02900 , 0.02900 , 0.02900 , 0.03300 , 0.03300 , 0.03300 , 0.03300 , 0.03300 ,  &

       0.00400 , 0.00400 , 0.00400 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &
       0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.19900 , 0.19900 , 0.19900 , 0.13900 , 0.13900 , 0.13900 , 0.13900 , 0.13900 ,  &
       0.17200 , 0.17200 , 0.17200 , 0.16400 , 0.16400 , 0.16400 , 0.16400 , 0.16400 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.01200 , 0.01200 , 0.01200 , 0.01400 , 0.01400 , 0.01400 , 0.01400 , 0.01400 ,  &
       0.01400 , 0.01400 , 0.01400 , 0.01600 , 0.01600 , 0.01600 , 0.01600 , 0.01600 ,  &

       0.04000 , 0.04000 , 0.04000 , 0.04600 , 0.04600 , 0.04600 , 0.04600 , 0.04600 ,  &
       0.04800 , 0.04800 , 0.04800 , 0.05500 , 0.05500 , 0.05500 , 0.05500 , 0.05500 ,  &

       0.08000 , 0.08000 , 0.08000 , 0.18600 , 0.18600 , 0.18600 , 0.18600 , 0.18600 ,  &
       0.11500 , 0.11500 , 0.11500 , 0.10900 , 0.10900 , 0.10900 , 0.10900 , 0.10900 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &
       0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 , 0.00500 ,  &

       0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 ,  &
       0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 ,  &

       0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 ,  &
       0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 , 0.00000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 ,  &
       0.40000 , 0.40000 , 0.40000 , 0.25000 , 0.25000 , 0.25000 , 0.25000 , 0.25000 ,  &

       0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 , 0.40000 ,  &
       0.40000 , 0.40000 , 0.40000 , 0.25000 , 0.25000 , 0.25000 , 0.25000 , 0.25000 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &
       0.06000 , 0.06000 , 0.06000 , 0.15000 , 0.15000 , 0.15000 , 0.15000 , 0.15000 ,  &

       0.08000 , 0.08000 , 0.08000 , 0.08000 , 0.08000 , 0.08000 , 0.08000 , 0.08000 ,  &
       0.08000 , 0.08000 , 0.08000 , 0.20000 , 0.20000 , 0.20000 , 0.20000 , 0.20000 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &
       0.06000 , 0.06000 , 0.06000 , 0.15000 , 0.15000 , 0.15000 , 0.15000 , 0.15000 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 ,  &
       0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 , 0.02500 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 ,  &
       0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &
       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &

       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &
       1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 , 1.00000 ,  &

       0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 ,  &
       0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 , 0.24000 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 ,  &
       0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 , 0.58000 ,  &

       0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 ,  &
       0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 , 0.01500 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 ,  &
       0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 , 0.48000 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &
       0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 ,  &
       0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 , 0.05000 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &
       0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 , 0.00300 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &
       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &

       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &
       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &
       0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 , 0.01000 ,  &

       0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 ,  &
       0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 , 0.10000 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &
       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &
       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &
       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &

       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &
       0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 , 0.06000 ,  &

       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &
       0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 , 0.03000 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &
       0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 , 0.00200 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00300 , 0.00300 , 0.00300 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &

       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 ,  &
       0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 , 0.00100 /),&
             SHAPE=(/16,N_EFFS_SPC/)           ,  &
                    ORDER=(/1,2/)                 )

end module CONSTS_MEGAN
