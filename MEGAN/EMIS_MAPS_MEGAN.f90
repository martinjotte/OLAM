module EMIS_MAPS_MEGAN

!=======================================================================
!  SPC_SPCAT.EXT
!  This include file contains MEGAN speciated species and their MW.
!
!  MEGAN v2.1.0
!  INPUT version 200
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
!=======================================================================

      INTEGER,PARAMETER :: N_SPCA_SPC = 150        ! Number of speciated species
      CHARACTER*20   SPCA_SPC( N_SPCA_SPC )   ! speciated species name
      REAL           SPCA_MWT( N_SPCA_SPC )   ! Mechanism species molecular weight

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! _a  = alpha, _b  = beta, _c  = cis, _al = allo,
! _g  = gamma, _d  = delta, _t  = trans, _m  = methyl,
! _p  = para, _o  = ortho, _e  = ene, _ol = ol ,
! met = methyl, 2met= dimethyl, MBO = methylbutenol        ,
! 2s  = disulfide, s   = sulfide, OXD = oxide, ACT = acetate,
! PPPP= propenylpropyl       , DCTT= decatetraene         ,
! COTHER= acetaldehyde
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Isoprene
      DATA  SPCA_SPC(  1), SPCA_MWT(  1) / 'isoprene', 68.12  /

! MTP
      DATA  SPCA_SPC(  2), SPCA_MWT(  2) / 'myrcene ', 136.23 /
      DATA  SPCA_SPC(  3), SPCA_MWT(  3) / 'sabinene', 136.23 /
      DATA  SPCA_SPC(  4), SPCA_MWT(  4) / 'limonene', 136.23 /
      DATA  SPCA_SPC(  5), SPCA_MWT(  5) / 'carene_3', 136.23 /
      DATA  SPCA_SPC(  6), SPCA_MWT(  6) / 'ocimene_t_b   ', 136.23 /
      DATA  SPCA_SPC(  7), SPCA_MWT(  7) / 'pinene_b', 136.23 /
      DATA  SPCA_SPC(  8), SPCA_MWT(  8) / 'pinene_a', 136.23 /
! Other MT
      DATA  SPCA_SPC(  9), SPCA_MWT(  9) / 'A_2met_styrene  ', 132.20 /
      DATA  SPCA_SPC( 10), SPCA_MWT( 10) / 'cymene_p', 134.22 /
      DATA  SPCA_SPC( 11), SPCA_MWT( 11) / 'cymene_o', 134.22 /
      DATA  SPCA_SPC( 12), SPCA_MWT( 12) / 'phellandrene_a', 136.23 /
      DATA  SPCA_SPC( 13), SPCA_MWT( 13) / 'thujene_a ', 136.23 /
      DATA  SPCA_SPC( 14), SPCA_MWT( 14) / 'terpinene_a   ', 136.23 /
      DATA  SPCA_SPC( 15), SPCA_MWT( 15) / 'terpinene_g   ', 136.23 /
      DATA  SPCA_SPC( 16), SPCA_MWT( 16) / 'terpinolene   ', 136.23 /
      DATA  SPCA_SPC( 17), SPCA_MWT( 17) / 'phellandrene_b', 136.23 /
      DATA  SPCA_SPC( 18), SPCA_MWT( 18) / 'camphene', 136.23 /
      DATA  SPCA_SPC( 19), SPCA_MWT( 19) / 'bornene ', 136.23 /
      DATA  SPCA_SPC( 20), SPCA_MWT( 20) / 'fenchene_a', 136.23 /
      DATA  SPCA_SPC( 21), SPCA_MWT( 21) / 'ocimene_al', 136.23 /
      DATA  SPCA_SPC( 22), SPCA_MWT( 22) / 'ocimene_c_b   ', 136.23 /
      DATA  SPCA_SPC( 23), SPCA_MWT( 23) / 'tricyclene', 136.23 /
      DATA  SPCA_SPC( 24), SPCA_MWT( 24) / 'estragole ', 148.20 /
      DATA  SPCA_SPC( 25), SPCA_MWT( 25) / 'camphor ', 152.23 /
      DATA  SPCA_SPC( 26), SPCA_MWT( 26) / 'fenchone', 152.23 /
      DATA  SPCA_SPC( 27), SPCA_MWT( 27) / 'piperitone', 152.23 /
      DATA  SPCA_SPC( 28), SPCA_MWT( 28) / 'thujone_a ', 152.23 /
      DATA  SPCA_SPC( 29), SPCA_MWT( 29) / 'thujone_b ', 152.23 /
      DATA  SPCA_SPC( 30), SPCA_MWT( 30) / 'cineole_1_8   ', 154.25 /
      DATA  SPCA_SPC( 31), SPCA_MWT( 31) / 'borneol ', 154.25 /
      DATA  SPCA_SPC( 32), SPCA_MWT( 32) / 'linalool', 154.25 /
      DATA  SPCA_SPC( 33), SPCA_MWT( 33) / 'terpineol_4   ', 154.25 /
      DATA  SPCA_SPC( 34), SPCA_MWT( 34) / 'terpineol_a   ', 154.25 /
      DATA  SPCA_SPC( 35), SPCA_MWT( 35) / 'linalool_OXD_c', 170.25 /
      DATA  SPCA_SPC( 36), SPCA_MWT( 36) / 'linalool_OXD_t', 170.25 /
      DATA  SPCA_SPC( 37), SPCA_MWT( 37) / 'ionone_b', 192.30 /
      DATA  SPCA_SPC( 38), SPCA_MWT( 38) / 'bornyl_ACT', 196.29 /

! SQT
      DATA  SPCA_SPC( 39), SPCA_MWT( 39) / 'farnescene_a  ', 204.35 /
      DATA  SPCA_SPC( 40), SPCA_MWT( 40) / 'caryophyllene_b ', 204.35 /
! Other SQT
      DATA  SPCA_SPC( 41), SPCA_MWT( 41) / 'acoradiene', 204.35 /
      DATA  SPCA_SPC( 42), SPCA_MWT( 42) / 'aromadendrene ', 204.35 /
      DATA  SPCA_SPC( 43), SPCA_MWT( 43) / 'bergamotene_a ', 204.35 /
      DATA  SPCA_SPC( 44), SPCA_MWT( 44) / 'bergamotene_b ', 204.35 /
      DATA  SPCA_SPC( 45), SPCA_MWT( 45) / 'bisabolene_a  ', 204.35 /
      DATA  SPCA_SPC( 46), SPCA_MWT( 46) / 'bisabolene_b  ', 204.35 /
      DATA  SPCA_SPC( 47), SPCA_MWT( 47) / 'bourbonene_b  ', 204.35 /
      DATA  SPCA_SPC( 48), SPCA_MWT( 48) / 'cadinene_d', 204.35 /
      DATA  SPCA_SPC( 49), SPCA_MWT( 49) / 'cadinene_g', 204.35 /
      DATA  SPCA_SPC( 50), SPCA_MWT( 50) / 'cedrene_a ', 204.35 /
      DATA  SPCA_SPC( 51), SPCA_MWT( 51) / 'copaene_a ', 204.35 /
      DATA  SPCA_SPC( 52), SPCA_MWT( 52) / 'cubebene_a', 204.35 /
      DATA  SPCA_SPC( 53), SPCA_MWT( 53) / 'cubebene_b', 204.35 /
      DATA  SPCA_SPC( 54), SPCA_MWT( 54) / 'elemene_b ', 204.35 /
      DATA  SPCA_SPC( 55), SPCA_MWT( 55) / 'farnescene_b  ', 204.35 /
      DATA  SPCA_SPC( 56), SPCA_MWT( 56) / 'germacrene_B  ', 204.35 /
      DATA  SPCA_SPC( 57), SPCA_MWT( 57) / 'germacrene_D  ', 204.35 /
      DATA  SPCA_SPC( 58), SPCA_MWT( 58) / 'gurjunene_b   ', 204.35 /
      DATA  SPCA_SPC( 59), SPCA_MWT( 59) / 'humulene_a', 204.35 /
      DATA  SPCA_SPC( 60), SPCA_MWT( 60) / 'humulene_g', 204.35 /
      DATA  SPCA_SPC( 61), SPCA_MWT( 61) / 'isolongifolene', 204.35 /
      DATA  SPCA_SPC( 62), SPCA_MWT( 62) / 'longifolene   ', 204.35 /
      DATA  SPCA_SPC( 63), SPCA_MWT( 63) / 'longipinene   ', 204.35 /
      DATA  SPCA_SPC( 64), SPCA_MWT( 64) / 'muurolene_a   ', 204.35 /
      DATA  SPCA_SPC( 65), SPCA_MWT( 65) / 'muurolene_g   ', 204.35 /
      DATA  SPCA_SPC( 66), SPCA_MWT( 66) / 'selinene_b', 204.35 /
      DATA  SPCA_SPC( 67), SPCA_MWT( 67) / 'selinene_d', 204.35 /
      DATA  SPCA_SPC( 68), SPCA_MWT( 68) / 'nerolidol_c   ', 222.37 /
      DATA  SPCA_SPC( 69), SPCA_MWT( 69) / 'nerolidol_t   ', 222.37 /
      DATA  SPCA_SPC( 70), SPCA_MWT( 70) / 'cedrol  ', 222.37 /

! VOC
      DATA  SPCA_SPC( 71), SPCA_MWT( 71) / 'MBO_2m3e2ol   ', 86.13  /
      DATA  SPCA_SPC( 72), SPCA_MWT( 72) / 'methanol', 32.04  /
      DATA  SPCA_SPC( 73), SPCA_MWT( 73) / 'acetone ', 58.08  /
      DATA  SPCA_SPC( 74), SPCA_MWT( 74) / 'methane ', 16.04  /
! Ammonia, NO2, and NO
      DATA  SPCA_SPC( 75), SPCA_MWT( 75) / 'ammonia ', 17.03  /
      DATA  SPCA_SPC( 76), SPCA_MWT( 76) / 'nitrous_OXD   ', 44.01  /
      DATA  SPCA_SPC( 77), SPCA_MWT( 77) / 'nitric_OXD', 30.01  /
! Acetaldehyde + ethanol
      DATA  SPCA_SPC( 78), SPCA_MWT( 78) / 'acetaldehyde  ', 44.05  /
      DATA  SPCA_SPC( 79), SPCA_MWT( 79) / 'ethanol ', 46.07  /
! Formic acid + formaldehyde + acetic acid
      DATA  SPCA_SPC( 80), SPCA_MWT( 80) / 'formic_acid   ', 46.03  /
      DATA  SPCA_SPC( 81), SPCA_MWT( 81) / 'formaldehyde  ', 30.03  /
      DATA  SPCA_SPC( 82), SPCA_MWT( 82) / 'acetic_acid   ', 60.05  /
! Other VC
      DATA  SPCA_SPC( 83), SPCA_MWT( 83) / 'MBO_3m2e1ol   ', 86.13  /
      DATA  SPCA_SPC( 84), SPCA_MWT( 84) / 'MBO_3m3e1ol   ', 86.13  /
      DATA  SPCA_SPC( 85), SPCA_MWT( 85) / 'benzaldehyde  ', 106.12 /
      DATA  SPCA_SPC( 86), SPCA_MWT( 86) / 'butanone_2', 72.11  /
      DATA  SPCA_SPC( 87), SPCA_MWT( 87) / 'decanal ', 156.27 /
      DATA  SPCA_SPC( 88), SPCA_MWT( 88) / 'dodecene_1', 168.32 /
      DATA  SPCA_SPC( 89), SPCA_MWT( 89) / 'geranyl_acetone ', 194.31 /
      DATA  SPCA_SPC( 90), SPCA_MWT( 90) / 'heptanal', 114.19 /
      DATA  SPCA_SPC( 91), SPCA_MWT( 91) / 'heptane ', 100.20 /
      DATA  SPCA_SPC( 92), SPCA_MWT( 92) / 'hexane  ', 86.18  /
      DATA  SPCA_SPC( 93), SPCA_MWT( 93) / 'met_benzoate  ', 136.15 /
      DATA  SPCA_SPC( 94), SPCA_MWT( 94) / 'met_heptenone ', 126.20 /
      DATA  SPCA_SPC( 95), SPCA_MWT( 95) / 'neryl_acetone ', 194.31 /
      DATA  SPCA_SPC( 96), SPCA_MWT( 96) / 'nonanal ', 142.24 /
      DATA  SPCA_SPC( 97), SPCA_MWT( 97) / 'nonenal ', 140.22 /
      DATA  SPCA_SPC( 98), SPCA_MWT( 98) / 'octanal ', 128.21 /
      DATA  SPCA_SPC( 99), SPCA_MWT( 99) / 'octanol ', 130.23 /
      DATA  SPCA_SPC(100), SPCA_MWT(100) / 'octenol_1e3ol ', 128.21 /
      DATA  SPCA_SPC(101), SPCA_MWT(101) / 'oxopentanal   ', 100.12 /
      DATA  SPCA_SPC(102), SPCA_MWT(102) / 'pentane ', 72.15  /
      DATA  SPCA_SPC(103), SPCA_MWT(103) / 'phenyl_CCO', 120.15 /
      DATA  SPCA_SPC(104), SPCA_MWT(104) / 'pyruvic_acid  ', 88.06  /
      DATA  SPCA_SPC(105), SPCA_MWT(105) / 'terpinyl_ACT_a', 196.29 /
      DATA  SPCA_SPC(106), SPCA_MWT(106) / 'tetradecene_1 ', 196.37 /
      DATA  SPCA_SPC(107), SPCA_MWT(107) / 'toluene ', 92.14  /
      DATA  SPCA_SPC(108), SPCA_MWT(108) / 'carbon_monoxide ', 28.01  /
      DATA  SPCA_SPC(109), SPCA_MWT(109) / 'butene  ', 56.11  /
      DATA  SPCA_SPC(110), SPCA_MWT(110) / 'ethane  ', 30.07  /
      DATA  SPCA_SPC(111), SPCA_MWT(111) / 'ethene  ', 28.05  /
      DATA  SPCA_SPC(112), SPCA_MWT(112) / 'hydrogen_cyanide', 27.03  /
      DATA  SPCA_SPC(113), SPCA_MWT(113) / 'propane ', 44.10  /
      DATA  SPCA_SPC(114), SPCA_MWT(114) / 'propene ', 42.08  /
      DATA  SPCA_SPC(115), SPCA_MWT(115) / 'carbon_2s ', 76.14  /
      DATA  SPCA_SPC(116), SPCA_MWT(116) / 'carbonyl_s', 60.08  /
      DATA  SPCA_SPC(117), SPCA_MWT(117) / 'diallyl_2s', 146.28 /
      DATA  SPCA_SPC(118), SPCA_MWT(118) / 'A_2met_2s ', 94.20  /
      DATA  SPCA_SPC(119), SPCA_MWT(119) / 'A_2met_s  ', 62.14  /
      DATA  SPCA_SPC(120), SPCA_MWT(120) / 'met_chloride  ', 50.49  /
      DATA  SPCA_SPC(121), SPCA_MWT(121) / 'met_bromide   ', 94.94  /
      DATA  SPCA_SPC(122), SPCA_MWT(122) / 'met_iodide', 141.94 /
      DATA  SPCA_SPC(123), SPCA_MWT(123) / 'hydrogen_s', 34.08  /
      DATA  SPCA_SPC(124), SPCA_MWT(124) / 'met_mercaptan ', 48.11  /
      DATA  SPCA_SPC(125), SPCA_MWT(125) / 'met_propenyl_2s ', 120.24 /
      DATA  SPCA_SPC(126), SPCA_MWT(126) / 'PPPP_2s ', 148.29 /
      DATA  SPCA_SPC(127), SPCA_MWT(127) / 'A_2met_nonatriene',150.26 /
      DATA  SPCA_SPC(128), SPCA_MWT(128) / 'met_salicylate', 152.15 /
      DATA  SPCA_SPC(129), SPCA_MWT(129) / 'indole  ', 117.15 /
      DATA  SPCA_SPC(130), SPCA_MWT(130) / 'jasmone ', 164.24 /
      DATA  SPCA_SPC(131), SPCA_MWT(131) / 'met_jasmonate ', 224.30 /
      DATA  SPCA_SPC(132), SPCA_MWT(132) / 'A_3met_3DCTT', 218.38 /
      DATA  SPCA_SPC(133), SPCA_MWT(133) / 'hexanal ', 100.16 /
      DATA  SPCA_SPC(134), SPCA_MWT(134) / 'hexanol_1 ', 102.17 /
      DATA  SPCA_SPC(135), SPCA_MWT(135) / 'hexenal_c3', 98.14  /
      DATA  SPCA_SPC(136), SPCA_MWT(136) / 'hexenal_t2', 98.14  /
      DATA  SPCA_SPC(137), SPCA_MWT(137) / 'hexenol_c3', 100.16 /
      DATA  SPCA_SPC(138), SPCA_MWT(138) / 'hexenyl_ACT_c3', 142.20 /
      DATA  SPCA_SPC(139), SPCA_MWT(139) / 'homosalate  ', 131 /
      DATA  SPCA_SPC(140), SPCA_MWT(140) / 'Ehsalate ', 131 /
      DATA  SPCA_SPC(141), SPCA_MWT(141) / 'pentanal ', 133 /
      DATA  SPCA_SPC(142), SPCA_MWT(142) / 'heptanone', 94 /
      DATA  SPCA_SPC(143), SPCA_MWT(143) / 'anisole ', 85 /
      DATA  SPCA_SPC(144), SPCA_MWT(144) / 'verbenene ', 10 /
      DATA  SPCA_SPC(145), SPCA_MWT(145) / 'benzyl-acetate', 85  /
      DATA  SPCA_SPC(146), SPCA_MWT(146) / 'myrtenal', 32  /
      DATA  SPCA_SPC(147), SPCA_MWT(147) / 'benzyl-alcohol', 85 /
      DATA  SPCA_SPC(148), SPCA_MWT(148) / 'meta-cymenene', 10 /
      DATA  SPCA_SPC(149), SPCA_MWT(149) / 'ipsenol  ', 32 /
      DATA  SPCA_SPC(150), SPCA_MWT(150) / 'Napthalene ', 129 /
!=======================================================================
!  MAP_MGN20T138.EXT
!  This include file contains conversion table for MEGAN species to
!  134 species
!
!
!  MEGAN v2.1.0
!  INPUT version 200
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
!=======================================================================

      INTEGER,PARAMETER :: N_SMAP_SPC = 150   ! Number of map species

      CHARACTER*16   SPCA_NAM( N_SMAP_SPC )   ! speciated species name
      INTEGER        SPCA_MAP( N_SMAP_SPC )   ! speciated species name
                                              ! mapped to SPCAT_SPC.EXT
      CHARACTER*16   MG20_NAM( N_SMAP_SPC )   ! MEGAN species
      INTEGER        MG20_MAP( N_SMAP_SPC )   ! MEGAN species mapped to
                                              ! MGN_SPC.EXT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! _a  = alpha, _b  = beta, _c  = cis, _al = allo,
! _g  = gamma, _d  = delta, _t  = trans, _m  = methyl,
! _p  = para, _o  = ortho, _e  = ene, _ol = ol ,
! met = methyl, 2met= dimethyl, MBO = methylbutenol        ,
! 2s  = disulfide, s   = sulfide, OXD = oxide, ACT = acetate,
! PPPP= propenylpropyl       , DCTT= decatetraene         ,
! COTHER= acetaldehyde
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      DATA  SPCA_NAM(  1), SPCA_MAP(  1), MG20_NAM(  1), MG20_MAP(  1)  &
          / 'isoprene', 1     , 'ISOP    ', 1             /
      DATA  SPCA_NAM(  2), SPCA_MAP(  2), MG20_NAM(  2), MG20_MAP(  2)  &
          / 'myrcene ', 2     , 'MYRC    ', 2             /
      DATA  SPCA_NAM(  3), SPCA_MAP(  3), MG20_NAM(  3), MG20_MAP(  3)  &
          / 'sabinene', 3     , 'SABI    ', 3             /
      DATA  SPCA_NAM(  4), SPCA_MAP(  4), MG20_NAM(  4), MG20_MAP(  4)  &
          / 'limonene', 4     , 'LIMO    ', 4             /
      DATA  SPCA_NAM(  5), SPCA_MAP(  5), MG20_NAM(  5), MG20_MAP(  5)  &
          / 'carene_3', 5     , '3CAR    ', 5             /
      DATA  SPCA_NAM(  6), SPCA_MAP(  6), MG20_NAM(  6), MG20_MAP(  6)  &
          / 'ocimene_t_b   ', 6     , 'OCIM    ', 6             /
      DATA  SPCA_NAM(  7), SPCA_MAP(  7), MG20_NAM(  7), MG20_MAP(  7)  &
          / 'pinene_b', 7     , 'BPIN    ', 7             /
      DATA  SPCA_NAM(  8), SPCA_MAP(  8), MG20_NAM(  8), MG20_MAP(  8)  &
          / 'pinene_a', 8     , 'APIN    ', 8             /
      DATA  SPCA_NAM(  9), SPCA_MAP(  9), MG20_NAM(  9), MG20_MAP(  9)  &
          / 'A_2met_styrene  ', 9     , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 10), SPCA_MAP( 10), MG20_NAM( 10), MG20_MAP( 10)  &
          / 'cymene_p', 10    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 11), SPCA_MAP( 11), MG20_NAM( 11), MG20_MAP( 11)  &
          / 'cymene_o', 11    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 12), SPCA_MAP( 12), MG20_NAM( 12), MG20_MAP( 12)  &
          / 'phellandrene_a', 12    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 13), SPCA_MAP( 13), MG20_NAM( 13), MG20_MAP( 13)  &
          / 'thujene_a ', 13    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 14), SPCA_MAP( 14), MG20_NAM( 14), MG20_MAP( 14)  &
          / 'terpinene_a   ', 14    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 15), SPCA_MAP( 15), MG20_NAM( 15), MG20_MAP( 15)  &
          / 'terpinene_g   ', 15    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 16), SPCA_MAP( 16), MG20_NAM( 16), MG20_MAP( 16)  &
          / 'terpinolene   ', 16    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 17), SPCA_MAP( 17), MG20_NAM( 17), MG20_MAP( 17)  &
          / 'phellandrene_b', 17    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 18), SPCA_MAP( 18), MG20_NAM( 18), MG20_MAP( 18)  &
          / 'camphene', 18    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 19), SPCA_MAP( 19), MG20_NAM( 19), MG20_MAP( 19)  &
          / 'bornene ', 19    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 20), SPCA_MAP( 20), MG20_NAM( 20), MG20_MAP( 20)  &
          / 'fenchene_a', 20    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 21), SPCA_MAP( 21), MG20_NAM( 21), MG20_MAP( 21)  &
          / 'ocimene_al', 21    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 22), SPCA_MAP( 22), MG20_NAM( 22), MG20_MAP( 22)  &
          / 'ocimene_c_b   ', 22    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 23), SPCA_MAP( 23), MG20_NAM( 23), MG20_MAP( 23)  &
          / 'tricyclene', 23    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 24), SPCA_MAP( 24), MG20_NAM( 24), MG20_MAP( 24)  &
          / 'estragole ', 24    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 25), SPCA_MAP( 25), MG20_NAM( 25), MG20_MAP( 25)  &
          / 'camphor ', 25    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 26), SPCA_MAP( 26), MG20_NAM( 26), MG20_MAP( 26)  &
          / 'fenchone', 26    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 27), SPCA_MAP( 27), MG20_NAM( 27), MG20_MAP( 27)  &
          / 'piperitone', 27    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 28), SPCA_MAP( 28), MG20_NAM( 28), MG20_MAP( 28)  &
          / 'thujone_a ', 28    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 29), SPCA_MAP( 29), MG20_NAM( 29), MG20_MAP( 29)  &
          / 'thujone_b ', 29    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 30), SPCA_MAP( 30), MG20_NAM( 30), MG20_MAP( 30)  &
          / 'cineole_1_8   ', 30    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 31), SPCA_MAP( 31), MG20_NAM( 31), MG20_MAP( 31)  &
          / 'borneol ', 31    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 32), SPCA_MAP( 32), MG20_NAM( 32), MG20_MAP( 32)  &
          / 'linalool', 32    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 33), SPCA_MAP( 33), MG20_NAM( 33), MG20_MAP( 33)  &
          / 'terpineol_4   ', 33    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 34), SPCA_MAP( 34), MG20_NAM( 34), MG20_MAP( 34)  &
          / 'terpineol_a   ', 34    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 35), SPCA_MAP( 35), MG20_NAM( 35), MG20_MAP( 35)  &
          / 'linalool_OXD_c', 35    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 36), SPCA_MAP( 36), MG20_NAM( 36), MG20_MAP( 36)  &
          / 'linalool_OXD_t', 36    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 37), SPCA_MAP( 37), MG20_NAM( 37), MG20_MAP( 37)  &
          / 'ionone_b', 37    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 38), SPCA_MAP( 38), MG20_NAM( 38), MG20_MAP( 38)  &
          / 'bornyl_ACT', 38    , 'OMTP    ', 9             /
      DATA  SPCA_NAM( 39), SPCA_MAP( 39), MG20_NAM( 39), MG20_MAP( 39)  &
          / 'farnescene_a  ', 39    , 'FARN    ', 10            /
      DATA  SPCA_NAM( 40), SPCA_MAP( 40), MG20_NAM( 40), MG20_MAP( 40)  &
          / 'caryophyllene_b ', 40    , 'BCAR    ', 11            /
      DATA  SPCA_NAM( 41), SPCA_MAP( 41), MG20_NAM( 41), MG20_MAP( 41)  &
          / 'acoradiene', 41    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 42), SPCA_MAP( 42), MG20_NAM( 42), MG20_MAP( 42)  &
          / 'aromadendrene ', 42    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 43), SPCA_MAP( 43), MG20_NAM( 43), MG20_MAP( 43)  &
          / 'bergamotene_a ', 43    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 44), SPCA_MAP( 44), MG20_NAM( 44), MG20_MAP( 44)  &
          / 'bergamotene_b ', 44    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 45), SPCA_MAP( 45), MG20_NAM( 45), MG20_MAP( 45)  &
          / 'bisabolene_a  ', 45    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 46), SPCA_MAP( 46), MG20_NAM( 46), MG20_MAP( 46)  &
          / 'bisabolene_b  ', 46    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 47), SPCA_MAP( 47), MG20_NAM( 47), MG20_MAP( 47)  &
          / 'bourbonene_b  ', 47    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 48), SPCA_MAP( 48), MG20_NAM( 48), MG20_MAP( 48)  &
          / 'cadinene_d', 48    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 49), SPCA_MAP( 49), MG20_NAM( 49), MG20_MAP( 49)  &
          / 'cadinene_g', 49    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 50), SPCA_MAP( 50), MG20_NAM( 50), MG20_MAP( 50)  &
          / 'cedrene_a ', 50    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 51), SPCA_MAP( 51), MG20_NAM( 51), MG20_MAP( 51)  &
          / 'copaene_a ', 51    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 52), SPCA_MAP( 52), MG20_NAM( 52), MG20_MAP( 52)  &
          / 'cubebene_a', 52    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 53), SPCA_MAP( 53), MG20_NAM( 53), MG20_MAP( 53)  &
          / 'cubebene_b', 53    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 54), SPCA_MAP( 54), MG20_NAM( 54), MG20_MAP( 54)  &
          / 'elemene_b ', 54    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 55), SPCA_MAP( 55), MG20_NAM( 55), MG20_MAP( 55)  &
          / 'farnescene_b  ', 55    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 56), SPCA_MAP( 56), MG20_NAM( 56), MG20_MAP( 56)  &
          / 'germacrene_B  ', 56    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 57), SPCA_MAP( 57), MG20_NAM( 57), MG20_MAP( 57)  &
          / 'germacrene_D  ', 57    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 58), SPCA_MAP( 58), MG20_NAM( 58), MG20_MAP( 58)  &
          / 'gurjunene_b   ', 58    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 59), SPCA_MAP( 59), MG20_NAM( 59), MG20_MAP( 59)  &
          / 'humulene_a', 59    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 60), SPCA_MAP( 60), MG20_NAM( 60), MG20_MAP( 60)  &
          / 'humulene_g', 60    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 61), SPCA_MAP( 61), MG20_NAM( 61), MG20_MAP( 61)  &
          / 'isolongifolene', 61    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 62), SPCA_MAP( 62), MG20_NAM( 62), MG20_MAP( 62)  &
          / 'longifolene   ', 62    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 63), SPCA_MAP( 63), MG20_NAM( 63), MG20_MAP( 63)  &
          / 'longipinene   ', 63    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 64), SPCA_MAP( 64), MG20_NAM( 64), MG20_MAP( 64)  &
          / 'muurolene_a   ', 64    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 65), SPCA_MAP( 65), MG20_NAM( 65), MG20_MAP( 65)  &
          / 'muurolene_g   ', 65    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 66), SPCA_MAP( 66), MG20_NAM( 66), MG20_MAP( 66)  &
          / 'selinene_b', 66    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 67), SPCA_MAP( 67), MG20_NAM( 67), MG20_MAP( 67)  &
          / 'selinene_d', 67    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 68), SPCA_MAP( 68), MG20_NAM( 68), MG20_MAP( 68)  &
          / 'nerolidol_c   ', 68    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 69), SPCA_MAP( 69), MG20_NAM( 69), MG20_MAP( 69)  &
          / 'nerolidol_t   ', 69    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 70), SPCA_MAP( 70), MG20_NAM( 70), MG20_MAP( 70)  &
          / 'cedrol  ', 70    , 'OSQT    ', 12            /
      DATA  SPCA_NAM( 71), SPCA_MAP( 71), MG20_NAM( 71), MG20_MAP( 71)  &
          / 'MBO_2m3e2ol   ', 71    , 'MBO     ', 13            /
      DATA  SPCA_NAM( 72), SPCA_MAP( 72), MG20_NAM( 72), MG20_MAP( 72)  &
          / 'methanol', 72    , 'MEOH    ', 14            /
      DATA  SPCA_NAM( 73), SPCA_MAP( 73), MG20_NAM( 73), MG20_MAP( 73)  &
          / 'acetone ', 73    , 'ACTO    ', 15            /
      DATA  SPCA_NAM( 74), SPCA_MAP( 74), MG20_NAM( 74), MG20_MAP( 74)  &
          / 'methane ', 74    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 75), SPCA_MAP( 75), MG20_NAM( 75), MG20_MAP( 75)  &
          / 'ammonia ', 75    , 'NO      ', 17            /
      DATA  SPCA_NAM( 76), SPCA_MAP( 76), MG20_NAM( 76), MG20_MAP( 76)  &
          / 'nitrous_OXD   ', 76    , 'NO      ', 17            /
      DATA  SPCA_NAM( 77), SPCA_MAP( 77), MG20_NAM( 77), MG20_MAP( 77)  &
          / 'nitric_OXD', 77    , 'NO      ', 17            /
      DATA  SPCA_NAM( 78), SPCA_MAP( 78), MG20_NAM( 78), MG20_MAP( 78)  &
          / 'acetaldehyde  ', 78    , 'BIDIR   ', 18            /
      DATA  SPCA_NAM( 79), SPCA_MAP( 79), MG20_NAM( 79), MG20_MAP( 79)  &
          / 'ethanol ', 79    , 'BIDIR   ', 18            /
      DATA  SPCA_NAM( 80), SPCA_MAP( 80), MG20_NAM( 80), MG20_MAP( 80)  &
          / 'formic_acid   ', 80    , 'BIDIR   ', 18           /
      DATA  SPCA_NAM( 81), SPCA_MAP( 81), MG20_NAM( 81), MG20_MAP( 81)  &
          / 'formaldehyde  ', 81    , 'BIDIR   ', 18            /
      DATA  SPCA_NAM( 82), SPCA_MAP( 82), MG20_NAM( 82), MG20_MAP( 82)  &
          / 'acetic_acid   ', 82    , 'BIDIR   ', 18            /
      DATA  SPCA_NAM( 83), SPCA_MAP( 83), MG20_NAM( 83), MG20_MAP( 83)  &
          / 'MBO_3m2e1ol   ', 83    , 'OTHER    ', 20            /
      DATA  SPCA_NAM( 84), SPCA_MAP( 84), MG20_NAM( 84), MG20_MAP( 84)  &
          / 'MBO_3m3e1ol   ', 84    , 'OTHER    ', 20            /
      DATA  SPCA_NAM( 85), SPCA_MAP( 85), MG20_NAM( 85), MG20_MAP( 85)  &
          / 'benzaldehyde  ', 85    , 'OTHER    ', 20            /
      DATA  SPCA_NAM( 86), SPCA_MAP( 86), MG20_NAM( 86), MG20_MAP( 86)  &
          / 'butanone_2', 86    , 'OTHER    ', 20            /
      DATA  SPCA_NAM( 87), SPCA_MAP( 87), MG20_NAM( 87), MG20_MAP( 87)  &
          / 'decanal ', 87    , 'OTHER    ', 20            /
      DATA  SPCA_NAM( 88), SPCA_MAP( 88), MG20_NAM( 88), MG20_MAP( 88)  &
          / 'dodecene_1', 88    , 'OTHER    ', 20            /
      DATA  SPCA_NAM( 89), SPCA_MAP( 89), MG20_NAM( 89), MG20_MAP( 89)  &
          / 'geranyl_acetone ', 89    , 'OTHER  ', 20            /
      DATA  SPCA_NAM( 90), SPCA_MAP( 90), MG20_NAM( 90), MG20_MAP( 90)  &
          / 'heptanal', 90    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 91), SPCA_MAP( 91), MG20_NAM( 91), MG20_MAP( 91)  &
          / 'heptane ', 91    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 92), SPCA_MAP( 92), MG20_NAM( 92), MG20_MAP( 92)  &
          / 'hexane  ', 92    , 'OTHER      ', 20            /
      DATA  SPCA_NAM( 93), SPCA_MAP( 93), MG20_NAM( 93), MG20_MAP( 93)  &
          / 'met_benzoate  ', 93    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 94), SPCA_MAP( 94), MG20_NAM( 94), MG20_MAP( 94)  &
          / 'met_heptenone ', 94    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 95), SPCA_MAP( 95), MG20_NAM( 95), MG20_MAP( 95)  &
          / 'neryl_acetone ', 95    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 96), SPCA_MAP( 96), MG20_NAM( 96), MG20_MAP( 96)  &
          / 'nonanal ', 96    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 97), SPCA_MAP( 97), MG20_NAM( 97), MG20_MAP( 97)  &
          / 'nonenal ', 97    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 98), SPCA_MAP( 98), MG20_NAM( 98), MG20_MAP( 98)  &
          / 'octanal ', 98    , 'OTHER     ', 20            /
      DATA  SPCA_NAM( 99), SPCA_MAP( 99), MG20_NAM( 99), MG20_MAP( 99)  &
          / 'octanol ', 99    , 'OTHER     ', 20            /
      DATA  SPCA_NAM(100), SPCA_MAP(100), MG20_NAM(100), MG20_MAP(100)  &
          / 'octenol_1e3ol ', 100   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(101), SPCA_MAP(101), MG20_NAM(101), MG20_MAP(101)  &
          / 'oxopentanal   ', 101   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(102), SPCA_MAP(102), MG20_NAM(102), MG20_MAP(102)  &
          / 'pentane ', 102   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(103), SPCA_MAP(103), MG20_NAM(103), MG20_MAP(103)  &
          / 'phenyl_CCO', 103   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(104), SPCA_MAP(104), MG20_NAM(104), MG20_MAP(104)  &
          / 'pyruvic_acid  ', 104   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(105), SPCA_MAP(105), MG20_NAM(105), MG20_MAP(105)  &
          / 'terpinyl_ACT_a', 105   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(106), SPCA_MAP(106), MG20_NAM(106), MG20_MAP(106)  &
          / 'tetradecene_1 ', 106   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(107), SPCA_MAP(107), MG20_NAM(107), MG20_MAP(107)  &
          / 'toluene ', 107   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(108), SPCA_MAP(108), MG20_NAM(108), MG20_MAP(108)  &
          / 'carbon_monoxide ', 108   , 'CO     ', 16            /
      DATA  SPCA_NAM(109), SPCA_MAP(109), MG20_NAM(109), MG20_MAP(109)  &
          / 'butene  ', 109   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(110), SPCA_MAP(110), MG20_NAM(110), MG20_MAP(110)  &
          / 'ethane  ', 110   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(111), SPCA_MAP(111), MG20_NAM(111), MG20_MAP(111)  &
          / 'ethene  ', 111   , 'STRESS     ', 20            /
      DATA  SPCA_NAM(112), SPCA_MAP(112), MG20_NAM(112), MG20_MAP(112)  &
          / 'hydrogen_cyanide', 112   , 'STRESS     ', 20            /
      DATA  SPCA_NAM(113), SPCA_MAP(113), MG20_NAM(113), MG20_MAP(113)  &
          / 'propane ', 113   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(114), SPCA_MAP(114), MG20_NAM(114), MG20_MAP(114)  &
          / 'propene ', 114   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(115), SPCA_MAP(115), MG20_NAM(115), MG20_MAP(115)  &
          / 'carbon_2s ', 115   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(116), SPCA_MAP(116), MG20_NAM(116), MG20_MAP(116)  &
          / 'carbonyl_s', 116   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(117), SPCA_MAP(117), MG20_NAM(117), MG20_MAP(117)  &
          / 'diallyl_2s', 117   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(118), SPCA_MAP(118), MG20_NAM(118), MG20_MAP(118)  &
          / 'A_2met_2s ', 118   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(119), SPCA_MAP(119), MG20_NAM(119), MG20_MAP(119)  &
          / 'A_2met_s  ', 119   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(120), SPCA_MAP(120), MG20_NAM(120), MG20_MAP(120)  &
          / 'met_chloride  ', 120   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(121), SPCA_MAP(121), MG20_NAM(121), MG20_MAP(121)  &
          / 'met_bromide   ', 121   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(122), SPCA_MAP(122), MG20_NAM(122), MG20_MAP(122)  &
          / 'met_iodide', 122   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(123), SPCA_MAP(123), MG20_NAM(123), MG20_MAP(123)  &
          / 'hydrogen_s', 123   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(124), SPCA_MAP(124), MG20_NAM(124), MG20_MAP(124)  &
          / 'met_mercaptan ', 124   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(125), SPCA_MAP(125), MG20_NAM(125), MG20_MAP(125)  &
          / 'met_propenyl_2s ', 125   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(126), SPCA_MAP(126), MG20_NAM(126), MG20_MAP(126)  &
          / 'PPPP_2s ', 126   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(127), SPCA_MAP(127), MG20_NAM(127), MG20_MAP(127)  &
          / 'A_2met_nonatriene ', 127   , 'STRESS     ', 20            /
      DATA  SPCA_NAM(128), SPCA_MAP(128), MG20_NAM(128), MG20_MAP(128)  &
          / 'met_salicylate', 128   , 'STRESS     ', 20            /
      DATA  SPCA_NAM(129), SPCA_MAP(129), MG20_NAM(129), MG20_MAP(129)  &
          / 'indole  ', 129   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(130), SPCA_MAP(130), MG20_NAM(130), MG20_MAP(130)  &
          / 'jasmone ', 130   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(131), SPCA_MAP(131), MG20_NAM(131), MG20_MAP(131)  &
          / 'met_jasmonate ', 131   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(132), SPCA_MAP(132), MG20_NAM(132), MG20_MAP(132)  &
          / '3met_3DCTT', 132   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(133), SPCA_MAP(133), MG20_NAM(133), MG20_MAP(133)  &
          / 'hexanal ', 133   , 'STRESS    ', 19            /
      DATA  SPCA_NAM(134), SPCA_MAP(134), MG20_NAM(134), MG20_MAP(134)  &
          / 'hexanol_1 ', 134   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(135), SPCA_MAP(135), MG20_NAM(135), MG20_MAP(135)  &
          / 'hexenal_c3', 135   , 'STRESS     ', 19            /
      DATA  SPCA_NAM(136), SPCA_MAP(136), MG20_NAM(136), MG20_MAP(136)  &
          / 'hexenal_t2', 136   , 'STRESS    ', 19            /
      DATA  SPCA_NAM(137), SPCA_MAP(137), MG20_NAM(137), MG20_MAP(137)  &
          / 'hexenol_c3', 137   , 'STRESS    ', 19            /
      DATA  SPCA_NAM(138), SPCA_MAP(138), MG20_NAM(138), MG20_MAP(138)  &
          / 'hexenyl_ACT_c3', 138   , 'STRESS     ', 19       /
      DATA  SPCA_NAM(139), SPCA_MAP(139), MG20_NAM(139), MG20_MAP(139)  &
          / 'homosalate  ', 139   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(140), SPCA_MAP(140), MG20_NAM(140), MG20_MAP(140)  &
          / 'Ehsalate ', 140   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(141), SPCA_MAP(141), MG20_NAM(141), MG20_MAP(141)  &
          / 'pentanal ', 141   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(142), SPCA_MAP(142), MG20_NAM(142), MG20_MAP(142)  &
          / 'heptanone', 142   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(143), SPCA_MAP(143), MG20_NAM(143), MG20_MAP(143)  &
          / 'anisole ', 143   , 'OTHER    ', 20            /
      DATA  SPCA_NAM(144), SPCA_MAP(144), MG20_NAM(144), MG20_MAP(144)  &
          / 'verbenene ', 144   , 'OMTP     ', 9            /
      DATA  SPCA_NAM(145), SPCA_MAP(145), MG20_NAM(145), MG20_MAP(145)  &
          / 'benzyl-acetate', 145   , 'OTHER     ', 20            /
      DATA  SPCA_NAM(146), SPCA_MAP(146), MG20_NAM(146), MG20_MAP(146)  &
          / 'myrtenal', 146   , 'OMTP    ',  9            /
      DATA  SPCA_NAM(147), SPCA_MAP(147), MG20_NAM(147), MG20_MAP(147)  &
          / 'benzyl-alcohol', 147   , 'OTHER    ', 20            /
      DATA  SPCA_NAM(148), SPCA_MAP(148), MG20_NAM(148), MG20_MAP(148)  &
          / 'meta-cymenene', 148   , 'OMTP    ',  9        /
      DATA  SPCA_NAM(149), SPCA_MAP(149), MG20_NAM(149), MG20_MAP(149)  &
          / 'ipsenol', 149   , 'OMTP    ',  9             /
      DATA  SPCA_NAM(150), SPCA_MAP(150), MG20_NAM(150), MG20_MAP(150)  &
          / 'Napthalene', 150   , 'OTHER     ', 20            /

!=======================================================================
!  SPC_CB05.EXT
!  This include file contains CB05 (CMAQ/CAMx) species and their MW.
!
!
!  Mechanism Name: CB05 (CMAQ/CAMx)
!  MEGAN v2.10
!  INPUT version x.x
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  bkoo         04/13/07 - Created
!=======================================================================

      CHARACTER*16   SPC_CB05MECH
      PARAMETER     (SPC_CB05MECH = 'CB05            ')

      INTEGER        N_CB05_SPC
      PARAMETER     (N_CB05_SPC = 20)

      CHARACTER*16, target :: MECH_SPC_CB05( N_CB05_SPC )  ! Mechanism species name
      REAL                    MECH_MWT_CB05( N_CB05_SPC )  ! Mechanism species molecular weight

      DATA  MECH_SPC_CB05(  1), MECH_MWT_CB05(  1) / 'ISOP            ',  80.00  /
      DATA  MECH_SPC_CB05(  2), MECH_MWT_CB05(  2) / 'TERP            ', 160.00  /
      DATA  MECH_SPC_CB05(  3), MECH_MWT_CB05(  3) / 'PAR             ',  16.00  /
      DATA  MECH_SPC_CB05(  4), MECH_MWT_CB05(  4) / 'XYL             ', 128.00  /
      DATA  MECH_SPC_CB05(  5), MECH_MWT_CB05(  5) / 'OLE             ',  32.00  /
      DATA  MECH_SPC_CB05(  6), MECH_MWT_CB05(  6) / 'NR              ',  16.00  /
      DATA  MECH_SPC_CB05(  7), MECH_MWT_CB05(  7) / 'MEOH            ',  16.00  /
      DATA  MECH_SPC_CB05(  8), MECH_MWT_CB05(  8) / 'CH4             ',  16.00  /
      DATA  MECH_SPC_CB05(  9), MECH_MWT_CB05(  9) / 'NH3             ',  17.00  /
      DATA  MECH_SPC_CB05( 10), MECH_MWT_CB05( 10) / 'NO              ',  46.00  /
      DATA  MECH_SPC_CB05( 11), MECH_MWT_CB05( 11) / 'ALD2            ',  32.00  /
      DATA  MECH_SPC_CB05( 12), MECH_MWT_CB05( 12) / 'ETOH            ',  32.00  /
      DATA  MECH_SPC_CB05( 13), MECH_MWT_CB05( 13) / 'FORM            ',  16.00  /
      DATA  MECH_SPC_CB05( 14), MECH_MWT_CB05( 14) / 'ALDX            ',  32.00  /
      DATA  MECH_SPC_CB05( 15), MECH_MWT_CB05( 15) / 'TOL             ', 112.00  /
      DATA  MECH_SPC_CB05( 16), MECH_MWT_CB05( 16) / 'IOLE            ',  64.00  /
      DATA  MECH_SPC_CB05( 17), MECH_MWT_CB05( 17) / 'CO              ',  28.00  /
      DATA  MECH_SPC_CB05( 18), MECH_MWT_CB05( 18) / 'ETHA            ',  32.00  /
      DATA  MECH_SPC_CB05( 19), MECH_MWT_CB05( 19) / 'ETH             ',  28.00  /
      DATA  MECH_SPC_CB05( 20), MECH_MWT_CB05( 20) / 'SESQ            ', 204.35  /

!=======================================================================
!  MAP_CB05_CV2CB05.EXT
!  This include file contains conversion table for 150 speciated species
!  to CB05 (CMAQ/CAMx) species
!
!
!  MEGAN v2.10
!  INPUT version x.x
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  bkoo         04/13/07 - Created
!  Tan          07/18/11 - Updated for MEGANv2.10
!=======================================================================

      CHARACTER*16  MAP_CB05MECH
      PARAMETER    (MAP_CB05MECH = 'CB05            ')

      INTEGER        N_CB05
      PARAMETER     (N_CB05 = 204)        ! Number of map species

      CHARACTER*16   SPMH_NAM_CB05( N_CB05 )   ! speciated species name
      INTEGER, target :: SPMH_MAP_CB05( N_CB05 )   ! speciated species name
                                              ! mapped to SPC_SPCAT.EXT
      CHARACTER*16   MECH_NAM_CB05( N_CB05 )   ! mechanism species
      INTEGER, target :: MECH_MAP_CB05( N_CB05 )   ! mechanism species mapped
                                              ! to SPC_CB4Q.EXT
      REAL, target    :: CONV_FAC_CB05( N_CB05 )   ! conversion factor


      DATA  SPMH_NAM_CB05(  1)     , SPMH_MAP_CB05(  1)  &
          / 'isoprene        ', 1             /
      DATA  MECH_NAM_CB05(  1)     , MECH_MAP_CB05(  1)  &
          / 'ISOP            ', 1             /
      DATA  CONV_FAC_CB05(  1)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  2)     , SPMH_MAP_CB05(  2)  &
          / 'myrcene         ', 2             /
      DATA  MECH_NAM_CB05(  2)     , MECH_MAP_CB05(  2)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  2)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  3)     , SPMH_MAP_CB05(  3)  &
          / 'sabinene        ', 3             /
      DATA  MECH_NAM_CB05(  3)     , MECH_MAP_CB05(  3)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  3)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  4)     , SPMH_MAP_CB05(  4)  &
          / 'limonene        ', 4             /
      DATA  MECH_NAM_CB05(  4)     , MECH_MAP_CB05(  4)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  4)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  5)     , SPMH_MAP_CB05(  5)  &
          / 'carene_3        ', 5             /
      DATA  MECH_NAM_CB05(  5)     , MECH_MAP_CB05(  5)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  5)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  6)     , SPMH_MAP_CB05(  6)  &
          / 'ocimene_t_b     ', 6             /
      DATA  MECH_NAM_CB05(  6)     , MECH_MAP_CB05(  6)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  6)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  7)     , SPMH_MAP_CB05(  7)  &
          / 'pinene_b        ', 7             /
      DATA  MECH_NAM_CB05(  7)     , MECH_MAP_CB05(  7)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  7)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  8)     , SPMH_MAP_CB05(  8)  &
          / 'pinene_a        ', 8             /
      DATA  MECH_NAM_CB05(  8)     , MECH_MAP_CB05(  8)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(  8)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(  9)     , SPMH_MAP_CB05(  9)  &
          / '2met_styrene    ', 9             /
      DATA  MECH_NAM_CB05(  9)     , MECH_MAP_CB05(  9)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(  9)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 10)     , SPMH_MAP_CB05( 10)  &
          / '2met_styrene    ', 9             /
      DATA  MECH_NAM_CB05( 10)     , MECH_MAP_CB05( 10)  &
          / 'XYL             ', 4             /
      DATA  CONV_FAC_CB05( 10)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 11)     , SPMH_MAP_CB05( 11)  &
          / '2met_styrene    ', 9             /
      DATA  MECH_NAM_CB05( 11)     , MECH_MAP_CB05( 11)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05( 11)  &
          / 0.50          /

      DATA  SPMH_NAM_CB05( 12)     , SPMH_MAP_CB05( 12)  &
          / 'cymene_p        ', 10            /
      DATA  MECH_NAM_CB05( 12)     , MECH_MAP_CB05( 12)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 12)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05( 13)     , SPMH_MAP_CB05( 13)  &
          / 'cymene_p        ', 10            /
      DATA  MECH_NAM_CB05( 13)     , MECH_MAP_CB05( 13)  &
          / 'XYL             ', 4             /
      DATA  CONV_FAC_CB05( 13)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 14)     , SPMH_MAP_CB05( 14)  &
          / 'cymene_o        ', 11            /
      DATA  MECH_NAM_CB05( 14)     , MECH_MAP_CB05( 14)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 14)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05( 15)     , SPMH_MAP_CB05( 15)  &
          / 'cymene_o        ', 11            /
      DATA  MECH_NAM_CB05( 15)     , MECH_MAP_CB05( 15)  &
          / 'XYL             ', 4             /
      DATA  CONV_FAC_CB05( 15)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 16)     , SPMH_MAP_CB05( 16)  &
          / 'phellandrene_a  ', 12            /
      DATA  MECH_NAM_CB05( 16)     , MECH_MAP_CB05( 16)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 16)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 17)     , SPMH_MAP_CB05( 17)  &
          / 'thujene_a       ', 13            /
      DATA  MECH_NAM_CB05( 17)     , MECH_MAP_CB05( 17)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 17)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 18)     , SPMH_MAP_CB05( 18)  &
          / 'terpinene_a     ', 14            /
      DATA  MECH_NAM_CB05( 18)     , MECH_MAP_CB05( 18)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 18)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 19)     , SPMH_MAP_CB05( 19)  &
          / 'terpinene_g     ', 15            /
      DATA  MECH_NAM_CB05( 19)     , MECH_MAP_CB05( 19)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 19)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 20)     , SPMH_MAP_CB05( 20)  &
          / 'terpinolene     ', 16            /
      DATA  MECH_NAM_CB05( 20)     , MECH_MAP_CB05( 20)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 20)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 21)     , SPMH_MAP_CB05( 21)  &
          / 'phellandrene_b  ', 17            /
      DATA  MECH_NAM_CB05( 21)     , MECH_MAP_CB05( 21)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 21)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 22)     , SPMH_MAP_CB05( 22)  &
          / 'camphene        ', 18            /
      DATA  MECH_NAM_CB05( 22)     , MECH_MAP_CB05( 22)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 22)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 23)     , SPMH_MAP_CB05( 23)  &
          / 'bornene         ', 19            /
      DATA  MECH_NAM_CB05( 23)     , MECH_MAP_CB05( 23)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 23)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 24)     , SPMH_MAP_CB05( 24)  &
          / 'fenchene_a      ', 20            /
      DATA  MECH_NAM_CB05( 24)     , MECH_MAP_CB05( 24)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 24)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 25)     , SPMH_MAP_CB05( 25)  &
          / 'ocimene_al      ', 21            /
      DATA  MECH_NAM_CB05( 25)     , MECH_MAP_CB05( 25)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 25)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 26)     , SPMH_MAP_CB05( 26)  &
          / 'ocimene_c_b     ', 22            /
      DATA  MECH_NAM_CB05( 26)     , MECH_MAP_CB05( 26)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 26)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 27)     , SPMH_MAP_CB05( 27)  &
          / 'tricyclene      ', 23            /
      DATA  MECH_NAM_CB05( 27)     , MECH_MAP_CB05( 27)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 27)  &
          / 9.00          /

      DATA  SPMH_NAM_CB05( 28)     , SPMH_MAP_CB05( 28)  &
          / 'tricyclene      ', 23            /
      DATA  MECH_NAM_CB05( 28)     , MECH_MAP_CB05( 28)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 28)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 29)     , SPMH_MAP_CB05( 29)  &
          / 'estragole       ', 24            /
      DATA  MECH_NAM_CB05( 29)     , MECH_MAP_CB05( 29)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 29)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 30)     , SPMH_MAP_CB05( 30)  &
          / 'camphor         ', 25            /
      DATA  MECH_NAM_CB05( 30)     , MECH_MAP_CB05( 30)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 30)  &
          / 8.00          /

      DATA  SPMH_NAM_CB05( 31)     , SPMH_MAP_CB05( 31)  &
          / 'camphor         ', 25            /
      DATA  MECH_NAM_CB05( 31)     , MECH_MAP_CB05( 31)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 31)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 32)     , SPMH_MAP_CB05( 32)  &
          / 'fenchone        ', 26            /
      DATA  MECH_NAM_CB05( 32)     , MECH_MAP_CB05( 32)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 32)  &
          / 8.00          /

      DATA  SPMH_NAM_CB05( 33)     , SPMH_MAP_CB05( 33)  &
          / 'fenchone        ', 26            /
      DATA  MECH_NAM_CB05( 33)     , MECH_MAP_CB05( 33)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 33)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05( 34)     , SPMH_MAP_CB05( 34)  &
          / 'piperitone      ', 27            /
      DATA  MECH_NAM_CB05( 34)     , MECH_MAP_CB05( 34)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 34)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 35)     , SPMH_MAP_CB05( 35)  &
          / 'thujone_a       ', 28            /
      DATA  MECH_NAM_CB05( 35)     , MECH_MAP_CB05( 35)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 35)  &
          / 9.00          /

      DATA  SPMH_NAM_CB05( 36)     , SPMH_MAP_CB05( 36)  &
          / 'thujone_a       ', 28            /
      DATA  MECH_NAM_CB05( 36)     , MECH_MAP_CB05( 36)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 36)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 37)     , SPMH_MAP_CB05( 37)  &
          / 'thujone_b       ', 29            /
      DATA  MECH_NAM_CB05( 37)     , MECH_MAP_CB05( 37)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 37)  &
          / 9.00          /

      DATA  SPMH_NAM_CB05( 38)     , SPMH_MAP_CB05( 38)  &
          / 'thujone_b       ', 29            /
      DATA  MECH_NAM_CB05( 38)     , MECH_MAP_CB05( 38)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 38)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 39)     , SPMH_MAP_CB05( 39)  &
          / 'cineole_1_8     ', 30            /
      DATA  MECH_NAM_CB05( 39)     , MECH_MAP_CB05( 39)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 39)  &
          / 9.00          /

      DATA  SPMH_NAM_CB05( 40)     , SPMH_MAP_CB05( 40)  &
          / 'cineole_1_8     ', 30            /
      DATA  MECH_NAM_CB05( 40)     , MECH_MAP_CB05( 40)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 40)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 41)     , SPMH_MAP_CB05( 41)  &
          / 'borneol         ', 31            /
      DATA  MECH_NAM_CB05( 41)     , MECH_MAP_CB05( 41)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 41)  &
          / 8.00          /

      DATA  SPMH_NAM_CB05( 42)     , SPMH_MAP_CB05( 42)  &
          / 'borneol         ', 31            /
      DATA  MECH_NAM_CB05( 42)     , MECH_MAP_CB05( 42)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 42)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05( 43)     , SPMH_MAP_CB05( 43)  &
          / 'linalool        ', 32            /
      DATA  MECH_NAM_CB05( 43)     , MECH_MAP_CB05( 43)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 43)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 44)     , SPMH_MAP_CB05( 44)  &
          / 'terpineol_4     ', 33            /
      DATA  MECH_NAM_CB05( 44)     , MECH_MAP_CB05( 44)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 44)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 45)     , SPMH_MAP_CB05( 45)  &
          / 'terpineol_a     ', 34            /
      DATA  MECH_NAM_CB05( 45)     , MECH_MAP_CB05( 45)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 45)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 46)     , SPMH_MAP_CB05( 46)  &
          / 'linalool_OXD_c  ', 35            /
      DATA  MECH_NAM_CB05( 46)     , MECH_MAP_CB05( 46)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 46)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 47)     , SPMH_MAP_CB05( 47)  &
          / 'linalool_OXD_t  ', 36            /
      DATA  MECH_NAM_CB05( 47)     , MECH_MAP_CB05( 47)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 47)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 48)     , SPMH_MAP_CB05( 48)  &
          / 'ionone_b        ', 37            /
      DATA  MECH_NAM_CB05( 48)     , MECH_MAP_CB05( 48)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05( 48)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 49)     , SPMH_MAP_CB05( 49)  &
          / 'ionone_b        ', 37            /
      DATA  MECH_NAM_CB05( 49)     , MECH_MAP_CB05( 49)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 49)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05( 50)     , SPMH_MAP_CB05( 50)  &
          / 'bornyl_ACT      ', 38            /
      DATA  MECH_NAM_CB05( 50)     , MECH_MAP_CB05( 50)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 50)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05( 51)     , SPMH_MAP_CB05( 51)  &
          / 'bornyl_ACT      ', 38            /
      DATA  MECH_NAM_CB05( 51)     , MECH_MAP_CB05( 51)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 51)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05( 52)     , SPMH_MAP_CB05( 52)  &
          / 'farnescene_a    ', 39            /
      DATA  MECH_NAM_CB05( 52)     , MECH_MAP_CB05( 52)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 52)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 53)     , SPMH_MAP_CB05( 53)  &
          / 'caryophyllene_b ', 40            /
      DATA  MECH_NAM_CB05( 53)     , MECH_MAP_CB05( 53)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 53)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 54)     , SPMH_MAP_CB05( 54)  &
          / 'acoradiene      ', 41            /
      DATA  MECH_NAM_CB05( 54)     , MECH_MAP_CB05( 54)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 54)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 55)     , SPMH_MAP_CB05( 55)  &
          / 'aromadendrene   ', 42            /
      DATA  MECH_NAM_CB05( 55)     , MECH_MAP_CB05( 55)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 55)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 56)     , SPMH_MAP_CB05( 56)  &
          / 'bergamotene_a   ', 43            /
      DATA  MECH_NAM_CB05( 56)     , MECH_MAP_CB05( 56)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 56)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 57)     , SPMH_MAP_CB05( 57)  &
          / 'bergamotene_b   ', 44            /
      DATA  MECH_NAM_CB05( 57)     , MECH_MAP_CB05( 57)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 57)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 58)     , SPMH_MAP_CB05( 58)  &
          / 'bisabolene_a    ', 45            /
      DATA  MECH_NAM_CB05( 58)     , MECH_MAP_CB05( 58)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 58)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 59)     , SPMH_MAP_CB05( 59)  &
          / 'bisabolene_b    ', 46            /
      DATA  MECH_NAM_CB05( 59)     , MECH_MAP_CB05( 59)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 59)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 60)     , SPMH_MAP_CB05( 60)  &
          / 'bourbonene_b    ', 47            /
      DATA  MECH_NAM_CB05( 60)     , MECH_MAP_CB05( 60)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 60)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 61)     , SPMH_MAP_CB05( 61)  &
          / 'cadinene_d      ', 48            /
      DATA  MECH_NAM_CB05( 61)     , MECH_MAP_CB05( 61)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 61)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 62)     , SPMH_MAP_CB05( 62)  &
          / 'cadinene_g      ', 49            /
      DATA  MECH_NAM_CB05( 62)     , MECH_MAP_CB05( 62)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 62)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 63)     , SPMH_MAP_CB05( 63)  &
          / 'cedrene_a       ', 50            /
      DATA  MECH_NAM_CB05( 63)     , MECH_MAP_CB05( 63)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 63)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 64)     , SPMH_MAP_CB05( 64)  &
          / 'copaene_a       ', 51            /
      DATA  MECH_NAM_CB05( 64)     , MECH_MAP_CB05( 64)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 64)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 65)     , SPMH_MAP_CB05( 65)  &
          / 'cubebene_a      ', 52            /
      DATA  MECH_NAM_CB05( 65)     , MECH_MAP_CB05( 65)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 65)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 66)     , SPMH_MAP_CB05( 66)  &
          / 'cubebene_b      ', 53            /
      DATA  MECH_NAM_CB05( 66)     , MECH_MAP_CB05( 66)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 66)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 67)     , SPMH_MAP_CB05( 67)  &
          / 'elemene_b       ', 54            /
      DATA  MECH_NAM_CB05( 67)     , MECH_MAP_CB05( 67)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 67)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 68)     , SPMH_MAP_CB05( 68)  &
          / 'farnescene_b    ', 55            /
      DATA  MECH_NAM_CB05( 68)     , MECH_MAP_CB05( 68)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 68)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 69)     , SPMH_MAP_CB05( 69)  &
          / 'germacrene_B    ', 56            /
      DATA  MECH_NAM_CB05( 69)     , MECH_MAP_CB05( 69)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 69)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 70)     , SPMH_MAP_CB05( 70)  &
          / 'germacrene_D    ', 57            /
      DATA  MECH_NAM_CB05( 70)     , MECH_MAP_CB05( 70)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 70)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 71)     , SPMH_MAP_CB05( 71)  &
          / 'gurjunene_b     ', 58            /
      DATA  MECH_NAM_CB05( 71)     , MECH_MAP_CB05( 71)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 71)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 72)     , SPMH_MAP_CB05( 72)  &
          / 'humulene_a      ', 59            /
      DATA  MECH_NAM_CB05( 72)     , MECH_MAP_CB05( 72)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 72)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 73)     , SPMH_MAP_CB05( 73)  &
          / 'humulene_g      ', 60            /
      DATA  MECH_NAM_CB05( 73)     , MECH_MAP_CB05( 73)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 73)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 74)     , SPMH_MAP_CB05( 74)  &
          / 'isolongifolene  ', 61            /
      DATA  MECH_NAM_CB05( 74)     , MECH_MAP_CB05( 74)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 74)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 75)     , SPMH_MAP_CB05( 75)  &
          / 'longifolene     ', 62            /
      DATA  MECH_NAM_CB05( 75)     , MECH_MAP_CB05( 75)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 75)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 76)     , SPMH_MAP_CB05( 76)  &
          / 'longipinene     ', 63            /
      DATA  MECH_NAM_CB05( 76)     , MECH_MAP_CB05( 76)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 76)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 77)     , SPMH_MAP_CB05( 77)  &
          / 'muurolene_a     ', 64            /
      DATA  MECH_NAM_CB05( 77)     , MECH_MAP_CB05( 77)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 77)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 78)     , SPMH_MAP_CB05( 78)  &
          / 'muurolene_g     ', 65            /
      DATA  MECH_NAM_CB05( 78)     , MECH_MAP_CB05( 78)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 78)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 79)     , SPMH_MAP_CB05( 79)  &
          / 'selinene_b      ', 66            /
      DATA  MECH_NAM_CB05( 79)     , MECH_MAP_CB05( 79)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 79)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 80)     , SPMH_MAP_CB05( 80)  &
          / 'selinene_d      ', 67            /
      DATA  MECH_NAM_CB05( 80)     , MECH_MAP_CB05( 80)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 80)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 81)     , SPMH_MAP_CB05( 81)  &
          / 'nerolidol_c     ', 68            /
      DATA  MECH_NAM_CB05( 81)     , MECH_MAP_CB05( 81)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 81)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 82)     , SPMH_MAP_CB05( 82)  &
          / 'nerolidol_t     ', 69            /
      DATA  MECH_NAM_CB05( 82)     , MECH_MAP_CB05( 82)  &
          / 'SESQ            ', 20            /
      DATA  CONV_FAC_CB05( 82)  &
          / 15.00         /

      DATA  SPMH_NAM_CB05( 83)     , SPMH_MAP_CB05( 83)  &
          / 'cedrol          ', 70            /
      DATA  MECH_NAM_CB05( 83)     , MECH_MAP_CB05( 83)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 83)  &
          / 13.00         /

      DATA  SPMH_NAM_CB05( 84)     , SPMH_MAP_CB05( 84)  &
          / 'cedrol          ', 70            /
      DATA  MECH_NAM_CB05( 84)     , MECH_MAP_CB05( 84)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 84)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05( 85)     , SPMH_MAP_CB05( 85)  &
          / 'MBO_2m3e2ol     ', 71            /
      DATA  MECH_NAM_CB05( 85)     , MECH_MAP_CB05( 85)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05( 85)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 86)     , SPMH_MAP_CB05( 86)  &
          / 'MBO_2m3e2ol     ', 71            /
      DATA  MECH_NAM_CB05( 86)     , MECH_MAP_CB05( 86)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 86)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05( 87)     , SPMH_MAP_CB05( 87)  &
          / 'methanol        ', 72            /
      DATA  MECH_NAM_CB05( 87)     , MECH_MAP_CB05( 87)  &
          / 'MEOH            ', 7             /
      DATA  CONV_FAC_CB05( 87)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 88)     , SPMH_MAP_CB05( 88)  &
          / 'acetone         ', 73            /
      DATA  MECH_NAM_CB05( 88)     , MECH_MAP_CB05( 88)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 88)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05( 89)     , SPMH_MAP_CB05( 89)  &
          / 'methane         ', 74            /
      DATA  MECH_NAM_CB05( 89)     , MECH_MAP_CB05( 89)  &
          / 'CH4             ', 8             /
      DATA  CONV_FAC_CB05( 89)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 90)     , SPMH_MAP_CB05( 90)  &
          / 'ammonia         ', 75            /
      DATA  MECH_NAM_CB05( 90)     , MECH_MAP_CB05( 90)  &
          / 'NH3             ', 9             /
      DATA  CONV_FAC_CB05( 90)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 91)     , SPMH_MAP_CB05( 91)  &
          / 'nitric_OXD      ', 77            /
      DATA  MECH_NAM_CB05( 91)     , MECH_MAP_CB05( 91)  &
          / 'NO              ', 10            /
      DATA  CONV_FAC_CB05( 91)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 92)     , SPMH_MAP_CB05( 92)  &
          / 'acetaldehyde    ', 78            /
      DATA  MECH_NAM_CB05( 92)     , MECH_MAP_CB05( 92)  &
          / 'ALD2            ', 11            /
      DATA  CONV_FAC_CB05( 92)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 93)     , SPMH_MAP_CB05( 93)  &
          / 'ethanol         ', 79            /
      DATA  MECH_NAM_CB05( 93)     , MECH_MAP_CB05( 93)  &
          / 'ETOH            ', 12            /
      DATA  CONV_FAC_CB05( 93)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 94)     , SPMH_MAP_CB05( 94)  &
          / 'formic_acid     ', 80            /
      DATA  MECH_NAM_CB05( 94)     , MECH_MAP_CB05( 94)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 94)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 95)     , SPMH_MAP_CB05( 95)  &
          / 'formaldehyde    ', 81            /
      DATA  MECH_NAM_CB05( 95)     , MECH_MAP_CB05( 95)  &
          / 'FORM            ', 13            /
      DATA  CONV_FAC_CB05( 95)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 96)     , SPMH_MAP_CB05( 96)  &
          / 'acetic_acid     ', 82            /
      DATA  MECH_NAM_CB05( 96)     , MECH_MAP_CB05( 96)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 96)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 97)     , SPMH_MAP_CB05( 97)  &
          / 'acetic_acid     ', 82            /
      DATA  MECH_NAM_CB05( 97)     , MECH_MAP_CB05( 97)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05( 97)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 98)     , SPMH_MAP_CB05( 98)  &
          / 'MBO_3m2e1ol     ', 83            /
      DATA  MECH_NAM_CB05( 98)     , MECH_MAP_CB05( 98)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05( 98)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05( 99)     , SPMH_MAP_CB05( 99)  &
          / 'MBO_3m2e1ol     ', 83            /
      DATA  MECH_NAM_CB05( 99)     , MECH_MAP_CB05( 99)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05( 99)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(100)     , SPMH_MAP_CB05(100)  &
          / 'MBO_3m3e1ol     ', 84            /
      DATA  MECH_NAM_CB05(100)     , MECH_MAP_CB05(100)  &
          / 'FORM            ', 13            /
      DATA  CONV_FAC_CB05(100)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(101)     , SPMH_MAP_CB05(101)  &
          / 'MBO_3m3e1ol     ', 84            /
      DATA  MECH_NAM_CB05(101)     , MECH_MAP_CB05(101)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(101)  &
          / 4.00          /

      DATA  SPMH_NAM_CB05(102)     , SPMH_MAP_CB05(102)  &
          / 'benzaldehyde    ', 85            /
      DATA  MECH_NAM_CB05(102)     , MECH_MAP_CB05(102)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(102)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(103)     , SPMH_MAP_CB05(103)  &
          / 'butanone_2      ', 86            /
      DATA  MECH_NAM_CB05(103)     , MECH_MAP_CB05(103)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(103)  &
          / 4.00          /

      DATA  SPMH_NAM_CB05(104)     , SPMH_MAP_CB05(104)  &
          / 'decanal         ', 87            /
      DATA  MECH_NAM_CB05(104)     , MECH_MAP_CB05(104)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(104)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(105)     , SPMH_MAP_CB05(105)  &
          / 'decanal         ', 87            /
      DATA  MECH_NAM_CB05(105)     , MECH_MAP_CB05(105)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(105)  &
          / 8.00          /

      DATA  SPMH_NAM_CB05(106)     , SPMH_MAP_CB05(106)  &
          / 'dodecene_1      ', 88            /
      DATA  MECH_NAM_CB05(106)     , MECH_MAP_CB05(106)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(106)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(107)     , SPMH_MAP_CB05(107)  &
          / 'dodecene_1      ', 88            /
      DATA  MECH_NAM_CB05(107)     , MECH_MAP_CB05(107)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(107)  &
          / 10.00         /

      DATA  SPMH_NAM_CB05(108)     , SPMH_MAP_CB05(108)  &
          / 'geranyl_acetone ', 89            /
      DATA  MECH_NAM_CB05(108)     , MECH_MAP_CB05(108)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(108)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(109)     , SPMH_MAP_CB05(109)  &
          / 'geranyl_acetone ', 89            /
      DATA  MECH_NAM_CB05(109)     , MECH_MAP_CB05(109)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(109)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(110)     , SPMH_MAP_CB05(110)  &
          / 'heptanal        ', 90            /
      DATA  MECH_NAM_CB05(110)     , MECH_MAP_CB05(110)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(110)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(111)     , SPMH_MAP_CB05(111)  &
          / 'heptanal        ', 90            /
      DATA  MECH_NAM_CB05(111)     , MECH_MAP_CB05(111)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(111)  &
          / 5.00          /

      DATA  SPMH_NAM_CB05(112)     , SPMH_MAP_CB05(112)  &
          / 'heptane         ', 91            /
      DATA  MECH_NAM_CB05(112)     , MECH_MAP_CB05(112)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(112)  &
          / 7.00          /

      DATA  SPMH_NAM_CB05(113)     , SPMH_MAP_CB05(113)  &
          / 'hexane          ', 92            /
      DATA  MECH_NAM_CB05(113)     , MECH_MAP_CB05(113)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(113)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05(114)     , SPMH_MAP_CB05(114)  &
          / 'met_benzoate    ', 93            /
      DATA  MECH_NAM_CB05(114)     , MECH_MAP_CB05(114)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(114)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(115)     , SPMH_MAP_CB05(115)  &
          / 'met_benzoate    ', 93            /
      DATA  MECH_NAM_CB05(115)     , MECH_MAP_CB05(115)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(115)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(116)     , SPMH_MAP_CB05(116)  &
          / 'met_heptenone   ', 94            /
      DATA  MECH_NAM_CB05(116)     , MECH_MAP_CB05(116)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(116)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05(117)     , SPMH_MAP_CB05(117)  &
          / 'met_heptenone   ', 94            /
      DATA  MECH_NAM_CB05(117)     , MECH_MAP_CB05(117)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(117)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(118)     , SPMH_MAP_CB05(118)  &
          / 'neryl_acetone   ', 95            /
      DATA  MECH_NAM_CB05(118)     , MECH_MAP_CB05(118)  &
          / 'IOLE            ', 16            /
      DATA  CONV_FAC_CB05(118)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(119)     , SPMH_MAP_CB05(119)  &
          / 'neryl_acetone   ', 95            /
      DATA  MECH_NAM_CB05(119)     , MECH_MAP_CB05(119)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(119)  &
          / 8.00          /

      DATA  SPMH_NAM_CB05(120)     , SPMH_MAP_CB05(120)  &
          / 'neryl_acetone   ', 95            /
      DATA  MECH_NAM_CB05(120)     , MECH_MAP_CB05(120)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(120)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(121)     , SPMH_MAP_CB05(121)  &
          / 'nonanal         ', 96            /
      DATA  MECH_NAM_CB05(121)     , MECH_MAP_CB05(121)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(121)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(122)     , SPMH_MAP_CB05(122)  &
          / 'nonanal         ', 96            /
      DATA  MECH_NAM_CB05(122)     , MECH_MAP_CB05(122)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(122)  &
          / 7.00          /

      DATA  SPMH_NAM_CB05(123)     , SPMH_MAP_CB05(123)  &
          / 'nonenal         ', 97            /
      DATA  MECH_NAM_CB05(123)     , MECH_MAP_CB05(123)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(123)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(124)     , SPMH_MAP_CB05(124)  &
          / 'nonenal         ', 97            /
      DATA  MECH_NAM_CB05(124)     , MECH_MAP_CB05(124)  &
          / 'IOLE            ', 16            /
      DATA  CONV_FAC_CB05(124)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(125)     , SPMH_MAP_CB05(125)  &
          / 'nonenal         ', 97            /
      DATA  MECH_NAM_CB05(125)     , MECH_MAP_CB05(125)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(125)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(126)     , SPMH_MAP_CB05(126)  &
          / 'octanal         ', 98            /
      DATA  MECH_NAM_CB05(126)     , MECH_MAP_CB05(126)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(126)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(127)     , SPMH_MAP_CB05(127)  &
          / 'octanal         ', 98            /
      DATA  MECH_NAM_CB05(127)     , MECH_MAP_CB05(127)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(127)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05(128)     , SPMH_MAP_CB05(128)  &
          / 'octanol         ', 99            /
      DATA  MECH_NAM_CB05(128)     , MECH_MAP_CB05(128)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(128)  &
          / 8.00          /

      DATA  SPMH_NAM_CB05(129)     , SPMH_MAP_CB05(129)  &
          / 'octenol_1e3ol   ', 100           /
      DATA  MECH_NAM_CB05(129)     , MECH_MAP_CB05(129)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(129)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05(130)     , SPMH_MAP_CB05(130)  &
          / 'octenol_1e3ol   ', 100           /
      DATA  MECH_NAM_CB05(130)     , MECH_MAP_CB05(130)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(130)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(131)     , SPMH_MAP_CB05(131)  &
          / 'oxopentanal     ', 101           /
      DATA  MECH_NAM_CB05(131)     , MECH_MAP_CB05(131)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(131)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(132)     , SPMH_MAP_CB05(132)  &
          / 'oxopentanal     ', 101           /
      DATA  MECH_NAM_CB05(132)     , MECH_MAP_CB05(132)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(132)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(133)     , SPMH_MAP_CB05(133)  &
          / 'pentane         ', 102           /
      DATA  MECH_NAM_CB05(133)     , MECH_MAP_CB05(133)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(133)  &
          / 5.00          /

      DATA  SPMH_NAM_CB05(134)     , SPMH_MAP_CB05(134)  &
          / 'phenyl_CCO      ', 103           /
      DATA  MECH_NAM_CB05(134)     , MECH_MAP_CB05(134)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(134)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(135)     , SPMH_MAP_CB05(135)  &
          / 'phenyl_CCO      ', 103           /
      DATA  MECH_NAM_CB05(135)     , MECH_MAP_CB05(135)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(135)  &
          / 0.50          /

      DATA  SPMH_NAM_CB05(136)     , SPMH_MAP_CB05(136)  &
          / 'pyruvic_acid    ', 104           /
      DATA  MECH_NAM_CB05(136)     , MECH_MAP_CB05(136)  &
          / 'FORM            ', 13            /
      DATA  CONV_FAC_CB05(136)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(137)     , SPMH_MAP_CB05(137)  &
          / 'pyruvic_acid    ', 104           /
      DATA  MECH_NAM_CB05(137)     , MECH_MAP_CB05(137)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(137)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(138)     , SPMH_MAP_CB05(138)  &
          / 'terpinyl_ACT_a  ', 105           /
      DATA  MECH_NAM_CB05(138)     , MECH_MAP_CB05(138)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(138)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(139)     , SPMH_MAP_CB05(139)  &
          / 'terpinyl_ACT_a  ', 105           /
      DATA  MECH_NAM_CB05(139)     , MECH_MAP_CB05(139)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(139)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(140)     , SPMH_MAP_CB05(140)  &
          / 'tetradecene_1   ', 106           /
      DATA  MECH_NAM_CB05(140)     , MECH_MAP_CB05(140)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(140)  &
          / 12.00         /

      DATA  SPMH_NAM_CB05(141)     , SPMH_MAP_CB05(141)  &
          / 'tetradecene_1   ', 106           /
      DATA  MECH_NAM_CB05(141)     , MECH_MAP_CB05(141)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(141)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(142)     , SPMH_MAP_CB05(142)  &
          / 'toluene         ', 107           /
      DATA  MECH_NAM_CB05(142)     , MECH_MAP_CB05(142)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(142)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(143)     , SPMH_MAP_CB05(143)  &
          / 'carbon_monoxide ', 108           /
      DATA  MECH_NAM_CB05(143)     , MECH_MAP_CB05(143)  &
          / 'CO              ', 17            /
      DATA  CONV_FAC_CB05(143)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(144)     , SPMH_MAP_CB05(144)  &
          / 'butene          ', 109           /
      DATA  MECH_NAM_CB05(144)     , MECH_MAP_CB05(144)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(144)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(145)     , SPMH_MAP_CB05(145)  &
          / 'butene          ', 109           /
      DATA  MECH_NAM_CB05(145)     , MECH_MAP_CB05(145)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(145)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(146)     , SPMH_MAP_CB05(146)  &
          / 'ethane          ', 110           /
      DATA  MECH_NAM_CB05(146)     , MECH_MAP_CB05(146)  &
          / 'ETHA            ', 18            /
      DATA  CONV_FAC_CB05(146)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(147)     , SPMH_MAP_CB05(147)  &
          / 'ethene          ', 111           /
      DATA  MECH_NAM_CB05(147)     , MECH_MAP_CB05(147)  &
          / 'ETH             ', 19            /
      DATA  CONV_FAC_CB05(147)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(148)     , SPMH_MAP_CB05(148)  &
          / 'propane         ', 113           /
      DATA  MECH_NAM_CB05(148)     , MECH_MAP_CB05(148)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(148)  &
          / 1.50          /

      DATA  SPMH_NAM_CB05(149)     , SPMH_MAP_CB05(149)  &
          / 'propane         ', 113           /
      DATA  MECH_NAM_CB05(149)     , MECH_MAP_CB05(149)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(149)  &
          / 1.50          /

      DATA  SPMH_NAM_CB05(150)     , SPMH_MAP_CB05(150)  &
          / 'propene         ', 114           /
      DATA  MECH_NAM_CB05(150)     , MECH_MAP_CB05(150)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(150)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(151)     , SPMH_MAP_CB05(151)  &
          / 'propene         ', 114           /
      DATA  MECH_NAM_CB05(151)     , MECH_MAP_CB05(151)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(151)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(152)     , SPMH_MAP_CB05(152)  &
          / 'diallyl_2s      ', 117           /
      DATA  MECH_NAM_CB05(152)     , MECH_MAP_CB05(152)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(152)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(153)     , SPMH_MAP_CB05(153)  &
          / 'diallyl_2s      ', 117           /
      DATA  MECH_NAM_CB05(153)     , MECH_MAP_CB05(153)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(153)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(154)     , SPMH_MAP_CB05(154)  &
          / '2met_2s         ', 118           /
      DATA  MECH_NAM_CB05(154)     , MECH_MAP_CB05(154)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(154)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(155)     , SPMH_MAP_CB05(155)  &
          / '2met_s          ', 119           /
      DATA  MECH_NAM_CB05(155)     , MECH_MAP_CB05(155)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(155)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(156)     , SPMH_MAP_CB05(156)  &
          / 'met_chloride    ', 120           /
      DATA  MECH_NAM_CB05(156)     , MECH_MAP_CB05(156)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(156)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(157)     , SPMH_MAP_CB05(157)  &
          / 'met_bromide     ', 121           /
      DATA  MECH_NAM_CB05(157)     , MECH_MAP_CB05(157)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(157)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(158)     , SPMH_MAP_CB05(158)  &
          / 'met_iodide      ', 122           /
      DATA  MECH_NAM_CB05(158)     , MECH_MAP_CB05(158)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(158)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(159)     , SPMH_MAP_CB05(159)  &
          / 'met_mercaptan   ', 124           /
      DATA  MECH_NAM_CB05(159)     , MECH_MAP_CB05(159)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(159)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(160)     , SPMH_MAP_CB05(160)  &
          / 'met_propenyl_2s ', 125           /
      DATA  MECH_NAM_CB05(160)     , MECH_MAP_CB05(160)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(160)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(161)     , SPMH_MAP_CB05(161)  &
          / 'met_propenyl_2s ', 125           /
      DATA  MECH_NAM_CB05(161)     , MECH_MAP_CB05(161)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(161)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(162)     , SPMH_MAP_CB05(162)  &
          / 'PPPP_2s         ', 126           /
      DATA  MECH_NAM_CB05(162)     , MECH_MAP_CB05(162)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(162)  &
          / 4.00          /

      DATA  SPMH_NAM_CB05(163)     , SPMH_MAP_CB05(163)  &
          / 'PPPP_2s         ', 126           /
      DATA  MECH_NAM_CB05(163)     , MECH_MAP_CB05(163)  &
          / 'OLE             ', 5             /
      DATA  CONV_FAC_CB05(163)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(164)     , SPMH_MAP_CB05(164)  &
          / '2met_nonatriene ', 127           /
      DATA  MECH_NAM_CB05(164)     , MECH_MAP_CB05(164)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(164)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(165)     , SPMH_MAP_CB05(165)  &
          / 'met_salicylate  ', 128           /
      DATA  MECH_NAM_CB05(165)     , MECH_MAP_CB05(165)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(165)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(166)     , SPMH_MAP_CB05(166)  &
          / 'met_salicylate  ', 128           /
      DATA  MECH_NAM_CB05(166)     , MECH_MAP_CB05(166)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(166)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(167)     , SPMH_MAP_CB05(167)  &
          / 'indole          ', 129           /
      DATA  MECH_NAM_CB05(167)     , MECH_MAP_CB05(167)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(167)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(168)     , SPMH_MAP_CB05(168)  &
          / 'indole          ', 129           /
      DATA  MECH_NAM_CB05(168)     , MECH_MAP_CB05(168)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(168)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(169)     , SPMH_MAP_CB05(169)  &
          / 'jasmone         ', 130           /
      DATA  MECH_NAM_CB05(169)     , MECH_MAP_CB05(169)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(169)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(170)     , SPMH_MAP_CB05(170)  &
          / 'jasmone         ', 130           /
      DATA  MECH_NAM_CB05(170)     , MECH_MAP_CB05(170)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(170)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(171)     , SPMH_MAP_CB05(171)  &
          / 'met_jasmonate   ', 131           /
      DATA  MECH_NAM_CB05(171)     , MECH_MAP_CB05(171)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(171)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(172)     , SPMH_MAP_CB05(172)  &
          / 'met_jasmonate   ', 131           /
      DATA  MECH_NAM_CB05(172)     , MECH_MAP_CB05(172)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(172)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(173)     , SPMH_MAP_CB05(173)  &
          / '3met_3DCTT      ', 132           /
      DATA  MECH_NAM_CB05(173)     , MECH_MAP_CB05(173)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(173)  &
          / 16.00         /

      DATA  SPMH_NAM_CB05(174)     , SPMH_MAP_CB05(174)  &
          / 'hexanal         ', 133           /
      DATA  MECH_NAM_CB05(174)     , MECH_MAP_CB05(174)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(174)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(175)     , SPMH_MAP_CB05(175)  &
          / 'hexanal         ', 133           /
      DATA  MECH_NAM_CB05(175)     , MECH_MAP_CB05(175)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(175)  &
          / 4.00          /

      DATA  SPMH_NAM_CB05(176)     , SPMH_MAP_CB05(176)  &
          / 'hexanol_1       ', 134           /
      DATA  MECH_NAM_CB05(176)     , MECH_MAP_CB05(176)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(176)  &
          / 6.00          /

      DATA  SPMH_NAM_CB05(177)     , SPMH_MAP_CB05(177)  &
          / 'hexenal_c3      ', 135           /
      DATA  MECH_NAM_CB05(177)     , MECH_MAP_CB05(177)  &
          / 'IOLE            ', 16            /
      DATA  CONV_FAC_CB05(177)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(178)     , SPMH_MAP_CB05(178)  &
          / 'hexenal_c3      ', 135           /
      DATA  MECH_NAM_CB05(178)     , MECH_MAP_CB05(178)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(178)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(179)     , SPMH_MAP_CB05(179)  &
          / 'hexenal_t2      ', 136           /
      DATA  MECH_NAM_CB05(179)     , MECH_MAP_CB05(179)  &
          / 'IOLE            ', 16            /
      DATA  CONV_FAC_CB05(179)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(180)     , SPMH_MAP_CB05(180)  &
          / 'hexenal_t2      ', 136           /
      DATA  MECH_NAM_CB05(180)     , MECH_MAP_CB05(180)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(180)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(181)     , SPMH_MAP_CB05(181)  &
          / 'hexenol_c3      ', 137           /
      DATA  MECH_NAM_CB05(181)     , MECH_MAP_CB05(181)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(181)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(182)     , SPMH_MAP_CB05(182)  &
          / 'hexenol_c3      ', 137           /
      DATA  MECH_NAM_CB05(182)     , MECH_MAP_CB05(182)  &
          / 'IOLE            ', 16            /
      DATA  CONV_FAC_CB05(182)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(183)     , SPMH_MAP_CB05(183)  &
          / 'hexenyl_ACT_c3  ', 138           /
      DATA  MECH_NAM_CB05(183)     , MECH_MAP_CB05(183)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(183)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(184)     , SPMH_MAP_CB05(184)  &
          / 'hexenyl_ACT_c3  ', 138           /
      DATA  MECH_NAM_CB05(184)     , MECH_MAP_CB05(184)  &
          / 'IOLE            ', 16            /
      DATA  CONV_FAC_CB05(184)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(185)     , SPMH_MAP_CB05(185)  &
          / 'hexenyl_ACT_c3  ', 138           /
      DATA  MECH_NAM_CB05(185)     , MECH_MAP_CB05(185)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(185)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(186)     , SPMH_MAP_CB05(186)  &
          / 'homosalate      ', 139           /
      DATA  MECH_NAM_CB05(186)     , MECH_MAP_CB05(186)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(186)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(187)     , SPMH_MAP_CB05(187)  &
          / 'homosalate      ', 139           /
      DATA  MECH_NAM_CB05(187)     , MECH_MAP_CB05(187)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(187)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(188)     , SPMH_MAP_CB05(188)  &
          / 'Ehsalate        ', 140           /
      DATA  MECH_NAM_CB05(188)     , MECH_MAP_CB05(188)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(188)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(189)     , SPMH_MAP_CB05(189)  &
          / 'Ehsalate        ', 140           /
      DATA  MECH_NAM_CB05(189)     , MECH_MAP_CB05(189)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(189)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(190)     , SPMH_MAP_CB05(190)  &
          / 'pentanal        ', 141           /
      DATA  MECH_NAM_CB05(190)     , MECH_MAP_CB05(190)  &
          / 'ALDX            ', 14            /
      DATA  CONV_FAC_CB05(190)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(191)     , SPMH_MAP_CB05(191)  &
          / 'pentanal        ', 141           /
      DATA  MECH_NAM_CB05(191)     , MECH_MAP_CB05(191)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(191)  &
          / 3.00          /

      DATA  SPMH_NAM_CB05(192)     , SPMH_MAP_CB05(192)  &
          / 'heptanone       ', 142           /
      DATA  MECH_NAM_CB05(192)     , MECH_MAP_CB05(192)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(192)  &
          / 7.00          /

      DATA  SPMH_NAM_CB05(193)     , SPMH_MAP_CB05(193)  &
          / 'anisole         ', 143           /
      DATA  MECH_NAM_CB05(193)     , MECH_MAP_CB05(193)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(193)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(194)     , SPMH_MAP_CB05(194)  &
          / 'verbenene       ', 144           /
      DATA  MECH_NAM_CB05(194)     , MECH_MAP_CB05(194)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(194)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(195)     , SPMH_MAP_CB05(195)  &
          / 'benzyl-acetate  ', 145           /
      DATA  MECH_NAM_CB05(195)     , MECH_MAP_CB05(195)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(195)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(196)     , SPMH_MAP_CB05(196)  &
          / 'benzyl-acetate  ', 145           /
      DATA  MECH_NAM_CB05(196)     , MECH_MAP_CB05(196)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(196)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(197)     , SPMH_MAP_CB05(197)  &
          / 'benzyl-acetate  ', 145           /
      DATA  MECH_NAM_CB05(197)     , MECH_MAP_CB05(197)  &
          / 'NR              ', 6             /
      DATA  CONV_FAC_CB05(197)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(198)     , SPMH_MAP_CB05(198)  &
          / 'myrtenal        ', 146           /
      DATA  MECH_NAM_CB05(198)     , MECH_MAP_CB05(198)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(198)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(199)     , SPMH_MAP_CB05(199)  &
          / 'benzyl-alcohol  ', 147           /
      DATA  MECH_NAM_CB05(199)     , MECH_MAP_CB05(199)  &
          / 'TOL             ', 15            /
      DATA  CONV_FAC_CB05(199)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(200)     , SPMH_MAP_CB05(200)  &
          / 'meta-cymenene   ', 148           /
      DATA  MECH_NAM_CB05(200)     , MECH_MAP_CB05(200)  &
          / 'XYL             ', 4             /
      DATA  CONV_FAC_CB05(200)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(201)     , SPMH_MAP_CB05(201)  &
          / 'meta-cymenene   ', 148           /
      DATA  MECH_NAM_CB05(201)     , MECH_MAP_CB05(201)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(201)  &
          / 2.00          /

      DATA  SPMH_NAM_CB05(202)     , SPMH_MAP_CB05(202)  &
          / 'ipsenol         ', 149           /
      DATA  MECH_NAM_CB05(202)     , MECH_MAP_CB05(202)  &
          / 'TERP            ', 2             /
      DATA  CONV_FAC_CB05(202)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(203)     , SPMH_MAP_CB05(203)  &
          / 'Napthalene      ', 150           /
      DATA  MECH_NAM_CB05(203)     , MECH_MAP_CB05(203)  &
          / 'XYL             ', 4             /
      DATA  CONV_FAC_CB05(203)  &
          / 1.00          /

      DATA  SPMH_NAM_CB05(204)     , SPMH_MAP_CB05(204)  &
          / 'Napthalene      ', 150           /
      DATA  MECH_NAM_CB05(204)     , MECH_MAP_CB05(204)  &
          / 'PAR             ', 3             /
      DATA  CONV_FAC_CB05(204)  &
          / 2.00          /

!=======================================================================
!  SPC_SAPRC99_SAPRC99.EXT
!  This include file contains SAPRC99 species and their MW.
!
!
!  Mechanism Name: SAPRC99
!  MEGAN v2.10
!  INPUT version x.x
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  bkoo         04/13/07 - Modified for new MECHANISM scheme
!=======================================================================

      CHARACTER*16   SPC_SAPRC99MECH
      PARAMETER     (SPC_SAPRC99MECH = 'SAPRC99         ')

      INTEGER        N_SAPRC99_SPC
      PARAMETER     (N_SAPRC99_SPC = 28)

      CHARACTER*16, target :: MECH_SPC_SAPRC99( N_SAPRC99_SPC ) ! Mechanism species name
      REAL                    MECH_MWT_SAPRC99( N_SAPRC99_SPC ) ! Mechanism species molecular weight

! Note conversion between 134 species and SAPRC99 is done by 1:1 mole

      DATA  MECH_SPC_SAPRC99(  1), MECH_MWT_SAPRC99(  1) / 'ISOPRENE        ', 68.00   /
      DATA  MECH_SPC_SAPRC99(  2), MECH_MWT_SAPRC99(  2) / 'TRP1            ', 136.00  /
      DATA  MECH_SPC_SAPRC99(  3), MECH_MWT_SAPRC99(  3) / 'MEOH            ', 32.00   /
      DATA  MECH_SPC_SAPRC99(  4), MECH_MWT_SAPRC99(  4) / 'ACET            ', 58.00   /
      DATA  MECH_SPC_SAPRC99(  5), MECH_MWT_SAPRC99(  5) / 'CH4             ', 16.00   /
      DATA  MECH_SPC_SAPRC99(  6), MECH_MWT_SAPRC99(  6) / 'NO              ', 30.00   /
      DATA  MECH_SPC_SAPRC99(  7), MECH_MWT_SAPRC99(  7) / 'NO2             ', 44.01   /
      DATA  MECH_SPC_SAPRC99(  8), MECH_MWT_SAPRC99(  8) / 'NH3             ', 17.00   /
      DATA  MECH_SPC_SAPRC99(  9), MECH_MWT_SAPRC99(  9) / 'CCHO            ', 44.00   /
      DATA  MECH_SPC_SAPRC99( 10), MECH_MWT_SAPRC99( 10) / 'HCOOH           ', 46.00   /
      DATA  MECH_SPC_SAPRC99( 11), MECH_MWT_SAPRC99( 11) / 'HCHO            ', 30.00   /
      DATA  MECH_SPC_SAPRC99( 12), MECH_MWT_SAPRC99( 12) / 'CCO_OH          ', 60.00   /
      DATA  MECH_SPC_SAPRC99( 13), MECH_MWT_SAPRC99( 13) / 'BALD            ', 106.00  /
      DATA  MECH_SPC_SAPRC99( 14), MECH_MWT_SAPRC99( 14) / 'MEK             ', 72.00   /
      DATA  MECH_SPC_SAPRC99( 15), MECH_MWT_SAPRC99( 15) / 'RCO_OH          ', 74.00   /
      DATA  MECH_SPC_SAPRC99( 16), MECH_MWT_SAPRC99( 16) / 'CO              ', 28.00   /
      DATA  MECH_SPC_SAPRC99( 17), MECH_MWT_SAPRC99( 17) / 'ETHENE          ', 28.00   /
      DATA  MECH_SPC_SAPRC99( 18), MECH_MWT_SAPRC99( 18) / 'ALK1            ', 30.10   /
      DATA  MECH_SPC_SAPRC99( 19), MECH_MWT_SAPRC99( 19) / 'ALK2            ', 36.70   /
      DATA  MECH_SPC_SAPRC99( 20), MECH_MWT_SAPRC99( 20) / 'ALK3            ', 58.60   /
      DATA  MECH_SPC_SAPRC99( 21), MECH_MWT_SAPRC99( 21) / 'ALK4            ', 77.60   /
      DATA  MECH_SPC_SAPRC99( 22), MECH_MWT_SAPRC99( 22) / 'ALK5            ', 118.90  /
      DATA  MECH_SPC_SAPRC99( 23), MECH_MWT_SAPRC99( 23) / 'ARO1            ', 98.60   /
      DATA  MECH_SPC_SAPRC99( 24), MECH_MWT_SAPRC99( 24) / 'ARO2            ', 118.70  /
      DATA  MECH_SPC_SAPRC99( 25), MECH_MWT_SAPRC99( 25) / 'OLE1            ', 72.30   /
      DATA  MECH_SPC_SAPRC99( 26), MECH_MWT_SAPRC99( 26) / 'OLE2            ', 75.80   /
      DATA  MECH_SPC_SAPRC99( 27), MECH_MWT_SAPRC99( 27) / 'RCHO            ', 58.00   /
      DATA  MECH_SPC_SAPRC99( 28), MECH_MWT_SAPRC99( 28) / 'NONR            ', 1.00    /

!=======================================================================
!  MAP_SAPRC99_CV2SAPRC99.EXT
!  This include file contains conversion table for 150 speciated species
!  to SAPRC99 species
!
!
!  MEGAN v2.10
!  INPUT version x.x
!
!  History:
!  Who          When       What
!  ---------------------------------------------------------------------
!  Tan          12/02/06 - Creates this file
!  bkoo         04/13/07 - Modified for new MECHANISM scheme
!  Tan          07/18/11 - Updated for MEGANv2.10
!=======================================================================

      CHARACTER*16   MAP_SAPRC99MECH
      PARAMETER     (MAP_SAPRC99MECH = 'SAPRC99         ')

      INTEGER        N_SAPRC99
      PARAMETER     (N_SAPRC99 = 150)        ! Number of map species

      CHARACTER*16       SPMH_NAM_SAPRC99( N_SAPRC99 )   ! speciated species name
      INTEGER, target :: SPMH_MAP_SAPRC99( N_SAPRC99 )   ! speciated species name
                                                         ! mapped to SPC_SPCAT.EXT
      CHARACTER*16       MECH_NAM_SAPRC99( N_SAPRC99 )   ! mechanism species
      INTEGER, target :: MECH_MAP_SAPRC99( N_SAPRC99 )   ! mechanism species mapped
                                                         ! to SPC_SAPRC99.EXT
      REAL, target ::    CONV_FAC_SAPRC99( N_SAPRC99 )   ! conversion factor


      DATA  SPMH_NAM_SAPRC99(  1)     , SPMH_MAP_SAPRC99(  1)  &
          / 'isoprene        ', 1              /
      DATA  MECH_NAM_SAPRC99(  1)     , MECH_MAP_SAPRC99(  1)  &
          / 'ISOPRENE        ', 1              /
      DATA  CONV_FAC_SAPRC99(  1)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  2)     , SPMH_MAP_SAPRC99(  2)  &
          / 'myrcene         ', 2              /
      DATA  MECH_NAM_SAPRC99(  2)     , MECH_MAP_SAPRC99(  2)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  2)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  3)     , SPMH_MAP_SAPRC99(  3)  &
          / 'sabinene        ', 3              /
      DATA  MECH_NAM_SAPRC99(  3)     , MECH_MAP_SAPRC99(  3)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  3)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  4)     , SPMH_MAP_SAPRC99(  4)  &
          / 'limonene        ', 4              /
      DATA  MECH_NAM_SAPRC99(  4)     , MECH_MAP_SAPRC99(  4)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  4)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  5)     , SPMH_MAP_SAPRC99(  5)  &
          / 'carene_3        ', 5              /
      DATA  MECH_NAM_SAPRC99(  5)     , MECH_MAP_SAPRC99(  5)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  5)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  6)     , SPMH_MAP_SAPRC99(  6)  &
          / 'ocimene_t_b     ', 6              /
      DATA  MECH_NAM_SAPRC99(  6)     , MECH_MAP_SAPRC99(  6)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  6)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  7)     , SPMH_MAP_SAPRC99(  7)  &
          / 'pinene_b        ', 7              /
      DATA  MECH_NAM_SAPRC99(  7)     , MECH_MAP_SAPRC99(  7)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  7)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  8)     , SPMH_MAP_SAPRC99(  8)  &
          / 'pinene_a        ', 8              /
      DATA  MECH_NAM_SAPRC99(  8)     , MECH_MAP_SAPRC99(  8)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(  8)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(  9)     , SPMH_MAP_SAPRC99(  9)  &
          / '2met_styrene    ', 9              /
      DATA  MECH_NAM_SAPRC99(  9)     , MECH_MAP_SAPRC99(  9)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99(  9)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 10)     , SPMH_MAP_SAPRC99( 10)  &
          / 'cymene_p        ', 10             /
      DATA  MECH_NAM_SAPRC99( 10)     , MECH_MAP_SAPRC99( 10)  &
          / 'ARO2            ', 24             /
      DATA  CONV_FAC_SAPRC99( 10)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 11)     , SPMH_MAP_SAPRC99( 11)  &
          / 'cymene_o        ', 11             /
      DATA  MECH_NAM_SAPRC99( 11)     , MECH_MAP_SAPRC99( 11)  &
          / 'ARO2            ', 24             /
      DATA  CONV_FAC_SAPRC99( 11)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 12)     , SPMH_MAP_SAPRC99( 12)  &
          / 'phellandrene_a  ', 12             /
      DATA  MECH_NAM_SAPRC99( 12)     , MECH_MAP_SAPRC99( 12)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 12)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 13)     , SPMH_MAP_SAPRC99( 13)  &
          / 'thujene_a       ', 13             /
      DATA  MECH_NAM_SAPRC99( 13)     , MECH_MAP_SAPRC99( 13)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 13)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 14)     , SPMH_MAP_SAPRC99( 14)  &
          / 'terpinene_a     ', 14             /
      DATA  MECH_NAM_SAPRC99( 14)     , MECH_MAP_SAPRC99( 14)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 14)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 15)     , SPMH_MAP_SAPRC99( 15)  &
          / 'terpinene_g     ', 15             /
      DATA  MECH_NAM_SAPRC99( 15)     , MECH_MAP_SAPRC99( 15)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 15)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 16)     , SPMH_MAP_SAPRC99( 16)  &
          / 'terpinolene     ', 16             /
      DATA  MECH_NAM_SAPRC99( 16)     , MECH_MAP_SAPRC99( 16)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 16)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 17)     , SPMH_MAP_SAPRC99( 17)  &
          / 'phellandrene_b  ', 17             /
      DATA  MECH_NAM_SAPRC99( 17)     , MECH_MAP_SAPRC99( 17)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 17)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 18)     , SPMH_MAP_SAPRC99( 18)  &
          / 'camphene        ', 18             /
      DATA  MECH_NAM_SAPRC99( 18)     , MECH_MAP_SAPRC99( 18)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 18)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 19)     , SPMH_MAP_SAPRC99( 19)  &
          / 'bornene         ', 19             /
      DATA  MECH_NAM_SAPRC99( 19)     , MECH_MAP_SAPRC99( 19)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 19)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 20)     , SPMH_MAP_SAPRC99( 20)  &
          / 'fenchene_a      ', 20             /
      DATA  MECH_NAM_SAPRC99( 20)     , MECH_MAP_SAPRC99( 20)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 20)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 21)     , SPMH_MAP_SAPRC99( 21)  &
          / 'ocimene_al      ', 21             /
      DATA  MECH_NAM_SAPRC99( 21)     , MECH_MAP_SAPRC99( 21)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 21)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 22)     , SPMH_MAP_SAPRC99( 22)  &
          / 'ocimene_c_b     ', 22             /
      DATA  MECH_NAM_SAPRC99( 22)     , MECH_MAP_SAPRC99( 22)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 22)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 23)     , SPMH_MAP_SAPRC99( 23)  &
          / 'tricyclene      ', 23             /
      DATA  MECH_NAM_SAPRC99( 23)     , MECH_MAP_SAPRC99( 23)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 23)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 24)     , SPMH_MAP_SAPRC99( 24)  &
          / 'estragole       ', 24             /
      DATA  MECH_NAM_SAPRC99( 24)     , MECH_MAP_SAPRC99( 24)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 24)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 25)     , SPMH_MAP_SAPRC99( 25)  &
          / 'camphor         ', 25             /
      DATA  MECH_NAM_SAPRC99( 25)     , MECH_MAP_SAPRC99( 25)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 25)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 26)     , SPMH_MAP_SAPRC99( 26)  &
          / 'fenchone        ', 26             /
      DATA  MECH_NAM_SAPRC99( 26)     , MECH_MAP_SAPRC99( 26)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 26)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 27)     , SPMH_MAP_SAPRC99( 27)  &
          / 'piperitone      ', 27             /
      DATA  MECH_NAM_SAPRC99( 27)     , MECH_MAP_SAPRC99( 27)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 27)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 28)     , SPMH_MAP_SAPRC99( 28)  &
          / 'thujone_a       ', 28             /
      DATA  MECH_NAM_SAPRC99( 28)     , MECH_MAP_SAPRC99( 28)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 28)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 29)     , SPMH_MAP_SAPRC99( 29)  &
          / 'thujone_b       ', 29             /
      DATA  MECH_NAM_SAPRC99( 29)     , MECH_MAP_SAPRC99( 29)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 29)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 30)     , SPMH_MAP_SAPRC99( 30)  &
          / 'cineole_1_8     ', 30             /
      DATA  MECH_NAM_SAPRC99( 30)     , MECH_MAP_SAPRC99( 30)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 30)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 31)     , SPMH_MAP_SAPRC99( 31)  &
          / 'borneol         ', 31             /
      DATA  MECH_NAM_SAPRC99( 31)     , MECH_MAP_SAPRC99( 31)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 31)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 32)     , SPMH_MAP_SAPRC99( 32)  &
          / 'linalool        ', 32             /
      DATA  MECH_NAM_SAPRC99( 32)     , MECH_MAP_SAPRC99( 32)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 32)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 33)     , SPMH_MAP_SAPRC99( 33)  &
          / 'terpineol_4     ', 33             /
      DATA  MECH_NAM_SAPRC99( 33)     , MECH_MAP_SAPRC99( 33)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 33)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 34)     , SPMH_MAP_SAPRC99( 34)  &
          / 'terpineol_a     ', 34             /
      DATA  MECH_NAM_SAPRC99( 34)     , MECH_MAP_SAPRC99( 34)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 34)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 35)     , SPMH_MAP_SAPRC99( 35)  &
          / 'linalool_OXD_c  ', 35             /
      DATA  MECH_NAM_SAPRC99( 35)     , MECH_MAP_SAPRC99( 35)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 35)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 36)     , SPMH_MAP_SAPRC99( 36)  &
          / 'linalool_OXD_t  ', 36             /
      DATA  MECH_NAM_SAPRC99( 36)     , MECH_MAP_SAPRC99( 36)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 36)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 37)     , SPMH_MAP_SAPRC99( 37)  &
          / 'ionone_b        ', 37             /
      DATA  MECH_NAM_SAPRC99( 37)     , MECH_MAP_SAPRC99( 37)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 37)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 38)     , SPMH_MAP_SAPRC99( 38)  &
          / 'bornyl_ACT      ', 38             /
      DATA  MECH_NAM_SAPRC99( 38)     , MECH_MAP_SAPRC99( 38)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 38)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 39)     , SPMH_MAP_SAPRC99( 39)  &
          / 'farnescene_a    ', 39             /
      DATA  MECH_NAM_SAPRC99( 39)     , MECH_MAP_SAPRC99( 39)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 39)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 40)     , SPMH_MAP_SAPRC99( 40)  &
          / 'caryophyllene_b ', 40             /
      DATA  MECH_NAM_SAPRC99( 40)     , MECH_MAP_SAPRC99( 40)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 40)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 41)     , SPMH_MAP_SAPRC99( 41)  &
          / 'acoradiene      ', 41             /
      DATA  MECH_NAM_SAPRC99( 41)     , MECH_MAP_SAPRC99( 41)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 41)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 42)     , SPMH_MAP_SAPRC99( 42)  &
          / 'aromadendrene   ', 42             /
      DATA  MECH_NAM_SAPRC99( 42)     , MECH_MAP_SAPRC99( 42)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 42)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 43)     , SPMH_MAP_SAPRC99( 43)  &
          / 'bergamotene_a   ', 43             /
      DATA  MECH_NAM_SAPRC99( 43)     , MECH_MAP_SAPRC99( 43)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 43)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 44)     , SPMH_MAP_SAPRC99( 44)  &
          / 'bergamotene_b   ', 44             /
      DATA  MECH_NAM_SAPRC99( 44)     , MECH_MAP_SAPRC99( 44)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 44)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 45)     , SPMH_MAP_SAPRC99( 45)  &
          / 'bisabolene_a    ', 45             /
      DATA  MECH_NAM_SAPRC99( 45)     , MECH_MAP_SAPRC99( 45)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 45)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 46)     , SPMH_MAP_SAPRC99( 46)  &
          / 'bisabolene_b    ', 46             /
      DATA  MECH_NAM_SAPRC99( 46)     , MECH_MAP_SAPRC99( 46)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 46)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 47)     , SPMH_MAP_SAPRC99( 47)  &
          / 'bourbonene_b    ', 47             /
      DATA  MECH_NAM_SAPRC99( 47)     , MECH_MAP_SAPRC99( 47)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 47)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 48)     , SPMH_MAP_SAPRC99( 48)  &
          / 'cadinene_d      ', 48             /
      DATA  MECH_NAM_SAPRC99( 48)     , MECH_MAP_SAPRC99( 48)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 48)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 49)     , SPMH_MAP_SAPRC99( 49)  &
          / 'cadinene_g      ', 49             /
      DATA  MECH_NAM_SAPRC99( 49)     , MECH_MAP_SAPRC99( 49)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 49)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 50)     , SPMH_MAP_SAPRC99( 50)  &
          / 'cedrene_a       ', 50             /
      DATA  MECH_NAM_SAPRC99( 50)     , MECH_MAP_SAPRC99( 50)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 50)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 51)     , SPMH_MAP_SAPRC99( 51)  &
          / 'copaene_a       ', 51             /
      DATA  MECH_NAM_SAPRC99( 51)     , MECH_MAP_SAPRC99( 51)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 51)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 52)     , SPMH_MAP_SAPRC99( 52)  &
          / 'cubebene_a      ', 52             /
      DATA  MECH_NAM_SAPRC99( 52)     , MECH_MAP_SAPRC99( 52)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 52)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 53)     , SPMH_MAP_SAPRC99( 53)  &
          / 'cubebene_b      ', 53             /
      DATA  MECH_NAM_SAPRC99( 53)     , MECH_MAP_SAPRC99( 53)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 53)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 54)     , SPMH_MAP_SAPRC99( 54)  &
          / 'elemene_b       ', 54             /
      DATA  MECH_NAM_SAPRC99( 54)     , MECH_MAP_SAPRC99( 54)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 54)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 55)     , SPMH_MAP_SAPRC99( 55)  &
          / 'farnescene_b    ', 55             /
      DATA  MECH_NAM_SAPRC99( 55)     , MECH_MAP_SAPRC99( 55)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 55)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 56)     , SPMH_MAP_SAPRC99( 56)  &
          / 'germacrene_B    ', 56             /
      DATA  MECH_NAM_SAPRC99( 56)     , MECH_MAP_SAPRC99( 56)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 56)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 57)     , SPMH_MAP_SAPRC99( 57)  &
          / 'germacrene_D    ', 57             /
      DATA  MECH_NAM_SAPRC99( 57)     , MECH_MAP_SAPRC99( 57)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 57)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 58)     , SPMH_MAP_SAPRC99( 58)  &
          / 'gurjunene_b     ', 58             /
      DATA  MECH_NAM_SAPRC99( 58)     , MECH_MAP_SAPRC99( 58)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 58)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 59)     , SPMH_MAP_SAPRC99( 59)  &
          / 'humulene_a      ', 59             /
      DATA  MECH_NAM_SAPRC99( 59)     , MECH_MAP_SAPRC99( 59)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 59)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 60)     , SPMH_MAP_SAPRC99( 60)  &
          / 'humulene_g      ', 60             /
      DATA  MECH_NAM_SAPRC99( 60)     , MECH_MAP_SAPRC99( 60)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 60)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 61)     , SPMH_MAP_SAPRC99( 61)  &
          / 'isolongifolene  ', 61             /
      DATA  MECH_NAM_SAPRC99( 61)     , MECH_MAP_SAPRC99( 61)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 61)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 62)     , SPMH_MAP_SAPRC99( 62)  &
          / 'longifolene     ', 62             /
      DATA  MECH_NAM_SAPRC99( 62)     , MECH_MAP_SAPRC99( 62)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 62)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 63)     , SPMH_MAP_SAPRC99( 63)  &
          / 'longipinene     ', 63             /
      DATA  MECH_NAM_SAPRC99( 63)     , MECH_MAP_SAPRC99( 63)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 63)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 64)     , SPMH_MAP_SAPRC99( 64)  &
          / 'muurolene_a     ', 64             /
      DATA  MECH_NAM_SAPRC99( 64)     , MECH_MAP_SAPRC99( 64)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 64)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 65)     , SPMH_MAP_SAPRC99( 65)  &
          / 'muurolene_g     ', 65             /
      DATA  MECH_NAM_SAPRC99( 65)     , MECH_MAP_SAPRC99( 65)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 65)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 66)     , SPMH_MAP_SAPRC99( 66)  &
          / 'selinene_b      ', 66             /
      DATA  MECH_NAM_SAPRC99( 66)     , MECH_MAP_SAPRC99( 66)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 66)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 67)     , SPMH_MAP_SAPRC99( 67)  &
          / 'selinene_d      ', 67             /
      DATA  MECH_NAM_SAPRC99( 67)     , MECH_MAP_SAPRC99( 67)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 67)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 68)     , SPMH_MAP_SAPRC99( 68)  &
          / 'nerolidol_c     ', 68             /
      DATA  MECH_NAM_SAPRC99( 68)     , MECH_MAP_SAPRC99( 68)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 68)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 69)     , SPMH_MAP_SAPRC99( 69)  &
          / 'nerolidol_t     ', 69             /
      DATA  MECH_NAM_SAPRC99( 69)     , MECH_MAP_SAPRC99( 69)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 69)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 70)     , SPMH_MAP_SAPRC99( 70)  &
          / 'cedrol          ', 70             /
      DATA  MECH_NAM_SAPRC99( 70)     , MECH_MAP_SAPRC99( 70)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 70)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 71)     , SPMH_MAP_SAPRC99( 71)  &
          / 'MBO_2m3e2ol     ', 71             /
      DATA  MECH_NAM_SAPRC99( 71)     , MECH_MAP_SAPRC99( 71)  &
          / 'ISOPRENE        ', 1              /
      DATA  CONV_FAC_SAPRC99( 71)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 72)     , SPMH_MAP_SAPRC99( 72)  &
          / 'methanol        ', 72             /
      DATA  MECH_NAM_SAPRC99( 72)     , MECH_MAP_SAPRC99( 72)  &
          / 'MEOH            ', 3              /
      DATA  CONV_FAC_SAPRC99( 72)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 73)     , SPMH_MAP_SAPRC99( 73)  &
          / 'acetone         ', 73             /
      DATA  MECH_NAM_SAPRC99( 73)     , MECH_MAP_SAPRC99( 73)  &
          / 'ACET            ', 4              /
      DATA  CONV_FAC_SAPRC99( 73)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 74)     , SPMH_MAP_SAPRC99( 74)  &
          / 'methane         ', 74             /
      DATA  MECH_NAM_SAPRC99( 74)     , MECH_MAP_SAPRC99( 74)  &
          / 'CH4             ', 5              /
      DATA  CONV_FAC_SAPRC99( 74)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 75)     , SPMH_MAP_SAPRC99( 75)  &
          / 'ammonia         ', 75             /
      DATA  MECH_NAM_SAPRC99( 75)     , MECH_MAP_SAPRC99( 75)  &
          / 'NH3             ', 8              /
      DATA  CONV_FAC_SAPRC99( 75)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 76)     , SPMH_MAP_SAPRC99( 76)  &
          / 'nitrous_OXD     ', 76             /
      DATA  MECH_NAM_SAPRC99( 76)     , MECH_MAP_SAPRC99( 76)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99( 76)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 77)     , SPMH_MAP_SAPRC99( 77)  &
          / 'nitric_OXD      ', 77             /
      DATA  MECH_NAM_SAPRC99( 77)     , MECH_MAP_SAPRC99( 77)  &
          / 'NO              ', 6              /
      DATA  CONV_FAC_SAPRC99( 77)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 78)     , SPMH_MAP_SAPRC99( 78)  &
          / 'acetaldehyde    ', 78             /
      DATA  MECH_NAM_SAPRC99( 78)     , MECH_MAP_SAPRC99( 78)  &
          / 'CCHO            ', 9              /
      DATA  CONV_FAC_SAPRC99( 78)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 79)     , SPMH_MAP_SAPRC99( 79)  &
          / 'ethanol         ', 79             /
      DATA  MECH_NAM_SAPRC99( 79)     , MECH_MAP_SAPRC99( 79)  &
          / 'ALK3            ', 20             /
      DATA  CONV_FAC_SAPRC99( 79)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 80)     , SPMH_MAP_SAPRC99( 80)  &
          / 'formic_acid     ', 80             /
      DATA  MECH_NAM_SAPRC99( 80)     , MECH_MAP_SAPRC99( 80)  &
          / 'HCOOH           ', 10             /
      DATA  CONV_FAC_SAPRC99( 80)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 81)     , SPMH_MAP_SAPRC99( 81)  &
          / 'formaldehyde    ', 81             /
      DATA  MECH_NAM_SAPRC99( 81)     , MECH_MAP_SAPRC99( 81)  &
          / 'HCHO            ', 11             /
      DATA  CONV_FAC_SAPRC99( 81)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 82)     , SPMH_MAP_SAPRC99( 82)  &
          / 'acetic_acid     ', 82             /
      DATA  MECH_NAM_SAPRC99( 82)     , MECH_MAP_SAPRC99( 82)  &
          / 'CCO_OH          ', 12             /
      DATA  CONV_FAC_SAPRC99( 82)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 83)     , SPMH_MAP_SAPRC99( 83)  &
          / 'MBO_3m2e1ol     ', 83             /
      DATA  MECH_NAM_SAPRC99( 83)     , MECH_MAP_SAPRC99( 83)  &
          / 'ISOPRENE        ', 1              /
      DATA  CONV_FAC_SAPRC99( 83)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 84)     , SPMH_MAP_SAPRC99( 84)  &
          / 'MBO_3m3e1ol     ', 84             /
      DATA  MECH_NAM_SAPRC99( 84)     , MECH_MAP_SAPRC99( 84)  &
          / 'ISOPRENE        ', 1              /
      DATA  CONV_FAC_SAPRC99( 84)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 85)     , SPMH_MAP_SAPRC99( 85)  &
          / 'benzaldehyde    ', 85             /
      DATA  MECH_NAM_SAPRC99( 85)     , MECH_MAP_SAPRC99( 85)  &
          / 'BALD            ', 13             /
      DATA  CONV_FAC_SAPRC99( 85)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 86)     , SPMH_MAP_SAPRC99( 86)  &
          / 'butanone_2      ', 86             /
      DATA  MECH_NAM_SAPRC99( 86)     , MECH_MAP_SAPRC99( 86)  &
          / 'MEK             ', 14             /
      DATA  CONV_FAC_SAPRC99( 86)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 87)     , SPMH_MAP_SAPRC99( 87)  &
          / 'decanal         ', 87             /
      DATA  MECH_NAM_SAPRC99( 87)     , MECH_MAP_SAPRC99( 87)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99( 87)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 88)     , SPMH_MAP_SAPRC99( 88)  &
          / 'dodecene_1      ', 88             /
      DATA  MECH_NAM_SAPRC99( 88)     , MECH_MAP_SAPRC99( 88)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99( 88)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 89)     , SPMH_MAP_SAPRC99( 89)  &
          / 'geranyl_acetone ', 89             /
      DATA  MECH_NAM_SAPRC99( 89)     , MECH_MAP_SAPRC99( 89)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99( 89)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 90)     , SPMH_MAP_SAPRC99( 90)  &
          / 'heptanal        ', 90             /
      DATA  MECH_NAM_SAPRC99( 90)     , MECH_MAP_SAPRC99( 90)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99( 90)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 91)     , SPMH_MAP_SAPRC99( 91)  &
          / 'heptane         ', 91             /
      DATA  MECH_NAM_SAPRC99( 91)     , MECH_MAP_SAPRC99( 91)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 91)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 92)     , SPMH_MAP_SAPRC99( 92)  &
          / 'hexane          ', 92             /
      DATA  MECH_NAM_SAPRC99( 92)     , MECH_MAP_SAPRC99( 92)  &
          / 'ALK4            ', 21             /
      DATA  CONV_FAC_SAPRC99( 92)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 93)     , SPMH_MAP_SAPRC99( 93)  &
          / 'met_benzoate    ', 93             /
      DATA  MECH_NAM_SAPRC99( 93)     , MECH_MAP_SAPRC99( 93)  &
          / 'ARO1            ', 23             /
      DATA  CONV_FAC_SAPRC99( 93)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 94)     , SPMH_MAP_SAPRC99( 94)  &
          / 'met_heptenone   ', 94             /
      DATA  MECH_NAM_SAPRC99( 94)     , MECH_MAP_SAPRC99( 94)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99( 94)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 95)     , SPMH_MAP_SAPRC99( 95)  &
          / 'neryl_acetone   ', 95             /
      DATA  MECH_NAM_SAPRC99( 95)     , MECH_MAP_SAPRC99( 95)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99( 95)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 96)     , SPMH_MAP_SAPRC99( 96)  &
          / 'nonanal         ', 96             /
      DATA  MECH_NAM_SAPRC99( 96)     , MECH_MAP_SAPRC99( 96)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99( 96)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 97)     , SPMH_MAP_SAPRC99( 97)  &
          / 'nonenal         ', 97             /
      DATA  MECH_NAM_SAPRC99( 97)     , MECH_MAP_SAPRC99( 97)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99( 97)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 98)     , SPMH_MAP_SAPRC99( 98)  &
          / 'octanal         ', 98             /
      DATA  MECH_NAM_SAPRC99( 98)     , MECH_MAP_SAPRC99( 98)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99( 98)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99( 99)     , SPMH_MAP_SAPRC99( 99)  &
          / 'octanol         ', 99             /
      DATA  MECH_NAM_SAPRC99( 99)     , MECH_MAP_SAPRC99( 99)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99( 99)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(100)     , SPMH_MAP_SAPRC99(100)  &
          / 'octenol_1e3ol   ', 100            /
      DATA  MECH_NAM_SAPRC99(100)     , MECH_MAP_SAPRC99(100)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(100)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(101)     , SPMH_MAP_SAPRC99(101)  &
          / 'oxopentanal     ', 101            /
      DATA  MECH_NAM_SAPRC99(101)     , MECH_MAP_SAPRC99(101)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99(101)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(102)     , SPMH_MAP_SAPRC99(102)  &
          / 'pentane         ', 102            /
      DATA  MECH_NAM_SAPRC99(102)     , MECH_MAP_SAPRC99(102)  &
          / 'ALK4            ', 21             /
      DATA  CONV_FAC_SAPRC99(102)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(103)     , SPMH_MAP_SAPRC99(103)  &
          / 'phenyl_CCO      ', 103            /
      DATA  MECH_NAM_SAPRC99(103)     , MECH_MAP_SAPRC99(103)  &
          / 'ARO1            ', 23             /
      DATA  CONV_FAC_SAPRC99(103)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(104)     , SPMH_MAP_SAPRC99(104)  &
          / 'pyruvic_acid    ', 104            /
      DATA  MECH_NAM_SAPRC99(104)     , MECH_MAP_SAPRC99(104)  &
          / 'RCO_OH          ', 15             /
      DATA  CONV_FAC_SAPRC99(104)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(105)     , SPMH_MAP_SAPRC99(105)  &
          / 'terpinyl_ACT_a  ', 105            /
      DATA  MECH_NAM_SAPRC99(105)     , MECH_MAP_SAPRC99(105)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(105)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(106)     , SPMH_MAP_SAPRC99(106)  &
          / 'tetradecene_1   ', 106            /
      DATA  MECH_NAM_SAPRC99(106)     , MECH_MAP_SAPRC99(106)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(106)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(107)     , SPMH_MAP_SAPRC99(107)  &
          / 'toluene         ', 107            /
      DATA  MECH_NAM_SAPRC99(107)     , MECH_MAP_SAPRC99(107)  &
          / 'ARO1            ', 23             /
      DATA  CONV_FAC_SAPRC99(107)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(108)     , SPMH_MAP_SAPRC99(108)  &
          / 'carbon_monoxide ', 108            /
      DATA  MECH_NAM_SAPRC99(108)     , MECH_MAP_SAPRC99(108)  &
          / 'CO              ', 16             /
      DATA  CONV_FAC_SAPRC99(108)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(109)     , SPMH_MAP_SAPRC99(109)  &
          / 'butene          ', 109            /
      DATA  MECH_NAM_SAPRC99(109)     , MECH_MAP_SAPRC99(109)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(109)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(110)     , SPMH_MAP_SAPRC99(110)  &
          / 'ethane          ', 110            /
      DATA  MECH_NAM_SAPRC99(110)     , MECH_MAP_SAPRC99(110)  &
          / 'ALK1            ', 18             /
      DATA  CONV_FAC_SAPRC99(110)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(111)     , SPMH_MAP_SAPRC99(111)  &
          / 'ethene          ', 111            /
      DATA  MECH_NAM_SAPRC99(111)     , MECH_MAP_SAPRC99(111)  &
          / 'ETHENE          ', 17             /
      DATA  CONV_FAC_SAPRC99(111)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(112)     , SPMH_MAP_SAPRC99(112)  &
          / 'hydrogen_cyanide', 112            /
      DATA  MECH_NAM_SAPRC99(112)     , MECH_MAP_SAPRC99(112)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(112)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(113)     , SPMH_MAP_SAPRC99(113)  &
          / 'propane         ', 113            /
      DATA  MECH_NAM_SAPRC99(113)     , MECH_MAP_SAPRC99(113)  &
          / 'ALK2            ', 19             /
      DATA  CONV_FAC_SAPRC99(113)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(114)     , SPMH_MAP_SAPRC99(114)  &
          / 'propene         ', 114            /
      DATA  MECH_NAM_SAPRC99(114)     , MECH_MAP_SAPRC99(114)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(114)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(115)     , SPMH_MAP_SAPRC99(115)  &
          / 'carbon_2s       ', 115            /
      DATA  MECH_NAM_SAPRC99(115)     , MECH_MAP_SAPRC99(115)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(115)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(116)     , SPMH_MAP_SAPRC99(116)  &
          / 'carbonyl_s      ', 116            /
      DATA  MECH_NAM_SAPRC99(116)     , MECH_MAP_SAPRC99(116)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(116)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(117)     , SPMH_MAP_SAPRC99(117)  &
          / 'diallyl_2s      ', 117            /
      DATA  MECH_NAM_SAPRC99(117)     , MECH_MAP_SAPRC99(117)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(117)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(118)     , SPMH_MAP_SAPRC99(118)  &
          / '2met_2s         ', 118            /
      DATA  MECH_NAM_SAPRC99(118)     , MECH_MAP_SAPRC99(118)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99(118)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(119)     , SPMH_MAP_SAPRC99(119)  &
          / '2met_s          ', 119            /
      DATA  MECH_NAM_SAPRC99(119)     , MECH_MAP_SAPRC99(119)  &
          / 'ALK4            ', 21             /
      DATA  CONV_FAC_SAPRC99(119)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(120)     , SPMH_MAP_SAPRC99(120)  &
          / 'met_chloride    ', 120            /
      DATA  MECH_NAM_SAPRC99(120)     , MECH_MAP_SAPRC99(120)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(120)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(121)     , SPMH_MAP_SAPRC99(121)  &
          / 'met_bromide     ', 121            /
      DATA  MECH_NAM_SAPRC99(121)     , MECH_MAP_SAPRC99(121)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(121)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(122)     , SPMH_MAP_SAPRC99(122)  &
          / 'met_iodide      ', 122            /
      DATA  MECH_NAM_SAPRC99(122)     , MECH_MAP_SAPRC99(122)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(122)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(123)     , SPMH_MAP_SAPRC99(123)  &
          / 'hydrogen_s      ', 123            /
      DATA  MECH_NAM_SAPRC99(123)     , MECH_MAP_SAPRC99(123)  &
          / 'NONR            ', 28             /
      DATA  CONV_FAC_SAPRC99(123)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(124)     , SPMH_MAP_SAPRC99(124)  &
          / 'met_mercaptan   ', 124            /
      DATA  MECH_NAM_SAPRC99(124)     , MECH_MAP_SAPRC99(124)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99(124)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(125)     , SPMH_MAP_SAPRC99(125)  &
          / 'met_propenyl_2s ', 125            /
      DATA  MECH_NAM_SAPRC99(125)     , MECH_MAP_SAPRC99(125)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(125)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(126)     , SPMH_MAP_SAPRC99(126)  &
          / 'PPPP_2s         ', 126            /
      DATA  MECH_NAM_SAPRC99(126)     , MECH_MAP_SAPRC99(126)  &
          / 'OLE1            ', 25             /
      DATA  CONV_FAC_SAPRC99(126)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(127)     , SPMH_MAP_SAPRC99(127)  &
          / '2met_nonatriene ', 127            /
      DATA  MECH_NAM_SAPRC99(127)     , MECH_MAP_SAPRC99(127)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(127)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(128)     , SPMH_MAP_SAPRC99(128)  &
          / 'met_salicylate  ', 128            /
      DATA  MECH_NAM_SAPRC99(128)     , MECH_MAP_SAPRC99(128)  &
          / 'ARO1            ', 23             /
      DATA  CONV_FAC_SAPRC99(128)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(129)     , SPMH_MAP_SAPRC99(129)  &
          / 'indole          ', 129            /
      DATA  MECH_NAM_SAPRC99(129)     , MECH_MAP_SAPRC99(129)  &
          / 'ARO2            ', 24             /
      DATA  CONV_FAC_SAPRC99(129)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(130)     , SPMH_MAP_SAPRC99(130)  &
          / 'jasmone         ', 130            /
      DATA  MECH_NAM_SAPRC99(130)     , MECH_MAP_SAPRC99(130)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(130)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(131)     , SPMH_MAP_SAPRC99(131)  &
          / 'met_jasmonate   ', 131            /
      DATA  MECH_NAM_SAPRC99(131)     , MECH_MAP_SAPRC99(131)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(131)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(132)     , SPMH_MAP_SAPRC99(132)  &
          / '3met_3DCTT      ', 132            /
      DATA  MECH_NAM_SAPRC99(132)     , MECH_MAP_SAPRC99(132)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(132)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(133)     , SPMH_MAP_SAPRC99(133)  &
          / 'hexanal         ', 133            /
      DATA  MECH_NAM_SAPRC99(133)     , MECH_MAP_SAPRC99(133)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99(133)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(134)     , SPMH_MAP_SAPRC99(134)  &
          / 'hexanol_1       ', 134            /
      DATA  MECH_NAM_SAPRC99(134)     , MECH_MAP_SAPRC99(134)  &
          / 'ALK5            ', 22             /
      DATA  CONV_FAC_SAPRC99(134)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(135)     , SPMH_MAP_SAPRC99(135)  &
          / 'hexenal_c3      ', 135            /
      DATA  MECH_NAM_SAPRC99(135)     , MECH_MAP_SAPRC99(135)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99(135)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(136)     , SPMH_MAP_SAPRC99(136)  &
          / 'hexenal_t2      ', 136            /
      DATA  MECH_NAM_SAPRC99(136)     , MECH_MAP_SAPRC99(136)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99(136)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(137)     , SPMH_MAP_SAPRC99(137)  &
          / 'hexenol_c3      ', 137            /
      DATA  MECH_NAM_SAPRC99(137)     , MECH_MAP_SAPRC99(137)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99(137)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(138)     , SPMH_MAP_SAPRC99(138)  &
          / 'hexenyl_ACT_c3  ', 138            /
      DATA  MECH_NAM_SAPRC99(138)     , MECH_MAP_SAPRC99(138)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99(138)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(139)     , SPMH_MAP_SAPRC99(139)  &
          / 'homosalate      ', 139            /
      DATA  MECH_NAM_SAPRC99(139)     , MECH_MAP_SAPRC99(139)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(139)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(140)     , SPMH_MAP_SAPRC99(140)  &
          / 'Ehsalate        ', 140            /
      DATA  MECH_NAM_SAPRC99(140)     , MECH_MAP_SAPRC99(140)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(140)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(141)     , SPMH_MAP_SAPRC99(141)  &
          / 'pentanal        ', 141            /
      DATA  MECH_NAM_SAPRC99(141)     , MECH_MAP_SAPRC99(141)  &
          / 'RCHO            ', 27             /
      DATA  CONV_FAC_SAPRC99(141)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(142)     , SPMH_MAP_SAPRC99(142)  &
          / 'heptanone       ', 142            /
      DATA  MECH_NAM_SAPRC99(142)     , MECH_MAP_SAPRC99(142)  &
          / 'OLE2            ', 26             /
      DATA  CONV_FAC_SAPRC99(142)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(143)     , SPMH_MAP_SAPRC99(143)  &
          / 'anisole         ',143             /
      DATA  MECH_NAM_SAPRC99(143)     , MECH_MAP_SAPRC99(143)  &
          / 'BALD            ', 13             /
      DATA  CONV_FAC_SAPRC99(143)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(144)     , SPMH_MAP_SAPRC99(144)  &
          / 'verbenene       ',144             /
      DATA  MECH_NAM_SAPRC99(144)     , MECH_MAP_SAPRC99(144)  &
          / 'ARO2            ', 24             /
      DATA  CONV_FAC_SAPRC99(144)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(145)     , SPMH_MAP_SAPRC99(145)  &
          / 'benzyl-acetate  ',145             /
      DATA  MECH_NAM_SAPRC99(145)     , MECH_MAP_SAPRC99(145)  &
          / 'ARO1            ', 23             /
      DATA  CONV_FAC_SAPRC99(145)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(146)     , SPMH_MAP_SAPRC99(146)  &
          / 'myrtenal        ',146             /
      DATA  MECH_NAM_SAPRC99(146)     , MECH_MAP_SAPRC99(146)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(146)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(147)     , SPMH_MAP_SAPRC99(147)  &
          / 'benzyl-alcohol  ',147             /
      DATA  MECH_NAM_SAPRC99(147)     , MECH_MAP_SAPRC99(147)  &
          / 'ARO1            ', 23             /
      DATA  CONV_FAC_SAPRC99(147)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(148)     , SPMH_MAP_SAPRC99(148)  &
          / 'meta-cymenene   ',148             /
      DATA  MECH_NAM_SAPRC99(148)     , MECH_MAP_SAPRC99(148)  &
          / 'ARO2            ', 24             /
      DATA  CONV_FAC_SAPRC99(148)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(149)     , SPMH_MAP_SAPRC99(149)  &
          / 'ipsenol         ',149             /
      DATA  MECH_NAM_SAPRC99(149)     , MECH_MAP_SAPRC99(149)  &
          / 'TRP1            ', 2              /
      DATA  CONV_FAC_SAPRC99(149)  &
          / 1.0            /
      DATA  SPMH_NAM_SAPRC99(150)     , SPMH_MAP_SAPRC99(150)  &
          / 'Napthalene      ', 150            /
      DATA  MECH_NAM_SAPRC99(150)     , MECH_MAP_SAPRC99(150)  &
          / 'ARO2            ', 24             /
      DATA  CONV_FAC_SAPRC99(150)  &
          / 1.0            /

end module EMIS_MAPS_MEGAN
