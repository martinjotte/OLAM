
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

! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header: /project/yoj/arc/CCTM/src/depv/m3dry/DEPVVARS.F,v 1.7 2012/01/19 14:21:45 yoj Exp $

MODULE DEPVVARS

!-----------------------------------------------------------------------
! Name:     Dry Deposition Variables
! Purpose:  Contains arrays specific to dry deposition species.
!           Initializes dry deposition arrays.
! Revised:  19 Aug 2005  Original version.  (T. Otte and W. Hutzell)
!           25 Nov 2006 J.Young: combine
!           30 Apr 2008  Removed references to RADMDRY.  Added five air
!                        toxic species to LTOTG.  (T. Otte and W. Hutzell)
!           16 Feb 2011 S.Roselle: replaced I/O API include files with
!                        UTILIO_DEFN
!           11 May 2011 J.Bash: Updated for NH3 bidirectional exchange
!           04 Jan 2012 J.Young: Initialize char variables with blank padding
!           28 Aug 2014  G. Sarwar: added deposition for CLNO2 
!           07 Nov 14 J.Bash: Updated for the ASX_DATA_MOD shared data module. 
!           05 May 2015: H.Pye: Added ISOPNN and MTNO3 with Donna Schwede
!           24 Aug 2015: H.Pye: Added IEPOX and HACET with Donna Schwede
!           26 Jan 2016: H. Pye added info for SOA SVOCs
!           24 Mar 2016: G. Sarwar: added iodine and bromine species

!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: LDDEP = 4
  INTEGER, PARAMETER :: LTOTG = 178
  INTEGER, PARAMETER :: LCMP  = 8

  INTEGER, PARAMETER :: N_SPC_M3DRY = LTOTG
  CHARACTER( 16 )    :: DEPV_METHOD

  LOGICAL            :: USE_DEPSPC( LTOTG ) = .FALSE.
  CHARACTER( 16 )    :: DEPSPC    ( LTOTG )

  integer, parameter :: l_sulf =  2
  integer, parameter :: l_no2  =  3
  integer, parameter :: l_o3   =  5
  integer, parameter :: l_nh3  = 13
  integer, parameter :: l_hono = 15

  REAL, PARAMETER    :: a0     = 8.0     ! [dim'less]
  REAL, PARAMETER    :: dwat   = 0.2178  ! [cm^2/s] at 273.15K
  REAL, PARAMETER    :: kvis   = 0.132   ! [cm^2 / s] at 273.15K
  REAL, PARAMETER    :: pr     = 0.709   ! [dim'less]
  REAL, PARAMETER    :: rsnow0 = 1000.0
  REAL, PARAMETER    :: rg0    = 1000.0  ! [s/m]
  REAL, PARAMETER    :: rcut0  = 3000.0  ! [s/m]

  REAL               :: ar       ( ltotg ) ! reactivity relative to HNO3
  REAL               :: dif0     ( ltotg ) ! molecular diffusivity [cm2/s]
  REAL               :: ddif0    ( ltotg ) ! molecular diffusivity / dwat
  REAL               :: lebas    ( ltotg ) ! Le Bas molar volume [cm3/mol ]
  REAL               :: meso     ( ltotg ) ! Exception for species that
                                                          ! react with cell walls.
                                                          ! Wesely 1989 eq 6.
  REAL               :: scc_pr_23( ltotg ) ! 0.2(PR/SCC)**2/3, fn of DIF0
  real               :: csnow    ( ltotg )
  real               :: cgnd     ( ltotg )
  real               :: ccut     ( ltotg )
  REAL               :: dw25     ( ltotg ) ! diffusivity of water at 25C

  CHARACTER( 16 )    :: subname  ( ltotg ) ! for subroutine HLCONST
  integer            :: hlspc    ( ltotg ) ! species index in hlconst
  logical            :: hleff    ( ltotg )

!-------------------------------------------------------------------------------
!-- Chemical-Dependent Parameters (Original Source: Modified ADOM - Padro)
!                                                        at 298.15 K
!     Species   Dif0(cm2/s) Alphastar   Reactivity     -DKHOR [K]  KH0 [M/atm]
!     _______   ___________ _________   __________     _________  __________
!  1   SO2      0.1089 ~    1000.       8.00            3100.0 ~     1.2e00 ~
!  2   SULF      --       --           --                --           --
!  3   NO2      0.1361 ~    1.00      8.00 2.00*        2500.0 ~     1.2e-2 ~
!  4   NO       0.1802 ~    1         2.0+              1500.0 ~     1.9e-3 ~
!  5   O3       0.1444 ~    10.00       15.0  8@        2700.0 ~     1.2e-2 ~
!  6   HNO3     0.1628 ~    1E9       18.0 800* 8000**  8700.0 ~     2.6e+6 ~
!  7   H2O2     0.2402 ~    1.00      12.0 30*          6600.0 ~     9.7e+4 ~
!  8   ACET ALD 0.1525 ~    1         10.+              5700.0 ~     1.3e+1 ~
!  9   HCHO     0.1877 ~    1.00      10.+              6400.0 ~     7.0e+3 ~
! 10   OP       0.1525 ~    1         10.0+             5200.0 ~     3.1e+2 ~
! 11   PAA      0.1220 ~    1         20+               5300.0 ~     8.4e+2 ~
! 12   ORA      0.1525 ~    1         20+               5700.0 ~     3.7e+3 ~
! 13   NH3      0.1978 ~    1E5       10.0              4100.0 ~     5.8e+1 ~
! 14   PAN      0.0938 ~    1         16.0~~            5900.0 ~     2.9e00 ~
! 15   HONO     0.1525 ~    1         20+               4800.0 ~     4.9e+1 ~
! 16   CO       0.1807 ~    1         5.0               1600.0 ~     9.5e-4 ~
! 17   METHANOL 0.1363 ~    1.0 ~       2.0 ~           5200.0 ~     2.2e+2 ~
! --   CO2      0.1381 ~                                2400.0 ~     3.4e-2 ~
! 18   N2O5     0.0808 ^^             5000.0**
! 19   NO3      0.1153 ^^             5000.0**
! 20   GEN ALD  0.1525 ##   1.0       10.0##
! 21   CL2      0.1080 %              10.0 %
! 22   HOCL     0.1300 %              10.0 %
! 23   HCL      0.1510 %              8000.0 %
! 24   FMCL     0.1094 %              10.0 %
! 27   HG       0.1194 $              0.1 $
! 28   HGIIGAS  0.0976 $              8000.0 $
! 50   HEXAMETHYLENE DIISOCYANATE       10.0 <>
! 51   HYDRAZINE                        20.0 <>
! 52   MALEIC ANHYDRIDE                 10.0 <>
! 53   TOLUENE DIISOCYANATE             10.0 <>
! 54   TRIETHYLAMINE                    10.0 <>
!
!---------Notes
!  * Updates based on literature review 7/96 JEP
!  # Diff and H based on Wesely (1988) same as RADM
!  + Estimated by JEP 2/97
!  @ Updated by JEP 9/01
!  ~ Added by YW 1/02.  Dif0 based on Massman (1998).  Henry's Law constant
!    is defined here as: h=cg/ca, where cg is the concentration of a species
!    in gas-phase, and ca is its aqueous-phase concentration.  The smaller h,
!    the larger solubility.  Henry's Law constant in another definition (KH):
!    KH = ca/pg [M/atm], KH = KH0 * exp(-DKH/R(1/T-1/T0)), where KH0 and -DKH
!    values are from Rolf Sander (1999).  h=1/(KH*R*T).
! ** Update by DBS based on estimates by JEP 1/03
! ^^ From Bill Massman, personal communication 4/03
! ## Diffusivity calculated by SPARC, reactivity = other aldehydes
! ++ Dif0 in Massman is diffusivity at temperature 0C and 1 atm (101.325kPa), so
!    chemicals that were not in Massman's paper need to be adjusted.  We assume
!    JEP's original values were for 25C and 1 atm.
!  % Added by G. Sarwar (10/04)
!  $ Added by R. Bullock (02/05) HG diffusivity is from Massman (1999).
!    HGIIGAS diffusivity calculated from the HG value and a mol. wt. scaling
!    factor of MW**(-2/3) from EPA/600/3-87/015. ORD, Athens, GA.  HGIIGAS
!    mol.wt. used is that of HgCl2.  Reactivity of HG is 1/20th of NO and NO2
!    values based on general atmospheric lifetimes of each species.  Reactivity
!    of HGIIGAS is based on HNO3 surrogate.
! @@ Mesophyll resistances for NO, NO2, and CO added by J. Pleim (07/07) based
!    on values in Pleim, Venkatram, and Yamartino, 1984:  ADOM/TADAP Model
!    Development Program, Volume 4, The Dry Deposition Module.  ERT, Inc.,
!    Concord, MA (peer reviewed).
! ~~ Reactivity for PAN changed from 4.0 to 16.0 by J. Pleim (07/07) based on
!    comparisons with Turnipseed et al., JGR, 2006.
! %% Species ICL1 and ICL2 are removed, not used in CB05.  G. Sarwar (07/07)
! <> Hazardous Air Pollutants that are believed to undergo significant dry
!    deposition. Hydrazine and triethylamine reactivities are based on analogies
!    to NH3. Maleic anhydride reactivity is assumed similar to aldehydes.
!    Toluene diisocyanate and hexamethylene diisocyanate reactivities are
!    assumed to be similar to SO2. Diffusivities are calculated with standard
!    formulas.  W. Hutzell (04/08)
! %% G. Sarwar: added data for iodine and bromine species (03/2016)
! %% B. Hutzell: added dry deposition data for methane, acrylic acid, methyl chloride,
!    and acetonitrile (09/2016)

  DATA subname(  1), dif0(  1), ar(  1), meso(  1), lebas(  1) / 'SO2             ', 0.1089,   10.0,      0.0,  35.0/
  DATA subname(  2), dif0(  2), ar(  2), meso(  2), lebas(  2) / 'H2SO4           ', 0.1091, 8000.0,      0.0,  49.0/
  DATA subname(  3), dif0(  3), ar(  3), meso(  3), lebas(  3) / 'NO2             ', 0.1361,    2.0,      0.1,  21.0/
  DATA subname(  4), dif0(  4), ar(  4), meso(  4), lebas(  4) / 'NO              ', 0.1802,    2.0,      0.0,  14.0/
  DATA subname(  5), dif0(  5), ar(  5), meso(  5), lebas(  5) / 'O3              ', 0.1444,   12.0,      1.0,  21.0/
  DATA subname(  6), dif0(  6), ar(  6), meso(  6), lebas(  6) / 'HNO3            ', 0.1067, 8000.0,      0.0,  35.0/
  !ar=34,000 such that r_cut=0.7 s/m as in Nguyen et al. 2015:
  DATA subname(  7), dif0(  7), ar(  7), meso(  7), lebas(  7) / 'H2O2            ', 0.1300,34000.0,      1.0,  28.0/
  DATA subname(  8), dif0(  8), ar(  8), meso(  8), lebas(  8) / 'ACETALDEHYDE    ', 0.1111,   10.0,      0.0,  56.0/
  DATA subname(  9), dif0(  9), ar(  9), meso(  9), lebas(  9) / 'FORMALDEHYDE    ', 0.1554,   10.0,      0.0,  35.0/
  !meso change from 0.1 to 0.3, Wolfe and Thornton 2011 ACP per J. Bash:
  DATA subname( 10), dif0( 10), ar( 10), meso( 10), lebas( 10) / 'METHYLHYDROPEROX', 0.1179,   10.0,      0.3,  49.0/
  DATA subname( 11), dif0( 11), ar( 11), meso( 11), lebas( 11) / 'PEROXYACETIC_ACI', 0.0868,   20.0,      0.1,  70.0/
  DATA subname( 12), dif0( 12), ar( 12), meso( 12), lebas( 12) / 'ACETIC_ACID     ', 0.0944,   20.0,      0.0,  63.0/
  DATA subname( 13), dif0( 13), ar( 13), meso( 13), lebas( 13) / 'NH3             ', 0.1978,   20.0,      0.0,  28.0/
  DATA subname( 14), dif0( 14), ar( 14), meso( 14), lebas( 14) / 'PAN             ', 0.0687,   16.0,      0.1,  91.0/
  DATA subname( 15), dif0( 15), ar( 15), meso( 15), lebas( 15) / 'HNO2            ', 0.1349,   20.0,      0.1,  28.0/
  DATA subname( 16), dif0( 16), ar( 16), meso( 16), lebas( 16) / 'CO              ', 0.1807,    5.0,      0.0,  14.0/
  DATA subname( 17), dif0( 17), ar( 17), meso( 17), lebas( 17) / 'METHANOL        ', 0.1329,    2.0,      0.0,  42.0/
  DATA subname( 18), dif0( 18), ar( 18), meso( 18), lebas( 18) / 'N2O5            ', 0.0808, 5000.0,      0.0,  49.0/
  DATA subname( 19), dif0( 19), ar( 19), meso( 19), lebas( 19) / 'NO3             ', 0.1153, 5000.0,      0.0,  28.0/
  DATA subname( 20), dif0( 20), ar( 20), meso( 20), lebas( 20) / 'GENERIC_ALDEHYDE', 0.0916,   10.0,      0.0,  56.0/
  DATA subname( 21), dif0( 21), ar( 21), meso( 21), lebas( 21) / 'CL2             ', 0.1080,   10.0,      0.0,  49.0/
  DATA subname( 22), dif0( 22), ar( 22), meso( 22), lebas( 22) / 'HOCL            ', 0.1300,   10.0,      0.0,  38.5/
  DATA subname( 23), dif0( 23), ar( 23), meso( 23), lebas( 23) / 'HCL             ', 0.1510, 8000.0,      0.0,  31.5/
  DATA subname( 24), dif0( 24), ar( 24), meso( 24), lebas( 24) / 'FMCL            ', 0.1094,   10.0,      0.0,  45.5/
  DATA subname( 25), dif0( 25), ar( 25), meso( 25), lebas( 25) / 'HG              ', 0.1194,    0.1,      0.0,  14.8/
  ! estimation from back calculating to get dw25 = 1.04e-5 (Garland et al, 1965):
  DATA subname( 26), dif0( 26), ar( 26), meso( 26), lebas( 26) / 'HGIIGAS         ', 0.0976, 8000.0,      0.0,  95.0/
  DATA subname( 27), dif0( 27), ar( 27), meso( 27), lebas( 27) / 'TECDD_2378      ', 0.0525,    2.0,      0.0, 217.0/
  DATA subname( 28), dif0( 28), ar( 28), meso( 28), lebas( 28) / 'PECDD_12378     ', 0.0508,    2.0,      0.0, 234.5/
  DATA subname( 29), dif0( 29), ar( 29), meso( 29), lebas( 29) / 'HXCDD_123478    ', 0.0494,    2.0,      0.0, 252.0/
  DATA subname( 30), dif0( 30), ar( 30), meso( 30), lebas( 30) / 'HXCDD_123678    ', 0.0494,    2.0,      0.0, 252.0/
  DATA subname( 31), dif0( 31), ar( 31), meso( 31), lebas( 31) / 'HXCDD_123478    ', 0.0494,    2.0,      0.0, 252.0/
  DATA subname( 32), dif0( 32), ar( 32), meso( 32), lebas( 32) / 'HPCDD_1234678   ', 0.0480,    2.0,      0.0, 269.5/
  DATA subname( 33), dif0( 33), ar( 33), meso( 33), lebas( 33) / 'OTCDD           ', 0.0474,    2.0,      0.0, 287.0/
  DATA subname( 34), dif0( 34), ar( 34), meso( 34), lebas( 34) / 'TECDF_2378      ', 0.0534,    2.0,      0.0, 210.0/
  DATA subname( 35), dif0( 35), ar( 35), meso( 35), lebas( 35) / 'PECDF_12378     ', 0.0517,    2.0,      0.0, 227.5/
  DATA subname( 36), dif0( 36), ar( 36), meso( 36), lebas( 36) / 'PECDF_23478     ', 0.0517,    2.0,      0.0, 227.5/
  DATA subname( 37), dif0( 37), ar( 37), meso( 37), lebas( 37) / 'HXCDF_123478    ', 0.0512,    2.0,      0.0, 245.0/
  DATA subname( 38), dif0( 38), ar( 38), meso( 38), lebas( 38) / 'HXCDF_123678    ', 0.0512,    2.0,      0.0, 245.0/
  DATA subname( 39), dif0( 39), ar( 39), meso( 39), lebas( 39) / 'HXCDF_234678    ', 0.0512,    2.0,      0.0, 245.0/
  DATA subname( 40), dif0( 40), ar( 40), meso( 40), lebas( 40) / 'HXCDF_123789    ', 0.0512,    2.0,      0.0, 245.0/
  DATA subname( 41), dif0( 41), ar( 41), meso( 41), lebas( 41) / 'HPCDF_1234678   ', 0.0487,    2.0,      0.0, 262.5/
  DATA subname( 42), dif0( 42), ar( 42), meso( 42), lebas( 42) / 'HPCDF_1234789   ', 0.0487,    2.0,      0.0, 262.5/
  DATA subname( 43), dif0( 43), ar( 43), meso( 43), lebas( 43) / 'OTCDF           ', 0.0474,    2.0,      0.0, 280.0/
  DATA subname( 44), dif0( 44), ar( 44), meso( 44), lebas( 44) / 'NAPHTHALENE     ', 0.0778,    4.0,      0.0, 119.0/
  DATA subname( 45), dif0( 45), ar( 45), meso( 45), lebas( 45) / '1NITRONAPHTHALEN', 0.0692,    4.0,      0.0, 133.0/
  DATA subname( 46), dif0( 46), ar( 46), meso( 46), lebas( 46) / '2NITRONAPHTHALEN', 0.0692,    4.0,      0.0, 133.0/
  DATA subname( 47), dif0( 47), ar( 47), meso( 47), lebas( 47) / '14NAPHTHOQUINONE', 0.0780,    4.0,      0.0, 119.0/
  DATA subname( 48), dif0( 48), ar( 48), meso( 48), lebas( 48) / 'HEXAMETHYLE_DIIS', 0.0380,   10.0,      0.0, 196.0/
  DATA subname( 49), dif0( 49), ar( 49), meso( 49), lebas( 49) / 'HYDRAZINE       ', 0.4164,   20.0,      0.0,  42.0/
  DATA subname( 50), dif0( 50), ar( 50), meso( 50), lebas( 50) / 'MALEIC_ANHYDRIDE', 0.0950,   10.0,      0.0,  70.0/
  DATA subname( 51), dif0( 51), ar( 51), meso( 51), lebas( 51) / '24-TOLUENE_DIIS ', 0.0610,   10.0,      0.0, 154.0/
  DATA subname( 52), dif0( 52), ar( 52), meso( 52), lebas( 52) / 'TRIETHYLAMINE   ', 0.0881,   20.0,      0.0, 154.0/
  ! assumes 58.2% C5H11O4N and 41.8% C5H11O3N:
  DATA subname( 53), dif0( 53), ar( 53), meso( 53), lebas( 53) / 'ORG_NTR         ', 0.0607,   16.0,      0.0, 160.0/
  DATA subname( 54), dif0( 54), ar( 54), meso( 54), lebas( 54) / 'HYDROXY_NITRATES', 0.0609,   16.0,      0.0, 156.1/
  DATA subname( 55), dif0( 55), ar( 55), meso( 55), lebas( 55) / 'MPAN            ', 0.0580,   16.0,      0.1, 133.0/
  DATA subname( 56), dif0( 56), ar( 56), meso( 56), lebas( 56) / 'PPN             ', 0.0631,   16.0,      0.1, 118.2/
  DATA subname( 57), dif0( 57), ar( 57), meso( 57), lebas( 57) / 'MVK             ', 0.0810,    8.0,      1.0,  88.8/
  DATA subname( 58), dif0( 58), ar( 58), meso( 58), lebas( 58) / 'DINTR           ', 0.0617,   16.0,      0.1, 169.8/
  DATA subname( 59), dif0( 59), ar( 59), meso( 59), lebas( 59) / 'NTR_ALK         ', 0.0688,   16.0,      0.1, 133.0/
  DATA subname( 60), dif0( 60), ar( 60), meso( 60), lebas( 60) / 'NTR_OH          ', 0.0665,   16.0,      0.1, 140.4/
  DATA subname( 61), dif0( 61), ar( 61), meso( 61), lebas( 61) / 'HYDROXY_NITRATES', 0.0646,   16.0,      0.0, 147.8/
  DATA subname( 62), dif0( 62), ar( 62), meso( 62), lebas( 62) / 'PROPNN          ', 0.0677,   16.0,      0.0, 133.0/
  ! dif0 estimated following Erickson III et al., JGR, 104, D7, 8347-8372, 1999:
  DATA subname( 63), dif0( 63), ar( 63), meso( 63), lebas( 63) / 'NITRYL_CHLORIDE ', 0.0888,    8.0,      0.0,  45.5/
  DATA subname( 64), dif0( 64), ar( 64), meso( 64), lebas( 64) / 'ISOPNN          ',0.0457,    8.0,      0.0,  206.8/
  DATA subname( 65), dif0( 65), ar( 65), meso( 65), lebas( 65) / 'MTNO3           ',0.0453,    8.0,      0.0,  251.2/
  DATA subname( 66), dif0( 66), ar( 66), meso( 66), lebas( 66) / 'IEPOX           ',0.0579,    8.0,      0.0,  110.8/
  ! dif0 from Nguyen 2015 PNAS:
  DATA subname( 67), dif0( 67), ar( 67), meso( 67), lebas( 67) / 'HACET           ',0.1060,    8.0,      0.0,   72.6/
  DATA subname( 68), dif0( 68), ar( 68), meso( 68), lebas( 68) / 'SVALK1          ',0.0514,   20.0,      0.0,  280.5/
  DATA subname( 69), dif0( 69), ar( 69), meso( 69), lebas( 69) / 'SVALK2          ',0.0546,   20.0,      0.0,  275.6/
  DATA subname( 70), dif0( 70), ar( 70), meso( 70), lebas( 70) / 'SVBNZ1          ',0.0642,   20.0,      0.0,  134.1/
  DATA subname( 71), dif0( 71), ar( 71), meso( 71), lebas( 71) / 'SVBNZ2          ',0.0726,   20.0,      0.0,  127.5/
  DATA subname( 72), dif0( 72), ar( 72), meso( 72), lebas( 72) / 'SVISO1          ',0.0733,   20.0,      0.0,  126.3/
  DATA subname( 73), dif0( 73), ar( 73), meso( 73), lebas( 73) / 'SVISO2          ',0.0729,   20.0,      0.0,  123.8/
  DATA subname( 74), dif0( 74), ar( 74), meso( 74), lebas( 74) / 'SVPAH1          ',0.0564,   20.0,      0.0,  235.7/
  DATA subname( 75), dif0( 75), ar( 75), meso( 75), lebas( 75) / 'SVPAH2          ',0.0599,   20.0,      0.0,  231.5/
  DATA subname( 76), dif0( 76), ar( 76), meso( 76), lebas( 76) / 'SVSQT           ',0.0451,   20.0,      0.0,  346.5/
  DATA subname( 77), dif0( 77), ar( 77), meso( 77), lebas( 77) / 'SVTOL1          ',0.0637,   20.0,      0.0,  153.7/
  DATA subname( 78), dif0( 78), ar( 78), meso( 78), lebas( 78) / 'SVTOL2          ',0.0607,   20.0,      0.0,  194.1/
  DATA subname( 79), dif0( 79), ar( 79), meso( 79), lebas( 79) / 'SVTRP1          ',0.0603,   20.0,      0.0,  194.9/
  DATA subname( 80), dif0( 80), ar( 80), meso( 80), lebas( 80) / 'SVTRP2          ',0.0559,   20.0,      0.0,  218.8/
  DATA subname( 81), dif0( 81), ar( 81), meso( 81), lebas( 81) / 'SVXYL1          ',0.0610,   20.0,      0.0,  154.6/
  DATA subname( 82), dif0( 82), ar( 82), meso( 82), lebas( 82) / 'SVXYL2          ',0.0585,   20.0,      0.0,  194.6/
  DATA subname( 83), dif0( 83), ar( 83), meso( 83), lebas( 83) / 'IO              ',0.1002,    8.0,      0.0,   44.4/
  DATA subname( 84), dif0( 84), ar( 84), meso( 84), lebas( 84) / 'OIO             ',0.0938,    8.0,      0.0,   51.8/
  DATA subname( 85), dif0( 85), ar( 85), meso( 85), lebas( 85) / 'I2O2            ',0.0732,    8.0,      0.0,   88.8/
  DATA subname( 86), dif0( 86), ar( 86), meso( 86), lebas( 86) / 'I2O3            ',0.0707,    8.0,      0.0,   96.2/
  DATA subname( 87), dif0( 87), ar( 87), meso( 87), lebas( 87) / 'I2O4            ',0.0684,    8.0,      0.0,  103.6/
  DATA subname( 88), dif0( 88), ar( 88), meso( 88), lebas( 88) / 'HI              ',0.1045,    8.0,      0.0,   40.7/
  DATA subname( 89), dif0( 89), ar( 89), meso( 89), lebas( 89) / 'HOI             ',0.0972,    8.0,      0.0,   48.1/
  DATA subname( 90), dif0( 90), ar( 90), meso( 90), lebas( 90) / 'INO             ',0.0882,    8.0,      0.0,   60.9/
  DATA subname( 91), dif0( 91), ar( 91), meso( 91), lebas( 91) / 'INO2            ',0.0883,   20.0,      0.0,   69.2/
  DATA subname( 92), dif0( 92), ar( 92), meso( 92), lebas( 92) / 'IONO2           ',0.0792,    8.0,      0.0,   77.5/
  DATA subname( 93), dif0( 93), ar( 93), meso( 93), lebas( 93) / 'BRO             ',0.1144,    1.0,      0.0,   34.4/
  DATA subname( 94), dif0( 94), ar( 94), meso( 94), lebas( 94) / 'HOBR            ',0.1101,    1.0,      0.0,   38.1/
  DATA subname( 95), dif0( 95), ar( 95), meso( 95), lebas( 95) / 'HBR             ',0.1216,    2.0,      0.0,   30.7/
  DATA subname( 96), dif0( 96), ar( 96), meso( 96), lebas( 96) / 'BRONO2          ',0.0855,    1.0,      0.0,   67.5/
  DATA subname( 97), dif0( 97), ar( 97), meso( 97), lebas( 97) / 'BRNO2           ',0.0909,    1.0,      0.0,   59.2/
  DATA subname( 98), dif0( 98), ar( 98), meso( 98), lebas( 98) / 'BRCL            ',0.0966,    1.0,      0.0,   51.6/
  DATA subname( 99), dif0( 99), ar( 99), meso( 99), lebas( 99) / 'DMS             ',0.0926,    2.0,      0.0,   77.4/
  DATA subname(100), dif0(100), ar(100), meso(100), lebas(100) / 'MSA             ',0.0896,    2.0,      0.0,   77.4/
  ! dif0, equation 9-22. Scwarzenbach et. (1993) Env. Org. Chem.:
  DATA subname(101), dif0(101), ar(101), meso(101), lebas(101) / 'METHANE         ',0.2107,    2.0,      0.0,   29.6/
  DATA subname(102), dif0(102), ar(102), meso(102), lebas(102) / 'ACRYACID        ',0.0908,    2.0,      0.0,   63.2/
  DATA subname(103), dif0(103), ar(103), meso(103), lebas(103) / 'CARBSULFIDE     ',0.1240,    5.0,      0.0,   51.5/
  DATA subname(104), dif0(104), ar(104), meso(104), lebas(104) / 'ACETONITRILE    ',0.1280,    5.0,      0.0,   52.3/
  ! dif0, equation 9-22. Scwarzenbach et. (1993) Env. Org. Chem.:
  DATA subname(105), dif0(105), ar(105), meso(105), lebas(105) / '6_NITRO_O_CRESOL',0.0664,   16.0,      0.0,  155.0/
  DATA subname(106), dif0(106), ar(106), meso(106), lebas(106) / 'GENERIC_ALDEHYDE',0.0646,   10.0, 0.0,   56.0 / ! PCVOC
  DATA subname(107), dif0(107), ar(107), meso(107), lebas(107) / 'NTR_OH          ',0.0722,   16.0, 0.1,  140.4 / ! INTR
  DATA subname(108), dif0(108), ar(108), meso(108), lebas(108) / 'METHYLHYDROPEROX',0.0853,   10.0, 0.3,   49.0 / ! ISPX diffusion should be ~ 0.0710 according to Wolfe and thornton 2011 ACP
  DATA subname(109), dif0(109), ar(109), meso(109), lebas(109) / 'METHYLHYDROPEROX',0.1371,   10.0, 0.3,   49.0 / ! ROOH diffusion should be ~ 0.0710 according to Wolfe and thornton 2011 ACP
  DATA subname(110), dif0(110), ar(110), meso(110), lebas(110) / 'ADIPIC_ACID     ',0.0646,90000.0, 0.0,   63.0 / ! LVPCSOG
  DATA subname(111), dif0(111), ar(111), meso(111), lebas(111) / 'ADIPIC_ACID     ',0.0456,    4.2, 0.0,   63.0 / ! VIVPO1
  DATA subname(112), dif0(112), ar(112), meso(112), lebas(112) / 'ADIPIC_ACID     ',0.0766,71624.8, 0.0,   63.0 / ! VLVOO1
  DATA subname(113), dif0(113), ar(113), meso(113), lebas(113) / 'ADIPIC_ACID     ',0.0766, 9042.0, 0.0,   63.0 / ! VLVOO2
  DATA subname(114), dif0(114), ar(114), meso(114), lebas(114) / 'ADIPIC_ACID     ',0.0533,13818.0, 0.0,   63.0 / ! VLVPO1
  DATA subname(115), dif0(115), ar(115), meso(115), lebas(115) / 'ADIPIC_ACID     ',0.0771, 1133.9, 0.0,   63.0 / ! VSVOO1
  DATA subname(116), dif0(116), ar(116), meso(116), lebas(116) / 'ADIPIC_ACID     ',0.0771,   18.1, 0.0,   63.0 / ! VSVOO2
  DATA subname(117), dif0(117), ar(117), meso(117), lebas(117) / 'ADIPIC_ACID     ',0.0775,    2.3, 0.0,   63.0 / ! VSVOO3
  DATA subname(118), dif0(118), ar(118), meso(118), lebas(118) / 'ADIPIC_ACID     ',0.0511, 1830.5, 0.0,   63.0 / ! VSVPO1
  DATA subname(119), dif0(119), ar(119), meso(119), lebas(119) / 'ADIPIC_ACID     ',0.0493,  241.0, 0.0,   63.0 / ! VSVPO2
  DATA subname(120), dif0(120), ar(120), meso(120), lebas(120) / 'ADIPIC_ACID     ',0.0474,   31.8, 0.0,   63.0 / ! VSVPO3
  DATA subname(121), dif0(121), ar(121), meso(121), lebas(121) / 'FORMIC_ACID     ',0.1411,   20.0, 0.0,   63.0 / ! FACD
  DATA subname(122), dif0(122), ar(122), meso(122), lebas(122) / 'MEK             ',0.0859,    1.0, 0.0,  108.2 / ! KET different in different mechanisms
  DATA subname(123), dif0(123), ar(123), meso(123), lebas(123) / 'ETHENE          ',0.1366,    1.0, 0.0,   58.1 / ! ETH
  DATA subname(124), dif0(124), ar(124), meso(124), lebas(124) / 'HNO4            ',0.1233,    1.0, 0.0,   45.2 / ! PNA
  DATA subname(125), dif0(125), ar(125), meso(125), lebas(125) / 'GLYOXAL         ',0.1188,    1.0, 0.0,   56.2 / ! GLY
  DATA subname(126), dif0(126), ar(126), meso(126), lebas(126) / 'GLYOXAL         ',0.1181,    1.0, 0.0,   56.4 / ! GLYD
  DATA subname(127), dif0(127), ar(127), meso(127), lebas(127) / 'METHYL_GLYOXAL  ',0.1038,    1.0, 0.0,   72.5 / ! MGLY
  DATA subname(128), dif0(128), ar(128), meso(128), lebas(128) / 'ETHANE          ',0.1312,    1.0, 0.0,   61.5 / ! ETHA
  DATA subname(129), dif0(129), ar(129), meso(129), lebas(129) / 'ETHANOL         ',0.1213,    1.0, 0.0,   59.1 / ! ETOH
  DATA subname(130), dif0(130), ar(130), meso(130), lebas(130) / 'ETHANE          ',0.0870,    1.0, 0.0,  111.1 / ! PAR as Pentane
  DATA subname(131), dif0(131), ar(131), meso(131), lebas(131) / 'ACETONE         ',0.1057,    1.0, 0.0,   75.2 / ! ACET
  DATA subname(132), dif0(132), ar(132), meso(132), lebas(132) / 'PROPANE         ',0.1095,    1.0, 0.0,   78.1 / ! PRPA
  DATA subname(133), dif0(133), ar(133), meso(133), lebas(133) / 'ACETYLENE       ',0.1523,    1.0, 0.0,   45.8 / ! ETHY
  DATA subname(134), dif0(134), ar(134), meso(134), lebas(134) / 'ETHENE          ',0.1135,    1.0, 0.0,   73.1 / ! OLE as Propene
  DATA subname(135), dif0(135), ar(135), meso(135), lebas(135) / 'ETHENE          ',0.0990,    1.0, 0.0,   89.5 / ! IOLE as Isobutene
  DATA subname(136), dif0(136), ar(136), meso(136), lebas(136) / 'MEK             ',0.0852,    1.0, 0.0,  101.2 / ! IEPOX different scavenging H in CB05 and
  DATA subname(137), dif0(137), ar(137), meso(137), lebas(137) / 'BENZENE         ',0.0942,    1.0, 0.0,   89.4 / ! BENZENE
  DATA subname(138), dif0(138), ar(138), meso(138), lebas(138) / '2-CRESOL        ',0.0850,    1.0, 0.0,  108.1 / ! CRES
  DATA subname(139), dif0(139), ar(139), meso(139), lebas(139) / 'TOLUENE         ',0.0860,    1.0, 0.0,  105.7 / ! TOL
  DATA subname(140), dif0(140), ar(140), meso(140), lebas(140) / 'O-XYLENE        ',0.0796,    1.0, 0.0,  122.0 / ! XYLMN
  DATA subname(141), dif0(141), ar(141), meso(141), lebas(141) / 'O-XYLENE        ',0.0777,    1.0, 0.0,  123.5 / ! NAPH
  DATA subname(142), dif0(142), ar(142), meso(142), lebas(142) / 'PHENOL          ',0.0844,    1.0, 0.0,  102.6 / ! CAT1
  DATA subname(143), dif0(143), ar(143), meso(143), lebas(143) / 'PINENE          ',0.0545,    1.0, 0.0,  251.5 / ! SESQ
  DATA subname(144), dif0(144), ar(144), meso(144), lebas(144) / 'PINENE          ',0.0700,    1.0, 0.0,  136.2 / ! TERP
  DATA subname(145), dif0(145), ar(145), meso(145), lebas(145) / 'ISOPRENE        ',0.0913,    1.0, 0.0,  136.2 / ! ISOP
  DATA subname(146), dif0(146), ar(146), meso(146), lebas(146) / 'METHACROLEIN    ',0.1033,    1.0, 0.0,   69.6 / ! OPEN C4H4O2
  DATA subname(147), dif0(147), ar(147), meso(147), lebas(147) / 'MEK             ',0.0950,    1.0, 0.0,   81.7 / ! XOPN C5H6O2
  DATA subname(148), dif0(148), ar(148), meso(148), lebas(148) / 'DECANE          ',0.0739,    1.0, 0.0,  142.8 / ! SOAALK as Propylcyclopentane
  DATA subname(149), dif0(149), ar(149), meso(149), lebas(149) / '13-BUTADIENE    ',0.1019,    1.0, 0.0,   84.8 / ! BUTADIENE13
  DATA subname(150), dif0(150), ar(150), meso(150), lebas(150) / 'ACROLEIN        ',0.1092,    1.0, 0.0,   70.5 /
  DATA subname(151), dif0(151), ar(151), meso(151), lebas(151) / 'SVMT1           ',0.0424,   20.0,      0.0, 355.2/ ! see Xu et al., 2018 ACPD: doi:10.5194/acp-2017-1109
  DATA subname(152), dif0(152), ar(152), meso(152), lebas(152) / 'SVMT2           ',0.0556,   20.0,      0.0, 236.8/
  DATA subname(153), dif0(153), ar(153), meso(153), lebas(153) / 'SVMT3           ',0.0583,   20.0,      0.0, 214.6/
  DATA subname(154), dif0(154), ar(154), meso(154), lebas(154) / 'SVMT4           ',0.0587,   20.0,      0.0, 229.4/
  DATA subname(155), dif0(155), ar(155), meso(155), lebas(155) / 'SVMT5           ',0.0619,   20.0,      0.0, 207.2/
  DATA subname(156), dif0(156), ar(156), meso(156), lebas(156) / 'SVMT6           ',0.0624,   20.0,      0.0, 222.0/
  DATA subname(157), dif0(157), ar(157), meso(157), lebas(157) / 'SVMT7           ',0.0661,   20.0,      0.0, 199.8/
  DATA subname(158), dif0(158), ar(158), meso(158), lebas(158) / 'SVAVB1          ',0.0560,100388.0,     0.0, 163.1/
  DATA subname(159), dif0(159), ar(159), meso(159), lebas(159) / 'SVAVB2          ',0.0600,  1461.2,     0.0, 163.2/
  DATA subname(160), dif0(160), ar(160), meso(160), lebas(160) / 'SVAVB3          ',0.0620,   175.2,     0.0, 163.0/
  DATA subname(161), dif0(161), ar(161), meso(161), lebas(161) / 'SVAVB4          ',0.0650,    20.8,     0.0, 162.7/
  DATA subname(162), dif0(162), ar(162), meso(162), lebas(162) / 'CLNO3           ',0.0902,    8.0,      0.0,  52.5/
  DATA subname(163), dif0(163), ar(163), meso(163), lebas(163) / 'FMBR            ',0.0965,   10.0,      0.0,  52.5/
  DATA subname(164), dif0(164), ar(164), meso(164), lebas(164) / 'I2              ',0.0795,    4.0,      0.0,  77.0/
  DATA subname(165), dif0(165), ar(165), meso(165), lebas(165) / 'CH3I            ',0.0881,    2.0,      0.0,  66.5/
  DATA subname(166), dif0(166), ar(166), meso(166), lebas(166) / 'ICL             ',0.0878,    4.0,      0.0,  63.0/
  DATA subname(167), dif0(167), ar(167), meso(167), lebas(167) / 'IBR             ',0.0851,    4.0,      0.0,  70.0/
  DATA subname(168), dif0(168), ar(168), meso(168), lebas(168) / 'MI2             ',0.0713,    2.0,      0.0,  98.0/
  DATA subname(169), dif0(169), ar(169), meso(169), lebas(169) / 'MIB             ',0.0753,    2.0,      0.0,  91.0/
  DATA subname(170), dif0(170), ar(170), meso(170), lebas(170) / 'MIC             ',0.0773,    2.0,      0.0,  84.0/
  DATA subname(171), dif0(171), ar(171), meso(171), lebas(171) / 'BR2             ',0.0925,    2.0,      0.0,  63.0/
  DATA subname(172), dif0(172), ar(172), meso(172), lebas(172) / 'MB3             ',0.0705,    2.0,      0.0, 108.5/
  DATA subname(173), dif0(173), ar(173), meso(173), lebas(173) / 'CH3BR           ',0.0980,    2.0,      0.0,  59.5/
  DATA subname(174), dif0(174), ar(174), meso(174), lebas(174) / 'MB2             ',0.0804,    2.0,      0.0,  84.0/
  DATA subname(175), dif0(175), ar(175), meso(175), lebas(175) / 'MB2C            ',0.0720,    2.0,      0.0, 101.5/
  DATA subname(176), dif0(176), ar(176), meso(176), lebas(176) / 'MBC2            ',0.0739,    2.0,      0.0,  94.5/
  DATA subname(177), dif0(177), ar(177), meso(177), lebas(177) / 'MBC             ',0.0834,    2.0,      0.0,  77.0/
  DATA subname(178), dif0(178), ar(178), meso(178), lebas(178) / 'CLO             ',0.1288,    8.0,      0.0,  31.5/


CONTAINS


  FUNCTION DEPVVARS_INIT( ) RESULT ( SUCCESS )

    IMPLICIT NONE

    LOGICAL         :: SUCCESS
    INTEGER         :: ALLOCSTAT
    INTEGER         :: L
    CHARACTER( 96 ) :: XMSG
    REAL, PARAMETER :: TWOTHIRD = 2.0 / 3.0

    SUCCESS = .TRUE.

    !-------------------------------------------------------------------------------
    ! For M3DRY, set up core species, and include toxic and chlorine compounds.
    !-------------------------------------------------------------------------------

    DEPV_METHOD  = 'M3DRY           '

    DEPSPC(  1 ) = 'SO2             '
    DEPSPC(  2 ) = 'SULF            '
    DEPSPC(  3 ) = 'NO2             '
    DEPSPC(  4 ) = 'NO              '
    DEPSPC(  5 ) = 'O3              '
    DEPSPC(  6 ) = 'HNO3            '
    DEPSPC(  7 ) = 'H2O2            '
    DEPSPC(  8 ) = 'ALD             '
    DEPSPC(  9 ) = 'HCHO            '
    DEPSPC( 10 ) = 'OP              '
    DEPSPC( 11 ) = 'PAA             '
    DEPSPC( 12 ) = 'ORA             '
    DEPSPC( 13 ) = 'NH3             '
    DEPSPC( 14 ) = 'PAN             '
    DEPSPC( 15 ) = 'HONO            '
    DEPSPC( 16 ) = 'CO              '
    DEPSPC( 17 ) = 'METHANOL        '
    DEPSPC( 18 ) = 'N2O5            '
    DEPSPC( 19 ) = 'NO3             '
    DEPSPC( 20 ) = 'GEN_ALD         '
    DEPSPC( 21 ) = 'CL2             '
    DEPSPC( 22 ) = 'HOCL            '
    DEPSPC( 23 ) = 'HCL             '
    DEPSPC( 24 ) = 'FMCL            '
    DEPSPC( 25 ) = 'HG              '
    DEPSPC( 26 ) = 'HGIIGAS         '
    DEPSPC( 27 ) = 'TECDD_2378      '
    DEPSPC( 28 ) = 'PECDD_12378     '
    DEPSPC( 29 ) = 'HXCDD_123478    '
    DEPSPC( 30 ) = 'HXCDD_123678    '
    DEPSPC( 31 ) = 'HXCDD_123789    '
    DEPSPC( 32 ) = 'HPCDD_1234678   '
    DEPSPC( 33 ) = 'OTCDD           '
    DEPSPC( 34 ) = 'TECDF_2378      '
    DEPSPC( 35 ) = 'PECDF_12378     '
    DEPSPC( 36 ) = 'PECDF_23478     '
    DEPSPC( 37 ) = 'HXCDF_123478    '
    DEPSPC( 38 ) = 'HXCDF_123678    '
    DEPSPC( 39 ) = 'HXCDF_234678    '
    DEPSPC( 40 ) = 'HXCDF_123789    '
    DEPSPC( 41 ) = 'HPCDF_1234678   '
    DEPSPC( 42 ) = 'HPCDF_1234789   '
    DEPSPC( 43 ) = 'OTCDF           '
    DEPSPC( 44 ) = 'NAPHTHALENE     '
    DEPSPC( 45 ) = '1NITRONAPHTHA   '
    DEPSPC( 46 ) = '2NITRONAPHTHA   '
    DEPSPC( 47 ) = '14NAPHTHOQUIN   '
    DEPSPC( 48 ) = 'HEXMETH_DIIS    '
    DEPSPC( 49 ) = 'HYDRAZINE       '
    DEPSPC( 50 ) = 'MAL_ANHYDRIDE   '
    DEPSPC( 51 ) = 'TOLUENE_DIIS    '
    DEPSPC( 52 ) = 'TRIETHYLAMINE   '
    DEPSPC( 53 ) = 'NTR             '
    DEPSPC( 54 ) = 'NTRM            '
    DEPSPC( 55 ) = 'MPAN            '
    DEPSPC( 56 ) = 'PPN             '
    DEPSPC( 57 ) = 'ISPD            '
    DEPSPC( 58 ) = 'NTRDN           '
    DEPSPC( 59 ) = 'NTRALK          '
    DEPSPC( 60 ) = 'NTROH           '
    DEPSPC( 61 ) = 'NTRPX           '
    DEPSPC( 62 ) = 'PROPNN          '
    DEPSPC( 63 ) = 'CLNO2           '
    DEPSPC( 64 ) = 'ISOPNN          '
    DEPSPC( 65 ) = 'MTNO3           '
    DEPSPC( 66 ) = 'IEPOX           '
    DEPSPC( 67 ) = 'HACET           '
    DEPSPC( 68 ) = 'SVALK1          '
    DEPSPC( 69 ) = 'SVALK2          '
    DEPSPC( 70 ) = 'SVBNZ1          '
    DEPSPC( 71 ) = 'SVBNZ2          '
    DEPSPC( 72 ) = 'SVISO1          '
    DEPSPC( 73 ) = 'SVISO2          '
    DEPSPC( 74 ) = 'SVPAH1          '
    DEPSPC( 75 ) = 'SVPAH2          '
    DEPSPC( 76 ) = 'SVSQT           '
    DEPSPC( 77 ) = 'SVTOL1          '
    DEPSPC( 78 ) = 'SVTOL2          '
    DEPSPC( 79 ) = 'SVTRP1          '
    DEPSPC( 80 ) = 'SVTRP2          '
    DEPSPC( 81 ) = 'SVXYL1          '
    DEPSPC( 82 ) = 'SVXYL2          '
    DEPSPC( 83 ) = 'IO              '
    DEPSPC( 84 ) = 'OIO             '
    DEPSPC( 85 ) = 'I2O2            '
    DEPSPC( 86 ) = 'I2O3            '
    DEPSPC( 87 ) = 'I2O4            '
    DEPSPC( 88 ) = 'HI              '
    DEPSPC( 89 ) = 'HOI             '
    DEPSPC( 90 ) = 'INO             '
    DEPSPC( 91 ) = 'INO2            '
    DEPSPC( 92 ) = 'IONO2           '
    DEPSPC( 93 ) = 'BRO             '
    DEPSPC( 94 ) = 'HOBR            '
    DEPSPC( 95 ) = 'HBR             '
    DEPSPC( 96 ) = 'BRONO2          '
    DEPSPC( 97 ) = 'BRNO2           '
    DEPSPC( 98 ) = 'BRCL            '
    DEPSPC( 99 ) = 'DMS             '
    DEPSPC( 100) = 'MSA             '
    DEPSPC( 101) = 'METHANE         '
    DEPSPC( 102) = 'ACRYACID        '
    DEPSPC( 103) = 'CARBSULFIDE     '
    DEPSPC( 104) = 'ACETONITRILE    '
    DEPSPC( 105) = 'METH_NIT_PHEN   ' ! 6-methyl-2-nitrophenol aka 6-nitro-o-cresol
    DEPSPC( 106) = 'PCVOC           '
    DEPSPC( 107) = 'INTR            '
    DEPSPC( 108) = 'ISPX            '
    DEPSPC( 109) = 'ROOH            '
    DEPSPC( 110) = 'LVPCSOG         '
    DEPSPC( 111) = 'VIVPO1          '
    DEPSPC( 112) = 'VLVOO1          '
    DEPSPC( 113) = 'VLVOO2          '
    DEPSPC( 114) = 'VLVPO1          '
    DEPSPC( 115) = 'VSVOO1          '
    DEPSPC( 116) = 'VSVOO2          '
    DEPSPC( 117) = 'VSVOO3          '
    DEPSPC( 118) = 'VSVPO1          '
    DEPSPC( 119) = 'VSVPO2          '
    DEPSPC( 120) = 'VSVPO3          '
    DEPSPC( 121) = 'FACD            '
    DEPSPC( 122) = 'KET             '
    DEPSPC( 123) = 'ETH             '
    DEPSPC( 124) = 'PNA             '
    DEPSPC( 125) = 'GLY             '
    DEPSPC( 126) = 'GLYD            '
    DEPSPC( 127) = 'MGLY            '
    DEPSPC( 128) = 'ETHA            '
    DEPSPC( 129) = 'ETOH            '
    DEPSPC( 130) = 'PAR             '
    DEPSPC( 131) = 'ACET            '
    DEPSPC( 132) = 'PRPA            '
    DEPSPC( 133) = 'ETHY            '
    DEPSPC( 134) = 'OLE             '
    DEPSPC( 135) = 'IOLE            '
    DEPSPC( 136) = 'IEPOX           '
    DEPSPC( 137) = 'BENZ            '
    DEPSPC( 138) = 'CRES            '
    DEPSPC( 139) = 'TOL             '
    DEPSPC( 140) = 'XYLMN           '
    DEPSPC( 141) = 'NAPH            '
    DEPSPC( 142) = 'CAT1            '
    DEPSPC( 143) = 'SESQ            '
    DEPSPC( 144) = 'TERP            '
    DEPSPC( 145) = 'ISOP            '
    DEPSPC( 146) = 'OPEN            '
    DEPSPC( 147) = 'XOPN            '
    DEPSPC( 148) = 'SOAALK          '
    DEPSPC( 149) = 'BUTADIENE13     '
    DEPSPC( 150) = 'ACROLEIN        '
    DEPSPC( 151) = 'SVMT1           '
    DEPSPC( 152) = 'SVMT2           '
    DEPSPC( 153) = 'SVMT3           '
    DEPSPC( 154) = 'SVMT4           '
    DEPSPC( 155) = 'SVMT5           '
    DEPSPC( 156) = 'SVMT6           '
    DEPSPC( 157) = 'SVMT7           '
    DEPSPC( 158) = 'SVAVB1          '
    DEPSPC( 159) = 'SVAVB2          '
    DEPSPC( 160) = 'SVAVB3          '
    DEPSPC( 161) = 'SVAVB4          '
    DEPSPC( 162) = 'CLNO3           '
    DEPSPC( 163) = 'FMBR            '
    DEPSPC( 164) = 'I2              '
    DEPSPC( 165) = 'CH3I            '
    DEPSPC( 166) = 'ICL             '
    DEPSPC( 167) = 'IBR             '
    DEPSPC( 168) = 'MI2             '
    DEPSPC( 169) = 'MIB             '
    DEPSPC( 170) = 'MIC             '
    DEPSPC( 171) = 'BR2             '
    DEPSPC( 172) = 'MB3             '
    DEPSPC( 173) = 'CH3BR           '
    DEPSPC( 174) = 'MB2             '
    DEPSPC( 175) = 'MB2C            '
    DEPSPC( 176) = 'MBC2            '
    DEPSPC( 177) = 'MBC             '
    DEPSPC( 178) = 'CLO             '

    DO L = 1, N_SPC_M3DRY
       SCC_PR_23( L ) = 0.2 * ( DIF0( L ) * PR / KVIS ) ** TWOTHIRD
       DDIF0    ( L ) = DIF0( L ) / DWAT
       CSNOW    ( L ) = AR( L ) / ( A0 * RSNOW0 )
       CGND     ( L ) = AR( L ) / ( A0 * RG0 )
       CCUT     ( L ) = AR( L ) / ( A0 * RCUT0 )
       DW25     ( L ) = 13.26E-5 / ( 0.8904**1.14 * LEBAS( L )**0.589 )
    END DO

  END FUNCTION DEPVVARS_INIT

END MODULE DEPVVARS
