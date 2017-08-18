
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
C $Header$

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

      SUBROUTINE HRINIT
C***********************************************************************
C
C  FUNCTION: To initialize species tolerances, arrays, and indices
C
C  PRECONDITIONS: For the CB05TUCL family of mechanisms
C
C  RETURN VALUES: None
C
C  KEY SUBROUTINES/FUNCTIONS CALLED: None
C
C  REVISION HISTORY: Created by EBI solver program, May 26, 2017
C   21 Jun 10 J.Young:   convert for Namelist redesign
C   18 Jul 14 B.Hutzell: revised to use the RXNS_DATA and RXNS_FUNCTION
C                        modules instead of include files
C***********************************************************************

      USE UTILIO_DEFN
      USE HRDATA
      USE RXNS_DATA
      USE RXNS_FUNCTION
      use mem_grid, only: mza

      IMPLICIT NONE

C.....INCLUDES:
      Include 'CONST.EXT'          ! common constants

C.....ARGUMENTS: NONE

C.....PARAMETERS: NONE


C.....SAVED VARIABLES:
      CHARACTER( 16 ), SAVE  ::  PNAME = 'HRINIT'   ! Program name


C.....LOCAL VARIABLES:
      CHARACTER( 132 ) ::  MSG       ! Log message

      INTEGER N               ! Loop index
      INTEGER IOS             ! Allocate error status

C***********************************************************************

      N_SPEC = NUMB_MECH_SPC
      N_RXNS = NRXNS
      N_SPCL_RKS =   0

      allocate( FORWARD_CONV(N_SPEC), 
     &          REVERSE_CONV(N_SPEC),
     &          RTOL        (N_SPEC) )
 
      !$omp parallel private(ios,msg)
      ALLOCATE( RKI( NRXNS ),
     &          RKI_SAV( MZA, NRXNS ),
     &          RXRAT( NRXNS ),
     &          YC(   N_SPEC),
     &          YC0(  N_SPEC),
     &          YCP(  N_SPEC),
     &          PROD( N_SPEC),
     &          LOSS( N_SPEC),
     &          PNEG( N_SPEC),
     &          STAT = IOS )

       IF ( IOS .NE. 0 ) THEN
          MSG = 'Error allocating RKI, RXRAT, RTOL, YC, YC0, YCP, FORWARD_CONV, '
     &       // 'REVERSE_CONV, PROD, LOSS, or PNEG'
          CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
       END IF
      !$omp end parallel

       IF( NSPECIAL .GT. 0 )THEN
          !$omp parallel private(ios,msg)
          ALLOCATE( SYC( MZA, N_SPEC ), STAT = IOS )
          IF ( IOS .NE. 0 ) THEN
               MSG = 'Error allocating SYC'
               CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
          END IF
          !$omp end parallel
          ALLOCATE( IOLD2NEW( N_SPEC, 1 ), STAT = IOS )
          IF ( IOS .NE. 0 ) THEN
               MSG = 'Error allocating IOLD2NEW'
               CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
          END IF
          DO N = 1, NUMB_MECH_SPC
             IOLD2NEW( N, 1 ) = N
          END DO
       END IF

       IF( .NOT. MAP_CHEMISTRY_SPECIES() )THEN
          MSG = 'Detected above error(s) when mapping Chemistry species from CGRID species'
          CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT2 )
       END IF

       DO N = 1, NUMB_MECH_SPC
          FORWARD_CONV( N ) = 1.0E-3 * MWAIR / SPECIES_MOLWT( N )
          REVERSE_CONV( N ) = 1.0E+3 / MWAIR * SPECIES_MOLWT( N )
       END DO

! initialize all reaction rates to zero so unused array elements are defined

       !$omp parallel
       RXRAT   = 0.0D0
       RKI_SAV = 0.0D0

! initialize all production and loss rates to zero so unused array elements are defined

       PROD  = 0.0D0
       LOSS  = 0.0D0
       !$omp end parallel

       NING1 = 4
       NING2 = 4



c..Allocate and define ebi species

      N_EBISP  = 108
      ALLOCATE( EBISP( N_EBISP ), STAT = IOS ) 
      IF ( IOS .NE. 0 ) THEN
         MSG = 'Error allocating EBISP' 
         CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
      END IF

      EBISP(   1 ) = HNO3
      EBISP(   2 ) = H2O2
      EBISP(   3 ) = XO2
      EBISP(   4 ) = XO2N
      EBISP(   5 ) = NTR
      EBISP(   6 ) = ROOH
      EBISP(   7 ) = FORM
      EBISP(   8 ) = ALD2
      EBISP(   9 ) = ALDX
      EBISP(  10 ) = PAR
      EBISP(  11 ) = CO
      EBISP(  12 ) = MEO2
      EBISP(  13 ) = MEPX
      EBISP(  14 ) = MEOH
      EBISP(  15 ) = HCO3
      EBISP(  16 ) = FACD
      EBISP(  17 ) = PACD
      EBISP(  18 ) = AACD
      EBISP(  19 ) = CXO3
      EBISP(  20 ) = PANX
      EBISP(  21 ) = ROR
      EBISP(  22 ) = OLE
      EBISP(  23 ) = ETH
      EBISP(  24 ) = IOLE
      EBISP(  25 ) = TOL
      EBISP(  26 ) = CRES
      EBISP(  27 ) = TO2
      EBISP(  28 ) = TOLRO2
      EBISP(  29 ) = OPEN
      EBISP(  30 ) = MGLY
      EBISP(  31 ) = CRO
      EBISP(  32 ) = CAT1
      EBISP(  33 ) = CRON
      EBISP(  34 ) = CRNO
      EBISP(  35 ) = CRN2
      EBISP(  36 ) = CRPX
      EBISP(  37 ) = OPO3
      EBISP(  38 ) = CAO2
      EBISP(  39 ) = OPAN
      EBISP(  40 ) = XYLMN
      EBISP(  41 ) = XYLRO2
      EBISP(  42 ) = NAPH
      EBISP(  43 ) = PAHRO2
      EBISP(  44 ) = ISOP
      EBISP(  45 ) = ISPD
      EBISP(  46 ) = ISOPRXN
      EBISP(  47 ) = TERP
      EBISP(  48 ) = TRPRXN
      EBISP(  49 ) = SO2
      EBISP(  50 ) = SULF
      EBISP(  51 ) = SULRXN
      EBISP(  52 ) = ETOH
      EBISP(  53 ) = ETHA
      EBISP(  54 ) = CL2
      EBISP(  55 ) = CL
      EBISP(  56 ) = HOCL
      EBISP(  57 ) = CLO
      EBISP(  58 ) = FMCL
      EBISP(  59 ) = HCL
      EBISP(  60 ) = CLNO2
      EBISP(  61 ) = TOLNRXN
      EBISP(  62 ) = TOLHRXN
      EBISP(  63 ) = XYLNRXN
      EBISP(  64 ) = XYLHRXN
      EBISP(  65 ) = BENZENE
      EBISP(  66 ) = BENZRO2
      EBISP(  67 ) = BNZNRXN
      EBISP(  68 ) = BNZHRXN
      EBISP(  69 ) = SESQ
      EBISP(  70 ) = SESQRXN
      EBISP(  71 ) = PAHNRXN
      EBISP(  72 ) = PAHHRXN
      EBISP(  73 ) = SOAALK
      EBISP(  74 ) = ALKRXN
      EBISP(  75 ) = H2NO3PIJ
      EBISP(  76 ) = H2NO3PK
      EBISP(  77 ) = ACLI
      EBISP(  78 ) = ACLJ
      EBISP(  79 ) = ACLK
      EBISP(  80 ) = AXYL1J
      EBISP(  81 ) = AOLGAJ
      EBISP(  82 ) = AXYL2J
      EBISP(  83 ) = ATOL1J
      EBISP(  84 ) = ATOL2J
      EBISP(  85 ) = ABNZ1J
      EBISP(  86 ) = ABNZ2J
      EBISP(  87 ) = ATRP1J
      EBISP(  88 ) = AOLGBJ
      EBISP(  89 ) = ATRP2J
      EBISP(  90 ) = AISO1J
      EBISP(  91 ) = AISO2J
      EBISP(  92 ) = ASQTJ
      EBISP(  93 ) = APAH1J
      EBISP(  94 ) = APAH2J
      EBISP(  95 ) = AALK1J
      EBISP(  96 ) = AALK2J
      EBISP(  97 ) = PCVOC
      EBISP(  98 ) = PCSOARXN
      EBISP(  99 ) = VLVPO1
      EBISP( 100 ) = VSVPO1
      EBISP( 101 ) = VSVPO2
      EBISP( 102 ) = VSVPO3
      EBISP( 103 ) = VIVPO1
      EBISP( 104 ) = VLVOO1
      EBISP( 105 ) = VLVOO2
      EBISP( 106 ) = VSVOO2
      EBISP( 107 ) = VSVOO3
      EBISP( 108 ) = VSVOO1


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Set species tolerances
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RTOL( NO2      ) = 1.0D-03
      RTOL( NO       ) = 1.0D-03
      RTOL( O        ) = 1.0D+00
      RTOL( O3       ) = 1.0D-03
      RTOL( NO3      ) = 1.0D-03
      RTOL( O1D      ) = 1.0D+00
      RTOL( OH       ) = 1.0D-03
      RTOL( HO2      ) = 1.0D-03
      RTOL( N2O5     ) = 1.0D-03
      RTOL( HNO3     ) = 1.0D-03
      RTOL( HONO     ) = 1.0D-03
      RTOL( PNA      ) = 1.0D-03
      RTOL( H2O2     ) = 1.0D-03
      RTOL( XO2      ) = 1.0D-03
      RTOL( XO2N     ) = 1.0D-03
      RTOL( NTR      ) = 1.0D+00
      RTOL( ROOH     ) = 1.0D-03
      RTOL( FORM     ) = 1.0D-03
      RTOL( ALD2     ) = 1.0D-03
      RTOL( ALDX     ) = 1.0D-03
      RTOL( PAR      ) = 1.0D-03
      RTOL( CO       ) = 1.0D-03
      RTOL( MEO2     ) = 1.0D-03
      RTOL( MEPX     ) = 1.0D-03
      RTOL( MEOH     ) = 1.0D-03
      RTOL( HCO3     ) = 1.0D+00
      RTOL( FACD     ) = 1.0D-03
      RTOL( C2O3     ) = 1.0D-03
      RTOL( PAN      ) = 1.0D-03
      RTOL( PACD     ) = 1.0D-03
      RTOL( AACD     ) = 1.0D-03
      RTOL( CXO3     ) = 1.0D-03
      RTOL( PANX     ) = 1.0D-03
      RTOL( ROR      ) = 1.0D-03
      RTOL( OLE      ) = 1.0D-03
      RTOL( ETH      ) = 1.0D-03
      RTOL( IOLE     ) = 1.0D-03
      RTOL( TOL      ) = 1.0D-03
      RTOL( CRES     ) = 1.0D-03
      RTOL( TO2      ) = 1.0D-03
      RTOL( TOLRO2   ) = 1.0D-03
      RTOL( OPEN     ) = 1.0D-03
      RTOL( MGLY     ) = 1.0D-03
      RTOL( CRO      ) = 1.0D-03
      RTOL( CAT1     ) = 1.0D-03
      RTOL( CRON     ) = 1.0D-03
      RTOL( CRNO     ) = 1.0D+00
      RTOL( CRN2     ) = 1.0D+00
      RTOL( CRPX     ) = 1.0D-01
      RTOL( OPO3     ) = 1.0D-03
      RTOL( CAO2     ) = 1.0D-03
      RTOL( OPAN     ) = 1.0D-03
      RTOL( XYLMN    ) = 1.0D-03
      RTOL( XYLRO2   ) = 1.0D-03
      RTOL( NAPH     ) = 1.0D-03
      RTOL( PAHRO2   ) = 1.0D+00
      RTOL( ISOP     ) = 1.0D-03
      RTOL( ISPD     ) = 1.0D-03
      RTOL( ISOPRXN  ) = 1.0D+00
      RTOL( TERP     ) = 1.0D-03
      RTOL( TRPRXN   ) = 1.0D+00
      RTOL( SO2      ) = 1.0D-03
      RTOL( SULF     ) = 1.0D+00
      RTOL( SULRXN   ) = 1.0D+00
      RTOL( ETOH     ) = 1.0D-03
      RTOL( ETHA     ) = 1.0D-03
      RTOL( CL2      ) = 1.0D+00
      RTOL( CL       ) = 1.0D+00
      RTOL( HOCL     ) = 1.0D+00
      RTOL( CLO      ) = 1.0D+00
      RTOL( FMCL     ) = 1.0D+00
      RTOL( HCL      ) = 1.0D+00
      RTOL( CLNO2    ) = 1.0D-01
      RTOL( TOLNRXN  ) = 1.0D+00
      RTOL( TOLHRXN  ) = 1.0D+00
      RTOL( XYLNRXN  ) = 1.0D+00
      RTOL( XYLHRXN  ) = 1.0D+00
      RTOL( BENZENE  ) = 1.0D-03
      RTOL( BENZRO2  ) = 1.0D-03
      RTOL( BNZNRXN  ) = 1.0D+00
      RTOL( BNZHRXN  ) = 1.0D+00
      RTOL( SESQ     ) = 1.0D-03
      RTOL( SESQRXN  ) = 1.0D+00
      RTOL( PAHNRXN  ) = 1.0D+00
      RTOL( PAHHRXN  ) = 1.0D+00
      RTOL( SOAALK   ) = 1.0D-03
      RTOL( ALKRXN   ) = 1.0D+00
      RTOL( H2NO3PIJ ) = 1.0D-03
      RTOL( H2NO3PK  ) = 1.0D-03
      RTOL( ACLI     ) = 1.0D-03
      RTOL( ACLJ     ) = 1.0D-03
      RTOL( ACLK     ) = 1.0D-03
      RTOL( AXYL1J   ) = 1.0D-03
      RTOL( AOLGAJ   ) = 1.0D-03
      RTOL( AXYL2J   ) = 1.0D-03
      RTOL( ATOL1J   ) = 1.0D-03
      RTOL( ATOL2J   ) = 1.0D-03
      RTOL( ABNZ1J   ) = 1.0D-03
      RTOL( ABNZ2J   ) = 1.0D-03
      RTOL( ATRP1J   ) = 1.0D-03
      RTOL( AOLGBJ   ) = 1.0D-03
      RTOL( ATRP2J   ) = 1.0D-03
      RTOL( AISO1J   ) = 1.0D-03
      RTOL( AISO2J   ) = 1.0D-03
      RTOL( ASQTJ    ) = 1.0D-03
      RTOL( APAH1J   ) = 1.0D-03
      RTOL( APAH2J   ) = 1.0D-03
      RTOL( AALK1J   ) = 1.0D-03
      RTOL( AALK2J   ) = 1.0D-03
      RTOL( PCVOC    ) = 1.0D-03
      RTOL( PCSOARXN ) = 1.0D-03
      RTOL( VLVPO1   ) = 1.0D-03
      RTOL( VSVPO1   ) = 1.0D-03
      RTOL( VSVPO2   ) = 1.0D-03
      RTOL( VSVPO3   ) = 1.0D-03
      RTOL( VIVPO1   ) = 1.0D-03
      RTOL( VLVOO1   ) = 1.0D-03
      RTOL( VLVOO2   ) = 1.0D-03
      RTOL( VSVOO2   ) = 1.0D-03
      RTOL( VSVOO3   ) = 1.0D-03
      RTOL( VSVOO1   ) = 1.0D-03


      RETURN

      END
