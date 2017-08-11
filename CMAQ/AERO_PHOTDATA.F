
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
      MODULE AERO_PHOTDATA

      IMPLICIT NONE

!     REAL, ALLOCATABLE :: SHELL_MOMENT3( :,: )  ! modal shell third moment conc.
!     REAL, ALLOCATABLE :: CORE_MOMENT3 ( :,: )  ! model core third moment conc.

!     REAL, ALLOCATABLE :: BLK_AE_VOL  ( :,: )     ! aerosol modal volumes [ m**3/m**3 ]
!     REAL, ALLOCATABLE :: BLK_AE_LSG  ( :,: )     ! aerosol modal log of geo. stand. dev
!     REAL, ALLOCATABLE :: BLK_AE_DGN  ( :,: )     ! aerosol modal mean geometeric diam. [ m ]
!     REAL, ALLOCATABLE :: AE_DGN_CORE ( :,: )   ! inner core modal mean geometeric diam. [ m ]
!     REAL, ALLOCATABLE :: AE_DGN_SHELL( :,: )   ! modal mean geometeric diam. [ m ]
!     REAL, ALLOCATABLE :: BLK_AE_CONC ( :,:,: ) ! aerosol species concentrations [ ug/m**3]

!     REAL, ALLOCATABLE :: AE_NR_SHELL ( :,: ) ! real part of refractive index for outer shell
!     REAL, ALLOCATABLE :: AE_NI_SHELL ( :,: ) ! imaginary part of refractive index for outer shell
!     REAL, ALLOCATABLE :: AE_NR_CORE  ( :,: ) ! real part of refractive index for inner core
!     REAL, ALLOCATABLE :: AE_NI_CORE  ( :,: ) ! imaginary part of refractive index for inner core
!     REAL, ALLOCATABLE :: BLK_AE_NR   ( :,:,: ) ! mean aerosol real part of refractive index
!     REAL, ALLOCATABLE :: BLK_AE_NI   ( :,:,: ) ! mean aerosol imaginary part of refractive index

      REAL, ALLOCATABLE :: AERO_ASYM_FAC ( :,: ) ! aerosol modal averaged asymmetry factor
      REAL, ALLOCATABLE :: AERO_EXTI_COEF( :,: ) ! aerosol modal averaged extinction coeff.
      REAL, ALLOCATABLE :: AERO_SCAT_COEF( :,: ) ! aerosol modal averaged scattering coeff.

      !$omp threadprivate( AERO_ASYM_FAC, AERO_EXTI_COEF, AERO_SCAT_COEF)

      real, allocatable :: ri_r( :,:,: )
      real, allocatable :: ri_i( :,:,: )

      LOGICAL           :: CORE_SHELL = .false.
      LOGICAL           :: MIE_CALC   = .false.

      PUBLIC            :: AERO_ASYM_FAC, AERO_EXTI_COEF, AERO_SCAT_COEF,
     &                     INIT_AERO_DATA, GET_AERO_DATA
      PRIVATE

      INTEGER :: LOGDEV                   ! unit number for the log file
      INTEGER, ALLOCATABLE :: REFRACT_INDX_MAP( : )  ! map array for refactive index

      REAL,    ALLOCATABLE :: VOL_MINS( : )          ! minmum volume permode
!     REAL,    ALLOCATABLE :: M3_FACTOR( : )         ! convert mass to portion of 3rd moment
      LOGICAL, ALLOCATABLE :: IS_CORE_SPCS( :,: )    ! flag denoting species is in aerosol core

      INTEGER, PARAMETER   :: NUMB_COR_SPCS = 3      ! number species in core

      REAL, allocatable :: INV_LAMBDA(:) ! reciprocal of wavelength [ 1/m ]
      REAL, allocatable ::  LAMDA_UM (:) ! wavelength  [ um ]

C *** Species in aerosol core

      CHARACTER( 16 ), PARAMETER :: CORE_SPCS( NUMB_COR_SPCS ) =
     &                              (/ 'AECI', 'AECJ', 'AECK' /)
      
      CONTAINS

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE INIT_AERO_DATA( nlays )

C Allocate and define arrays used to calculate aerosol optical properties

         USE UTILIO_DEFN
         USE AERO_DATA           ! aero variable data
         USE AEROMET_DATA        ! constants and met data
         USE CSQY_DATA,   ONLY: NWL, NUM_REFRACTIVE, REFRACTIVE_INDEX,
     &                          WAVELENGTH => EFFWL_REF
         use oname_coms,   only: nl

         IMPLICIT NONE

C Arguments:

         integer, intent(in) :: nlays

C Local:
         CHARACTER(  32 ) :: PNAME             = 'INIT_AERO_DATA'
!        CHARACTER(  16 ) :: CORE_SHELL_OPTICS = 'CORESHELL_OPTICS'
!        CHARACTER(  16 ) :: OPTICS_MIE_CALC   = 'MIE_OPTICS'
         CHARACTER( 120 ) :: XMSG 

         INTEGER          :: ALLOCSTAT
!        INTEGER          :: ESTAT          ! status from environment var check
         INTEGER          :: JDATE, JTIME
         INTEGER          :: M, N, SPC, V   ! loop counters
      
         LOGICAL, SAVE    :: INITIALIZED  = .FALSE.

#ifdef verbose_phot
         character( 26 ), allocatable :: lambda_list( : )
#endif

         IF ( INITIALIZED ) THEN
            RETURN
         END IF
          
         LOGDEV = INIT3()
         JDATE  = 0
         JTIME  = 0

!...Get flag to use core-shell mixing model for aerosol optical properties

         core_shell = nl%core_shell_photolysis
!         CORE_SHELL = .false.         ! default
!         XMSG = 'Use core-shell mixing model for aerosol optical properties'
!         CORE_SHELL = ENVYN( CORE_SHELL_OPTICS, XMSG, CORE_SHELL, ESTAT )
!         IF ( ESTAT .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) XMSG
!         IF ( ESTAT .EQ. 1 ) THEN
!            XMSG = 'Environment variable improperly formatted'
!            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!         ELSE IF ( ESTAT .EQ. -1 ) THEN
!            XMSG =
!     &          'Environment variable set, but empty ... Using default:'
!            WRITE( LOGDEV, '(5X, A, L5)' ) XMSG, CORE_SHELL
!         ELSE IF ( ESTAT .EQ. -2 ) THEN
!            XMSG = 'Environment variable not set ... Using default:'
!            WRITE( LOGDEV, '(5X, A, L5)' ) XMSG, CORE_SHELL
!         END IF

!...Get flag to use fast optics for volume mixing model for aerosol optical properties

         mie_calc = nl%mie_calc_photolysis
!         MIE_CALC = .FALSE.         ! default
!         XMSG = 'Use Mie theory for aerosol optical properties of Internal mixing model'
!         MIE_CALC = ENVYN( OPTICS_MIE_CALC, XMSG, MIE_CALC, ESTAT )
!         IF ( ESTAT .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) XMSG
!         IF ( ESTAT .EQ. 1 ) THEN
!            XMSG = 'Environment variable improperly formatted'
!            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
!         ELSE IF ( ESTAT .EQ. -1 ) THEN
!            XMSG =
!     &          'Environment variable set, but empty ... Using default:'
!            WRITE( LOGDEV, '(5X, A, L5)' ) XMSG, MIE_CALC
!         ELSE IF ( ESTAT .EQ. -2 ) THEN
!            XMSG = 'Environment variable not set ... Using default:'
!            WRITE( LOGDEV, '(5X, A, L5)' ) XMSG, MIE_CALC
!         END IF

C...Allocate needed arrays

         ALLOCATE ( VOL_MINS( N_MODE ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating VOL_MINS'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

!         ALLOCATE ( BLK_AE_CONC( N_AEROSPC,N_MODE,NLAYS ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating BLK_AE_CONC'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( SHELL_MOMENT3( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating SHELL_MOMENT3'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( BLK_AE_VOL( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating BLK_AE_VOL'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( BLK_AE_LSG( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating BLK_AE_LSG'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( BLK_AE_DGN( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating BLK_AE_DGN'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( AE_DGN_SHELL( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating AE_DGN_SHELL'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( AE_NR_SHELL( nwl,N_MODE ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating  AE_NR_SHELL'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

!         ALLOCATE ( AE_NI_SHELL( nwl,N_MODE ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating  AE_NI_SHELL'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF

         !$omp parallel private(allocstat, xmsg)
         ALLOCATE ( AERO_ASYM_FAC( NLAYS,NWL ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating AERO_ASYM_FAC'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE ( AERO_EXTI_COEF( NLAYS,NWL ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating AERO_EXTI_FAC'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ALLOCATE ( AERO_SCAT_COEF( NLAYS,NWL ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating AERO_EXTI_FAC'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         AERO_SCAT_COEF = 0.0
         AERO_EXTI_COEF = 0.0
         AERO_ASYM_FAC  = 0.0
         !$omp end parallel

         CALL MAP_AERO()
                  
!         ALLOCATE( M3_FACTOR( N_AEROSPC ), STAT = ALLOCSTAT )
!         IF ( ALLOCSTAT .NE. 0 ) THEN
!            XMSG = 'Failure allocating  M3_FACTOR'
!            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!         END IF
         
         ALLOCATE( REFRACT_INDX_MAP( N_AEROSPC ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating REFRACT_INDX_MAP'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         
!        IF( CORE_SHELL )THEN
            ALLOCATE( IS_CORE_SPCS( N_AEROSPC,N_MODE ), STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'Failure allocating IS_CORE_SPCS'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         
            IS_CORE_SPCS = .FALSE.

!            ALLOCATE ( CORE_MOMENT3( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!            IF ( ALLOCSTAT .NE. 0 ) THEN
!               XMSG = 'Failure allocating CORE_MOMENT3'
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
      
!            ALLOCATE ( AE_DGN_CORE( N_MODE,NLAYS ), STAT = ALLOCSTAT )
!            IF ( ALLOCSTAT .NE. 0 ) THEN
!               XMSG = 'Failure allocating AE_DGN_CORE'
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
      
!            ALLOCATE ( AE_NR_CORE( nwl,N_MODE ), STAT = ALLOCSTAT )
!            IF ( ALLOCSTAT .NE. 0 ) THEN
!               XMSG = 'Failure allocating AE_NR_CORE'
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
      
!            ALLOCATE ( AE_NI_CORE( nwl,N_MODE ), STAT = ALLOCSTAT )
!            IF ( ALLOCSTAT .NE. 0 ) THEN
!               XMSG = 'Failure allocating AE_NI_CORE'
!               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
!            END IF
!        END IF
C ***    Create map for refactive index for each species

#ifdef verbose_phot
         allocate( lambda_list( nwl ) )
         do v = 1, nwl
            write( lambda_list( v ), 5000)int(wavelength( v )), 
     &      int(wavelength( v ))
         end do
         write( logdev, 5001)(lambda_list( v ), v = 1, nwl)
#endif 

         IF ( MIE_CALC .AND. .NOT. CORE_SHELL ) THEN
            CALL M3MESG('For optics of internally mixed aerosol, Mie calculation used')
            CALL M3MESG('Core-shell mixing model option is set to off')
         ELSE IF ( MIE_CALC .AND. CORE_SHELL ) THEN
            CALL M3MESG('For AERO mixing models, core-shell model is possible')
            CALL M3MESG('For optics of aerosols with no core, Mie calculation used')
         ELSE IF ( .NOT. MIE_CALC .AND. CORE_SHELL ) THEN
            CALL M3MESG('For AERO mixing models, core-shell model is possible')
            CALL M3MESG('For optics of aerosols with no core, FASTER_OPTICS used')
         ELSE
            CALL M3MESG('Using internal mixing model for aerosol and FASTER_OPTICS')
         END IF

         FORALL ( SPC = 1:N_MODE ) VOL_MINS( SPC ) = PI6 * AEROMODE( SPC )%MIN_M3CONC

         DO SPC = 1, N_AEROSPC

!           M3_FACTOR( SPC ) = 1.0E-9 * F6DPI / AEROSPC( SPC )%DENSITY

            DO M = 1, NUM_REFRACTIVE
               IF ( AEROSPC( SPC )%OPTIC_SURR .EQ. REFRACTIVE_INDEX( M )%NAME ) THEN
                  REFRACT_INDX_MAP( SPC ) = M
                  DO N = 1, N_MODE

                     IF ( CORE_SHELL ) THEN
                        DO V = 1, NUMB_COR_SPCS
                           IF ( TRIM(AEROSPC( SPC )%NAME( N )) .EQ. TRIM(CORE_SPCS( V )) ) THEN
                               IS_CORE_SPCS( SPC, N ) = .TRUE.
                           END IF
                        END DO
                     END IF
#ifdef verbose_phot
                     if ( aero_missing(spc,n) ) cycle
                     write( logdev, 5002)trim( aerospc( spc )%name( n ) ), 
     &               trim( refractive_index( m )%name ),is_core_spcs( spc,n ),
     &               (refractive_index( m )%real_part( n,v ), 
     &                refractive_index( m )%imag_part( n,v ), v = 1, nwl)
#endif                  
                  END DO 
                  EXIT
               END IF
            END DO
            IF ( REFRACT_INDX_MAP( SPC ) .LT. 1 ) THEN
               DO M = 1, N_MODE
                  IF ( aero_missing(spc,n) ) CYCLE
                  XMSG = 'Species ' // TRIM( AEROSPC( SPC )%NAME( M ) )
     &                 // ' does not use defined refractive index.'
                  CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT3 )
               END DO
            END IF
         END DO
         
         allocate(lamda_um  (nwl))
         allocate(inv_lambda(nwl))

         do v = 1, nwl
            LAMDA_UM  (v) = 1.0E-3 * WAVELENGTH( V ) ! convert to micro-meters
            INV_LAMBDA(v) = 1.0E+9 / WAVELENGTH( V )
         enddo

         allocate(ri_r(nwl,NUM_REFRACTIVE,n_mode))
         allocate(ri_i(nwl,NUM_REFRACTIVE,n_mode))
         do m = 1, n_mode
            do n = 1, NUM_REFRACTIVE
               do v = 1, nwl
                  ri_r(v,n,m) = REFRACTIVE_INDEX( n )%REAL_PART( m,v )
                  ri_i(v,n,m) = REFRACTIVE_INDEX( n )%IMAG_PART( m,v )
               enddo
            enddo
         enddo

      RETURN

5000  FORMAT(2X,'NR_',I3.3,7X,'NI_',I3.3,5X)
5002  FORMAT(2(A16,1X),L12,1X,200(ES12.4,1X))                    
5001  FORMAT('     AE_SPC     ',1X,'  REFRACT_INDX  ',1X,'PART_OF_CORE',1X,200(A26))

      END SUBROUTINE INIT_AERO_DATA

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE GET_AERO_DATA ( iw, NLAYS, koff )

C-----------------------------------------------------------------------
C  FUNCTION:  This subroutine calculates the volume, the natural logs of 
C     geometric standard of all 3 modes. It assumes that aerosols have inner 
C     core with an outer shell to compute a core and total diameter per mode.  
C     The code also calculates the average modal aerosol refractive indices 
C     for the inner core and outer shell.
C
C     Internally mixed wet aerosols are assumed.
C
C     NOTE
C       2nd moment concentrations (M2) are passed into GETPAR in the
C       CBLK array.  The GETPAR calculations implicitly assume
C       that the input value of M2 is up to date.  
C       If, for example, the input M2 value was calculated
C       for a "dry" aerosol and the Wet_Moments_FLAG is .TRUE., GETPAR would
C       incorrectly adjust the M2 concentrations!
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  EXTRACT_AERO, EXTRACT_SOA, GETPAR
C                                     AERO_OPTICAL2, AERO_OPTICAL_CS  
C  REVISION HISTORY
C     06/13/13 B.Hutzell - initial version created for aero_photdata.F 
C                          from version 5.01
C     05/2016 H. Pye and B. Murphy - updated treatment of aerosol
C                          moments
C
C  REFERENCES:
C     Horvath, H., Size segrated light absorption coefficient for
C     the atmospheric aerosol, Atmos. Environ., Vol. 29, No. 8,
C     pp 875-883, 1995.
C
C     Beyer, K.D. A.R. Ravishankara, and E. R. Lovejoy, Meanurements
C     of H2SO4/H2O and H2SO4/HNO3/H2O solutions, J. Geophys. Res. Vol.
C     101, D9, pp 14519-14524, 1996.
C
C-----------------------------------------------------------------------

      USE RXNS_DATA
      USE CGRID_SPCS          ! CGRID mechanism species
      USE UTILIO_DEFN

      USE AERO_DATA           ! aero variable data
      USE AEROMET_DATA        ! constants and met data
      USE CSQY_DATA,  ONLY : NWL, NUM_REFRACTIVE, WAVELENGTH => EFFWL_REF

      USE RRTMG_AERO_OPTICAL_UTIL_MODULE, ONLY : AERO_OPTICAL2, AERO_OPTICAL_CS
      use cgrid_defn, only: cgrid

      IMPLICIT NONE

C Arguments:

      INTEGER, INTENT( IN ) :: iw               ! specified column index
      INTEGER, INTENT( IN ) :: NLAYS            ! # of vertical layers
      INTEGER, INTENT( IN ) :: koff             ! virtical grid offset

C Parameters:

      REAL, PARAMETER :: NEGLECT_FRAC = 1.0e-6      ! skip if below this volume fraction
      REAL, PARAMETER :: ONE_THIRD    = 1.0 / 3.0
      REAL, PARAMETER :: ONE_OVER_PI  = 1.0 / PI
      real, parameter :: lsg_max      = log(2.05)
      real, parameter :: lsg_min      = log(1.05)

C Local variables:

      REAL :: FACTOR
      REAL :: FRACTION_CORE
      REAL :: MOMENT3_PORTION

      REAL DGN_CORE       ! geometric mean diameters of aerosol core
      REAL VOL            ! Total Modal aerosol volume
      REAL VOL_CORE       ! Modal aerosol volume of aerosol core
      REAL DGN_SHELL      ! geometric mean diameters of aerosol core
      REAL VOL_SHELL      ! Modal aerosol volume of aerosol core
      REAL XXLSG          ! geometric standard deviation
      REAL LSG            ! log of geo. stand. dev
      REAL MIE_PARAMETER  ! Modal Mie parameter
      REAL CORE_MOMENT3
      REAL SHELL_MOMENT3
      real expfac2
      real expfac3

C***  variables for calculating modal averaged properties per layer and lambda

      COMPLEX NR_SHELL
      COMPLEX NR_CORE

      REAL BEXT_MODE   ! total aerosol extinction coefficient [ 1/m ]
      REAL BSCAT_MODE  ! total aerosol scattering coefficient [ 1/m ]
      REAL GBAR_MODE   ! total aerosol asymmetry factor
      REAL BEXT(nwl)   ! total aerosol extinction coefficient [ 1/m ]
      REAL VFAC        ! unit correction factors
      REAL BSCAT(nwl)  ! total aerosol scattering coefficient [ 1/m ]
      REAL G_BAR(nwl)  ! total aerosol asymmetry factor

      real AE_NR_SHELL(nwl)
      real AE_NI_SHELL(nwl)
      real AE_NR_CORE(nwl)
      real AE_NI_CORE(nwl)

      INTEGER L, V, M, k            ! loop counters
      INTEGER SPC                   ! species loop counter
      INTEGER MODE                  ! aerosol mode loop counter
      LOGICAL SUCCESS
      
      real, parameter :: asym_max = 0.990
      real, parameter :: exti_min = 1.e-9
      real, parameter :: scat_min = 1.e-9

C-----------------------------------------------------------------------

      ! *** initialize internal and external dimensions of aerosols

      AERO_SCAT_COEF = scat_min
      AERO_EXTI_COEF = exti_min
      AERO_ASYM_FAC  = asym_max

      IF ( N_AE_SPC .LE. 0 ) RETURN

      SUCCESS       = .TRUE.

      LAY_LOOP: DO L = 1, NLAYS

         k = l + koff

         ! *** extract grid cell concentrations of aero species from CGRID into aerospc_conc 
         ! *** in aero_data module
         !     Also converts dry surface area to wet second moment
         CALL EXTRACT_AERO ( CGRID( k,IW,: ), .TRUE. )

         ! *** Calculate geometric mean diameters and standard deviations of "wet" size distribution
         CALL GETPAR ( .FALSE. )

         G_BAR(1:nwl) = 0.0
         BEXT (1:nwl) = 0.0
         BSCAT(1:nwl) = 0.0

         MODE_LOOP: DO MODE = 1, N_MODE

           if (moment3_conc(mode) < 1.e5 * aeromode( mode )%min_m3conc ) cycle mode_loop
           if (moment2_conc(mode) < 1.e5 * aeromode( mode )%min_m2conc ) cycle mode_loop
           if (moment0_conc(mode) < 1.e5 * aeromode( mode )%min_numconc) cycle mode_loop

            AE_NR_SHELL(1:NWL) = 0.0
            AE_NI_SHELL(1:NWL) = 0.0
            AE_NR_CORE (1:NWL) = 0.0
            AE_NI_CORE (1:NWL) = 0.0

            CORE_MOMENT3  = 0.0
            SHELL_MOMENT3 = 0.0

            SPC_LOOP: DO SPC = 1, N_AEROSPC

               ! neglect tracer and undefined species
               IF ( AEROSPC( SPC )%TRACER ) CYCLE spc_loop
               if ( AERO_MISSING( SPC,MODE ) ) CYCLE spc_loop

               ! neglect species with small concentration
               if ( AEROSPC_CONC( SPC,MODE ) < 1.e-23 ) cycle spc_loop

               MOMENT3_PORTION = AEROSPC_M3FAC( SPC ) * AEROSPC_CONC( SPC,MODE )

               M = REFRACT_INDX_MAP( SPC )

               IF( IS_CORE_SPCS( SPC,MODE ) ) THEN  

                  CORE_MOMENT3 = CORE_MOMENT3 + MOMENT3_PORTION
                  DO V = 1, NWL
                     AE_NR_CORE(v) = AE_NR_CORE(v) + MOMENT3_PORTION * ri_r(v,m,mode)
                     AE_NI_CORE(v) = AE_NI_CORE(v) + MOMENT3_PORTION * ri_i(v,m,mode)
                  ENDDO

               ELSE

                  SHELL_MOMENT3 = SHELL_MOMENT3 + MOMENT3_PORTION
                  DO V = 1, NWL
                     AE_NR_SHELL(v) = AE_NR_SHELL(v) + MOMENT3_PORTION * ri_r(v,m,mode)
                     AE_NI_SHELL(v) = AE_NI_SHELL(v) + MOMENT3_PORTION * ri_i(v,m,mode)
                  ENDDO

               END IF 

            END DO SPC_LOOP

            moment3_portion = CORE_MOMENT3 + SHELL_MOMENT3

            if (moment3_portion < 1.e5 * aeromode( mode )%min_m3conc) cycle mode_loop

            vol = PI6 * moment3_portion
            LSG = AEROMODE_LNSG( MODE )

            dgn_shell = AEROMODE_DIAM( MODE )
            dgn_core  = 0.0

            IF (CORE_SHELL) THEN
               FRACTION_CORE = CORE_MOMENT3 / MOMENT3_PORTION
               IF ( FRACTION_CORE > NEGLECT_FRAC ) THEN
                  FACTOR   = min(FRACTION_CORE ** ONE_THIRD, 0.999)
                  DGN_CORE = FACTOR * AEROMODE_DIAM( MODE )
               ELSE
                  DGN_CORE = 0.0
               ENDIF
            ENDIF

            do v = 1, nwl

               if (core_shell) then
                  AE_NR_CORE(v) = AE_NR_CORE(v) / max(CORE_MOMENT3, 1.e-36)
                  AE_NI_CORE(v) = AE_NI_CORE(v) / max(CORE_MOMENT3, 1.e-36)
               END IF
               
               AE_NR_SHELL(v) = AE_NR_SHELL(v) / max(SHELL_MOMENT3, 1.e-36)
               AE_NI_SHELL(v) = AE_NI_SHELL(v) / max(SHELL_MOMENT3, 1.e-36)

            enddo
               
            if (CORE_SHELL .and. DGN_CORE > 1.e-20) then

               VOL_CORE  = PI6 * CORE_MOMENT3
               VOL_SHELL = PI6 * SHELL_MOMENT3

               XXLSG = EXP( lsg )

               do v = 1, nwl

                  success = .true.

                  NR_CORE  = CMPLX( AE_NR_CORE ( v ), AE_NI_CORE ( v ) )
                  NR_SHELL = CMPLX( AE_NR_SHELL( v ), AE_NI_SHELL( v ) )

                  IF ( AE_NR_CORE( v ) .LT. 1.e-30 ) cycle
                  IF ( AE_NI_CORE( v ) .LT. 1.e-30 ) cycle

                  IF ( AE_NR_SHELL( v ) .LT. 1.e-30 ) cycle
                  IF ( AE_NI_SHELL( v ) .LT. 1.e-30 ) cycle

                  ! *** Core-Shell case for optics
                  CALL AERO_OPTICAL_CS( LAMDA_UM(v), NR_CORE, NR_SHELL,  
     &                                   VOL_CORE, VOL_SHELL, DGN_CORE,      
     &                                   DGN_SHELL, XXLSG,
     &                                   BEXT_MODE, BSCAT_MODE, GBAR_MODE, SUCCESS )

                  !***sum to get total extinction and scattering
                  !***  and contribution to the overall asymmetry factor

                  BEXT (v) = BEXT (v) + BEXT_MODE
                  BSCAT(v) = BSCAT(v) + BSCAT_MODE
                  G_BAR(v) = G_BAR(v) + BSCAT_MODE * GBAR_MODE

               enddo

            ELSE   ! Not core shell

               if (mie_calc) then
                  XXLSG = EXP( lsg )
               else
                  lsg = min(lsg, lsg_max)
                  lsg = max(lsg, lsg_min)
                  xxlsg = exp( lsg )
                  MIE_PARAMETER = dgn_shell * exp( 3.0 * lsg * lsg )
                  EXPFAC2 = EXP( 2.0 * lsg * lsg )
                  EXPFAC3 = EXP( 4.5 * lsg * lsg )
               endif

               do v = 1, nwl

                  IF ( AE_NR_SHELL( V ) .LT. 1.e-30 ) cycle
                  IF ( AE_NI_SHELL( V ) .LT. 1.e-30 ) cycle

                  IF ( MIE_CALC ) THEN

                     !*** FSB Internal Volume mixture case for optics 
                     NR_SHELL = CMPLX( AE_NR_SHELL( v ), AE_NI_SHELL( v ) )
                     CALL AERO_OPTICAL2( LAMDA_UM(v), NR_SHELL, VOL, 
     &                                   DGN_SHELL, XXLSG,
     &                                   BEXT_MODE, BSCAT_MODE, GBAR_MODE, SUCCESS )
                  ELSE 

                     ! use old fast optics method
                     CALL FAST_OPTICS( AE_NR_SHELL(v), AE_NI_SHELL(v), 
     &                                 WAVELENGTH(V), XXLSG, lsg, mie_parameter, 
     &                                 expfac2, expfac3, 
     &                                 BEXT_MODE, BSCAT_MODE, GBAR_MODE )

                     VFAC       = VOL * INV_LAMBDA(v)
                     BEXT_MODE  = BEXT_MODE  * VFAC
                     BSCAT_MODE = BSCAT_MODE * VFAC

                  ENDIF

                  !***sum to get total extinction and scattering
                  !***  and contribution to the overall asymmetry factor

                  BEXT (v) = BEXT (v) + BEXT_MODE
                  BSCAT(v) = BSCAT(v) + BSCAT_MODE
                  G_BAR(v) = G_BAR(v) + BSCAT_MODE * GBAR_MODE

               enddo

            END IF

         END DO MODE_LOOP

         G_BAR(1:nwl) = G_BAR(1:nwl) / MAX( BSCAT(1:nwl), 1.0E-28 )

         AERO_SCAT_COEF( L,1:nwl ) = max( BSCAT(1:nwl), scat_min)
         AERO_EXTI_COEF( L,1:nwl ) = max( BEXT (1:nwl), exti_min)
         AERO_ASYM_FAC ( L,1:nwl ) = min( G_BAR(1:nwl), asym_max)

      END DO LAY_LOOP

      END SUBROUTINE GET_AERO_DATA

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE FAST_OPTICS( NR, NI, LAMBDA, SIGMA_G, XLNSIG, MIE_PAR, 
     &                        expfac2, expfac3, BETA_EXT, BETA_SCAT, G )
C-----------------------------------------------------------------------
C  A parameterization of the aerosol extinction and scattering code
C     Coded by Dr. Francis S. Binkowski
C     CEMP, The Institute for the Environment
C     The University of North Carolina at Chapel Hill
C     email: frank_binkowski@unc.edu
C     Code dates. Began February 25, 2005, current March 18, 2005
C     modified by FSB May 10, 2005 to calculate asymmetry factor by the
C     semi-empirical method of Hanna and Mathur. Note that
C     the normalized scattering coefficient (non-dimensional) is
C     interpreted as an average of  scattering efficiency over the
C     size distribution
C     Modified 05/25/2005 to ensure that if ni = 0.0, that
C     beta_scat = beta_ext.
C     Modified 03/17/2008 to ensure match with reference calculation for
C     all cases of sigma_g.
C     Modified 10/25/2009 to be consistent with research version of code.
C     The modifications include better representation of the
C     asymmetry factor.
C     The form of the mathematical expression for extinction
C     and scattering is based upon that of
C     Heintzenberg & Baker (1976), and further inspired by the work
C     of Willeke and Brockmann(1977) who showed the basic shape of
C     the extinction resulting from an integration of the Mie extinction
C     efficiency over a log normal particle size distribution.
C     This current code has been developed from a direct calculation of
C     extinction and scattering using a Mie code integrated over
C     a log normal distribution. The results of that calculation
C     replicate the Willee and Brockmann (1977) results.
C     The parameterization is based upon a fit to the standard CMAQ
C     radiative transfer code for visibility. to which a module for
C     absorption was added so that absorption and scattering would also be
C     available.
C
C     The  shape has been modified for absorption by small particles
C     by using the formulae of Penndorf (1962). For large particles,
C     the shape has been altered by using the large sphere values from
C     Table I of  Irvine (1965)
C
C     Results are consistent with those of Willeke and Brockmann (1977).
C
C     The asymmetry factor is calculated by a new algorithm FSB - 04/04/2008
C
C     Oct  2014: B.Hutzell-Added mininum values for scatttering efficiency for small
C                highly absorbing particles such as Aitken particle that
C                are mostly composed of elemental carbon
C
C  References:
C
C     Heintzenberg, J. and M. Baker, Spherical particle populations:
C     approximate analytic relationship between size distribution
C     parameters and integral optical properties, Applied Optics,
C     Volume 15, Number 5 pp 1178-1181, May 1976.
C
C     Irvine, W. M. Light Scattering by spherical particles: Radiation
C     pressure, asymmetry factor, and extinction cross section,
C     Journal of eh Optical Society of Amerioca, Vol. 55, NO. 1,
C     January ,1965
C
C     Penndorf, R. Scattering and Extinction Coefficeints for small
C     spherical particles, J. Atmospheric Sciences, Volume 19, p 193,
C     March 1962.
C
C     Willeke, K. and  J. E. Brockmann, Extinction coefficients for
C     multmodal atmospheric particle size distributions, Atmospheric
C     Environment, vol. 11, pp 95-999, 1977.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C***include files

      INCLUDE 'CONST.EXT'      ! physical constants

C***arguments

      REAL, INTENT( IN )     :: NR, NI     ! indices of refraction
      REAL, INTENT( IN )     :: LAMBDA     ! wavelength, nm
      REAL, INTENT( IN )     :: SIGMA_G    ! geometric standard deviation
      REAL, INTENT( IN )     :: XLNSIG     ! log of geometric standard deviation
      REAL, INTENT( IN )     :: mie_par
      REAL, INTENT( IN )     :: expfac2
      REAL, INTENT( IN )     :: expfac3

      REAL, INTENT( OUT )    :: BETA_EXT   ! normalized extinction coefficient
      REAL, INTENT( OUT )    :: BETA_SCAT  ! normalized scattering coefficient
      REAL, INTENT( OUT )    :: G          ! asymmetry factor

C***internal variables

      REAL NR1, NR2, TAU
      REAL C, CEXT, CSCAT
      REAL B, BEXT, BSCAT
      REAL BBFAC
      REAL ALPHV
      REAL ALPHA_I
      REAL A, LOGX2, XLNSIG2, MM1

      REAL, PARAMETER :: SIGMA_G_TWO = 2.0
      REAL, PARAMETER :: LOGSIG_TWO = 0.693147
      REAL, PARAMETER :: LOGSIG_TWO_SQD = LOGSIG_TWO * LOGSIG_TWO
      REAL, PARAMETER :: LOGSIG_105 = 1.102500
      REAL, PARAMETER :: DIFFSIG = SIGMA_G_TWO - 1.05
      REAL, PARAMETER :: A_TWO = 0.5 / LOGSIG_TWO_SQD
      REAL LARGESCAT  ! large sphere limit for scattering
      REAL LARGEEXT   ! large sphere limit for extinction
      REAL SMALL_G, LARGE_G

      REAL ALPHV2, ALPHV3
      REAL X_ALPHA, X_ALPHA2, X_ALPHA3
      REAL FCORR
!     REAL EXPFAC2, EXPFAC3
      REAL EFAC_EXT, EFAC_SCAT
      REAL FAC1
      REAL F1, F2, F3
      REAL G1, G2, G3, G4, G5
      REAL NN, TT

      REAL T1F1, T1F2, T2F2, T1F3, T2F3
      REAL T1G1, T2G1, T1G2, T2G2, T1G3, T2G3, T1G4, T2G4
      REAL T1G5, T2G5, T1P1, T2P1

C***the following are for calculating the Penndorff Coefficients

      REAL PENN1, A1, A2, A3, PENN2
      REAL XNR, XNI, XNR2, XNI2, XNRI, XNRI2, XNRMI
      REAL XRI, XRI2, XRI36, XNX, XNX2
      REAL Z1, Z12, Z2, XC1

C***coefficients for polynomials

      REAL, PARAMETER :: F1A0 = -1.309193E-1
      REAL, PARAMETER :: F1A1 =  2.565668E+0

      REAL, PARAMETER :: F2A0 =  1.489233E+1
      REAL, PARAMETER :: F2A1 = -5.311351E+1
      REAL, PARAMETER :: F2A2 =  8.180334E+1
      REAL, PARAMETER :: F2A3 = -4.549854E+1

      REAL, PARAMETER :: F3A0 =  5.549359E-1
      REAL, PARAMETER :: F3A1 =  2.571002E-1
      REAL, PARAMETER :: F3A2 =  9.212703E-1
      REAL, PARAMETER :: F3A3 =  5.155047E-2

      REAL, PARAMETER :: G1A0 =  1.000000E+0
      REAL, PARAMETER :: G1A1 = -7.636121E-1
      REAL, PARAMETER :: G1A2 =  1.285532E+0
      REAL, PARAMETER :: G1A3 =  1.629161E-2

      REAL, PARAMETER :: G2A0 =  1.000000E+0
      REAL, PARAMETER :: G2A1 = -2.507954E+0
      REAL, PARAMETER :: G2A2 =  5.170246E+0
      REAL, PARAMETER :: G2A3 = -4.945515E+0

      REAL, PARAMETER :: G3A0 =  1.000000E+0
      REAL, PARAMETER :: G3A1 = -2.650134E+0
      REAL, PARAMETER :: G3A2 =  4.515327E+0
      REAL, PARAMETER :: G3A3 = -4.027442E+0

      REAL, PARAMETER :: G4A0 =  1.000000E+0
      REAL, PARAMETER :: G4A1 = -8.532647E-1
      REAL, PARAMETER :: G4A2 =  6.766073E-1

      REAL, PARAMETER :: G5A0 =  1.0
      REAL, PARAMETER :: G5A1 = -6.674946E-1
      REAL, PARAMETER :: G5A2 = -5.426378E-1

C *** Coefficients for quadratic fit for <cos> large particles
C     from Irvine (1965) Table I
      REAL, PARAMETER :: IRVA0 =  1.076232E+0
      REAL, PARAMETER :: IRVA1 = -4.891972E-2
      REAL, PARAMETER :: IRVA2 = -2.243449E-2


C FSB Coefficients for SMALL_G  calculation
      REAL, PARAMETER :: XXA0 =  3.392224E-2
      REAL, PARAMETER :: XXA1 =  8.276532E-1
      REAL, PARAMETER :: XXA2 = -3.784926E-1
      REAL, PARAMETER :: XXA3 =  5.853108E-2
      REAL, PARAMETER :: G_ALPHA_I = XXA0 + XXA1 + XXA2 + XXA3

C FSB Coefficients for LARGE_G calculation
      REAL, PARAMETER :: YYA0 =  6.776173E-1
      REAL, PARAMETER :: YYA1 =  2.385339E-2
      REAL, PARAMETER :: YYA2 = -7.952263E-4
      REAL, PARAMETER :: YYA3 =  9.410370E-6

C FSB Coefficientf for adjusting asymmetry factor calculation
      REAL, PARAMETER :: CXG0 = -5.845196E-2
      REAL, PARAMETER :: CXG1 =  3.229183E-1
      REAL, PARAMETER :: CXG2 =  2.764716E-1
      REAL, PARAMETER :: CXG3 = -8.790411E-2

      REAL, PARAMETER :: AA0 =  88.466
      REAL, PARAMETER :: AA1 = -61.628
      REAL, PARAMETER :: AA2 =  34.483
      REAL, PARAMETER :: AA3 =  -8.543
      REAL, PARAMETER :: AA4 =   0.77435

      REAL QF1, QF2, QF3
      REAL QQSUM, QQF1,QQF2, QQF3, QQCORR

      REAL, PARAMETER :: DEGTORAD = PI180
      REAL, PARAMETER :: THREE_PI_TWO = 3.0 * PI / 2.0

      REAL, PARAMETER :: SCALE = 1.00E+9

C FSB check range of SIGMA_G
C *** Maximum value of SIGMA_G allowed is 2.0
!     SIGMA_G = MIN( 2.0, SIGMA_G )
C *** Minimum allowed value fo SIGMA_F is 1.05
!     SIGMA_G = MAX( 1.05, SIGMA_G )
C***FSB start calculation
!     XLNSIG = LOG( SIGMA_G ) 

      ALPHV  = SCALE * PI * MIE_PAR / LAMBDA
      ALPHV2 = ALPHV * ALPHV
      ALPHV3 = ALPHV * ALPHV * ALPHV

      XLNSIG2 = XLNSIG * XLNSIG
      A = 0.5 / XLNSIG2

      NR1 = NR - 1.0
      NR2 = NR * NR

C***evaluate polynomials
C***  optimize for piplined microprocessor
C***  appproach recommended by Dr. Carlie Coats

      NN = NR1 * NR1

      T1F1 = F1A0 + F1A1 * NR1
      T1F2 = F2A0 + F2A1 * NR1
      T2F2 = F2A2 + F2A3 * NR1
      T1F3 = F3A0 + F3A1 * NR1
      T2F3 = F3A2 + F3A3 * NR1

      F1 = T1F1                 ! linear
      F2 = T1F2 + T2F2 * NN     ! cubic
      F3 = T1F3 + T2F3 * NN     ! cubic

      C       = F1

C FSB correct for values of SIGMA_ less than 2.0.
C     The smallest value of SIGMA_G allowed is 1.05
C     2.0 - 1.05 = 0.95
C      FCORR   = (1.0 - 0.13 * (SIGMA_G - 1.05)  / 0.95 )
C     1.0 / 0.95 = 1.052632.  1.0 / 0.95 * 0.13 = 0.136842

      FCORR   = 1.0 - 0.136842 * ( SIGMA_G - 1.05 )

      CEXT    = C * FCORR
      CSCAT   = CEXT
      B       = F3 * A_TWO
      ALPHA_I = F2
      BEXT    = B
      BSCAT   = B
      PENN1   = 0.0
      PENN2   = 0.0

      IF ( NI .GT. 0.0 ) THEN

         TAU = NI / NR1

C***evaluate more polynomials
C***  optimize for piplined microprocessor
C***  appproach recommended by Dr. Carlie Coats

         TT = TAU * TAU
         T1G1 = G1A0 + G1A1 * TAU
         T2G1 = G1A2 + G1A3 * TAU
         T1G2 = G2A0 + G2A1 * TAU
         T2G2 = G2A2 + G2A3 * TAU
         T1G3 = G3A0 + G3A1 * TAU
         T2G3 = G3A2 + G3A3 * TAU
         T1G4 = G4A0 + G4A1 * TAU
         T2G4 = G4A2            ! quadratic
         T1G5 = G5A0 + G5A1 * TAU
         T2G5 = G5A2            ! quadratic
         G1   = T1G1 + T2G1 * TT
         G2   = T1G2 + T2G2 * TT
         G3   = T1G3 + T2G3 * TT
         G4   = T1G4 + T2G4 * TT
         G5   = T1G5 + T2G5 * TT

C *** adjust the variables

         CEXT    = CEXT    * G1
         CSCAT   = CSCAT   * G2
         BEXT    = BEXT    * G3
         BSCAT   = BSCAT   * G4
         ALPHA_I = ALPHA_I * G5

C*** Calculate the Penndorf Coefficients for the small particle limit

         XNR   = NR
         XNI   = NI
         XNR2  = XNR   * XNR
         XNI2  = XNI   * XNI
         XNRI  = XNR2  + XNI2
         XNRI2 = XNRI  * XNRI
         XNRMI = XNR2  - XNI2
         XRI   = XNR   * XNI
         XRI2  = XRI   * XRI
         XRI36 = 36.0  * XRI2
         XNX   = XNRI2 + XNRMI - 2.0
         XNX2  = XNX   * XNX

         Z1    = XNRI2 + 4.0 * XNRMI + 4.0
         Z12   = Z1    * Z1
         Z2    = 4.0   * XNRI2 + 12.0 * XNRMI + 9.0
         XC1   = 8.0   / ( 3.0 * Z12 )
         A1    = 24.0  * XRI / Z1

         A2    = 4.0   * XRI / 15.0 + 20.0 * XRI / ( 3.0 * Z2 ) +
     &           4.8   * XRI * (  7.0 * XNRI2 +
     &           4.0   * ( XNRMI - 5.0 ) ) / Z12

         A3    = XC1   * ( XNX2 - XRI36 )

!        EXPFAC2 = EXP( 2.0 * XLNSIG2 )
!        EXPFAC3 = EXP( 4.5 * XLNSIG2 )

         T1P1 = A1 + A2 * ALPHV2 * EXPFAC2
         T2P1 = A3 * ALPHV3 * EXPFAC3

C***PENN1 is the analytic integral of the Pendorff formulae over
C***   a log normal particle size distribution.

         PENN1 = THREE_PI_TWO * ( T1P1 + T2P1 )
         PENN2 = THREE_PI_TWO * T2P1

      END IF                    ! test of ni > 0.0

      X_ALPHA  = ALPHV / ALPHA_I

      LOGX2    = LOG( X_ALPHA )** 2

      BBFAC    = BEXT * A  / ( BEXT + A )

      EFAC_EXT = EXP( -BBFAC * LOGX2 )

C***FSB calculate normalized extinction and scattering coefficients

      BETA_EXT = THREE_PI_TWO * CEXT *
     &                      SQRT( A / ( BEXT + A ) ) * EFAC_EXT

      BETA_SCAT = BETA_EXT      ! NI = 0.0 case

C *** Check for 0.0 < NI

      IF ( NI .GT. 0.0 ) THEN

!        BBFAC = 1.5 * BBFAC ! reset BBFAC for scattering
         BBFAC = 1.1 * BBFAC ! reset BBFAC for scattering

         EFAC_SCAT = EXP( -BBFAC * LOGX2 )

C *** recalculate the normalized scattering coefficient

         BETA_SCAT = THREE_PI_TWO * CSCAT *
     &                       SQRT( A / ( BSCAT + A ) ) * EFAC_SCAT

C *** Adjust beta_ext for small particle absorption:

         IF ( X_ALPHA .LT. 0.13 ) THEN

            FAC1 = ALPHV + 0.6 * ALPHV2

            BETA_EXT = MAX( BETA_EXT,
     &                 ( BETA_EXT  * FAC1 + ( 1.0 - FAC1 ) * PENN1 ) )

            BETA_SCAT = MAX( BETA_SCAT,
     &                 ( BETA_SCAT * FAC1 + ( 1.0 - FAC1 ) * PENN2 ) )

         END IF ! test for   XALPHA < 0.13
        
! small highly absorbing  particles such as a Tyndall scatterers such as
! described in Jacobson (1999) in Fundmentals of Atmospheric Modeling. This
! is approximately consistent with Penndorf (1962).

         BETA_SCAT = MAX( 0.001 * BETA_EXT, BETA_SCAT )
          
      END IF  ! test for 0.0  < NI

C *** Calculate and adjust for large sphere limits(Irvine, 1965)

      IF ( ALPHV .GE. 10.0 ) THEN

         MM1 = EXP( 0.5 * XLNSIG2 ) / ALPHV

C FSB large sphere limit - scattering

         LARGESCAT = THREE_PI_TWO * ( 8.652439E-1 + 1.501772E-1 * NR ) * MM1
         BETA_SCAT = MAX( BETA_SCAT, LARGESCAT )

C FSB large spnere limit for extinction ( no edge effectss )

         LARGEEXT = THREE_PI_TWO * ( 2.0 * MM1 )
         BETA_EXT = MAX( BETA_EXT, LARGEEXT )

      END IF ! test for 10.0 < ALPHV

C***FSB now calculate the asymmetry factor
C***  using the Hanna-Mathur quasi-empirical method
!     GMAX = 0.0

!     QSCAT_AVG = BETA_SCAT
!     QANGLE = AA0

!     IF ( QSCAT_AVG .GT. 1.0E-4 ) THEN
!        QQ  = QSCAT_AVG * QSCAT_AVG
!        QF1 = AA0 + AA1 * QSCAT_AVG
!        QF2 = AA2 + AA3 * QSCAT_AVG
!        QF3 = AA4 * QSCAT_AVG * QSCAT_AVG
!        QANGLE = QF1 + QF2 * QQ + QF3 * QQ
!     END IF

!     QANGLE = DEGTORAD * QANGLE
!     G      = COS( QANGLE)     ! asymmetry factor

C***FSB now calculate the asymmetry factor
C***   using a parametric fit to Mie calculations

      X_ALPHA2 = X_ALPHA * X_ALPHA
      X_ALPHA3 = X_ALPHA * X_ALPHA * X_ALPHA

      SMALL_G = 0.0
      LARGE_G = 0.0

      IF ( ALPHV .LT. 2.6 ) THEN

C FSB Calculate SMALL_G for sigma_g = 2.0

         QF1 = XXA0 + XXA1 * ALPHV
         QF2 = XXA2 * ALPHV2
         QF3 = XXA3 * ALPHV3
         SMALL_G = QF1 + QF2 + QF3 ! valid at sigma_g = 2.0

C FSB adjust SMALL_G for values less than 2.0
C     The smallest sigma_g allowed is 1.05

         QQF1 = CXG0 + CXG1 * ALPHV
         QQF2 = CXG2 * ALPHV2
         QQF3 = CXG3 * ALPHV3
         QQSUM = QQF1 + QQF2 + QQF3
         FAC1 = 1.052632 * ( 2.0 - SIGMA_G)
         QQSUM = MAX( 0.0, QQSUM )
         QQCORR =   ( 1.0 - FAC1 ) + FAC1 * QQSUM
         SMALL_G = SMALL_G * QQCORR
         G = SMALL_G

      ELSE IF ( ALPHV .LT. 100.0 ) THEN

         QF1 =  YYA0 + YYA1 * X_ALPHA
         QF2 = YYA2 * X_ALPHA2
         QF3 = YYA3 * X_ALPHA3
         LARGE_G = QF1 + QF2 + QF3
         G = LARGE_G

      ELSE

C     Very large sphere limit from fit to Irvine (1965).

         G = 1.124484 - 1.153869E-1 * NR

      END IF

C Simple and gross fix if appromation exceed realistic bounds
      G = MAX( -0.9999, MIN( G, 0.9999))

      END SUBROUTINE FAST_OPTICS

      END MODULE AERO_PHOTDATA
