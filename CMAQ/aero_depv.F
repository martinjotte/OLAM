
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
C $Header: /project/yoj/arc/CCTM/src/aero/aero5/aero_depv.F,v 1.12 2012/01/19 13:12:14 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)aero_depv.F       1.3 /project/mod3/CMAQ/src/ae_depv/aero_depv/SCCS/s.aero_depv.F 18 Jun 1997 12:55:48

      SUBROUTINE AERO_DEPV ( iw )

C aerosol dry deposition routine
C   written 4/9/97 by Dr. Francis S. Binkowski
C   uses code from modpar and vdvg from the aerosol module.
C   This routine uses a single block to hold information
C   for the lowest layer.
C NOTES: This version assumes that RA is available on the met file.
c        Array structure for vector optimization
C 26 Apr 97 Jeff - many mods
C 13 Dec 97 Jeff - expect uncoupled CGRID, concs as micro-g/m**3, #/m**3
C
C 1/11/99 David Wong at LM - change NUMCELLS to CELLNUM in the loop index
C FSB 3/17/99 changed to accommodate surface area/second moment and
C    encapsulated the actual drydep calculation into a subroutine which
C    is attached to this code
C Jeff - Dec 00 - move CGRID_MAP into f90 module
C FSB 12/11/2000. Logic added to allow deposition of particles at their
C     "wet" diameters; that is, accounting for the water on the particles.
C     This is done by adjusting the third and second moments for the
C     presence of water assuming that the geometric standard deviations
C     are not changed by this process. This appears to be a very good
C     assumption.
C 30 Aug 01 J.Young: Dyn alloc; Use HGRD_DEFN
C    Jan 03 J.Young: Change CGRID dimensions, eliminate re-allocations
C  6 Mar 03 J.Young: eliminate a lot of allocate/deallocates
C  7 Aug 03 S.Roselle: updated code for loading the min aero conc array
C 17 Dec 03 S.Roselle: Adjust 2nd and 3rd moments to include SOA,
C     without affecting the geometric standard deviations.
C 31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C     domain specifications in one module
C 07 Jun 05 P.Bhave: Added code to handle new species in the AE4
C     mechanism: ANAI, ANAJ, ANAK, ACLI, ACLJ, ACLK, ASO4K, AH2OK,
C     and ANO3K; look for ASEAS only when using AE3 mechanism
C 30 Jan 08 S.Napelenok & P.Bhave: Added code to handle new SOA species
C     in AE5; defined DRY aerosol to include nonvolatile SOA spcs
C 14 Apr 08 J.Kelly: Added code to handle new species ANH4K and SRFCOR.
C     Also added code to handle variable coarse mode standard deviation
C     in AE5 (no longer fixed at 2.2).
C 08 Sep 08 P.Bhave: Backward compatibility with AE4 mechanisms
C     standardized names of all coarse-mode variables
C 19 Apr 10 S.Howard: aero re-engineering for modularity
C 23 Apr 10 J.Young: replace chem mechanism include files with namelists
C 10 Mar 11 S.Howard: Renamed met_data to aeromet_data
C 25 Mar 11 S.Roselle: Replaced I/O API include files with UTILIO_DEFN
C 20 May 11 D.Schwede: Modified for mosaic
C 31 Aug 11 J.Bash: Moved shared mosaic variables to MOSAIC_MOD
C 27 Sep 11 David Wong: replaced all run time dynamic arrays with allocatable
C                       arrays to avoid run time memory issue  
C 08 Jun 12 J.Young: remove full character blank padding for GNU Fortran (GCC) 4.1.2
C 07 Jul 14 B.Hutzell: replaced mechanism include file(s) with fortran module
C 07 Nov 14 J.Bash: Updated for the ASX_DATA_MOD shared data module.
C    May 16 B. Murphy, H. Pye: Updated treatment of aerosol moments

      USE RXNS_DATA           ! chemical mechanism data
      USE CGRID_SPCS          ! CGRID mechanism species
      USE UTILIO_DEFN      
      USE AERO_DATA           ! aero variable data

      use const_data, only: grav, stdatmpa
      use getpar_mod, only: getpar
      use cgrid_conv, only: rev_cgrid_sfc_iw
      use depv_defn
      use mem_ijtabs, only: itab_w
      use mem_leaf,   only: land, itab_wl
      use mem_sea,    only: sea,  itab_ws
      use mem_turb,   only: pblh
      use mem_basic,  only: rho, press, theta
      use mem_grid,   only: lsw, lpw

      IMPLICIT NONE

C Arguments
      integer, intent(in) :: iw

C Parameters
      REAL,    PARAMETER :: T0 = 288.15      ! [ K ] ! starting standard surface temp.
      REAL,    PARAMETER :: ONE3 = 1.0 / 3.0
      REAL,    PARAMETER :: TWO3 = 2.0 / 3.0
      INTEGER, PARAMETER :: N_AE_DEP_SPC = 9 ! no. of surrogates for aerosol dry dep velocities

C set up species dimension and indices for deposition velocity internal array VDEP
      INTEGER, PARAMETER :: VDNATK = 1,  ! Aitken mode number
     &                      VDNACC = 2,  ! accumulation mode number
     &                      VDNCOR = 3,  ! coarse mode number
     &                      VDMATK = 4,  ! Aitken mode mass
     &                      VDMACC = 5,  ! accumulation mode mass
     &                      VDMCOR = 6,  ! coarse mode mass
     &                      VDSATK = 7,  ! Aitken mode surface area
     &                      VDSACC = 8,  ! accumulation mode surface area
     &                      VDSCOR = 9   ! coarse mode surface area

C Local variables:

      CHARACTER( 16 ) :: VDAE_NAME( N_AE_DEP_SPC ) ! dep vel surrogate name table
      DATA         VDAE_NAME( 1 ) / 'VNUMATKN' /
      DATA         VDAE_NAME( 2 ) / 'VNUMACC ' /
      DATA         VDAE_NAME( 3 ) / 'VNUMCOR ' /
      DATA         VDAE_NAME( 4 ) / 'VMASSI  ' /
      DATA         VDAE_NAME( 5 ) / 'VMASSJ  ' /
      DATA         VDAE_NAME( 6 ) / 'VMASSC  ' /
      DATA         VDAE_NAME( 7 ) / 'VSRFATKN' /
      DATA         VDAE_NAME( 8 ) / 'VSRFACC ' /
      DATA         VDAE_NAME( 9 ) / 'VSRFCOR ' /

      INTEGER, ALLOCATABLE, SAVE :: DEPV_SUR( : )   ! pointer to surrogate

C Meteorological variables

      REAL :: XXLSGAT(lsw(iw))           ! log of standard deviation
      REAL :: XXLSGAC(lsw(iw))
      REAL :: XXLSGCO(lsw(iw))
 
      REAL :: DGATK(lsw(iw))             ! geometric mean diameter
      REAL :: DGACC(lsw(iw))
      REAL :: DGCOR(lsw(iw))

      REAL :: PDENSAT(lsw(iw))           ! particle density         
      REAL :: PDENSAC(lsw(iw))
      REAL :: PDENSCO(lsw(iw))

      real :: airdens
      real :: airpres
      real :: airtemp

      REAL :: XLM               ! mean free path [ m ]
      REAL :: AMU               ! dynamic viscosity [ kg m**-1 s**-1 ]

      real :: vdep(n_ae_dep_spc)  ! deposition  velocity [ m/s ]

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: LOGDEV                   ! unit number for the log file
      CHARACTER( 16 ), SAVE :: PNAME = 'AERO_DEPV'
      CHARACTER( 96 ) :: XMSG = ' '

      INTEGER  V, N                  ! loop counters
      INTEGER  ALLOCSTAT

      INTEGER :: JDATE               ! current model date , coded YYYYDDD
      INTEGER :: JTIME               ! current model time , coded HHMMSS

      real    :: ustar, wstar, ra, wtv0
      integer :: jwls, iwls, ks, kw

      real :: cngrd(lsw(iw),nspcsd) 

      jdate = 0
      jtime = 0

      IF ( FIRSTIME ) THEN

         !$omp barrier
         !$omp single
         FIRSTIME = .FALSE.
         LOGDEV = INIT3()

         ALLOCATE ( DEPV_SUR( N_AE_DEPV ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DEPV_SUR'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         ! Set the dep vel surrogate pointers
         DO V = 1, N_AE_DEPV
            N = INDEX1( AE_DEPV( V ), N_AE_DEP_SPC, VDAE_NAME )
            IF ( N .NE. 0 ) THEN
               DEPV_SUR( V ) = N
            ELSE
               XMSG = 'Could not find ' // AE_DEPV( V ) // ' in aerosol' //
     &                ' surrogate table. >>> Dep vel set to zero <<< '
               CALL M3WARN( PNAME, JDATE, JTIME, XMSG )
               DEPV_SUR( V ) = 0
            END IF
         END DO
         !$omp end single

      END IF

      IF ( N_AE_SPC .LE. 0 ) RETURN

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Put the grid cell physical data in the block arrays
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call rev_cgrid_sfc_iw( iw, cngrd )

      do ks = 1, lsw(iw)

         ! Extract grid cell concentrations of aero species from CGRID
         ! into aerospc_conc in aero_data module. Also determines
         ! second moment from surface area and adds wet species.

         call extract_aero( cngrd( ks,: ), .true. )

         ! Calculate geometric mean diameters and standard deviations of
         ! the "wet" size distribution

         call getpar( .false., noM3=.true. )

         ! Save getpar values to arrays

         XXLSGAT(ks) = AEROMODE_LNSG( 1 )
         XXLSGAC(ks) = AEROMODE_LNSG( 2 )
         XXLSGCO(ks) = AEROMODE_LNSG( 3 )

         DGATK(ks)   = AEROMODE_DIAM( 1 )
         DGACC(ks)   = AEROMODE_DIAM( 2 )
         DGCOR(ks)   = AEROMODE_DIAM( 3 )

         PDENSAT(ks) = AEROMODE_DENS( 1 )
         PDENSAC(ks) = AEROMODE_DENS( 2 )
         PDENSCO(ks) = AEROMODE_DENS( 3 )

      enddo

      ! Land Cells
      do jwls = 1, itab_w(iw)%nland

         iwls = itab_w(iw)%iland(jwls)
         kw   = itab_wl(iwls)%kw
         ks   = kw - lpw(iw) + 1

         ! Set meteorological data for this surface cell
         AIRDENS = rho(kw,iw)
         AIRPRES = press(kw,iw)
         AIRTEMP = land%cantemp(iwls)
         ra      = 1.0 / land%ggaer(iwls)
         ustar   = land%ustar(iwls)
         wtv0    = land%wthv(iwls)

         if (wtv0 > 0.0) then
            wstar = (grav * pblh(iw) * wtv0 / theta(kw,iw)) ** one3
         else
            wstar = 0.0
         endif

         ! Calculate mean free path [ m ]:
         XLM = 6.6328E-8 * STDATMPA * AIRTEMP / ( T0 * AIRPRES )

         ! Calcualte dynamic viscosity [ kg m**-1 s**-1 ]:
         AMU = 1.458E-6 * AIRTEMP * SQRT( AIRTEMP )
     &       / ( AIRTEMP + 110.4 )

         ! get dry deposition velocities:

         CALL GETDEP_V ( XLM, AMU, DGATK(ks), DGACC(ks), DGCOR(ks),
     &                   XXLSGAT(ks), XXLSGAC(ks), XXLSGCO(ks),
     &                   PDENSAT(ks), PDENSAC(ks), PDENSCO(ks), VDEP, 
     &                   AIRTEMP, AIRDENS, ustar, wstar, ra )

         DO V = 1, N_AE_DEPV
            IF ( DEPV_SUR( V ) .GT. 0 ) THEN
               adepv_land( iwls,V ) = VDEP( DEPV_SUR( V ) )
            ELSE
               adepv_land( iwls,V ) = 0.0
            END IF
         END DO

      enddo

      ! SEA Cells
      do jwls = 1, itab_w(iw)%nsea
         iwls = itab_w(iw)%isea(jwls)
         kw   = itab_ws(iwls)%kw
         ks   = kw - lpw(iw) + 1

         ! Set meteorological data for this surface cell
         AIRDENS = rho(kw,iw)
         AIRPRES = press(kw,iw)
         AIRTEMP = sea%cantemp(iwls)
         ra      = 1.0 / sea%ggaer(iwls)
         ustar   = sea%ustar(iwls)
         wtv0    = sea%wthv(iwls)

         if (wtv0 > 0.0) then
            wstar = (grav * pblh(iw) * wtv0 / theta(kw,iw)) ** one3
         else
            wstar = 0.0
         endif

         ! Calculate mean free path [ m ]:
         XLM = 6.6328E-8 * STDATMPA * AIRTEMP / ( T0 * AIRPRES )

         ! Calcualte dynamic viscosity [ kg m**-1 s**-1 ]:
         AMU = 1.458E-6 * AIRTEMP * SQRT( AIRTEMP ) / ( AIRTEMP + 110.4 )

         ! get dry deposition velocities:

         CALL GETDEP_V ( XLM, AMU, DGATK(ks), DGACC(ks), DGCOR(ks),
     &                   XXLSGAT(ks), XXLSGAC(ks), XXLSGCO(ks),
     &                   PDENSAT(ks), PDENSAC(ks), PDENSCO(ks), VDEP, 
     &                   AIRTEMP, AIRDENS, ustar, wstar, ra )

         DO V = 1, N_AE_DEPV
            IF ( DEPV_SUR( V ) .GT. 0 ) THEN
               adepv_sea( iwls,V ) = VDEP( DEPV_SUR( V ) )
            ELSE
               adepv_sea( iwls,V ) = 0.0
            END IF
         END DO
            
      enddo

      contains

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE GETDEP_V ( XLM, AMU, DGATK, DGACC, DGCOR,
     &                      XXLSGAT, XXLSGAC, XXLSGCO,
     &                      PDENSAT, PDENSAC, PDENSCO,VDEP, 
     &                      TEMP2, DENS1, USTAR, WSTAR, RA )

C *** Calculate deposition velocity for Aitken, accumulation, and
C     coarse modes.
C     Reference:
C     Binkowski F. S., and U. Shankar, The regional particulate
C     model 1. Model description and preliminary results.
C     J. Geophys. Res., 100, D12, 26191-26209, 1995.
 
C    May 05 D.Schwede: added impaction term to coarse mode dry deposition
C 25 May 05 J.Pleim:  Updated dry dep velocity calculation for aerosols
C                     to Venkatram and Pleim (1999)
C 20 Jul 05 J.Pleim:  Changed impaction term using modal integration of
C                     Stokes**2 / 400 (Giorgi, 1986, JGR)
C 14 Apr 08 J.Kelly:  Added code to calculate deposition velocity of
C                     coarse surface area and to account for variable
C                     standard deviation of the coarse mode.
C 08 Sep 08 P.Bhave:  Backward compatibility with AE4 mechanisms
C                     standardized names of all coarse-mode variables
C-----------------------------------------------------------------------

      use const_data, only: pi, boltzmann, grav

      IMPLICIT NONE

C *** input arguments

C     atmospheric properties
      REAL, INTENT( IN ) :: XLM      ! atmospheric mean free path [ m ]
      REAL, INTENT( IN ) :: AMU      ! atmospheric dynamic viscosity [ kg/(m s) ]

C     aerosol properties:

C     modal diameters: [ m ]
      REAL, INTENT( IN ) :: DGATK    ! nuclei mode geometric mean diameter  [ m ]
      REAL, INTENT( IN ) :: DGACC    ! accumulation geometric mean diameter [ m ]
      REAL, INTENT( IN ) :: DGCOR    ! coarse mode geometric mean diameter  [ m ]

C     log of modal geometric standard deviations
      REAL, INTENT( IN ) :: XXLSGAT  ! Aitken mode
      REAL, INTENT( IN ) :: XXLSGAC  ! accumulation mode
      REAL, INTENT( IN ) :: XXLSGCO  ! coarse mode

C     average modal particle densities  [ kg/m**3 ]
      REAL, INTENT( IN ) :: PDENSAT  ! average particle density in nuclei mode
      REAL, INTENT( IN ) :: PDENSAC  ! average particle density in accumulation mode
      REAL, INTENT( IN ) :: PDENSCO  ! average particle density in coarse mode

c     atmospheric variables
      real, intent( in ) :: temp2
      real, intent( in ) :: dens1
      real, intent( in ) :: ustar
      real, intent( in ) :: wstar
      real, intent( in ) :: ra

C *** output arguments

      REAL, INTENT( OUT ) :: VDEP( n_ae_dep_spc ) ! deposition  velocity [ m/s ]

C *** array indices hardcoded to match SUBROUTINE AERO_DEPV
      INTEGER, PARAMETER :: VDNATK = 1,  ! Aitken mode number
     &                      VDNACC = 2,  ! accumulation mode number
     &                      VDNCOR = 3,  ! coarse mode number
     &                      VDMATK = 4,  ! Aitken mode mass
     &                      VDMACC = 5,  ! accumulation mode mass
     &                      VDMCOR = 6,  ! coarse mode mass
     &                      VDSATK = 7,  ! Aitken mode surface area
     &                      VDSACC = 8,  ! accumulation mode surface area
     &                      VDSCOR = 9   ! coarse mode surface area

C modal Knudsen numbers
      REAL KNATK   ! Aitken mode Knudsen number
      REAL KNACC   ! accumulation "
      REAL KNCOR   ! coarse mode

C modal particle diffusivities for number, 2nd, and 3rd moment, or mass:
      REAL DCHAT0N, DCHAT0A, DCHAT0C
      REAL DCHAT2N, DCHAT2A, DCHAT2C
      REAL DCHAT3N, DCHAT3A, DCHAT3C

C modal sedimentation velocities for number, 2nd, and 3rd moment, or mass:
!     REAL VGHAT0N, VGHAT0A, VGHAT0C
!     REAL VGHAT2N, VGHAT2A, VGHAT2C
!     REAL VGHAT3N, VGHAT3A, VGHAT3C

      REAL DCONST1, DCONST1N, DCONST1A, DCONST1C
      REAL DCONST2, DCONST3N, DCONST3A, DCONST3C
      REAL SC0NI2, SC0AI2, SC0CI2    ! Inv Schmidt number squared for number
      REAL SC2NI2, SC2AI2, SC2CI2    ! Inv Schmidt number squared for 2ND MOMENT
      REAL SC3NI2, SC3AI2, SC3CI2    ! Inv Schmidt number squared for 3rd moment
      REAL STOKENF, STOKEAF, STOKECF ! Stokes number factors for each mode
      REAL RD0N, RD0A, RD0C     ! canopy resistance for number
      REAL RD2N, RD2A, RD2C     ! canopy resistance for 2nd moment
      REAL RD3N, RD3A, RD3C     ! canopy resisteance for 3rd moment
      REAL UTSCALI              ! scratch function of USTAR and WSTAR
      REAL NUI                  ! kinematic viscosity [ m**2 s**-1 ]
      REAL USTFAC               ! scratch function of USTAR, NU, and GRAV
      REAL TWOXLM               ! 2 X atmospheric mean free path

      REAL, PARAMETER :: BHAT    = 1.246 ! Constant from Cunningham slip correction
      REAL, PARAMETER :: THREEPI = 3.0 * PI
      REAL, PARAMETER :: ONE3    = 1.0 / 3.0
      REAL, PARAMETER :: TWO3    = 2.0 / 3.0
      real, parameter :: gravi   = 1.0 / grav

C Scalar variables for VARIABLE standard deviations.

      REAL    L2SGAT, L2SGAC   ! see usage
      REAL    L2SGCO

!     REAL    EAT1             ! Aitken mode exp( log^2( sigmag )/8 )
!     REAL    EAC1             ! accumulation mode exp( log^2( sigmag )/8 )
!     REAL    ECO1             ! coarse mode exp( log^2( sigmag )/8 )

      REAL    ESAT04           ! Aitken       " **4
      REAL    ESAC04           ! accumulation "
      REAL    ESCO04           ! coarse       "

      REAL    ESAT08           ! Aitken       " **8
      REAL    ESAC08           ! accumulation "
      REAL    ESCO08           ! coarse       "

      REAL    ESAT12           ! Aitken       " **12
      REAL    ESAC12           ! accumulation "     
      REAL    ESCO12           ! coarse       "     

      REAL    ESAT16           ! Aitken       " **16
      REAL    ESAC16           ! accumulation "
      REAL    ESCO16           ! coarse       "

      REAL    ESAT20           ! Aitken       " **20
      REAL    ESAC20           ! accumulation "
      REAL    ESCO20           ! coarse       "

      REAL    ESAT28           ! Aitken       " **28
      REAL    ESAC28           ! accumulation "
      REAL    ESCO28           ! coarse       "

      REAL    ESAT32           ! Aitken       " **32
      REAL    ESAC32           ! accumulation "
      REAL    ESCO32           ! coarse       "

      REAL    ESAT36           ! Aitken       " **36
      REAL    ESAC36           ! accumulation "
      REAL    ESCO36           ! coarse       "

      REAL    ESAT48           ! Aitken       " **48
      REAL    ESAC48           ! accumulation "     
      REAL    ESCO48           ! coarse       "     

      REAL    ESAT64           ! Aitken       " **64
      REAL    ESAC64           ! accumulation "
      REAL    ESCO64           ! coarse       "

      REAL    ESAT128          ! Aitken       " **128
      REAL    ESAC128          ! accumulation "
      REAL    ESCO128          ! coarse       "

      REAL    ESAT160          ! Aitken       " **160
      REAL    ESAC160          ! accumulation "
      REAL    ESCO160          ! coarse       "

      REAL    ESATM12          ! Aitken       " **(-12)
      REAL    ESACM12          ! accumulation "        
      REAL    ESCOM12          ! coarse       "        

      REAL    ESATM16          ! Aitken       " **(-16)
      REAL    ESACM16          ! accumulation "        
      REAL    ESCOM16          ! coarse       "        

      REAL    ESATM20          ! Aitken       " **(-20)
      REAL    ESACM20          ! accumulation "
      REAL    ESCOM20          ! coarse       "

      REAL    ESATM32          ! Aitken       " **(-32)
      REAL    ESACM32          ! accumulation "
      REAL    ESCOM32          ! coarse       "

      REAL    EIM              ! Impaction efficiency

C-----------------------------------------------------------------------

C *** Calculate Knudsen numbers

      TWOXLM = XLM + XLM
      KNATK = TWOXLM / DGATK
      KNACC = TWOXLM / DGACC
      KNCOR = TWOXLM / DGCOR

C *** Calculate functions of variable standard deviation.

      L2SGAT = XXLSGAT * XXLSGAT
      L2SGAC = XXLSGAC * XXLSGAC
      L2SGCO = XXLSGCO * XXLSGCO

!     EAT1   = EXP( 0.125 * L2SGAT )
!     EAC1   = EXP( 0.125 * L2SGAC )
!     ECO1   = EXP( 0.125 * L2SGCO )

      ESAT04  = EXP( 0.5 * L2SGAT )
      ESAC04  = EXP( 0.5 * L2SGAC )
      ESCO04  = EXP( 0.5 * L2SGCO )

      ESAT08  = ESAT04 * ESAT04
      ESAC08  = ESAC04 * ESAC04
      ESCO08  = ESCO04 * ESCO04

      ESAT12  = ESAT04 * ESAT08
      ESAC12  = ESAC04 * ESAC08
      ESCO12  = ESCO04 * ESCO08

      ESAT16  = ESAT08 * ESAT08
      ESAC16  = ESAC08 * ESAC08
      ESCO16  = ESCO08 * ESCO08

      ESAT20  = ESAT16 * ESAT04
      ESAC20  = ESAC16 * ESAC04
      ESCO20  = ESCO16 * ESCO04

      ESAT28  = ESAT20 * ESAT08
      ESAC28  = ESAC20 * ESAC08
      ESCO28  = ESCO20 * ESCO08

      ESAT32  = ESAT16 * ESAT16
      ESAC32  = ESAC16 * ESAC16
      ESCO32  = ESCO16 * ESCO16

      ESAT36  = ESAT16 * ESAT20
      ESAC36  = ESAC16 * ESAC20
      ESCO36  = ESCO16 * ESCO20

      ESAT48  = ESAT36 * ESAT12
      ESAC48  = ESAC36 * ESAC12
      ESCO48  = ESCO36 * ESCO12

      ESAT64  = ESAT32 * ESAT32
      ESAC64  = ESAC32 * ESAC32
      ESCO64  = ESCO32 * ESCO32

      ESAT128 = ESAT64 * ESAT64
      ESAC128 = ESAC64 * ESAC64
      ESCO128 = ESCO64 * ESCO64

      ESAT160 = ESAT128* ESAT32
      ESAC160 = ESAC128* ESAC32
      ESCO160 = ESCO128* ESCO32

C *** calculate inverses:

      ESATM12 = 1.0 / ESAT12
      ESACM12 = 1.0 / ESAC12
      ESCOM12 = 1.0 / ESCO12

      ESATM16 = 1.0 / ESAT16
      ESACM16 = 1.0 / ESAC16
      ESCOM16 = 1.0 / ESCO16

      ESATM20 = 1.0 / ESAT20
      ESACM20 = 1.0 / ESAC20
      ESCOM20 = 1.0 / ESCO20

      ESATM32 = 1.0 / ESAT32
      ESACM32 = 1.0 / ESAC32
      ESCOM32 = 1.0 / ESCO32

      DCONST1  = BOLTZMANN * TEMP2 / ( THREEPI * AMU )
      DCONST1N = DCONST1 / DGATK
      DCONST1A = DCONST1 / DGACC
      DCONST1C = DCONST1 / DGCOR
      DCONST2  = GRAV / ( 18.0 * AMU )
      DCONST3N = DCONST2 * PDENSAT * DGATK * DGATK
      DCONST3A = DCONST2 * PDENSAC * DGACC * DGACC
      DCONST3C = DCONST2 * PDENSCO * DGCOR * DGCOR

C i-mode
      DCHAT0N  = DCONST1N * ( ESAT04  + BHAT * KNATK * ESAT16 )
      DCHAT2N  = DCONST1N * ( ESATM12 + BHAT * KNATK * ESATM16 )
      DCHAT3N  = DCONST1N * ( ESATM20 + BHAT * KNATK * ESATM32 )
!     VGHAT0N  = DCONST3N * ( ESAT16  + BHAT * KNATK * ESAT04 )
!     VGHAT2N  = DCONST3N * ( ESAT48  + BHAT * KNATK * ESAT20 )
!     VGHAT3N  = DCONST3N * ( ESAT64  + BHAT * KNATK * ESAT28 )

C j-mode
      DCHAT0A  = DCONST1A * ( ESAC04  + BHAT * KNACC * ESAC16 )
      DCHAT2A  = DCONST1A * ( ESACM12 + BHAT * KNACC * ESACM16 )
      DCHAT3A  = DCONST1A * ( ESACM20 + BHAT * KNACC * ESACM32 )
!     VGHAT0A  = DCONST3A * ( ESAC16  + BHAT * KNACC * ESAC04 )
!     VGHAT2A  = DCONST3A * ( ESAC48  + BHAT * KNACC * ESAC20 )
!     VGHAT3A  = DCONST3A * ( ESAC64  + BHAT * KNACC * ESAC28 )

C coarse mode
      DCHAT0C  = DCONST1C * ( ESCO04  + BHAT * KNCOR * ESCO16 )
      DCHAT2C  = DCONST1C * ( ESCOM12 + BHAT * KNCOR * ESCOM16 )
      DCHAT3C  = DCONST1C * ( ESCOM20 + BHAT * KNCOR * ESCOM32 )
!     VGHAT0C  = DCONST3C * ( ESCO16  + BHAT * KNCOR * ESCO04 )
!     VGHAT2C  = DCONST3C * ( ESCO48  + BHAT * KNCOR * ESCO20 )
!     VGHAT3C  = DCONST3C * ( ESCO64  + BHAT * KNCOR * ESCO28 )

C now calculate the deposition velocities

      NUI     = DENS1 / AMU
      USTFAC  = NUI * GRAVI * USTAR ** 2
      STOKENF = 0.0025 * (DCONST3N * USTFAC) ** 2
      STOKEAF = 0.0025 * (DCONST3A * USTFAC) ** 2
      STOKECF = 0.0025 * (DCONST3C * USTFAC) ** 2
      UTSCALI = ustar / ( ustar**2 + 0.24 * wstar**2 )

C first do 0th moment for the deposition of number

C  Aitken mode
      SC0NI2 = (NUI * DCHAT0N) ** 2
      EIM    = MIN( STOKENF * ESAT64, 1.0 )
      RD0N   = UTSCALI / ( SC0NI2**ONE3 + EIM )

!     VDEP( VDNATK ) = VGHAT0N
!    &               / ( 1.0 - EXP( -VGHAT0N * ( RA + RD0N ) ) )
      VDEP( VDNATK ) = 1.0 / ( RA + RD0N )

C accumulation mode
      SC0AI2 = (NUI * DCHAT0A) ** 2
      EIM    = MIN( STOKEAF * ESAC64, 1.0 )
      RD0A   = UTSCALI / ( SC0AI2**ONE3 + EIM )
      
!     VDEP( VDNACC ) = VGHAT0A
!    &               / ( 1.0 - EXP( -VGHAT0A * ( RA + RD0A ) ) )
      VDEP( VDNACC ) = 1.0 / ( RA + RD0A )

C coarse mode
      SC0CI2 = (NUI * DCHAT0C) ** 2
      EIM    = MIN( STOKECF * ESCO64, 1.0 )
      RD0C   = UTSCALI / ( SC0CI2**ONE3 + EIM )

!     VDEP( VDNCOR ) = VGHAT0C
!    &               / ( 1.0 - EXP( -VGHAT0C * ( RA + RD0C ) ) )
      VDEP( VDNCOR ) = 1.0 / ( RA + RD0C )

C now do 2nd moment for the deposition of surface area

C  Aitken mode
      SC2NI2 = (NUI * DCHAT2N) ** 2
      EIM    = MIN( STOKENF * ESAT128, 1.0 )
      RD2N   = UTSCALI / ( SC2NI2**ONE3 + EIM )

!     VDEP( VDSATK ) = VGHAT2N
!    &               / ( 1.0 - EXP( -VGHAT2N * ( RA + RD2N ) ) )
      VDEP( VDSATK ) = 1.0 / ( RA + RD2N )

C accumulation mode
      SC2AI2 = (NUI * DCHAT2A) ** 2
      EIM    = MIN( STOKEAF * ESAC128, 1.0 )
      RD2A   = UTSCALI / ( SC2AI2**ONE3 + EIM )

!     VDEP( VDSACC ) = VGHAT2A
!    &               / ( 1.0 - EXP( -VGHAT2A * ( RA + RD2A ) ) )
      VDEP( VDSACC ) =  1.0 / ( RA + RD2A )

C coarse mode
      SC2CI2 = (NUI * DCHAT2C) ** 2
      EIM    = MIN( STOKECF * ESCO128, 1.0 )
      RD2C   = UTSCALI / ( SC2CI2**ONE3 + EIM )

!     VDEP( VDSCOR ) = VGHAT2C
!    &               / ( 1.0 - EXP( -VGHAT2C * ( RA + RD2C ) ) )
      VDEP( VDSCOR ) = 1.0 / ( RA + RD2C )

C now do 3rd moment for the deposition of mass

C  Aitken mode
      SC3NI2 = (NUI * DCHAT3N) ** 2
      EIM    = MIN( STOKENF * ESAT160, 1.0 )
      RD3N   = UTSCALI / ( SC3NI2**ONE3 + EIM )

!     VDEP( VDMATK ) = VGHAT3N
!    &               / ( 1.0 - EXP( -VGHAT3N * ( RA + RD3N ) ) )
      VDEP( VDMATK ) = 1.0 / ( RA + RD3N )

C accumulation mode
      SC3AI2 = (NUI * DCHAT3A) ** 2
      EIM    = MIN( STOKEAF * ESAC160, 1.0 )
      RD3A   = UTSCALI / ( SC3AI2**ONE3 + EIM )

!     VDEP( VDMACC ) = VGHAT3A
!    &               / ( 1.0 - EXP( -VGHAT3A * ( RA + RD3A ) ) )
      VDEP( VDMACC ) = 1.0 / ( RA + RD3A )

C coarse mode
      SC3CI2 = (NUI * DCHAT3C) ** 2
      EIM    = MIN( STOKECF * ESCO160, 1.0 )
      RD3C   = UTSCALI / ( SC3CI2**ONE3 + EIM )

!     VDEP( VDMCOR ) = VGHAT3C
!    &               / ( 1.0 - EXP( -VGHAT3C * ( RA + RD3C ) ) )
      VDEP( VDMCOR ) = 1.0 / ( RA + RD3C )

      END SUBROUTINE GETDEP_V

      end subroutine aero_depv
