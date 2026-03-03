
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
! $Header: /project/yoj/arc/CCTM/src/emis/emis/DUST_EMIS.F,v 1.7 2012/01/19 15:21:13 yoj Exp $

      module dust_emis

! Description:
! * Extracts selected landuse types from BELD01 and BELD03 and merges
! * the selections into a dust-related landuse array (ULAND).

!   Optionally, reads 3 gridded crop calendar file and calculates an
!   erodible agriculture land fraction. (cropcal)

! * Applies a predetermined removal fraction in and below canopy to
! * ULAND and determines a transport factor (TFB) for this regime.
! * = applies to tfbelow

! Function: 3d point source emissions interface to the chemistry-transport model

! Revision History:
! 16 Dec 10 J.Young: Adapting Daniel Tong`s work on windblown dust
! 21 Apr 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
! 11 May 11 D.Wong: incorporated twoway model implementation
!  8 Jul 11 J.Young: unified string lengths in character lists for compiler
!                    compatibility
!  8 Jun 12 J.Young: remove full character blank padding for GNU Fortran (GCC) 4.1.2
! 30 Sep 13 J.Young: corrected diag file units description; added snow cover adjustment;
!                    adjusted F/G (vertical to horizontal flux) ratio to be continuous
!                    for clay content > 0.2; convert volumetric soil moisture to
!                    gravimetric water content; corrected soil moisture factor (fmoit);
!                    use lwmask>0 rather than sltyp>0 (non-existent) for over water test

      use aero_data,  only: n_mode
      use const_data, only: pi
      use leaf_coms,  only: nvtyp

      implicit none

!Bob: Hardwire nstyp for now, but should remove references to soil classes
      integer, parameter :: nstyp = 12

! Sea-Salt Emissions Rates
      real, allocatable, save :: dustoutm( :,:,: )  ! mass emission rates [g/m**3/s]
      real, allocatable, save :: dustoutn( :,: )    ! number emission rates [1/m**3/s]
      real, allocatable, save :: dustouts( :,: )    ! surface-area emisrates [m2/m**3/s]

      public dustoutm, dustoutn, dustouts,
     &       dust_emis_init, get_dust_emis, dust_em
      private

      real, allocatable, save :: dust_em( : )  ! total dust emissions [g/m**3/s]

! updated values of mass fraction for "freshly emitted dust"
! based on Kok [PNAS, 2011] and Nabat et al. [ACP, 2012]
      real, parameter :: fracmj = 0.07  ! mass fraction assigned to accum mode
      real, parameter :: fracmk = 0.93  ! mass fraction assigned to coarse mode

! diam`s from FRACMJ,FRACMK-weighted 2 2-bin averages of geom means
      real, parameter :: dgvj = 1.3914  ! geom mean diam of accum mode [um]
      real, parameter :: dgvk = 5.2590  ! geom mean diam of coarse mode [um]
      real, parameter :: sigj = 2.0000  ! geom std deviation of accum mode flux
      real, parameter :: sigk = 2.0000  ! geom std deviation of coarse mode flux

! Factors for converting 3rd moment emission rates into number and 2nd moment
! emission rates.  (Diameters in [um] changed to [m] ) See Equations 7b and 7c
! of Binkowski & Roselle (2003)
      real, parameter :: l2sgj    = log( sigj ) ** 2
      real, parameter :: l2sgk    = log( sigk ) ** 2
      real, parameter :: factnumj = exp( 4.5 * l2sgj ) / dgvj ** 3 * 1.0e18
      real, parameter :: factnumk = exp( 4.5 * l2sgk ) / dgvk ** 3 * 1.0e18
      real, parameter :: factm2j  = exp( 0.5 * l2sgj ) / dgvj * 1.0e6
      real, parameter :: factm2k  = exp( 0.5 * l2sgk ) / dgvk * 1.0e6
      real, parameter :: factsrfj = pi * factm2j
      real, parameter :: factsrfk = pi * factm2k

      integer, parameter :: ndp = 4   ! number of soil texture type particle sizes

!---Mean mass median particle diameter (m) for each soil texture type
!   Chatenet et al. [Sedimentology,1996] and Menut et al. [JGR,2013]
      real, parameter :: dp( ndp ) = (/ 690.0E-6,   ! Coarse sand
     &                                  210.0E-6,   ! Fine-medium sand
     &                                  125.0E-6,   ! Silt
     &                                    2.0E-6 /) ! Clay

C---Soil texture: the amount of
C   1: Coarse sand, 2: Fine-medium sand, 3: Silt, 4: Clay
C   in each soil type [Kg/Kg]. from Menut et al. [JGR,2013]
         real, parameter :: soiltxt( ndp,nstyp ) = reshape (
     &           (/ 0.46,   0.46,   0.05,   0.03,     ! Sand
     &              0.41,   0.41,   0.18,   0.00,     ! Loamy Sand
     &              0.29,   0.29,   0.32,   0.10,     ! Sandy Loam
     &              0.00,   0.17,   0.70,   0.13,     ! Silt Loam
     &              0.00,   0.43,   0.39,   0.18,     ! Loam
     &              0.29,   0.29,   0.15,   0.27,     ! Sandy Clay Loam
     &              0.00,   0.10,   0.56,   0.34,     ! Silty Clay Loam
     &              0.00,   0.32,   0.34,   0.34,     ! Clay Loam
     &              0.00,   0.52,   0.06,   0.42,     ! Sandy Clay
     &              0.00,   0.06,   0.47,   0.47,     ! Silty Clay
     &              0.00,   0.22,   0.20,   0.58,     ! Clay
     &              0.00,   0.00,   0.00,   0.00 /),  ! Other
     &           (/ ndp,nstyp /), order = (/ 1,2 /) ) ! fill rows first

         real, parameter :: flxfac1( nstyp ) =
     &           (/ 5.8860e-05,     ! Sand
     &              5.8860e-05,     ! Loamy Sand
     &              5.2974e-05,     ! Sandy Loam
     &              5.2974e-05,     ! Silt Loam
     &              5.2974e-05,     ! Loam
     &              9.4176e-05,     ! Sandy Clay Loam
     &              5.2974e-05,     ! Silty Clay Loam
     &              5.2974e-05,     ! Clay Loam
     &              2.3544e-05,     ! Sandy Clay
     &              2.3544e-05,     ! Silty Clay
     &              2.3544e-05,     ! Clay
     &              0.00      /)    ! Other

         real, parameter :: flxfac2( nstyp ) =
     &           (/ 1.5215430,      ! Sand
     &              1.5215430,      ! Loamy Sand
     &              1.0758933,      ! Sandy Loam
     &              1.0758933,      ! Silt Loam
     &              1.0758933,      ! Loam
     &              1.0758933,      ! Sandy Clay Loam
     &              1.0758933,      ! Silty Clay Loam
     &              1.0758933,      ! Clay Loam
     &              0.1964303,      ! Sandy Clay
     &              0.1964303,      ! Silty Clay
     &              0.1964303,      ! Clay
     &              0.3402273 /)    ! Other

! Local Variables:

      real, save :: dustmode_dens( n_mode ) ! average modal density [kg/m**3]

      integer, parameter :: ndltyp = 3   ! number of desert land type beld1
      character( 16 ) :: vnmld( ndltyp ) = (/
     &                   'USGS_shrubland ',
     &                   'USGS_shrubgrass',
     &                   'USGS_sprsbarren' /)

! Number of soil types: For both WRF and MM5-PX met models, there are 16 types;
! the first 12 soil types are used and the rest lumped into Other.
      integer, parameter :: nsltyp = 13

! Variables for the windblown dust diagnostic file:
      logical, save      :: dustem_diag   ! flag for dustemis diagnostic file
      integer, parameter :: ndust_diag = 51 ! number of diagnostic output variables
!     real, allocatable, save :: dustbf( :,:,: ) ! diagnostic accumulate buffer

#ifdef verbose
      real, save         :: sdiagv( ndust_diag )  ! global sum of each diag output var
#endif

      character( 16 ) :: diagnm( ndust_diag ) = (/  ! output variable names
     &                   'QSHLD  ', 'QSHGR  ', 'QSPBA  ', 'QCROP  ',
     &                   'ETOTL  ',
     &                   'LSHLD  ', 'LSHGR  ', 'LSPBA  ', 'LCROP  ',
     &                   'TLDST  ',
     &                   'TFABV  ', 'TFINF  ',
     &                   'FMOIT  ', 'SD_EP  ', 'KV2HZ  ', 'WTMAX  ',
     &                   'ULDBR  ', 'ULDVG  ', 'ULDUB  ', 'ULDFR  ',
     &                   'USTAR  ' ,
     &                   'USHLD  ', 'USHGR  ', 'USPBA  ', 'UCROP  ',
     &                   'ASO4J  ', 'ANO3J  ', 'ACLJ   ', 'ANH4J  ',
     &                   'ANAJ   ', 'ACAJ   ', 'AMGJ   ', 'AKJ    ',
     &                   'APOCJ  ', 'APNCOMJ', 'AECJ   ',
     &                   'AFEJ   ', 'AALJ   ', 'ASIJ   ', 'ATIJ   ', 'AMNJ   ',
     &                   'AH2OJ  ', 'AOTHRJ ',
     &                   'ASO4K  ', 'ANO3K  ', 'ACLK   ', 'ASOIL  ',
     &                   'ANUMJ  ', 'ANUMK  ', 'ASRFJ  ', 'ASRFK  ' /)
!    &                   'STDEVJ ', 'STDEVK '  /)

! diagnostic output from duster (in PRMTR( :,:,1:29 ):
!  1 'QSHLD'  QAM( C,R,1 ) * 1.E6  ! g/s->ug/s   emissions for shrub landuse type
!  2 'QSHGR'  QAM( C,R,2 ) * 1.E6  ! g/s->ug/s      "      "  shrubgrass landuse type
!  3 'QSPBA'  QAM( C,R,3 ) * 1.E6  ! g/s->ug/s      "      "  barren landuse type
!  4 'QCROP'  QAM( C,R,4 ) * 1.E6  ! g/s->ug/s      "      "  cropland landuse type
!  5 'ETOTL'  QAT( C,R )           ! g/s         total for the 3 desert types and cropland
!  6 'LSHLD'  ELUS( C,R,1 )        ! shrub       erodible landuse fraction (%)
!  7 'LSHGR'  ELUS( C,R,2 )        ! shrubgrass     "        "       "
!  8 'LSPBA'  ELUS( C,R,3 )        ! barren         "        "       "
!  9 'LCROP'  ELUS( C,R,4 )        ! cropland       "        "       "
! 10 'TLDST'                       ! total desert fraction
!  =  LADUT( C,R,1 ) + LADUT( C,R,2 ) + LADUT( C,R,3 ) + AGLAND( C,R )
! 11 'TFABV'  TFA ( C,R )          ! transport fraction above canopy
! 12 'TFINF'  TFB ( C,R )          ! transport fraction in and below canopy
! 13 'FMOIT'  FMOIT( C,R )         ! soil moisture factor for threshold u*
! 14 'SD_EP'  SD_EP ( C,R )        ! soil->dust erodiblity potential
! 15 'KV2HZ'  KVH  ( C,R )         ! ratio of vertical to horizontal flux
! 16 'WTMAX'  WMAX ( C,R )         ! max adsorbed water fraction
! 17 'ULDBR'  ULAND( C,R,1 )       ! water & barren land coverage
! 18 'ULDVG'  ULAND( C,R,2 )       ! vege land coverage
! 19 'ULDUB'  ULAND( C,R,3 )       ! urban land coverage
! 20 'ULDFR'  ULAND( C,R,4 )       ! forest land coverage
! 21 'USTAR'  USTAR ( C,R )        ! u* (friction velocity)
! 22 'USHLD'  UTSTAR ( C,R,1 ) or 0.0  ! threshold u* for shrub
! 23 'USHGR'  UTSTAR ( C,R,2 ) or 0.0  !     "     "      shrubgrass
! 24 'USPBA'  UTSTAR ( C,R,3 ) or 0.0  !     "     "      barren
! 25 'UCROP'  UTSTAR ( C,R,4 ) or 0.0  !     "     "      cropland
! 26 'ASO4J'                       ! accumulation mode Sulfate
! 27 'ANO3J'                       ! accumulation mode Nitrate
! 28 'ACLJ'                        ! accumulation mode Chlorine
! 29 'ANH4J'                       ! accumulation mode Ammonium
! 30 'ANAJ'                        ! accumulation mode Sodium
! 31 'ACAJ'                        ! accumulation mode Calcium
! 32 'AMGJ'                        ! accumulation mode Magnesium
! 33 'AKJ'                         ! accumulation mode Potassium
! 34 'APOCJ'                       ! accumulation mode Organic Carbon
! 35 'APNCOMJ'                     ! accumulation mode Non-Carbon Organic Matter
! 36 'AECJ'                        ! accumulation mode Elemental Carbon
! 37 'AFEJ'                        ! accumulation mode Iron
! 38 'AALJ'                        ! accumulation mode Aluminum
! 39 'ASIJ'                        ! accumulation mode Silicon
! 40 'ATIJ'                        ! accumulation mode Titanium
! 41 'AMNJ'                        ! accumulation mode Manganese
! 42 'AH2OJ'                       ! accumulation mode Water
! 43 'AOTHRJ'                      ! accumulation mode Other
! 44 'ASO4K'                       ! coarse mode Sulfate
! 45 'ANO3K'                       ! coarse mode Nitrate
! 46 'ACLK'                        ! coarse mode Chlorine
! 47 'ASOIL'                       ! coarse mode total
! 48 'ANUMJ'                       ! accumulation mode number   [#/s]
! 49 'ANUMK'                       ! coarse mode number         [#/s]
! 50 'ASRFJ'                       ! accumulation mode surface area  [m**2/s]
! 51 'ASRFK'                       ! coarse mode surface area        [m**2/s]
!X 52 'STDEVJ'                     ! accumulation mode geometric standard deviation
!X 53 'STDEVK'                     ! coarse mode geometric standard deviation

! Module shared variables:
!     real,    allocatable, save :: sltyp ( :,: ) ! soil type (13 USDA categories)
!     real,    allocatable, save :: lwmask( :,: ) ! land-water mask (1=land, 0=water)
!     real,    allocatable, save :: agland( :,: ) ! agriculture land fraction
!      real,    allocatable, save :: tfb   ( : )   ! transport fraction in and below canopy
      real,    allocatable, save :: wmax  ( : )   ! max adsorb water fraction [%]
!      real,    allocatable, save :: kvh   ( : )   ! ratio of vertical flux / horizontal (k factor)
!     real,    allocatable, save :: sd_ep ( : )   ! soil->dust erodiblity potential
!     integer, allocatable, save :: uland ( : )   ! fraction of 4 landtype for dust removal [%]
!      real,    allocatable, save :: ladut ( : )   ! % of desertland (BELD01)

!      real, allocatable :: erod( : )
!      real, allocatable :: ht  ( : )

      real    :: sd_ep(  nstyp )
!      real    :: wmax (  nstyp )   ! max adsorb water fraction [%]
      real    :: z0s  (  nstyp )
      real    :: frd  (  nstyp )
      integer :: uland( 0:nvtyp)
      integer :: erod ( 0:nvtyp)
      real    :: hv   ( 0:nvtyp)

!     integer, save :: sdate, stime  ! scenario start date & time

      real, parameter :: eropot( 3 ) =  ! erodible potential of soil components
     &        (/ 0.08,   ! clay
     &           1.00,   ! silt
     &           0.12 /) ! sand

      integer, save :: logdev

      real, parameter :: mv       = 0.16
      real, parameter :: sigv     = 1.45
      real, parameter :: betav    = 202.0
      real, parameter :: sigv_mv  = sigv  * mv
      real, parameter :: betav_mv = betav * mv

      real, parameter :: mb       = 0.5
      real, parameter :: sigb     = 1.0
      real, parameter :: betab    = 90.0
      real, parameter :: sigb_mb  = sigb  * mb
      real, parameter :: betab_mb = betab * mb

      ! Removal fraction in and below canopy
      real, parameter :: fw( 4 ) = (/ 0.00,   ! water/barren
     &                                0.10,   ! grass&crop
     &                                0.50,   ! building
     &                                0.95 /) ! forest

!---Roughness density for solid elements
! from Darmenova et al. [JGR,2009] and Xi and Sokolik [JGR,2015]
      real, parameter :: lambdab( 4 ) =
     &           (/ 0.03,     ! shrubland
     &              0.04,     ! shrubgrass
     &              0.0001,   ! barrenland
     &              0.15 /)   ! cropland

!---Compound for computational efficiency
      real, parameter :: hb_lambdab( 4 ) =
     &           (/ 6.0e-04,   ! shrubland
     &              8.0e-04,   ! shrubgrass
     &              2.0e-06,   ! barrenland
     &              3.0e-03 /) ! cropland

      CONTAINS

!=======================================================================

         function dust_emis_init( jdate, jtime ) result( success )

         use aero_data        ! aerosol species definitions
         use const_data
         use utilio_defn

         use mem_land,   only: land, mland, nzg

! Arguments:
         integer, intent( in ) :: jdate   ! current model date, coded YYYYDDD
         integer, intent( in ) :: jtime   ! current model time, coded HHMMSS
         logical               :: success

! External Functions:
         integer, external :: setup_logdev

! Local variables:
!        character( 16 ) :: ctm_dustem_diag = 'CTM_DUSTEM_DIAG' ! env var for
                                                            ! diagnostic file
!        character( 16 ) :: ctm_erode_agland = 'CTM_ERODE_AGLAND' ! env var to
                                                       ! use erodible cropland
          character( 16 ) :: pname = 'DUST_EMIS_INIT'
!         character( 16 ) :: vname
!         character( 80 ) :: vardesc
          character( 120 ) :: xmsg = ' '

!        logical :: erode_agland = .true. ! default
         integer :: status
         integer :: i, j, n, idx, iland
         real    :: sumsplit, sumfrac

!         interface

!            subroutine cropcal ( jdate, jtime, agland )
!            integer, intent( in )  :: jdate, jtime
!            real,    intent( out ) :: agland( :,: )
!            end subroutine cropcal
!
!            subroutine tfbelow ( jdate, jtime, uland, tfb )
!            integer, intent( in )  :: jdate, jtime
!            real,    intent( out ) :: uland( :,:,: )
!            real,    intent( out ) :: tfb( :,: )
!            end subroutine tfbelow

!         end interface

!-----------------------------------------------------------------------

         logdev = setup_logdev()
         success = .true.

         allocate ( dustoutm( ndust_spc,n_mode,mland ),
     &              dustoutn(           n_mode,mland ),
     &              dustouts(           n_mode,mland ), stat = status )
         if ( status .ne. 0 ) then
            xmsg = '*** DUSTOUTM, DUSTOUTN or DUSTOUTS memory allocation failed'
            call m3warn ( pname, jdate, jtime, xmsg )
            success = .false.; return
         end if

! Allocate emissions array
         allocate( dust_em( mland ), stat = status )
         if ( status .ne. 0 ) then
            xmsg = '*** Failure allocating DUST_EM'
            call m3warn( pname, jdate, jtime, xmsg )
            success = .false.; return
         end if

! Allocate private arrays
!         allocate( agland( ncols,nrows ),
!     &             wmax  ( ncols,nrows ),
!     &             kvh   ( ncols,nrows ),
!     &             sd_ep ( ncols,nrows ),
!     &             sltyp ( ncols,nrows ),
!     &             lwmask( ncols,nrows ),
!     &             tfb   ( ncols,nrows ),
!     &             uland ( ncols,nrows,4 ),
!     &             ladut ( ncols,nrows,ndltyp ), stat = status )

         allocate( wmax ( mland ))

         ! Convert residual soil moisture from volumetric (m**3/m**3) to
         ! gravimetric (Kg/Kg)

         do iland = 2, mland
            wmax(iland) = land%wresid_vg(nzg,iland) * 1000.0 / land%bulkdens_drysoil(nzg,iland)
         enddo

!        allocate( kvh  ( mland ), stat = status )
!     &             tfb  ( mland ),
!     &             ladut( mland ), stat = status )

!         allocate( erod( mland) )
!         allocate( ht  ( mland) )

!         if ( status .ne. 0 ) then
!            xmsg = '*** Failure allocating AGLAND, WMAX, KVH, SD_EP,'
!     &           // ' SLTYP, LWMASK, TFB, ULAND, or LADUT'
!            call m3warn( pname, jdate, jtime, xmsg )
!            success = .false.; return
!         end if

!        agland = 0.0   ! array assignment
!         kvh    = 0.0   ! array assignment
!        sd_ep  = 0.0   ! array assignment

! Get env var for diagnostic output
         dustem_diag = .false.         ! default
!         vardesc = 'Flag for writing the windblown dust emission diagnostic file'
!         dustem_diag = envyn( ctm_dustem_diag, vardesc, dustem_diag, status )
!         if ( status .ne. 0 ) write( logdev,'( 5x, a )' ) vardesc
!         if ( status .eq. 1 ) then
!            xmsg = 'Environment variable improperly formatted'
!            call m3warn( pname, jdate, jtime, xmsg )
!            success = .false.; return
!         else if ( status .eq. -1 ) then
!            xmsg = 'Environment variable set, but empty ... Using default:'
!            write( logdev,'( 5x, a, i9 )' ) xmsg, jtime
!         else if ( status .eq. -2 ) then
!            xmsg = 'Environment variable not set ... Using default:'
!            write( logdev,'( 5x, a, i9 )' ) xmsg, jtime
!         end if

!         if ( dustem_diag ) then    ! Open the sea-salt emission dignostic file
!
!C Allocate diagnostic emissions arrays
!            allocate( dustbf( ndust_diag,ncols,nrows ), stat = status )
!            if ( status .ne. 0 ) then
!               xmsg = '*** Failure allocating DUSTBF'
!               call m3warn( pname, jdate, jtime, xmsg )
!               success = .false.; return
!            end if
!
!            sdate = envint( 'CTM_STDATE', 'Scenario Start (YYYYJJJ)', 0, status )
!            stime = envint( 'CTM_STTIME', 'Scenario Start (HHMMSS)', 0, status )
!
!            if ( mype .eq. 0 )
!     &         call opdust_emis ( sdate, stime, tstep, ndust_diag, diagnm )
!
!         end if

!         call subhfile ( dust_lu_1, gxoff, gyoff,
!     &                   strtcol_01, endcol_01, strtrow_01, endrow_01 )

! Get desert land (fraction) data

!         do i = 1, ndltyp
!            if ( .not. xtract3( dust_lu_1, vnmld( i ),
!     &                          1,1, strtrow_01,endrow_01, strtcol_01,endcol_01,
!     &                          0, 0, ladut( 1,1,i ) ) ) then
!               xmsg = 'Could not read ' // trim( vnmld( i ) ) // ' from DUST_LU_1'
!               call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!            end if
!         end do
!
!         call subhfile ( met_cro_2d, gxoff, gyoff,
!     &                   strtcolmc2, endcolmc2, strtrowmc2, endrowmc2 )
!
!         call subhfile ( grid_cro_2d, gxoff, gyoff,
!     &                   strtcolgc2, endcolgc2, strtrowgc2, endrowgc2 )

! Get soil type data
!        if ( metmod .eq. 'mm5' ) then
!            vname = 'SLTYP'
!        else if (metmod == 'NAM') then
!           vname = 'sotyp'
!        end if
!         if ( .not. interpx( met_cro_2d, vname, pname,
!     &                       strtcolmc2,endcolmc2, strtrowmc2,endrowmc2, 1,1,
!     &                       jdate, jtime, sltyp ) ) then
!            xmsg = 'Could not read ' // trim( vname ) // ' from MET_CRO_2D'
!            call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!         end if
!
!         vname = 'LWMASK'
!         if ( .not. interpx( grid_cro_2d, vname, pname,
!     &                       strtcolgc2,endcolgc2, strtrowgc2,endrowgc2, 1,1,
!     &                       jdate, jtime, lwmask ) ) then
!            xmsg = 'Could not read ' // trim( vname ) // ' from GRID_CRO_2D'
!            call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!         end if

! Get grid cell area
! note: crude estimate is made for LAT/LONG coordinate systems
!         if ( gdtyp_gd .eq. latgrd3 ) then
!            dx1 = dg2m * real( xcell_gd, 4 ) ! in m
!            dx2 = dg2m * real( ycell_gd, 4 )
!     &          * cos( pi180*( yorig_gd + real( ycell_gd, 4 )
!     &          * float( gl_nrows/2 ) ) ) ! in m
!         else
!            dx1 = real( xcell_gd, 4 )        ! in m
!            dx2 = real( ycell_gd, 4 )        ! in m
!         end if
!         cellarea = dx1 * dx2

! Get env var for erodible agriculture land fraction
!         erode_agland = .false.         ! default
!         vardesc = 'Flag for calculating erodible agriculture land fraction'
!         erode_agland = envyn( ctm_erode_agland, vardesc, erode_agland, status )
!         if ( status .ne. 0 ) write( logdev,'( 5x, a )' ) vardesc
!         if ( status .eq. 1 ) then
!            xmsg = 'Environment variable improperly formatted'
!            call m3warn( pname, jdate, jtime, xmsg )
!            success = .false.; return
!         else if ( status .eq. -1 ) then
!            xmsg = 'Environment variable set, but empty ... Using default:'
!            write( logdev,'( 5x, a, i9 )' ) xmsg, jtime
!         else if ( status .eq. -2 ) then
!            xmsg = 'Environment variable not set ... Using default:'
!            write( logdev,'( 5x, a, i9 )' ) xmsg, jtime
!         end if
!
!         if ( erode_agland ) then
!            call cropcal ( sdate, stime, agland )
!            do r = 1, my_nrows
!            do c = 1, my_ncols
!               if ( agland( c,r ) .lt. 0.0 .or. agland( c,r ) .gt. 100.0 ) then
!                  xmsg = '*** ERROR in AGLAND'
!                  call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!               end if
!            end do
!            end do
!         end if

         ! Flags for dust land use category:
         ! 0 = nodust; 1 = barren, 2 = grass, 3 = urban, 4 = forest
         uland( 0) = 0   ! Ocean
         uland( 1) = 0   ! Lakes, rivers, streams
         uland( 2) = 0   ! Ice cap/glacier
         uland( 3) = 1   ! Desert, bare soil
         uland( 4) = 0   ! Evergreen needleleaf tree
         uland( 5) = 0   ! Deciduous needleleaf tree
         uland( 6) = 0   ! Deciduous broadleaf tree
         uland( 7) = 0   ! Evergreen broadleaf tree
         uland( 8) = 2   ! Short grass
         uland( 9) = 2   ! Tall grass
         uland(10) = 1   ! Semi-desert
         uland(11) = 0   ! Tundra
         uland(12) = 2   ! Evergreen shrub
         uland(13) = 2   ! Deciduous shrub
         uland(14) = 4   ! Mixed woodland
         uland(15) = 2   ! Crop/mixed farming, C3 grassland
         uland(16) = 0   ! Irrigated crop
         uland(17) = 0   ! Bog or marsh
         uland(18) = 4   ! Wooded grassland
         uland(19) = 3   ! Urban and built up
         uland(20) = 0   ! Wetland evergreen broadleaf tree
         uland(21) = 2   ! Deforested (Amazon - Medvigy)

         ! Flags for erodable vegetation type
         ! 0 = nodust; 1 = shrub, 2 = grass, 3 = barren, 4 = ag
         erod( 0) = 3   ! Ocean
         erod( 1) = 3   ! Lakes, rivers, streams
         erod( 2) = 3   ! Ice cap/glacier
         erod( 3) = 3   ! Desert, bare soil
         erod( 4) = 2   ! Evergreen needleleaf tree
         erod( 5) = 2   ! Deciduous needleleaf tree
         erod( 6) = 2   ! Deciduous broadleaf tree
         erod( 7) = 2   ! Evergreen broadleaf tree
         erod( 8) = 2   ! Short grass
         erod( 9) = 2   ! Tall grass
         erod(10) = 1   ! Semi-desert
         erod(11) = 1   ! Tundra
         erod(12) = 1   ! Evergreen shrub
         erod(13) = 1   ! Deciduous shrub
         erod(14) = 2   ! Mixed woodland
         erod(15) = 4   ! Crop/mixed farming, C3 grassland
         erod(16) = 4   ! Irrigated crop
         erod(17) = 2   ! Bog or marsh
         erod(18) = 2   ! Wooded grassland
         erod(19) = 2   ! Urban and built up
         erod(20) = 2   ! Wetland evergreen broadleaf tree
         erod(21) = 2   ! Deforested (Amazon - Medvigy)

         ! Flags for erodable vegetation height
         hv( 0) = 0.10  ! Ocean
         hv( 1) = 0.10  ! Lakes, rivers, streams
         hv( 2) = 0.10  ! Ice cap/glacier
         hv( 3) = 0.10  ! Desert, bare soil
         hv( 4) = 1.00  ! Evergreen needleleaf tree
         hv( 5) = 1.00  ! Deciduous needleleaf tree
         hv( 6) = 1.00  ! Deciduous broadleaf tree
         hv( 7) = 1.00  ! Evergreen broadleaf tree
         hv( 8) = 0.08  ! Short grass
         hv( 9) = 0.15  ! Tall grass
         hv(10) = 0.10  ! Semi-desert
         hv(11) = 0.10  ! Tundra
         hv(12) = 0.15  ! Evergreen shrub
         hv(13) = 0.15  ! Deciduous shrub
         hv(14) = 1.00  ! Mixed woodland
         hv(15) = 0.20  ! Crop/mixed farming, C3 grassland
         hv(16) = 0.20  ! Irrigated crop
         hv(17) = 0.20  ! Bog or marsh
         hv(18) = 1.00  ! Wooded grassland
         hv(19) = 0.50  ! Urban and built up
         hv(20) = 1.00  ! Wetland evergreen broadleaf tree
         hv(21) = 0.10  ! Deforested (Amazon - Medvigy)

      do j = 1, nstyp-1
         sd_ep( j ) = eropot( 1 ) *  soiltxt( 4,j )
     &              + eropot( 2 ) *  soiltxt( 3,j )
     &              + eropot( 3 ) * (soiltxt( 1,j ) + soiltxt( 2,j ))

         z0s  ( j ) = sum( dp(:) * soiltxt(:,j) ) / 30.

         frd  ( j ) = 1.0 / log( 0.7 * (122.55 / z0s(j))**0.8 )

      enddo

! Calculate modal average dust particle densities (accum and coarse modes) [ kg/m**3 ]
! The following works because the dust_spc`s are a fixed split of the total emitted
! mass.
         dustmode_dens( 1 ) = 0.0
         do n = 2, n_mode
            sumsplit = 0.0; sumfrac = 0.0
            do i = 1, ndust_spc
               sumsplit = sumsplit + dust_spc( i )%spcfac( n ) ! should = 1.0
               idx = findAero( dust_spc( i )%name( n ), .true. )
               sumfrac = sumfrac + dust_spc( i )%spcfac( n ) * aerospc_spvol( idx )
            end do
            dustmode_dens( n ) = sumsplit / sumfrac
         end do

#ifdef verbose
         write( logdev,* ) ' '
         write( logdev,* ) '    l2sgj,l2sgk:         ', l2sgj, l2sgk
         write( logdev,* ) '    factnumj,factnumk:   ', factnumj, factnumk
         write( logdev,* ) '    factm2j,factm2k:     ', factm2j, factm2k
         write( logdev,* ) '    factsrfj,factsrfk:   ', factsrfj, factsrfk
         write( logdev,* ) '    modal avg dens(j/k): ', dustmode_dens( 2 ), dustmode_dens( 3 )
         write( logdev,* ) ' '
#endif

         end function dust_emis_init

!=======================================================================

!         subroutine opdust_emis ( jdate, jtime, tstep, ndust_diag, wrdust_spc )
!
!C   27 Dec 10 J.Young: initial
!
!         use grid_conf           ! horizontal & vertical domain specifications
!         use utilio_defn
!
!         implicit none
!
!         include SUBST_FILES_ID  ! file name parameters
!
!C Arguments:
!         integer      jdate      ! current model date, coded YYYYDDD
!         integer      jtime      ! current model time, coded HHMMSS
!         integer      tstep      ! output time step
!         integer      ndust_diag
!         character( 16 ) :: wrdust_spc( ndust_diag )
!
!C Local variables:
!         character( 16 ) :: pname = 'OPDUST_EMIS     '
!         character( 96 ) :: xmsg = ' '
!
!         integer      v, l       ! loop induction variables
!
!C-----------------------------------------------------------------------
!
!C Try to open existing file for update
!         if ( .not. open3( ctm_dust_emis_1, fsrdwr3, pname ) ) then
!            xmsg = 'Could not open CTM_DUST_EMIS_1 for update - '
!     &           // 'try to open new'
!            call m3mesg( xmsg )
!
!C Set output file characteristics based on COORD.EXT and open diagnostic file
!            ftype3d = grdded3
!            sdate3d = jdate
!            stime3d = jtime
!            tstep3d = tstep
!            call nextime( sdate3d, stime3d, tstep3d ) !  start the next hour
!
!            nvars3d = ndust_diag
!            ncols3d = gl_ncols
!            nrows3d = gl_nrows
!            nlays3d = 1
!            nthik3d = 1
!            gdtyp3d = gdtyp_gd
!            p_alp3d = p_alp_gd
!            p_bet3d = p_bet_gd
!            p_gam3d = p_gam_gd
!            xorig3d = xorig_gd
!            yorig3d = yorig_gd
!            xcent3d = xcent_gd
!            ycent3d = ycent_gd
!            xcell3d = xcell_gd
!            ycell3d = ycell_gd
!            vgtyp3d = vgtyp_gd
!            vgtop3d = vgtop_gd
!!           vgtpun3d = vgtpun_gd ! currently, not defined
!            do l = 1, nlays3d + 1
!               vglvs3d( l ) = vglvs_gd( l )
!            end do
!            gdnam3d = grid_name  ! from HGRD_DEFN
!
!            do v = 1, nvars3d
!               vtype3d( v ) = m3real
!               vname3d( v ) = wrdust_spc( v )
!               if ( vname3d( v )(1:4) .eq. 'ANUM' ) then
!                  units3d( v ) = 'number/m**3/s'
!               else if ( vname3d( v )(1:4) .eq. 'ASRF' ) then
!                  units3d( v ) = 'm**2/m**3/s'
!               else
!                  units3d( v ) = 'g/m**3/s'
!               end if
!               vdesc3d( v ) = 'hourly avg ' // trim( vname3d( v ) )
!     &                     // ' wbdust emis rate per grid cell volume'
!            end do
!
!            fdesc3d( 1 ) = 'hourly layer-1 windblown dust emission rates'
!            do l = 2, mxdesc3
!               fdesc3d( l ) = ' '
!            end do
!
!C Open windblown dust emissions diagnostic file
!            if ( .not. open3( ctm_dust_emis_1, fsnew3, pname ) ) then
!               xmsg = 'Could not create the CTM_DUST_EMIS_1 file'
!               call m3exit( pname, sdate3d, stime3d, xmsg, xstat1 )
!            end if
!
!         end if
!
!         return
!
!         end subroutine opdust_emis

!=======================================================================

         subroutine get_dust_emis( )

         use mem_land,  only: mland, omland
         use mem_para,  only: myrank
         use mem_sfcg,  only: itab_wsfc
         use misc_coms, only: iparallel

         implicit none

         integer :: iland, iwsfc

! Initialize windblown dust diagnostics output buffer
!        if ( dustem_diag .and. wstep .eq. 0 ) then
!            dustbf = 0.0   ! array assignment
!#ifdef verbose
!            sdiagv = 0.0   ! array assignment
!#endif
!         end if

         !$omp parallel do
         do iland = 2, mland
            iwsfc = iland + omland

            ! Skip this cell if running in parallel and cell rank is not MYRANK
            if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

            call get_dust_emis1(iland)
         enddo
         !$omp end parallel do

         end subroutine get_dust_emis

!=======================================================================

         subroutine get_dust_emis1( iland )

         use utilio_defn
         use aero_data
         use const_data
         use mem_land,   only: land, omland, nzg
         use leaf_coms,  only: wcap_min
         use mem_grid,   only: volti, gdz_abov
         use mem_sfcg,   only: itab_wsfc, sfcg
         use mem_megan,  only: pfts
         use consts_coms, only: vonk

! Arguments:

         integer, intent(in) :: iland

! Parameters:
         real, parameter :: sc  = 32.0     ! SENSIT-based scaling parameter
!        real, parameter :: z0s = 0.00001  ! Surface roughness for dust [m],
                                           ! Ref. Marticorena et al. JGR 1997

         real, parameter :: amen = 2.8     ! values are from 0 to 3.5, mean 2.8
!        real, parameter :: amen = 32.0    ! values adjusted based on Dale's MNORT measurements
         real, parameter :: radg = 122.87  ! g*s**2/m**4 122.87->120.8
!        real, parameter :: fdsa = 0.5     ! factor of dust source area 1/.1/.01
         real, parameter :: fdsa = 1.0     ! factor of dust source area 1/.1/.01
         real, parameter :: gmad = 0.78    ! Factor to convert geometric to aerodynamic diameter
         real, parameter :: cfac = amen * radg * fdsa * gmad  ! combined factor for qam
         real, parameter :: gpkg = 1.0e03  ! g/kg

!        character( 16 ) :: pname = 'GET_DUST_EMIS   '
!        character( 16 ) :: vname
!        character( 96 ) :: xmsg
         integer :: n, v, kw, indx, iw, j, lc, iwsfc, ja

!         integer, save   :: wstep = 0      ! local write counter
!         integer         :: mdate, mtime   ! diagnostic file write date&time

         real :: ku
         real :: fmoit  ! factor of soil moisture on u*t
         real :: ustar  ! friction velocity [m/s]
!        real :: utstar ! threshold U* [m/s]
         real :: soimt  ! soil moisture (m**3/m**3)
!        real :: utem   ! U term = U*(U**2-Ut**2) [(m/s)**3]
         real :: tfa    ! transport fraction above canopy
         real :: tfb    ! transport fraction below canopy
         real :: ferod  ! fraction of cell area that is erodable
         real :: vflux   ! vertical flux
         real :: kvh    ! vertical to horizontal flux ratio
         real :: slmsts_top
         real :: vels, vels10, dens

! automatic arrays:
!        real       :: fmoit ( ncols,nrows )
!        real       :: rn    ( ncols,nrows )        ! nonconvective preciptation
!        real       :: rc    ( ncols,nrows )        ! convective preciptation
!        real       :: snow  ( ncols,nrows )        ! grid cell fraction snow-covered
!        real       :: ustar ( ncols,nrows )        ! friction velocity [m/s]
!        real       :: utstar( ncols,nrows,ndltyp+1 ) ! threshold U* [m/s]
!        real       :: soimt ( ncols,nrows )        ! soil moisture (m**3/m**3)
!        real       :: utem  ( ncols,nrows,ndltyp+1 ) ! U term = U*(U**2-Ut**2) [(m/s)**3]
!        real       :: qam   ( ncols,nrows,ndltyp+1 ) ! emis for landuse type [g/m**/s]
!        real       :: elus  ( ncols,nrows,ndltyp+1 ) ! Erodible landuse fraction (0~100)
!        real       :: wind10( ncols,nrows )        ! wind speed at 10m altitude [m/s]
!        real       :: zruf  ( ncols,nrows )        ! surface roughness length [m]
!        real       :: tfa   ( ncols,nrows )        ! transport fraction above canopy
!        real       :: wrbuf ( ncols,nrows )        ! diagnositc write buffer

!        real       :: diagv ( ndust_diag )         ! diagnostic output variables

         real       :: edust( n_mode )  ! mass emis rate [g/s] per mode (only accum & coarse)
         real       :: m3j              ! 3rd moment accumulation (J) mode emis rates [m3/m3/s]
         real       :: m3k              ! 3rd moment coarse mode (K) emis rates [m3/m3/s]


         real       :: vegfree, lambda, lambdav, fruf, ftot, vegheight, z0

! Threshold Friction Velocity [m/s] for 13 soil types and 3 land use types, which are:
!    1 -shrubland,
!    2 -shrubgrass land, and
!    3 -barren land
! Local (surface) Threshold Friction Velocity from field measurements
!   Ref. Gillette et al., JGR, 1980 for desert soils (J=2).
!   Ref. Gillette, JGR, 1988 for Loose Agr. Soils (J=3).
! The following soil types were not measured for desert land (we chose to use
! agr. data): Sandy Clay Loam, Clay Loam, Sandy Clay, and Silty Clay.
! Modified values compatiable with both MM5 & NAM.
! There is no measurement of this value for Silt. The values for Silt are chosen from
! Silty Loam since the soil composition is close.
! Other includes all types higher than 12. The values of Other are too high to allow
! any dust emission.
         real, parameter :: uthfv( nsltyp,3 ) = reshape (
     &           (/ 0.80,   0.42,   0.28,     ! Sand
     &              1.00,   0.51,   0.34,     ! Loamy Sand
     &              1.40,   0.66,   0.29,     ! Sandy Loam
     &              1.70,   0.34,   1.08,     ! Silt Loam
     &              1.70,   0.34,   1.08,     ! Silt
     &              1.70,   0.49,   0.78,     ! Loam
     &              1.70,   0.78,   0.78,     ! Sandy Clay Loam
     &              1.70,   0.33,   0.64,     ! Silty Clay Loam
     &              1.70,   0.71,   0.71,     ! Clay Loam
     &              1.70,   0.71,   0.71,     ! Sandy Clay
     &              1.70,   0.56,   0.56,     ! Silty Clay
     &              1.70,   0.78,   0.54,     ! Clay
     &              9.99,   9.99,   9.99 /),  ! Other
     &           (/ nsltyp,3 /), order = (/ 2,1 /) )

         real, parameter :: dmsk( 4 ) = ! land use type desert fraction
     &           (/ 0.50,    ! shrubland
     &              0.25,    ! shrubgrass
     &              0.75,    ! barrenland
     &              0.75 /)  ! cropland

! Based on Slinn 1982, fugitive dust in coarse mode, Vd average from
! PM2~PM10; 3 wind velocities: 10, 5, 1 m/s; 4 canopy characteristics:
! water and lamda = 1.0, 3.5, 5.0

         real, parameter :: vd( 4,3 ) = reshape ( ! deposition velocity [m/s]
     &           (/ 0.0051,    ! natural water Vd(1,1), WindSpeed = 10 m/s
     &              0.0152,    ! lamda = 2.0   Vd(2,1), WindSpeed = 10 m/s
     &              0.0268,    ! lamda = 3.5   Vd(3,1), WindSpeed = 10 m/s
     &              0.0382,    ! lamda = 5.0   Vd(4,1), WindSpeed = 10 m/s
     &              0.0018,    ! natural water Vd(1,2), WindSpeed =  5 m/s
     &              0.0056,    ! lamda = 2.0   Vd(2,2), WindSpeed =  5 m/s
     &              0.0099,    ! lamda = 3.5   Vd(3,2), WindSpeed =  5 m/s
     &              0.0141,    ! lamda = 5.0   Vd(4,2), WindSpeed =  5 m/s
     &              0.0018,    ! natural water Vd(1,3), WindSpeed =  1 m/s
     &              0.0020,    ! lamda = 2.0   Vd(2,3), WindSpeed =  1 m/s
     &              0.0021,    ! lamda = 3.5   Vd(3,3), WindSpeed =  1 m/s
     &              0.0030 /), ! lamda = 5.0   Vd(4,3), WindSpeed =  1 m/s
     &           (/ 4,3 /), order = (/ 1,2 /) )

!         interface
!            subroutine tfabove ( ustar, swind, uland, tfa )
!            real,    intent( in )  :: ustar( :,: )
!            real,    intent( in )  :: swind( :,: )
!            real,    intent( in )  :: uland( :,:,: )
!            real,    intent( out ) :: tfa( :,: )
!            end subroutine tfabove
!         end interface
!
!#ifdef verbose
!         integer dryhit
!         integer dusthit
!#endif

! Convert layer-average U* into surface local U*
!    -> U*t has to be surface values from measurements.

!         do r = 1, my_nrows
!         do c = 1, my_ncols
!!           ustar( c,r ) = ustar( c,r ) * log( zi ) / log( zi * zruf( c,r ) / z0s )
!            zi = sc * zruf( c,r )
!            ustar( c,r ) = ustar( c,r ) * log( zi / zruf( c,r ) ) / log( zi / z0s )
!         end do
!         end do
!
! calculate transport factor above the canopy
!         call tfabove ( ustar, wind10, uland, tfa )

!         vname = 'SOIM1'
!!        vname = 'SOILW'   ??? 'NAM' ???
!         if ( .not. interpx( met_cro_2d, vname, pname,
!     &                       strtcolmc2,endcolmc2, strtrowmc2,endrowmc2, 1,1,
!     &                       jdate, jtime, soimt ) ) then
!            xmsg = 'Could not read '// trim( vname ) // ' from MET_CRO_2D '
!            call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!         end if

!-----------------------------------------------------------------------
!  formula of Draxler & Gillette (2001) Atmos. Environ.
!      F   =  K A (r/g) U * ( U*^2 - Ut*^2 )
!  where
!      F   = emission flux                         [g/m**2-s]
!      K   = constant 2.0E-04                      [1/m]
!      A   = 0~3.5  mean = 2.8
!      r/g = ratio of air density to gravity=120.8 [g-s**2/m**4]
!      U*  = friction velocity                     [m/s]
!      Ut* = threshold friction velocity           [m/s]
!-----------------------------------------------------------------------

!#ifdef verbose
!         dryhit = 0
!         dusthit = 0
!#endif

            iwsfc = iland + omland

            dust_em     (iland) = 0.0
            dustoutn  (:,iland) = 0.0
            dustouts  (:,iland) = 0.0
            dustoutm(:,:,iland) = 0.0

            ! Skip if no erodable land use type

            lc = sfcg%leaf_class(iwsfc)
            if (uland(lc) == 0) return

            ! Skip if grid cell is snow or water covered

            if ( land%sfcwater_mass(1,iland) >= wcap_min ) return

            ! Skip if top of soil is near saturation

            slmsts_top = land%wsat_vg(nzg,iland)
            if ( land%soil_water(nzg,iland) > 0.6 * slmsts_top ) return

            ! Skip if precipitating

            if ( sfcg%pcpg(iwsfc) > wcap_min ) return

            ! Skip if no bare ground

            ferod =  pfts( 0,iland)                        ! barren
     &            + (pfts(15,iland) + pfts(16,iland)) * 0.7  ! crop

            if (ferod < 0.01) return

            ! Convert layer-average U* into surface local U*
            ! -> U*t has to be surface values from measurements.

            ! zi = sc * sfcg%rough(iwsfc)
            ! ustar = sfcg%ustar(iwsfc) * log( zi / sfcg%rough(iwsfc) ) / log( zi / z0s )

            dens  = sfcg%rhos (iwsfc)
            ustar = sfcg%ustar(iwsfc)
            vels  = sfcg%vels (iwsfc)

! ???? possiblity if rain < 1/100 in. (1 in. = 2.540 cm) and not over water (SLTYP = 0)
! Dust possiblity if rain < 1/100 in. (1 in. = 2.540 cm) and not over water (SLTYP = 14)
!           if ( ( ( rn( c,r ) + rc( c,r ) ) .le. 0.0254 ) .and. ! rn, rc = [cm]
!??  &           ( sltyp( c,r ) .gt. 0.0 ) ) then
!    &           ( sltyp( c,r ) .lt. 14.0 ) ) then
! Dust possiblity if not over water and rain < 1/100 in. (1 in. = 2.540 cm)
!            if ( ( lwmask( c,r ) .gt. 0.0 ) .and.
!     &           ( rn( c,r ) + rc( c,r ) ) .le. 0.0254 ) then ! rn, rc = [cm]
!
!               j = nint( sltyp( c,r ) )

            ! Soil textural class
            ! PX and OLAM combine "silt" with "silt loam"

            !j = land%ntext_soil(nzg,iland)

!           ! Convert volumetric soil moisture (m**3/m**3) to
            ! gravimetric water content (Kg/Kg)

            soimt = land%soil_water(nzg,iland) * 1000.0 / land%bulkdens_drysoil(nzg,iland)

            ! Soil moisture effect on U*t

            if ( soimt <= wmax(iland) ) then
               fmoit = 1.0
            else
               fmoit = sqrt( 1.0 + 27.7195 * ( soimt - wmax(iland) )**0.68 )
            end if

            ! Adjust U*t based on soil moisture ---

!
!                  do m = 1, 3   ! landuse type
!                     elus( c,r,m ) = ladut( c,r,m ) * dmsk( m )  ! desert land [%]
!                  end do
!                  elus( c,r,4 ) = agland( c,r ) * dmsk( 4 )      ! crop land [%]

!                  do m = 1, 4   ! landuse type

!                     if ( elus( c,r,m ) .gt. 100.0 .or. elus( c,r,m ) .lt. 0.0 ) then
!                        write( xmsg,2009 ) elus( c,r,m ), c, r, m
!                        call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!                     end if


            n = erod(lc)

            vegfree   = max(0.001, min(pfts(0,iland), 1.0))
            lambdav   = -0.35 * log( vegfree ) ! Shao et al. [Aus. J. Soil Res.,1996]

!            fruf2  = ( 1.0 - sigv_mv  * lambdav )
!     &             * ( 1.0 + betav_mv * lambdav )
!     &             * ( 1.0 - sigb_mb  * lambdab( n ) / vegfree )
!     &             * ( 1.0 + betab_mb * lambdab( n ) / vegfree )
!
!            fruf = sqrt(max(fruf2, 1.0))

            lambda    = lambdab( n ) + lambdav
            vegheight = ( hb_lambdab( n ) + hv(lc) * lambdav ) / lambda
!           vegheight = ( hb_lambdab( n ) + veg_ht(lc) * lambdav ) / lambda

!            if ( lambda .le. 0.2 ) then
!               z0 = 0.96  * ( lambda ** 1.07 )  * vegheight
!            else
!               z0 = 0.083 * ( lambda ** ( -0.46 ) ) * vegheight
!            end if

            if ( lambda < 0.0408 ) then
!              z0 = vegheight * 10 ** ( 1.31 * log10(lambda) + 0.66 )
               z0 = vegheight * exp( 1.31 * log(lambda) + 1.52 )
            else
!              z0 = vegheight * 0.3135
               z0 = vegheight * 0.0692
            endif

            z0 = max(z0, z0s(j))

            fruf = 1.0 / ( 1.0 - min( frd(j) * log( z0 / z0s(j) ), 0.99 ) )

 ! BOB: This section of code references volti(kw,iw), so it needs to know the ATM cell.
 ! However, everything else about this subroutine is defined by the land cell itself.
 ! I recommend that dust_em be defined as kg/m^2 instead of kg/m^3 so that it would be
 ! independent of the ATM cell.

       !!   do ja = 1,itab_wsfc(iwsfc)%nwatm

            ja = 1  ! Temporary setting (while question about volti is being answered)

            iw = itab_wsfc(iwsfc)%iwatm(ja) ! local index
            kw = itab_wsfc(iwsfc)%kwatm(ja)


!            utstar =  uthfv( j,3 ) * fmoit

            ! wind erosion occurs only when U* > U*t

!            if ( ustar .le. utstar ) return

!            utem = ustar * ( ustar * ustar - utstar * utstar )
!           qam  = cfac * kvh(iland) * utem * sd_ep(iland) * dmsk( ladut(iland) ) * land%area(iland) * volti(kw,iw)

            ! Dust removal by surrounding vegetation <-???
            ! Adjust dust emissions for transport factors

            vels10 = vels * log(10.          / sfcg%rough(iwsfc))
     &                    / log(gdz_abov(kw) / sfcg%rough(iwsfc))

!BOB  ustar was accessed above as sfcg%ustar(iwsfc); should that one be used, or
!     should the following value be used?

            ! Calculate friction velocity (U*) at the surface applicable to dust emission
            ustar = vonk * vels10 / log( 10.0 / z0 )

!            ustar = land%ustar(iland) * log(gdz_abov(kw) / land%rough(iland))
!     &                              / log(gdz_abov(kw) / z0)
!!     &                              / log(gdz_abov(kw) / z0s(j))

            if     ( vels10 .le. 2.0 ) then  ! surface wind [0-2] range
               indx = 3
            elseif ( vels10 .lt. 7.0 ) then  ! surface wind (2-7) range
               indx = 2
            else                             ! surface wind [7-^) range
               indx = 1
            endif

            ku  = 0.08 * ustar
            tfa = ku / ( ku + vd( uland(lc),indx ) ) ! ratio dqup/dqrd

            tfb = 1.0 - fw( uland(lc) )

            ! Vert-to-Horiz dust flux ratio : Kang et al. [JGR, 2011] : Eq. (12)
            kvh   = flxfac1(j) * ( 0.24 + flxfac2(j) * ustar )

            ftot  = fmoit * fruf

            vflux = kvh * dust_hflux(j, ustar, dens, ftot)                           ! [g/m^2/s]

!            if (vflux > 1.e-20) then
!            if (ustar > 0.2) then
!               write(*,'(2I6,10g13.5)') iland, lc, fruf, fmoit, ustar, z0, vflux
!            endif

!            write(*,'(2I6,10g13.5)') iland, uland(lc), vegfree, lambda, land%rough(iland), z0,
!     &                land%ustar(iland), ustar, vonk * vels10 / log( 10.0 / z0 ), fruf, fmoit, vflux

            dust_em(iland) = vflux * tfa * tfb * ferod * sfcg%area(iwsfc) * volti(kw,iw) ! [g/m^3/s]

!                     if ( elus( c,r,m ) .gt. 0.0 ) then
!                        utstar( c,r,m ) = uthfv( j,3 ) * fmoit( c,r )
!                        if ( ustar( c,r ) .gt. utstar( c,r,m ) ) then
!                           ! wind erosion occurs only when U* > U*t
!                           utem( c,r,m ) = ustar( c,r )
!!                          utem( c,r,m ) = ( ustar( c,r ) + utstar( c,r,m ) )
!     &                                   * ( ustar( c,r ) * ustar( c,r )
!     &                                   -   utstar( c,r,m ) * utstar( c,r,m ) )
!                           qam( c,r,m ) = cfac * kvh( c,r ) * utem( c,r,m ) * sd_ep( c,r )
!                                                                          ! [g/m**2/s]
!                           qam( c,r,m ) = qam( c,r,m ) * rlay1hgt
!     &                                  * ( elus( c,r,m ) * 0.01 )        ! [g/m**3/s]
!                        end if
!                     end if
!
!                     if ( elus( c,r,m ) .eq. 0.0 .and. qam( c,r,m ) .ne. 0.0 ) then
!                        xmsg = 'ERROR: erodible land use = 0, but emissions .ne. 0'
!                        call m3exit( pname, jdate, jtime, xmsg, xstat1 )
!                     end if
!
!                     dust_em( c,r ) = dust_em( c,r ) + qam( c,r,m )
!
!                  end do   ! m landuse

! Dust removal by surrounding vegetation <-???
! Adjust dust emissions for transport factors

! J/K mass emis rate [g/m**3/s]  (edust( 1 ) not used)
            edust( 2 ) = fracmj * dust_em( iland )
            edust( 3 ) = fracmk * dust_em( iland )

!           do v = 1, ndust_spc
!              dustoutm( v,1,iland ) = 0.0
!           end do

            do n = 2, n_mode
               do v = 1, ndust_spc
                  dustoutm( v,n,iland ) = edust( n ) * dust_spc( v )%spcfac( n )
               end do
            end do

! J/K 3rd moment emis rate [1/s] (needed for number and surface)
            m3j = edust( 2 ) * f6dpi / ( gpkg * dustmode_dens( 2 ) )
            m3k = edust( 3 ) * f6dpi / ( gpkg * dustmode_dens( 3 ) )

! Mode-specific emission rates of particle number [1/m**3/s]
!           dustoutn( 1,iland ) = 0.0
            dustoutn( 2,iland ) = m3j * factnumj
            dustoutn( 3,iland ) = m3k * factnumk

! Mode-specific dry surface area emission rates [m2/m**3/s].
! 2nd moment multiplied by PI to obtain the surface area emissions rate.
!           dustouts( 1,iland ) = 0.0
            dustouts( 2,iland ) = m3j * factsrfj
            dustouts( 3,iland ) = m3k * factsrfk

!#ifdef verbose
!            if ( m3j .ne. 0.0 ) dusthit = dusthit + 1
!#endif
!
!            if ( dustem_diag ) then
!               diagv(  1 ) = qam( c,r,1 )            ! g/m**3/s
!               diagv(  2 ) = qam( c,r,2 )            ! g/m**3/s
!               diagv(  3 ) = qam( c,r,3 )            ! g/m**3/s
!               diagv(  4 ) = qam( c,r,4 )            ! g/m**3/s
!               diagv(  5 ) = dust_em( c,r )          ! g/m**3/s
!               diagv(  6 ) = elus( c,r,1 )           ! shrub
!               diagv(  7 ) = elus( c,r,2 )           ! shrubgrass
!               diagv(  8 ) = elus( c,r,3 )           ! barren
!               diagv(  9 ) = elus( c,r,4 )           ! cropland
!               diagv( 10 ) = ladut( c,r,1 ) + ladut( c,r,2 )
!     &                     + ladut( c,r,3 ) + agland( c,r )
!               diagv( 11 ) = tfa ( c,r )
!               diagv( 12 ) = tfb ( c,r )
!               diagv( 13 ) = fmoit( c,r )
!               diagv( 14 ) = sd_ep( c,r )
!               diagv( 15 ) = kvh  ( c,r )
!               diagv( 16 ) = wmax ( c,r )
!               diagv( 17 ) = uland( c,r,1 )
!               diagv( 18 ) = uland( c,r,2 )
!               diagv( 19 ) = uland( c,r,3 )
!               diagv( 20 ) = uland( c,r,4 )
!               diagv( 21 ) = ustar( c,r )
!
!               do v = 1, 4
!                  if ( utstar( c,r,v ) .ne. 999.9 ) then
!                     diagv( 21+v ) = utstar( c,r,v )
!                  else
!                     diagv( 21+v ) = 0.0
!                  end if
!               end do   ! diagv(25)
!
!               do v = 1, ndust_spc-1
!                  diagv( 25+v ) = dustoutm( v,2,c,r )
!               end do   ! diagv(25+18) = diagv(43)
!
!               do v = 1, 3
!                  diagv( 43+v ) = dustoutm( v,3,c,r )
!               end do   ! diagv(43+3) = diagv(46)
!               diagv( 47 ) = dustoutm( ndust_spc,3,c,r ) ! 'ASOIL'
!
!               diagv( 48 ) = dustoutn( 2,c,r )
!               diagv( 49 ) = dustoutn( 3,c,r )
!               diagv( 50 ) = dustouts( 2,c,r )
!               diagv( 51 ) = dustouts( 3,c,r )
!
!!              if ( stdevj( c,r ) .ne. 0.0 ) then
!!                 diagv( 52 ) = stdevj( c,r )
!!                 diagv( 53 ) = stdevk( c,r )
!!              end if
!
!               do v = 1, ndust_diag
!                  dustbf( v,c,r ) = dustbf( v,c,r ) + diagv( v )
!     &                                              * float( time2sec ( tstep( 2 ) ) )
!#ifdef verbose
!                  sdiagv( v ) = sdiagv( v ) + diagv( v )
!     &                                      * float( time2sec ( tstep( 2 ) ) )
!#endif
!               end do
!            end if   ! dustem_diag

!#ifdef verbose
!         write( logdev,'( /5x, a, 2i8 / )' ) 'dust hit count, out of total:',
!     &                                                dusthit, (c-1)*(r-1)
!#endif
!
!         if ( dustem_diag ) then
!C If last call this hour, write out the windblown dust emissions dignostics.
!C Then reset the sea-salt emissions array and local write counter.
!
!            wstep = wstep + time2sec( tstep( 2 ) )
!
!            if ( wstep .ge. time2sec( tstep( 1 ) ) ) then
!               if ( .not. currstep( jdate, jtime, sdate, stime, tstep( 1 ),
!     &                              mdate, mtime ) ) then
!                  xmsg = 'Cannot get step date and time'
!                  call m3exit( pname, jdate, jtime, xmsg, xstat3 )
!               end if
!               call nextime( mdate, mtime, tstep( 1 ) )
!               do v = 1, ndust_diag
!                  do r = 1, my_nrows
!                     do c = 1, my_ncols
!                        wrbuf( c,r ) = dustbf( v,c,r ) / float( wstep )
!                     end do
!                  end do
!                  if ( .not. WRITE3( ctm_dust_emis_1, diagnm( v ),
!     &                       mdate, mtime, wrbuf ) ) then
!                     xmsg = 'Could not write ' // trim( diagnm( v ) )
!     &                    // ' to CTM_DUST_EMIS_1'
!                     call m3exit( pname, mdate, mtime, xmsg, xstat1 )
!                  end if
!               end do
!               write( logdev, '( /5x, 2( a, 1x ), i8, ":", i6.6 )' )
!     &               'Timestep written to CTM_DUST_EMIS_1',
!     &               'for date and time', mdate, mtime
!#ifdef verbose
!               write( logdev,2015 ) jdate, jtime
!               do v = 1, ndust_diag
!                  if ( diagnm( v )(1:4) .ne. 'ANUM' ) then
!                     write( logdev,2019 ) v, diagnm( v ), sdiagv( v ) / float( wstep )
!                  else
!                     write( logdev,2023 ) v, diagnm( v ), sdiagv( v ) / float( wstep )
!                  end if
!               end do
!               sdiagv = 0.0   ! array assignment
!#endif
!               wstep = 0
!               dustbf = 0.0   ! array assignment
!            end if   !  time to write
!         end if   ! dustem_diag
!
!2009     Format( '*** Erodible landuse incorrect ', 1pe13.5, 1x, 'at: ', 3i4 )
!2015     format( /5x, 'Total grid time-avg sum of dust emis variables at:',
!     &            1x, i8, ":", I6.6  )
!2019     format( i10, 1x, a6, f20.5 )
!2023     format( i10, 1x, a6, e20.3 )


         end subroutine get_dust_emis1


         real function dust_hflux(j, ustr, dens, ftot)

         use const_data, only: grav

         implicit none

         integer, intent(in) :: j
         real,    intent(in) :: ustr
         real,    intent(in) :: dens
         real,    intent(in) :: ftot

         real, parameter :: amen = 1.0          ! Marticorena and Bergametti [JGR,1997]
         real, parameter :: cfac = 1000.0 * amen / grav

         real, parameter :: A = 260.60061       ! 0.0123 * 2650.0 * 9.81 / 1.227
         real, parameter :: B = 1.6540342e-06   ! 0.0123 * 0.000165 / 1.227
         real, parameter :: utfac(ndp) = sqrt( A * dp + B / dp )

         real    :: utstar                      ! threshold U* [m/s]
         real    :: hflux
         integer :: n

         hflux = 0.0

         do n = 1, ndp
            utstar = utfac(n) * ftot   ! Shao and Lu [JGR,2000]

            ! wind erosion occurs only if U* > U*t
            if (ustr > utstar) then

               !---Horiz. Flux from White (1979)
               hflux = hflux + soiltxt(n,j) * ( ustr + utstar )
     &                       * ( ustr * ustr - utstar * utstar )

            endif
         enddo

         dust_hflux = cfac * hflux * dens * sd_ep(j)

         end function dust_hflux

      end module dust_emis
