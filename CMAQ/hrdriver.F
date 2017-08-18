
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

      subroutine chem( mrl )

      USE RXNS_DATA
      USE RXNS_FUNCTION
      USE UTILIO_DEFN 
      use hrdata
      use mem_ijtabs, only: jtab_w, jtw_prog
      use misc_coms,  only: io6

      implicit none

      integer, intent(in) :: mrl
      LOGICAL, SAVE       :: LFIRST = .TRUE.
      integer             :: j, iw
      CHARACTER( 132 )    :: MSG

      INTEGER :: MXEBI
      REAL    :: SMEBI
      INTEGER :: NCALL
!
      integer :: maxebi
      real    :: sumebi
      integer :: ncall_ebi
      integer :: avgebi

      IF ( NUMB_MECH_SPC .EQ. 0 ) RETURN

!  On first call, call routines to set-up for EBI solver

      IF ( LFIRST ) THEN

         LOGDEV = INIT3( )

         IF ( MECHNAME .NE. 'CB05TUCL_AE6_AQ' ) THEN
            MSG = 'This version of the EBI solver can only be used with'
     &           // ' the CB05TUCL_AE6_AQ chemical mechanism'
            CALL M3EXIT( 'HRDRIVER', 0, 0, MSG, XSTAT1 )
         END IF

         IF( INDEX( MECHNAME, 'AE' ) .NE. 0 ) THEN
            L_AE_VRSN = .TRUE.
         ELSE
            L_AE_VRSN = .FALSE.
         END IF

         IF( INDEX( MECHNAME, 'AQ' ) .NE. 0 ) THEN
            L_AQ_VRSN = .TRUE.
         ELSE
            L_AQ_VRSN = .FALSE.
         END IF

         LFIRST = .FALSE.

      END IF      ! First time

      maxebi = 0
      sumebi = 0.
      ncall_ebi = 0

      !$omp parallel do private(iw,mxebi,smebi,ncall) reduction(+:sumebi,ncall_ebi) reduction(max:maxebi)
      do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

         call chem1(iw,mxebi,smebi,ncall)
         
         maxebi    = max(maxebi,mxebi)
         sumebi    = sumebi + smebi
         ncall_ebi = ncall_ebi + ncall

      enddo
      !$omp end parallel do

#ifdef EBI_STATS
      avgebi = nint(sumebi / real(ncall_ebi))
      write(io6,'(1x,A,2(I0,1x))') " EBI mn, mx, rates: ", avgebi, maxebi
#endif

      end subroutine chem


      SUBROUTINE CHEM1( iw, mxebi, smebi, ncall )

C**********************************************************************
C
C  FUNCTION: Driver subroutine for Euler Backward Iterative solver
C
C  PRECONDITIONS: For the CB05TUCL family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  HRINIT, PHOT, HRCALCKS, HRSOLVER
C
C  REVISION HISTORY: Created by EBI solver program, May 26, 2017
C                       Based on the algorithm in "Test of Two Numerical
C                       Schemes for Use in Atmospheric Transport-Chemistry
C                       Models", O. Hertel, R. Berkowicz, J. Christensen,
C                       and O. Hov, Atm Env., Vol. 27A, No. 16, 1993.
C                       Original MEBI code developed by Ho-Chun Huang,
C                       SUNY, Albany -- "On the performance of numerical
C                       solvers for a chemistry submodel in three-dimensional
C                       air quality models 1. Box model simulations",
C                       H. Huang and J.S. Chang, JGR, Vol 106, No. D17, 2001.
C                       This version replaces Huang and Chang use of numerical
C                       solutions with analytical solutions derived in
C                       Hertel et al.
C   21 Jun 10 J.Young: convert for Namelist redesign
C   11 May 11 D.Wong: incorporated twoway model implementation
C   27 Sep 11 B.Hutzell: revised method for defining CALL_DEG variable
C   18 Jul 14 B.Hutzell: revised: 1) to use the RXNS_DATA and RXNS_FUNCTION
C                        modules instead of include files, 2) to enable
C                        reactions between species types, 3) to calculate
C                        heterogeneous rate constants via AEROSOL_CHEMISTRY
C                        module, 4) to replace call to HRCALCKS with
C                        subroutine in RXNS_FUNCTION module and 5) to change
C                        how degrade routines are used, if present
C   02 Dec 14 B.Hutzell: 1) added terrestrial data to conduct surface
C                        dependent reactions and 2) modified the call CALC_RCONST
C                        routine
C**********************************************************************

!     USE HGRD_DEFN             ! horizontal domain specifications
!     USE VGRD_DEFN             ! vertical layer specifications
      USE CGRID_SPCS            ! CGRID mechanism species
      USE UTILIO_DEFN           ! IOAPI parameters and functions declarations
      USE RXNS_DATA
      USE AEROSOL_CHEMISTRY
      USE RXNS_FUNCTION
      USE HRDATA
      USE AERO_DATA,  ONLY: N_MODE
      use cgrid_defn, only: cgrid
      use misc_coms,  only: dtlm
      use mem_grid,   only: lpw, mza
      use mem_ijtabs, only: itab_w
      use mem_basic,  only: tair, sh_v, press, rho
      use mem_sea,    only: sea
      use mem_turb,   only: frac_land
      use cgrid_defn, only: cgrid

      IMPLICIT NONE

C..Includes:
!     INCLUDE 'FILES_CTM.EXT'  ! CMAQ files
      INCLUDE 'CONST.EXT'      ! CMAQ constants

!     INCLUDE SUBST_EMISPRM   ! Emissions processing control parameters
!     INCLUDE SUBST_PACTL_ID  ! Process analysis control parameters

C..Arguments:
!     REAL, POINTER :: CGRID( :,:,:,: )  ! Species concentrations
!     INTEGER JDATE           ! Current date (YYYYDDD)
!     INTEGER JTIME           ! Current time (HHMMSS)
!     INTEGER TSTEP( 3 )      ! Time step vector (HHMMSS)
      integer, intent(in) :: iw

C..Parameters:
      REAL( 8 ), PARAMETER :: DCONMIN = 1.0D-30               ! minimum species concentration allowed
      REAL,      PARAMETER :: CONCMIN = 1.0E-30               ! minimum species concentration allowed
      REAL,      PARAMETER :: MAOMV   = 1.0E6 * MWAIR / MWWAT ! Mol Wt of air over Mol Wt of water times 1.0E6

C..External Functions:


C..Saved Local Variables:

!     CHARACTER( 16 ), SAVE :: PNAME = 'HRDRIVER'     ! Program name

!      INTEGER, SAVE :: ISTFL            ! Unit no. of iteration stat output file
!      LOGICAL, SAVE :: LFIRST = .TRUE.  ! Flag for first call to this subroutine

!     REAL( 8 ), SAVE :: PA2ATM   ! Pascal to atm conv fac
      real( 8 ), parameter :: pa2atm = 1.d0 / real(STDATMPA,8)


C..Scratch Local Variables:
!      CHARACTER( 132 ) :: MSG           ! Message text
!      CHARACTER(  16 ) :: VNAME         ! Name of I/O API data variable

      integer :: s, k, jws, iws, mrlw
      real    :: area, frac_ice

!     INTEGER DELT_SEC        ! EBI max time step in seconds
!     INTEGER ESP             ! Loop index for emissions species
!     INTEGER ITMSTEP         ! Chemistry integration interval (sec)
!     INTEGER LEV             ! Layer index
!     INTEGER MIDDATE         ! Date at time step midpoint
!     INTEGER MIDTIME         ! Time at time step midpoint
      INTEGER NPH             ! Index for number of phot. rxns in PHOT
      INTEGER SPC             ! Species loop index

      INTEGER, intent(inout) :: MXEBI           ! Max no. of EBI iterations
      REAL,    intent(inout) :: SMEBI           ! Sum of EBI iterations
      INTEGER, intent(inout) :: NCALL
!     real :: rjmax, hemax, rkmax

      LOGICAL LSUNLIGHT       ! Flag for sunlight

      REAL    ::      INV_DENS(mza)     ! reciprocal of air mass density, m3/Kg
      REAL(8) ::      TEMP    (mza)
      REAL(8) ::      ATMPRES (mza)
      REAL(8) ::      H2O     (mza)

      real(8) :: het            ( mza, nhetero )
      real    :: chem_m3dry_init( n_mode,  mza )
      integer :: kmax
      logical :: land

      INTERFACE
        SUBROUTINE HRSOLVER( JDATE, JTIME, C, R, L )
           IMPLICIT NONE
           INTEGER, INTENT( IN ) :: JDATE    ! Current date (YYYYDDD)
           INTEGER, INTENT( IN ) :: JTIME    ! Current time (HHMMSS)
           INTEGER, INTENT( IN ) :: C, R, L  ! Cell col, row, lev
        END SUBROUTINE HRSOLVER
      END INTERFACE

#ifdef hrstats
!     MNEBI = 1000
      hemax = 0.
      rjmax = 0.
      rkmax = 0.
#endif
      MXEBI = 0
      SMEBI = 0.0
      NCALL = 0

      mrlw = itab_w(iw)%mrlw

      EBI_TMSTEP = dtlm(mrlw) / 60.d0

         do k = mza-1, lpw(iw), -1
            if (press(k,iw) > 100.d2) exit
         enddo

         kmax = k + 1

         ! Don't do tropospheric chem above 100 mb

         do spc = 1, numb_mech_spc
            s = cgrid_index( spc )
            do k = kmax+1, mza
               cgrid( k, iw, s ) = max( cgrid( k, iw, s ), concmin)
            enddo
         enddo

         ! Set physical quantities

         do k = lpw(iw), kmax
            TEMP    ( k ) = tair(k,iw)
            ATMPRES ( k ) = PA2ATM * press(k,iw)
            H2O     ( k ) = MAX( MAOMV * sh_v(k,iw), 0.0 )
            inv_dens( k ) = 1.0 / rho(k,iw)
         enddo

         ! Is column mostly land or open sea (map seaice to land)

         frac_ice = 0.0

         if (itab_w(iw)%nsea > 0) then
            area = 0.0
            do jws = 1, itab_w(iw)%nsea
               iws = itab_w(iw)%isea(jws)
               frac_ice = frac_ice + sea%area(iws) * sea%seaicec(iws)
               area     = area     + sea%area(iws)
            enddo
            frac_ice = frac_ice / area
         endif

         if ((frac_land(iw) + (1.0-frac_land(iw)) * frac_ice) > 0.5) then
            land = .true.
         else
            land = .false.
         end if

C.. Get Heteorogeneous rates and Update Aerosol Distribution Properties
         call hetchem_rates( iw, kmax, het, chem_m3dry_init )

         if ( any( rj(lpw(iw),:,iw) > 1.d-7 ) ) then
            lsunlight = .true.
         else
            lsunlight = .false.
         endif

c..Get rate constants

         CALL CALC_RCONST( TEMP, ATMPRES, H2O, RJ(:,:,iw), HET, LSUNLIGHT, LAND, RKI_SAV, lpw(iw), kmax )

#ifdef hrstats
         hemax = max( maxval(het    (lpw(iw):kmax,:   )), hemax )
         rjmax = max( maxval(rj     (lpw(iw):kmax,:,iw)), rjmax )
         rkmax = max( maxval(rki_sav(lpw(iw):kmax,:   )), rkmax )
#endif

c..Call EBI solver

         do k = lpw(iw), kmax

            ! Load ICs
            do spc = 1, numb_mech_spc
               s = cgrid_index( spc )

               if ( convert_conc( spc ) )then
                  yc( spc ) = max( concmin, inv_dens( k ) * forward_conv( spc ) * cgrid( k,iw,s ) )
               else
                  yc( spc ) = max( concmin, cgrid( k,iw,s ) )
               endif
            enddo

            ! Load reaction rates
            do nph = 1, nrxns
               rki(nph) = rki_sav(k,nph)
            enddo

            CALL HRSOLVER( 0, 0, iw, itab_w(iw)%iwglobe, k )

            MXEBI = MAX( MXEBI, N_EBI_IT )
            SMEBI = SMEBI + REAL( N_EBI_IT )
            NCALL = NCALL + 1

c..Update concentration array
            DO SPC = 1, NUMB_MECH_SPC
               S = CGRID_INDEX( SPC )
               IF ( CONVERT_CONC( SPC ) )THEN
                  CGRID( k, iw, S ) = REVERSE_CONV( SPC ) * rho( k,iw ) * MAX( YC( SPC ), DCONMIN )
               ELSE
                  CGRID( k, iw, S ) = MAX( YC( SPC ), DCONMIN )
               END IF
            END DO

         enddo

         ! Update Aerosol Surface Area
         CALL HETCHEM_UPDATE_AERO( chem_m3dry_init, iw, lpw(iw), kmax )

#ifdef hrstats
!      AVGEBI = nint( SUMEBI / REAL( NCALL ) )
!      WRITE( ISTFL, 94040 ) JDATE, JTIME, MNEBI, AVGEBI, MXEBI
!      write(io6,'(1x,I0.2,A,2(I0,1x),3(g12.5,1x))') myrank, " EBI mn, mx, rates: ", avgebi, mxebi, hemax, rjmax, rkmax
#endif

C*********************** FORMAT STATEMENTS ****************************

92000 FORMAT( / 10X, 'Euler Backward Iterative Parameters -'
     &        / 10X, 'Chemistry Integration Time Interval (min):', F12.4,
     &        / 10X, 'EBI maximum time step (min):              ', F12.4 )

92020 FORMAT( //10X, 'Species convergence tolerances:' )

92040 FORMAT(   10X, A16, 2X, 1PE12.2 )

92060 FORMAT( / 10X, 'Emissions Processing in Chemistry ...'
     &        / 10X, 'Number of Emissions Layers:         ', I3
     &        / 10X, 'out of total Number of Model Layers:', I3 )


94020 FORMAT( 'DATE      TIME ', 'MNEBI AVEBI MXEBI' )

94040 FORMAT( I7, 1X, I6, 1X, 3( I5, 1X ) )

      END
