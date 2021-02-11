
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
      USE UTILIO_DEFN
      USE AEROSOL_CHEMISTRY
      use hrdata
      use mem_ijtabs, only: jtab_w, jtw_prog
      use aero_data,  only: aerospc_m3fac, aer_str, spec_map
      use cgrid_spcs, only: n_gc_spc, gc_spc
      use cgrid_defn, only: cgrid

#ifdef EBI_STATS
      use misc_coms, only: io6
#endif

      implicit none

      integer, intent(in) :: mrl
      LOGICAL, SAVE       :: LFIRST = .TRUE.
      integer             :: j, iw, n, v, s, l
      CHARACTER( 132 )    :: MSG

      INTEGER :: MXEBI
      REAL    :: SMEBI
      INTEGER :: NCALL

      integer :: maxebi, avgebi
      real    :: sumebi
      integer :: ncall_ebi

      IF ( NUMB_MECH_SPC .EQ. 0 ) RETURN

!  On first call, call routines to set-up for EBI solver

      IF ( LFIRST ) THEN

         LOGDEV = INIT3( )

         IF( MECHNAME .NE. 'CB6R3_AE6_AQ' ) THEN
             MSG = 'This version of the EBI solver can only be used with'
     &            // ' the CB6R3_AE6_AQ chemical mechanism'
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

         do n = 1, n_spec
            if (convert_conc(n)) then
               v = cgrid_index(n) - aer_str + 1
               s = spec_map(v)
               m3_conv(n) = reverse_conv(n) * aerospc_m3fac(s)
            else
               m3_conv(n) = 1.0
            endif
         enddo

         CALL HETCHEM_INIT()

         nrxn = 0
         do n = 1, n_gc_spc
            l = len_trim(gc_spc(n))
            if (l > 2) then
               if (gc_spc(n)(l-2:l) == 'RXN') nrxn = nrxn + 1
            endif
         enddo

         allocate(irxns(nrxn))

         s = 0
         do n = 1, n_gc_spc
            l = len_trim(gc_spc(n))
            if (l > 2) then
               if (gc_spc(n)(l-2:l) == 'RXN') then
                  s = s + 1
                  irxns(s) = n
               endif
            endif
         enddo

         LFIRST = .FALSE.

      END IF      ! First time

      maxebi = 0
      sumebi = 0.
      ncall_ebi = 0

!$omp parallel do private(iw,s,mxebi,smebi,ncall) schedule(guided)
!$omp&            reduction(+:sumebi,ncall_ebi) reduction(max:maxebi)
      do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

         do s = 1, nrxn
            cgrid(:,iw,irxns(s)) = 1.e-30
         enddo

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
C  PRECONDITIONS: For the CB6R3 family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  HRINIT, PHOT, HRCALCKS, HRSOLVER
C
C  REVISION HISTORY: Created by EBI solver program, Mar 10, 2017
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

      USE CGRID_SPCS            ! CGRID mechanism species
      USE UTILIO_DEFN           ! IOAPI parameters and functions declarations
      USE RXNS_DATA
      USE AEROSOL_CHEMISTRY
      USE RXNS_FUNCTION
      USE HRDATA
      USE CONST_DATA
      USE AERO_DATA,  ONLY: N_MODE, SRF_STR
      USE CGRID_DEFN, ONLY: CGRID

      use misc_coms,  only: dtlm
      use mem_grid,   only: lpw, mza
      use mem_ijtabs, only: itab_w
      use mem_sfcg,   only: sfcg, itab_wsfc
      use mem_basic,  only: tair, rr_v, press, rho
      use mem_sea,    only: sea, omsea
      use mem_turb,   only: frac_land

      IMPLICIT NONE

C..Arguments:
      integer, intent(in) :: iw

C..Parameters:
      REAL( 8 ), PARAMETER :: DCONMIN = 1.0D-30               ! minimum species concentration allowed
      REAL,      PARAMETER :: CONCMIN = 1.0E-30               ! minimum species concentration allowed
      REAL,      PARAMETER :: MAOMV   = 1.0E6 * MWAIR / MWWAT ! Mol Wt of air over Mol Wt of water times 1.0E6
      REAL,      PARAMETER :: PA2ATM  = 1.0 / STDATMPA
      REAL( 8 ), PARAMETER :: MAXPRED = 1.0D+03    ! Upper limit on predicted conc

C..Saved Local Variables:

!     CHARACTER( 16 ), SAVE :: PNAME = 'HRDRIVER'     ! Program name

!     INTEGER, SAVE :: ISTFL            ! Unit no. of iteration stat output file
!     LOGICAL, SAVE :: LFIRST = .TRUE.  ! Flag for first call to this subroutine



C..Scratch Local Variables:
!      CHARACTER( 132 ) :: MSG           ! Message text
!      CHARACTER(  16 ) :: VNAME         ! Name of I/O API data variable

      integer :: s, k, j, jsfc, iwsfc, jasfc, isea, mrlw
      real    :: area, frac_ice

!     INTEGER DELT_SEC        ! EBI max time step in seconds
!     INTEGER ESP             ! Loop index for emissions species
!     INTEGER ITMSTEP         ! Chemistry integration interval (sec)
!     INTEGER LEV             ! Layer index
!     INTEGER MIDDATE         ! Date at time step midpoint
!     INTEGER MIDTIME         ! Time at time step midpoint
!     INTEGER NPH             ! Index for number of phot. rxns in PHOT
      INTEGER SPC             ! Species loop index

      INTEGER, intent(inout) :: MXEBI           ! Max no. of EBI iterations
      REAL,    intent(inout) :: SMEBI           ! Sum of EBI iterations
      INTEGER, intent(inout) :: NCALL
!     real :: rjmax, hemax, rkmax

      LOGICAL LSUNLIGHT       ! Flag for sunlight

      REAL ::      DENS    (mza)
      REAL ::      INV_DENS(mza)     ! reciprocal of air mass density, m3/Kg
      REAL ::      TEMP    (mza)
      REAL ::      ATMPRES (mza)
      REAL ::      H2O     (mza)

      real(8) :: het            ( mza, nhetero )
      real    :: chem_m3_init( n_mode,  mza )
      integer :: kmax, m
      logical :: land
      real    :: moment3_conc(n_mode)

      real, parameter :: onethird = 1. / 3.

      real(8) :: rki_sav (mza, nrxns)
      real(8) :: rki_trns(nrxns, mza)

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
            ATMPRES ( k ) = PA2ATM * real( press(k,iw) )
            H2O     ( k ) = MAX( MAOMV * rr_v(k,iw), 0.0 )
            DENS    ( k ) = real( rho(k,iw) )
            inv_dens( k ) = 1.0 / dens(k)
         enddo

         ! Is column mostly land or open sea (map seaice to land)

         frac_ice = 0.0

         if (itab_w(iw)%jsea2 > 0) then
            area = 0.0
            do jsfc = itab_w(iw)%jsea1, itab_w(iw)%jsea2
               iwsfc = itab_w(iw)%iwsfc(jsfc)
               jasfc = itab_w(iw)%jasfc(jsfc)

               isea = iwsfc - omsea

               frac_ice = frac_ice + itab_wsfc(iwsfc)%arc(jasfc) * sea%seaicec(isea)
               area     = area     + itab_wsfc(iwsfc)%arc(jasfc)
            enddo
            frac_ice = frac_ice / area
         endif

         if ((frac_land(iw) + (1.0-frac_land(iw)) * frac_ice) > 0.5) then
            land = .true.
         else
            land = .false.
         end if

C.. Get Heteorogeneous rates and Update Aerosol Distribution Properties
         call hetchem_rates( iw, kmax, het, chem_m3_init )

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

         rki_trns = transpose( rki_sav )

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

            ! Special for methane
            yc(ech4) = max(yc(ech4) - ATM_CH4, dconmin)

            yc  = min(yc, maxpred)
            yc0 = yc

            ! Load reaction rates
            rki = rki_trns(:,k)

            CALL HRSOLVER( 0, 0, iw, itab_w(iw)%iwglobe, k )

            MXEBI = MAX( MXEBI, N_EBI_IT )
            SMEBI = SMEBI + REAL( N_EBI_IT )
            NCALL = NCALL + 1

            ! Special for methane

            yc(ech4) = cgrid(k,iw,cgrid_index(ech4)) + yc(ech4) - yc0(ech4)

            ! Update concentration array

            DO SPC = 1, NUMB_MECH_SPC
               S = CGRID_INDEX( SPC )
               IF ( CONVERT_CONC( SPC ) )THEN
                  CGRID( k, iw, S ) = REVERSE_CONV( SPC ) * rho( k,iw ) * MAX( YC( SPC ), DCONMIN )
               ELSE
                  CGRID( k, iw, S ) = MAX( YC( SPC ), DCONMIN )
               END IF
            END DO

            ! Update Aerosol Surface Area. These need to follow the species
            ! in RXNS_DATA_MODULE

            moment3_conc(1) = chem_m3_init(1,k)
     &                      + real(yc(107) - yc0(107)) * m3_conv(107) * dens(k)

            moment3_conc(2) = chem_m3_init(2,k)
     &                      + real(yc(108) - yc0(108)) * m3_conv(108)

            moment3_conc(3) = chem_m3_init(3,k)
     &                      + real(yc(109) - yc0(109)) * m3_conv(109) * dens(k)

            do spc = 110, 128
               moment3_conc(2) = moment3_conc(2)
     &                         + real(yc(spc) - yc0(spc)) * m3_conv(spc)
            enddo

            Do m = 1, n_mode
               s = srf_str + m - 1
               CGRID( K,IW,s ) = CGRID( K,IW,s )
     &              * ( (MOMENT3_CONC( m ) / chem_m3_init(m,k))**2 ) ** ONETHIRD
            End Do

         End Do

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
