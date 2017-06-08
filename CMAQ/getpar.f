
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
C $Header: /project/yoj/arc/CCTM/src/aero/aero5/getpar.f,v 1.7 2012/01/19 13:13:27 yoj Exp $

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Subroutine getpar( limit_sg  )

C  Calculates the 3rd moments (M3), masses, aerosol densities, and
C  geometric mean diameters (Dg) of all 3 modes, and the natural logs of
C  geometric standard deviations (Sg) of the Aitken and accumulation modes.

C  The logical variable, WET_MOMENTS_FLAG, dictates whether the
C  calculations in GETPAR are to assume that the aerosol is "wet" or
C  "dry."  In the present context, a "wet" aerosol consists of all
C  chemical components of the aerosol.  A "dry" aerosol excludes
C  particle-bound water and also excludes semivol secondary organic aerosol.

C  NOTE! 2nd moment concentrations (M2) are passed into GETPAR in the
C  CBLK array and are modified within GETPAR only in the event that
C  the Sg value of a given mode has gone outside of the acceptable
C  range (1.05 to 2.50).  The GETPAR calculations implicitly assume
C  that the input value of M2 is consistent with the input value of
C  WET_MOMENTS_FLAG.  If, for example, the input M2 value was calculated
C  for a "dry" aerosol and the WET_MOMENTS_FLAG is .TRUE., GETPAR would
C  incorrectly adjust the M2 concentrations!
C
C  
C  Outputs: 
C    moment3_conc
C    moment2_conc (adjusted if standard dev. hits limit)
C    aeromode_dens
C    aeromode_lnsg
C    aeromode_diam
C    aeromode_mass
C
C SH  03/10/11 Renamed met_data to aeromet_data
C HP and BM 4/2016: Updated use of wet_moments_flag which is now
C    available through AERO_DATA consistent with the moments it refers to

C-----------------------------------------------------------------------

      Use aero_data
      Use aeromet_data   ! Includes CONST.EXT

      Implicit None

C Arguments:
      Logical, Intent( In ) :: limit_sg  ! fix coarse and accum Sg's to the input value?

C Output variables:
C  updates arrays in aero_data module
C  moment3_conc   3rd moment concentration [ ug /m**3 ]
C  aeromode_mass  mass concentration: [ ug / m**3 ]
C  aeromode_dens  avg particle density [ kg / m**3 ]
C  aeromode_diam  geometric mean diameter [ m ]
C  aeromode_lnsg  log of geometric standard deviation

C Local Variables:
!     Real :: xfsum       ! (ln(M0)+2ln(M3))/3; used in Sg calcs
!     Real :: lxfm2       ! ln(M2); used in Sg calcs
      Real :: l2sg        ! square of ln(Sg); used in diameter calcs
!     Real :: es36        ! exp(4.5*l2sg); used in diameter calcs
      real :: es36_one3
      real :: xm0_one3
      real :: xm3_one3

      Real, Parameter :: one3  = 1.0 / 3.0
      Real, Parameter :: dgmin = 1.0E-09   ! minimum particle diameter [ m ]
      Real, Parameter :: densmin = 1.0E03  ! minimum particle density [ kg/m**3 ]
      real, parameter :: dgmax(n_mode) = min(def_diam(1:n_mode) * 100., 5.e-6)
      Real, parameter :: minl2sg = Log( min_sigma_g ) ** 2
      Real, parameter :: maxl2sg = Log( max_sigma_g ) ** 2
      Real, parameter :: minel2sg = exp( minl2sg )
      Real, parameter :: maxel2sg = exp( maxl2sg )

      Real( 8 ) :: sumM3  ( n_mode )
      Real( 8 ) :: sumMass( n_mode )
      Integer   :: n, spc   ! loop counters

      real :: exfsum, el2sg

C-----------------------------------------------------------------------

C *** Calculate aerosol 3rd moment concentrations [ m**3 / m**3 ]

      Do n = 1, n_mode
         sumM3  (n) = 0.0_8
         sumMass(n) = 0.0_8
      End Do

      Do spc = 1, n_aerospc
         If ( aerospc( spc )%tracer) cycle
         If ( aerospc( spc )%no_M2Wet .AND. .Not. wet_moments_flag ) Cycle
         do n = 1, n_mode
            if ( aero_missing(spc,n) ) cycle
            sumM3  (n) = sumM3  (n) + aerospc_conc(spc,n) * aerospc_m3fac(spc)
            sumMass(n) = sumMass(n) + aerospc_conc(spc,n)
         End Do
      End Do

C *** Calculate modal average particle densities [ kg/m**3 ]
      do n = 1, n_mode
         moment3_conc (n) = Max( Real(sumM3(n)), aeromode( n )%min_m3conc )
         aeromode_mass(n) = sumMass(n)
         aeromode_dens(n) = max(f6piove9 * aeromode_mass(n) / moment3_conc(n), densmin)
      enddo

C *** Calculate geometric standard deviations as follows:
c        ln^2(Sg) = 1/3*ln(M0) + 2/3*ln(M3) - ln(M2)
c     NOTES:
c      1. Equation 10-5a of [Binkowski:1999] and Equation 5a of
c         Binkowski&Roselle(2003) contain typographical errors.
c      2. If the square of the logarithm of the geometric standard
c         deviation is out of an acceptable range, reset this value and
c         adjust the second moments to be consistent with this value.
c         In this manner, M2 is artificially increased when Sg exceeds
c         the maximum limit.  M2 is artificially decreased when Sg falls
c         below the minimum limit.

C *** Aitken Mode:

      if ( limit_sg ) then

         Do n = 1, n_mode
            l2sg = aeromode_lnsg(n) ** 2
            
            xm0_one3 = moment0_conc(n) ** one3
            xm3_one3 = moment3_conc(n) ** one3

            moment2_conc(n) = xm0_one3 * xm3_one3 ** 2 * Exp( -l2sg )

            ES36_one3 = Exp( 1.5 * l2sg )

            aeromode_diam( n ) = max( dgmin, ( xm3_one3 / ( xm0_one3 * es36_one3 ) ) )

            aeromode_diam( n ) = min( dgmax(n), aeromode_diam(n) )

         End Do

      else

         Do n = 1, n_mode

            xm0_one3 = moment0_conc(n) ** one3
            xm3_one3 = moment3_conc(n) ** one3

            exfsum = xm0_one3 * xm3_one3 ** 2

            el2sg = exfsum / moment2_conc( n ) 

            el2sg = Max( el2sg, minel2sg )
            el2sg = Min( el2sg, maxel2sg )

            moment2_conc( n )  = exfsum / el2sg

            l2sg = log( el2sg )

            aeromode_lnsg( n ) = sqrt( l2sg )

!           ES36_one3 = Exp( 1.5 * l2sg )
            ES36_one3 = el2sg * sqrt(el2sg)

            aeromode_diam( n ) = max( dgmin, ( xm3_one3 / ( xm0_one3 * es36_one3 ) ) )

            aeromode_diam( n ) = min( dgmax(n), aeromode_diam(n) )

         End Do

      endif

      End Subroutine getpar
