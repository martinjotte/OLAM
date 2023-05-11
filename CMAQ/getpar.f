
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

      module getpar_mod

      contains

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Subroutine getpar( fixed_sg, moment0_conc, moment2_conc, moment3_conc,
     &                   aeromode_lnsg, aeromode_diam, i1, i2 )

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
C  Outputs:
C    moment3_conc  third moment, porportional to volume [ m3/m3 ]
C    moment2_conc  second moment, prop. to surface area [ m2/m3 ]
C                     (adjusted if standard dev. hits limit)
C    aeromode_dens [ kg/m3 ]
C    aeromode_lnsg log of geometric standard deviation
C    aeromode_diam geometric mean diameter [ m ]
C    aeromode_mass mass concentration: [ ug / m**3 ]
C
C SH  03/10/11 Renamed met_data to aeromet_data
C HP and BM 4/2016: Updated use of wet_moments_flag which is now
C    available through AERO_DATA consistent with the moments it refers to

C-----------------------------------------------------------------------

      Use aero_data,  only: min_diam_g, min_sigma_g, max_sigma_g, n_mode,
     &                      min_l2sg, max_l2sg

      Implicit None

C Arguments:

      Logical, Intent( In ) :: fixed_sg  ! If TRUE, then the second moment is modified
                                         ! during each call in order to preserve the
                                         ! standard deviaiton at the current value.
                                         !
                                         ! If FALSE, then the standard deviation is
                                         ! recalculated to be consistent with the current
                                         ! combination of the 0th, 2nd and 3rd moments.
                                         ! During this calculation, standard deviaiton is
                                         ! limited by parameters in AERO_DATA (min_sigma_g
                                         ! and max_sigma_g)

      Real,    Intent( In    ) :: moment0_conc ( n_mode )
      Real,    Intent( In    ) :: moment3_conc ( n_mode )
      Real,    Intent( InOut ) :: moment2_conc ( n_mode )
      Real,    Intent( InOut ) :: aeromode_lnsg( n_mode ) ! log of geometric standard deviation
      Real,    Intent( Out   ) :: aeromode_diam( n_mode ) ! geometric mean diameter [ m ]
      Integer, Intent( In    ) :: i1, i2

C Local Variables:

      Real :: l2sg        ! square of ln(Sg); used in diameter calcs
      real :: es36_one3
      real :: xm0_one3
      real :: xm3_one3

      Real, Parameter :: one3     = 1.0 / 3.0
      Real, parameter :: minel2sg = exp( min_l2sg )
      Real, parameter :: maxel2sg = exp( max_l2sg )

      Integer :: n, ns, ne
      real    :: exfsum, el2sg

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

      ns = max(i1,1)
      ne = min(i2,n_mode)

      if ( fixed_sg ) then

         Do n = ns, ne
            l2sg = aeromode_lnsg(n) ** 2

            xm0_one3 = moment0_conc(n) ** one3
            xm3_one3 = moment3_conc(n) ** one3

            moment2_conc(n) = xm0_one3 * xm3_one3 ** 2 * Exp( -l2sg )

            ES36_one3 = Exp( 1.5 * l2sg )

            aeromode_diam( n ) = max( min_diam_g(n), ( xm3_one3 / (xm0_one3 * es36_one3) ) )
         End Do

      else

         Do n = ns, ne
            xm0_one3 = moment0_conc(n) ** one3
            xm3_one3 = moment3_conc(n) ** one3

            exfsum = xm0_one3 * xm3_one3 ** 2

            el2sg = max( minel2sg, exfsum / moment2_conc( n ) )
            el2sg = Min( maxel2sg, el2sg )

            moment2_conc( n )  = exfsum / el2sg

            l2sg = log( el2sg )

            aeromode_lnsg( n ) = sqrt( l2sg )

            ES36_one3 = el2sg * sqrt(el2sg)

            aeromode_diam( n ) = max( min_diam_g(n), ( xm3_one3 / (xm0_one3 * es36_one3) ) )
         End Do

      endif

      End Subroutine getpar

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Subroutine getdens(aeromode_mass, aeromode_dens,
     &                   aerospc_conc, moment3_conc, i1, i2)

!     Calculates total aerosol mass and average particle density
!     for each mode specified.

      use const_data, only: f6piove9
      use aero_data,  only: n_mode, idry_str, idry_end, iwet_end, n_aerospc

      implicit None

      real,    intent(in ) :: moment3_conc (n_mode)
      real,    intent(in ) :: aerospc_conc (n_aerospc, n_mode)
      integer, intent(in ) :: i1, i2
      real,    intent(out) :: aeromode_mass(n_mode)
      real,    intent(out) :: aeromode_dens(n_mode)

      integer              :: n, ns, ne
      real,      parameter :: densmin = 1.e3  ! minimum particle density [ kg/m**3 ]

      ns = max(i1,1)
      ne = min(i2,n_mode)

      do n = ns, ne

         ! dry species mass [ ug / m**3 ]
         aeromode_mass(n) = sum( aerospc_conc( idry_str:idry_end(n), n ) )

         ! add wet and volotile species [ ug / m**3 ]
         aeromode_mass(n) = aeromode_mass(n) + sum( aerospc_conc( 1:iwet_end(n), n ) )

         ! Calculate modal average particle densities [ kg/m**3 ]
         aeromode_dens(n) = max( f6piove9 * aeromode_mass(n) / moment3_conc(n), densmin )

      enddo

      end subroutine getdens

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      Subroutine getdens_conc(aeromode_mass, aeromode_dens,
     &                        conc, moment3_conc, i1, i2)

!     Calculates total aerosol mass and average particle density
!     for each mode specified.

      use const_data, only: f6piove9
      use aero_data,  only: n_mode, mode_map, aer_trac, aer_str, aer_end
      use cgrid_spcs, only: nspcsd

      implicit None

      real,    intent(in ) :: moment3_conc(n_mode)
      real,    intent(in ) :: conc(nspcsd)
      integer, intent(in ) :: i1, i2
      real,    intent(out) :: aeromode_mass(n_mode)
      real,    intent(out) :: aeromode_dens(n_mode)

      integer              :: n, ns, ne
      real,      parameter :: densmin = 1.e3  ! minimum particle density [ kg/m**3 ]

      ns = max(i1,1)
      ne = min(i2,n_mode)

      do n = ns, ne
         ! modal mass [ ug / m**3 ]
         aeromode_mass(n) = sum( conc(aer_str:aer_end), mask=(mode_map==n .and. .not. aer_trac) )

         ! modal mean particle densities [ kg/m**3 ]
         aeromode_dens(n) = max( f6piove9 * aeromode_mass(n) / moment3_conc(n), densmin )
      enddo

      end subroutine getdens_conc

      end module getpar_mod
