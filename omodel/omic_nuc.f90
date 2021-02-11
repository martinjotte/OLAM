!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

!===============================================================================
Module ccnbin_coms

  ! Physical parameters used in CCN nucleation calculations

  real, parameter :: pi1  = 3.141592654
  real, parameter :: pi2  = pi1 * 2.0
  real, parameter :: pio2 = pi1 * 0.5
  real, parameter :: pio4 = pi1 * 0.25
  real, parameter :: pio6 = pi1 / 6.0
  real, parameter :: grav = 9.8

  real, parameter :: cp     = 1004.   ! dry air spec heat at const P [J/(kg K)]
  real, parameter :: rdry   = 287.    ! dry air gas constant [J/(kg K)]
  real, parameter :: rvap   = 461.    ! water vapor gas constant [J/(kg K)]
  real, parameter :: alvl   = 2.50e6  ! latent heat of evaporation [J/kg]
  real, parameter :: rhow   = 1.0e3   ! density of water [kg/m^3]
  real, parameter :: afhh   = 2.25    ! Kumar et al. (2011) Eq. (4) parameter
  real, parameter :: bfhh   = 1.2     ! Kumar et al. (2011) Eq. (4) parameter
  real, parameter :: mbfhh  = -bfhh
  real, parameter :: d2moli = 1.818e9 ! 1 / (2 * diam of h2o molecule) [m^-1]

  real, parameter :: xalphac = 1. ! condens accomodation coeff {was 0.042 in GNF}

  ! Petters and Kreidenweis (2007) used the following fixed values of temperature
  ! and surface tension in determining kappa values for various aerosol types.
  ! Kappa values may need to be re-calibrated if ever these values of tempkap
  ! and sfctension are changed.

  real, parameter :: tempkap = 298.15   ! [K]
  real, parameter :: sfctension = 0.072 ! [N/m]
  real, parameter :: curvcof = 4.0 * sfctension / (rvap * tempkap * rhow)

  real, parameter :: dminprog  =  3.e-6 ! minimum droplet diameter for prognosis
  real, parameter :: dmaxinit  = 10.e-6 ! maximum droplet initial diameter (unless large drydiam)
  real, parameter :: dactivate = 20.e-6 ! minimum droplet diameter to guarantee activation

  ! Memory for CCN bins

  integer :: nccntyp    ! number of CCN types
  integer :: nnuc       ! number of prognosed CCN + GCCN + IFN types
  integer :: nbins      ! total number of CCN bins, summed over all CCN types
  integer :: idust1 = 0 ! index for Fine Dust in ccntyp and nucx arrays
  integer :: idust2 = 0 ! index for Coarse Dust size 2 in ccntyp and nucx arrays
  integer :: idust3 = 0 ! index for Coarse Dust size 3 in ccntyp and nucx arrays
  integer :: idust4 = 0 ! index for Coarse Dust size 4 in ccntyp and nucx arrays
  integer :: isalt  = 0 ! index for sea salt (film mode) in ccntyp and nucx arrays
  integer :: igccnx = 0 ! index for GCCN in nucx arrays
  integer :: iifnx  = 0 ! index for IFN in nucx arrays

  ! Arrays dimensioned by nccntype

  character(40), allocatable :: ccntyp_name(:) ! ccn type name

  integer, allocatable :: nbins_ccntyp(:) ! # of bins to allocate for ccn type

  real, allocatable :: ccntyp_dmin (:) ! ccn type dry diam [m] in smallest-size bin
  real, allocatable :: ccntyp_dmax (:) ! ccn type dry diam [m] in largest-size bin
  real, allocatable :: ccntyp_dmed (:) ! ccn type median dry diam [m]
  real, allocatable :: ccntyp_dsig (:) ! ccn type dry diam shape parameter []
  real, allocatable :: ccntyp_kappa(:) ! ccn type kappa param []
  real, allocatable :: ccntyp_alpha(:) ! Fraction of each CCN type that has IFN properties

  ! Arrays dimensioned by nnuc = nccntyp + 1 (if GCCN prognosed) + 1 (if IFN prognosed)

  real, allocatable ::    rho_nucx(:) ! CCN, GCCN, or IFN dry density [kg/m^3]
  real, allocatable :: bkappa_nucx(:) ! CCN, GCCN, or IFN kappa parameter []
  real, allocatable :: diam_nucx(:,:) ! CCN, GCCN, or IFN wet diameter [m] at 7 relhum
                                      ! values, but diam_nucx(1,:) = dry diameter

  ! Arrays dimensioned by nbins

  character(40), allocatable :: ccnbin_name(:)
  integer,       allocatable :: iccntyp    (:)
  integer,       allocatable :: iccnjbin   (:)
  integer,       allocatable :: ihyg       (:)

  real, allocatable :: bkappa    (:)
  real, allocatable :: relcon_bin(:)
  real, allocatable :: drydiam   (:)
  real, allocatable :: drydiam3  (:)
  real, allocatable :: drydiam11 (:)
  real, allocatable :: drydiam3bk(:)
  real, allocatable :: wetdiam099(:)
  real, allocatable :: wetdiam100(:)
  real, allocatable :: wetdiammin(:)
  real, allocatable :: critdiam  (:)
  real, allocatable :: critsat   (:)
  real, allocatable :: d01sat    (:)
  real, allocatable :: d10sat    (:)
  real, allocatable :: d30sat    (:)
  real, allocatable :: d50sat    (:)
  real, allocatable :: crits     (:)
  real, allocatable :: dmaxbin   (:)

  real, allocatable :: dfrac01   (:)
  real, allocatable :: dfrac10   (:)
  real, allocatable :: dfrac30   (:)
  real, allocatable :: dfrac50   (:)
  real, allocatable :: dfraccrit (:)
end module ccnbin_coms

!=============================================================================

subroutine ccnbin_init()

  use ccnbin_coms, only: pi2, pio6, rdry, rvap, cp, &
                         rhow, tempkap, sfctension, xalphac, &
                         nccntyp, nnuc, nbins, idust1, idust2, idust3, idust4, &
                         isalt, igccnx, iifnx, dmaxbin, wetdiammin, drydiam11, &
                         ccntyp_name, nbins_ccntyp, ccntyp_dmin, ccntyp_dmax, &
                         ccntyp_dmed, ccntyp_dsig, ccntyp_kappa, ccntyp_alpha, &
                         bkappa, relcon_bin, drydiam, drydiam3, wetdiam099, &
                         wetdiam100, critdiam, critsat, d01sat, d10sat, drydiam3bk, &
                         d30sat, d50sat, crits, ccnbin_name, iccntyp, iccnjbin, ihyg, &
                         dminprog, dmaxinit, dactivate, diam_nucx, rho_nucx, bkappa_nucx, &
                         dfrac01, dfrac10, dfrac30, dfrac50, dfraccrit

  use micro_coms,  only: iccn, igccn, iifn
  use misc_coms,   only: io6

  implicit none

  integer :: jbins, ic, jbin, ibin, inewt, i0, inuc
  real :: d01, d10, d30, d50
  real :: ax, dsiglog, dmedlog, diamfac
  real :: aux1, aux2, aux3, aux4, scal
  real :: dbnd1, dbnd2
  real :: ddlog, wetd, dwetd, satk, satkp, satkpp
  real :: tot,tot_ifn, f0
  real, allocatable :: ccn_mass(:), ccn_diam2(:)

  ! The value of ICCN determines a unique set of CCN types to be prognosed.
  ! New sets may be added here by adding ELSEIF sections in the next two
  ! IF-ENDIF blocks and filling with appropriate values.  The first IF-ENDIF
  ! block defines the number of CCN types so that arrays can be dimensioned
  ! to the proper size.

  ! For ICCN = 1, in which case CCN are specified rather than diagnosed,
  ! set nccntyp to 1 so that ccntyp_alpha can be used.

  if (iccn == 1) then
     nccntyp = 1
  elseif (iccn == 2) then
     nccntyp = 1
  elseif (iccn == 3) then
     nccntyp = 1
  elseif (iccn == 4) then
     nccntyp = 3
  elseif (iccn == 5) then
     nccntyp = 8
  elseif (iccn == 6) then
     nccntyp = 9
  endif

  allocate (ccntyp_name (nccntyp))
  allocate (nbins_ccntyp(nccntyp))
  allocate (ccntyp_dmin (nccntyp))
  allocate (ccntyp_dmax (nccntyp))
  allocate (ccntyp_dmed (nccntyp))
  allocate (ccntyp_dsig (nccntyp))
  allocate (ccntyp_kappa(nccntyp))
  allocate (ccntyp_alpha(nccntyp))

  allocate (ccn_mass (nccntyp))
  allocate (ccn_diam2(nccntyp))

  ! Sum number of prognosed CCN + GCCN + IFN types and allocate arrays to include all

  nnuc = 0

  if (iccn >= 2) nnuc = nccntyp
  if (igccn == 2) then
     nnuc = nnuc + 1
     igccnx = nnuc
  endif
  if (iifn == 2) then
     nnuc = nnuc + 1
     iifnx = nnuc
  endif

  allocate (   diam_nucx(7,nnuc))
  allocate (    rho_nucx  (nnuc))
  allocate ( bkappa_nucx  (nnuc))

  ! Fill information for each CCN type.  Kappa should be at least 0.005 for
  ! kappa formulation to apply, and should not exceed 1.28, the highest known
  ! value for any CCN (NaCl).  If kappa is below the minimum threshold, the
  ! particle is assumed to be governed by the adsorption formula.
  ! ccntyp_dmin must not be less than 0.001 micron or greater than 100 microns.

  ! Set zero default values for idust1, idust2, idust3, idust4, isalt, and ccntyp_alpha(:)

  idust1 = 0
  idust2 = 0
  idust3 = 0
  idust4 = 0
  isalt = 0
  ccntyp_alpha(:) = 0.

  if (iccn == 1) then  ! ICCN denotes diagnostic CCN type; ccntyp_alpha(1) is required

     nbins_ccntyp(1) = 0
     ccntyp_name (1) = 'Dust'
     ccntyp_alpha(1) = 0.088

  elseif (iccn == 2) then ! 1 prognostic CCN type: Fine Dust

     idust1 = 1
     ccntyp_name (1) = 'Fine Dust'
     nbins_ccntyp(1) = 20
     ccntyp_dmed (1) = 0.2e-6
     ccntyp_dsig (1) = 2.0
     ccntyp_kappa(1) = 0.00
     diam_nucx (1,1) = 0.2e-6
      rho_nucx   (1) = 2.5e3

  elseif (iccn == 3) then ! 1 prognostic CCN type: Sea Salt.

     isalt = 1
     ccntyp_name (1) = 'Sea Salt'
     nbins_ccntyp(1) = 20
     ccntyp_dmed (1) = 0.15e-6
     ccntyp_dsig (1) = 1.5
     ccntyp_kappa(1) = 1.20
     diam_nucx (1,1) = 0.15e-6
      rho_nucx   (1) = 2.165e3

  elseif (iccn == 4) then ! 3 prognostic CCN types: Fine Dust, Coarse Dust, Sea Salt

     idust1 = 1
     ccntyp_name (1) = 'Fine Dust'
     nbins_ccntyp(1) = 20
     ccntyp_dmed (1) = 0.2e-6
     ccntyp_dsig (1) = 5.0
     ccntyp_kappa(1) = 0.00
     diam_nucx (1,1) = 0.2e-6
      rho_nucx   (1) = 2.5e3

     idust2 = 2
     idust3 = 2
     idust4 = 2
     ccntyp_name (2) = 'Coarse Dust'
     nbins_ccntyp(2) = 20
     ccntyp_dmed (2) = 2.0e-6
     ccntyp_dsig (2) = 2.0
     ccntyp_kappa(2) = 0.00
     diam_nucx (1,2) = 3.0e-6
      rho_nucx   (2) = 2.65e3

     isalt = 3
     ccntyp_name (3) = 'Sea Salt'
     nbins_ccntyp(3) = 20
     ccntyp_dmed (3) = 0.15e-6
     ccntyp_dsig (3) = 1.5
     ccntyp_kappa(3) = 1.20
     diam_nucx (1,3) = 0.15e-6
      rho_nucx   (3) = 2.165e3

  elseif (iccn == 5) then

     ! ICCN = 5 denotes a set of the following 8 CCN types that were selected
     ! for a study of aerosol impacts on tropical cyclones.  Source functions
     ! for all of these CCN are derived from the GEOS-Chem model, while standard
     ! source functions in OLAM for dust and sea salt are applied as well.  This
     ! is ORIGINAL set with 3 coarse dust modes combined, and SO4 not combined
     ! with NO3_NH4.

     idust1 = 1
     ccntyp_name (1) = 'Fine Dust'
     nbins_ccntyp(1) = 20
     ccntyp_dmed (1) = 0.2e-6
     ccntyp_dsig (1) = 2.0
     ccntyp_kappa(1) = 0.00
     diam_nucx (1,1) = 0.2e-6
      rho_nucx   (1) = 2.5e3

     idust2 = 2
     idust3 = 2
     idust4 = 2
     ccntyp_name (2) = 'Coarse Dust'
     nbins_ccntyp(2) = 20
     ccntyp_dmed (2) = 3.0e-6
     ccntyp_dsig (2) = 2.0
     ccntyp_kappa(2) = 0.00
     diam_nucx (1,2) = 3.0e-6
      rho_nucx   (2) = 2.65e3

     isalt = 3
     ccntyp_name (3) = 'Sea Salt'
     nbins_ccntyp(3) = 20
     ccntyp_dmed (3) = 0.15e-6
     ccntyp_dsig (3) = 1.5
     ccntyp_kappa(3) = 1.20
     diam_nucx (1,3) = 0.15e-6
      rho_nucx   (3) = 2.165e3

     ccntyp_name (4) = 'SO4'
     nbins_ccntyp(4) = 20
     ccntyp_dmed (4) = 0.11e-6
     ccntyp_dsig (4) = 1.6
     ccntyp_kappa(4) = 1.0
     diam_nucx (1,4) = 0.11e-6
      rho_nucx   (4) = 2.5e3 !?

     ccntyp_name (5) = 'NO3_NH4'
     nbins_ccntyp(5) = 20
     ccntyp_dmed (5) = 0.11e-6
     ccntyp_dsig (5) = 1.6
     ccntyp_kappa(5) = 1.00
     diam_nucx (1,5) = 0.11e-6
      rho_nucx   (5) = 2.5e3 !?

     ccntyp_name (6) = 'Black Carbon' ! Combined Hydrophobic and Hydrophylic carbon from GEOS-Chem
     nbins_ccntyp(6) = 20
     ccntyp_dmed (6) = 0.02e-6
     ccntyp_dsig (6) = 1.6
     ccntyp_kappa(6) = 0.00
     diam_nucx (1,6) = 0.02e-6
      rho_nucx   (6) = 2.5e3 !?

     ccntyp_name (7) = 'Hydrophilic Organic Carbon'
     nbins_ccntyp(7) = 20
     ccntyp_dmed (7) = 0.09e-6
     ccntyp_dsig (7) = 1.6
     ccntyp_kappa(7) = 0.10
     diam_nucx (1,7) = 0.09e-6
      rho_nucx   (7) = 2.5e3 !?

     ccntyp_name (8) = 'Hydrophobic Organic Carbon'
     nbins_ccntyp(8) = 20
     ccntyp_dmed (8) = 0.09e-6
     ccntyp_dsig (8) = 1.6
     ccntyp_kappa(8) = 0.01
     diam_nucx (1,8) = 0.09e-6
      rho_nucx   (8) = 2.5e3 !?

  elseif (iccn == 6) then

     ! ICCN = 6 denotes a set of the following 9 CCN types that were selected
     ! for a study of aerosol impacts on tropical cyclones.  Source functions
     ! for all of these CCN are derived from the GEOS-Chem model, while standard
     ! source functions in OLAM for dust and sea salt are applied as well.  This
     ! is NEW set with 3 coarse dust modes treated as separate prognostic species,
     ! and SO4 is combined with NO3_NH4 because they have the same mean diameter,
     ! size spectral width, and kappa.

     idust1 = 1
     ccntyp_name (1) = 'Dust1'
     nbins_ccntyp(1) = 20
     ccntyp_dmed (1) = 0.3e-6 * 2. ! estimated, revise later
     ccntyp_dsig (1) = 1.5
     ccntyp_kappa(1) = 0.00
     diam_nucx (1,1) = ccntyp_dmed(1)
      rho_nucx   (1) = 2.5e3

     idust2 = 2
     ccntyp_name (2) = 'Dust2'
     nbins_ccntyp(2) = 20
     ccntyp_dmed (2) = 1.3e-6 * 2. ! estimated, revise later
     ccntyp_dsig (2) = 1.3
     ccntyp_kappa(2) = 0.00
     diam_nucx (1,2) = ccntyp_dmed(2)
      rho_nucx   (2) = 2.65e3

     idust3 = 3
     ccntyp_name (3) = 'Dust3'
     nbins_ccntyp(3) = 20
     ccntyp_dmed (3) = 2.3e-6 * 2. ! estimated, revise later
     ccntyp_dsig (3) = 1.3
     ccntyp_kappa(3) = 0.00
     diam_nucx (1,3) = ccntyp_dmed(3)
      rho_nucx   (3) = 2.65e3

     idust4 = 4
     ccntyp_name (4) = 'Dust4'
     nbins_ccntyp(4) = 20
     ccntyp_dmed (4) = 4.0e-6 * 2.  ! estimated, revise later
     ccntyp_dsig (4) = 1.3
     ccntyp_kappa(4) = 0.00
     diam_nucx (1,4) = ccntyp_dmed(4)
      rho_nucx   (4) = 2.65e3

     isalt = 5
     ccntyp_name (5) = 'Sea Salt'
     nbins_ccntyp(5) = 20
     ccntyp_dmed (5) = 0.15e-6 * 2.
     ccntyp_dsig (5) = 1.5
     ccntyp_kappa(5) = 1.20
     diam_nucx (1,5) = ccntyp_dmed(5)
      rho_nucx   (5) = 2.165e3

     ccntyp_name (6) = 'SO4_NO3_NH4'
     nbins_ccntyp(6) = 20
     ccntyp_dmed (6) = 0.11e-6 * 2.
     ccntyp_dsig (6) = 1.6
     ccntyp_kappa(6) = 1.0
     diam_nucx (1,6) = ccntyp_dmed(6)
      rho_nucx   (6) = 2.0e3

     ccntyp_name (7) = 'Black Carbon' ! Combined Hydrophilic and Hydrophobic carbon from GEOS-Chem
     nbins_ccntyp(7) = 20
     ccntyp_dmed (7) = 0.02e-6 * 2.
     ccntyp_dsig (7) = 1.6
     ccntyp_kappa(7) = 0.00
     diam_nucx (1,7) = ccntyp_dmed(7)
      rho_nucx   (7) = 0.8e3

     ccntyp_name (8) = 'OCPI'
     nbins_ccntyp(8) = 20
     ccntyp_dmed (8) = 0.09e-6 * 2.
     ccntyp_dsig (8) = 1.6
     ccntyp_kappa(8) = 0.10
     diam_nucx (1,8) = ccntyp_dmed(8)
      rho_nucx   (8) = 0.4e3

     ccntyp_name (9) = 'OCPO'
     nbins_ccntyp(9) = 20
     ccntyp_dmed (9) = 0.09e-6 * 2.
     ccntyp_dsig (9) = 1.6
     ccntyp_kappa(9) = 0.0
     diam_nucx (1,9) = ccntyp_dmed(9)
      rho_nucx   (9) = 0.4e3

  endif

  if (iccn >= 2) bkappa_nucx(1:nccntyp) = ccntyp_kappa(1:nccntyp)

  if (igccn == 2) then
       diam_nucx(1,igccnx) = 1.0e-6
        rho_nucx(  igccnx) = 2.5e3
     bkappa_nucx(  igccnx) = 1.20 ! Value for salt
  endif

  if (iifn == 2) then
       diam_nucx(1,iifnx) = 1.0e-6 ! Only used for sedimentation velocity
        rho_nucx(  iifnx) = 2.5e3
     bkappa_nucx(  iifnx) = 0. ! Value for dust
  endif

! Count up total number of bins over all CCN types

  nbins = sum(nbins_ccntyp(1:nccntyp))

  if (nbins < 1) return

  ! Loop over active CCN types

  do ic = 1, nccntyp

     ! The minimum and maximum bin sizes, ccntyp_dmin and ccntyp_dmax, are set to
     ! ccntyp_dmed/diamfac and ccntyp*diamfac, respectively, where diamfac is
     ! given by the following empirical formulas.  Using 20 bins, this results in
     ! the first and last bins each holding about 0.5% or less of the population.

     if     (ccntyp_dsig(ic) < 2.0) then
        diamfac =   0.8 +  5.2 * (ccntyp_dsig(ic) - 1.0)
     elseif (ccntyp_dsig(ic) < 3.0) then
        diamfac =   6.0 + 11.0 * (ccntyp_dsig(ic) - 2.0)
     elseif (ccntyp_dsig(ic) < 4.0) then
        diamfac =  17.0 + 16.0 * (ccntyp_dsig(ic) - 3.0)
     elseif (ccntyp_dsig(ic) < 5.0) then
        diamfac =  33.0 + 24.0 * (ccntyp_dsig(ic) - 4.0)
     elseif (ccntyp_dsig(ic) < 6.0) then
        diamfac =  57.0 + 35.0 * (ccntyp_dsig(ic) - 5.0)
     elseif (ccntyp_dsig(ic) < 7.0) then
        diamfac =  92.0 + 52.0 * (ccntyp_dsig(ic) - 6.0)
     else
        diamfac = 144.0 + 80.0 * (ccntyp_dsig(ic) - 7.0)
     endif

     ccntyp_dmin(ic) = ccntyp_dmed(ic) / diamfac
     ccntyp_dmax(ic) = ccntyp_dmed(ic) * diamfac

     ! Limit ccntyp_dmax to values no larger than 2 microns since GCCN are treated
     ! separately in the model (as the GCCN category).  However, allow dust2 category
     ! to have larger (unlimited) sizes to represent large mineral dust without 
     ! significant solute content.

     ! DISABLE THIS ACTION FOR NOW...EXPLORING ALTERNATIVE (AND MORE CONTROLLED)
     ! SOURCE FUNCTION FOR GCCN.

     ! DISABLED  if (ic /= idust2) then
     ! DISABLED    ccntyp_dmax(ic) = min(ccntyp_dmax(ic), 2.0e-6)
     ! DISABLED  endif

     print*, 'diamfac ',ic,diamfac,ccntyp_dsig(ic)

  enddo

  ! Allocate arrays for all size-bins of all CCN types

  allocate (ccnbin_name(nbins))
  allocate (iccntyp    (nbins))
  allocate (iccnjbin   (nbins))
  allocate (ihyg       (nbins))
  allocate (bkappa     (nbins))
  allocate (relcon_bin (nbins))
  allocate (drydiam    (nbins))
  allocate (drydiam3   (nbins))
  allocate (drydiam11  (nbins))
  allocate (drydiam3bk (nbins))
  allocate (wetdiam099 (nbins))
  allocate (wetdiam100 (nbins))
  allocate (wetdiammin (nbins))
  allocate (critdiam   (nbins))
  allocate (critsat    (nbins))
  allocate (d01sat     (nbins))
  allocate (d10sat     (nbins))
  allocate (d30sat     (nbins))
  allocate (d50sat     (nbins))
  allocate (crits      (nbins))
  allocate (dmaxbin    (nbins))

  allocate (dfrac01    (nbins))
  allocate (dfrac10    (nbins))
  allocate (dfrac30    (nbins))
  allocate (dfrac50    (nbins))
  allocate (dfraccrit  (nbins))

  ! Initialize bins from input "ccntyp" arrays

  jbins = 0
  do ic = 1, nccntyp

     write(io6,'(a,i5)') 'ccntyp ',ic

     if (ccntyp_dsig(ic) < 1.1) then
        stop 'stop: ccntyp_dsig less than 1.1'
     endif

!    ax = log(ccntyp_dmax(ic) / ccntyp_dmin(ic)) / (nbins_ccntyp(ic) * log(2.))
     ax = log(ccntyp_dmax(ic) / ccntyp_dmin(ic)) / real(nbins_ccntyp(ic))

     dbnd1 = ccntyp_dmin(ic)
     tot = 0.
     tot_ifn = 0.

     do jbin = 1,nbins_ccntyp(ic)
        ibin = jbins + jbin

        ccnbin_name(ibin) = ccntyp_name(ic)
        iccntyp    (ibin) = ic
        iccnjbin   (ibin) = jbin

        ! Set up the CCN size spectrum -- geometric progression of diameter

!       dbnd2 = ccntyp_dmin(ic) * 2.0**(jbin * ax)
        dbnd2 = ccntyp_dmin(ic) * exp(jbin * ax)

        drydiam   (ibin) = 0.5 * (dbnd1 + dbnd2)
        drydiam3  (ibin) = drydiam(ibin)**3
        drydiam11 (ibin) = drydiam(ibin) * 1.1
        bkappa    (ibin) = ccntyp_kappa(ic)
        drydiam3bk(ibin) = drydiam3(ibin) * (1. - bkappa(ibin))
        dmaxbin   (ibin) = max(drydiam11(ibin), dmaxinit)
        wetdiammin(ibin) = 1.0025 * drydiam(ibin)

        ! Set up initial aerosol distribution

        dsiglog = log(ccntyp_dsig(ic))
        dmedlog = log(ccntyp_dmed(ic))

        aux1 = (log(dbnd1) - dmedlog)**2
        aux2 = (log(dbnd2) - dmedlog)**2
        aux1 = aux1 / (2. * dsiglog**2)
        aux2 = aux2 / (2. * dsiglog**2)
        aux3 = exp (-aux1)
        aux4 = exp (-aux2)
        aux3 = aux3 / (sqrt(pi2) * dbnd1 * dsiglog)
        aux4 = aux4 / (sqrt(pi2) * dbnd2 * dsiglog)
        relcon_bin(ibin) = 0.5 * (aux4 + aux3) * (dbnd2 - dbnd1)
        tot = tot + relcon_bin(ibin)

        if (dbnd1 >= 0.5e-6) then
           tot_ifn = tot_ifn + relcon_bin(ibin)
        elseif (dbnd2 > 0.5e-6) then
           tot_ifn = tot_ifn + relcon_bin(ibin) &
                   * (dbnd2 - 0.5e-6) / (dbnd2 - dbnd1)
        endif

        dbnd1 = dbnd2

     enddo

     ! For dust, determine ccntyp_alpha parameter based on fraction of CCN in
     ! bins for sizes 0.5 microns and greater.

     if (ic == idust1) then
        ccntyp_alpha(ic) = tot_ifn / tot
        print*, 'ccntyp_alpha for idust1 ',ic,ccntyp_alpha(ic)
     elseif (ic == idust2) then
        ccntyp_alpha(ic) = tot_ifn / tot
        print*, 'ccntyp_alpha for idust2 ',ic,ccntyp_alpha(ic)
     elseif (ic == idust3) then
        ccntyp_alpha(ic) = tot_ifn / tot
        print*, 'ccntyp_alpha for idust3 ',ic,ccntyp_alpha(ic)
     elseif (ic == idust4) then
        ccntyp_alpha(ic) = tot_ifn / tot
        print*, 'ccntyp_alpha for idust4 ',ic,ccntyp_alpha(ic)
     endif

     ! Rescale distribution to compensate for any accumulated errors

     scal = 1. / tot

     tot = 0.

     ccn_mass (ic) = 0.
     ccn_diam2(ic) = 0.

     do jbin = 1,nbins_ccntyp(ic)
        ibin = jbins + jbin

        relcon_bin(ibin) = relcon_bin(ibin) * scal
        tot = tot + relcon_bin(ibin)

        ccn_mass (ic) = ccn_mass (ic) + relcon_bin(ibin) * rho_nucx(ic) * pio6 * drydiam3(ibin)
        ccn_diam2(ic) = ccn_diam2(ic) + relcon_bin(ibin) * drydiam(ibin)**2

        print*, 'filledb ',ic,jbin,ibin,drydiam(ibin),relcon_bin(ibin)

     enddo

     jbins = jbins + nbins_ccntyp(ic)

     print*, 'ic,ccn_mass_ug,ccn_reff ',ic,ccn_mass(ic)*1.e9,0.5 * sqrt(ccn_diam2(ic))

  enddo

  ! Loop over bins and for each one diagnose selected points on the Kohler curve
  ! using calls to subroutine satkap

  do ibin = 1,nbins

     ! Choose initial wetd greater than drydiam to ensure and accelerate
     ! convergence.  For absorption, factor of 1.01 is ok with bkappa in the
     ! range (.001, 1.28) and drydiam in the range (.001, 100.) microns.
     ! For adsorption (flagged by bkappa < .001), base wetd on drydiam.

     if (bkappa(ibin) > 1.e-3) then  ! Absorption (kappa formula)
        wetd = drydiam(ibin) * 1.01
     else                            ! Adsorption (Kumar et al. formula)
        ddlog = log(drydiam(ibin)) + 6.6
        wetd = drydiam(ibin) * (1.33 + .0075 * ddlog**2)
     endif

     ! Diagnose critical diameter and saturation using Newton-Raphson method

     do inewt = 1,50

        call satkap(wetd, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                    satk, satkp, satkpp, 1)

        if (satkpp < 0.) then
           dwetd = -satkp / satkpp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-negative satkpp ',ibin,inewt, &
           drydiam(ibin),wetd,satk,satkp,satkpp
           stop
        endif

     enddo
     critdiam(ibin) = wetd
     critsat (ibin) = satk

     ! Diagnose equilibrium diameter at satenv = 0.99 using Newton-Raphson method

     if (bkappa(ibin) > 1.e-3) then
        wetd = drydiam(ibin) * 1.0001
     else
        ddlog = log(drydiam(ibin)) + 8.6
        wetd = drydiam(ibin) * (1.0005 + .00004 * ddlog**4)
     endif

     do inewt = 1,50
        call satkap(wetd, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 0.99) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp99 ',ibin,inewt,drydiam(ibin),wetd,satk,satkp
           stop
        endif
     enddo
     wetdiam099(ibin) = wetd

     ! Diagnose equilibrium diameter at satenv = 1.00 using Newton-Raphson method

     if (bkappa(ibin) > 1.e-3) then
        wetd = drydiam(ibin) * 1.0001
     else
        ddlog = log(drydiam(ibin)) + 6.6
        wetd = drydiam(ibin) * (1.07 + .003 * ddlog**2)
     endif

     do inewt = 1,50
        call satkap(wetd, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 1.00) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp100 ',ibin,inewt,drydiam(ibin),wetd,satk,satkp
           stop
        endif
     enddo
     wetdiam100(ibin) = wetd

     ! Call to satkap for equilibrium saturation at d01, d10, d30, and d50
     ! between wetdiam100 and critdiam

     d01 = 0.99 * wetdiam100(ibin) + 0.01 * critdiam(ibin)
     d10 = 0.90 * wetdiam100(ibin) + 0.10 * critdiam(ibin)
     d30 = 0.70 * wetdiam100(ibin) + 0.30 * critdiam(ibin)
     d50 = 0.50 * wetdiam100(ibin) + 0.50 * critdiam(ibin)

     call satkap(d01, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                 satk, satkp, satkpp, 0)
     d01sat(ibin) = satk

     call satkap(d10, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                 satk, satkp, satkpp, 0)
     d10sat(ibin) = satk

     call satkap(d30, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                 satk, satkp, satkpp, 0)
     d30sat(ibin) = satk

     call satkap(d50, drydiam(ibin), drydiam3(ibin), bkappa(ibin), &
                 satk, satkp, satkpp, 0)
     d50sat(ibin) = satk

  enddo

  do ibin = 1, nbins
     dfrac01  (ibin) = .01 / (d01sat (ibin) - 1.00)
     dfrac10  (ibin) = .09 / (d10sat (ibin) - d01sat(ibin))
     dfrac30  (ibin) = .20 / (d30sat (ibin) - d10sat(ibin))
     dfrac50  (ibin) = .20 / (d50sat (ibin) - d30sat(ibin))
     dfraccrit(ibin) = .50 / (critsat(ibin) - d50sat(ibin))
  enddo

  ! Rank all bins from lowest to highest critical saturation values.
  ! Save this ordering in ihyg array.

  do ibin = 1,nbins
     crits(ibin) = critsat(ibin)
     ihyg(ibin) = ibin
  enddo

  do ibin = 1,nbins-1
     do jbin = ibin+1,nbins
        if (crits(jbin) < crits(ibin)) then
           f0 = crits(ibin)
           crits(ibin) = crits(jbin)
           crits(jbin) = f0

           i0 = ihyg(ibin)
           ihyg(ibin) = ihyg(jbin)
           ihyg(jbin) = i0
        endif
     enddo
  enddo

! Special print to examine relative hygroscopicity of TC-project aerosols

   do ibin = 1,nbins
      jbin = ihyg(ibin)

      write(6,'(a,3i7,3f10.6)') 'ibin,jbin,iccntyp(jbin) ', &
            ibin,jbin,iccntyp(jbin),crits(ibin),critsat(jbin),drydiam(jbin)*1.e6
  enddo

  ! Loop over nnuc (number of CCN + GCCN + IFN types) and for a chosen particle
  ! diameter characteristic of each type, diagnose selected points on the
  ! Kohler curve using calls to subroutine satkap.  The diagnosed wet diameters
  ! for different relative humidity values are used to compute gravitational
  ! settling velocity for deposition. 

  do inuc = 1,nnuc

     ! Choose initial wetd greater than dry diameter, stored in diam_nucx(1,inuc),
     ! to ensure and accelerate convergence.  For absorption, factor of 1.01
     ! is ok with bkappa_nucx in the range (.001, 1.28) and dry diameter in the
     ! range (.001, 100.) microns.  For adsorption (flagged by bkappa_nucx < .001),
     ! base wetd on dry diameter.

     ! Diagnose equilibrium diameter at satenv = 0.90 using Newton-Raphson method

     if (bkappa_nucx(inuc) > 1.e-3) then
        wetd = diam_nucx(1,inuc) * 1.0001
     else
        wetd = diam_nucx(1,inuc)
        go to 10
     endif

     do inewt = 1,50
        call satkap(wetd, diam_nucx(1,inuc), diam_nucx(1,inuc)**3, bkappa_nucx(inuc), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 0.90) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp90 nucx ',inuc,inewt,diam_nucx(1,inuc),wetd,satk,satkp
           stop
        endif
     enddo
     10 continue
     diam_nucx(2,inuc) = wetd

     ! Diagnose equilibrium diameter at satenv = 0.97 using Newton-Raphson method

     if (bkappa_nucx(inuc) > 1.e-3) then
        wetd = diam_nucx(1,inuc) * 1.0001
     else
        ddlog = log(diam_nucx(1,inuc)) + 8.6
        wetd = diam_nucx(1,inuc) * (1.0005 + .00004 * ddlog**4)
     endif

     do inewt = 1,50
        call satkap(wetd, diam_nucx(1,inuc), diam_nucx(1,inuc)**3, bkappa_nucx(inuc), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 0.97) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp97 nucx ',inuc,inewt,diam_nucx(1,inuc),wetd,satk,satkp
           stop
        endif
     enddo
     diam_nucx(3,inuc) = wetd

     ! Diagnose equilibrium diameter at satenv = 0.99 using Newton-Raphson method

     if (bkappa_nucx(inuc) > 1.e-3) then
        wetd = diam_nucx(1,inuc) * 1.0001
     else
        ddlog = log(diam_nucx(1,inuc)) + 8.6
        wetd = diam_nucx(1,inuc) * (1.0005 + .00004 * ddlog**4)
     endif

     do inewt = 1,50
        call satkap(wetd, diam_nucx(1,inuc), diam_nucx(1,inuc)**3, bkappa_nucx(inuc), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 0.99) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp99 nucx ',inuc,inewt,diam_nucx(1,inuc),wetd,satk,satkp
           stop
        endif
     enddo
     diam_nucx(4,inuc) = wetd

     ! Diagnose equilibrium diameter at satenv = 0.997 using Newton-Raphson method

     if (bkappa_nucx(inuc) > 1.e-3) then
        wetd = diam_nucx(1,inuc) * 1.0001
     else
        ddlog = log(diam_nucx(1,inuc)) + 8.6
        wetd = diam_nucx(1,inuc) * (1.0005 + .00004 * ddlog**4)
     endif

     do inewt = 1,50
        call satkap(wetd, diam_nucx(1,inuc), diam_nucx(1,inuc)**3, bkappa_nucx(inuc), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 0.997) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp997 nucx ',inuc,inewt,diam_nucx(1,inuc),wetd,satk,satkp
           stop
        endif
     enddo
     diam_nucx(5,inuc) = wetd

     ! Diagnose equilibrium diameter at satenv = 0.999 using Newton-Raphson method

     if (bkappa_nucx(inuc) > 1.e-3) then
        wetd = diam_nucx(1,inuc) * 1.0001
     else
        ddlog = log(diam_nucx(1,inuc)) + 8.6
        wetd = diam_nucx(1,inuc) * (1.0005 + .00004 * ddlog**4)
     endif

     do inewt = 1,50
        call satkap(wetd, diam_nucx(1,inuc), diam_nucx(1,inuc)**3, bkappa_nucx(inuc), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 0.999) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp999 nucx ',inuc,inewt,diam_nucx(1,inuc),wetd,satk,satkp
           stop
        endif
     enddo
     diam_nucx(6,inuc) = wetd

     ! Diagnose equilibrium diameter at satenv = 1.00 using Newton-Raphson method

     if (bkappa_nucx(inuc) > 1.e-3) then
        wetd = diam_nucx(1,inuc) * 1.0001
     else
        ddlog = log(diam_nucx(1,inuc)) + 6.6
        wetd = diam_nucx(1,inuc) * (1.07 + .003 * ddlog**2)
     endif

     do inewt = 1,50
        call satkap(wetd, diam_nucx(1,inuc), diam_nucx(1,inuc)**3, bkappa_nucx(inuc), &
                    satk, satkp, satkpp, 1)

        if (satkp > 0.) then
           dwetd = - (satk - 1.00) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp100 nucx ',inuc,inewt,diam_nucx(1,inuc),wetd,satk,satkp
           stop
        endif
     enddo
     diam_nucx(7,inuc) = wetd

  enddo

  do inuc = 1,nnuc
     write(6,'(a,i5,7f10.2)') 'diam_nucx ',inuc,1.e6 * diam_nucx(1,inuc), &
                                   diam_nucx(2,inuc) / diam_nucx(1,inuc), &
                                   diam_nucx(3,inuc) / diam_nucx(1,inuc), &
                                   diam_nucx(4,inuc) / diam_nucx(1,inuc), &
                                   diam_nucx(5,inuc) / diam_nucx(1,inuc), &
                                   diam_nucx(6,inuc) / diam_nucx(1,inuc), &
                                   diam_nucx(7,inuc) / diam_nucx(1,inuc)
  enddo

end subroutine ccnbin_init

!=============================================================================

! Subroutine CCNBIN simulates growth of water droplets on CCN in an environment
! of increasing saturation.  CCN are sorted into bins, each of which represents
! a particular hygroscopicity (represented by the kappa parameter; Petters and
! Kreidenweis, 2007) and diameter.  For bins whose critical saturation is
! exceeded by environmental saturation, the wetted CCN grow freely as small
! cloud droplets, while most of the remainder (those that have not crossed the
! energy barrier represented by the Kohler curve) remain in stable equilibrium
! as unactivated haze particles.  An exception is made for very large CCN that
! can take on substantial water before reaching the critical diameter at the
! peak of the Kohler curve; these droplets are allowed to grow freely in the
! manner of activated cloud droplets once they have exceeded a specified
! threshold diameter.  The particular set of bins that develop into cloud
! droplets is determined by the maximum supersaturation attained and is the
! principal result of the simulation, along with the amount of water condensed
! onto the CCN.  The simulation is stopped automatically when the point is
! reached where enough new cloud drops have developed and grown to sufficient
! size that they take on vapor at a rate that exceeds environmental production
! of supersaturation.  However, a maximum duration of the simulation is
! specified externally based on the particular application (such as the length
! of a host model timestep) and need not be long enough for saturation to reach
! a maximum and begin to decrease.

! This code is adapted from earlier versions of program/subroutine PARCEL:
! Michael Sabin, NCAR (5/20/88)
! Larry Miloshevich, NCAR (10/16/90)
! Graham Feingold, CIRES, CIRA (1991-2000)
! Bob Walko, CSU (2000)
! Steve Saleeby, CSU (2000-2011)
! Dan Ward, CSU (2009)
! Dave Lerach, CSU (2010)

! The CSU developers used PARCEL to generate nucleation tables for RAMS (and
! later, OLAM).  The Ward/Lerach developments introduced the kappa parameter as
! a measure of hygroscopicity (Petters and Kreidenweis, 2007), while the
! Walko/Saleeby version retained the original form without kappa.

! The present CCNBIN code uses the kappa form and was adapted by Bob Walko
! (2014) from the earlier Walko, Saleeby, and Ward/Lerach versions.  Major
! restructuring and modification were done to make the code callable as
! subroutines in OLAM.  This enables CCN activation to be represented during
! simulation time without resorting to lookup tables.  The added computational
! cost of a simulation is sometimes outweighed by the flexibility of
! representing any mixture of different CCN types (kappas) and diameters,
! including multiple size modes.  The present form of the code could also be
! used to generate lookup tables, given an appropriate driver program, although
! it is impractical for lookup tables to exploit the full range of CCN size and
! kappa mixtures that can occur.

! CCNBIN differs substantially from PARCEL in several ways:
! (1) Prognostic equations for atmospheric vapor and droplet water content are
!     reformulated as specific densities (relative to combined mass of all
!     atmospheric constituents) instead of absolute densities (relative to
!     volume), which simplifies their form and reduces computation.
!     (In PARCEL, some terms that involved atmospheric density tendency
!     cancelled anyway and thus did not need to be evaluated.)
! (2) Subroutine ZZBREN, which used the bisection rule to find solutions
!     involving the Kohler curve, was replaced by the Newton-Raphson method,
!     which converges faster.
! (3) Subroutine vode, which solved a system of ODEs in PARCEL, was found to
!     be too slow, especially with a large number of bins, so it was replaced
!     by a more direct in-line solver that exploits the known behavior of the
!     physical system to gain efficiency.  Vode also would have required 
!     restructuring (e.g., replacement of common blocks) for multi-threaded
!     applications of OLAM.
! (4) Very small (sub-micron size) unactivated droplets are in highly stable
!     equilibrium with environmental saturation and respond extremely rapidly
!     to changes in saturation with a small change in diameter.  It is more
!     computationally efficient to diagnose, rather than prognose, droplet
!     diameters in this regime in order to avoid the need for extremely short
!     prognostic timesteps (or frequent evaluations of droplet surface
!     saturation as was performed in PARCEL via subroutine vode).  Droplets of
!     larger diameter are more efficiently prognosed.  We thus use a threshold
!     wet diameter (of a few microns) to determine whether to diagnose or
!     prognose a droplet.  In addition, diagnosis of unactivated droplet
!     diameter is approximated by linear interpolation between selected droplet
!     wet diameters where equilibrium saturation is first pre-diagnosed
!     iteratively.  Interpolation is sufficiently accurate and is far more
!     efficient than more precise diagnosis arrived at iteratively.

subroutine ccnbin(iw0, k, ntim, timespan, &
                  tempk0, press0, rhoa0, rhov0, satenv0, &
                  tempk1, press1, rhoa1, rhov1, satenv1, &
                  con_cp, con_ccny, cactivated, rho_wbc)

  use ccnbin_coms, only: pi2, pio2, pio4, pio6, rdry, rvap, cp, &
                         rhow, tempkap, sfctension, xalphac, bkappa, &
                         nccntyp, nbins, nbins_ccntyp, drydiam, &
                         relcon_bin, drydiam3, drydiam3bk, wetdiammin, &
                         wetdiam099, wetdiam100, critdiam, critsat, d01sat, &
                         d10sat, d30sat, d50sat, ihyg, d2moli, drydiam11, &
                         dminprog, dmaxbin, dactivate, curvcof, mbfhh, afhh, &
                         dfrac01, dfrac10, dfrac30, dfrac50, dfraccrit
  use therm_lib, only: rhovsl_inv

  implicit none

  integer, intent(in) :: iw0         ! OLAM grid horizontal index
  integer, intent(in) :: k           ! OLAM grid vertical index
  integer, intent(in) :: ntim        ! # of timesteps to be run in parcel model

  real, intent(in ) :: timespan
  real, intent(in ) ::  tempk0, tempk1   ! initial & final air temperatures [K]
  real, intent(in ) ::  press0, press1   ! initial & final pressures [Pa]
  real, intent(in ) ::   rhoa0,  rhoa1   ! initial & final air densities   [kg_a/m^3]
  real, intent(in ) ::   rhov0,  rhov1   ! initial & final vapor densities [kg_v/m^3]
  real, intent(in ) :: satenv0, satenv1  ! initial & final saturation (env) [ ]
  real, intent(in ) :: con_cp            ! concentration of cloud + pristine ice [#/m^3]
  real, intent(in ) :: con_ccny(nccntyp) ! prognosed concentration of each CCN type [#/m^3]
  real, intent(out) :: cactivated, rho_wbc

! Local arrays

  logical :: iprognose(nbins)  ! bin prognose flag
  logical :: skipbin  (nbins)  ! bin skip flag

  real :: con_bin (nbins)
  real :: wetdiam (nbins)
  real :: wetdiam3(nbins)

  ! Local variables

  integer :: ibin, ic, jbin, itim, jbins
  real :: dtim, tim, fractim, diagdiam, frac
  real :: tot, dtim4, dtim12
  real :: rho_wb
  real :: satenv, satenvmax
  real :: tempk, press, rhoa, rhov, rhov_remain
  real :: satk, esat
  real :: rlv, pratio, rlambda, dvn, rkan, fdcfac, fkafac
  real :: rad, rad2, dsvi, fksai, t1, t2, t3, tn, td, t1a, t3a
  real :: curvterm, actterm

  real, parameter :: afhh0    = afhh * d2moli ** mbfhh
  real, parameter :: onethird = 1. / 3.
  real, parameter :: dminprogeps = nearest(dminprog, -1.0)
  real, parameter :: fdc0 = sqrt( pi2 / rvap ) / xalphac
  real, parameter :: fka0 = sqrt( pi2 / rdry ) / (0.96 * cp)
  real, external  :: few, flhv

  cactivated = 0.
  rho_wbc    = 0.

  ! Return if we are not near saturation

  if (satenv1 < 0.99 .and. satenv0 < 0.99) return

  ! If critical saturation of all bins exceeds ending environmental saturation, return

  if (all(critsat(1:nbins) > satenv1)) return

  ! Fill con_bin values for current CCNBIN integration

  jbins = 0
  do ic = 1, nccntyp
     do jbin = 1,nbins_ccntyp(ic)
        ibin = jbins + jbin
        con_bin(ibin) = relcon_bin(ibin) * con_ccny(ic)
     enddo
     jbins = jbins + nbins_ccntyp(ic)
  enddo

  ! If cloud droplets and/or pristine ice already exist in this parcel, assume
  ! that they have nucleated from the most hygroscopic bins of CCN and that
  ! those bins are therefore not currently available for nucleation.  Set bin 
  ! concentrations to zero, beginning with the most hygroscopic, until reaching
  ! the number of already-activated hydrometeors input to this routine.  (If 
  ! pristine ice has nucleated from IN, its concentrations will be insignificant
  ! compared to CCN, so the zeroing will have negligible effect.  If pristine
  ! ice was homogeneously frozen from cloud droplets, it will be at high
  ! concentration but zeroing of CCN bins is appropriate because of the cloud
  ! droplet origin of the pristine ice.)

  tot = 0.

  do ibin = 1,nbins
     jbin = ihyg(ibin)
     tot = tot + con_bin(jbin)

     if (tot < con_cp) then
        con_bin(jbin) = 0.
     else
        con_bin(jbin) = max(0.,tot - con_cp)
        exit
     endif
  enddo

  ! Skip bins that are already nucleated or have negligible CCN concentrations

  skipbin(1:nbins) = ( con_bin(1:nbins) < 1.e2 )

  if (all(skipbin(1:nbins))) return

  ! Initialize quantities before time loop

  iprognose(:) = .false.
  satenvmax = 0.

  rlv    = FLHV(tempk0)
  rho_wb = 0.0

  dtim   = timespan / ntim
  dtim4  = dtim *  4.
  dtim12 = dtim * 12.
  tim    = 0.

  wetdiam (1:nbins) = drydiam (1:nbins)
  wetdiam3(1:nbins) = drydiam3(1:nbins)

  ! Start the time loop

  do itim = 1, ntim

     tim = tim + dtim
     fractim = tim / timespan

     ! Sum water condensate density over all bins [kg_wb/m^3]

     if (itim > 1) then
        rho_wb = rhow * pio6 &
              * sum(con_bin(1:nbins) * (wetdiam3(1:nbins) - drydiam3(1:nbins)))
     endif

     ! Update environmental variables

     press = press0 + fractim * (press1 - press0)
     rhoa  = rhoa0  + fractim * ( rhoa1 -  rhoa0)
     rhov  = rhov0  + fractim * ( rhov1 -  rhov0)
     tempk = tempk0 + fractim * (tempk1 - tempk0) + rlv * rho_wb / (rhoa * cp)

     ! Vapor density remaining after bin water uptake [kg_v/m^3]

     rhov_remain = rhov - rho_wb

     ! Saturation ratio

     satenv = rhov_remain * rhovsl_inv(tempk-273.15)

     ! Limit satenv on first 4 iterations to smooth sudden onset of latent heating

     if (itim == 1) then
        satenv = min(satenv,0.99001)
     elseif (itim == 2) then
        satenv = min(satenv,0.995)
     elseif (itim == 3) then
        satenv = min(satenv,0.998)
     elseif (itim == 4) then
        satenv = min(satenv,1.000)
     elseif (itim == 5) then
        satenv = min(satenv,1.002)
     endif

     ! If running ccnbin as part of OLAM (or other model), halt integration of
     ! ccnbin once satenv has begun to decrease in time (due to newly activated
     ! droplets consuming excess vapor faster than environmental production).  This
     ! prevents newly activated droplets from consuming excess moisture (aided by
     ! their solute effect) before they are introduced to the cloud droplet
     ! population which has no solute effect.

     satenvmax = max(satenvmax,satenv)
     if (satenv <= nearest(satenvmax,-1.)) exit

     ! Latent heat of vaporization, saturated vapor pressure

     rlv  = FLHV(tempk)
     esat = few (tempk)

     ! Atmospheric properties used for computing droplet vapor diffusional growth

     pratio  = 1.01325e5 / press
     rlambda = 6.6e-8 * pratio * (tempk / 293.15)   ! mean free path P&K 10-140
     dvn     = 0.211e-4 * pratio * (tempk / 273.15)**1.94   ! diffusivity P&K 13-3 but MKS
     rkan    = 4.38e-3 + 7.118e-5 * tempk   ! therm cond P&K 13-16 but in MKS with K

     fdcfac  = fdc0 * dvn  /  SQRT(tempk)
     fkafac  = fka0 * rkan / (SQRT(tempk) * rhoa)

     t1a = rlv * rhow / tempk
     t2  = rlv / (tempk * rvap) - 1.   ! P&K 13-28 denom term 2B
     t3a = rhow * tempk * rvap / esat

     if (satenv >= 0.99) then

        ! For bins that have not yet been prognosed, diagnose diameter by
        ! interpolation between pre-determined diameters

        do ibin = 1,nbins
           if ((.not. iprognose(ibin)) .and. (.not. skipbin(ibin))) then

              if (satenv < 1.00) then
                 frac = (satenv - 0.99) / 0.01
                 diagdiam = (       frac) * wetdiam100(ibin) &
                          + (1.00 - frac) * wetdiam099(ibin)
              elseif (satenv < d01sat(ibin)) then
                 frac = (satenv - 1.00) * dfrac01(ibin)
                 diagdiam = (       frac) * critdiam(ibin) &
                          + (1.00 - frac) * wetdiam100(ibin)
              elseif (satenv < d10sat(ibin)) then
                 frac = (satenv - d01sat(ibin)) * dfrac10(ibin)
                 diagdiam = (.01 + frac) * critdiam(ibin) &
                          + (.99 - frac) * wetdiam100(ibin)
              elseif (satenv < d30sat(ibin)) then
                 frac = (satenv - d10sat(ibin)) * dfrac30(ibin)
                 diagdiam = (.10 + frac) * critdiam(ibin) &
                          + (.90 - frac) * wetdiam100(ibin)
              elseif (satenv < d50sat(ibin)) then
                 frac = (satenv - d30sat(ibin)) * dfrac50(ibin)
                 diagdiam = (.30 + frac) * critdiam(ibin) &
                          + (.70 - frac) * wetdiam100(ibin)
              elseif (satenv < critsat(ibin)) then
                 frac = (satenv - d50sat(ibin)) * dfraccrit(ibin)
                 diagdiam = (.50 + frac) * critdiam(ibin) &
                          + (.50 - frac) * wetdiam100(ibin)
              else
                 diagdiam = max(dminprog,critdiam(ibin))
              endif

              diagdiam = max(diagdiam, drydiam11(ibin))

              if ( diagdiam >= dminprogeps ) then
                 iprognose(ibin) = .true.
                 wetdiam  (ibin) = min(diagdiam, dmaxbin(ibin))
                 wetdiam3 (ibin) = wetdiam(ibin)**3
              endif

           endif  ! iprognose(ibin)
        enddo

     endif  ! satenv > 0.99

     ! Prognose wet diameter for bins flagged for prognosis

     do ibin = 1,nbins
        if (iprognose(ibin)) then

           ! Time derivative of droplet diameter due to vapor diffusion.

           rad  = wetdiam(ibin) * 0.5
           rad2 = rad * rad / (rlambda + rad)

           dsvi  = (rad2 + fdcfac) / (rad * dvn )  ! inv of modified diffusivity P&K 13-13
           fksai = (rad2 + fkafac) / (rad * rkan)  ! inv of modified therm cond P&K 13-20

           t1 = t1a * fksai ! P&K 13-28 denom term 2A
           t3 = t3a * dsvi  ! P&K 13-28 denom term 1

           if (bkappa(ibin) > 1.e-4) then

              curvterm = exp(curvcof / wetdiam(ibin))
              actterm  = (wetdiam3(ibin) - drydiam3  (ibin)) &
                       / (wetdiam3(ibin) - drydiam3bk(ibin))
              satk     = curvterm * actterm

           else

              curvterm = curvcof / wetdiam(ibin)
              ! exp/log is often faster than ** for repetitive tasks
              ! actterm  = afhh0 * (wetdiam(ibin) - drydiam(ibin)) ** mbfhh
              actterm  = afhh0 * exp( log(wetdiam(ibin) - drydiam(ibin)) * mbfhh)
              satk     = exp(curvterm - actterm)

           endif

           tn = satenv - satk
           td = t1 * t2 + t3

           ! Update of wet diameter.  To improve computational stability,
           ! use volume tendency of growing drop and diameter tendency of
           ! shrinking drop.

           if (tn > 0.) then
              wetdiam(ibin) = ( wetdiam3(ibin) &
                              + dtim12 * wetdiam(ibin) * tn / td ) ** onethird
           else

              wetdiam(ibin) = max( wetdiammin(ibin), &
                                   wetdiam(ibin) + dtim4 * tn / (wetdiam(ibin) * td) )
           endif
           wetdiam3(ibin) = wetdiam(ibin)**3

        endif  ! iprognose(ibin)
     enddo  ! ibin

  enddo    ! itim

  ! Sum number concentration [#_droplets/kg_a] and bulk specific density
  ! [kg_wb/kg_a] of droplets over all ACTIVATED bins.  Define "activated" ccn
  ! as those whose critical saturation has been exceeded by environmental
  ! saturation and/or whose wet diameter exceeds dactivate.  Even though very
  ! large CCN are often not truly activated (because there was insufficient
  ! time for them to grow to their critical diameters), they are still more
  ! amenable to continued growth than smaller truly-activated CCN and thus
  ! qualify to become cloud droplets.

  do ibin = 1, nbins

     if (critsat(ibin) < satenvmax .or. wetdiam(ibin) > dactivate) then
        cactivated = cactivated + con_bin(ibin)
        rho_wbc = rho_wbc + rhow * pio6 &
                          * con_bin(ibin) * (wetdiam3(ibin) - drydiam3(ibin))
     endif

  enddo

end subroutine ccnbin

!=============================================================================

real function FEW (t)

! Function FEW calculates saturation vapor pressure [Pa] over water for
! a given temperature [K].

  implicit none

  real, intent(in) :: t
  real             :: ar, ar1, ari, br, cr, dw, er

  real, parameter  :: ts  = 373.16
  real, parameter  :: ps  = 101324.6
  real, parameter  :: tsi = 1. / ts
  real, parameter  :: l10 = log(10.)
  real, parameter  :: br0 = -7.90298   * l10
  real, parameter  :: cr0 =  5.02808
  real, parameter  :: dw0 = -1.3816e-7 * l10
  real, parameter  :: dw1 = 11.344     * l10
  real, parameter  :: er0 =  8.1328e-3 * l10
  real, parameter  :: er1 = -3.49149   * l10
  real, parameter  :: fw0 = ps * exp(-dw0-er0)

  ar  = ts / t
  ar1 = ar - 1.
  ari = tsi * t

  br  = br0 * ar1
  cr  = cr0 * log( ar )
  dw  = dw0 * exp( dw1 * ar1 * ari)
  er  = er0 * exp( er1 * ar1      )

  FEW = fw0 * exp( br + cr + dw + er )

!MJO:  Terms of the form 10.**x were replaced with exp(log(10) * x)
!MJO:  for speed.

end function few

!=============================================================================

real function FLHV (t)

!  Purpose:  This function calculates the latent heat of vaporization for
!            water as per P&K eqn. 4-85a. The result is then converted to
!            J/kg.  T is deg K.

  implicit none

  real, intent(in) :: t
  real             :: rga, rlhv

  rga  = 0.167 + (3.67e-4 * t)
  rlhv  = 597.3 * (273.15 / t)**rga
  FLHV = rlhv / 2.38844e-04

!Bob:  A reasonable quadratic fit (in J/kg) is:
!Bob:  flhv = 2500795. + tempc * (2.3 * tempc - 2453.)
!Bob:  (Result is within 300 J/kg of above formula between +/- 40 deg C)

end function flhv

!=============================================================================

! Subroutine satkap evaluates the saturation vapor mixing ratio at the surface
! of an aerosol particle and its first and second derivatives with respect to
! wet diameter.  For bkappa > 1.e-4, the aerosol is assumed to be sufficiently
! soluble that wet growth is described by Eqn (6) of Petters and Kreidenweis
! (2007) using the kappa parameter.  For bkappa < 1.e-4, the aerosol is assumed
! to be insoluble with wet growth occuring by adsorption and described by
! Eqn (4) of Kumar et al. (2011).

subroutine satkap(wetdiam, drydiam, drydiam3, bkappa, satk, satkp, satkpp, ids)

  use ccnbin_coms, only: rvap, rhow, tempkap, sfctension, &
                         afhh, bfhh, mbfhh, d2moli, curvcof

  implicit none

  real, intent(in)  :: drydiam, wetdiam, drydiam3, bkappa
  real, intent(out) :: satk, satkp, satkpp

  integer, intent(in) :: ids

  real(kind=8) :: wetdiam2, wetdiam3, curvterm
  real(kind=8) :: actdenom, actterm, term1, term2, term3, term4, term5
  real(kind=8) :: curv1, curv2, act1, act2, act3

  if (bkappa > 1.e-4) then

     wetdiam2 = wetdiam * wetdiam
     wetdiam3 = wetdiam2 * wetdiam

     curvterm = exp(curvcof / wetdiam)
     actdenom = wetdiam3 - drydiam3 * (1. - bkappa)
     actterm = (wetdiam3 - drydiam3) / actdenom

     satk = real(actterm * curvterm)

     if (ids == 0) return

     term1 = 3. * curvterm * wetdiam2 / actdenom
     term2 = term1 * actterm
     term3 = curvcof * satk / wetdiam2

     satkp = real(term1 - term2 - term3)

     term4 = (6. / actdenom) * (1. - actterm) &
           * (wetdiam - curvcof - 3.0 * wetdiam2**2 / actdenom)

     term5 = (actterm * curvcof / wetdiam3) * (2. + curvcof / wetdiam)

     satkpp = curvterm * real(term4 + term5)

  else

     curv1 = curvcof / wetdiam
     act1 = (wetdiam - drydiam) * d2moli
     act2 = afhh * act1 ** mbfhh

     satk = exp(real(curv1 - act2))

     if (ids == 0) return

     curv2 = curv1 / wetdiam
     act3 = act2 * bfhh * d2moli / act1

     satkp = satk * (-curv2 + act3)

     satkpp = satkp * (-curv2 + act3) &
            + satk * (2. * curv2 / wetdiam - act3 * (bfhh + 1.) * d2moli / act1)

  endif

end subroutine satkap

!=============================================================================

subroutine cldnuc(iw0,lpw0,dtli0,nbincall, &
                  rx,cx,qr,qx,con_ccnx,con_gccnx,rhov,rhoi,rhoa,press0, &
                  tair,tairc,wc0,rhovslair,rnuc_vc,rnuc_vd,cnuc_vc,cnuc_vd)

use micro_coms,  only: jnmb, emb0, emb0i, emb1, mza0, ncat, &
                       igccn, rxmin
use ccnbin_coms, only: nccntyp
use misc_coms,   only: io6, dtlm
use consts_coms, only: rvap, grav, alvl, cp, cliq, alli, eps_vapi
use therm_lib,   only: rhovsl

implicit none

integer, intent(in) :: iw0
integer, intent(in) :: lpw0

real, intent(in) :: dtli0

integer, intent(inout) :: nbincall

real, intent(inout) :: rx(mza0,ncat)
real, intent(inout) :: cx(mza0,ncat)
real, intent(inout) :: qr(mza0,ncat)
real, intent(inout) :: qx(mza0,ncat)

real, intent(inout) :: con_ccnx (mza0,nccntyp)
real, intent(inout) :: con_gccnx(mza0)

real, intent(inout) :: rhov  (mza0)
real, intent(in) :: rhoi     (mza0)
real, intent(in) :: press0   (mza0)
real, intent(in) :: tair     (mza0)
real, intent(in) :: tairc    (mza0)
real, intent(in) :: wc0      (mza0)
real, intent(in) :: rhovslair(mza0)

real, intent(inout) :: cnuc_vc(mza0)
real, intent(inout) :: cnuc_vd(mza0)
real, intent(inout) :: rnuc_vc(mza0)
real, intent(inout) :: rnuc_vd(mza0)

real, intent(in) :: rhoa(mza0)

integer :: k

real :: excessrhov, excessnum, rxnuc
real :: rnuc_vc_min, rnuc_vc_max
!real :: a1inv,w_pseudo
real :: cactivated, rx_wbc

! CCNBIN variables

  real ::  tempkA,  tempkB ! initial & final air temperatures [K]
  real ::  pressA,  pressB ! initial & final pressures [Pa]
  real ::   rhoaA,   rhoaB ! initial & final air densities   [kg_a/m^3]
  real ::   rhovA,   rhovB ! initial & final vapor densities [kg_v/m^3]
  real :: satenvA, satenvB
  real :: timespan, con_cp
  real :: con_ccny(nccntyp)

  integer, parameter :: ntim = 50

  ! Loop over all grid levels in column

  nbincall = 0

  do k = lpw0,mza0

     ! Diagnose excess of vapor density over saturation and number of drizzle
     ! particles it can nucleate.  Cycle if not > 0.

     excessrhov = rhov(k) - 1.00001 * rhovslair(k) ! x rhoa
     if (excessrhov <= 0.) cycle                   ! x rhoa
     excessnum = excessrhov * emb0i(8)

     ! Cotton (2016) pointed out that with GCCN, the solute effect accelerates
     ! droplet growth sufficiently (even at zero supersaturation) to produce
     ! fairly large droplets in a fairly short amount of time.  Accordingly, we
     ! choose to nucleate GCCN directly to the drizzle category, and to do so
     ! before CCN nucleation to cloud droplets.  With this choice, GCCN are NOT
     ! included in the CCNBIN computation of nucleation (of cloud).

     ! If GCCN are not prognosed, repeated nucleation over successive timesteps
     ! could lead to an over-abundance of drizzle particles unless the number
     ! nucleated are limited in some manner based on the nucleation that
     ! already occurred on previous time steps.  The number concentration of
     ! drizzle [cx(k,8)] is used as a proxy for the number of GCCN previously
     ! nucleated in the limiter below.

     ! if GCCN are prognosed, they are immediately scavenged upon nucleation to
     ! drizzle.  This is necessary because drizzle droplets can be far more
     ! numerous than GCCN when they arise from cloud droplet collisions, GCCN
     ! scavenging from drizzle collisions and sedimentation is therefore not
     ! feasible.

     if (igccn == 2) then ! GCCN are prognosed; nucleate and scavenge GCCN
        cnuc_vd(k) = min(excessnum, con_gccnx(k))
        con_gccnx(k) = con_gccnx(k) - cnuc_vd(k)
     else                 ! GCCN NOT prognosed; limit nucleation by drizzle present
        cnuc_vd(k) = min(excessnum, max(0., con_gccnx(k) - cx(k,8)))
     endif

     rnuc_vd(k) = cnuc_vd(k) * emb0(8)

     ! Update vapor mass and drizzle concentration and mass from nucleation

     rhov(k) = rhov(k) - rnuc_vd(k)
     cx(k,8) = cx(k,8) + cnuc_vd(k) ! (#/m^3)  x rhoa
     rx(k,8) = rx(k,8) + rnuc_vd(k) ! (kg/m^3) x rhoa
     qr(k,8) = qr(k,8) + rnuc_vd(k) * (tairc(k) * cliq + alli) ! (x rhoa)

     if (rx(k,8) >= rxmin(8)) qx(k,8) = qr(k,8) / rx(k,8)

     ! Now that GCCN/DRIZZLE nucleation is complete, diagnose remaining excess
     ! of vapor density over saturation; cycle if not > 0.

     excessrhov = rhov(k) - 1.00001 * rhovslair(k) ! x rhoa
     if (excessrhov <= 0.) cycle                   ! x rhoa

     ! Check whether cloud droplet number (and CCN number) are prognosed

     if (jnmb(1) == 4) then

        ! If NOT prognosing number concentration of cloud, use con_ccnx(k,:) as
        ! the potential number concentration of cloud droplets.  Limit number
        ! of newly nucleated cloud droplets by number of droplets already
        ! present and by available supersaturation vapor (subject to the
        ! minimum cloud droplet mass).

        if (nccntyp == 1) then
           cnuc_vc(k) = max(0., con_ccnx(k,1) - cx(k,1))
        else
           cnuc_vc(k) = max(0., sum(con_ccnx(k,1:nccntyp)) - cx(k,1))
        endif

        cnuc_vc(k) = min(cnuc_vc(k), excessrhov * emb0i(1))

        ! Assume that half of available supersaturation vapor is transferred to
        ! newly-nucleated cloud droplets, subject however to limits on cloud 
        ! droplet size and number nucleated.

        rnuc_vc_min = cnuc_vc(k) * emb0(1)
        rnuc_vc_max = cnuc_vc(k) * emb1(1)

        rnuc_vc(k) = max(rnuc_vc_min, min(rnuc_vc_max, 0.5 * excessrhov))

        ! Update vapor mass and cloud concentration and mass from nucleation

        rhov(k) = rhov(k) - rnuc_vc(k)
        cx(k,1) = cx(k,1) + cnuc_vc(k) ! (#/m^3)  x rhoa
        rx(k,1) = rx(k,1) + rnuc_vc(k) ! (kg/m^3) x rhoa
        qr(k,1) = qr(k,1) + rnuc_vc(k) * (tairc(k) * cliq + alli) ! (x rhoa)

        if (rx(k,1) >= rxmin(1)) qx(k,1) = qr(k,1) / rx(k,1)

     else   ! (jnmb(1) == 5)

        ! If prognosing number concentration of cloud (and CCN), prepare CCN
        ! fields and initial and final environmental fields for call to CCNBIN.

        ! There are many options to explore regarding how to apply CCNBIN.
        ! Ideally, the CCNBIN integration should follow a Lagrangian parcel
        ! through its entire CCN nucleation phase.  In the OLAM Eulerian
        ! formulation, Lagrangian parcels must be constructed, and it is an
        ! open question whether to reconstruct them each timestep (as in a
        ! semi-Lagriangian advection method) or whether to carry them through
        ! multiple timesteps that span the entire CCN nucleation episode.  For
        ! the present, we take the practical approach of performing the CCN
        ! integration over time periods equal to the OLAM model long timestep
        ! (DTLONG or DTL), in which case it is to be hoped that a lengthy
        ! nucleation episode will be adequately represented over the time
        ! series of integrations in Eulerian grid cells.  A true Lagrangian
        ! approach remains to be investigated.

        ! Use current conditions in grid cell as "final environmental state"
        ! in ccnbin integration.

        con_cp = cx(k,1) + cx(k,3)
        con_ccny(1:nccntyp) = con_ccnx(k,1:nccntyp)

        tempkB  = tair(k)
        pressB  = press0(k)
        rhoaB   = rhoa(k)
        rhovB   = rhov(k)
        satenvB = rhov(k) / rhovslair(k)

        ! For the present, it will be assumed that "initial environmental
        ! state" conditions are identical with final environmental state
        ! conditions, with the exception of vapor density, which is taken
        ! to be the saturation value.  A more accurate representation of
        ! the initial environmental state would account for parcel movement
        ! and existing microphysical tendencies over the past timestep or,
        ! as discussed above for a Lagrangian parcel, a longer preceding period.

        ! In the future when Lagrangian methods are considered, the following
        ! code may be useful:

        !Lag   ! Compute pseudo vertical velocity based on supersaturation and its
        !Lag   ! assumed production over 1 timestep, using inverse of A1 coefficient
        !Lag   ! in Pruppacher and Klett Eqn. 13-31.  w_nuc is the larger of w_pseudo
        !Lag   ! and wc0(k).  Eventually replace this with Lagrangian parcel
        !Lag   ! supersaturation tendency.

        !Lag   supsat = excessrhov / (1.00001 * rhovslair(k))
        !Lag   supsat_rate = supsat * dtli0
        !Lag   a1inv = (rvap * tair(k)) / (grav * (alvl / (cp * tair(k)) - eps_vapi))
        !Lag   w_pseudo = a1inv * supsat_rate
        !Lag   w_nuc = max(wc0(k),w_pseudo)

        tempkA  = tempkB
        pressA  = pressB
        rhoaA   = rhoaB
        rhovA   = rhovslair(k)
        satenvA = 1.0

        timespan = dtlm(1)
        nbincall = nbincall + 1

        call ccnbin(iw0, k, ntim, timespan, &
                    tempkA, pressA, rhoaA, rhovA, satenvA, &
                    tempkB, pressB, rhoaB, rhovB, satenvB, &
                    con_cp, con_ccny, cactivated, rx_wbc)

        rxnuc = max(rx_wbc, cactivated * emb0(1))

        if (rxnuc > excessrhov) then
           cactivated = cactivated * excessrhov / rxnuc
           rxnuc      = excessrhov
        endif

        cnuc_vc(k) = cactivated  ! [#/m^3]
        rnuc_vc(k) = rxnuc       ! [kg_w/m^3]

        rhov(k) = rhov(k) - rnuc_vc(k)

        cx(k,1) = cx(k,1) + cnuc_vc(k) ! [#/m^3]  x rhoa
        rx(k,1) = rx(k,1) + rnuc_vc(k) ! [kg/m^3] x rhoa
        qr(k,1) = qr(k,1) + rnuc_vc(k) * (tairc(k) * cliq + alli) ! (x rhoa)

        if (rx(k,1) >= rxmin(1)) qx(k,1) = qr(k,1) / rx(k,1)

! Save current environmental variables for next timestep (parcel run only)

        !par    tempkA  = tempkB
        !par    pressA  = pressB
        !par    rhoaA   = rhoaB
        !par    rhovA   = rhovB
        !par    rr_vA   = rr_vB
        !par    satenvA = satenvB

     endif ! (jnmb(1) == 5)

  enddo ! k

end subroutine cldnuc

!===============================================================================

subroutine icenuc(k1,k2,lpw0,mrl0,iw0, &
   rx,cx,qr,qx,emb,vap,tx,rhov,rhoa,press0,dynvisc,thrmcon, &
   tair,tairc,rhovslair,rhovsiair,con_ccnx,con_ifnx,dtl0, &
   rnuc_cp_hom,rnuc_dp_hom,rnuc_vp_haze,rnuc_vp_immers, &
   cnuc_cp_hom,cnuc_dp_hom,cnuc_vp_haze,cnuc_vp_immers)

  use micro_coms,  only: mza0, ncat, dnfac, pwmasi, rxmin, ndnc, ddnc, dtc, &
                         drhhz, dthz, frachz, emb0, emb0i, iccn, igccn, jnmb, &
                         fracc
  use ccnbin_coms, only: nccntyp, ccntyp_alpha
  use consts_coms, only: cice, cliq, alli
  use misc_coms,   only: io6

  implicit none

  integer, intent(inout) :: k1(11)
  integer, intent(inout) :: k2(11)

  integer, intent(in) :: lpw0
  integer, intent(in) :: mrl0
  integer, intent(in) :: iw0

  real, intent(inout) :: rx (mza0,ncat)
  real, intent(inout) :: cx (mza0,ncat)
  real, intent(inout) :: qr (mza0,ncat)
  real, intent(inout) :: qx (mza0,ncat)
  real, intent(in)    :: emb(mza0,ncat)
  real, intent(in)    :: vap(mza0,ncat)
  real, intent(in)    :: tx (mza0,ncat)

  real, intent(inout) :: rhov  (mza0)
  real, intent(in) :: press0   (mza0)
  real, intent(in) :: dynvisc  (mza0)
  real, intent(in) :: thrmcon  (mza0)
  real, intent(in) :: tair     (mza0)
  real, intent(in) :: tairc    (mza0)
  real, intent(in) :: rhovslair(mza0)
  real, intent(in) :: rhovsiair(mza0)
  real, intent(in) :: con_ccnx (mza0,nccntyp)
  real, intent(in) :: con_ifnx (mza0)

  real, intent(inout) :: rnuc_cp_hom   (mza0)
  real, intent(inout) :: rnuc_dp_hom   (mza0)
  real, intent(inout) :: rnuc_vp_haze  (mza0)
  real, intent(inout) :: rnuc_vp_immers(mza0)

  real, intent(inout) :: cnuc_cp_hom   (mza0)
  real, intent(inout) :: cnuc_dp_hom   (mza0)
  real, intent(inout) :: cnuc_vp_haze  (mza0)
  real, intent(inout) :: cnuc_vp_immers(mza0)

  real, intent(in) :: rhoa(mza0)

  real, intent(in) :: dtl0

  integer :: k,idnc,itc,irhhz,ithz

  real :: dn1,dn8,fraccld,ridnc,wdnc2,tc,ritc,wtc2,rhhz
  real :: rirhhz,wrhhz2,thz,rithz,wthz2,frachaz,excessrhov
  real :: con_ifnx_sum

  !---> D14 parameters
  real, parameter :: d14a = 0., d14b = 1., d14c = 0.46, d14d = -11.6
  real, parameter :: cf = 3., cfo1000 = cf / 1000.

  ! FREEZING OF CLOUD DROPLETS: Loop from lowest to highest grid level
  ! that may contain cloud water

  do k = k1(1),k2(1)

     ! Skip current grid level if cloud water is not present

     if (rx(k,1) < rxmin(1)) cycle

     ! HOMOGENEOUS FREEZING OF CLOUD DROPLETS - needs air temperature below -30 C

     if (tairc(k) <= -30.01) then

        dn1 = dnfac(1) * emb(k,1) ** pwmasi(1) ! characteristic diameter
        ridnc = max(1.,min(real(ndnc-1),dn1 / ddnc))
        idnc = int(ridnc)
        wdnc2 = ridnc - real(idnc)

        tc = max(-49.99,tairc(k))
        ritc = (tc + 50.00) / dtc + 1.0
        itc = int(ritc)
        wtc2 = ritc - real(itc)

        fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,mrl0) &
                +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,mrl0) &
                + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,mrl0) &
                +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,mrl0)

        ! cnuc_cp_hom is the number of cloud droplets that are diagnosed to
        ! homogeneously freeze at the given temperature of tairc.  The number
        ! represents a fraction (fraccld) of the number concentration of cloud
        ! droplets [cx(k,1)] present and does not depend on model time step
        ! length.  Repeated application of homogeneous freezing over successive
        ! timesteps would lead to freezing scavenging of (essentially) all cloud
        ! droplets unless the number nucleated were limited in some manner based
        ! on the freezing that already occurred on previous time steps.  The 
        ! number concentration of pristine ice [cx(k,3)] is used as a proxy for
        ! the number of cloud droplets previously frozen in the limiter below.
        ! Saleeby (2009, RAMS) argued that this limit should not be used when
        ! cloud number concentration and CCN are prognosed. 

        if (jnmb(1) == 5 .and. iccn > 0) then
           cnuc_cp_hom(k) = max(0.,fraccld * cx(k,1)) ! x rhoa
        else
           cnuc_cp_hom(k) = max(0.,fraccld * cx(k,1) - cx(k,3)) ! x rhoa
        endif

        rnuc_cp_hom(k) = cnuc_cp_hom(k) * emb(k,1) ! x rhoa

        ! Subtract homogeneous freezing number and mass from cloud water

        cx(k,1) = cx(k,1) - cnuc_cp_hom(k) ! x rhoa
        rx(k,1) = rx(k,1) - rnuc_cp_hom(k) ! x rhoa
        qr(k,1) = qr(k,1) - rnuc_cp_hom(k) * qx(k,1) ! (x rhoa)

        if (rx(k,1) >= rxmin(1)) qx(k,1) = qr(k,1) / rx(k,1)

        ! Add homogeneous freezing number and mass to pristine ice

        cx(k,3) = cx(k,3) + cnuc_cp_hom(k) ! x rhoa
        rx(k,3) = rx(k,3) + rnuc_cp_hom(k) ! x rhoa
        qr(k,3) = qr(k,3) + rnuc_cp_hom(k) * cice * tairc(k) ! (x rhoa)

        if (rx(k,3) >= rxmin(3)) qx(k,3) = qr(k,3) / rx(k,3)

        ! (Changes in qr1 and qr3 have different magnitudes; energy difference will be
        ! accounted for in next theta diagnosis from theta_il and ice/liquid contents.)

        ! if homogeneous freezing has depleted the cloud drops, skip immersion freezing
        if (rx(k,1) < rxmin(1)) cycle

     endif

!--------------- PARCEL EXPTS: turn off hom frz cld
!   cnuc_cp_hom(k) = 0.
!   rnuc_cp_hom(k) = 0.
!---------------------------------------------------

     ! HETEROGENEOUS IMMERSION FREEZING (OF CLOUD DROPLETS) BASED ON DEMOTT (2014)

     !----------------------------------------------------------------------
     ! Code from Gustavo:
     !----------------------------------------------------------------------

     !if (cifnx(k) > 0.) .and. tairc(k) < 0.) then
     !   ! saved as #/cm3 
     !   N_0 = cifnx(k)
  
     !   aux1 = (cf * N_0) ** (-d14a*tairc(k) + d14b)    
     !   aux2 = exp(-d14c*tairc(k) + d14d)
     !   ! D14 result is  STP #/L 
     !   N_a = aux1 * aux2 
     !   N_a = N_a/1000.*(dn0(k)/dn0d10) !ambient per cm3
  
     !   N_a = min(N_0,N_a)

     !   !from #/cc to #/kg
     !   N_a = N_a * 1.0e6/dn0(k)  
     !   if (N_a>1e-4) diagni = N_a
     !endif

     !---------------------------------------------------------------------------
     ! Bob's code version of DeMott 2014 heterogeneous nucleation by immersion
     ! freezing.  Gustavo provided code with parameter values d14a = 0 and
     ! d14b = 1 (although d14b is 1.25 in parts of DeMott 2014, which seems
     ! unjustified physically).  With Gustavo's values, the DeMott formula
     ! simplifies and nucleated number concentration (cnuc_vp_immers) varies
     ! linearly with ambient IFN number concentration (cifnx).  The DeMott
     ! formula expresses ambient IFN number concentration in #/cm^3 and nucleated
     ! number concentration in #/L, but both are #/m^3 in OLAM microphysics, so
     ! a factor of .001 (part of cfo1000) is included here.  Based on discussions
     ! with Bill Cotton, this immersion freezing parameterization requires that
     ! cloud water be present, so this computation is grouped with homogeneous
     ! freezing of cloud water for which the presence of cloud water is already
     ! checked.  Although cloud water is required to initiate immersion freezing,
     ! we assume that nucleated ice particles acquire (the majority of) their
     ! mass from vapor and not from cloud water.  This is consistent with the
     ! minimum ice crystal mass emb0(3) being much larger than the minimum cloud
     ! droplet mass emb0(1).  In addition, we do not deplete cloud droplet number
     ! when immersion freezing occurs.  This is justified from the fact that
     ! concentrations of heterogeneous-nucleated ice crystals are typically
     ! 3 or more orders of magnitude lower than cloud droplet concentrations.
     !---------------------------------------------------------------------------

     if (tairc(k) < 0.) then

        ! Regardless of whether CCN are prognosed or specified, assume that they
        ! can contribute to immersion freezing according to ccntyp_alpha factor

        con_ifnx_sum = con_ifnx(k) &
                     + sum(con_ccnx(k,1:nccntyp) * ccntyp_alpha(1:nccntyp))

        cnuc_vp_immers(k) = con_ifnx_sum &
                          * min(1.,cfo1000 * exp(-d14c*tairc(k) + d14d))

        ! cnuc_vp_immers is the number of IFN that are diagnosed to nucleate to
        ! ice at the given temperature of tairc.  The number represents a
        ! fraction of the IFN present [con_ifnx_sum] and does not depend on model
        ! time step length.  Repeated application of heterogeneous freezing over
        ! successive timesteps would lead to freezing scavenging of (essentially)
        ! all IFN unless the number nucleated were limited in some manner based
        ! on the freezing that already occurred on previous time steps.  The
        ! number concentration of pristine ice [cx(k,3)] is used as a proxy for
        ! the number of IFN previously nucleated in the limiter below.

        cnuc_vp_immers(k) = max(0.,cnuc_vp_immers(k) - cx(k,3)) ! x rhoa
        rnuc_vp_immers(k) = cnuc_vp_immers(k) * emb0(3) ! x rhoa

        ! Prevent heterogeneous nucleation from using up too much available water

        excessrhov = max(rhov(k) - rhovsiair(k), 0.0)

        if (rnuc_vp_immers(k) > excessrhov) then
           cnuc_vp_immers(k) = excessrhov * emb0i(3)
           rnuc_vp_immers(k) = excessrhov
        endif

        ! Add heterogeneous freezing number and mass to pristine ice

        cx(k,3) = cx(k,3) + cnuc_vp_immers(k) ! x rhoa
        rx(k,3) = rx(k,3) + rnuc_vp_immers(k) ! x rhoa
        qr(k,3) = qr(k,3) + rnuc_vp_immers(k) * cice * tairc(k) ! (x rhoa)

        if (rx(k,3) >= rxmin(3)) qx(k,3) = qr(k,3) / rx(k,3)

        ! Subtract heterogeneous freezing mass from water vapor
        ! to conserver water
        rhov(k) = rhov(k) - rnuc_vp_immers(k)

     endif

  enddo

  !  Homogeneous nucleation of haze: Loop over all grid levels

  do k = lpw0,mza0
     rhhz = rhov(k) / rhovslair(k)  ! stays the same

     ! Skip current grid level if r.h. is less than 82% or temp is above -35.01 C

     if (rhhz < 0.82 .or. tairc(k) > -35.01) cycle

     rirhhz = min(0.1799,rhhz - 0.82) / drhhz + 1.0
     irhhz = int(rirhhz)
     wrhhz2 = rirhhz - real(irhhz)

     thz = max(-59.99,tairc(k))
     rithz = (thz + 60.00) / dthz + 1.0
     ithz = int(rithz)
     wthz2 = rithz - real(ithz)

     frachaz = (1.-wrhhz2) * (1.-wthz2) * frachz(irhhz  ,ithz  ) &
             +     wrhhz2  * (1.-wthz2) * frachz(irhhz+1,ithz  ) &
             + (1.-wrhhz2) *     wthz2  * frachz(irhhz  ,ithz+1) &
             +     wrhhz2  *     wthz2  * frachz(irhhz+1,ithz+1)
     frachaz = 1. - exp(-frachaz * dtl0)

     ! Saleeby (2009, RAMS) and Walko (2016) update: Scale the haze nuclei to
     ! the CCN concentration.

     ! cnuc_vp_haz is the number of haze particles that are diagnosed to
     ! homogeneously freeze at the given temperature of tairc.  The number
     ! represents a fraction (frachaz) of the number concentration of haze
     ! particles present and does not depend on model time step length.
     ! Repeated application of homogeneous freezing over successive time steps
     ! would lead to freezing scavenging of (essentially) all haze particles
     ! unless the number nucleated were limited in some manner based on the
     ! freezing that already occurred on previous time steps.  The number
     ! concentration of pristine ice [cx(k,3)] is used as a proxy for the
     ! number of cloud droplets previously frozen in the limiter below.
     ! Saleeby (2009, RAMS) argued that this limit should not be used when
     ! cloud number concentration and CCN are prognosed.

     if (jnmb(1) == 5 .and. iccn > 1) then
        cnuc_vp_haze(k) = max(0.,frachaz * sum(con_ccnx(k,1:nccntyp))) ! x rhoa
     else
        cnuc_vp_haze(k) = max(0.,frachaz * sum(con_ccnx(k,1:nccntyp)) - cx(k,3)) ! x rhoa
     endif

     rnuc_vp_haze(k) = cnuc_vp_haze(k) * emb0(3) ! x rhoa

     ! Prevent haze nucleation from using up too much available water

     excessrhov = max(rhov(k) - rhovsiair(k), 0.0)

     if (rnuc_vp_haze(k) > excessrhov) then
        cnuc_vp_haze(k) = excessrhov * emb0i(3)
        rnuc_vp_haze(k) = excessrhov
     endif

!--------------- PARCEL EXPTS: turn off hom nuc haze
!   cnuc_vp_haze(k) = 0.
!   rnuc_vp_haze(k) = 0.
!---------------------------------------------------

     rx(k,3) = rx(k,3) + rnuc_vp_haze(k)
     qr(k,3) = qr(k,3) + rnuc_vp_haze(k) * cice * tairc(k)
     cx(k,3) = cx(k,3) + cnuc_vp_haze(k)

     if (rx(k,3) >= rxmin(3)) qx(k,3) = qr(k,3) / rx(k,3)

     ! Subtract haze nucleation mass from water vapor to conserve water
     rhov(k) = rhov(k) - rnuc_vp_haze(k)

  enddo

  ! HOMOGENEOUS FREEZING OF DRIZZLE (CLOUD2): Loop from lowest to highest grid level
  ! that may contain drizzle

  do k = k1(8),k2(8)

     ! Skip current grid cell if drizzle is not present or if temperature is above -30.01 C

     if (rx(k,8) < rxmin(8) .or. tairc(k) > -30.01) cycle

     dn8 = dnfac(8) * emb(k,8) ** pwmasi(8) ! characteristic diameter
     ridnc = max(1.,min(real(ndnc-1),dn8 / ddnc))
     idnc = int(ridnc)
     wdnc2 = ridnc - real(idnc)

     tc = max(-49.99,tairc(k))
     ritc = (tc + 50.00) / dtc + 1.0
     itc = int(ritc)
     wtc2 = ritc - real(itc)

     fraccld = (1.-wdnc2) * (1.-wtc2) * fracc(idnc  ,itc  ,mrl0) &
             +     wdnc2  * (1.-wtc2) * fracc(idnc+1,itc  ,mrl0) &
             + (1.-wdnc2) *     wtc2  * fracc(idnc  ,itc+1,mrl0) &
             +     wdnc2  *     wtc2  * fracc(idnc+1,itc+1,mrl0)

     ! cnuc_dp_hom is the number of drizzle droplets that are diagnosed to
     ! homogeneously freeze at the given temperature of tairc.  The number
     ! represents a fraction (fraccld) of the number concentration of drizzle
     ! droplets [cx(k,8)] present and does not depend on model time step
     ! length.  Repeated application of homogeneous freezing over successive
     ! timesteps would lead to freezing scavenging of (essentially) all
     ! drizzle unless the number nucleated were limited in some manner based
     ! on the freezing that already occurred on previous time steps.  The 
     ! number concentration of pristine ice [cx(k,3)] is used as a proxy for
     ! the number of cloud droplets previously frozen in the limiter below.
     ! Saleeby (2009, RAMS) argued that this limit should not be used when
     ! GCCN are prognosed. 

     if (igccn == 2) then
        cnuc_dp_hom(k) = max(0.,fraccld * cx(k,8)) ! x rhoa
     else
        cnuc_dp_hom(k) = max(0.,fraccld * cx(k,8) - cx(k,3)) ! x rhoa
     endif

     rnuc_dp_hom(k) = cnuc_dp_hom(k) * emb(k,8) ! x rhoa

!--------------- PARCEL EXPTS: turn off hom frz drizzle
!   cnuc_dp_hom(k) = 0.
!   rnuc_dp_hom(k) = 0.
!---------------------------------------------------

     ! Subtract homogeneous freezing number and mass from drizzle and add
     ! to pristine ice

     rx(k,3) = rx(k,3) + rnuc_dp_hom(k) ! x rhoa
     cx(k,3) = cx(k,3) + cnuc_dp_hom(k) ! x rhoa
     qr(k,3) = qr(k,3) + rnuc_dp_hom(k) * cice * tairc(k) ! (x rhoa)

     if (rx(k,3) >= rxmin(3)) qx(k,3) = qr(k,3) / rx(k,3)

     rx(k,8) = rx(k,8) - rnuc_dp_hom(k) ! x rhoa
     cx(k,8) = cx(k,8) - cnuc_dp_hom(k) ! x rhoa
     qr(k,8) = qr(k,8) - rnuc_dp_hom(k) * qx(k,8) ! (x rhoa)

     if (rx(k,8) >= rxmin(8)) qx(k,8) = qr(k,8) / rx(k,8)

  enddo

  ! Meyers did ice habit diagnosis here.  Option 1 is to use habit at cloud top.
  ! Option 2 is to use new habit at each level.  Need to consider other options.

end subroutine icenuc

