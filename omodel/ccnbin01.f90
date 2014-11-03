
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Bob notes 9/24/2014: Subroutine ccnbin and related routines below planed for
! "on-line" incorporation into OLAM microphysics for nucleating CCN during a
! model run.  This effort suspended until further considering Lagrangian parcel
! formulation, which is the correct way to consider nucleation and
! environmental forcing rates.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


! Driver to set up arguments for subroutine ccnbin and to call subroutine ccnbin
! for an individual grid point or parcel.  

  program ccnbin_driver

  implicit none

  integer :: nccntyp   ! # of ccn types (e.g. prognosed ccn species in host model)

  integer, allocatable :: nbins_ccntyp(:) ! # of bins to allocate for ccn type

  real, allocatable :: ccntyp_dmin (:) ! ccn type minimum dry diam [cm]
  real, allocatable :: ccntyp_dmax (:) ! ccn type maximum dry diam [cm]
  real, allocatable :: ccntyp_dmed (:) ! ccn type median dry diam [cm]
  real, allocatable :: ccntyp_dsig (:) ! ccn type dry diam shape param []
  real, allocatable :: ccntyp_conc (:) ! ccn type concentration [#/g_a]
  real, allocatable :: ccntyp_kappa(:) ! ccn type kappa param []

  integer :: nkiw, kiw
  integer :: ntim, nccntyp, nbins, ic

  real :: tempk0, tempk1 ! initial & final air temperatures [K]
  real :: press0, press1 ! initial & final pressures [dyn/cm^2]
  real ::  rhoa0,  rhoa1 ! initial & final air densities   [g_a/cm^3]
  real ::  rhov0,  rhov1 ! initial & final vapor densities [g_v/cm^3]
  real ::  sh_v0,  sh_v1 ! initial & final vapor spec dens [g_v/g_a]

  real :: timespan

  real, external :: rhovsl

  nkiw = 1000

  ntim = 500

  timespan = 100.

! Loop over "grid points" or "parcels"

  do kiw = 1,nkiw

     nccntyp = 1

     allocate (nbins_ccntyp(nccntyp))

     allocate (ccntyp_dmin (nccntyp))
     allocate (ccntyp_dmax (nccntyp))
     allocate (ccntyp_dmed (nccntyp))
     allocate (ccntyp_dsig (nccntyp))
     allocate (ccntyp_conc (nccntyp))
     allocate (ccntyp_kappa(nccntyp))

! Fill information for each ccntyp.  Kappa should be (at least slightly)
! greater than 0.0 for particle to be a CCN, and should not exceed 1.28, the
! highest known value for any CCN (NaCl).

!     nbins_ccntyp(1) = 100
!     ccntyp_dmin (1) =  0.001e-4
!     ccntyp_dmax (1) = 60.000e-4
!     ccntyp_dmed (1) =  0.100e-4
!     ccntyp_dsig (1) = 3.162
!     ccntyp_conc (1) = 300.e3
!     ccntyp_kappa(1) = 0.06 

!     nbins_ccntyp(2) = 100
!     ccntyp_dmin (2) =  0.001e-4
!     ccntyp_dmax (2) = 60.000e-4
!     ccntyp_dmed (2) =  0.100e-4
!     ccntyp_dsig (2) = 3.162
!     ccntyp_conc (2) = 300.e3
!     ccntyp_kappa(2) = 0.06 

     nbins_ccntyp(1) = 100
     ccntyp_dmin (1) =  0.001e-4
     ccntyp_dmax (1) = 60.000e-4
     ccntyp_dmed (1) =  0.100e-4
     ccntyp_dsig (1) = 2.0
     ccntyp_kappa(1) = 0.001 + 1.279 * real(kiw-1)/real(nkiw-1) 
     ccntyp_conc (1) = 300.e3

! Count up total number of bins for CCNBIN simulation

     nbins = sum(nbins_ccntyp(1:nccntyp))

! Fill environmental values

     tempk0 = 300.
     press0 = 700.e3
     rhoa0  = 0.8e-3
     rhov0  = 0.999 * (1.0e-3 * rhovsl(tempk0-273.15))
     sh_v0  = rhov0 / rhoa0

     tempk1 = 300.
     press1 = 700.e3
     rhoa1  = 0.8e-3
     rhov1  = 1.11 * (1.0e-3 * rhovsl(tempk1-273.15))
     sh_v1  = rhov1 / rhoa1

     print*, 'kiw ', kiw, ccntyp_kappa(1)
     print*, ' '

! Call CCNBIN subroutine

     if (nbins >= 1) then
        call ccnbin(ntim, nbins, nccntyp, &
                    nbins_ccntyp, ccntyp_dmin, ccntyp_dmax, &
                    ccntyp_dmed, ccntyp_dsig, ccntyp_conc, ccntyp_kappa, &  
                    tempk0, press0, rhoa0, rhov0, sh_v0, &
                    tempk1, press1, rhoa1, rhov1, sh_v1, timespan)
     endif

     deallocate (nbins_ccntyp)

     deallocate (ccntyp_dmin)
     deallocate (ccntyp_dmax)
     deallocate (ccntyp_dmed)
     deallocate (ccntyp_dsig)
     deallocate (ccntyp_conc)
     deallocate (ccntyp_kappa)

  enddo  ! kiw

  end program ccnbin_driver

!=============================================================================

Module ccnbin_coms

  real, parameter :: pi1  = 3.141592654
  real, parameter :: pi2  = pi1 * 2.0
  real, parameter :: pio2 = pi1 * 0.5
  real, parameter :: pio4 = pi1 * 0.25
  real, parameter :: pio6 = pi1 / 6.0

  real, parameter :: cpcgs = 1.0048e7   ! Specific heat of air at const P
  real, parameter :: rdrycgs = 2.8704e6 ! Gas constant of air
  real, parameter :: rvapcgs = 4.615e6  ! Gas constant of vapor
  real, parameter :: rhowcgs = 1.0      ! density of water [g/cm^3]

  real, parameter :: xalphac = 1. ! condens accomodation coeff {was 0.042 in GNF}

! Petters and Kreidenweis (2007) used the following fixed values of temperature
! and surface tension in determining kappa values for various aerosol types.
! Kappa values may need to be re-calibrated if ever these values of tempkap
! and sfctension are changed.

  real, parameter :: tempkap = 298.15
  real, parameter :: sfctension = 72.0

end module ccnbin_coms

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
!     prognose a droplet's diameter.  In addition, diagnosis of unactivated
!     droplet diameter is approximated by linear interpolation between selected
!     droplet wet diameters where equilibrium saturation is first pre-diagnosed
!     iteratively.  Interpolation is sufficiently accurate and is far more
!     efficient than more precise diagnosis arrived at iteratively.

  subroutine ccnbin(ntim, nbins, nccntyp, &
                    nbins_ccntyp, ccntyp_dmin, ccntyp_dmax, &
                    ccntyp_dmed, ccntyp_dsig, ccntyp_conc, ccntyp_kappa, &
                    tempk0, press0, rhoa0, rhov0, sh_v0, &
                    tempk1, press1, rhoa1, rhov1, sh_v1, timespan)

  use ccnbin_coms, only: pi2, pio2, pio4, pio6, rdrycgs, rvapcgs, cpcgs, &
                         rhowcgs, tempkap, sfctension, xalphac

  implicit none

  integer, intent(in) :: ntim        ! # of timesteps to be run in parcel model
  integer, intent(in) :: nbins       ! total # of bins for parcel model
  integer, intent(in) :: nccntyp     ! # of ccn types (chemical species)

  integer, intent(in) :: nbins_ccntyp(nccntyp) ! # of bins for each ccntype

  real, intent(inout) :: ccntyp_dmin (nccntyp) ! ccn type minimum dry diam [cm]
  real, intent(inout) :: ccntyp_dmax (nccntyp) ! ccn type maximum dry diam [cm]
  real, intent(in)    :: ccntyp_dmed (nccntyp) ! ccn type median dry diam [cm]
  real, intent(in)    :: ccntyp_dsig (nccntyp) ! ccn type dry diam shape param []
  real, intent(in)    :: ccntyp_conc (nccntyp) ! ccn type concentration [#/g_a]
  real, intent(in)    :: ccntyp_kappa(nccntyp) ! ccn type kappa param []

  real, intent(in) :: timespan
  real, intent(in) :: tempk0, tempk1 ! initial & final air temperatures [K]
  real, intent(in) :: press0, press1 ! initial & final pressures [dyn/cm^2]
  real, intent(in) ::  rhoa0,  rhoa1 ! initial & final air densities   [g_a/cm^3]
  real, intent(in) ::  rhov0,  rhov1 ! initial & final vapor densities [g_v/cm^3]
  real, intent(in) ::  sh_v0,  sh_v1 ! initial & final vapor spec dens [g_v/g_a]

! Local arrays

  integer :: iprognose (nbins) ! bin prognose flag  (0 = no, 1 = yes)
  integer :: iactivated(nbins) ! bin activated flag (0 = no, 1 = yes)

  real :: drydiam   (nbins)
  real :: drydiam3  (nbins)
  real :: bkappa    (nbins)
  real :: conc      (nbins)
  real :: wetdiam   (nbins)
  real :: wetdiam099(nbins)
  real :: wetdiam100(nbins)
  real :: critdiam  (nbins)
  real :: critsat   (nbins)
  real :: d01sat    (nbins)
  real :: d10sat    (nbins)
  real :: d30sat    (nbins)
  real :: d50sat    (nbins)

! Local variables

  integer :: ibin, ic, jbin, maxsat, itim
  integer :: itim_maxsat, nactivated, inewt

  real :: smax, dtim, tim, fractim, diagdiam, frac
  real :: d01, d10, d30, d50
  real :: tot, xerr1, xerr2, err1, err2

  real :: satenv0, satenv1, satenv, satenvmax
  real :: ax, sqrt2, dsiglog, dmedlog
  real :: tempk, press, rhoa, rhov, sh_v, sh_wb, dsh_wb, rhov_remain
  real :: wetd, dwetd, satk, satkp, satkpp, wetdiam3, dminprog, dmaxinit

  real :: dbnd1, dbnd2, dsig5

  real :: rlvcgs, pratio, rlambda, dvn, rkan, fdcfac, fkafac
  real :: rad, rad2, dsv, fksa, t1, t2, t3, t4, t5, tn, td

  real, external :: erf2, fewcgs, flhvcgs, rhovsl

  dminprog = 3.e-4  ! minimum droplet diameter for prognosis
  dmaxinit = 10.e-4 ! maximum droplet initial diameter (unless large drydiam)

! Initialize bins from input "ccntyp" arrays

  sqrt2 = sqrt(2.)
  ibin = 0
  do ic = 1,nccntyp

     ax = log(ccntyp_dmax(ic) / ccntyp_dmin(ic)) / (nbins_ccntyp(ic) * log(2.))
     dbnd1 = ccntyp_dmin(ic)
     tot = 0.
 
     do jbin = 1,nbins_ccntyp(ic)
        ibin = ibin + 1

! Set up the CCN size spectrum -- geometric progression of diameter

        dbnd2 = ccntyp_dmin(ic) * 2.0**(jbin * ax)
        drydiam (ibin) = 0.5 * (dbnd1 + dbnd2)
        drydiam3(ibin) = drydiam(ibin)**3
        bkappa  (ibin) = ccntyp_kappa(ic)
        wetdiam (ibin) = drydiam(ibin)

! Set up initial aerosol distribution

        dsiglog = log(ccntyp_dsig(ic))
        dmedlog = log(ccntyp_dmed(ic))

        xerr1 = (log(dbnd1) - dmedlog) / (sqrt2 * dsiglog)
        xerr2 = (log(dbnd2) - dmedlog) / (sqrt2 * dsiglog)
        err1  = erf2(xerr1)
        err2  = erf2(xerr2)
        conc(ibin) = abs((err2 - err1) * ccntyp_conc(ic) * .5) ! [#/g_a]
        tot = tot + conc(ibin)

        dbnd1 = dbnd2
     enddo

  enddo

! Initial and final environmental saturation values (actual final value will
! be less than satenv1 due to condensation onto CCN)

  satenv0 = rhov0 / (1.0e-3 * rhovsl(tempk0-273.15))
  satenv1 = rhov1 / (1.0e-3 * rhovsl(tempk1-273.15))

! Loop over bins and for each one diagnose selected points on the Kohler curve
! using calls to subroutine satkap

  do ibin = 1,nbins

! Diagnose critical diameter and saturation using Newton-Raphson method
! (start wetd greater than drydiam for faster convergence; factor of
! 1.01 is ok for bkappa >= 0.001 and drydiam >= 0.001 micron.)

     wetd = drydiam(ibin) * 1.01
     do inewt = 1,50
        call satkap(wetd, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 1)

        if (satkpp < 0.) then
           dwetd = - satkp / satkpp
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

     wetd = drydiam(ibin) * 1.0001
     do inewt = 1,50
        call satkap(wetd, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 0)

        if (satkp > 0.) then
           dwetd = - (satk - 0.99) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp ',ibin,inewt,drydiam,wetd,satk,satkp
           stop
        endif
     enddo
     wetdiam099(ibin) = wetd

! Diagnose equilibrium diameter at satenv = 1.00 using Newton-Raphson method

     wetd = drydiam(ibin) * 1.0001
     do inewt = 1,50
        call satkap(wetd, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 0)

        if (satkp > 0.) then
           dwetd = - (satk - 1.00) / satkp
           if (abs(dwetd/wetd) < 1.e-6) exit
           wetd = wetd + dwetd
        else
           print*, 'non-positive satkp ',ibin,inewt,drydiam,wetd,satk,satkp
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

     call satkap(d01, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 0)
     d01sat(ibin) = satk

     call satkap(d10, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 0)
     d10sat(ibin) = satk

     call satkap(d30, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 0)
     d30sat(ibin) = satk

     call satkap(d50, drydiam3(ibin), bkappa(ibin), satk, satkp, satkpp, 0)
     d50sat(ibin) = satk

  enddo

! If critical saturation of all bins exceeds ending environmental saturation, return

  if (all(critsat(1:nbins) > satenv1)) return

! Initialize quantities before time loop

  satenvmax = satenv0
  rlvcgs = FLHVcgs(tempk0)
  maxsat = 0
  iprognose(:) = 0
  iactivated(:) = 0

  dtim = timespan / ntim
  tim = 0.

! Start the time loop

  do itim = 1,ntim 

     tim = tim + dtim
     fractim = tim / timespan

! Sum water mixing ratio over all bins [g_wb/g_a]

     sh_wb = rhowcgs * pio6 &
           * sum(conc(1:nbins) * (wetdiam(1:nbins)**3 - drydiam3(1:nbins)))

! Update environmental variables

     tempk = tempk0 + fractim * (tempk1 - tempk0) + rlvcgs * sh_wb / cpcgs
     press = press0 + fractim * (press1 - press0)
     rhoa  = rhoa0  + fractim * ( rhoa1 -  rhoa0)
     rhov  = rhov0  + fractim * ( rhov1 -  rhov0)
     sh_v  = sh_v0  + fractim * ( sh_v1 -  sh_v0)

! Vapor density remaining after bin water uptake [g_v/cm^3]

     rhov_remain = rhov * (1. - sh_wb / sh_v)

! Saturation ratio; Latent heat of vaporization

     satenv = rhov_remain / (1.0e-3 * rhovsl(tempk-273.15))
     rlvcgs = FLHVcgs(tempk)

! Atmospheric properties used for computing droplet vapor diffusional growth

     pratio  = 1.01325e6 / press
     rlambda = 6.6e-6 * pratio * (tempk / 293.15)   ! mean free path P&K 10-106
     dvn     = 0.211 * pratio * (tempk / 273.15)**1.94   ! diffusivity P&K 13-3
     rkan    = 4326. + 7.118 * tempk   ! therm cond P&K 13-16 but in cgs with K
     fdcfac  = dvn * SQRT(pi2 / (rvapcgs * tempk)) / xalphac
     fkafac  = rkan * SQRT(pi2 / (rdrycgs * tempk)) / (0.96 * cpcgs * rhoa)

! Loop over bins to update droplet diameters

     do ibin = 1,nbins

! For bins that have not yet been prognosed, diagnose diameter by interpolation
! between pre-determined diameters

        if (iprognose(ibin) == 0) then

           if (satenv >= critsat(ibin)) then
              diagdiam = max(dminprog,critdiam(ibin))
           elseif (satenv >= d50sat(ibin)) then
              frac = (satenv - d50sat(ibin)) / (critsat(ibin) - d50sat(ibin))
              diagdiam = (.50 + .50 * frac) * critdiam(ibin) &
                       + (.50 - .50 * frac) * wetdiam100(ibin)
           elseif (satenv >= d30sat(ibin)) then
              frac = (satenv - d30sat(ibin)) / (d50sat(ibin) - d30sat(ibin))
              diagdiam = (.30 + .20 * frac) * critdiam(ibin) &
                       + (.70 - .20 * frac) * wetdiam100(ibin)
           elseif (satenv >= d10sat(ibin)) then
              frac = (satenv - d10sat(ibin)) / (d30sat(ibin) - d10sat(ibin))
              diagdiam = (.10 + .20 * frac) * critdiam(ibin) &
                       + (.90 - .20 * frac) * wetdiam100(ibin)
           elseif (satenv >= d01sat(ibin)) then
              frac = (satenv - d01sat(ibin)) / (d10sat(ibin) - d01sat(ibin))
              diagdiam = (.01 + .09 * frac) * critdiam(ibin) &
                       + (.99 - .09 * frac) * wetdiam100(ibin)
           elseif (satenv >= 1.00) then
              frac = (satenv - 1.00) / (d01sat(ibin) - 1.00)
              diagdiam = (       .01 * frac) * critdiam(ibin) &
                       + (1.00 - .01 * frac) * wetdiam100(ibin)
           elseif (satenv >= 0.99) then
              frac = (satenv - 0.99) / 0.01
              diagdiam = (       frac) * wetdiam100(ibin) &
                       + (1.00 - frac) * wetdiam099(ibin)
           else
              diagdiam = wetdiam099(ibin)
           endif

           if (satenv >= 0.99) then
              if (diagdiam >= dminprog * 0.999 .or. &
                  drydiam(ibin)*1.1 >= dminprog) then

                 wetdiam(ibin) = max(drydiam(ibin)*1.1, min(dmaxinit,diagdiam))
                 iprognose(ibin) = 1
              endif
           endif

        endif  ! iprognose(ibin) == 0

! Prognose wet diameter for bins flagged for prognosis

        if (iprognose(ibin) == 1) then

! Time derivative of droplet diameter due to vapor diffusion.  The quantity
! (t5 * exp(t4)) below is the rhs of Eq. 6 of Petters and Kreidenweis (2007).

           rad  = wetdiam(ibin) * 0.5
           rad2 = rad * rad / (rlambda + rad)
           dsv  = dvn * rad / (rad2 + fdcfac) ! modified diffusivity P&K 13-13
           fksa = rkan * rad / (rad2 + fkafac) ! modified therm cond P&K 13-20

           t1 = rlvcgs * rhowcgs / (tempk * fksa) ! P&K 13-28 denom term 2A
           t2 = rlvcgs / (tempk * rvapcgs) - 1.   ! P&K 13-28 denom term 2B
           t3 = rhowcgs * tempk * rvapcgs / (dsv * FEWcgs(tempk)) ! P&K 13-28 denom term 1

           call satkap(wetdiam(ibin), drydiam3(ibin), bkappa(ibin), &
                       satk, satkp, satkpp, 0)

           tn = satenv - satk
           td = t3 + (t1 * t2)

! Update of wet diameter squared.  To improve computational stability, use
! volume tendency of growing drop and diameter tendency of shrinking drop.

           if (tn > 0.) then
              wetdiam(ibin) = (wetdiam(ibin)**3 &
                            + dtim * 12. * wetdiam(ibin) * tn / td)**(1./3.)
           else
              wetdiam(ibin) = max(drydiam(ibin) * 1.0001, wetdiam(ibin) &
                            + dtim * 4. * tn / (wetdiam(ibin) * td))
           endif

        endif  ! iprognose(ibin) == 1

     enddo  ! ibin

     satenvmax = max(satenvmax,satenv)

! If running ccnbin as part of OLAM (or other model), halt integration of
! ccnbin once satenv has begun to decrease in time (due to newly activated
! droplets consuming excess vapor faster than environmental production).  This
! prevents newly activated droplets from consuming excess moisture (aided by
! their solute effect) before they are introduced to the cloud droplet
! population which has no solute effect. 

     if (satenv < satenvmax - 0.00001) exit

  enddo    ! itim

! Identify "activated" particles as those whose critical saturation was
! exceeded by environmental saturation during the time integration.  Even
! though very large CCN are often not truly activated (because there was
! insufficient time for them to grow to their critical diameters), they
! are still more amenable to continued growth than smaller truly-activated
! CCN and thus qualify to become cloud droplets.

  nactivated = 0

  do ibin = 1,nbins
     if (critsat(ibin) < satenvmax) then
        iactivated(ibin) = 1
        nactivated = nactivated + 1
     endif
  enddo

  end subroutine ccnbin

!=============================================================================

      real function FEWcgs (t)

! Function FEWcgs calculates saturation vapor pressure [dyn/cm^2] over water for
! a given temperature [K].

      data ts,sr /373.16, 3.0057166/ !log(1013.246)=3.0057166
      ar  = ts/t
      ar1 = ar - 1.
      br  = 7.90298 * ar1
      cr  = 5.02808 * LOG10(ar)
      dw  = 1.3816e-7 * (10.**(11.344 * (ar1/ar)) - 1.)
      er  = 8.1328e-3 * (10.**(-3.49149 * ar1) - 1.)
      answer = 10.**(cr-dw+er+sr-br)
      FEWcgs = answer * 1000.

      end function fewcgs

!===============================================================================

real function rhovsl(tc)

! Function rhovsl calculates the density of water vapor at saturation
! over liquid as a function of Celsius temperature

implicit none

real, intent(in) :: tc

real, parameter :: c0 = .6105851e+03 ,c1 = .4440316e+02 ,c2 =  .1430341e+01
real, parameter :: c3 = .2641412e-01 ,c4 = .2995057e-03 ,c5 =  .2031998e-05
real, parameter :: c6 = .6936113e-08 ,c7 = .2564861e-11 ,c8 = -.3704404e-13
real, parameter :: rvap = 461.
real :: esl,x

x = max(-80.,tc)
esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rhovsl = esl / (rvap * (tc + 273.15))

end function rhovsl

!=============================================================================

      real function FLHVcgs (t)

!  purpose:  This function calculates the latent heat of vaporization for     
!    water as per P&K eqn. 4-85a. The result is then converted to     
!           ergs/g.  T is deg K.                                             

      rga  = 0.167 + (3.67e-4 * t)
      rlhv  = 597.3 * (273.15 / t)**rga
      FLHVcgs = rlhv / 2.38844e-08

!Bob:  A reasonable quadratic fit (in J/kg) is:
!Bob:  flhv = 2500795. + tempc * (2.3 * tempc - 2453.)
!Bob:  (Result is within 300 J/kg of above formula between +/- 40 deg C)

      end function flhvcgs

!=============================================================================

      real FUNCTION erf2(x)
      REAL x, gammp
!     USES gammp
!      Returns the error function erf(x).
      if(x.lt.0.)then
         erf2=-gammp(.5,x**2)
      else
         erf2=gammp(.5,x**2)
      endif
      END function erf2

!=============================================================================

      real FUNCTION gammp(a,x)
      REAL a,x
!     USES gcf,gser
!      Returns the incomplete gamma function P(a, x).
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'fbad arguments in gammp'
      if(x.lt.a+1.)then         !Use the series representation.
         call gser(gamser,a,x,gln)
      gammp=gamser
      else                      !Use the continued fraction representation
         call gcf(gammcf,a,x,gln)
         gammp=1.-gammcf        !and take its complement.
      endif
      END function gammp

!=============================================================================

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
!     USES gammln
!     Returns the incomplete gamma function P(a, x) evaluated by its 
!     series representation as gamser. Also returns ln(a) as gln.
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
         if(x.lt.0.)pause 'x < 0 in gser'
         gamser=0.
         return
      endif
      ap=a
      sum=1./a
      del=sum
      do n=1,ITMAX
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if(abs(del).lt.abs(sum)*EPS)goto 1
      enddo 
      pause 'a too large, ITMAX too small in gser'
 1    gamser=sum*exp(-x+a*log(x)-gln)      
      END subroutine gser

!=============================================================================

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!     USES gammln
!     Returns the incomplete gamma function Q(a, x) evaluated by its 
!     continued fraction representation as gammcf. Also returns ln (a) 
!     as gln.
!     Parameters: ITMAX is the maximum allowed number of iterations; 
!     EPS is the relative accuracy; FPMIN is a number near the smallest 
!     representable floating-point number.
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a                  !Set up for evaluating continued 
                                !fraction by modified
      c=1./FPMIN                !Lentzs method (5.2) with b0 = 0.
      d=1./b
      h=d

      do i=1,ITMAX              !Iterate to convergence.
         an=-i*(i-a)
         b=b+2.
         d=an*d+b
         if(abs(d).lt.FPMIN)d=FPMIN
         c=b+an/c
         if(abs(c).lt.FPMIN)c=FPMIN
         d=1./d
         del=d*c
         h=h*del
         if(abs(del-1.).lt.EPS)goto 1
      enddo 
      pause 'a too large, ITMAX too small in gcf'
 1    gammcf=exp(-x+a*log(x)-gln)*h !Put factors in front.
      END subroutine gcf

!=============================================================================

      DOUBLE PRECISION FUNCTION gammln(xx)
!     Returns the value of Lanczos' approximation to ln(gamma(xx)), for n > 0.
!     Taken from Numerical Recipes in Fortran, 2nd ed., Cambridge University Press.
!
      INTEGER j
      REAL xx
      DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
      SAVE cof, stp
      DATA cof, stp/76.18009172947146d0, -86.50532032941677d0, &
      24.01409824083091d0, -1.231739572450155d0, &
      0.1208650973866179d-2, -.5395239384953d-5, &
      2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 10 j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y
 10   continue
      gammln=tmp+dlog(stp*ser/x)
      END function gammln

!=============================================================================

! Subroutine satkap evaluates the saturation vapor mixing ratio at the surface
! of a wetted hygroscopic aerosol [based on Eqn (6) of Petters and Kreidenweis
! (2007)] and its first and second derivatives with respect to wet diameter.  

  subroutine satkap(wetdiam, drydiam3, bkappa, satk, satkp, satkpp, ids)

  use ccnbin_coms, only: rvapcgs, rhowcgs, tempkap, sfctension

  implicit none

  real, intent(in)  :: wetdiam, drydiam3, bkappa
  real, intent(out) :: satk, satkp, satkpp

  integer, intent(in) :: ids

  real(kind=8) :: wetdiam2, wetdiam3, curvcof, curvterm
  real(kind=8) :: actdenom, actterm, term1, term2, term3, term4, term5

  wetdiam2 = wetdiam * wetdiam
  wetdiam3 = wetdiam2 * wetdiam

  curvcof  = 4.0 * sfctension / (rvapcgs * tempkap * rhowcgs)
  curvterm = exp(curvcof / wetdiam)
  actdenom = wetdiam3 - drydiam3 * (1 - bkappa)
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

  end subroutine satkap
