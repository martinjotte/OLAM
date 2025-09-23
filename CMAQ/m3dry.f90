
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
! $Header: /project/yoj/arc/CCTM/src/depv/m3dry/m3dry.F,v 1.12 2012/01/19 14:19:43 yoj Exp $

SUBROUTINE m3dry ( iw, abflux, sfc_hono, nsfc, depvel_gas )

!-------------------------------------------------------------------------------
! Name:     Models-3 Dry Deposition
! Purpose:  Computes dry deposition velocities using Rst and Ra, and
!           elements of ADOM DD model.
! Revised:  21 Jan 1998  Original version.  (J. Pleim and A. Bourgeois)
!           18 Sep 2001  Made general for USGS 24-category system.
!                        (T. Otte, J. Pleim, and W. Hutzell)
!           14 Jan 2002  Added temperature dependence to Henry's Law
!                        constants.  Added temperature and pressure
!                        dependence to diffusivity.  Added new dry
!                        deposition species, methanol.  (Y. Wu and T. Otte)
!           18 Jan 2002  Changed the reference wet cuticle resistance.
!                        (J. Pleim)
!           09 Jun 2003  Added logic for modeling snow covered surfaces.
!                        Changed the reactivities for SO2, HNO3 and NH3.
!                        Changed pH values to have an east-west variation.
!                        Using the Henry's law constant function from CMAQ in
!                        place of local code.  Also changed the code for
!                        deposition to water to use a pH of 8.1 and the
!                        temperature of water in calculating the Henry's law
!                        constant.  Adjusted values of RSNOW0 = 1000 and
!                        A(NH3) = 20.  Added new dry deposition species: N2O5,
!                        NO3, Generic_aldehyde.  Corrected diffusivities of
!                        chemicals and water and viscosity of air to all be at
!                        the same temperature (273.15K).  Temperature and
!                        pressure adjustments to the values are not needed
!                        because the diffusivities and viscosity are always used
!                        as ratios, so the temperature-pressure dependence was
!                        removed.  Removed dry deposition species, ATRA and
!                        ATRAP, from output.  (D. Schwede, J. Pleim, and
!                        T. Otte)
!           28 Feb 2005  Added optional dry deposition species for chlorine
!                        and mercury.  (G. Sarwar, R. Bullock, and T. Otte)
!           02 Feb 2006  Added mesophyll resistance to dry deposition velocity
!                        calculation, and defined non-zero value for mercury.
!                        (D. Schwede, J. Pleim, and R. Bullock)
!           01 Aug 2007  Added a non-zero mesophyll resistance for NO, NO2, and
!                        CO.  Restored wet cuticle resistance for O3 based on
!                        field study measurements.  Added wet ground resistance.
!                        Changed ground resistance to include partitioning of
!                        wet and dry ground.  Updated pH of rain water for
!                        eastern United States and outside of North America.
!                        Changed reactivity for PAN.  Removed dry deposition
!                        velocity calculations for obsolete chlorine species
!                        ICL1 and ICL2.  Corrected error in the calculation of
!                        surface resistance over water where (Sc/Pr)**(2/3) had
!                        been inadvertently omitted from the numerator.
!                        Surface resistance over water is now a function of
!                        species.  Surface resistance over water now uses wet
!                        bulb temperature rather than ground (water) temperature
!                        in the calculation of the effective Henry's law
!                        constant, and the algorithm has been updated.  Changed
!                        (Sc/Pr)**(2/3) over water to a species-dependent,
!                        meteorologically dependent variable.  Effective Henry's
!                        law constant over land now uses 2-m temperature rather
!                        than layer 1 temperature. Changed ES
!                        into ES_AIR and ES_GRND, and changed QSS into QSS_AIR
!                        and QSS_GRND to clarify usage.  (J. Pleim, E. Cooter,
!                        J. Bash, T. Otte, and G. Sarwar)
!           07 Dec 2007  Add into CMAQ for in-line deposition velocities.
!                        (W. Hutzell, J. Young and T. Otte)
!           07 Jan 2008  Changed the value of d3, the scaling parameter used to
!                        estimate the friction velocity in surface waters from
!                        the atmospheric friction velocity to a value following
!                        Slinn et al. (1978) and Fairall et al. (2007).
!                        (J. Bash)
!           01 Feb 2008  Added bidirectional NH3 flux calculations. (J. Pleim and
!                        J. Young)
!           20 Mar 2008  Added a trap for undefined dry deposition velocities
!                        (e.g., NaN's).  (T. Otte)
!           21 Mar 2008  Added heterogeneous reaction for HONO. It affects HONO, NO2
!                        and HNO3 (G. Sarwar)
!           30 Apr 2008  Added five air toxic species to output.  (W. Hutzell
!                        and T. Otte)
!           August 2008  Applied a minimum value to ustg (0.001 m/s) to prevent
!                        negative rbg values in the bidi calculation (J. Pleim)
!           05 Oct 2009  Added condition that vegetation fraction must be
!                        greater than zero to be considered a land point.  This
!                        works around intermittent inconsistencies in surface
!                        fields in some WRF data sets.  (T. Otte)
!           Dec 2009     Revised bidirectional NH3 flux calculations to use soil
!                        Gamma values read from gridded file.  Bidi flux calcs
!                        based on comparisons with Lillington, NC corn data 2007 (J. Pleim)
!           30 Mar 2010  Modified to output the NH3 bidi stomtal, cuticle and soil component
!                        fluxes and chaged NH3 bidi variables used in estimating the compensation
!                        point double. (J. Bash)
!           16 Feb 2011 S.Roselle: replaced I/O API include files with UTILIO_DEFN
!           20 May 2011  D.Schwede: add MOSAIC processing
!           14 Jul 2011  Replaced dw25 calculation with Hayduk and Laudie method.
!                        LeBas molar volumes are from the Schroeder additive method
!                        with the exception of HGIIGAS (modeled as HgCl2) which was
!                        obtained using the Tyn and Calus method. Also, ICL1 and ICL2
!                        were removed. (D. Schwede)
!           27 Jul 2011 J.Bash: Parmaterized the mesophyll resistance as a function
!                        of solubility following Wesely 1989 Atmos Environ.
!           15 Aug 2011  Modified HONO calculation so that deposition velocity for NO2
!                        that is output in DEPV file does not include the loss due to
!                        the heterogeneous reaction. This additional loss is now
!                        accounted for in vdiff.F (D. Schwede and G. Sarwar)
!           29 Aug 2011 Added NH3 bidirectional flux variables and modules and integrated
!                       NH3 bidi algorithms with MOSAIC algorithms. NH3 bidi routines now
!                       read in foratted EPIC output and maintain a soil NH4 budget and
!                       fluxes are calculated for individual land cover types.
!                       (J. Bash and D. Schwede)
!           22 Sep 2011 -- incorporated twoway model implementation
!                       -- removed non-use dluse array
!                          (David Wong)
!           26 Sep 2011 -- made the number of actual and dummy arguments the same in
!                          calling subroutine Init_ABFlux
!                          (David Wong)
!-------------------------------------------------------------------------------

  use depvvars
  use depv_defn,   only: ie_hono, ic_no2, n_unique_gdepv
  use utilio_defn
  use const_data
  use cgrid_defn,  only: cgrid, vdemis_gc
  use mem_ijtabs,  only: itab_w
  use mem_sfcg,    only: itab_wsfc, sfcg
  use mem_land,    only: land, omland, nzg
  use mem_sea,     only: sea, omsea
  use mem_grid,    only: dzt, volti, lpw
  use mem_basic,   only: theta, rho
  use mem_turb,    only: pblh
  use leaf_coms,   only: tai_max, veg_frac, wcap_min
  use hlconst_mod, only: hlconst
  use micro_coms,  only: ccnparm
  use mem_micro,   only: cldnum
  use therm_lib,   only: rhovsl_inv, qtk, qtk_sea, qwtk
  use mem_megan,   only: pfts
  use oname_coms,  only: nl

  implicit none

! Arguments:

  integer, intent( in  ) :: iw
  logical, intent( in  ) :: abflux
  logical, intent( in  ) :: sfc_hono
  integer, intent( in  ) :: nsfc
  real,    intent( out ) :: depvel_gas(n_unique_gdepv, nsfc)  ! deposition  velocity [ m/s ]

! Local Variables:

  REAL,            PARAMETER :: d3         = 1.38564e-2 ! [dim'less]
                                               ! k*sqrt(rhoair/rhowater) from Slinn 78
  REAL                       :: dw
  LOGICAL                    :: effective  ! true=compute effective Henry's Law const
  REAL                       :: heff                 ! effective Henry's Law constant
  REAL                       :: heff_ap              ! Henry's Law constant for leaf apoplast M/atm
  REAL                       :: hplus
  REAL,            PARAMETER :: hplus_ap    = 1.0e-6     ! pH=6.0 leaf apoplast solution Ph (Massad et al 2008)
  REAL,            PARAMETER :: hplus_dirty = 1.0e-5     ! pH=5.0
  REAL,            PARAMETER :: hplus_clean = 3.16228e-6 ! 10.0**(-5.5)
  REAL,            PARAMETER :: hplus_h2o   = 7.94328e-9 ! 10.0**(-8.1)
  REAL,            PARAMETER :: hplus_def   = 1.0e-5     ! pH=5.0
  real                       :: hveg
  LOGICAL                    :: isurban
  LOGICAL                    :: isocean
  integer                    :: iland, isea, jsfc, iwsfc, jasfc
  REAL                       :: kviswi      ! 1 / kinematic viscosity of water [s/cm^2]
  INTEGER                    :: kw
  INTEGER                    :: ks
  REAL                       :: arcoarkw
  INTEGER                    :: l
  REAL                       :: lai, laiv          ! leaf area index
  REAL                       :: tai, taiv          ! total area index
  INTEGER                    :: n
  CHARACTER( 16 ), PARAMETER :: pname      = 'M3DRY'
  real                       :: pvd
  real                       :: rai                  ! cell aerodynamic resistance
  real                       :: raicei               ! aerodynamic resistance over ice
  real                       :: rawi                 ! aerodynamic resistance over water
  REAL                       :: rbci
  REAL                       :: rbsulfi
  REAL                       :: rci
  REAL,            PARAMETER :: rcw0       = 125000.0 ! adapted from Slinn 78
  REAL,            PARAMETER :: rcw0i      = 1.0 / rcw0

  REAL                       :: rgndi
  REAL                       :: rgwi                 ! resist for water-covered sfc

  real,            parameter :: cwetO3 = 1.0 / 1250.0 ! O3 wet cuticle conductance, Altimir et al 2006

  REAL,            PARAMETER :: rgwet0     = 25000.0      ! [s/m]
  REAL,            PARAMETER :: rgwet0i    = 1.0 / rgwet0 ! m/s]

  REAL                       :: rh_air               ! rel humidity (air)
  REAL                       :: rinci                ! conductance through canopy to ground
  REAL                       :: rsi                  ! stomatal resistance [s/m]
  REAL,            PARAMETER :: rsndiff    = 10.0    ! snow diffusivity fac
  REAL                       :: rsnowi
  REAL                       :: rstomi
  REAL                       :: rsurfi
  REAL                       :: rweti                ! wet sfc resist (cuticle or grnd)
  REAL                       :: rwetsfci
  REAL                       :: scw_pr_23            ! (scw/pr)**2/3
  real                       :: sst
  REAL,            PARAMETER :: svp2       = 17.67   ! from MM5 and WRF
  REAL,            PARAMETER :: svp3       = 29.65   ! from MM5 and WRF
  REAL,            PARAMETER :: rt25inK    = 1.0/(stdtemp + 25.0) ! 298.15K = 25C
  REAL                       :: tempcr               ! cell canopy temp
  REAL                       :: tempgr               ! cell ground temp
  REAL                       :: tice
  REAL                       :: tw
  REAL,            PARAMETER :: onethird   = 1.0 / 3.0
  REAL,            PARAMETER :: twothird   = 2.0 / 3.0
  REAL                       :: ustar                ! cell friction velocity
  REAL                       :: ustari               ! friction velocity over ice
  REAL                       :: ustarw               ! friction velocity over water
  REAL                       :: vegfr                ! cell veg coverage fraction
  real                       :: wrcr                 ! vegetation water (kg/m^2)
  REAL                       :: wrmax
  REAL                       :: wstar
  REAL                       :: wtv0
  REAL                       :: xm                   ! liquid water mass frac
  REAL                       :: xt                   ! liquid water mass frac
  real                       :: rh_func              ! RH function for the development of a water
                                                         !      film on leaf cuticles
  real                       :: wsat_vg
  real                       :: soil_water1
  logical                    :: hasveg
  integer                    :: lc
  real                       :: taimax
  real                       :: laimax

  real :: tempvg
  real :: deltag, liqfrg
  real :: deltav, liqfrv
  real :: fact, sfacm
  real :: rgliqi
  real :: rdryi
  real :: rcuti
  real :: rgdi
  real :: rgndci
  real :: rmesoi

!-------------------------------------------------------------------------------
! For ozone exchange over the ocean
  REAL                       :: pchang      ! p value used in equation (12) in Chang et al( 2004)
  REAL                       :: kwchang     ! kw value used in Chang et al (2004)
  REAL                       :: ciodide     ! iodide concentration (umol/L)
  REAL                       :: qiodide     ! q in Chang et al (2004)

!-------------------------------------------------------------------------------
! For heterogenous hono on surfaces
  REAL                       :: kno2        ! first order rate constant for the heterogenous reaction [1/s]
  REAL                       :: surf_bldg   ! Surface area of bldgs to volume of air [-]
  REAL                       :: surf_leaf   ! Surface area of leaves to volume of air [-]
  REAL                       :: zf          ! layer height [m]

!-------------------------------------------------------------------------------
! For Ammonia bi-directional flux
! REAL                       :: sltyp     ! soil type
! REAL                       :: wg        ! soil moisture in top 1 cm (vol frc)
! REAL                       :: w2        ! soil moisture in top 1 m (vol frc)
! REAL                       :: cnh3      ! Layer 1 NH3 concentration from CGRID [ppm]
! REAL                       :: tpvd      ! temp pvd variable
! REAL                       :: f_stom    ! stomatal flux component
! REAL                       :: f_cut     ! cuticular flux component
! REAL                       :: f_soil    ! soil flux component
! REAL                       :: f_emis    ! Emission flux component
! REAL                       :: f_dep     ! Deposition flux component
! REAL                       :: f_ag      ! flux from agriculture
! REAL                       :: f_nat     ! flux from natural systems
! REAL                       :: f_wat     ! direct flux to surface water
! REAL                       :: lnh3      ! loss part of ammonia flux
                                          ! Note that lnh3 is stored in DEPVEL_GAS(nh3) for use in vdiff

!-------------------------------------------------------------------------------
! Loop over land cells and calculate dry deposition.
!-------------------------------------------------------------------------------

  effective = .TRUE.

  depvel_gas = 0.0

  do jsfc = itab_w(iw)%jland1, itab_w(iw)%jland2
     iwsfc = itab_w(iw)%iwsfc(jsfc)
     jasfc = itab_w(iw)%jasfc(jsfc)

     kw = itab_wsfc(iwsfc)%kwatm(jasfc)

     iland = iwsfc - omland

     ks = kw - lpw(iw) + 1

     arcoarkw = itab_wsfc(iwsfc)%arcoarkw(jasfc)

     ustar = sfcg%ustar(iwsfc)
     zf    = dzt(kw)
     rai   = sfcg%ggaer(iwsfc)
     wtv0  = sfcg%wthv(iwsfc)
     lc    = sfcg%leaf_class(iwsfc)

     if (tai_max(lc) < 0.1) then
        hveg   = 0.0
        vegfr  = 0.0
        lai    = 0.0
        tai    = 0.0
        tempvg = land%veg_temp(iland)
        hasveg = .false.
     else
        hveg   = land%veg_height(iland)
        vegfr  = max(1.0 - pfts(0,iland), 0.1)
        lai    = land%veg_lai(iland)
        tai    = land%veg_tai(iland)
        tempvg = land%veg_temp(iland)
        hasveg = .true.

        taimax = tai_max(lc) / veg_frac(lc)
        laimax = taimax * lai / tai

        laiv = min(lai / vegfr, laimax)
        taiv = max(tai / vegfr, taimax)
     endif

     tempcr = sfcg%cantemp(iwsfc)
     rh_air = 100.0 * sfcg%canrrv(iwsfc) * real(rho(kw,iw)) * rhovsl_inv( sfcg%cantemp(iwsfc) - 273.15 )
     rh_air = min( 100.0, max( rh_air, 0.0 ) )

     arcoarkw = itab_wsfc(iwsfc)%arcoarkw(jasfc)
     isurban = ( lc == 19 )

     if (wtv0 > 0.0) then
        wstar  = (grav * pblh(iw) * wtv0 / theta(kw,iw)) ** onethird
     else
        wstar  = 0.0
     endif

     wsat_vg = land%wsat_vg(nzg,iland)
     soil_water1 = min(land%soil_water(nzg,iland), wsat_vg)

     if (land%sfcwater_mass(1,iland) >= wcap_min) then
        if (land%sfcwater_mass(2,iland) >= wcap_min) then
           ! Surface wetness factor
           deltag = min( 1.0, (land%sfcwater_mass(2,iland) / 0.2)**2 ) ** onethird
           ! Sfcwater temperature and liquid water mass fraction
           call qtk( land%sfcwater_energy(2,iland), tempgr, liqfrg )
        else
           deltag = min( 1.0, (land%sfcwater_mass(1,iland) / 0.2)**2 ) ** onethird
           call qtk( land%sfcwater_energy(1,iland), tempgr, liqfrg )
        endif

        if (hasveg) then
           sfacm = 1.0 - land%snowfac(iland)

           tai   = tai   * sfacm
           lai   = lai   * sfacm
           hveg  = hveg  * sfacm
           vegfr = vegfr * sfacm
           taiv  = taiv  * sfacm
           laiv  = laiv  * sfacm
        endif

     else

        liqfrg  = 0.0
        deltag  = 0.0

        ! Determine the bare ground temperature
        call qwtk( land%soil_energy(nzg,iland),land%soil_water(nzg,iland)*1.e3, &
                   land%specifheat_drysoil(nzg,iland), tempgr, xm )
     endif

     ! Canopy Wetness

     hasveg = (lai > 0.01 .and. vegfr > 0.01)

     if (hasveg) then
        wrcr   = land%veg_water(iland)
        wrmax  = 0.2 * max(tai,1.e-10)
        deltav = min( 1.0, (wrcr / wrmax)**2 ) ** onethird
        liqfrv = merge(0.0, 1.0, tempvg < 273.16)
     else
        deltav = 0.0
        liqfrv = 1.0
     endif

     ! PH of canopy/ground water

     if (ccnparm < 1.e6) then

        fact = 0.05 * (cldnum(iw) * 1.e-7 - 10.0)
        fact = min( max(fact, 0.0), 1.0 )
        hplus = fact * hplus_dirty + (1.0 - fact) * hplus_clean

     else

        hplus = hplus_def

     endif

     IF ( sfc_hono ) THEN

! Calculate A/V for leaves.
! LAI was multiplied by 2 to account for the fact that surface area
! is provided by both sides of the leaves.
! Matthews Jones, Ammonia deposition to semi-natural vegetation,
! PhD dissertation, University of Dundee, Scotland, 2006

        surf_leaf = 2.0 * tai / zf

! Calculate A/V for buildings and other structures.
! Buildings and other structures can provide additional surfaces in
! urban areas for the heterogeneous reaction to occur. However, such
! information is not readily available; in the absence of such information,
! it is scaled to purb(c,r). Svensson et al., (1987) suggests a typical value
! of 0.2 for A/V for buildings in urban environments. A maximum value of 0.2
! for A/V for buildings is assigned to the grid cell containing the highest
! purb(c,r) i.e., 100.0. A/V for buildings for other grid-cell is calculated
! as purb(c,r)*(0.2/100.0); Cai et al. (2006) used a value of 1.0 for their
! study at New York (total A/V)

        if ( isurban ) then
           surf_bldg = 0.2
        else
           surf_bldg = 0.0
        endif

! Calculate rate constant for the reaction (psudeo-first order reaction,
! unit per second). Calculate pseudo-first order rate constant using Eq 1
! of Vogel et al. (2003).  Unit of KNO2 is in 1/min in the paper; divide it
! by 60 to convert it into 1/sec.

        kno2 = MAX( 0.0, 5.0E-5 * (surf_leaf + surf_bldg) )

     endif

     ! Loop over species to calculate dry deposition velocities.

     n = 0
     dloopl: DO l = 1, n_spc_m3dry

        IF ( .NOT. use_depspc( l ) ) CYCLE dloopl

        n = n + 1

        IF ( l .EQ. l_sulf ) THEN  ! Sulfate (SULF)

           ! Sulfate calculation follows Wesely (1985), Eqn. 11.

           rbsulfi = 0.002 * (ustar**2 + 0.24 * wstar**2) / ustar
           depvel_gas( n,ks ) = depvel_gas( n,ks ) + arcoarkw * &
                                rai * rbsulfi / ( rai + rbsulfi )
        else

           ! Use CMAQ function for calculating the effective Henry's Law
           ! constant.  Note that original M3DRY wants inverse,
           ! non-dimensional Henry's Law (caq/cg).

           heff = hlconst( tempcr, hleff( l ), hplus, hlspc( l ) )

           if ( hleff( l ) .and. vegfr > 0.01 .and. lai > 0.01) then
              heff_ap = hlconst( tempcr, hleff( l ), hplus_ap, hlspc( l ) )
           else
              heff_ap = heff
           endif

           ! Make Henry's Law constant non-dimensional.

           heff  = heff * 0.08205 * tempcr

           ! Wet cuticle resistance.

           IF ( l .NE. l_o3 ) THEN
              rweti = heff * rcw0i
           ELSE
              ! Canopy level wet resistence Rwet to ozone was found to be about
              ! 200 s/m on basis of Keysburg exp. Using LAI(1-sided) of about
              ! 6.25 measured at Keysburg gives leaf level rwet about 1250 s/m
              ! Leaf level rwet from Altimir et al 2006 gives about 1350 s/m
              rweti = cweto3
           END IF

           ! wet (liquid) ground conductance

           IF ( l .Ne. l_O3 ) THEN
              rgliqi = heff * rgwet0i
           else
              rgliqi = 0.002  ! resistance of 500 s/m
           endif

           ! snow / frozen ground conductance

           rsnowi = csnow( l )

           ! wet ground conductance

           rgwi = (1.0 - liqfrg) * rsnowi + liqfrg * rgliqi

           ! dry ground conductance

           IF ( l .Ne. l_O3 ) THEN
              rgdi = cgnd( l )
           else
              rgdi = wsat_vg / (200.0 * wsat_vg + 300.0 * soil_water1)
           endif

           ! Bare ground conductance

           rgndi  = ( 1.0 - deltag ) * rgdi + deltag * rgwi

           ! Vegetation calculations

           if (hasveg) then

              ! Stomatal conductance

              rsi = 1.0 / land%stom_resist(iland)

              ! Canopy conductance

              rinci = ustar / max(14.0 * taiv * hveg, 1.e-4)

              ! Ground conductance under canopy vegetation

              rgndci = rgndi * rinci / (rgndi + rinci)

              ! Wet cuticle conductance with snow/ice

              rwetsfci = ( 1.0 - liqfrv) * rsnowi + liqfrv * rweti

              ! Dry cuticle conductance

              if ( l == l_nh3 ) then
                 rdryi = 0.00025 * exp( 0.054 * rh_air )
              elseif ( l == l_o3 ) then
                 rh_func = max( 0.0, (rh_air - 70.0) / 30.0 )
                 rdryi = ( 1.0 - rh_func ) * ccut( l ) + rh_func * rweti
              else
                 rdryi = ccut( l )
              endif

              ! Total cuticle conductance

              rcuti = taiv * ( (1.0 - deltav) * rdryi + deltav * rwetsfci )

              ! Bulk stomatal and mesophyll conductance (added in series)

              rstomi = ddif0( l ) * rsi
              rmesoi = heff_ap / 3000. + meso( l )

              rstomi = laiv * rmesoi * rstomi / (rstomi + rmesoi)

              ! Bulk surface conductance with vegetation (added in parallel)

              rci =  vegfr * (rstomi + rcuti + rgndci) + (1.0 - vegfr) * rgndi

           else

              rci = rgndi

           endif

           ! Compute dry deposition velocity (aerodynamic, surface, and laminar
           ! conductances added in series)

           rbci = ustar * scc_pr_23( l )

           depvel_gas( n,ks ) = depvel_gas( n,ks ) + arcoarkw * &
                                rci * rai * rbci / (rci*rai + rci * rbci + rai * rbci)

! HONO production via heterogeneous reaction on ground surfaces,
! 2NO2 = HONO + HNO3
! Rate constant for the reaction = (3.0E-3/60)* (A/V),
! where A/V is surface area/volume ratio
! HONO is produced and released into the atmosphere
! NO2 is lost via chemical reaction
! HNO3 is sticky and stays on the surfaces

! Determine NO2 concentration needed for HONO production term.

           IF ( sfc_hono .and. nl%do_emis /= 0 .and. l .EQ. l_no2 ) THEN

! Loss of NO2 via the heterogeneous reaction is accounted as additional
! depositional loss. Add the loss of NO2 via the heterogeneous reaction
! to the regular deposition velocity (increased dep. vel.).  This will
! reduce the NO2 conc. in the atmosphere. Dep vel is adjusted back to the
! original value in vdiffacm2 after NO2 conc is reduced but before calculating
! depositional loss.

              depvel_gas( n,ks ) = depvel_gas( n,ks ) + arcoarkw * 2.0 * kno2 * zf

! Calculate production (pvd) for HONO; unit = ppm * m/s

              pvd = kno2 * cgrid(kw,iw,ic_no2) * zf

! Add this to the HONO emissions for now when we apply it in vertical diffusion

              vdemis_gc(kw,iw,ie_hono) = vdemis_gc(kw,iw,ie_hono) &
                                       + pvd * sfcg%area(iwsfc) * volti(kw,iw)

!BOB: Does the line above needs replacement for sfcg%area?

           ENDIF

        END IF   ! special condition for sulfate (SULF)

     END DO dloopl   ! (l = 1, n_spc_m3dry)

  enddo   ! land

!-------------------------------------------------------------------------------
! Loop over sea cells and calculate dry deposition.
!-------------------------------------------------------------------------------

  do jsfc = itab_w(iw)%jsea1, itab_w(iw)%jsea2
     iwsfc = itab_w(iw)%iwsfc(jsfc)
     jasfc = itab_w(iw)%jasfc(jsfc)

     kw = itab_wsfc(iwsfc)%kwatm(jasfc)

     isea = iwsfc - omsea

     ks = kw - lpw(iw) + 1

     arcoarkw = itab_wsfc(iwsfc)%arcoarkw(jasfc)

     isocean = (sfcg%leaf_class(iwsfc) == 0)

     ustar  = sfcg%ustar(iwsfc)
     ustarw = sea%sea_ustar(isea)
     rawi   = sea%sea_ggaer(isea)

     sst   = sea%seatc(isea)
     tw    = sea%sea_cantemp(isea) ! water surface film temperature
     wtv0  = sfcg%wthv(iwsfc)

     if (wtv0 > 0.0) then
        wstar = (grav * pblh(iw) * wtv0 / theta(kw,iw)) ** onethird
     else
        wstar = 0.0
     endif

     kviswi = EXP( 0.025 * ( tw - stdtemp ) ) / 0.017

     if  (sea%nlev_seaice(isea) > 0 ) then
        ustari = sea%ice_ustar(isea)
        raicei = sea%ice_ggaer(isea)
        tice   = sea%ice_cantemp(isea) ! ice surface film temperature

        ! Determine the seaice liquid water mass fraction
        call qtk_sea( sea%seaice_energy(sea%nlev_seaice(isea),isea), xt, xm)
     endif

     ! Loop over species to calculate dry deposition velocities.

     n = 0

     dloops: DO l = 1, n_spc_m3dry

        IF ( .NOT. use_depspc( l ) ) CYCLE dloops

        n = n + 1

        IF ( l .EQ. l_sulf ) THEN  ! Sulfate (SULF)

           ! Sulfate calculation follows Wesely (1985), Eqn. 11.

           rbsulfi = 0.002 * (ustar**2 + 0.24 * wstar**2) / ustar

           depvel_gas( n,ks ) = depvel_gas( n,ks ) + arcoarkw * &
                                rawi * rbsulfi / ( rawi + rbsulfi )

        else

           ! Use CMAQ function for calculating the effective Henry's Law
           ! constant.  Note that original M3DRY wants inverse, non-
           ! dimensional Henry's Law (caq/cg).   Water pH is different
           ! than rain, and we need to use the water temperature.

           heff  = hlconst( tw, hleff( l ), hplus_h2o, hlspc( l ) )

           ! Make Henry's Law constant non-dimensional.

           heff  = heff * 0.08205 * tw

           ! from Hayduk and Laudie

           dw        = dw25( l ) * ( tw * rt25inK ) * ( 0.009025 * kviswi )
           scw_pr_23 = d3 * ( (dw * pr * kviswi)**2 ) ** onethird

           if (l == l_o3) then  ! implement Chang et al(2004)

              ! pChang is a/H or alpha/H which would be 1/H in current model
              ! note that in Chang et al (2004) and Garland et al (1980), their
              ! H is Cair/Cwater wich is the inverse of heff

              pChang  = 1.75
              kwChang = ustarw * scw_pr_23

              ! If a file of chlorophyll concentrations is provided, Iodide concentration
              ! are estimated from a fit to the Rebello et al 1990 data. The slope and
              ! correlation are given in the paper but not the intercept, so the data in
              ! Tables 3 & 4 were fit to get the relationship below. The regression gives
              ! the concentration in umol/L and is converted to mol/L for use in
              ! Chang et al eq. The slope and correlation are a slightly different than
              ! in Table 5. If chlorophyll concs are not available, a constant value
              ! for [I-] of 100e-9 mol/l is used Use ocean file variables to determine if
              ! the water cell is ocean or lake; method is only for ocean cells

              if (isocean) then
                 ! Iodide in sea-water based on SST  (mol /dm-3)
                 ciodide = 1.46E6 * EXP( -9134.0 / sst )
                 qiodide = sqrt ( 2.0e9 * ciodide * dw * 1.e-4 ) * heff
              else
                 qiodide = 0.0
              endif

              rsurfi = pChang * kwchang + qiodide

           else  ! not ozone

              rsurfi = heff * ustarw * scw_pr_23

           endif

           ! Compute dry deposition velocity over sea

           rbci = ustarw * scc_pr_23( l )

           depvel_gas( n,ks ) = depvel_gas( n,ks ) + arcoarkw * (1.0 - sea%seaicec(isea)) * &
                                rsurfi * rawi * rbci / ( rsurfi * rawi + rsurfi * rbci + rawi * rbci )

           ! Include fractional contribution of seaice if present

           if  (sea%nlev_seaice(isea) > 0 ) then

              ! Use CMAQ function for calculating the effective Henry's Law
              ! constant.  Note that original M3DRY wants inverse,
              ! non-dimensional Henry's Law (caq/cg).

              heff  = hlconst( tice, hleff( l ), hplus_h2o, hlspc( l ) )

              ! Make Henry's Law constant non-dimensional

              heff  = heff * 0.08205 * tice

              ! Wet surface resistance

              rgwi  = rgwet0i * heff

              ! Assume seaice behaves similarly to snowcover

              rsurfi = (1.0 - xm) * csnow(l) + xm * rgwi

              ! Compute dry deposition velocity contribution from seaice

              rbci  = ustari * scc_pr_23( l )

              depvel_gas( n,ks ) = depvel_gas( n,ks ) + arcoarkw * sea%seaicec(isea) * &
                   rsurfi * raicei * rbci / ( rsurfi * raicei + rsurfi * rbci + raicei * rbci )

           endif  ! seaice

        endif  ! special condition for sulfate (SULF)

     END DO dloops   ! (l = 1, n_spc_m3dry)

  enddo   ! sea

END SUBROUTINE m3dry
