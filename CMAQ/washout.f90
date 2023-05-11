subroutine precip_scav( conc, dconc, remov, airtemp, airdens, airpres, &
                        dz, f_act, hplus, pcpflx, trflx, f_liq, wcbar )

  use aero_data
  use aq_data
  use const_data,  only: inv_mwairkg, molvol, stdtemp, mwairkg
  use cgrid_spcs,  only: nspcsd
  use hlconst_mod, only: hlconst

  implicit none

  real,    intent(in ) :: conc (nspcsd+1)
  real,    intent(in ) :: trflx(nspcsd+1)
  real,    intent(out) :: dconc(nspcsd+1)
  real,    intent(out) :: remov(nspcsd+1)
  real,    intent(in ) :: airdens, airpres, airtemp, f_act, hplus
  real,    intent(in ) :: pcpflx, dz, f_liq, wcbar

  real    :: scavfrac(n_mode,3), fact, facti, gasfac, conv
  real    :: rtch, rate, cwat
  integer :: spc, sp2, i, pntr
  real    :: lwcv, conc_air
  real    :: kh(n_gas_scav)
  real    :: zfact, dconc_air
  logical :: doaer, dogas

  real, parameter :: rholiq  = 1.e+3
  real, parameter :: rholiqi = 1.e-3
  real, parameter :: rgas = molvol / stdtemp

  dconc = 0.0
  remov = 0.0

  ! skip if precip rate is small
  if (pcpflx < 1.e-12) return

  doaer = (f_act < .99)
  dogas = (f_liq > .01)

  ! converts 1/mol_air to 1/m^2
  fact  = airdens * dz * inv_mwairkg
  facti = 1.0 / fact

  ! skip if no washout (no interstitial aerosol and no liquid precip for gas uptake)
  if ((.not. doaer) .and. (.not. dogas)) return

  call calc_scav(conc, scavfrac, gasfac, airtemp, airdens, airpres, pcpflx, doaer, dogas)

  ! Washout of aerosols

  if (doaer) then

     if (f_act > .001) then
        scavfrac(2:3,:) = scavfrac(2:3,:) * (1.0 - f_act)
     endif

     do i = 1, naitkn
        spc = iaitkn_map    (i)
        sp2 = iaitkn_2_accum(i)
        dconc(spc) = - conc (spc) * scavfrac(1,3)
        remov(sp2) = - dconc(spc) * fact
     enddo

     do i = 1, naccum
        spc = iaccum_map(i)
        dconc(spc) = - conc (spc) * scavfrac(2,3)
        remov(spc) = remov(spc) - dconc(spc) * fact
     enddo

     do i = 1, ncoars
        spc = icoars_map(i)
        dconc(spc) = - conc (spc) * scavfrac(3,3)
        remov(spc) = - dconc(spc) * fact
     enddo

     do i = 1, 3
        spc = num_str + i - 1
                    dconc(spc) = - conc (spc) * scavfrac(i,1)
        if (i /= 1) remov(spc) = - dconc(spc) * fact

        spc = srf_str + i - 1
                    dconc(spc) = - conc (spc) * scavfrac(i,2)
        if (i /= 1) remov(spc) = - dconc(spc) * fact
     enddo

     spc = srf_str
     sp2 = srf_str + 1
     remov(sp2) = remov(sp2) - dconc(spc) * fact

  endif

  ! Washout of gases

  if (dogas) then

     rate = gasfac * f_liq * dz

     ! get "nondimensional" Henry's law constants for gases (m^3 air / m^3 water)
     RTCH = RGAS * AIRTEMP
     DO SPC = 1, N_GAS_SCAV
        KH(SPC) = RTCH * HLCONST( AIRTEMP, HLEFF_MAP( SPC ), HPLUS, HLSPC_MAP( SPC ) )
     ENDDO

     ! per mol air to per vol
     conv = airdens * inv_mwairkg

     ! liquid water flux (m^3 h2o / m^3 air) (m/s)
     lwcv = pcpflx * f_liq * rholiqi

     ! cloud water
     cwat = wcbar * rholiqi

     DO SPC = 1, N_GAS_SCAV
        PNTR = GAS_SCAV_MAP( SPC )

        ! ambient concentration converted to mol/m^3
        conc_air = conc(pntr) * conv

        ! remove fraction of tracer in cloud water
        if (wcbar > 1.e-12) then
           conc_air = conc_air / (1.0 + kh(spc) * cwat)
        endif

!       zfact = 1.0 - exp(-gasfac * kh(spc))
        zfact = rate / kh(spc)
        if (zfact > 2.e-5) zfact = 1.0 - exp(-zfact)

!       ! exponential function...
!       tfact = kh(spc) * lwcv * zfact * dt / dz
!       if (tfact > 2.e-5) tfact = 1.0 - exp(-tfact)
!
!       ! equilibrium concentration from Henry's Law
!       conc_eq = trflx(pntr) / (kh(spc) * lwcv)
!
!       ! change in ambient concentration (mol/m^2)
!       dconc_air = -(conc_eq - conc_air) * tfact * dz

        ! change in ambient concentration (mol/m^2/s)
        dconc_air = (trflx(pntr) - conc_air * kh(spc) * lwcv) * zfact

        remov(pntr) = - dconc_air
        dconc(pntr) = dconc_air * facti
     ENDDO

  endif

end subroutine precip_scav




subroutine calc_scav(conc, scav, gasfac, airtemp, airdens, airpres, pcpflx, &
                     doaer, dogas)

  use cgrid_spcs, only: nspcsd, ae_strt, ae_fini
  use const_data
  use getpar_mod, only: getpar, getdens_conc
  use cgrid_conv, only: rev_cgrid_one
  use aq_data,    only: convb
  use aero_data,  only: aer_str, aer_end, aer_num, aer_m3, mode_map, &
                        aer_trac, n_mode, aeromode, m0_min, m2_min, &
                        num_str, srf_str

  implicit none

  real,    intent(in ) :: conc(nspcsd+1)
  real,    intent(in ) :: airdens, airpres, airtemp, pcpflx
  real,    intent(out) :: scav(n_mode,3), gasfac
  logical, intent(in ) :: doaer, dogas

  real :: m3(aer_num)
  real :: cngrd(nspcsd)

  real :: aeromode_mass(n_mode)
  real :: aeromode_dens(n_mode)
  real :: aeromode_lnsg(n_mode)
  real :: aeromode_diam(n_mode)
  real :: moment0_conc (n_mode)
  real :: moment2_conc (n_mode)
  real :: moment3_conc (n_mode)
  real :: sqrt2_lnsg   (n_mode)
  real :: volm1        (n_mode)
  real :: sfcm1        (n_mode)
  real :: effic        (3)

  integer :: i, n, s

  real :: xlm, amu, wmu, pcpodmb, pcpvel
  real :: dmb, re, a1r, sstar, dmbi, muratio, tau, st, dratio, stcorr
  real :: collect1, collect2, dconst1, tau0, cc, densratio, sre

! integer, parameter :: nn = 10
!
! real :: x(nn,2) = reshape ( [ &
!          -3.436159119,   7.6404329E-6, &
!          -2.532731674,   0.00134364575, &
!          -1.756683649,   0.03387439446, &
!          -1.03661083,    0.2401386111, &
!          -0.3429013272,  0.6108626337, &
!           0.3429013272,  0.610862634, &
!           1.03661083,    0.2401386111, &
!           1.756683649,   0.03387439446, &
!           2.532731674,   0.001343645747, &
!           3.436159119,   7.64043286E-6], [nn,2], order=[2,1] )

  integer, parameter :: nn = 16

  real :: x(nn,2) = reshape ( [ &
           -4.688738939,    2.65480747E-10, &
           -3.869447905,    2.32098084E-7, &
           -3.176999162,    2.711860093E-5, &
           -2.546202158,    9.32284009E-4, &
           -1.951787991,    0.0128803115, &
           -1.380258539,    0.0838100414, &
           -0.822951449,    0.280647459, &
           -0.2734810461,   0.507929479, &
            0.273481046,    0.507929479, &
            0.8229514491,   0.2806474585, &
            1.380258539,    0.0838100414, &
            1.951787991,    0.01288031154, &
            2.546202158,    9.32284009E-4, &
            3.176999162,    2.71186009E-5, &
            3.869447905,    2.320980845E-7, &
            4.688738939,    2.65480747E-10], [nn,2], order=[2,1] )

! integer, parameter :: nn = 20
!
! real :: x(nn,2) = reshape ( [ &
!          -5.38748089,     2.229393646E-13, &
!          -4.60368245,     4.39934099E-10, &
!          -3.94476404,     1.086069371E-7, &
!          -3.347854567,    7.80255648E-6, &
!          -2.788806058,    2.28338636E-4, &
!          -2.254974002,    0.003243773342, &
!          -1.738537712,    0.02481052089, &
!          -1.234076215,    0.109017206, &
!          -0.7374737286,   0.2866755054, &
!          -0.2453407083,   0.4622436696, &
!           0.245340708,    0.4622436696, &
!           0.7374737286,   0.2866755054, &
!           1.234076215,    0.109017206, &
!           1.738537712,    0.0248105209, &
!           2.254974002,    0.003243773342, &
!           2.788806058,    2.28338636E-4, &
!           3.347854567,    7.80255648E-6, &
!           3.94476404,     1.08606937E-7, &
!           4.60368245,     4.399340992E-10, &
!           5.38748089,     2.22939365E-13], [nn,2], order=[2,1] )

  real :: dp(nn)
  real :: ecollect(nn)

  real :: vsed, schm, dg, fv

  real, parameter :: t0      = 288.15      ! [ K ] ! starting standard surface temp.
  real, parameter :: one3    = 1./ 3.
  real, parameter :: two3    = 2./ 3.
  real, parameter :: sqrt2   = sqrt(2.0)
  real, parameter :: threepi = 3.0 * pi
  real, parameter :: sqrtpi  = sqrt(pi)
  real, parameter :: spim    = 1.0 / sqrtpi
  real, parameter :: dg0     = 1.e-5
  real, parameter :: dg00    = dg0 * stdatmpa * stdtemp ** (-1.75)

  ! Calculate air dynamic viscosity [ kg m**-1 s**-1 ]
  AMU = 1.458E-6 * AIRTEMP * SQRT( AIRTEMP ) / ( AIRTEMP + 110.4 )

  dmb    = 0.0025                ! rain diameter
  pcpvel = 7.33 / sqrt(airdens)  ! rain fall speed
  dmbi   = 1.0 / dmb

  re  = 0.5 * dmb * pcpvel * airdens / amu
  sre = sqrt(re)

  ! Compute gas scavenging factor

  if (dogas) then

     ! diffusivity of gases in air
     dg = dg00 * airtemp**2 * airtemp**(-0.25) / airpres

     ! ventilation factor (Pruppacher and Klett)
     fv = 0.78 + 0.308 * sre * (amu / (airdens * dg)) ** one3

     ! gas scavenging factor
     gasfac = 12. * dg * dmbi**2 * fv / pcpvel

  endif

  ! Compute erosol scavenging factor

  if (doaer) then

     ! convert CGRID aerosols to densities
     cngrd(ae_strt:ae_fini) = conc(ae_strt:ae_fini) * convb(ae_strt:ae_fini)
     call rev_cgrid_one(cngrd, airdens)

     ! compute aerosol moments directly from concentration array

     m3 = cngrd(aer_str:aer_end) * aer_m3

     do i = 1, 3
        n = num_str + i - 1
        s = srf_str + i - 1

        moment3_conc( i ) = Max( sum(m3, mask=(mode_map==i .and. .not. aer_trac) ), &
                                   aeromode( i )%min_m3conc )
        moment0_conc( i ) = Max( cngrd( n ), m0_min( i ) )
        moment2_conc( i ) = Max( cngrd( s ) * oneovpi, m2_min( i ) )
     enddo

     ! compute aerosol diameters and geometric standard deviation

     call getpar( .false., moment0_conc, moment2_conc, moment3_conc, &
                  aeromode_lnsg, aeromode_diam, 1, 3 )

     do i = 1, 3
        sqrt2_lnsg(i) = aeromode_lnsg(i) * sqrt2
        volm1     (i) = 1.0 / (sqrtpi * aeromode_diam(i)**3  * exp(4.5 * aeromode_lnsg(i)**2))
        sfcm1     (i) = 1.0 / (sqrtpi * aeromode_diam(i)**2  * exp(2.0 * aeromode_lnsg(i)**2))
     enddo

     ! Compute mean particle densities
     call getdens_conc(aeromode_mass, aeromode_dens, cngrd, moment3_conc, 1, 3)

     ! Calculate mean free path [ m ]:
     XLM = 6.6328E-8 * STDATMPA * AIRTEMP / ( T0 * AIRPRES )

     ! Dynamic viscosity of liquid water
     wmu = 0.509528 * (max(airtemp,260.) - 229.7) ** (-1.5)

     ! Collection of terms in particle Schmidt number
     dconst1 = threepi * amu **2 / (boltzmann * airtemp * airdens)

     a1r   = log(1. + re)
     sstar = (1.2 + 0.0833333 * a1r) / (1. + a1r)

     pcpodmb = 1.5e-3 * pcpflx * dmbi

     muratio = amu / wmu

     do i = 1, n_mode

        tau0 = (aeromode_dens(i) - airdens) / (18. * amu)

        densratio = sqrt(aeromode_dens(i) * 1.e-3)

        do n = 1, nn

           ! diameter of particles at the Gaussian weighting points
           dp(n) = aeromode_diam(i) * exp( sqrt2_lnsg(i) * x(n,1) )

           ! Cunningham slip correction
           cc = 1.0 + 2.0 * xlm / dp(n) * (1.257 + 0.4 * exp( -0.55 * dp(n) / xlm ))

           ! particle characteristic fall relaxation time
           tau = tau0 * dp(n) * dp(n) * cc

           ! particle fall velocity
           vsed = tau * grav

           ! Stokes number of particle
           st = 2. * tau * (pcpvel - vsed) * dmbi

           ! Schmidt number of particle
           schm = dconst1 * dp(n) / cc

           ! Collection efficiency due to Brownian diffusion
           collect1 = (4.0 + sre * (1.6 * schm**one3 + 0.64 * sqrt(schm))) / (re * schm)

           ! Collection efficiency due to interception
           dratio = dp(n) * dmbi
           collect2 = 4. * dratio * muratio + (4. + 8. * sre) * dratio**2

           ! Collection efficiency due to impaction
           stcorr = 0.0
           if (st > sstar) then
              stcorr = ((st - sstar) / (st - sstar + two3))**1.5 * densratio
           endif

           ! Total collection efficiency
           ecollect(n) = min(collect1 + collect2 + stcorr, 1.0)
        enddo

        ! Sum number, area, and mass efficiencies over all Gaussian integration points
        effic(1) = spim     * sum( x(:,2) * ecollect(:) )
        effic(2) = sfcm1(i) * sum( x(:,2) * ecollect(:) * dp(:)**2 )
        effic(3) = volm1(i) * sum( x(:,2) * ecollect(:) * dp(:)**3 )

        ! Below-cloud scavenging factors
        scav(i,1:3) = min(effic(1:3), 1.0) * pcpodmb
     enddo

  endif

END SUBROUTINE calc_scav
