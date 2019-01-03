!------------------Subroutine deposition_driver--------------
!This subroutine calculates dry deposition, rain and in-cloud
!scavenging. Originally written for RAMS 4.3.0 by Nair et al.
!Land surface types updated to RAMS 6.0 by Rob Seigel (2010). 
!Steve Saleeby (2011) checked code and updated wet scavenging to 
!match collection efficiency and scavenging rates documented in
!Seinfeld and Pandis (2006). Terms were adjusted in the collection
!efficiency for wet scavenging to be appropriate for multiple 
!aerosol types. Particle sizes and fall speeds were derived from
!RAMS parameters. 

!==============================================================================
subroutine nuclei_deposition(iw, k1, k2, dtl0, voa, rhoa, press, tair, tairc, &
   dynvisc, rhov, rhovslair, dmb, pcpvel, pcpfluxr, con_ccnx, con_gccnx, con_ifnx)

use micro_coms,  only: ncat, iccn, igccn, iifn
use ccnbin_coms, only: nccntyp, nnuc, rhow, rho_nucx, diam_nucx
use mem_grid,    only: mza, lpw, lsw, zm, zt, dzt_bot, dzit, zfacm2, zfacim2, &
                       arw, volti
use mem_ijtabs,  only: itab_w
use mem_sea,     only: sea, itab_ws
use mem_leaf,    only: land, itab_wl
use consts_coms, only: r8, grav, vonk, gravi
use nuclei_coms, only: jsfcinx
use mem_basic,   only: vxe, vye, vze
use nuclei_coms, only: e0, acoll, alpha, ggamma, jsfcinx, &
                       pi3, boltz, t_suth, t_ref, p_ref, dva_ref, airmfp_ref
use mem_cuparm,  only: conprr, kcutop, kcubot

implicit none

  integer, intent(in) :: iw      ! Grid horizontal index

  integer, intent(in) :: k1(11) ! Lowest grid level of each precipitation type
  integer, intent(in) :: k2(11) ! Highest grid level of each precipitation type

  real, intent(in) :: dtl0   ! Model long timestep [s]

  real, intent(in) :: voa      (mza) ! Ratio of source cell volume to top horizontal area
  real(r8), intent(in) :: rhoa (mza) ! Density of air [kg/m^3]
  real, intent(in) :: press    (mza) ! Air pressure [Pa]
  real, intent(in) :: tair     (mza) ! Air temperature [K]
  real, intent(in) :: tairc    (mza) ! Air temperature [deg C]
  real, intent(in) :: dynvisc  (mza) ! Air dynamic viscosity [kg/(m s)]
  real, intent(in) :: rhov     (mza) ! Water vapor density [kg/m^3]
  real, intent(in) :: rhovslair(mza) ! Saturation water vapor density [kg/m^3]

  real, intent(inout) :: dmb      (mza,ncat) ! Diameter of mean-mass precipitation [m]
  real, intent(inout) :: pcpvel   (mza,ncat) ! Precipitation fall velocity [m/s]
  real, intent(inout) :: pcpfluxr (mza,ncat) ! Precipitation flux rate [kg/(m^2 s)]

  real, intent(inout) :: con_ccnx (mza,nccntyp) ! CCN  number concentration [#/m^3]
  real, intent(inout) :: con_gccnx(mza)         ! GCCN number concentration [#/m^3]
  real, intent(inout) :: con_ifnx (mza)         ! IFN  number concentration [#/m^3]

  integer :: lcat ! Precipitation category index
  integer :: k    ! grid vertical index
  integer :: inuc ! Nuclei category index
  integer :: ic   ! CCN category index
  integer :: nsea
  integer :: jws
  integer :: iws
  integer :: kw
  integer :: nland
  integer :: jwl
  integer :: iwl
  integer :: leaf_class
  integer :: isfcinx    ! sfc type index
  real    :: wwet(mza)  ! Weight for diam_nucx table interpolation
  integer :: iwet(mza)  ! Index for diam_nucx table interpolation

  integer :: ks1(11), ks2(11) ! Copy of k1, k2

  real :: re       (mza,ncat) ! Reynolds number of precipitating hydrometeors [ ]
  real :: sqrt_re  (mza,ncat) ! Square root of Reynolds number [ ]
  real :: sstar    (mza,ncat) ! parameter used to calculate collection efficient
  real :: pcpodmb  (mza,ncat) ! Ratio of precipitation flux rate to hydrometeor diameter

  real :: sqrt_rhoi_nuc (mza) ! Square root of (rhow / rho_nucx)

  real :: vels          (mza) !  wind speed for all levels
  real :: airmfp        (mza) ! air free path [m]
  real :: dynvisc_liqwat(mza) ! dynamic viscosity of liquid water [kg/(m s)] 
  real :: sqrt_schm     (mza) ! Square root of Schmidt number [ ]
  real :: cbrt_schm     (mza) ! Cube root of Schmidt number [ ]
  real :: scavfrac      (mza) ! scavenging fraction of nuclei this timestep [ ]

  real :: con_nucx   (mza,nnuc) ! aerosol concentration [#/m^3]
  real :: dwetnuc    (mza,nnuc) ! nuclei wet diameter [m]
  real :: rho_wetnuc (mza,nnuc) ! nuclei wet density [kg/m^3]
  real :: slipc      (mza,nnuc) ! slip correction factor
  real :: diff_nuc   (mza,nnuc) ! Nuclei diffusion coefficient due to Brownian motion
  real :: schm       (mza,nnuc) ! Schmidt number of nuclei [ ]
  real :: tau        (mza,nnuc) ! nuclei inertial relaxation time; tau*g = grav settling vel
  real :: vwetnuc    (mza,nnuc) ! nuclei gravitional settling velocity
  real :: xfernuc    (mza,nnuc) ! number transferred across arw by dry deposition
  real :: xfernuc_sfc(mza,nnuc) ! number transferred to land & sea cells by dry deposition
  real :: sourcec    (mza,nnuc) ! number of nuclei in source cell per m^2 of source cell TOP area

  real :: vels10    ! wind speed at 10 m reference height
  real :: st        ! Stokes number
  real :: stcorr    ! Stokes correction for collection efficiency
  real :: collect1  ! Nuclei collection efficiency
  real :: collect2  ! Nuclei collection efficiency
  real :: ecollect  ! Nuclei collection efficiency
  real :: ra        ! aerodynamic resistance [m/s]
  real :: eb        ! contribution to collection efficiency from Brownian diffusion
  real :: ein       ! " " from interception
  real :: eim       ! " " particle inertia/impaction
  real :: r1        ! Rebound fraction
  real :: rs        ! surface resistance
  real :: vnucsfc   ! surface deposition velocity [m/s]
  real :: dratio    ! ratio of nuclei to precip diameter [ ]
  real :: slipc_dry ! Slip correction factor of dry nuclei
  real :: vdrynuc   ! gravitational settling velocity of dry nuclei
  real :: cd        ! drag coefficient over water (Baron & Willeke (2001) 4-23
  real :: rvd       ! auxiliary variable
  real :: a1r       ! auxiliary variable
  real :: dzicld    ! Inverse of cloud vertical thickness in convective param [1/m] 
  real :: xkcprime
  real :: xkdprime
  real :: xkc
  real :: xkd
  real :: relhum    ! relative humidity
  real :: fracwkk   ! nuclei fractional fall distance out of source cell

  real, parameter :: m3l10 = -3.0 * log(10.0)

  ! Copy nuclei concentrations (con_nucx) from microphysics column arrays

  inuc = 0

  if (iccn >= 2) then
     do ic = 1,nccntyp
        inuc = inuc + 1
        con_nucx(:,inuc) = con_ccnx(:,ic)
     enddo
  endif

  if (igccn == 2) then
     inuc = inuc + 1
     con_nucx(:,inuc) = con_gccnx(:)
  endif

  if (iifn == 2) then
     inuc = inuc + 1
     con_nucx(:,inuc) = con_ifnx(:)
  endif

  ! Loop over T levels in grid column

  do k = lpw(iw), mza

     ! Get interpolation index and weight for diam_nucx table

     relhum = min(1.,rhov(k) / rhovslair(k))

     if (relhum < 0.7) then
        iwet(k) = 1
        wwet(k) = 0.
     elseif (relhum < 0.9) then
        iwet(k) = 1
        wwet(k) = 5. * (relhum - 0.8)
     elseif (relhum < 0.97) then
        iwet(k) = 2
        wwet(k) = 14.2857 * (relhum - 0.9)
     elseif (relhum < 0.99) then
        iwet(k) = 3
        wwet(k) = 50. * (relhum - 0.97)
     elseif (relhum < 0.997) then
        iwet(k) = 4
        wwet(k) = 142.857 * (relhum - 0.99)
     elseif (relhum < 0.999) then
        iwet(k) = 5
        wwet(k) = 500. * (relhum - 0.997)
     elseif (relhum < 1.0) then
        iwet(k) = 6
        wwet(k) = 1000. * (relhum - 0.999)
     else
        iwet(k) = 6
        wwet(k) = 1.
     endif

     ! Wind speed

     vels(k) = sqrt(vxe(k,iw)**2 + vye(k,iw)**2 + vze(k,iw)**2)

     ! Air mean free path [micron], Baron and Willeke(2001) eq.(4-6)

     airmfp(k) = airmfp_ref * (p_ref / press(k)) * (tair(k) / t_ref) &
             * ((1. + t_suth / t_ref) / (1. + t_suth / tair(k)))

     ! Dynamic viscosity of liquid water: linear fit valid at 0 C and 20 C

     dynvisc_liqwat(k) = 1.787e-3 - 0.03925e-3 * max(0.,min(20.,tairc(k)))

  enddo

  ! Loop over microphysics aerosol types

  do inuc = 1,nnuc

     ! Loop over T levels in grid column

     do k = lpw(iw), mza

        ! Diagnose aerosol wet diameter by interpolation from diam_nucx table based on r.h.

        dwetnuc(k,inuc) = diam_nucx(iwet(k)  ,inuc) * (1. - wwet(k)) &
                        + diam_nucx(iwet(k)+1,inuc) *       wwet(k)

        rho_wetnuc(k,inuc) = (diam_nucx(1,inuc) / dwetnuc(k,inuc))**3 * (rho_nucx(inuc) - rhow) + rhow

        ! Slip correction factor from Seinfeld and Pandis (2006) eq.(9.34)

        slipc(k,inuc) = 1.0 + 2.0 * airmfp(k) / dwetnuc(k,inuc) &
              * (1.257 + 0.4 * exp(-0.55 * dwetnuc(k,inuc) / airmfp(k)))

        ! Aerosol diffusivity [m2/s], Seinfeld and Pandis (2006) eq. (9.73)
        ! note dp is in unit of microns so we correct with a 1.e6 factor

        diff_nuc(k,inuc) = boltz * tair(k) * slipc(k,inuc) / (pi3 * dynvisc(k) * dwetnuc(k,inuc))

        schm(k,inuc) = dynvisc(k) / (real(rhoa(k)) * diff_nuc(k,inuc))

        ! Characteristic relaxation time of particle

        tau(k,inuc) = rho_wetnuc(k,inuc) * dwetnuc(k,inuc)**2 * slipc(k,inuc) / (18. * dynvisc(k))

        vwetnuc(k,inuc) = grav * tau(k,inuc)

        ! Nuclei number in source grid cell per m^2 of source grid cell BOTTOM horiz area

        sourcec(k,inuc) = con_nucx(k,inuc) * voa(k) * zfacim2(k) * zfacm2(k-1)

        ! Fractional fall distance of source grid cell nuclei [across arw(k-1)]

        fracwkk = min(1.0, vwetnuc(k,inuc) * dtl0 * dzit(k))

        ! Nuclei number transfer [#] across arw(k-1,iw) this timestep (positive downward)

        xfernuc(k-1,inuc) = sourcec(k,inuc) * fracwkk * arw(k-1,iw)

     enddo ! k

  enddo ! inuc

  ! Initialize surface deposition number to zero for all grid levels and nuc types

  xfernuc_sfc(:,:) = 0.

  ! Check for sea area beneath this atmospheric grid column

  nsea = itab_w(iw)%nsea
  if (nsea > 0) then

     ! Loop over sea cells beneath this atmospheric grid column

     do jws = 1,nsea
        iws = itab_w(iw)%isea(jws)
        kw = itab_ws(iws)%kw

        ! Diagnose wind speed at 10 m height

        vels10 = vels(kw) * log(10.         / sea%sea_rough(iws)) &
                          / log(dzt_bot(kw) / sea%sea_rough(iws))

        ! Loop over microphysics aerosol types

        do inuc = 1,nnuc

           if (vels10 <= 1.e-3) then
              vnucsfc = vwetnuc(kw,inuc)
           else

              ! Deposition velocity based on the Slinn and Slinn (1980)

              slipc_dry = 1.0 + 2.0 * airmfp(kw) / diam_nucx(1,inuc) &
                 * (1.257 + 0.4 * exp(-0.55 * diam_nucx(1,inuc) / airmfp(kw)))

              vdrynuc = rho_nucx(inuc) * diam_nucx(1,inuc)**2 * grav * slipc_dry &
                      / (18. * dynvisc(kw))

              rvd = rho_wetnuc(kw,inuc) * vwetnuc(kw,inuc) * dwetnuc(kw,inuc)

              cd = 24. * dynvisc(kw) / rvd + .4704

              st = vwetnuc(kw,inuc) * gravi * sea%ustar(iws)**2 * real(rhoa(kw)) / dynvisc(kw)

              xkcprime = cd * vels10 / (1. - vonk)
              xkdprime = cd * vels10 * (1. / sqrt(schm(kw,inuc)) + exp( m3l10 / st)) / vonk

              xkc = xkcprime + vdrynuc
              xkd = xkdprime + vwetnuc(kw,inuc)
              vnucsfc = xkc * xkd / (xkd + xkc - vdrynuc)

           endif

           ! Fractional fall distance of source grid cell nuclei [to sea cell]

           fracwkk = min(1.0, vnucsfc * dtl0 * dzit(kw))

           ! Nuclei number transfer [#] to this sea cell this timestep (positive downward)

           xfernuc_sfc(kw-1,inuc) = xfernuc_sfc(kw-1,inuc) &
                                  + sourcec(kw,inuc) * fracwkk * sea%area(iws)

        enddo ! inuc

     enddo ! jws

  endif ! nsea > 0

  ! Check for land area beneath this atmospheric grid column

  nland = itab_w(iw)%nland
  if (nland > 0) then

     ! Loop over land cells beneath this atmospheric grid column

     do jwl = 1,nland
        iwl = itab_w(iw)%iland(jwl)
        kw = itab_wl(iwl)%kw

        ! Diagnose wind speed at 10 m height

        vels10 = vels(kw) * log(10.         / land%rough(iwl)) &
                          / log(dzt_bot(kw) / land%rough(iwl))

        leaf_class = land%leaf_class(iwl)
        isfcinx = jsfcinx(leaf_class)

        ! aerodynamic resistance
        Ra = vels10 / land%ustar(iwl)**2

        ! Loop over microphysics aerosol types

        do inuc = 1,nnuc

           ! Dry deposition velocity on different surfaces
           !   vnucsfc = Vg + 1/(Ra+Rs), Vg: gratitational settling,
           !                             Ra: aerodynamic resistances
           !                             Rs: surface resistance
           ! Both surface and aerodynamic resistance are factored into the net
           ! fallspeed.  Based on Slinn (1981)

           Eb = schm(kw,inuc) ** (-ggamma(isfcinx))

           if (isfcinx /= 8  .and. isfcinx /= 9 .and. isfcinx /= 12) then
              St = vwetnuc(kw,inuc) * land%ustar(iwl) / (acoll(isfcinx) * grav * 1.0e-3)
              Ein = 0.5 * (dwetnuc(kw,inuc) * 1.0e-6 / (acoll(isfcinx) * 1.0e-3))**2
           else
              St = vwetnuc(kw,inuc) * land%ustar(iwl)**2 * real(rhoa(kw)) * gravi / dynvisc(kw)
              Ein = 0.0
           endif

           Eim = (St / (alpha(isfcinx) + St))**2

           if (isfcinx == 0) then
              R1 = 1. ! No rebound from wet surfaces
           else
              R1 = exp ( -sqrt(St) )
           endif
           Rs = 1. / (e0 * land%ustar(iwl) * (Eb + Eim + Ein) * R1)

           vnucsfc = vwetnuc(kw,inuc) + 1. / (Ra + Rs)

           ! Fractional fall distance of source grid cell nuclei [to land cell]

           fracwkk = min(1.0, vnucsfc * dtl0 * dzit(kw))

           ! Nuclei number transfer [#] to this land cell this timestep (positive downward)

           xfernuc_sfc(kw-1,inuc) = xfernuc_sfc(kw-1,inuc) &
                                  + sourcec(kw,inuc) * fracwkk * land%area(iwl)

        enddo ! inuc

     enddo ! jwl

  endif ! nland > 0

  ! Loop over all grid W levels in this column and apply fluxes to nucx concentrations

  do k = lpw(iw),mza-1
     do inuc = 1,nnuc
        con_nucx(k+1,inuc) = con_nucx(k+1,inuc) - xfernuc(k,inuc) * volti(k+1,inuc)
        con_nucx(k  ,inuc) = con_nucx(k  ,inuc) + xfernuc(k,inuc) * volti(k  ,inuc)

        if (k <= lsw(iw)) then
           con_nucx(k,inuc) = con_nucx(k,inuc) - xfernuc_sfc(k-1,inuc) * volti(k,inuc)
        endif
     enddo
  enddo

  ! Precip Scavenging for nuclei (Saleeby 2011) Slinn (1983) method detailed in
  ! (X.Wang et al. 2010; Seinfeld & Pandis 2006).  Currently assuming spheres
  ! for snow/aggr scavenging.  Literature suggests snow scavenging coefficients
  ! are higher than for rain, so the spherical assumption likely sets a lower
  ! bound for snow scavenging for now.

  ! If there is convective (rain) precipitation, assign a mean raindrop
  ! diameter, fall speed, and precipitation flux profile so that scavenging can
  ! be evaluated in the same manner as for resolved microphysics precipitation.
  ! Enter the information as the cloud water category (lcat = 1), which is
  ! otherwise not used for scavenging.

  ks1(1) = 0
  ks2(1) = 0

  ks1(2:11) = k1(2:11)
  ks2(2:11) = k2(2:11)

  if (conprr(iw) > 1.e-12) then
     ks1(1)  = lpw(iw)
     ks2(1)  = kcutop(iw)
     ks2(11) = max(ks2(11),ks2(1))
     dzicld = 1. / (zm(kcutop(iw)) - zm(kcubot(iw)-1))

     do k = lpw(iw),ks2(1)
        dmb(k,1) = 0.002
        pcpvel(k,1) = 10.

        if (k < kcubot(iw)) then
           pcpfluxr(k,1) = conprr(iw) * dtl0
        else
           pcpfluxr(k,1) = conprr(iw) * dtl0 * (zm(kcutop(iw)) - zt(k)) * dzicld
        endif
     enddo
  endif

  ! No scavenging if there is no precipitation in this IW column

  if (lpw(iw) > ks2(11)) go to 50

  ! Loop over precipitation categories

  do lcat = 1,ncat
     if (lpw(iw) > ks2(lcat)) cycle

     ! Loop over precipation levels for this category

     do k = ks1(lcat),ks2(lcat)

        ! Compute quantities dependent on precipitation but not aerosol nuclei

        re       (k,lcat) = 0.5 * dmb(k,lcat) * pcpvel(k,lcat) * real(rhoa(k)) / dynvisc(k)
        sqrt_re  (k,lcat) = sqrt(re(k,lcat))
        a1r               = log(1. + re(k,lcat))
        sstar    (k,lcat) = (1.2 + 0.083333 * a1r) / (1. + a1r)
        pcpodmb  (k,lcat) = pcpfluxr(k,lcat) * 1.e-3 / dmb(k,lcat)
     enddo
  enddo ! lcat

! Loop over aerosol-nuclei categories

  do inuc = 1,nnuc

     ! Loop over grid levels with possibility of precipitation

     do k = lpw(iw),ks2(11)

        ! Compute quantities dependent on aerosol nuclei but not precipitation

        sqrt_schm(k) = sqrt(schm(k,inuc))
        cbrt_schm(k) = schm(k,inuc)**0.33333333
        sqrt_rhoi_nuc(k) = sqrt(rhow / rho_wetnuc(k,inuc))

        scavfrac(k) = 0.

     enddo

     ! Loop over precipitation categories

     do lcat = 1,ncat
        if (lpw(iw) > ks2(lcat)) cycle

       ! Loop over precipation levels for this category

        do k = ks1(lcat),ks2(lcat)

           st = 2 * tau(k,inuc) * (pcpvel(k,lcat) - vwetnuc(k,inuc)) / dmb(k,lcat)

           dratio = dwetnuc(k,inuc) / dmb(k,lcat)

           collect1 = (4. / (re(k,lcat) * schm(k,inuc))) &
                    * (1.0 + sqrt_re(k,lcat) * (0.4 * cbrt_schm(k) + 0.16 * sqrt_schm(k)))

           collect2 = (4. * dratio) &
                    * (dynvisc(k) / dynvisc_liqwat(k) + (1. + 2. * sqrt_re(k,lcat)) * dratio)

           ecollect = collect1 + collect2

           if (st > sstar(k,lcat)) then
              stcorr = ((st - sstar(k,lcat)) / (st - sstar(k,lcat) + 2. / 3.))**1.5
              ecollect = ecollect + stcorr * sqrt_rhoi_nuc(k)
           endif
           if (ecollect > 1.0) ecollect = 1.0

           scavfrac(k) = scavfrac(k) + 1.5 * ecollect * pcpodmb(k,lcat)

        enddo ! k

     enddo ! lcat

     ! Loop over grid levels with possibility of precipitation

     do k = lpw(iw),ks2(11)

        con_nucx(k,inuc) = con_nucx(k,inuc) - scavfrac(k) * con_nucx(k,inuc)

        if (con_nucx(k,inuc) < 0.) con_nucx(k,inuc) = 0.

     enddo

  enddo ! inuc

  50 continue ! For skipping precipitation scavenging

  ! Copy con_nucx back to separate con_ccnx, con_gccnx, and con_ifnx arrays

  inuc = 0
  if (iccn >= 2) then
     do ic = 1,nccntyp
        inuc = inuc + 1
        con_ccnx(:,ic) = con_nucx(:,inuc)
     enddo
  endif

  if (igccn == 2) then
     inuc = inuc + 1
     con_gccnx(:) = con_nucx(:,inuc)
  endif

  if (iifn == 2) then
     inuc = inuc + 1
     con_ifnx(:) = con_nucx(:,inuc)
  endif

end subroutine nuclei_deposition
