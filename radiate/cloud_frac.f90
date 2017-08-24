subroutine get_cloud_frac(iw, ka, frac)

  use mem_grid,    only: mza, zm, zt, glatw, glonw, dzt, dzim
  use misc_coms,   only: icfrac, cfracrh1, cfracrh2, cfraccup
  use consts_coms, only: t00
  use mem_basic,   only: rho, tair, sh_v, sh_w
  use misc_coms,   only: nqparm
  use mem_ijtabs,  only: itab_w
  use mem_cuparm,  only: iactcu, kcubot, kcutop, qwcon, conprr
  use mem_turb,    only: frac_land
  use mem_micro,   only: sh_c, sh_p
  use clouds_gno,  only: cu_cldfrac
  use therm_lib,   only: rhovsl, rhovsi

  implicit none

  integer, intent(in)    :: iw
  integer, intent(in)    :: ka
  real   , intent(inout) :: frac(mza)

  real    :: rhov (mza)       ! vapor density [kg_vap/m^3]
  real    :: rhoc (mza)       ! bulk cloud water density [kg_cld/m^3]
  real    :: rhop (mza)       ! bulk pristine ice density [kg_pris/m^3]
  real    :: tc, rh
  real    :: rhl, rhi, fracl, fraci
  real    :: fland
  integer :: k, mrlw
  real    :: rh00
  real    :: abslat, wt20, wt60, cfrh1, cfrh2, dcfrhi
  real    :: rh_tot(mza)
  real    :: qc_sub(mza)
  real    :: cu_cldf(mza)
  real    :: qsat, qw

  if (allocated(frac_land)) then
     fland = frac_land(iw)
  else
     fland = 1.0
  endif

  do k = ka, mza

     rhov(k) = max(0.,sh_v(k,iw)) * rho(k,iw)

     if (allocated(sh_c)) then
        rhoc(k) = max(0.,sh_c(k,iw)) * rho(k,iw)
     else
        rhoc(k) = 0.
     endif

     if (allocated(sh_p)) then
        rhop(k) = max(0.,sh_p(k,iw)) * rho(k,iw)
     else
        rhop(k) = 0.
     endif

  enddo

  if (icfrac == 1) then

! Fractional cloudiness based on RH from Mocko and Cotton (1995).

     if (fland > 0.5) then
        rh00 = 0.85
     else
        rh00 = 0.75
     endif

     do k = ka, mza
        tc = tair(k,iw) - t00
        if (tc > -10.0) then
           rh = rhov(k) / rhovsl(tc)
        else
           rh = rhov(k) / rhovsi(tc)
        endif
        rh = min(rh, 1.0)
        frac(k) = max( 1.0 - sqrt(( 1.0 - rh ) / ( 1.0 - rh00)), 0.0)
     enddo

! Walko's linear forms with inclusion of cloud and pristine ice condensate... 

  else

! Use adjustable lower and upward RH thresholds from namelist

     if (icfrac == 2) then

        cfrh1 = cfracrh1
        cfrh2 = cfracrh2

     else

! Latitudinal and land/sea variation of cloud fraction parameters

        abslat = abs(glatw(iw))

        if     (abslat < 20.) then
           wt60 = 0.
        elseif (abslat > 60.) then
           wt60 = 1.
        else
           wt60 = (abslat - 20.) / 40.
        endif

        wt20 = 1. - wt60

! Select set of parameters with namelist flag icfrac

        if (icfrac == 3) then

           if (fland > 0.5) then
              cfrh1 = wt20 * 0.90 + wt60 * 0.90  ! land set 1
              cfrh2 = wt20 * 1.40 + wt60 * 1.40  ! land set 1
           else
              cfrh1 = wt20 * 1.00 + wt60 * 0.95  ! sea set 1
              cfrh2 = wt20 * 1.20 + wt60 * 1.20  ! sea set 1
           endif

        elseif (icfrac == 4) then

           if (fland > 0.5) then
              cfrh1 = wt20 * 0.80 + wt60 * 0.90  ! land set 2
              cfrh2 = wt20 * 1.05 + wt60 * 1.40  ! land set 2
           else
              cfrh1 = wt20 * 1.00 + wt60 * 0.95  ! sea set 2
              cfrh2 = wt20 * 1.20 + wt60 * 1.20  ! sea set 2
           endif

        elseif (icfrac == 5) then

           if (fland > 0.5) then
              cfrh1 = wt20 * 0.85 + wt60 * 0.90  ! land set 3
              cfrh2 = wt20 * 1.00 + wt60 * 1.40  ! land set 3
           else
              cfrh1 = wt20 * 1.00 + wt60 * 0.95  ! sea set 3
              cfrh2 = wt20 * 1.20 + wt60 * 1.20  ! sea set 3
           endif

        elseif (icfrac == 6) then

           if (fland > 0.5) then
              cfrh1 = wt20 * 0.80 + wt60 * 0.90  ! land set 4
              cfrh2 = wt20 * 1.00 + wt60 * 1.40  ! land set 4
           else
              cfrh1 = wt20 * 1.00 + wt60 * 0.95  ! sea set 4
              cfrh2 = wt20 * 1.20 + wt60 * 1.20  ! sea set 4
           endif

        endif

     endif

     dcfrhi = 1. / max(1.e-6, cfrh2-cfrh1)

     do k = ka, mza
        tc = tair(k,iw) - t00
        rhl = (rhov(k) + rhoc(k) + rhop(k)) / rhovsl(tc)
        rhi = (rhov(k) +           rhop(k)) / rhovsi(tc)

        fracl = (rhl - cfrh1) * dcfrhi
        fraci = (rhi - cfrh1) * dcfrhi

        frac(k) = min(1.0, max(0.0, fracl, fraci))
     enddo

  endif

  ! If there is subgrid convection, modify the estimated cloud fraction to include
  ! the convective clouds from the cumulus scheme

  if (iactcu(iw) == 1) then

     ! This section estimates the cloud fraction from subgrid cumulus and 
     ! any resolved clouds based on a lookup table of the scheme of 
     ! Bony and Emanuel (2001, JAS)

     do k = kcubot(iw), kcutop(iw)
        tc = tair(k,iw) - t00
        if (tc > -10.0) then
           qsat = rhovsl(tc) / rho(k,iw)
        else
           qsat = rhovsi(tc) / rho(k,iw)
        endif

        qw = max( sh_w(k,iw), 1.e-8 )

        rh_tot(k) = qw / qsat
        qc_sub(k) = sqrt( qwcon(k,iw) / qw )
     enddo

     call cu_cldfrac(kcubot(iw), kcutop(iw), rh_tot, qc_sub, cu_cldf)

     ! Do we want to overwrite the resolved cloud fraction or merge the two?
     do k = kcubot(iw), kcutop(iw)
        frac(k) = max( min(cu_cldf(k), 1.0), 0.0 )
     enddo

     if (conprr(iw) > 1.e-12) then

        ! If there is deep (precipitating) convection, 
        ! include subgrid clouds below convective cloud base

        do k = ka, kcubot(iw) -1
           frac(k) = max(frac(k), frac(kcubot(iw)))
        enddo

        ! If there is deep convection, limit the resolved cloud fraction 
        ! below and just above the cumulus to create some breaks

        do k = ka, kcubot(iw) - 1
           frac(k) = min(frac(k), cfraccup)
        enddo
        frac(kcutop(iw)+1) = min(frac(k), cfraccup)

     endif

  endif

end subroutine get_cloud_frac
