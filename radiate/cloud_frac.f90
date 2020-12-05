subroutine calc_3d_cloud_fraction(mrl)

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_grid,    only: mza, lpw
  use mem_radiate, only: cloud_frac
  use mem_cuparm,  only: iactcu, qwcon
  use mem_basic,   only: rr_w, rr_v
  use misc_coms,   only: icfrac

  implicit none

  integer, intent(in) :: mrl
  integer             :: j, iw, ka, k
  real                :: cond
  real                :: frac(mza)

  real,     parameter :: cond_min = 1.e-8
  real,     parameter :: frac_min = 0.10

  if (mrl == 0) return

  !$omp parallel private(frac)
  !$omp do private(iw, ka, k, cond) schedule(guided)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     ka = lpw(iw)

     if (icfrac > 0) then
        call get_cloud_frac(iw,ka,frac)
     else
        frac(ka:mza) = 1.0
     endif

     do k = ka, mza

        if (iactcu(iw) > 0) then
           cond = rr_w(k,iw) - rr_v(k,iw) + qwcon(k,iw)
        else
           cond = rr_w(k,iw) - rr_v(k,iw)
        endif

        if (cond > cond_min) then
           cloud_frac(k,iw) = max(frac(k), frac_min)
        else
           cloud_frac(k,iw) = 0.0
        endif

     enddo

  enddo
  !$omp end do
  !$omp end parallel

end subroutine calc_3d_cloud_fraction




subroutine get_cloud_frac(iw, ka, frac)

  use mem_grid,    only: mza, glatw
  use misc_coms,   only: icfrac, cfracrh1, cfracrh2, cfraccup
  use consts_coms, only: t00
  use mem_basic,   only: rho, tair, rr_v, rr_w
  use mem_cuparm,  only: iactcu, kcubot, kcutop, qwcon, conprr
  use mem_turb,    only: frac_land
  use mem_micro,   only: rr_c, rr_p
  use clouds_gno,  only: cu_cldfrac
  use therm_lib,   only: rhovsl_inv, rhovsi_inv

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
  integer :: k
  real    :: rh00
  real    :: abslat, wt20, wt60, cfrh1, cfrh2, dcfrhi
  real    :: rh_tot(mza)
  real    :: qc_sub(mza)
  real    :: cu_cldf(mza)
  real    :: qw, rhow

  if (allocated(frac_land)) then
     fland = frac_land(iw)
  else
     fland = 1.0
  endif

  do k = ka, mza

     rhov(k) = max(0.,rr_v(k,iw)) * real(rho(k,iw))

     if (allocated(rr_c)) then
        rhoc(k) = max(0.,rr_c(k,iw)) * real(rho(k,iw))
     else
        rhoc(k) = 0.
     endif

     if (allocated(rr_p)) then
        rhop(k) = max(0.,rr_p(k,iw)) * real(rho(k,iw))
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
        rh = rhov(k) * rhovsl_inv(tc)

        if (rhop(k) > 1.e-20) then
           rhi = rhov(k) * rhovsi_inv(tc)
           rh  = (rh * rhoc(k) + rhi * rhop(k)) / (rhoc(k) + rhop(k))
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
        rhl = (rhov(k) + rhoc(k) + rhop(k)) * rhovsl_inv(tc)
        rhi = (rhov(k) +           rhop(k)) * rhovsi_inv(tc)

        fracl = (rhl - cfrh1) * dcfrhi
        fraci = (rhi - cfrh1) * dcfrhi

        frac(k) = min(1.0, max(0.0, fracl, fraci))
     enddo

  endif

  ! If there is subgrid convection, modify the estimated cloud fraction to include
  ! the convective clouds from the cumulus scheme

  if (iactcu(iw) > 1) then

     ! This section estimates the cloud fraction from subgrid cumulus and
     ! any resolved clouds based on a lookup table of the scheme of
     ! Bony and Emanuel (2001, JAS)

     do k = kcubot(iw), kcutop(iw)

        tc   = tair(k,iw) - t00
        qw   = max( rr_w(k,iw), 1.e-8 )
        rhow = qw * real(rho(k,iw))

        if (tc > -10.0) then
           rh_tot(k) = rhow * rhovsl_inv(tc)
        else
           rh_tot(k) = rhow * rhovsi_inv(tc)
        endif

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
