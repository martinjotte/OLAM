subroutine timeseries_plots()

  use mem_basic,  only: vze
  use mem_micro,  only: sh_c, sh_d, sh_r, sh_p, sh_s, sh_a, sh_g, sh_h, &
                        accpd, accpr, accpp, accps, accpa, accpg, accph
  use mem_ijtabs, only: jtab_w, jtw_prog
  use misc_coms,  only: io6, time_istp8, timmax8, dtlong, runtype, time8
  use oplot_coms, only: op
  use oname_coms, only: nl
  use mem_grid,   only: zt
  use mem_plot,   only: latheat_liq_accum_prev0, latheat_liq_accum_prev1, &
                        latheat_ice_accum_prev0, latheat_ice_accum_prev1, &
                        time8_prev0, time8_prev1

  integer, save :: ncall = 0, ncall_tot
  integer :: ifill = 1
  integer, parameter :: ny = 64  ! Number of vertical levels to plot

  integer :: iw, j, kh, k, nc

  real, save, allocatable :: ge1(:), vctr18(:)
  real, save, allocatable :: he1(:,:), he2(:,:), he3(:,:), he4(:,:), &
                             he5(:,:), he6(:,:), he7(:,:), he8(:,:)
  real, save, allocatable :: height(:), enw(:)

  real, save :: aspect = .7
  real, save :: scalelab = .014
  real, save :: timebeg,timeend,timedif,timeinc
  real :: enprog, lwc, pcp

  ncall = ncall + 1

  ! On the first call to this subroutine, compute the plotting time increment
  ! and allocate arrays

  if (ncall == 1) then

     if (runtype == 'PLOTONLY') then
        ncall_tot = nl%nplt_files
     else
        ncall_tot = int(timmax8 / dtlong) + 2
     endif

     allocate (ge1(ncall_tot), vctr18(ncall_tot))
     allocate (he1(ny,ncall_tot), he2(ny,ncall_tot), he3(ny,ncall_tot), &
               he4(ny,ncall_tot), he5(ny,ncall_tot), he6(ny,ncall_tot), &
               he7(ny,ncall_tot), he8(ny,ncall_tot))

     allocate (enw(ny), height(ny))
     height(1:ny) = zt(2:ny+1) * 1.e-3

     ge1(:)   = 0.
     he1(:,:) = 0.
     he2(:,:) = 0.
     he3(:,:) = 0.
     he4(:,:) = 0.
     he5(:,:) = 0.
     he6(:,:) = 0.
     he7(:,:) = 0.
     he8(:,:) = 0.

 !   timebeg = real(time_istp8)   / 3600.
     timebeg = 0.
     timeend = real(timmax8) / 3600.
     timedif = timeend - timebeg

     if (timedif < .03) then
        timeinc = .001
     elseif (timedif < .06) then
        timeinc = .002
     elseif (timedif < .1) then
        timeinc = .004

     elseif (timedif < .3) then
        timeinc = .01
     elseif (timedif < .6) then
        timeinc = .02
     elseif (timedif < 1.) then
        timeinc = .04

     elseif (timedif < 3.) then
        timeinc = .1
     elseif (timedif < 6.) then
        timeinc = .2
     elseif (timedif < 10.) then
 !      timeinc = .4
        timeinc = .2

     elseif (timedif < 30.) then
        timeinc = 1.
     elseif (timedif < 60.) then
        timeinc = 2.
     elseif (timedif < 100.) then
        timeinc = 4.

     elseif (timedif < 300.) then
        timeinc = 10.
    elseif (timedif < 600.) then
        timeinc = 20.
     elseif (timedif < 1000.) then
        timeinc = 40.
     endif
      
  endif

  vctr18(ncall) = real(time_istp8) / 3600.

  enprog = 0.
  enw(:) = 0.

  ! Horizontal loop over all IW points

  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

    enprog = enprog + 1.

    ge1(ncall) = ge1(ncall) + accpd(iw) + accpr(iw) + accpp(iw) + accps(iw) &
                            + accpa(iw) + accpg(iw) + accph(iw)

    do kh = 1,ny
       k = kh + 1

       if (he1(kh,ncall) < vze(k,iw)) he1(kh,ncall) = vze(k,iw)

       if (vze(k,iw) > 0.) then
          enw(kh)       = enw(kh) + 1.
          he2(kh,ncall) = he2(kh,ncall) + vze(k,iw)**2
       endif

       lwc = 1000. * (sh_c(k,iw) + sh_d(k,iw))

       if (he3(kh,ncall) < lwc) &
           he3(kh,ncall) = lwc

       he4(kh,ncall) = he4(kh,ncall) + lwc**2

       pcp = 1000. * (sh_r(k,iw) + sh_h(k,iw))

       if (he5(kh,ncall) < pcp) &
           he5(kh,ncall) = pcp

       he6(kh,ncall) = he6(kh,ncall) + pcp**2

       he7(kh,ncall) = he7(kh,ncall) &
                     + latheat_liq_accum_prev0(k,iw) - latheat_liq_accum_prev1(k,iw)

       he8(kh,ncall) = he8(kh,ncall) &
                     + latheat_ice_accum_prev0(k,iw) - latheat_ice_accum_prev1(k,iw)
    enddo

  enddo

  ge1(ncall) = ge1(ncall) / enprog

  do kh = 1,ny
     he2(kh,ncall) = sqrt(he2(kh,ncall) / max(1.,enw(kh)))
     he4(kh,ncall) = sqrt(he4(kh,ncall) / enprog)
     he6(kh,ncall) = sqrt(he6(kh,ncall) / enprog)
     he7(kh,ncall) = he7(kh,ncall) * 86400. / ((time8_prev0 - time8_prev1) * enprog)
     he8(kh,ncall) = he8(kh,ncall) * 86400. / ((time8_prev0 - time8_prev1) * enprog)
  enddo

  if (real(time_istp8) + .5 * real(dtlong) < real(timmax8)) return

  ! Reopen the current graphics output workstation if it is closed

  call o_reopnwk()

  ! Plot time series

  !-------------------------------------------------------------------
  call plotback()
  call oplot_xy2('0','N','N',aspect,scalelab,10,0,               &
                 ncall_tot,  vctr18,ge1,                         &
                 'time(hours)','accumulated precipitation (mm)', &
                 timebeg,timeend,timeinc,5,                      &
                 -0.1,1.0,.02,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,              &
                  'max W ',' (m/s)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he1, 77,ifill, &
                  timebeg,timeend,timeinc,5,                &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,              &
                  'rms +W ',' (m/s)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he2,443,ifill, &
                  timebeg,timeend,timeinc,5,                &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,                 &
                  'max LWC ',' (g/kg)', ncall_tot,ny,vctr18,height, & 
                  'time(hours)','height (km)',he3,418,ifill,    &
                  timebeg,timeend,timeinc,5,                   &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,                 &
                  'rms LWC ',' (g/kg)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he4,443,ifill,    &
                  timebeg,timeend,timeinc,5,                   &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,                 &
                  'max PCP ',' (g/kg)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he5,418,ifill,    &
                  timebeg,timeend,timeinc,5,                   &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,                 &
                  'rms PCP ',' (g/kg)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he6,443,ifill,    &
                  timebeg,timeend,timeinc,5,                   &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,                 &
                  'latent_heating_liq ',' (K/day)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he7,445,ifill,    &
                  timebeg,timeend,timeinc,5,                   &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------
  call plotback()
  call oplot_zxy2('0','N','c',aspect,scalelab,                 &
                  'latent_heating_ice ',' (K/day)', ncall_tot,ny,vctr18,height, &
                  'time(hours)','height (km)',he8,445,ifill,    &
                  timebeg,timeend,timeinc,5,                   &
                  0.0,height(ny),1.0,5  )
  call o_frame()
  !-------------------------------------------------------------------

  ! Close the current workstation if output
  ! is to a NCAR graphics meta file. This allows viewing the complete
  ! meta file (including the last frame) during a run and in case the
  ! simulation crashes.

  if ((trim(runtype) /= 'PLOTONLY') .and. (op%plttype == 0)) then
     call o_clswk()
  endif

end subroutine timeseries_plots
