subroutine timeseries_plots(type)

  use mem_basic,  only: wc, vxe, vye, vze, tair, press, rho, rr_v
  use mem_micro,  only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, &
                        accpd, accpr, accpp, accps, accpa, accpg, accph
  use mem_ijtabs, only: jtab_w, jtw_prog
  use misc_coms,  only: io6, time_istp8, timmax8, dtlong, runtype, time8
  use oplot_coms, only: op
  use oname_coms, only: nl
  use mem_grid,   only: zt, dzt, xew, yew, zew, lpw
  use mem_plot,   only: latheat_liq_accum_prev0, latheat_liq_accum_prev1, &
                        latheat_ice_accum_prev0, latheat_ice_accum_prev1, &
                        time8_prev0, time8_prev1
  use consts_coms, only: cp, alvlocp, r8
  use therm_lib,  only: rhovsl

  implicit none

  character(1), intent(in) :: type

  integer, save :: ncall = 0, ncall_tot, numtimes, numsims = 1
  integer :: ifill = 1
  integer, parameter :: ny = 38  ! Number of vertical levels to plot (Harvey: k = 38 is at 19.8 km)

  integer, parameter :: nrad_lh = 50

  integer, save :: nlh(nrad_lh)

  integer :: iw, j, kh, k, nc, ka, jtime, jsim, jft, irad, isim, irad_lh, icolor

  integer :: khe, ihe

  integer, save, allocatable :: iw0(:,:)

  real, save, allocatable :: fe1(:,:), fe2(:,:), fe3(:,:)
  real, save, allocatable :: ge1(:,:), ge2(:,:), ge3(:,:), ge4(:,:), ge5(:,:), ge6(:,:)
  real, save, allocatable :: he1(:,:,:),  he2(:,:,:), he3(:,:,:), he4(:,:,:), &
                             he5(:,:,:),  he6(:,:,:), he7(:,:,:), he8(:,:,:), &
                             he9(:,:,:), he10(:,:,:), he11(:,:,:), he12(:,:,:), &
                             he13(:,:,:),he14(:,:,:), he15(:,:,:)
  real, save, allocatable :: ae6(:), vctr17(:), vctr18(:,:), vctr19(:)
  real, save, allocatable :: height(:), enw(:), ensat(:)

  real, save :: aspect = .7
  real, save :: scalelab = .014
  real, save :: timebeg,timeend,timedif,timeinc
  real :: ennow, enall, radnow, radall, lwc, pcp, vels2, velh, wt
  real :: pmsl, ss_liq, ss_thetadif

  character(len=1) :: ipanel
  character(len=60) :: xlab, ylab

  RETURN  ! Comment out this line to execute this subroutine

  ncall = ncall + 1

  ! On the first call to this subroutine, compute the plotting time increment
  ! and allocate arrays

  if (ncall == 1) then

     if (runtype == 'PLOTONLY' .and. type == 'H') then
        ncall_tot = nl%nplt_files
        numsims = 3
     elseif (runtype == 'PLOTONLY' .and. type == 'L') then  ! For this option, need to comment out
        ncall_tot = nl%nlite_files                          ! things that are not on lite files 
        numsims = 1
     else
        ncall_tot = int(timmax8 / dtlong) + 2
        numsims = 1
     endif

     numtimes = ncall_tot / numsims

     if (numtimes * numsims /= ncall_tot) then
        print*, 'numtimes, numsims, ncall_tot ',numtimes, numsims, ncall_tot
        stop 'stopping: numtimes'
     endif

     allocate (fe1(nrad_lh,numsims),fe2(nrad_lh,numsims),fe3(nrad_lh,numsims))
     allocate (ge1(numtimes,numsims), ge2(numtimes,numsims), ge3(numtimes,numsims), ge4(numtimes,numsims), &
               ge5(numtimes,numsims), ge6(numtimes,numsims), vctr18(numtimes,numsims))
     allocate (he1 (ny,numtimes,numsims), he2(ny,numtimes,numsims), he3(ny,numtimes,numsims), &
               he4 (ny,numtimes,numsims), he5(ny,numtimes,numsims), he6(ny,numtimes,numsims), &
               he7 (ny,numtimes,numsims), he8(ny,numtimes,numsims), he9(ny,numtimes,numsims), &
               he10(ny,numtimes,numsims), he11(ny,numtimes,numsims), he12(ny,numtimes,numsims), &
               he13(ny,numtimes,numsims), he14(ny,numtimes,numsims), he15(ny,numtimes,numsims))
     allocate (ae6(17), vctr17(17), vctr19(nrad_lh))
     allocate (iw0(49,3))

     allocate (enw(ny), height(ny), ensat(ny))
     height(1:ny) = zt(2:ny+1) * 1.e-3

     ge1 (:,:)   = 0.
     ge2 (:,:)   = 0.
     ge3 (:,:)   = 0.
     ge4 (:,:)   = 0.
     ge5 (:,:)   = 0.
     ge6 (:,:)   = 2000.

     he1 (:,:,:) = 0.
     he2 (:,:,:) = 0.
     he3 (:,:,:) = 0.
     he4 (:,:,:) = 0.
     he5 (:,:,:) = 0.
     he6 (:,:,:) = 0.
     he7 (:,:,:) = 0.
     he8 (:,:,:) = 0.
     he9 (:,:,:) = 0.
     he10(:,:,:) = 0.
     he11(:,:,:) = -100.
     he12(:,:,:) = 0.
     he13(:,:,:) = 0.
     he14(:,:,:) = 0.
     he15(:,:,:) = 0.

     fe1(:,:) = 0.
     fe2(:,:) = 0.
     fe3(:,:) = 0.

     ae6(:) = 0.

 !   timebeg = real(time_istp8)   / 3600.
     timebeg = 0.
     timeend = real(timmax8) / 3600.

     timeend = 48.  ! (hours)  [For incomplete simulation]

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
      
     ! iw0 is location of minimum MSLP in N, NA, and NA10X simulations

     iw0( 1,1) = 197494
     iw0( 2,1) = 440589
     iw0( 3,1) = 363891
     iw0( 4,1) = 177718
     iw0( 5,1) = 226159
     iw0( 6,1) = 364690
     iw0( 7,1) = 164705
     iw0( 8,1) = 454718
     iw0( 9,1) = 364600
     iw0(10,1) = 215274
     iw0(11,1) = 365545
     iw0(12,1) = 297927
     iw0(13,1) = 391030
     iw0(14,1) = 297954
     iw0(15,1) = 292944
     iw0(16,1) = 172457
     iw0(17,1) = 229390
     iw0(18,1) = 443344
     iw0(19,1) = 318441
     iw0(20,1) = 266537
     iw0(21,1) = 460675
     iw0(22,1) = 335921
     iw0(23,1) = 377728
     iw0(24,1) = 293593
     iw0(25,1) = 377774
     iw0(26,1) = 378736
     iw0(27,1) = 378737
     iw0(28,1) = 125517
     iw0(29,1) = 326627
     iw0(30,1) = 300882
     iw0(31,1) = 402487
     iw0(32,1) = 337951
     iw0(33,1) = 446231
     iw0(34,1) = 324407
     iw0(35,1) = 348380
     iw0(36,1) = 402443
     iw0(37,1) = 337928
     iw0(38,1) = 405162
     iw0(39,1) = 428283
     iw0(40,1) = 405186
     iw0(41,1) = 322756
     iw0(42,1) = 431027
     iw0(43,1) = 370302
     iw0(44,1) = 458052
     iw0(45,1) = 407729
     iw0(46,1) = 407723
     iw0(47,1) = 459776
     iw0(48,1) = 408095
     iw0(49,1) = 330930
     iw0( 1,2) = 197494
     iw0( 2,2) = 332316
     iw0( 3,2) = 288972
     iw0( 4,2) = 385580
     iw0( 5,2) = 297030
     iw0( 6,2) = 325339
     iw0( 7,2) = 441326
     iw0( 8,2) = 418613
     iw0( 9,2) = 296965
     iw0(10,2) = 226087
     iw0(11,2) = 226090
     iw0(12,2) = 417580
     iw0(13,2) = 278857
     iw0(14,2) = 429436
     iw0(15,2) = 289551
     iw0(16,2) = 198514
     iw0(17,2) = 437642
     iw0(18,2) = 281154
     iw0(19,2) = 208920
     iw0(20,2) = 305727
     iw0(21,2) = 367470
     iw0(22,2) = 396765
     iw0(23,2) = 427958
     iw0(24,2) = 437962
     iw0(25,2) = 269923
     iw0(26,2) = 428073
     iw0(27,2) = 423906
     iw0(28,2) = 293923
     iw0(29,2) = 326624
     iw0(30,2) = 399589
     iw0(31,2) = 269365
     iw0(32,2) = 402158
     iw0(33,2) = 282671
     iw0(34,2) = 421688
     iw0(35,2) = 419150
     iw0(36,2) = 330336
     iw0(37,2) = 355932
     iw0(38,2) = 405266
     iw0(39,2) = 432789
     iw0(40,2) = 417206
     iw0(41,2) = 260511
     iw0(42,2) = 405325
     iw0(43,2) = 168237
     iw0(44,2) = 404907
     iw0(45,2) = 277588
     iw0(46,2) = 407731
     iw0(47,2) = 237282
     iw0(48,2) = 371350
     iw0(49,2) = 357848
     iw0( 1,3) = 197494
     iw0( 2,3) = 385602
     iw0( 3,3) = 207373
     iw0( 4,3) = 350353
     iw0( 5,3) = 429099
     iw0( 6,3) = 433618
     iw0( 7,3) = 374792
     iw0( 8,3) = 436986
     iw0( 9,3) = 226096
     iw0(10,3) = 124354
     iw0(11,3) = 333043
     iw0(12,3) = 460451
     iw0(13,3) = 278857
     iw0(14,3) = 108075
     iw0(15,3) = 257802
     iw0(16,3) = 307075
     iw0(17,3) = 319586
     iw0(18,3) = 276350
     iw0(19,3) = 393855
     iw0(20,3) = 298846
     iw0(21,3) = 299828
     iw0(22,3) = 299832
     iw0(23,3) = 453909
     iw0(24,3) = 377735
     iw0(25,3) = 418942
     iw0(26,3) = 377782
     iw0(27,3) = 270217
     iw0(28,3) = 348068
     iw0(29,3) = 232688
     iw0(30,3) = 232668
     iw0(31,3) = 445303
     iw0(32,3) = 402560
     iw0(33,3) = 234324
     iw0(34,3) = 234322
     iw0(35,3) = 348377
     iw0(36,3) = 294227
     iw0(37,3) = 402341
     iw0(38,3) = 180910
     iw0(39,3) = 405265
     iw0(40,3) = 380610
     iw0(41,3) = 317963
     iw0(42,3) = 438915
     iw0(43,3) = 405327
     iw0(44,3) = 405296
     iw0(45,3) = 404902
     iw0(46,3) = 235716
     iw0(47,3) = 357726
     iw0(48,3) = 416900
     iw0(49,3) = 408096

     do irad = 1,nrad_lh
        vctr19(irad) = 10.0 * real(irad) - 5.0
     enddo
 
  endif

  jtime = mod(ncall-1,numtimes) + 1

  jsim = (ncall - 1) / numtimes + 1

  jft = 25 + 10 * (jtime - 1)

  write (6,'(a,9i6)') 'numsims, numtimes, ncall_tot, ncall, jtime, jsim ', &
                       numsims, numtimes, ncall_tot, ncall, jtime, jsim

  vctr18(jtime,jsim) = real(time_istp8) / 3600.

  ennow = 0.
  enall = 0.
  enw(:) = 0.
  ensat(:) = 0.

  nlh(:) = 0

  ! Horizontal loop over all IW points

  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

    radnow = sqrt((xew(iw)-xew(iw0(jft,jsim)))**2 + (yew(iw)-yew(iw0(jft,jsim)))**2 + (zew(iw)-zew(iw0(jft,jsim)))**2)
    radall = sqrt((xew(iw)-xew(396574))**2 + (yew(iw)-yew(396574))**2 + (zew(iw)-zew(396574))**2)

    if (radall <= 800.e3) then
       enall = enall + 1.
       ge1(jtime,jsim) = ge1(jtime,jsim) + accpd(iw) + accpr(iw) + accpp(iw) + accps(iw) &
                       + accpa(iw) + accpg(iw) + accph(iw)
    endif

    if (radnow <= 500.e3) then
       irad_lh = int(radnow / 10.e3) + 1

       nlh(irad_lh) = nlh(irad_lh) + 1

       do kh = 1,ny
          k = kh + 1

          if (jtime /= 1) then
             fe1(irad_lh,jsim) = fe1(irad_lh,jsim) &
                + cp * rho(k,iw) * dzt(k) * (latheat_liq_accum_prev0(k,iw) - latheat_liq_accum_prev1(k,iw)) &
                / max(1._r8,time8_prev0 - time8_prev1)
             fe2(irad_lh,jsim) = fe2(irad_lh,jsim) &
                + cp * rho(k,iw) * dzt(k) * (latheat_ice_accum_prev0(k,iw) - latheat_ice_accum_prev1(k,iw)) &
                / max(1._r8,time8_prev0 - time8_prev1)
          endif
       enddo

    endif


    if (radnow <= 100.e3) then
       ennow = ennow + 1.

       ka = lpw(iw)

       pmsl = press(ka,iw) * (1. - .0065 * zt(ka) / (tair(ka,iw) + .0065 * zt(ka)))**(-5.257) * .01  ! hydrostatic eqn.

       if (iw == iw0(jft,jsim)) then
          print*, ' '
          print*, 'ncall, pmsl ',ncall,pmsl,iw,ka
          print*, press(ka,iw),zt(ka),tair(ka,iw),pmsl
          print*, ge6(jtime,jsim)
          print*, ' '
       endif


       if (ge6(jtime,jsim) > pmsl) ge6(jtime,jsim) = pmsl

       do kh = 1,ny
          k = kh + 1

          lwc = 1000. * (rr_c(k,iw) + rr_d(k,iw))
          pcp = 1000. * (rr_r(k,iw) + rr_h(k,iw))
          ss_liq = (rr_v(k,iw) * rho(k,iw) / rhovsl(tair(k,iw)-273.15) - 1.0) * 1.e2

          ss_thetadif = min(0.,(rhovsl(tair(k,iw)-273.15) - rr_v(k,iw) * real(rho(k,iw)))) &
                      * alvlocp / real(rho(k,iw))

          wt = 0.5 * (wc(k,iw) + wc(k-1,iw))

          vels2 = vxe(k,iw)**2 + vye(k,iw)**2 + vze(k,iw)**2

          velh = sqrt(vels2 - wt**2)

          if (ge2(jtime,jsim) < lwc .and. tair(k,iw) < 273.15) &
              ge2(jtime,jsim) = lwc

          if (ge3(jtime,jsim) < wt) ge3(jtime,jsim) = wt

          ! Max VELH below 100 m

          if (zt(k) < 100.) then
             if (ge4(jtime,jsim) < velh) ge4(jtime,jsim) = velh
          endif

          if (ge5(jtime,jsim) < lwc) ge5(jtime,jsim) = lwc

          if (he1(kh,jtime,jsim) < wt) he1(kh,jtime,jsim) = wt

          if (wt > 0.) then
             enw(kh)       = enw(kh) + 1.
             he2(kh,jtime,jsim) = he2(kh,jtime,jsim) + wt**2
          endif


          if (he3(kh,jtime,jsim) < lwc) &
              he3(kh,jtime,jsim) = lwc

          if (he9(kh,jtime,jsim) < lwc .and. tair(k,iw) < 273.15) &
              he9(kh,jtime,jsim) = lwc

          if (he11(kh,jtime,jsim) < ss_liq) &
              he11(kh,jtime,jsim) = ss_liq

          if (he15(kh,jtime,jsim) > ss_thetadif) &
              he15(kh,jtime,jsim) = ss_thetadif

!D          if (ss_liq > 0.) then
             ensat(kh) = ensat(kh) + 1.

             if (ss_liq > 5.) then
                he12(kh,jtime,jsim) = he12(kh,jtime,jsim) + 1.
             endif

             if (ss_liq > 10.) then
                he13(kh,jtime,jsim) = he13(kh,jtime,jsim) + 1.
             endif

             if (ss_liq > 30.) then
                he14(kh,jtime,jsim) = he14(kh,jtime,jsim) + 1.
             endif
!D          endif

          he4(kh,jtime,jsim) = he4(kh,jtime,jsim) + lwc**2

          if (he5(kh,jtime,jsim) < pcp) &
              he5(kh,jtime,jsim) = pcp

          he6(kh,jtime,jsim) = he6(kh,jtime,jsim) + pcp**2

          if (jtime == 1) then
             he7(kh,jtime,jsim) = he7(kh,jtime,jsim) &
                           + latheat_liq_accum_prev0(k,iw)

             he8(kh,jtime,jsim) = he8(kh,jtime,jsim) &
                           + latheat_ice_accum_prev0(k,iw)

          else

             he7(kh,jtime,jsim) = he7(kh,jtime,jsim) &
                           + latheat_liq_accum_prev0(k,iw) - latheat_liq_accum_prev1(k,iw)

             he8(kh,jtime,jsim) = he8(kh,jtime,jsim) &
                           + latheat_ice_accum_prev0(k,iw) - latheat_ice_accum_prev1(k,iw)

          endif

       enddo

    endif

  enddo

  ge1(jtime,jsim) = ge1(jtime,jsim) / enall

  do kh = 1,ny
     he2 (kh,jtime,jsim) = sqrt(he2(kh,jtime,jsim) / max(1.,enw(kh)))
     he4 (kh,jtime,jsim) = sqrt(he4(kh,jtime,jsim) / ennow)
     he6 (kh,jtime,jsim) = sqrt(he6(kh,jtime,jsim) / ennow)
     he7 (kh,jtime,jsim) = he7 (kh,jtime,jsim) * 86400. / (max(1._r8,time8_prev0 - time8_prev1) * ennow)
     he8 (kh,jtime,jsim) = he8 (kh,jtime,jsim) * 86400. / (max(1._r8,time8_prev0 - time8_prev1) * ennow)
     he12(kh,jtime,jsim) = he12(kh,jtime,jsim) * 100. / max(1.,ensat(kh))
     he13(kh,jtime,jsim) = he13(kh,jtime,jsim) * 100. / max(1.,ensat(kh))
     he14(kh,jtime,jsim) = he14(kh,jtime,jsim) * 100. / max(1.,ensat(kh))
  enddo

  if (ncall < ncall_tot) return

  do irad_lh = 1,nrad_lh
     do jsim = 1,numsims
        fe1(irad_lh,jsim) = fe1(irad_lh,jsim) / ((numtimes-1) * nlh(irad_lh))
        fe2(irad_lh,jsim) = fe2(irad_lh,jsim) / ((numtimes-1) * nlh(irad_lh))
        fe3(irad_lh,jsim) = fe1(irad_lh,jsim) + fe2(irad_lh,jsim)

write(6,'(a,4i12,3e12.3)') 'fes ',irad_lh,jsim,nlh(irad_lh),numtimes, &
   fe1(irad_lh,jsim),fe2(irad_lh,jsim),fe3(irad_lh,jsim)

     enddo
  enddo

  ! Reopen the current graphics output workstation if it is closed

  call o_reopnwk()

  ! Plot time series

  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'radius (km)' ; ylab = 'latent heating liquid (W/m^2)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''            ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''            ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''            ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0, &
                    nrad_lh,  vctr19, &
                    fe1(: ,isim), &
                    xlab,ylab, &
                    0.0,500.0,10.0,5, &
                    -1500.0,10000.0,400.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'radius (km)' ; ylab = 'latent heating ice (W/m^2)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''            ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''            ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''            ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0, &
                    nrad_lh,  vctr19, &
                    fe2(: ,isim), &
                    xlab,ylab, &
                    0.0,500.0,10.0,5, &
                    -1500.0,10000.0,400.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'radius (km)' ; ylab = 'latent heating tot (W/m^2)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''            ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''            ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''            ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0, &
                    nrad_lh,  vctr19, &
                    fe3(: ,isim), &
                    xlab,ylab, &
                    0.0,500.0,10.0,5, &
                    -1500.0,10000.0,400.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'time (hours)' ; ylab = 'accumulated precipitation (mm)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''             ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''             ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''             ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0, &
                    numtimes,  vctr18, &
                    ge1(: ,isim), &
                    xlab,ylab, &
                    timebeg,timeend,timeinc,5, &
                    -0.1,60.0,2.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'time (hours)' ; ylab = 'max supercooled LWC (g/kg)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''             ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''             ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''             ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0,  &
                    numtimes,  vctr18, &
                    ge2(: ,isim), &
                    xlab,ylab, &
                    timebeg,timeend,timeinc,5, &
                    -.2,10.0,.2,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'time (hours)' ; ylab = 'max W (m/s)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''             ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''             ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''             ; ylab = ''
     endif

  call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0,  &
                 numtimes,  vctr18, &
                 ge3(: ,isim), &
                 xlab,ylab, &
                 timebeg,timeend,timeinc,5, &
                 -.2,40.0,1.,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'time (hours)' ; ylab = 'max V (m/s)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''             ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''             ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''             ; ylab = ''
     endif

  call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0,  &
                 numtimes,  vctr18, &
                 ge4(: ,isim), &
                 xlab,ylab, &
                 timebeg,timeend,timeinc,5, &
                 -.2,70.0,1.,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'time (hours)' ; ylab = 'max LWC (g/kg)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''             ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''             ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''             ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0,  &
                    numtimes,  vctr18, &
                    ge5(: ,isim), &
                    xlab,ylab, &
                    timebeg,timeend,timeinc,5, &
                    -.2,10.0,.2,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

     ! SPECIAL PLOT OF OBSERVED PMSL_MIN

ae6( 1) = 986.0    ; vctr17( 1) =  0.
ae6( 2) =  982.0   ; vctr17( 2) =  3.
ae6( 3) = 978.0    ; vctr17( 3) =  6.
ae6( 4) =  975.5   ; vctr17( 4) =  9.
ae6( 5) = 973.0    ; vctr17( 5) = 12.
ae6( 6) =  969.5   ; vctr17( 6) = 15.
ae6( 7) = 966.0    ; vctr17( 7) = 18.
ae6( 8) =  957.5   ; vctr17( 8) = 21.
ae6( 9) = 949.0    ; vctr17( 9) = 24.
ae6(10) =  946.0   ; vctr17(10) = 27.
ae6(11) = 943.0    ; vctr17(11) = 30.
ae6(12) =  942.0   ; vctr17(12) = 33.
ae6(13) = 941.0    ; vctr17(13) = 36.
ae6(14) = 937.0    ; vctr17(14) = 39.
ae6(15) = 948.0    ; vctr17(15) = 42.
ae6(16) =  963.0   ; vctr17(16) = 45.
ae6(17) = 978.0    ; vctr17(17) = 48.

     ! Set color to green

     call oplot_xy2('1','N','a','N',aspect,scalelab,9,0,  &
                    17,  vctr17, &
                    ae6(1:17),      &
                    '','', &
                    timebeg,timeend,timeinc,5, &
                    930.,1000.,2.,5  )

  do isim = 1,numsims

     if     (isim == 1) then
        icolor = 10 ; xlab = 'time (hours)' ; ylab = 'min PMSL (mb)'
     elseif (isim == 2) then
        icolor =  1 ; xlab = ''             ; ylab = ''
     elseif (isim == 3) then
        icolor =  8 ; xlab = ''             ; ylab = ''
     elseif (isim == 4) then
        icolor =  9 ; xlab = ''             ; ylab = ''
     endif

     call oplot_xy2('1','N','a','N',aspect,scalelab,icolor,0,  &
                    numtimes,  vctr18, &
                    ge6(: ,isim), &
                    xlab,ylab, &
                    timebeg,timeend,timeinc,5, &
                    930.,1000.,2.,5  )

  enddo

  call o_frame()

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------

  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'max W ',' (m/s)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he1(:,: ,isim), &
                     81,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
GOTO 11

  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'rms +W ',' (m/s)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he2(:,: ,isim), &
                     443,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
11 CONTINUE

  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'max LWC ',' (g/kg)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he3(:,: ,isim), &
                     456,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )
  enddo

  call o_frame()
  !-------------------------------------------------------------------
GOTO 12
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'rms LWC ',' (g/kg)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he4(:,: ,isim), &
                     443,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
12 CONTINUE
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'max PCP ',' (g/kg)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he5(:,: ,isim), &
                     456,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
GOTO 13
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'rms PCP ',' (g/kg)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he6(:,: ,isim), &
                     443,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
13 CONTINUE
  !-------------------------------------------------------------------
!! GOTO 14
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'latent_heating_liq ',' (K/day)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he7(:,: ,isim), &
                     467,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  ! SPECIAL DIFFERENCE PLOT (NA10X MINUS N) 

     do ihe = 1,numtimes
        do khe = 1,ny
           he7(khe,ihe,3) = he7(khe,ihe,3) &
                          - he7(khe,ihe,1)
        enddo
     enddo

     call oplot_zxy2('4','N','a','c',aspect,scalelab, &
                     'diff latent_heating_liq ',' (K/day)', numtimes,ny,vctr18,height, &
                     '','', &
                     he7(:,: ,3), &
                     466,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  call o_frame()
14 CONTINUE
  !-------------------------------------------------------------------
!! GOTO 15
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'latent_heating_ice ',' (K/day)', numtimes,ny,vctr18,height, &
                     xlab,ylab, &
                     he8(:,: ,isim), &
                     468,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  ! SPECIAL DIFFERENCE PLOT (NA10X MINUS N) 

     do ihe = 1,numtimes
        do khe = 1,ny
           he8(khe,ihe,3) = he8(khe,ihe,3) &
                          - he8(khe,ihe,1)
        enddo
     enddo

     call oplot_zxy2('4','N','a','c',aspect,scalelab, &
                     'diff latent_heating_ice ',' (K/day)', numtimes,ny,vctr18,height, &
                     '','', &
                      he8(:,: ,3), &
                     466,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )


  call o_frame()
15 CONTINUE
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'max supercooled LWC ',' (g/kg)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he9(:,: ,isim), &
                     418,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------

  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'max SS ',' (%)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he11(:,: ,isim), &
                     442,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

44 continue

  call o_frame()
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'SS thetadif ',' (K)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he15(:,: ,isim), &
                     462,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  ! SPECIAL DIFFERENCE PLOT (N MINUS NA10X) 

     do ihe = 1,numtimes
        do khe = 1,ny
           he15(khe,ihe,3) = -he15(khe,ihe,3) &
                           +  he15(khe,ihe,1)
        enddo
     enddo

     call oplot_zxy2('4','N','a','c',aspect,scalelab, &
                     'diff SS thetadif',' (K)', numtimes,ny,vctr18,height, & 
                     '','', &
                     he15(:,: ,3), &
                     465,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'area SS > 5% ',' (%)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he12(:,: ,isim), &
                     463,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'area SS > 10% ',' (%)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he13(:,: ,isim), &
                     463,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------
  call plotback()

  do isim = 1,numsims

     if     (isim == 1) then
        ipanel = '3' ; xlab = ''             ; ylab = 'height (km)'
     elseif (isim == 2) then
        ipanel = '1' ; xlab = 'time (hours)' ; ylab = 'height (km)'
     elseif (isim == 3) then
        ipanel = '2' ; xlab = 'time (hours)' ; ylab = ''
     elseif (isim == 4) then
        ipanel = '4' ; xlab = ''             ; ylab = ''
     endif

     call oplot_zxy2(ipanel,'N','a','c',aspect,scalelab, &
                     'area SS > 30% ',' (%)', numtimes,ny,vctr18,height, & 
                     xlab,ylab, &
                     he14(:,: ,isim), &
                     463,ifill,timebeg,timeend,timeinc,5, &
                     0.0,height(ny),1.0,5  )

  enddo

  call o_frame()
  !-------------------------------------------------------------------

6 continue

  ! Close the current workstation if output
  ! is to a NCAR graphics meta file. This allows viewing the complete
  ! meta file (including the last frame) during a run and in case the
  ! simulation crashes.

  if ((trim(runtype) /= 'PLOTONLY') .and. (op%plttype == 0)) then
     call o_clswk()
  endif

end subroutine timeseries_plots
