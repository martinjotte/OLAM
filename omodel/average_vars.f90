subroutine inc_mavg_vars()

  use mem_average_vars, only: &
     npoints_mavg, npoints_avg24, nz_avg, &
     press_mavg, rho_mavg, tempk_mavg, sh_v_mavg, sh_w_mavg, wc_mavg, &
     vxe_mavg, vye_mavg, vze_mavg, &
     rshort_mavg, rshort_top_mavg, rshortup_mavg, rshortup_top_mavg, &
     rlong_mavg, rlongup_mavg, rlongup_top_mavg, &
     latflux_mavg, sensflux_mavg, windspeed_mavg, &
     accpmic_mtot, accpcon_mtot, &
     press_avg24, tempk_avg24, shum_avg24, vxe_avg24, vye_avg24, vze_avg24, &
     rshort_avg24, rshort_top_avg24, rshortup_avg24, rshortup_top_avg24, &
     rlong_avg24, rlongup_avg24, rlongup_top_avg24, &
     latflux_avg24, sensflux_avg24, &
     accpmic_tot24, accpcon_tot24

  use mem_grid,    only: mwa, lpw, zt, zm, mza, dzim, mua, lpu, lpv, mva
  use mem_basic,   only: wc, sh_v, sh_w, tair, press, rho, vxe, vye, vze
  use consts_coms, only: alvl, cp
  use mem_turb,    only: sflux_r, sflux_t
  use misc_coms,   only: io6, dtlong, time8
  use mem_cuparm,  only: conprr
  use mem_micro,   only: pcpgr
  use mem_radiate, only: rshort, rshort_top, rshortup_top, rlong, rlongup, &
                         rlongup_top, albedt
  use mem_ijtabs,  only: jtab_w, jtw_prog

  implicit none

  integer :: iw, k, dhr, kll, j
  real(kind=8) :: dsec
  real :: windspeed
  integer, dimension(nz_avg-1) :: level_list = (/ 8, 13, 18, 26, 30, 34 /)

  npoints_mavg = npoints_mavg + 1

  dsec = mod(time8,86400.d0)
  dhr  = int(dsec / 3600.d0) + 1

  npoints_avg24(dhr) = npoints_avg24(dhr) + 1

! Horizontal loop over all prognostic W/T points

  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

! Loop over k levels from level_list table

     do k = 1, nz_avg
        if (k == 1) kll = lpw(iw)
        if (k > 1)  kll = level_list(k-1)

! Accumulate sums at lpw and at k levels from level_list table

        press_mavg(k,iw) = press_mavg(k,iw) + real(press(kll,iw))
          rho_mavg(k,iw) =   rho_mavg(k,iw) + real(rho  (kll,iw))
        tempk_mavg(k,iw) = tempk_mavg(k,iw) + tair(kll,iw)
         sh_v_mavg(k,iw) =  sh_v_mavg(k,iw) + sh_v(kll,iw)
         sh_w_mavg(k,iw) =  sh_w_mavg(k,iw) + sh_w(kll,iw)
           wc_mavg(k,iw) =    wc_mavg(k,iw) + wc  (kll,iw) 
          vxe_mavg(k,iw) =   vxe_mavg(k,iw) + vxe (kll,iw) 
          vye_mavg(k,iw) =   vye_mavg(k,iw) + vye (kll,iw) 
          vze_mavg(k,iw) =   vze_mavg(k,iw) + vze (kll,iw) 
     enddo

! Accumulate of sums for 2D surface quantities

     k = lpw(iw)

     windspeed = sqrt(vxe(k,iw)**2 + vye(k,iw)**2 + vze(k,iw)**2)

           rshort_mavg(iw) =       rshort_mavg(iw) + rshort      (iw) 
       rshort_top_mavg(iw) =   rshort_top_mavg(iw) + rshort_top  (iw) 
         rshortup_mavg(iw) =     rshortup_mavg(iw) + rshort(iw) * albedt(iw)
     rshortup_top_mavg(iw) = rshortup_top_mavg(iw) + rshortup_top(iw) 
            rlong_mavg(iw) =        rlong_mavg(iw) + rlong       (iw) 
          rlongup_mavg(iw) =      rlongup_mavg(iw) + rlongup     (iw) 
      rlongup_top_mavg(iw) =  rlongup_top_mavg(iw) + rlongup_top (iw) 

       latflux_mavg(iw) =   latflux_mavg(iw) + sflux_r(iw) * alvl
      sensflux_mavg(iw) =  sensflux_mavg(iw) + sflux_t(iw) * cp
     windspeed_mavg(iw) = windspeed_mavg(iw) + windspeed

     accpmic_mtot(iw) = accpmic_mtot(iw) + pcpgr(iw) * dtlong
     if (allocated(conprr)) then
        accpcon_mtot(iw) = accpcon_mtot(iw) + conprr(iw) * dtlong
     endif

! Accumulate sums for individual hours of diurnal cycle

     press_avg24(dhr,iw) = press_avg24(dhr,iw) + real(press(k,iw))
     tempk_avg24(dhr,iw) = tempk_avg24(dhr,iw) + tair(k,iw)
      shum_avg24(dhr,iw) =  shum_avg24(dhr,iw) + sh_v(k,iw)
       vxe_avg24(dhr,iw) =   vxe_avg24(dhr,iw) + vxe(k,iw)
       vye_avg24(dhr,iw) =   vye_avg24(dhr,iw) + vye(k,iw)
       vze_avg24(dhr,iw) =   vze_avg24(dhr,iw) + vze(k,iw)

           rshort_avg24(dhr,iw) =       rshort_avg24(dhr,iw) + rshort(iw)
       rshort_top_avg24(dhr,iw) =   rshort_top_avg24(dhr,iw) + rshort_top(iw)
         rshortup_avg24(dhr,iw) =     rshortup_avg24(dhr,iw) + rshort(iw) * albedt(iw)
     rshortup_top_avg24(dhr,iw) = rshortup_top_avg24(dhr,iw) + rshortup_top(iw)
            rlong_avg24(dhr,iw) =        rlong_avg24(dhr,iw) + rlong(iw)
          rlongup_avg24(dhr,iw) =      rlongup_avg24(dhr,iw) + rlongup(iw)
      rlongup_top_avg24(dhr,iw) =  rlongup_top_avg24(dhr,iw) + rlongup_top(iw)
          latflux_avg24(dhr,iw) =      latflux_avg24(dhr,iw) + sflux_r(iw) * alvl
         sensflux_avg24(dhr,iw) =     sensflux_avg24(dhr,iw) + sflux_t(iw) * cp
          accpmic_tot24(dhr,iw) =      accpmic_tot24(dhr,iw) + pcpgr (iw) * dtlong
     if (allocated(conprr)) then
          accpcon_tot24(dhr,iw) =      accpcon_tot24(dhr,iw) + conprr(iw) * dtlong
     endif

  enddo

  return
end subroutine inc_mavg_vars

!========================================================================

subroutine inc_davg_vars()

  use mem_average_vars, only: &
     npoints_davg, nz_avg, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
       canltempk_davg,   canltempk_dmin,   canltempk_dmax, &
        vegtempk_davg,    vegtempk_dmin,    vegtempk_dmax, &
       soiltempk_davg,   soiltempk_dmin,   soiltempk_dmax, &
       canstempk_davg,   canstempk_dmin,   canstempk_dmax, &
        tempk_lf_davg,    tempk_lf_dmin,    tempk_lf_dmax, &
     cantempk_lf_davg, cantempk_lf_dmin, cantempk_lf_dmax, &
       sfluxt_lf_davg,   sfluxr_lf_davg,                   &
        tempk_sf_davg,    tempk_sf_dmin,    tempk_sf_dmax, &
     cantempk_sf_davg, cantempk_sf_dmin, cantempk_sf_dmax, &
       sfluxt_sf_davg,   sfluxr_sf_davg

  use mem_ijtabs,  only: itabg_w, jtab_w, jtw_prog
  use mem_grid,    only: mwa, lpw, zt, zm, mza, dzim
  use mem_basic,   only: wc, sh_v, sh_w, tair, press, rho, vxe, vye, vze
  use consts_coms, only: alvl, cp
  use mem_turb,    only: sflux_r, sflux_t
  use misc_coms,   only: io6, dtlong, time8, isubdomain
  use mem_cuparm,  only: conprr
  use mem_micro,   only: pcpgr
  use mem_radiate, only: rshort, rshort_top, rshortup_top, rlong, rlongup, &
                         rlongup_top, albedt

  use leaf_coms, only: nzg, mwl, slcpd
  use mem_leaf,  only: land, itabg_wl, itab_wl
  use sea_coms,  only: mws
  use mem_sea,   only: sea, itabg_ws, itab_ws
  use mem_sflux, only: landflux, seaflux, jlandflux, jseaflux
  use mem_para,  only: myrank

  implicit none
  
  integer :: j, iw, k, iwl, iws, kw, ilf, isf
  real :: tempk, cantempk, vegtempk, soiltempk, fracliq
  integer, dimension(nz_avg-1) :: level_list = (/ 8, 13, 18, 26, 30, 34 /)

  npoints_davg = npoints_davg + 1

! Horizontal loop over all prognostic W/T points

  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

! Accumulate of sums for 2D surface quantities

     k = lpw(iw)

      press_davg(iw) =  press_davg(iw) + real(press(k,iw))
        vxe_davg(iw) =    vxe_davg(iw) + vxe(k,iw)
        vye_davg(iw) =    vye_davg(iw) + vye(k,iw)
        vze_davg(iw) =    vze_davg(iw) + vze(k,iw)
     rshort_davg(iw) = rshort_davg(iw) + rshort(iw) 
      tempk_davg(iw) =  tempk_davg(iw) + tair(k,iw)

     if (tempk_dmin(iw) > tair(k,iw)) tempk_dmin(iw) = tair(k,iw)
     if (tempk_dmax(iw) < tair(k,iw)) tempk_dmax(iw) = tair(k,iw)

     accpmic_dtot(iw) = accpmic_dtot(iw) + pcpgr(iw) * dtlong
     if (allocated(conprr)) then
        accpcon_dtot(iw) = accpcon_dtot(iw) + conprr(iw) * dtlong
     endif

     press_ul_davg(iw) = press_ul_davg(iw) + real(press(level_list(5),iw))
       vxe_ul_davg(iw) =   vxe_ul_davg(iw) + vxe(level_list(5),iw)
       vye_ul_davg(iw) =   vye_ul_davg(iw) + vye(level_list(5),iw)
       vze_ul_davg(iw) =   vze_ul_davg(iw) + vze(level_list(5),iw)

  enddo

! Horizontal loop over all land cells

  do iwl = 2, mwl

     if (isubdomain == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

     cantempk = land%can_temp(iwl)
     vegtempk = land%veg_temp(iwl)

     call qwtk(land%soil_energy(nzg,iwl),       &
               land%soil_water(nzg,iwl)*1.e3,   &
               slcpd(land%ntext_soil(nzg,iwl)), &
               soiltempk, fracliq)

     canltempk_davg(iwl) = canltempk_davg(iwl) +  cantempk
      vegtempk_davg(iwl) =  vegtempk_davg(iwl) +  vegtempk
     soiltempk_davg(iwl) = soiltempk_davg(iwl) + soiltempk

     if (canltempk_dmin(iwl) > cantempk) canltempk_dmin(iwl) = cantempk
     if (canltempk_dmax(iwl) < cantempk) canltempk_dmax(iwl) = cantempk

     if (vegtempk_dmin(iwl) > vegtempk) vegtempk_dmin(iwl) = vegtempk
     if (vegtempk_dmax(iwl) < vegtempk) vegtempk_dmax(iwl) = vegtempk

     if (soiltempk_dmin(iwl) > soiltempk) soiltempk_dmin(iwl) = soiltempk
     if (soiltempk_dmax(iwl) < soiltempk) soiltempk_dmax(iwl) = soiltempk
  enddo

! Horizontal loop over all sea cells

  do iws = 2, mws

     if (isubdomain == 1 .and. itab_ws(iws)%irank /= myrank) cycle

     cantempk = sea%can_temp(iws)

     canstempk_davg(iws) = canstempk_davg(iws) + cantempk

     if (canstempk_dmin(iws) > cantempk) canstempk_dmin(iws) = cantempk
     if (canstempk_dmax(iws) < cantempk) canstempk_dmax(iws) = cantempk
  enddo

! Horizontal loop over all landflux cells

  do j = 1,jlandflux(1)%jend(1)
     ilf = jlandflux(1)%ilandflux(j)
     iw  = landflux(ilf)%iw         ! global index
     kw  = landflux(ilf)%kw
     iwl = landflux(ilf)%iwls       ! global index

     if (isubdomain == 1) then
        iw  = itabg_w (iw )%iw_myrank  ! Local rank index (if MPI parallel)
        iwl = itabg_wl(iwl)%iwl_myrank ! Local rank index (if MPI parallel)
     endif

     cantempk = land%can_temp(iwl)

        tempk_lf_davg(ilf) =    tempk_lf_davg(ilf) + tair(kw,iw)
     cantempk_lf_davg(ilf) = cantempk_lf_davg(ilf) + cantempk
       sfluxt_lf_davg(ilf) =   sfluxt_lf_davg(ilf) + landflux(ilf)%sfluxt
       sfluxr_lf_davg(ilf) =   sfluxr_lf_davg(ilf) + landflux(ilf)%sfluxr

     if (tempk_lf_dmin(ilf) > tair(kw,iw)) tempk_lf_dmin(ilf) = tair(kw,iw)
     if (tempk_lf_dmax(ilf) < tair(kw,iw)) tempk_lf_dmax(ilf) = tair(kw,iw)

     if (cantempk_lf_dmin(ilf) > cantempk) cantempk_lf_dmin(ilf) = cantempk
     if (cantempk_lf_dmax(ilf) < cantempk) cantempk_lf_dmax(ilf) = cantempk
  enddo

! Horizontal loop over all seaflux cells

  do j = 1,jseaflux(1)%jend(1)
     isf = jseaflux(1)%iseaflux(j)
     iw  = seaflux(isf)%iw         ! global index
     kw  = seaflux(isf)%kw
     iws = seaflux(isf)%iwls           ! global index

     if (isubdomain == 1) then
        iw  = itabg_w (iw )%iw_myrank  ! Local rank index (if MPI parallel)
        iws = itabg_ws(iws)%iws_myrank ! Local rank index (if MPI parallel)
     endif

     cantempk = sea%can_temp(iws)

        tempk_sf_davg(isf) =    tempk_sf_davg(isf) + tair(kw,iw)
     cantempk_sf_davg(isf) = cantempk_sf_davg(isf) + cantempk
       sfluxt_sf_davg(isf) =   sfluxt_sf_davg(isf) + seaflux(isf)%sfluxt
       sfluxr_sf_davg(isf) =   sfluxr_sf_davg(isf) + seaflux(isf)%sfluxr

     if (tempk_sf_dmin(isf) > tair(kw,iw)) tempk_sf_dmin(isf) = tair(kw,iw)
     if (tempk_sf_dmax(isf) < tair(kw,iw)) tempk_sf_dmax(isf) = tair(kw,iw)

     if (cantempk_sf_dmin(isf) > cantempk) cantempk_sf_dmin(isf) = cantempk
     if (cantempk_sf_dmax(isf) < cantempk) cantempk_sf_dmax(isf) = cantempk
  enddo

  return
end subroutine inc_davg_vars

!========================================================================

subroutine norm_mavg_vars()

  use mem_average_vars, only: &
     npoints_mavg, npoints_avg24, &
     press_mavg, rho_mavg, tempk_mavg, sh_v_mavg, sh_w_mavg, wc_mavg, &
     vxe_mavg, vye_mavg, vze_mavg, &
     rshort_mavg, rshort_top_mavg, rshortup_mavg, rshortup_top_mavg, &
     rlong_mavg, rlongup_mavg, rlongup_top_mavg, &
     latflux_mavg, sensflux_mavg, windspeed_mavg, &
     press_avg24, tempk_avg24, shum_avg24, vxe_avg24, vye_avg24, vze_avg24, &
     rshort_avg24, rshort_top_avg24, rshortup_avg24, rshortup_top_avg24, &
     rlong_avg24, rlongup_avg24, rlongup_top_avg24, &
     latflux_avg24, sensflux_avg24

  use mem_grid,   only: mwa
  use mem_ijtabs, only: jtab_w, jtw_prog

  implicit none

  integer :: iw, ih, j
  real :: n24i(24), ni

  ni = 1.0 / npoints_mavg
  do ih = 1, 24
     n24i(ih) = 1.0 / npoints_avg24(ih)
  enddo

  press_mavg = press_mavg * ni
    rho_mavg =   rho_mavg * ni
  tempk_mavg = tempk_mavg * ni
   sh_v_mavg =  sh_v_mavg * ni
   sh_w_mavg =  sh_w_mavg * ni
     wc_mavg =    wc_mavg * ni
    vxe_mavg =   vxe_mavg * ni
    vye_mavg =   vye_mavg * ni
    vze_mavg =   vze_mavg * ni

    do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

           rshort_mavg(iw) =       rshort_mavg(iw) * ni
       rshort_top_mavg(iw) =   rshort_top_mavg(iw) * ni
         rshortup_mavg(iw) =     rshortup_mavg(iw) * ni
     rshortup_top_mavg(iw) = rshortup_top_mavg(iw) * ni
            rlong_mavg(iw) =        rlong_mavg(iw) * ni
          rlongup_mavg(iw) =      rlongup_mavg(iw) * ni
      rlongup_top_mavg(iw) =  rlongup_top_mavg(iw) * ni
          latflux_mavg(iw) =      latflux_mavg(iw) * ni
         sensflux_mavg(iw) =     sensflux_mavg(iw) * ni
        windspeed_mavg(iw) =    windspeed_mavg(iw) * ni
     
     do ih = 1, 24
        press_avg24(ih,iw) = press_avg24(ih,iw) * n24i(ih)
        tempk_avg24(ih,iw) = tempk_avg24(ih,iw) * n24i(ih)
         shum_avg24(ih,iw) =  shum_avg24(ih,iw) * n24i(ih)
          vxe_avg24(ih,iw) =   vxe_avg24(ih,iw) * n24i(ih)
          vye_avg24(ih,iw) =   vye_avg24(ih,iw) * n24i(ih)
          vze_avg24(ih,iw) =   vze_avg24(ih,iw) * n24i(ih)

              rshort_avg24(ih,iw) =       rshort_avg24(ih,iw) * n24i(ih)
          rshort_top_avg24(ih,iw) =   rshort_top_avg24(ih,iw) * n24i(ih)
            rshortup_avg24(ih,iw) =     rshortup_avg24(ih,iw) * n24i(ih)
        rshortup_top_avg24(ih,iw) = rshortup_top_avg24(ih,iw) * n24i(ih)
               rlong_avg24(ih,iw) =        rlong_avg24(ih,iw) * n24i(ih)
             rlongup_avg24(ih,iw) =      rlongup_avg24(ih,iw) * n24i(ih)
         rlongup_top_avg24(ih,iw) =  rlongup_top_avg24(ih,iw) * n24i(ih)
             latflux_avg24(ih,iw) =      latflux_avg24(ih,iw) * n24i(ih)
            sensflux_avg24(ih,iw) =     sensflux_avg24(ih,iw) * n24i(ih)
     enddo

  enddo

  return
end subroutine norm_mavg_vars

!========================================================================

subroutine norm_davg_vars()

  use mem_average_vars, only: &
     npoints_davg, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, tempk_davg, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     canltempk_davg,    vegtempk_davg, soiltempk_davg, canstempk_davg, &
      tempk_lf_davg, cantempk_lf_davg, sfluxt_lf_davg, sfluxr_lf_davg, &
      tempk_sf_davg, cantempk_sf_davg, sfluxt_sf_davg, sfluxr_sf_davg

  use mem_grid,  only: mwa
  use leaf_coms, only: mwl
  use mem_leaf,  only: itab_wl
  use sea_coms,  only: mws
  use mem_sea,   only: itab_ws
  use mem_sflux, only: jlandflux, jseaflux
  use mem_ijtabs,only: jtab_w, jtw_prog
  use misc_coms, only: isubdomain
  use mem_para,  only: myrank

  implicit none

  integer :: iw, iwl, iws, j, ilf, isf
  real :: ni

  ni = 1. / npoints_davg

  do j = 1, jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

        press_davg(iw) =    press_davg(iw) * ni
          vxe_davg(iw) =      vxe_davg(iw) * ni
          vye_davg(iw) =      vye_davg(iw) * ni
          vze_davg(iw) =      vze_davg(iw) * ni
       rshort_davg(iw) =   rshort_davg(iw) * ni
        tempk_davg(iw) =    tempk_davg(iw) * ni
     press_ul_davg(iw) = press_ul_davg(iw) * ni
       vxe_ul_davg(iw) =   vxe_ul_davg(iw) * ni
       vye_ul_davg(iw) =   vye_ul_davg(iw) * ni
       vze_ul_davg(iw) =   vze_ul_davg(iw) * ni
  enddo

  do iwl = 2, mwl
     if (isubdomain == 1 .and. itab_wl(iwl)%irank /= myrank) cycle
     canltempk_davg(iwl) = canltempk_davg(iwl) * ni
      vegtempk_davg(iwl) =  vegtempk_davg(iwl) * ni
     soiltempk_davg(iwl) = soiltempk_davg(iwl) * ni
  enddo

  do iws = 2, mws
     if (isubdomain == 1 .and. itab_ws(iws)%irank /= myrank) cycle
     canstempk_davg(iws) = canstempk_davg(iws) * ni
  enddo

  do j = 1,jlandflux(1)%jend(1)
     ilf = jlandflux(1)%ilandflux(j)

        tempk_lf_davg(ilf) =    tempk_lf_davg(ilf) * ni
     cantempk_lf_davg(ilf) = cantempk_lf_davg(ilf) * ni
       sfluxt_lf_davg(ilf) =   sfluxt_lf_davg(ilf) * ni
       sfluxr_lf_davg(ilf) =   sfluxr_lf_davg(ilf) * ni
  enddo

! Horizontal loop over all seaflux cells

  do j = 1,jseaflux(1)%jend(1)
     isf = jseaflux(1)%iseaflux(j)

        tempk_sf_davg(isf) =    tempk_sf_davg(isf) * ni
     cantempk_sf_davg(isf) = cantempk_sf_davg(isf) * ni
       sfluxt_sf_davg(isf) =   sfluxt_sf_davg(isf) * ni
       sfluxr_sf_davg(isf) =   sfluxr_sf_davg(isf) * ni
  enddo

  return
end subroutine norm_davg_vars

!===============================================================

subroutine write_mavg_vars(outyear,outmonth)
  
  use mem_average_vars, only: &
     nz_avg, &
     press_mavg, rho_mavg, tempk_mavg, sh_v_mavg, sh_w_mavg, wc_mavg, &
     vxe_mavg, vye_mavg, vze_mavg, &
     rshort_mavg, rshort_top_mavg, rshortup_mavg, rshortup_top_mavg, &
     rlong_mavg, rlongup_mavg, rlongup_top_mavg, &
     latflux_mavg, sensflux_mavg, windspeed_mavg, &
     accpmic_mtot, accpcon_mtot, &
     press_avg24, tempk_avg24, shum_avg24, vxe_avg24, vye_avg24, vze_avg24, &
     rshort_avg24, rshort_top_avg24, rshortup_avg24, rshortup_top_avg24, &
     rlong_avg24, rlongup_avg24, rlongup_top_avg24, &
     latflux_avg24, sensflux_avg24, &
     accpmic_tot24, accpcon_tot24

  use misc_coms,  only: iparallel, isubdomain, hfilepref, iclobber, io6, time8, iyear1, &
                        imonth1, idate1, itime1, ipar_out
  use mem_para,   only: myrank, iwa_globe_primary, iwa_local_primary
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa, nwa
  use leaf_coms,  only: mwl, nzg, nzs, dslz
  use mem_leaf,   only: land, itab_wl

  implicit none

  integer, intent(in) :: outyear, outmonth
  character(len=10) :: post
  character(len=128) :: hnamel
  integer :: lenl, iwl, k
  logical exans
  character(len=32) :: varn
  integer :: ndims, idims(2), month_use, year_use
  real :: wstorage(mwl)
  real :: tot_soil_water
  integer, pointer :: ilpts(:), igpts(:)
  integer :: nglobe
  
  ! Compute total water

  do iwl = 2, mwl
     if (isubdomain == 1 .and. itab_wl(iwl)%irank /= myrank) cycle
     wstorage(iwl) = sum(land%sfcwater_mass(1:nzs,iwl)) +   &
          land%veg_water(iwl) +   &
          land%rhos(iwl) * land%can_depth(iwl) * land%can_shv(iwl)
     tot_soil_water = 0.0
     do k = 1, nzg
        tot_soil_water = tot_soil_water +   &
             dslz(k) * 1000. * land%soil_water(k,iwl)
     enddo
     wstorage(iwl) = wstorage(iwl) + tot_soil_water
  enddo

! Determine output file name and open file

  if(iparallel == 0 .or. ipar_out == 1) then
     post = '$'
  else
     write(post,'(i10)')myrank
     post = 'r'//trim(adjustl(post))
  endif

  month_use = outmonth - 1
  year_use = outyear
  if(month_use == 0)then
     month_use = 12
     year_use = year_use - 1
  endif

  call makefnam8(hnamel,hfilepref,0.d0,year_use,month_use,1,  &
       0,'M',post,'h5')  ! '$' argument suppresses grid# encoding

  lenl = len_trim(hnamel)
  
  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*)'!!!   trying to open file name :'
     write(io6,*)'!!!       ',hnamel
     write(io6,*)'!!!   but it already exists. run is ended.'
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'history_write'
  endif
  
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  write(io6,*)'history_write:open file:',trim(hnamel)
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  call shdf5_open(hnamel,'W',iclobber)
  
  !  Write month-average variables
  
  if (iparallel == 1 .and. ipar_out == 1) then
     ilpts => iwa_local_primary
     igpts => iwa_globe_primary
     nglobe = nwa

     ndims    = 2
     idims(1) = nz_avg
     idims(2) = mwa

     call shdf5_orec(ndims,idims,'PRESS_MAVG',rvara=press_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'RHO_MAVG',rvara=  rho_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'TEMPK_MAVG',rvara=tempk_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'SH_V_MAVG',rvara= sh_v_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'SH_W_MAVG',rvara= sh_w_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,   'WC_MAVG',rvara=   wc_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'VXE_MAVG',rvara=  vxe_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'VYE_MAVG',rvara=  vye_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'VZE_MAVG',rvara=  vze_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     ndims    = 1
     idims(1) = mwa
     ilpts => iwa_local_primary
     igpts => iwa_globe_primary
     nglobe = nwa

     call shdf5_orec(ndims,idims,      'RSHORT_MAVG',rvara=      rshort_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'RSHORT_TOP_MAVG',rvara=  rshort_top_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,    'RSHORTUP_MAVG',rvara=    rshortup_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'RSHORTUP_TOP_MAVG',rvara=rshortup_top_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,       'RLONG_MAVG',rvara=       rlong_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,     'RLONGUP_MAVG',rvara=     rlongup_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'RLONGUP_TOP_MAVG',rvara= rlongup_top_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,     'LATFLUX_MAVG',rvara=     latflux_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,    'SENSFLUX_MAVG',rvara=    sensflux_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,   'WINDSPEED_MAVG',rvara=   windspeed_mavg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,     'ACCPMIC_MTOT',rvara=     accpmic_mtot, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,     'ACCPCON_MTOT',rvara=     accpcon_mtot, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     call shdf5_orec(ndims,idims,        'TOT_WATER',rvara=         wstorage, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  else

     ndims      = 2
     idims(1) = nz_avg
     idims(2) = mwa

     call shdf5_orec(ndims,idims,'PRESS_MAVG',rvara=press_mavg)
     call shdf5_orec(ndims,idims,  'RHO_MAVG',rvara=  rho_mavg)
     call shdf5_orec(ndims,idims,'TEMPK_MAVG',rvara=tempk_mavg)
     call shdf5_orec(ndims,idims, 'SH_V_MAVG',rvara= sh_v_mavg)
     call shdf5_orec(ndims,idims, 'SH_W_MAVG',rvara= sh_w_mavg)
     call shdf5_orec(ndims,idims,   'WC_MAVG',rvara=   wc_mavg)
     call shdf5_orec(ndims,idims,  'VXE_MAVG',rvara=  vxe_mavg)
     call shdf5_orec(ndims,idims,  'VYE_MAVG',rvara=  vye_mavg)
     call shdf5_orec(ndims,idims,  'VZE_MAVG',rvara=  vze_mavg)

     ndims      = 1
     idims(1) = mwa

     call shdf5_orec(ndims,idims,      'RSHORT_MAVG',rvara=      rshort_mavg)
     call shdf5_orec(ndims,idims,  'RSHORT_TOP_MAVG',rvara=  rshort_top_mavg)
     call shdf5_orec(ndims,idims,    'RSHORTUP_MAVG',rvara=    rshortup_mavg)
     call shdf5_orec(ndims,idims,'RSHORTUP_TOP_MAVG',rvara=rshortup_top_mavg)
     call shdf5_orec(ndims,idims,       'RLONG_MAVG',rvara=       rlong_mavg)
     call shdf5_orec(ndims,idims,     'RLONGUP_MAVG',rvara=     rlongup_mavg)
     call shdf5_orec(ndims,idims, 'RLONGUP_TOP_MAVG',rvara= rlongup_top_mavg)
     call shdf5_orec(ndims,idims,     'LATFLUX_MAVG',rvara=     latflux_mavg)
     call shdf5_orec(ndims,idims,    'SENSFLUX_MAVG',rvara=    sensflux_mavg)
     call shdf5_orec(ndims,idims,   'WINDSPEED_MAVG',rvara=   windspeed_mavg)
     call shdf5_orec(ndims,idims,     'ACCPMIC_MTOT',rvara=     accpmic_mtot)
     call shdf5_orec(ndims,idims,     'ACCPCON_MTOT',rvara=     accpcon_mtot)
     call shdf5_orec(ndims,idims,        'TOT_WATER',rvara=         wstorage)

  endif

!  ndims = 1
!  idims(1) = mwl
!  call shdf5_orec(ndims,idims,'TOT_RUNOFF',rvara=tot_runoff)

!  ndims = 2
!  idims(1) = 24
!  idims(2) = mwa

!  call shdf5_orec(ndims,idims,'PRESS_AVG24',rvara=press_avg24)
!  call shdf5_orec(ndims,idims,'TEMPK_AVG24',rvara=tempk_avg24)
!  call shdf5_orec(ndims,idims, 'SHUM_AVG24',rvara= shum_avg24)
!  call shdf5_orec(ndims,idims,  'VXE_AVG24',rvara=  vxe_avg24)
!  call shdf5_orec(ndims,idims,  'VYE_AVG24',rvara=  vye_avg24)
!  call shdf5_orec(ndims,idims,  'VZE_AVG24',rvara=  vze_avg24)

!  call shdf5_orec(ndims,idims,      'RSHORT_AVG24',rvara=      rshort_avg24)
!  call shdf5_orec(ndims,idims,  'RSHORT_TOP_AVG24',rvara=  rshort_top_avg24)
!  call shdf5_orec(ndims,idims,    'RSHORTUP_AVG24',rvara=    rshortup_avg24)
!  call shdf5_orec(ndims,idims,'RSHORTUP_TOP_AVG24',rvara=rshortup_top_avg24)
!  call shdf5_orec(ndims,idims,       'RLONG_AVG24',rvara=       rlong_avg24)
!  call shdf5_orec(ndims,idims,     'RLONGUP_AVG24',rvara=     rlongup_avg24)
!  call shdf5_orec(ndims,idims, 'RLONGUP_TOP_AVG24',rvara= rlongup_top_avg24)
!  call shdf5_orec(ndims,idims,     'LATFLUX_AVG24',rvara=     latflux_avg24)
!  call shdf5_orec(ndims,idims,    'SENSFLUX_AVG24',rvara=    sensflux_avg24)
!  call shdf5_orec(ndims,idims,     'ACCPMIC_TOT24',rvara=     accpmic_tot24)
!  call shdf5_orec(ndims,idims,     'ACCPCON_TOT24',rvara=     accpcon_tot24)

  call shdf5_close()
  
  return
end subroutine write_mavg_vars

!===============================================================

subroutine write_davg_vars(outyear,outmonth,outdate)
  
  use mem_average_vars, only: &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     canltempk_davg, canltempk_dmin, canltempk_dmax, &
      vegtempk_davg,  vegtempk_dmin,  vegtempk_dmax, &
     soiltempk_davg, soiltempk_dmin, soiltempk_dmax, &
     canstempk_davg, canstempk_dmin, canstempk_dmax, &
        tempk_lf_davg,    tempk_lf_dmin,    tempk_lf_dmax, &
     cantempk_lf_davg, cantempk_lf_dmin, cantempk_lf_dmax, &
       sfluxt_lf_davg,   sfluxr_lf_davg,                   &
        tempk_sf_davg,    tempk_sf_dmin,    tempk_sf_dmax, &
     cantempk_sf_davg, cantempk_sf_dmin, cantempk_sf_dmax, &
       sfluxt_sf_davg,   sfluxr_sf_davg

  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa, nwa
  use misc_coms,  only: iparallel, hfilepref, iclobber, io6, time8, iyear1, &
                        imonth1, idate1, itime1, ipar_out
  use mem_para,   only: myrank, iwa_globe_primary, iwa_local_primary, &
                                iwl_globe_primary, iwl_local_primary, &
                                iws_globe_primary, iws_local_primary, &
                                ifl_globe_primary, ifl_local_primary, &
                                ifs_globe_primary, ifs_local_primary

  use leaf_coms,  only: mwl, nwl
  use sea_coms,   only: mws, nws
  use mem_sflux,  only: nseaflux, mseaflux, nlandflux, mlandflux

  implicit none

  integer, intent(in) :: outyear, outmonth, outdate
  character(len=10) :: post
  character(len=128) :: hnamel
  integer :: lenl
  logical exans
  character(len=32) :: varn
  integer :: ndims,idims(2), month_use, year_use, date_use
  integer, pointer :: ilpts(:), igpts(:)
  integer :: nglobe

! Determine output file name and open file

  if(iparallel == 0 .or. ipar_out == 1) then
     post = '$'
  else
     write(post,'(i10)')myrank
     post = 'r'//trim(adjustl(post))
  endif

  date_use = outdate - 1
  month_use = outmonth
  year_use = outyear
  if(date_use == 0)then
     month_use = outmonth - 1
     if(month_use == 0)then
        month_use = 12
        year_use = year_use - 1
     endif
     date_use = 31
     if(month_use == 4 .or. month_use == 6 .or. month_use == 9   &
          .or. month_use == 11)then
        date_use = 30
     elseif(month_use == 2)then
        date_use = 28
        if(mod(year_use,4) == 0)date_use = 29
     endif
  endif

  call makefnam8(hnamel,hfilepref,0.d0,year_use,month_use,date_use,  &
       0,'D',post,'h5')  ! '$' argument suppresses grid# encoding

  lenl = len_trim(hnamel)
  
  inquire(file=hnamel,exist=exans)
  if (exans .and. iclobber == 0) then
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*)'!!!   trying to open file name :'
     write(io6,*)'!!!       ',hnamel
     write(io6,*)'!!!   but it already exists. run is ended.'
     write(io6,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'history_write'
  endif

  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  write(io6,*)'history_write:open file:',trim(hnamel)
  write(io6,*)'+++++++++++++++++++++++++++++++++++++'
  call shdf5_open(hnamel,'W',iclobber)

  !  Write day-average variables

  if (iparallel .and. ipar_out == 1) then

     ndims    = 1
     idims(1) = mwa

     ilpts => iwa_local_primary
     igpts => iwa_globe_primary
     nglobe = nwa

     call shdf5_orec(ndims,idims,  'PRESS_DAVG',rvara=  press_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,    'VXE_DAVG',rvara=    vxe_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,    'VYE_DAVG',rvara=    vye_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,    'VZE_DAVG',rvara=    vze_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'RSHORT_DAVG',rvara= rshort_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'TEMPK_DAVG',rvara=  tempk_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'TEMPK_DMIN',  rvara=tempk_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'TEMPK_DMAX',  rvara=tempk_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'ACCPMIC_DTOT',rvara=accpmic_dtot, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'ACCPCON_DTOT',rvara=accpcon_dtot, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     call shdf5_orec(ndims,idims,'PRESS_UL_DAVG',rvara=press_ul_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'VXE_UL_DAVG',rvara=  vxe_ul_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'VYE_UL_DAVG',rvara=  vye_ul_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'VZE_UL_DAVG',rvara=  vze_ul_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     idims(1) = mwl

     ilpts => iwl_local_primary
     igpts => iwl_globe_primary
     nglobe = nwl

     call shdf5_orec(ndims,idims,'CANLTEMPK_DAVG',rvara=canltempk_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANLTEMPK_DMIN',rvara=canltempk_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANLTEMPK_DMAX',rvara=canltempk_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'VEGTEMPK_DAVG',rvara= vegtempk_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'VEGTEMPK_DMIN',rvara= vegtempk_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims, 'VEGTEMPK_DMAX',rvara= vegtempk_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'SOILTEMPK_DAVG',rvara=soiltempk_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'SOILTEMPK_DMIN',rvara=soiltempk_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'SOILTEMPK_DMAX',rvara=soiltempk_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     idims(1) = mws

     ilpts => iws_local_primary
     igpts => iws_globe_primary
     nglobe = nws

     call shdf5_orec(ndims,idims,'CANSTEMPK_DAVG',rvara=canstempk_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANSTEMPK_DMIN',rvara=canstempk_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANSTEMPK_DMAX',rvara=canstempk_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     idims(1) = mlandflux

     ilpts => ifl_local_primary
     igpts => ifl_globe_primary
     nglobe = nlandflux

     call shdf5_orec(ndims,idims,   'TEMPK_LF_DAVG',rvara=   tempk_lf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,   'TEMPK_LF_DMIN',rvara=   tempk_lf_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,   'TEMPK_LF_DMAX',rvara=   tempk_lf_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANTEMPK_LF_DAVG',rvara=cantempk_lf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANTEMPK_LF_DMIN',rvara=cantempk_lf_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANTEMPK_LF_DMAX',rvara=cantempk_lf_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'SFLUXT_LF_DAVG',rvara=  sfluxt_lf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'SFLUXR_LF_DAVG',rvara=  sfluxr_lf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

     idims(1) = mseaflux

     ilpts => ifs_local_primary
     igpts => ifs_globe_primary
     nglobe = nseaflux

     call shdf5_orec(ndims,idims,   'TEMPK_SF_DAVG',rvara=   tempk_sf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,   'TEMPK_SF_DMIN',rvara=   tempk_sf_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,   'TEMPK_SF_DMAX',rvara=   tempk_sf_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANTEMPK_SF_DAVG',rvara=cantempk_sf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANTEMPK_SF_DMIN',rvara=cantempk_sf_dmin, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,'CANTEMPK_SF_DMAX',rvara=cantempk_sf_dmax, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'SFLUXT_SF_DAVG',rvara=  sfluxt_sf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
     call shdf5_orec(ndims,idims,  'SFLUXR_SF_DAVG',rvara=  sfluxr_sf_davg, &
          lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  else

     ndims    = 1
     idims(1) = mwa

     call shdf5_orec(ndims,idims,  'PRESS_DAVG',rvara=  press_davg)
     call shdf5_orec(ndims,idims,    'VXE_DAVG',rvara=    vxe_davg)
     call shdf5_orec(ndims,idims,    'VYE_DAVG',rvara=    vye_davg)
     call shdf5_orec(ndims,idims,    'VZE_DAVG',rvara=    vze_davg)
     call shdf5_orec(ndims,idims, 'RSHORT_DAVG',rvara= rshort_davg)
     call shdf5_orec(ndims,idims,  'TEMPK_DAVG',rvara=  tempk_davg)
     call shdf5_orec(ndims,idims,  'TEMPK_DMIN',rvara=  tempk_dmin)
     call shdf5_orec(ndims,idims,  'TEMPK_DMAX',rvara=  tempk_dmax)
     call shdf5_orec(ndims,idims,'ACCPMIC_DTOT',rvara=accpmic_dtot)
     call shdf5_orec(ndims,idims,'ACCPCON_DTOT',rvara=accpcon_dtot)

     call shdf5_orec(ndims,idims,'PRESS_UL_DAVG',rvara=press_ul_davg)
     call shdf5_orec(ndims,idims,  'VXE_UL_DAVG',  rvara=vxe_ul_davg)
     call shdf5_orec(ndims,idims,  'VYE_UL_DAVG',  rvara=vye_ul_davg)
     call shdf5_orec(ndims,idims,  'VZE_UL_DAVG',  rvara=vze_ul_davg)

     idims(1) = mwl

     call shdf5_orec(ndims,idims,'CANLTEMPK_DAVG',rvara=canltempk_davg)
     call shdf5_orec(ndims,idims, 'VEGTEMPK_DAVG',rvara= vegtempk_davg)
     call shdf5_orec(ndims,idims,'SOILTEMPK_DAVG',rvara=soiltempk_davg)
     call shdf5_orec(ndims,idims,'CANLTEMPK_DMIN',rvara=canltempk_dmin)
     call shdf5_orec(ndims,idims,'CANLTEMPK_DMAX',rvara=canltempk_dmax)
     call shdf5_orec(ndims,idims, 'VEGTEMPK_DMIN',rvara= vegtempk_dmin)
     call shdf5_orec(ndims,idims, 'VEGTEMPK_DMAX',rvara= vegtempk_dmax)
     call shdf5_orec(ndims,idims,'SOILTEMPK_DMIN',rvara=soiltempk_dmin)
     call shdf5_orec(ndims,idims,'SOILTEMPK_DMAX',rvara=soiltempk_dmax)

     idims(1) = mws

     call shdf5_orec(ndims,idims,'CANSTEMPK_DAVG',rvara=canstempk_davg)
     call shdf5_orec(ndims,idims,'CANSTEMPK_DMIN',rvara=canstempk_dmin)
     call shdf5_orec(ndims,idims,'CANSTEMPK_DMAX',rvara=canstempk_dmax)

     idims(1) = mlandflux

     call shdf5_orec(ndims,idims,   'TEMPK_LF_DAVG',rvara=   tempk_lf_davg)
     call shdf5_orec(ndims,idims,   'TEMPK_LF_DMIN',rvara=   tempk_lf_dmin)
     call shdf5_orec(ndims,idims,   'TEMPK_LF_DMAX',rvara=   tempk_lf_dmax)
     call shdf5_orec(ndims,idims,'CANTEMPK_LF_DAVG',rvara=cantempk_lf_davg)
     call shdf5_orec(ndims,idims,'CANTEMPK_LF_DMIN',rvara=cantempk_lf_dmin)
     call shdf5_orec(ndims,idims,'CANTEMPK_LF_DMAX',rvara=cantempk_lf_dmax)
     call shdf5_orec(ndims,idims,  'SFLUXT_LF_DAVG',rvara=  sfluxt_lf_davg)
     call shdf5_orec(ndims,idims,  'SFLUXR_LF_DAVG',rvara=  sfluxr_lf_davg)

     idims(1) = mseaflux

     call shdf5_orec(ndims,idims,   'TEMPK_SF_DAVG',rvara=   tempk_sf_davg)
     call shdf5_orec(ndims,idims,   'TEMPK_SF_DMIN',rvara=   tempk_sf_dmin)
     call shdf5_orec(ndims,idims,   'TEMPK_SF_DMAX',rvara=   tempk_sf_dmax)
     call shdf5_orec(ndims,idims,'CANTEMPK_SF_DAVG',rvara=cantempk_sf_davg)
     call shdf5_orec(ndims,idims,'CANTEMPK_SF_DMIN',rvara=cantempk_sf_dmin)
     call shdf5_orec(ndims,idims,'CANTEMPK_SF_DMAX',rvara=cantempk_sf_dmax)
     call shdf5_orec(ndims,idims,  'SFLUXT_SF_DAVG',rvara=  sfluxt_sf_davg)
     call shdf5_orec(ndims,idims,  'SFLUXR_SF_DAVG',rvara=  sfluxr_sf_davg)

  endif

  call shdf5_close()

  return
end subroutine write_davg_vars

!===============================================================

subroutine read_mavg_vars(mavgfile)

  use mem_average_vars, only: &
     nz_avg, &
     press_mavg, rho_mavg, tempk_mavg, sh_v_mavg, sh_w_mavg, wc_mavg, &
     vxe_mavg, vye_mavg, vze_mavg, &
     rshort_mavg, rshort_top_mavg, rshortup_mavg, rshortup_top_mavg, &
     rlong_mavg, rlongup_mavg, rlongup_top_mavg, &
     latflux_mavg, sensflux_mavg, windspeed_mavg, &
     accpmic_mtot, accpcon_mtot, &
     press_avg24, tempk_avg24, shum_avg24, vxe_avg24, vye_avg24, vze_avg24, &
     rshort_avg24, rshort_top_avg24, rshortup_avg24, rshortup_top_avg24, &
     rlong_avg24, rlongup_avg24, rlongup_top_avg24, &
     latflux_avg24, sensflux_avg24, &
     accpmic_tot24, accpcon_tot24

  use misc_coms, only: io6
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_grid, only: mwa, nwa

  implicit none

  character(*), intent(in) :: mavgfile

  logical exans
  integer :: ndims, idims(2)
  character(len=32) :: varn
  
  inquire(file=mavgfile,exist=exans)

  if (exans) then
  
! Month-average file exists.  Open, read, and close file.

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening month-average file ', trim(mavgfile)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(trim(mavgfile),'R')

 ! Read month-average variables
  
     ndims    = 2
     idims(1) = nz_avg
     idims(2) = mwa

     call shdf5_irec(ndims,idims,'PRESS_MAVG',rvara=press_mavg)
     call shdf5_irec(ndims,idims,  'RHO_MAVG',rvara=  rho_mavg)
     call shdf5_irec(ndims,idims,'TEMPK_MAVG',rvara=tempk_mavg)
     call shdf5_irec(ndims,idims, 'SH_V_MAVG',rvara= sh_v_mavg)
     call shdf5_irec(ndims,idims, 'SH_W_MAVG',rvara= sh_w_mavg)
     call shdf5_irec(ndims,idims,   'WC_MAVG',rvara=   wc_mavg)
     call shdf5_irec(ndims,idims,  'VXE_MAVG',rvara=  vxe_mavg)
     call shdf5_irec(ndims,idims,  'VYE_MAVG',rvara=  vye_mavg)
     call shdf5_irec(ndims,idims,  'VZE_MAVG',rvara=  vze_mavg)

     ndims    = 1
     idims(1) = mwa

     call shdf5_irec(ndims,idims,      'RSHORT_MAVG',rvara=      rshort_mavg)
     call shdf5_irec(ndims,idims,  'RSHORT_TOP_MAVG',rvara=  rshort_top_mavg)
     call shdf5_irec(ndims,idims,    'RSHORTUP_MAVG',rvara=    rshortup_mavg)
     call shdf5_irec(ndims,idims,'RSHORTUP_TOP_MAVG',rvara=rshortup_top_mavg)
     call shdf5_irec(ndims,idims,       'RLONG_MAVG',rvara=       rlong_mavg)
     call shdf5_irec(ndims,idims,     'RLONGUP_MAVG',rvara=     rlongup_mavg)
     call shdf5_irec(ndims,idims, 'RLONGUP_TOP_MAVG',rvara= rlongup_top_mavg)
     call shdf5_irec(ndims,idims,     'LATFLUX_MAVG',rvara=     latflux_mavg)
     call shdf5_irec(ndims,idims,    'SENSFLUX_MAVG',rvara=    sensflux_mavg)
     call shdf5_irec(ndims,idims,   'WINDSPEED_MAVG',rvara=   windspeed_mavg)
     call shdf5_irec(ndims,idims,     'ACCPMIC_MTOT',rvara=     accpmic_mtot)
     call shdf5_irec(ndims,idims,     'ACCPCON_MTOT',rvara=     accpcon_mtot)
!     call shdf5_irec(ndims,idims,       'TOT_WATER',rvara=         wstorage)

!     ndims = 1
!     idims(1) = mwl
!     call shdf5_irec(ndims,idims,'TOT_RUNOFF',rvara=tot_runoff)

!     ndims = 2
!     idims(1) = 24
!     idims(2) = mwa

!     call shdf5_irec(ndims,idims,'PRESS_AVG24',rvara=press_avg24)
!     call shdf5_irec(ndims,idims,'TEMPK_AVG24',rvara=tempk_avg24)
!     call shdf5_irec(ndims,idims, 'SHUM_AVG24',rvara= shum_avg24)
!     call shdf5_irec(ndims,idims,  'VXE_AVG24',rvara=  vxe_avg24)
!     call shdf5_irec(ndims,idims,  'VYE_AVG24',rvara=  vye_avg24)
!     call shdf5_irec(ndims,idims,  'VZE_AVG24',rvara=  vze_avg24)

!     call shdf5_irec(ndims,idims,      'RSHORT_AVG24',rvara=      rshort_avg24)
!     call shdf5_irec(ndims,idims,  'RSHORT_TOP_AVG24',rvara=  rshort_top_avg24)
!     call shdf5_irec(ndims,idims,    'RSHORTUP_AVG24',rvara=    rshortup_avg24)
!     call shdf5_irec(ndims,idims,'RSHORTUP_TOP_AVG24',rvara=rshortup_top_avg24)
!     call shdf5_irec(ndims,idims,       'RLONG_AVG24',rvara=       rlong_avg24)
!     call shdf5_irec(ndims,idims,     'RLONGUP_AVG24',rvara=     rlongup_avg24)
!     call shdf5_irec(ndims,idims, 'RLONGUP_TOP_AVG24',rvara= rlongup_top_avg24)
!     call shdf5_irec(ndims,idims,     'LATFLUX_AVG24',rvara=     latflux_avg24)
!     call shdf5_irec(ndims,idims,    'SENSFLUX_AVG24',rvara=    sensflux_avg24)
!     call shdf5_irec(ndims,idims,     'ACCPMIC_TOT24',rvara=     accpmic_tot24)
!     call shdf5_irec(ndims,idims,     'ACCPCON_TOT24',rvara=     accpcon_tot24)

     call shdf5_close()
  
  else

   ! Month-average file does not exist, stop model.
   
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open month-average file:'
     write(io6,*) '!!!   '//trim(mavgfile)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'in read_mavg_vars'
   
  endif

  return
end subroutine read_mavg_vars

!===============================================================

subroutine read_davg_vars(davgfile)
  
  use mem_average_vars, only: &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     canltempk_davg, canltempk_dmin, canltempk_dmax, &
      vegtempk_davg,  vegtempk_dmin,  vegtempk_dmax, &
     soiltempk_davg, soiltempk_dmin, soiltempk_dmax, &
     canstempk_davg, canstempk_dmin, canstempk_dmax, &
        tempk_lf_davg,    tempk_lf_dmin,    tempk_lf_dmax, &
     cantempk_lf_davg, cantempk_lf_dmin, cantempk_lf_dmax, &
       sfluxt_lf_davg,   sfluxr_lf_davg,                   &
        tempk_sf_davg,    tempk_sf_dmin,    tempk_sf_dmax, &
     cantempk_sf_davg, cantempk_sf_dmin, cantempk_sf_dmax, &
       sfluxt_sf_davg,   sfluxr_sf_davg

  use misc_coms,  only: io6
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa
  use leaf_coms,  only: mwl
  use sea_coms,   only: mws
  use mem_sflux,  only: mseaflux, mlandflux

  implicit none

  character(*), intent(in) :: davgfile

  logical exans
  integer :: ndims, idims(2)
  character(len=32) :: varn

  inquire(file=davgfile,exist=exans)

  if (exans) then
  
! Day-average file exists.  Open, read, and close file.

     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(io6,*) 'Opening day-average file ', trim(davgfile)
     write(io6,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     call shdf5_open(trim(davgfile),'R')

 ! Read day-average variables
  
     ndims    = 1
     idims(1) = mwa

     call shdf5_irec(ndims,idims,  'PRESS_DAVG',rvara=  press_davg)
     call shdf5_irec(ndims,idims,    'VXE_DAVG',rvara=    vxe_davg)
     call shdf5_irec(ndims,idims,    'VYE_DAVG',rvara=    vye_davg)
     call shdf5_irec(ndims,idims,    'VZE_DAVG',rvara=    vze_davg)
     call shdf5_irec(ndims,idims, 'RSHORT_DAVG',rvara= rshort_davg)
     call shdf5_irec(ndims,idims,  'TEMPK_DAVG',rvara=  tempk_davg)
     call shdf5_irec(ndims,idims,  'TEMPK_DMIN',rvara=  tempk_dmin)
     call shdf5_irec(ndims,idims,  'TEMPK_DMAX',rvara=  tempk_dmax)
     call shdf5_irec(ndims,idims,'ACCPMIC_DTOT',rvara=accpmic_dtot)
     call shdf5_irec(ndims,idims,'ACCPCON_DTOT',rvara=accpcon_dtot)

     call shdf5_irec(ndims,idims,'PRESS_UL_DAVG',rvara=press_ul_davg)
     call shdf5_irec(ndims,idims,  'VXE_UL_DAVG',rvara=  vxe_ul_davg)
     call shdf5_irec(ndims,idims,  'VYE_UL_DAVG',rvara=  vye_ul_davg)
     call shdf5_irec(ndims,idims,  'VZE_UL_DAVG',rvara=  vze_ul_davg)

     idims(1) = mwl

     call shdf5_irec(ndims,idims,'CANLTEMPK_DAVG',rvara=canltempk_davg)
     call shdf5_irec(ndims,idims, 'VEGTEMPK_DAVG',rvara= vegtempk_davg)
     call shdf5_irec(ndims,idims,'SOILTEMPK_DAVG',rvara=soiltempk_davg)
     call shdf5_irec(ndims,idims,'CANLTEMPK_DMIN',rvara=canltempk_dmin)
     call shdf5_irec(ndims,idims,'CANLTEMPK_DMAX',rvara=canltempk_dmax)
     call shdf5_irec(ndims,idims, 'VEGTEMPK_DMIN',rvara= vegtempk_dmin)
     call shdf5_irec(ndims,idims, 'VEGTEMPK_DMAX',rvara= vegtempk_dmax)
     call shdf5_irec(ndims,idims,'SOILTEMPK_DMIN',rvara=soiltempk_dmin)
     call shdf5_irec(ndims,idims,'SOILTEMPK_DMAX',rvara=soiltempk_dmax)

     idims(1) = mws

     call shdf5_irec(ndims,idims,'CANSTEMPK_DAVG',rvara=canstempk_davg)
     call shdf5_irec(ndims,idims,'CANSTEMPK_DMIN',rvara=canstempk_dmin)
     call shdf5_irec(ndims,idims,'CANSTEMPK_DMAX',rvara=canstempk_dmax)
  
     idims(1) = mlandflux

     call shdf5_irec(ndims,idims,   'TEMPK_LF_DAVG',rvara=   tempk_lf_davg)
     call shdf5_irec(ndims,idims,   'TEMPK_LF_DMIN',rvara=   tempk_lf_dmin)
     call shdf5_irec(ndims,idims,   'TEMPK_LF_DMAX',rvara=   tempk_lf_dmax)
     call shdf5_irec(ndims,idims,'CANTEMPK_LF_DAVG',rvara=cantempk_lf_davg)
     call shdf5_irec(ndims,idims,'CANTEMPK_LF_DMIN',rvara=cantempk_lf_dmin)
     call shdf5_irec(ndims,idims,'CANTEMPK_LF_DMAX',rvara=cantempk_lf_dmax)
     call shdf5_irec(ndims,idims,  'SFLUXT_LF_DAVG',rvara=  sfluxt_lf_davg)
     call shdf5_irec(ndims,idims,  'SFLUXR_LF_DAVG',rvara=  sfluxr_lf_davg)

     idims(1) = mseaflux

     call shdf5_irec(ndims,idims,   'TEMPK_SF_DAVG',rvara=   tempk_sf_davg)
     call shdf5_irec(ndims,idims,   'TEMPK_SF_DMIN',rvara=   tempk_sf_dmin)
     call shdf5_irec(ndims,idims,   'TEMPK_SF_DMAX',rvara=   tempk_sf_dmax)
     call shdf5_irec(ndims,idims,'CANTEMPK_SF_DAVG',rvara=cantempk_sf_davg)
     call shdf5_irec(ndims,idims,'CANTEMPK_SF_DMIN',rvara=cantempk_sf_dmin)
     call shdf5_irec(ndims,idims,'CANTEMPK_SF_DMAX',rvara=cantempk_sf_dmax)
     call shdf5_irec(ndims,idims,  'SFLUXT_SF_DAVG',rvara=  sfluxt_sf_davg)
     call shdf5_irec(ndims,idims,  'SFLUXR_SF_DAVG',rvara=  sfluxr_sf_davg)

     call shdf5_close()
  else

   ! Day-average file does not exist, stop model.
   
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(io6,*) '!!!   Trying to open day-average file:'
     write(io6,*) '!!!   '//trim(davgfile)
     write(io6,*) '!!!   but it does not exist. The run is ended.'
     write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop 'in read_davg_vars'
   
  endif

  return
end subroutine read_davg_vars

!=============================================================================

subroutine makefnam8(fname,prefix,tinc8,iyr,imn,idy,itm,ftype,post,fmt)

! creates standard timestamped filename

  implicit none

  character(len=*), intent(out) :: fname
  character(len=*), intent(in) :: prefix,post,fmt
  real(kind=8), intent(in) :: tinc8
  integer, intent(in) :: iyr, imn, idy, itm
  character*(*), intent(in) :: ftype

  integer :: oyr, omn, ody, otm

  character(len=40) :: dstring

  if (tinc8 == 0.) then
     oyr = iyr ; omn = imn ; ody = idy ; otm = itm
  else
     call date_add_to8(iyr,imn,idy,itm,tinc8,'s',oyr,omn,ody,otm)
  endif

  write(dstring,100) '-',trim(ftype),'-',oyr,'-',omn,'-',ody,'-',otm
  100 format(3a,i4.4,a1,i2.2,a1,i2.2,a1,i6.6)

  fname = trim(prefix)//trim(dstring)
  if (post(1:1) /= '$') then
     fname = trim(fname)//'_'//trim(post)
  endif
  fname = trim(fname)//'.'//trim(fmt)

  return
end subroutine makefnam8
