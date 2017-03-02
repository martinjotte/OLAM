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

  use mem_grid,    only: mwa, lpw, zt, zm, mza, dzim, lpv, mva
  use mem_basic,   only: wc, sh_v, sh_w, tair, press, rho, vxe, vye, vze
  use consts_coms, only: alvl, cp
  use mem_turb,    only: sfluxr, sfluxt
  use misc_coms,   only: io6, dtlong, time8, isubdomain
  use mem_cuparm,  only: conprr
!obs  use mem_micro,   only: pcpgr
  use mem_radiate, only: rshort, rshort_top, rshortup_top, rlong, rlongup, &
                         rlongup_top, albedt
  use mem_ijtabs,  only: jtab_w, jtw_prog

  implicit none

  integer :: iw, k, dhr, kll, j, kw
  real(kind=8) :: dsec
  real :: windspeed
  integer, dimension(nz_avg-1) :: level_list = (/ 8, 13, 18, 26, 30, 34 /)

  npoints_mavg = npoints_mavg + 1

  dsec = mod(time8,86400.d0)
  dhr  = int(dsec / 3600.d0) + 1

  npoints_avg24(dhr) = npoints_avg24(dhr) + 1

! Horizontal loop over all prognostic W/T points

  !$omp parallel do private(iw,k,kll,windspeed)
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

       latflux_mavg(iw) =   latflux_mavg(iw) + sfluxr(iw) * alvl
      sensflux_mavg(iw) =  sensflux_mavg(iw) + sfluxt(iw) * cp
     windspeed_mavg(iw) = windspeed_mavg(iw) + windspeed

!obs     accpmic_mtot(iw) = accpmic_mtot(iw) + pcpgr(iw) * dtlong
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
          latflux_avg24(dhr,iw) =      latflux_avg24(dhr,iw) + sfluxr(iw) * alvl
         sensflux_avg24(dhr,iw) =     sensflux_avg24(dhr,iw) + sfluxt(iw) * cp
!obs          accpmic_tot24(dhr,iw) =      accpmic_tot24(dhr,iw) + pcpgr (iw) * dtlong
     if (allocated(conprr)) then
          accpcon_tot24(dhr,iw) =      accpcon_tot24(dhr,iw) + conprr(iw) * dtlong
     endif

  enddo
  !$omp end parallel do

end subroutine inc_mavg_vars

!========================================================================

subroutine inc_davg_vars()

  use mem_average_vars, only: &
     npoints_davg, nz_avg, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     airtempk_l_davg,  airtempk_l_dmin,  airtempk_l_dmax, &
     cantempk_l_davg,  cantempk_l_dmin,  cantempk_l_dmax, &
       vegtempk_davg,    vegtempk_dmin,    vegtempk_dmax, &
      soiltempk_davg,   soiltempk_dmin,   soiltempk_dmax, &
       sfluxt_l_davg,    sfluxr_l_davg,                   &
     airtempk_s_davg,  airtempk_s_dmin,  airtempk_s_dmax, &
     cantempk_s_davg,  cantempk_s_dmin,  cantempk_s_dmax, &
       sfluxt_s_davg,    sfluxr_s_davg

  use mem_ijtabs,  only: itabg_w, jtab_w, jtw_prog
  use mem_grid,    only: mwa, lpw, zt, zm, mza, dzim
  use mem_basic,   only: wc, sh_v, sh_w, tair, press, rho, vxe, vye, vze
  use consts_coms, only: alvl, cp
  use mem_turb,    only: sfluxr, sfluxt
  use misc_coms,   only: io6, dtlong, time8, isubdomain
  use mem_cuparm,  only: conprr
!obs  use mem_micro,   only: pcpgr
  use mem_radiate, only: rshort, rshort_top, rshortup_top, rlong, rlongup, &
                         rlongup_top, albedt

  use leaf_coms, only: nzg, mwl, slcpd
  use mem_leaf,  only: land, itab_wl
  use sea_coms,  only: mws
  use mem_sea,   only: sea, itab_ws
  use mem_para,  only: myrank

  implicit none

  integer :: j, iw, k, iwl, iws, kw
  real :: airtempk, cantempk, vegtempk, soiltempk, fracliq
  integer, dimension(nz_avg-1) :: level_list = (/ 8, 13, 18, 26, 30, 34 /)

  npoints_davg = npoints_davg + 1

  !$omp parallel

! Horizontal loop over all prognostic W/T points

  !$omp do private (iw,k)
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

!obs     accpmic_dtot(iw) = accpmic_dtot(iw) + pcpgr(iw) * dtlong
     if (allocated(conprr)) then
        accpcon_dtot(iw) = accpcon_dtot(iw) + conprr(iw) * dtlong
     endif

     press_ul_davg(iw) = press_ul_davg(iw) + real(press(level_list(5),iw))
       vxe_ul_davg(iw) =   vxe_ul_davg(iw) + vxe(level_list(5),iw)
       vye_ul_davg(iw) =   vye_ul_davg(iw) + vye(level_list(5),iw)
       vze_ul_davg(iw) =   vze_ul_davg(iw) + vze(level_list(5),iw)

  enddo
  !$omp end do

! Horizontal loop over all land cells

  !$omp do private(iw,kw,airtempk,cantempk,vegtempk,soiltempk,fracliq)
  do iwl = 2, mwl
     iw = itab_wl(iwl)%iw         ! global index
     kw = itab_wl(iwl)%kw

     ! If run is parallel, get local rank indices
     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

     airtempk = tair(kw,iw)
     cantempk = land%cantemp(iwl)
     vegtempk = land%veg_temp(iwl)

     call qwtk(land%soil_energy(nzg,iwl),       &
               land%soil_water(nzg,iwl)*1.e3,   &
               slcpd(land%ntext_soil(nzg,iwl)), &
               soiltempk, fracliq)

     airtempk_l_davg(iwl) = airtempk_l_davg(iwl) +  airtempk
     cantempk_l_davg(iwl) = cantempk_l_davg(iwl) +  cantempk
       vegtempk_davg(iwl) =   vegtempk_davg(iwl) +  vegtempk
      soiltempk_davg(iwl) =  soiltempk_davg(iwl) + soiltempk
       sfluxt_l_davg(iwl) =   sfluxt_l_davg(iwl) + land%sfluxt(iwl)
       sfluxr_l_davg(iwl) =   sfluxr_l_davg(iwl) + land%sfluxr(iwl)

     if (airtempk_l_dmin(iwl) > airtempk) airtempk_l_dmin(iwl) = airtempk
     if (airtempk_l_dmax(iwl) < airtempk) airtempk_l_dmax(iwl) = airtempk

     if (cantempk_l_dmin(iwl) > cantempk) cantempk_l_dmin(iwl) = cantempk
     if (cantempk_l_dmax(iwl) < cantempk) cantempk_l_dmax(iwl) = cantempk

     if (vegtempk_dmin(iwl) > vegtempk) vegtempk_dmin(iwl) = vegtempk
     if (vegtempk_dmax(iwl) < vegtempk) vegtempk_dmax(iwl) = vegtempk

     if (soiltempk_dmin(iwl) > soiltempk) soiltempk_dmin(iwl) = soiltempk
     if (soiltempk_dmax(iwl) < soiltempk) soiltempk_dmax(iwl) = soiltempk
  enddo
  !$omp end do nowait

! Horizontal loop over all sea cells

  !$omp do private(iw,kw,airtempk,cantempk)
  do iws = 2, mws
     iw = itab_ws(iws)%iw         ! global index
     kw = itab_ws(iws)%kw

     ! If run is parallel, get local rank indices
     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

     airtempk = tair(kw,iw)
     cantempk = sea%cantemp(iws)

     airtempk_s_davg(iws) = airtempk_s_davg(iws) + airtempk
     cantempk_s_davg(iws) = cantempk_s_davg(iws) + cantempk
       sfluxt_s_davg(iws) =   sfluxt_s_davg(iws) + sea%sfluxt(iws)
       sfluxr_s_davg(iws) =   sfluxr_s_davg(iws) + sea%sfluxr(iws)

     if (airtempk_s_dmin(iws) > airtempk) airtempk_s_dmin(iws) = airtempk
     if (airtempk_s_dmax(iws) < airtempk) airtempk_s_dmax(iws) = airtempk

     if (cantempk_s_dmin(iws) > cantempk) cantempk_s_dmin(iws) = cantempk
     if (cantempk_s_dmax(iws) < cantempk) cantempk_s_dmax(iws) = cantempk
  enddo
  !$omp end do nowait

  !$omp end parallel
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

  !$omp parallel do private(iw,ih)
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
  !$omp end parallel do

end subroutine norm_mavg_vars

!========================================================================

subroutine norm_davg_vars()

  use mem_average_vars, only: &
     npoints_davg, &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, tempk_davg, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
         vegtempk_davg, soiltempk_davg, &
     cantempk_l_davg, airtempk_l_davg, sfluxt_l_davg, sfluxr_l_davg, &
     cantempk_s_davg, airtempk_s_davg, sfluxt_s_davg, sfluxr_s_davg

  use mem_grid,  only: mwa
  use leaf_coms, only: mwl
  use mem_leaf,  only: itab_wl
  use sea_coms,  only: mws
  use mem_sea,   only: itab_ws
  use mem_ijtabs,only: jtab_w, jtw_prog
  use misc_coms, only: isubdomain
  use mem_para,  only: myrank

  implicit none

  integer :: iw, iwl, iws, j
  real :: ni

  ni = 1. / npoints_davg

  !$omp parallel

  !$omp do private(iw)
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
  !$omp end do nowait

  !$omp do
  do iwl = 2, mwl
     airtempk_l_davg(iwl) = airtempk_l_davg(iwl) * ni
     cantempk_l_davg(iwl) = cantempk_l_davg(iwl) * ni
       vegtempk_davg(iwl) =   vegtempk_davg(iwl) * ni
      soiltempk_davg(iwl) =  soiltempk_davg(iwl) * ni
       sfluxt_l_davg(iwl) =   sfluxt_l_davg(iwl) * ni
       sfluxr_l_davg(iwl) =   sfluxr_l_davg(iwl) * ni
  enddo
  !$omp end do nowait

  !$omp do
  do iws = 2, mws
     airtempk_s_davg(iws) = airtempk_s_davg(iws) * ni
     cantempk_s_davg(iws) = cantempk_s_davg(iws) * ni
       sfluxt_s_davg(iws) =   sfluxt_s_davg(iws) * ni
       sfluxr_s_davg(iws) =   sfluxr_s_davg(iws) * ni
  enddo
  !$omp end do nowait

  !$omp end parallel
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
                        imonth1, idate1, itime1
  use mem_para,   only: myrank, iwa_globe_primary, iwa_local_primary
  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa, nwa
  use leaf_coms,  only: mwl, nzg, nzs, dslz
  use mem_leaf,   only: land, itab_wl
  use mem_basic,  only: rho
  use mem_ijtabs, only: itabg_w

  implicit none

  integer, intent(in) :: outyear, outmonth
  character(len=128) :: hnamel
  integer :: lenl, iwl, k, kw, iw
  logical exans
  character(len=32) :: varn
  integer :: ndims, idims(2), month_use, year_use
  real :: wstorage(mwl)
  real :: tot_soil_water
  integer, pointer :: ilpts(:), igpts(:)
  integer :: nglobe
  
  ! Compute total water

  do iwl = 2, mwl
     iw = itab_wl(iwl)%iw         ! global index
     kw = itab_wl(iwl)%kw

     ! If run is parallel, get local rank indices
     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

     wstorage(iwl) = sum(land%sfcwater_mass(1:nzs,iwl)) +   &
          land%veg_water(iwl) +   &
          rho(kw,iw) * land%can_depth(iwl) * land%canshv(iwl)
     tot_soil_water = 0.0
     do k = 1, nzg
        tot_soil_water = tot_soil_water +   &
             dslz(k) * 1000. * land%soil_water(k,iwl)
     enddo
     wstorage(iwl) = wstorage(iwl) + tot_soil_water
  enddo

! Determine output file name and open file

  month_use = outmonth - 1
  year_use = outyear
  if(month_use == 0)then
     month_use = 12
     year_use = year_use - 1
  endif

  call makefnam8(hnamel,hfilepref,0.d0,year_use,month_use,1,  &
       0,'M','$','h5')  ! '$' argument suppresses grid# encoding

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
  
  ilpts => iwa_local_primary
  igpts => iwa_globe_primary
  nglobe = nwa

  ndims    = 2
  idims(1) = nz_avg
  idims(2) = mwa

  call shdf5_orec(ndims,idims,'PRESS_MAVG',rvar2=press_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'RHO_MAVG',rvar2=  rho_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'TEMPK_MAVG',rvar2=tempk_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SH_V_MAVG',rvar2= sh_v_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SH_W_MAVG',rvar2= sh_w_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,   'WC_MAVG',rvar2=   wc_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VXE_MAVG',rvar2=  vxe_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VYE_MAVG',rvar2=  vye_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VZE_MAVG',rvar2=  vze_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  ndims    = 1
  idims(1) = mwa
  ilpts => iwa_local_primary
  igpts => iwa_globe_primary
  nglobe = nwa

  call shdf5_orec(ndims,idims,      'RSHORT_MAVG',rvar1=      rshort_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'RSHORT_TOP_MAVG',rvar1=  rshort_top_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'RSHORTUP_MAVG',rvar1=    rshortup_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'RSHORTUP_TOP_MAVG',rvar1=rshortup_top_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,       'RLONG_MAVG',rvar1=       rlong_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,     'RLONGUP_MAVG',rvar1=     rlongup_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'RLONGUP_TOP_MAVG',rvar1= rlongup_top_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,     'LATFLUX_MAVG',rvar1=     latflux_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'SENSFLUX_MAVG',rvar1=    sensflux_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,   'WINDSPEED_MAVG',rvar1=   windspeed_mavg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,     'ACCPMIC_MTOT',rvar1=     accpmic_mtot, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,     'ACCPCON_MTOT',rvar1=     accpcon_mtot, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,        'TOT_WATER',rvar1=         wstorage, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)


!  ndims = 1
!  idims(1) = mwl
!  ilpts => iwl_local_primary
!  igpts => iwl_globe_primary
!  nglobe = nwl
!
!  call shdf5_orec(ndims,idims,'TOT_RUNOFF',rvar1=tot_runoff)

!  ndims = 2
!  idims(1) = 24
!  idims(2) = mwa
!  ilpts => iwa_local_primary
!  igpts => iwa_globe_primary
!  nglobe = nwa
!
!  call shdf5_orec(ndims,idims,'PRESS_AVG24',rvar2=press_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,'TEMPK_AVG24',rvar2=tempk_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims, 'SHUM_AVG24',rvar2= shum_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,  'VXE_AVG24',rvar2=  vxe_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,  'VYE_AVG24',rvar2=  vye_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,  'VZE_AVG24',rvar2=  vze_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

!  call shdf5_orec(ndims,idims,      'RSHORT_AVG24',rvar2=      rshort_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,  'RSHORT_TOP_AVG24',rvar2=  rshort_top_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,    'RSHORTUP_AVG24',rvar2=    rshortup_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,'RSHORTUP_TOP_AVG24',rvar2=rshortup_top_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,       'RLONG_AVG24',rvar2=       rlong_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,     'RLONGUP_AVG24',rvar2=     rlongup_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims, 'RLONGUP_TOP_AVG24',rvar2= rlongup_top_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,     'LATFLUX_AVG24',rvar2=     latflux_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,    'SENSFLUX_AVG24',rvar2=    sensflux_avg24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,     'ACCPMIC_TOT24',rvar2=     accpmic_tot24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
!  call shdf5_orec(ndims,idims,     'ACCPCON_TOT24',rvar2=     accpcon_tot24, &
!       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  call shdf5_close()
  
end subroutine write_mavg_vars

!===============================================================

subroutine write_davg_vars(outyear,outmonth,outdate)
  
  use mem_average_vars, only: &
     press_davg, vxe_davg, vye_davg, vze_davg, rshort_davg, &
     tempk_davg, tempk_dmin, tempk_dmax, accpmic_dtot, accpcon_dtot, &
     press_ul_davg, vxe_ul_davg, vye_ul_davg, vze_ul_davg, &
     airtempk_l_davg, airtempk_l_dmin, airtempk_l_dmax, &
     cantempk_l_davg, cantempk_l_dmin, cantempk_l_dmax, &
       vegtempk_davg,   vegtempk_dmin,   vegtempk_dmax, &
      soiltempk_davg,  soiltempk_dmin,  soiltempk_dmax, &
       sfluxt_l_davg,   sfluxr_l_davg,                  &
     airtempk_s_davg, airtempk_s_dmin, airtempk_s_dmax, &
     cantempk_s_davg, cantempk_s_dmin, cantempk_s_dmax, &
       sfluxt_s_davg,   sfluxr_s_davg

  use hdf5_utils, only: shdf5_orec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa, nwa
  use misc_coms,  only: iparallel, hfilepref, iclobber, io6, time8, iyear1, &
                        imonth1, idate1, itime1
  use mem_para,   only: myrank, iwa_globe_primary, iwa_local_primary, &
                                iwl_globe_primary, iwl_local_primary, &
                                iws_globe_primary, iws_local_primary

  use leaf_coms,  only: mwl, nwl
  use sea_coms,   only: mws, nws

  implicit none

  integer, intent(in) :: outyear, outmonth, outdate
  character(len=128) :: hnamel
  integer :: lenl
  logical exans
  character(len=32) :: varn
  integer :: ndims,idims(2), month_use, year_use, date_use
  integer, pointer :: ilpts(:), igpts(:)
  integer :: nglobe

! Determine output file name and open file

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
       0,'D','$','h5')  ! '$' argument suppresses grid# encoding

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

  ndims    = 1
  idims(1) = mwa

  ilpts => iwa_local_primary
  igpts => iwa_globe_primary
  nglobe = nwa

  call shdf5_orec(ndims,idims,  'PRESS_DAVG',rvar1=  press_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'VXE_DAVG',rvar1=    vxe_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'VYE_DAVG',rvar1=    vye_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,    'VZE_DAVG',rvar1=    vze_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'RSHORT_DAVG',rvar1= rshort_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'TEMPK_DAVG',rvar1=  tempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'TEMPK_DMIN',  rvar1=tempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'TEMPK_DMAX',  rvar1=tempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'ACCPMIC_DTOT',rvar1=accpmic_dtot, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'ACCPCON_DTOT',rvar1=accpcon_dtot, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  call shdf5_orec(ndims,idims,'PRESS_UL_DAVG',rvar1=press_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VXE_UL_DAVG',rvar1=  vxe_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VYE_UL_DAVG',rvar1=  vye_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VZE_UL_DAVG',rvar1=  vze_ul_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  idims(1) = mwl

  ilpts => iwl_local_primary
  igpts => iwl_globe_primary
  nglobe = nwl

  call shdf5_orec(ndims,idims,'AIRTEMPK_L_DAVG',rvar1=airtempk_l_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'AIRTEMPK_L_DMIN',rvar1=airtempk_l_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'AIRTEMPK_L_DMAX',rvar1=airtempk_l_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_L_DAVG',rvar1=cantempk_l_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_L_DMIN',rvar1=cantempk_l_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_L_DMAX',rvar1=cantempk_l_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VEGTEMPK_DAVG',rvar1=  vegtempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VEGTEMPK_DMIN',rvar1=  vegtempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'VEGTEMPK_DMAX',rvar1=  vegtempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SOILTEMPK_DAVG',rvar1= soiltempk_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SOILTEMPK_DMIN',rvar1= soiltempk_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims, 'SOILTEMPK_DMAX',rvar1= soiltempk_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'SFLUXT_L_DAVG',rvar1=  sfluxt_l_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'SFLUXR_L_DAVG',rvar1=  sfluxr_l_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)

  idims(1) = mws

  ilpts => iws_local_primary
  igpts => iws_globe_primary
  nglobe = nws

  call shdf5_orec(ndims,idims,'AIRTEMPK_S_DAVG',rvar1=airtempk_s_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'AIRTEMPK_S_DMIN',rvar1=airtempk_s_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'AIRTEMPK_S_DMAX',rvar1=airtempk_s_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_S_DAVG',rvar1=cantempk_s_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_S_DMIN',rvar1=cantempk_s_dmin, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,'CANTEMPK_S_DMAX',rvar1=cantempk_s_dmax, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'SFLUXT_S_DAVG',rvar1=  sfluxt_s_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  call shdf5_orec(ndims,idims,  'SFLUXR_S_DAVG',rvar1=  sfluxr_s_davg, &
       lpoints=ilpts, gpoints = igpts, nglobe=nglobe)
  
  call shdf5_close()

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

  use misc_coms,  only: io6
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa, nwa
  use mem_ijtabs, only: itab_w
  use mem_leaf,   only: itab_wl

  implicit none

  character(*), intent(in) :: mavgfile

  logical :: exans
  integer :: ndims, idims(2)
  character(len=32) :: varn
  integer, pointer :: ilocal(:)

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
     ilocal => itab_w(:)%iwglobe

     call shdf5_irec(ndims,idims,'PRESS_MAVG',rvar2=press_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,  'RHO_MAVG',rvar2=  rho_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,'TEMPK_MAVG',rvar2=tempk_mavg, points=ilocal)
     call shdf5_irec(ndims,idims, 'SH_V_MAVG',rvar2= sh_v_mavg, points=ilocal)
     call shdf5_irec(ndims,idims, 'SH_W_MAVG',rvar2= sh_w_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,   'WC_MAVG',rvar2=   wc_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VXE_MAVG',rvar2=  vxe_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VYE_MAVG',rvar2=  vye_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VZE_MAVG',rvar2=  vze_mavg, points=ilocal)

     ndims    = 1
     idims(1) = mwa
     ilocal => itab_w(:)%iwglobe

     call shdf5_irec(ndims,idims,      'RSHORT_MAVG',rvar1=      rshort_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,  'RSHORT_TOP_MAVG',rvar1=  rshort_top_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,    'RSHORTUP_MAVG',rvar1=    rshortup_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,'RSHORTUP_TOP_MAVG',rvar1=rshortup_top_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,       'RLONG_MAVG',rvar1=       rlong_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,     'RLONGUP_MAVG',rvar1=     rlongup_mavg, points=ilocal)
     call shdf5_irec(ndims,idims, 'RLONGUP_TOP_MAVG',rvar1= rlongup_top_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,     'LATFLUX_MAVG',rvar1=     latflux_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,    'SENSFLUX_MAVG',rvar1=    sensflux_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,   'WINDSPEED_MAVG',rvar1=   windspeed_mavg, points=ilocal)
     call shdf5_irec(ndims,idims,     'ACCPMIC_MTOT',rvar1=     accpmic_mtot, points=ilocal)
     call shdf5_irec(ndims,idims,     'ACCPCON_MTOT',rvar1=     accpcon_mtot, points=ilocal)
!    call shdf5_irec(ndims,idims,        'TOT_WATER',rvar1=         wstorage, points=ilocal)

!     ndims = 1
!     idims(1) = mwl
!     ilocal => itab_wl(:)%iwglobe
!
!     call shdf5_irec(ndims,idims,'TOT_RUNOFF',rvar1=tot_runoff, points=ilocal)
!
!     ndims = 2
!     idims(1) = 24
!     idims(2) = mwa
!     ilocal => itab_w(:)%iwglobe
!
!     call shdf5_irec(ndims,idims,'PRESS_AVG24',rvar2=press_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,'TEMPK_AVG24',rvar2=tempk_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims, 'SHUM_AVG24',rvar2= shum_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,  'VXE_AVG24',rvar2=  vxe_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,  'VYE_AVG24',rvar2=  vye_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,  'VZE_AVG24',rvar2=  vze_avg24, points=ilocal)
!
!     call shdf5_irec(ndims,idims,      'RSHORT_AVG24',rvar2=      rshort_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,  'RSHORT_TOP_AVG24',rvar2=  rshort_top_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,    'RSHORTUP_AVG24',rvar2=    rshortup_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,'RSHORTUP_TOP_AVG24',rvar2=rshortup_top_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,       'RLONG_AVG24',rvar2=       rlong_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,     'RLONGUP_AVG24',rvar2=     rlongup_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims, 'RLONGUP_TOP_AVG24',rvar2= rlongup_top_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,     'LATFLUX_AVG24',rvar2=     latflux_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,    'SENSFLUX_AVG24',rvar2=    sensflux_avg24, points=ilocal)
!     call shdf5_irec(ndims,idims,     'ACCPMIC_TOT24',rvar2=     accpmic_tot24, points=ilocal)
!     call shdf5_irec(ndims,idims,     'ACCPCON_TOT24',rvar2=     accpcon_tot24, points=ilocal)

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
     airtempk_l_davg, airtempk_l_dmin, airtempk_l_dmax, &
     cantempk_l_davg, cantempk_l_dmin, cantempk_l_dmax, &
       vegtempk_davg,   vegtempk_dmin,   vegtempk_dmax, &
      soiltempk_davg,  soiltempk_dmin,  soiltempk_dmax, &
       sfluxt_l_davg,   sfluxr_l_davg,                  &
     airtempk_s_davg, airtempk_s_dmin, airtempk_s_dmax, &
     cantempk_s_davg, cantempk_s_dmin, cantempk_s_dmax, &
       sfluxt_s_davg,   sfluxr_s_davg

  use misc_coms,  only: io6
  use hdf5_utils, only: shdf5_irec, shdf5_open, shdf5_close
  use mem_grid,   only: mwa
  use leaf_coms,  only: mwl
  use sea_coms,   only: mws
  use mem_ijtabs, only: itabg_w, itab_w, itab_v, itab_m
  use mem_leaf,   only: itab_wl
  use mem_sea,    only: itab_ws

  implicit none

  character(*), intent(in) :: davgfile

  logical :: exans
  integer :: ndims, idims(2)
  character(len=32) :: varn
  integer, pointer :: ilocal(:)

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
     ilocal => itab_w(:)%iwglobe

     call shdf5_irec(ndims,idims,  'PRESS_DAVG',rvar1=  press_davg, points=ilocal)
     call shdf5_irec(ndims,idims,    'VXE_DAVG',rvar1=    vxe_davg, points=ilocal)
     call shdf5_irec(ndims,idims,    'VYE_DAVG',rvar1=    vye_davg, points=ilocal)
     call shdf5_irec(ndims,idims,    'VZE_DAVG',rvar1=    vze_davg, points=ilocal)
     call shdf5_irec(ndims,idims, 'RSHORT_DAVG',rvar1= rshort_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'TEMPK_DAVG',rvar1=  tempk_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'TEMPK_DMIN',rvar1=  tempk_dmin, points=ilocal)
     call shdf5_irec(ndims,idims,  'TEMPK_DMAX',rvar1=  tempk_dmax, points=ilocal)
     call shdf5_irec(ndims,idims,'ACCPMIC_DTOT',rvar1=accpmic_dtot, points=ilocal)
     call shdf5_irec(ndims,idims,'ACCPCON_DTOT',rvar1=accpcon_dtot, points=ilocal)

     call shdf5_irec(ndims,idims,'PRESS_UL_DAVG',rvar1=press_ul_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VXE_UL_DAVG',rvar1=  vxe_ul_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VYE_UL_DAVG',rvar1=  vye_ul_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VZE_UL_DAVG',rvar1=  vze_ul_davg, points=ilocal)

     idims(1) = mwl
     ilocal => itab_wl(:)%iwglobe

     call shdf5_irec(ndims,idims,'AIRTEMPK_L_DAVG',rvar1=airtempk_l_davg, points=ilocal)
     call shdf5_irec(ndims,idims,'AIRTEMPK_L_DMIN',rvar1=airtempk_l_dmin, points=ilocal)
     call shdf5_irec(ndims,idims,'AIRTEMPK_L_DMAX',rvar1=airtempk_l_dmax, points=ilocal)
     call shdf5_irec(ndims,idims,'CANTEMPK_L_DAVG',rvar1=cantempk_l_davg, points=ilocal)
     call shdf5_irec(ndims,idims,'CANTEMPK_L_DMIN',rvar1=cantempk_l_dmin, points=ilocal)
     call shdf5_irec(ndims,idims,'CANTEMPK_L_DMAX',rvar1=cantempk_l_dmax, points=ilocal)
     call shdf5_irec(ndims,idims,  'VEGTEMPK_DAVG',rvar1=  vegtempk_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'VEGTEMPK_DMIN',rvar1=  vegtempk_dmin, points=ilocal)
     call shdf5_irec(ndims,idims,  'VEGTEMPK_DMAX',rvar1=  vegtempk_dmax, points=ilocal)
     call shdf5_irec(ndims,idims, 'SOILTEMPK_DAVG',rvar1= soiltempk_davg, points=ilocal)
     call shdf5_irec(ndims,idims, 'SOILTEMPK_DMIN',rvar1= soiltempk_dmin, points=ilocal)
     call shdf5_irec(ndims,idims, 'SOILTEMPK_DMAX',rvar1= soiltempk_dmax, points=ilocal)
     call shdf5_irec(ndims,idims,  'SFLUXT_L_DAVG',rvar1=  sfluxt_l_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'SFLUXR_L_DAVG',rvar1=  sfluxr_l_davg, points=ilocal)

     idims(1) = mws
     ilocal => itab_ws(:)%iwglobe

     call shdf5_irec(ndims,idims,'AIRTEMPK_S_DAVG',rvar1=airtempk_s_davg, points=ilocal)
     call shdf5_irec(ndims,idims,'AIRTEMPK_S_DMIN',rvar1=airtempk_s_dmin, points=ilocal)
     call shdf5_irec(ndims,idims,'AIRTEMPK_S_DMAX',rvar1=airtempk_s_dmax, points=ilocal)
     call shdf5_irec(ndims,idims,'CANTEMPK_S_DAVG',rvar1=cantempk_s_davg, points=ilocal)
     call shdf5_irec(ndims,idims,'CANTEMPK_S_DMIN',rvar1=cantempk_s_dmin, points=ilocal)
     call shdf5_irec(ndims,idims,'CANTEMPK_S_DMAX',rvar1=cantempk_s_dmax, points=ilocal)
     call shdf5_irec(ndims,idims,  'SFLUXT_S_DAVG',rvar1=  sfluxt_s_davg, points=ilocal)
     call shdf5_irec(ndims,idims,  'SFLUXR_S_DAVG',rvar1=  sfluxr_s_davg, points=ilocal)

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
