Module mem_average_vars

  use mem_ijtabs, only: itab_w
  use misc_coms,  only: ngrids

  ! Number of time points to average over in whole month

  integer :: npoints_mavg

  ! Number of hourly-monthly time points

  integer :: npoints_avg24(24)

  ! Number of time points to average over in a day

  integer :: npoints_davg

  ! Number of levels on which to average the 3-d variables

  integer, parameter :: nz_avg = 7

  ! Monthly variables

  real, allocatable :: press_mavg(:,:)
  real, allocatable ::   rho_mavg(:,:)
  real, allocatable :: tempk_mavg(:,:)
  real, allocatable ::  sh_v_mavg(:,:)
  real, allocatable ::  sh_w_mavg(:,:)
  real, allocatable ::    wc_mavg(:,:)
  real, allocatable ::   vxe_mavg(:,:)
  real, allocatable ::   vye_mavg(:,:)
  real, allocatable ::   vze_mavg(:,:)

  real, allocatable ::       rshort_mavg(:)
  real, allocatable ::   rshort_top_mavg(:)
  real, allocatable ::     rshortup_mavg(:)
  real, allocatable :: rshortup_top_mavg(:)
  real, allocatable ::        rlong_mavg(:)
  real, allocatable ::      rlongup_mavg(:)
  real, allocatable ::  rlongup_top_mavg(:)
  real, allocatable ::      latflux_mavg(:)
  real, allocatable ::     sensflux_mavg(:)
  real, allocatable ::    windspeed_mavg(:)

  real, allocatable :: accpmic_mtot(:)
  real, allocatable :: accpcon_mtot(:)

  ! Monthly-average diurnal cycle

  real, allocatable :: press_avg24(:,:)
  real, allocatable :: tempk_avg24(:,:)
  real, allocatable ::  shum_avg24(:,:)
  real, allocatable ::   vxe_avg24(:,:)
  real, allocatable ::   vye_avg24(:,:)
  real, allocatable ::   vze_avg24(:,:)

  real, allocatable ::       rshort_avg24(:,:)
  real, allocatable ::   rshort_top_avg24(:,:)
  real, allocatable ::     rshortup_avg24(:,:)
  real, allocatable :: rshortup_top_avg24(:,:)
  real, allocatable ::        rlong_avg24(:,:)
  real, allocatable ::      rlongup_avg24(:,:)
  real, allocatable ::  rlongup_top_avg24(:,:)
  real, allocatable ::      latflux_avg24(:,:)
  real, allocatable ::     sensflux_avg24(:,:)

  real, allocatable :: accpmic_tot24(:,:)
  real, allocatable :: accpcon_tot24(:,:)

  ! Daily values - atmosphere

  real, allocatable ::  press_davg(:)
  real, allocatable ::    vxe_davg(:)
  real, allocatable ::    vye_davg(:)
  real, allocatable ::    vze_davg(:)
  real, allocatable :: rshort_davg(:)

  real, allocatable :: tempk_davg(:)
  real, allocatable :: tempk_dmin(:)
  real, allocatable :: tempk_dmax(:)

  real, allocatable :: accpmic_dtot(:)
  real, allocatable :: accpcon_dtot(:)

  real, allocatable :: press_ul_davg(:)
  real, allocatable ::   vxe_ul_davg(:)
  real, allocatable ::   vye_ul_davg(:)
  real, allocatable ::   vze_ul_davg(:)

  ! Daily values - land

  real, allocatable :: canltempk_davg(:)
  real, allocatable :: canltempk_dmin(:)
  real, allocatable :: canltempk_dmax(:)
  real, allocatable ::  vegtempk_davg(:)
  real, allocatable ::  vegtempk_dmin(:)
  real, allocatable ::  vegtempk_dmax(:)
  real, allocatable :: soiltempk_davg(:)
  real, allocatable :: soiltempk_dmin(:)
  real, allocatable :: soiltempk_dmax(:)

  ! Daily values - sea

  real, allocatable :: canstempk_davg(:)
  real, allocatable :: canstempk_dmin(:)
  real, allocatable :: canstempk_dmax(:)

  ! Daily values - landflux

  real, allocatable :: tempk_lf_davg(:)
  real, allocatable :: tempk_lf_dmin(:)
  real, allocatable :: tempk_lf_dmax(:)

  real, allocatable :: cantempk_lf_davg(:)
  real, allocatable :: cantempk_lf_dmin(:)
  real, allocatable :: cantempk_lf_dmax(:)

  real, allocatable :: sfluxt_lf_davg(:)
  real, allocatable :: sfluxr_lf_davg(:)

  ! Daily values - seaflux

  real, allocatable :: tempk_sf_davg(:)
  real, allocatable :: tempk_sf_dmin(:)
  real, allocatable :: tempk_sf_dmax(:)

  real, allocatable :: cantempk_sf_davg(:)
  real, allocatable :: cantempk_sf_dmin(:)
  real, allocatable :: cantempk_sf_dmax(:)

  real, allocatable :: sfluxt_sf_davg(:)
  real, allocatable :: sfluxr_sf_davg(:)

Contains

!===============================================================================

  subroutine alloc_average_vars(mwa, mwl, mws, mlandflux, mseaflux)

    use oname_coms,  only: nl

    implicit none

    integer, intent(in) :: mwa, mwl, mws, mlandflux, mseaflux

    if (nl%ioutput_mavg == 1) then

       allocate(press_mavg(nz_avg,mwa))
       allocate(rho_mavg  (nz_avg,mwa))
       allocate(tempk_mavg(nz_avg,mwa))
       allocate(sh_v_mavg (nz_avg,mwa))
       allocate(sh_w_mavg (nz_avg,mwa))
       allocate(wc_mavg   (nz_avg,mwa))
       allocate(vxe_mavg  (nz_avg,mwa))
       allocate(vye_mavg  (nz_avg,mwa))
       allocate(vze_mavg  (nz_avg,mwa))

       allocate(      rshort_mavg(mwa))
       allocate(  rshort_top_mavg(mwa))
       allocate(    rshortup_mavg(mwa))
       allocate(rshortup_top_mavg(mwa))
       allocate(       rlong_mavg(mwa))
       allocate(     rlongup_mavg(mwa))
       allocate( rlongup_top_mavg(mwa))
       allocate(     latflux_mavg(mwa))
       allocate(    sensflux_mavg(mwa))
       allocate(   windspeed_mavg(mwa))

       allocate(accpmic_mtot(mwa))
       allocate(accpcon_mtot(mwa))

       allocate(press_avg24(24,mwa))
       allocate(tempk_avg24(24,mwa))
       allocate( shum_avg24(24,mwa))
       allocate(  vxe_avg24(24,mwa))
       allocate(  vye_avg24(24,mwa))
       allocate(  vze_avg24(24,mwa))

       allocate(      rshort_avg24(24,mwa))
       allocate(  rshort_top_avg24(24,mwa))
       allocate(    rshortup_avg24(24,mwa))
       allocate(rshortup_top_avg24(24,mwa))
       allocate(       rlong_avg24(24,mwa))
       allocate(     rlongup_avg24(24,mwa))
       allocate( rlongup_top_avg24(24,mwa))
       allocate(     latflux_avg24(24,mwa))
       allocate(    sensflux_avg24(24,mwa))

       allocate(accpmic_tot24(24,mwa))
       allocate(accpcon_tot24(24,mwa))

    endif

  ! Daily values - atmosphere

    if (nl%ioutput_davg == 1) then

       allocate( press_davg(mwa))
       allocate(   vxe_davg(mwa))
       allocate(   vye_davg(mwa))
       allocate(   vze_davg(mwa))
       allocate(rshort_davg(mwa))

       allocate(tempk_davg(mwa))
       allocate(tempk_dmin(mwa))
       allocate(tempk_dmax(mwa))

       allocate(accpmic_dtot(mwa))
       allocate(accpcon_dtot(mwa))

       allocate(press_ul_davg(mwa))
       allocate(  vxe_ul_davg(mwa))
       allocate(  vye_ul_davg(mwa))
       allocate(  vze_ul_davg(mwa))

  ! Daily values - land

       allocate(canltempk_davg(mwl))
       allocate(canltempk_dmin(mwl))
       allocate(canltempk_dmax(mwl))
       allocate( vegtempk_davg(mwl))
       allocate( vegtempk_dmin(mwl))
       allocate( vegtempk_dmax(mwl))
       allocate(soiltempk_davg(mwl))
       allocate(soiltempk_dmin(mwl))
       allocate(soiltempk_dmax(mwl))

  ! Daily values - sea

       allocate(canstempk_davg(mws))
       allocate(canstempk_dmin(mws))
       allocate(canstempk_dmax(mws))

  ! Daily values - landflux

       allocate(tempk_lf_davg(mlandflux))
       allocate(tempk_lf_dmin(mlandflux))
       allocate(tempk_lf_dmax(mlandflux))

       allocate(cantempk_lf_davg(mlandflux))
       allocate(cantempk_lf_dmin(mlandflux))
       allocate(cantempk_lf_dmax(mlandflux))

       allocate(sfluxt_lf_davg(mlandflux))
       allocate(sfluxr_lf_davg(mlandflux))

  ! Daily values - seaflux

       allocate(tempk_sf_davg(mseaflux))
       allocate(tempk_sf_dmin(mseaflux))
       allocate(tempk_sf_dmax(mseaflux))

       allocate(cantempk_sf_davg(mseaflux))
       allocate(cantempk_sf_dmin(mseaflux))
       allocate(cantempk_sf_dmax(mseaflux))

       allocate(sfluxt_sf_davg(mseaflux))
       allocate(sfluxr_sf_davg(mseaflux))

    endif

    return
  end subroutine alloc_average_vars

!===============================================================================

  subroutine reset_mavg_vars()
    
    implicit none

    npoints_mavg = 0
    npoints_avg24(:) = 0

    press_mavg(:,:) = 0.0
      rho_mavg(:,:) = 0.0
    tempk_mavg(:,:) = 0.0
     sh_v_mavg(:,:) = 0.0
     sh_w_mavg(:,:) = 0.0
       wc_mavg(:,:) = 0.0
      vxe_mavg(:,:) = 0.0
      vye_mavg(:,:) = 0.0
      vze_mavg(:,:) = 0.0

          rshort_mavg(:) = 0.0
      rshort_top_mavg(:) = 0.0
        rshortup_mavg(:) = 0.0
    rshortup_top_mavg(:) = 0.0
           rlong_mavg(:) = 0.0
         rlongup_mavg(:) = 0.0
     rlongup_top_mavg(:) = 0.0
         latflux_mavg(:) = 0.0
        sensflux_mavg(:) = 0.0
       windspeed_mavg(:) = 0.0

    accpmic_mtot(:) = 0.0
    accpcon_mtot(:) = 0.0

  ! Monthly-average diurnal cycle

    press_avg24(:,:) = 0.0
    tempk_avg24(:,:) = 0.0
     shum_avg24(:,:) = 0.0
      vxe_avg24(:,:) = 0.0
      vye_avg24(:,:) = 0.0
      vze_avg24(:,:) = 0.0

          rshort_avg24(:,:) = 0.0
      rshort_top_avg24(:,:) = 0.0
        rshortup_avg24(:,:) = 0.0
    rshortup_top_avg24(:,:) = 0.0
           rlong_avg24(:,:) = 0.0
         rlongup_avg24(:,:) = 0.0
     rlongup_top_avg24(:,:) = 0.0
         latflux_avg24(:,:) = 0.0
        sensflux_avg24(:,:) = 0.0

    accpmic_tot24(:,:) = 0.0
    accpcon_tot24(:,:) = 0.0

    return
  end subroutine reset_mavg_vars

!===============================================================================

  subroutine reset_davg_vars()
    
    implicit none

    npoints_davg = 0

     press_davg(:) = 0.0
       vxe_davg(:) = 0.0
       vye_davg(:) = 0.0
       vze_davg(:) = 0.0
    rshort_davg(:) = 0.0

    tempk_davg(:) = 0.0
    tempk_dmin(:) = 10000.0
    tempk_dmax(:) = 0.0

    accpmic_dtot(:) = 0.0
    accpcon_dtot(:) = 0.0

    press_ul_davg(:) = 0.0
      vxe_ul_davg(:) = 0.0
      vye_ul_davg(:) = 0.0
      vze_ul_davg(:) = 0.0

  ! Daily values - land

    canltempk_davg(:) = 0.0
    canltempk_dmin(:) = 10000.0
    canltempk_dmax(:) = 0.0
     vegtempk_davg(:) = 0.0
     vegtempk_dmin(:) = 10000.0
     vegtempk_dmax(:) = 0.0
    soiltempk_davg(:) = 0.0
    soiltempk_dmin(:) = 10000.0
    soiltempk_dmax(:) = 0.0

  ! Daily values - sea

    canstempk_davg(:) = 0.0
    canstempk_dmin(:) = 10000.0
    canstempk_dmax(:) = 0.0

  ! Daily values - landflux

  tempk_lf_davg(:) = 0.0
  tempk_lf_dmin(:) = 10000.0
  tempk_lf_dmax(:) = 0.0

  cantempk_lf_davg(:) = 0.0
  cantempk_lf_dmin(:) = 10000.0
  cantempk_lf_dmax(:) = 0.0

  sfluxt_lf_davg(:) = 0.0
  sfluxr_lf_davg(:) = 0.0

  ! Daily values - seaflux

  tempk_sf_davg(:) = 0.0
  tempk_sf_dmin(:) = 10000.0
  tempk_sf_dmax(:) = 0.0

  cantempk_sf_davg(:) = 0.0
  cantempk_sf_dmin(:) = 10000.0
  cantempk_sf_dmax(:) = 0.0

  sfluxt_sf_davg(:) = 0.0
  sfluxr_sf_davg(:) = 0.0

    return
  end subroutine reset_davg_vars

End Module mem_average_vars
