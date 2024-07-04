Module mem_flux_accum

  use consts_coms, only: r8
  implicit none

  private :: r8

! Atmosphere grid arrays

  real(r8), allocatable ::       rshort_accum(:)
  real(r8), allocatable ::     rshortup_accum(:)
  real(r8), allocatable ::        rlong_accum(:)
  real(r8), allocatable ::      rlongup_accum(:)
  real(r8), allocatable ::   rshort_top_accum(:)
  real(r8), allocatable :: rshortup_top_accum(:)
  real(r8), allocatable ::  rlongup_top_accum(:)

  real(r8), allocatable ::          vc_accum(:,:)
  real(r8), allocatable ::          wc_accum(:,:)
  real(r8), allocatable ::       press_accum(:,:)
  real(r8), allocatable ::        tair_accum(:,:)
  real(r8), allocatable ::        rr_v_accum(:,:)
  real(r8), allocatable :: latheat_liq_accum(:,:) ! filled in omic_driv
  real(r8), allocatable :: latheat_ice_accum(:,:) ! filled in omic_driv

! Surface grid arrays

  real(r8), allocatable ::     vels_accum(:)
  real(r8), allocatable ::  airtemp_accum(:)
  real(r8), allocatable ::   airrrv_accum(:)
  real(r8), allocatable ::  cantemp_accum(:)
  real(r8), allocatable ::   canrrv_accum(:)
  real(r8), allocatable :: skintemp_accum(:)
  real(r8), allocatable ::   sfluxt_accum(:)
  real(r8), allocatable ::   sfluxr_accum(:) ! fast can nud
  real(r8), allocatable ::      pcp_accum(:) ! fast can nud
  real(r8), allocatable ::  sfctemp_accum(:) ! fast can nud
  real(r8), allocatable ::  fracliq_accum(:) ! fast can nud
  real(r8), allocatable ::   runoff_accum(:)

  real(r8), allocatable :: airtemp_dmin_accum (:)
  real(r8), allocatable :: airtemp_dmax_accum (:)
  real(r8), allocatable :: cantemp_dmin_accum (:)
  real(r8), allocatable :: cantemp_dmax_accum (:)

! Land arrays

  real(r8), allocatable :: wxferi_accum (:) ! filled in leaf4_soil
  real(r8), allocatable :: wxferp_accum (:) ! filled in leaf4_soil
  real(r8), allocatable :: wxfer1_accum (:) ! filled in leaf4_soil

  real(r8), allocatable ::  vegtemp_dmin_accum (:)
  real(r8), allocatable ::  vegtemp_dmax_accum (:)
  real(r8), allocatable :: soiltemp_dmin_accum (:)
  real(r8), allocatable :: soiltemp_dmax_accum (:)

Contains

!===============================================================================

subroutine alloc_flux_accum(mza,mva,mwa,mwsfc,mland)

  use misc_coms,   only: ilwrtyp, iswrtyp
  use leaf_coms,   only: isfcl
  use consts_coms, only: r8
  use oname_coms,  only: nl

  implicit none

  integer, intent(in) :: mza, mva, mwa, mwsfc, mland

  if (.not. nl%do_accum) return

! Allocate arrays for accumulated flux quantities
! Initialize arrays to zero

  if (iswrtyp > 0) then
     allocate (      rshort_accum(mwa)) ;       rshort_accum = 0._r8
     allocate (    rshortup_accum(mwa)) ;     rshortup_accum = 0._r8
     allocate (  rshort_top_accum(mwa)) ;   rshort_top_accum = 0._r8
     allocate (rshortup_top_accum(mwa)) ; rshortup_top_accum = 0._r8
  endif

  if (ilwrtyp > 0) then
     allocate (          rlong_accum(mwa)) ;           rlong_accum = 0._r8
     allocate (        rlongup_accum(mwa)) ;         rlongup_accum = 0._r8
     allocate (    rlongup_top_accum(mwa)) ;     rlongup_top_accum = 0._r8
  endif

  if (.false.) then
     allocate (   vc_accum(mza,mva)) ;    vc_accum = 0._r8
     allocate    (wc_accum(mza,mwa)) ;    wc_accum = 0._r8
     allocate (press_accum(mza,mwa)) ; press_accum = 0._r8
     allocate ( tair_accum(mza,mwa)) ;  tair_accum = 0._r8
     allocate ( rr_v_accum(mza,mwa)) ;  rr_v_accum = 0._r8
  endif

  if (.true.) then
     allocate (latheat_liq_accum(mza,mwa)) ; latheat_liq_accum = 0._r8
     allocate (latheat_ice_accum(mza,mwa)) ; latheat_ice_accum = 0._r8
  endif

  if (isfcl > 0) then
     allocate (    vels_accum(mwsfc)) ;     vels_accum = 0._r8
     allocate ( airtemp_accum(mwsfc)) ;  airtemp_accum = 0._r8
     allocate (  airrrv_accum(mwsfc)) ;   airrrv_accum = 0._r8
     allocate ( cantemp_accum(mwsfc)) ;  cantemp_accum = 0._r8
     allocate (  canrrv_accum(mwsfc)) ;   canrrv_accum = 0._r8
     allocate (skintemp_accum(mwsfc)) ; skintemp_accum = 0._r8
     allocate (  sfluxt_accum(mwsfc)) ;   sfluxt_accum = 0._r8
     allocate (  sfluxr_accum(mwsfc)) ;   sfluxr_accum = 0._r8 ! fast can nud
     allocate (     pcp_accum(mwsfc)) ;      pcp_accum = 0._r8 ! fast can nud
     allocate ( sfctemp_accum(mwsfc)) ;  sfctemp_accum = 0._r8 ! fast can nud
     allocate ( fracliq_accum(mwsfc)) ;  fracliq_accum = 0._r8 ! fast can nud
     allocate (  runoff_accum(mwsfc)) ;   runoff_accum = 0._r8

     allocate (wxferi_accum  (mland)) ; wxferi_accum   = 0._r8
     allocate (wxferp_accum  (mland)) ; wxferp_accum   = 0._r8
     allocate (wxfer1_accum  (mland)) ; wxfer1_accum   = 0._r8

     allocate ( airtemp_dmin_accum (mwsfc)) ;  airtemp_dmin_accum = 0._r8
     allocate ( airtemp_dmax_accum (mwsfc)) ;  airtemp_dmax_accum = 0._r8
     allocate ( cantemp_dmin_accum (mwsfc)) ;  cantemp_dmin_accum = 0._r8
     allocate ( cantemp_dmax_accum (mwsfc)) ;  cantemp_dmax_accum = 0._r8

     allocate ( vegtemp_dmin_accum (mland)) ;  vegtemp_dmin_accum = 0._r8
     allocate ( vegtemp_dmax_accum (mland)) ;  vegtemp_dmax_accum = 0._r8
     allocate (soiltemp_dmin_accum (mland)) ; soiltemp_dmin_accum = 0._r8
     allocate (soiltemp_dmax_accum (mland)) ; soiltemp_dmax_accum = 0._r8
  endif

end subroutine alloc_flux_accum

!===============================================================================

subroutine filltab_flux_accum()

  use var_tables, only: increment_vtable
  use oname_coms, only: nl

  implicit none

  if (.not. nl%do_accum) return

  if (allocated(      rshort_accum)) call increment_vtable(      'RSHORT_ACCUM','AW', dvar1=      rshort_accum)
  if (allocated(    rshortup_accum)) call increment_vtable(    'RSHORTUP_ACCUM','AW', dvar1=    rshortup_accum)
  if (allocated(       rlong_accum)) call increment_vtable(       'RLONG_ACCUM','AW', dvar1=       rlong_accum)
  if (allocated(     rlongup_accum)) call increment_vtable(     'RLONGUP_ACCUM','AW', dvar1=     rlongup_accum)
  if (allocated(  rshort_top_accum)) call increment_vtable(  'RSHORT_TOP_ACCUM','AW', dvar1=  rshort_top_accum)
  if (allocated(rshortup_top_accum)) call increment_vtable('RSHORTUP_TOP_ACCUM','AW', dvar1=rshortup_top_accum)
  if (allocated( rlongup_top_accum)) call increment_vtable( 'RLONGUP_TOP_ACCUM','AW', dvar1= rlongup_top_accum)

  if (allocated(          vc_accum)) call increment_vtable(          'VC_ACCUM','AV', dvar2=          vc_accum)
  if (allocated(          wc_accum)) call increment_vtable(          'WC_ACCUM','AW', dvar2=          wc_accum)
  if (allocated(       press_accum)) call increment_vtable(       'PRESS_ACCUM','AW', dvar2=       press_accum)
  if (allocated(        tair_accum)) call increment_vtable(        'TAIR_ACCUM','AW', dvar2=        tair_accum)
  if (allocated(        rr_v_accum)) call increment_vtable(        'RR_V_ACCUM','AW', dvar2=        rr_v_accum)
  if (allocated( latheat_liq_accum)) call increment_vtable( 'LATHEAT_LIQ_ACCUM','AW', dvar2= latheat_liq_accum)
  if (allocated( latheat_ice_accum)) call increment_vtable( 'LATHEAT_ICE_ACCUM','AW', dvar2= latheat_ice_accum)

  if (allocated(      vels_accum)) call increment_vtable(    'VELS_ACCUM','CW', dvar1=    vels_accum)
  if (allocated(   airtemp_accum)) call increment_vtable( 'AIRTEMP_ACCUM','CW', dvar1= airtemp_accum)
  if (allocated(    airrrv_accum)) call increment_vtable(  'AIRRRV_ACCUM','CW', dvar1=  airrrv_accum)
  if (allocated(   cantemp_accum)) call increment_vtable( 'CANTEMP_ACCUM','CW', dvar1= cantemp_accum)
  if (allocated(    canrrv_accum)) call increment_vtable(  'CANRRV_ACCUM','CW', dvar1=  canrrv_accum)
  if (allocated(  skintemp_accum)) call increment_vtable('SKINTEMP_ACCUM','CW', dvar1=skintemp_accum)
  if (allocated(    sfluxt_accum)) call increment_vtable(  'SFLUXT_ACCUM','CW', dvar1=  sfluxt_accum)
  if (allocated(    sfluxr_accum)) call increment_vtable(  'SFLUXR_ACCUM','CW', dvar1=  sfluxr_accum) ! fast can nud
  if (allocated(       pcp_accum)) call increment_vtable(     'PCP_ACCUM','CW', dvar1=     pcp_accum) ! fast can nud
  if (allocated(   sfctemp_accum)) call increment_vtable( 'SFCTEMP_ACCUM','CW', dvar1= sfctemp_accum) ! fast can nud
  if (allocated(   fracliq_accum)) call increment_vtable( 'FRACLIQ_ACCUM','CW', dvar1= fracliq_accum) ! fast can nud
  if (allocated(    runoff_accum)) call increment_vtable(  'RUNOFF_ACCUM','CW', dvar1=  runoff_accum)

  if (allocated(    wxferi_accum)) call increment_vtable(  'WXFERI_ACCUM','LW', dvar1=  wxferi_accum)
  if (allocated(    wxferp_accum)) call increment_vtable(  'WXFERP_ACCUM','LW', dvar1=  wxferp_accum)
  if (allocated(    wxfer1_accum)) call increment_vtable(  'WXFER1_ACCUM','LW', dvar1=  wxfer1_accum)

  if (allocated( airtemp_dmin_accum)) call increment_vtable( 'AIRTEMP_DMIN_ACCUM','CW', dvar1= airtemp_dmin_accum)
  if (allocated( airtemp_dmax_accum)) call increment_vtable( 'AIRTEMP_DMAX_ACCUM','CW', dvar1= airtemp_dmax_accum)
  if (allocated( cantemp_dmin_accum)) call increment_vtable( 'CANTEMP_DMIN_ACCUM','CW', dvar1= cantemp_dmin_accum)
  if (allocated( cantemp_dmax_accum)) call increment_vtable( 'CANTEMP_DMAX_ACCUM','CW', dvar1= cantemp_dmax_accum)

  if (allocated( vegtemp_dmin_accum)) call increment_vtable( 'VEGTEMP_DMIN_ACCUM','LW', dvar1= vegtemp_dmin_accum)
  if (allocated( vegtemp_dmax_accum)) call increment_vtable( 'VEGTEMP_DMAX_ACCUM','LW', dvar1= vegtemp_dmax_accum)
  if (allocated(soiltemp_dmin_accum)) call increment_vtable('SOILTEMP_DMIN_ACCUM','LW', dvar1=soiltemp_dmin_accum)
  if (allocated(soiltemp_dmax_accum)) call increment_vtable('SOILTEMP_DMAX_ACCUM','LW', dvar1=soiltemp_dmax_accum)

end subroutine filltab_flux_accum

!===============================================================================

subroutine flux_accum()

  use misc_coms,   only: time_istp8p, dtlm, dtsm, ilwrtyp, iswrtyp, iparallel

  use mem_ijtabs,  only: istp, jtab_v, jtab_w, jtv_prog, &
                         jtw_prog, mrl_begl, mrl_endl

  use mem_basic,   only: vc, wc, press, tair, rr_v

  use mem_radiate, only: albedt, rshort, rlong, rlongup, &
                         rshort_top, rshortup_top, rlongup_top

  use mem_average_vars, only:  airtempk_dmin,  airtempk_dmax, &
                               cantempk_dmin,  cantempk_dmax, &
                               vegtempk_dmin,  vegtempk_dmax, &
                              soiltempk_dmin, soiltempk_dmax

  use leaf_coms,   only: isfcl

  use consts_coms, only: r8

  use mem_sfcg,    only: sfcg, mwsfc, itab_wsfc
  use mem_land,    only: land, mland, omland, nzg, dslz
  use mem_lake,    only: lake, mlake, omlake
  use mem_sea,     only: sea, msea, omsea
  use mem_grid,    only: mza, lpv, lpw
  use therm_lib,   only: qtk, qwtk
  use mem_para,    only: myrank
  use oname_coms,  only: nl

  implicit none

  integer :: j, iv, iw, iwsfc, iland, ilake, isea, k, nls

  real :: soiltempk, tempk, fracliq, fldval

  real :: w_comb        ! (sfcwater + soil) water mass [kg/m^2]
  real :: qw_comb       ! (sfcwater + soil) energy [J/m^2]
  real :: hcapsoil      ! soil heat capacity [J/(m^2 K)]

  ! DO_ACCUM is a namelist option to turn on/off computing and writing the
  ! accumulation fields to the history file. Currently it defaults to .false.
  ! because it is causing crashes on large parallel or OpenMP runs. The issue
  ! is with the _dmin and _dmax variables from mem_average that are used here;
  ! I need to verify that they are being computed properly in parallel runs.
  ! They will also need to be added to the vtables and history files for
  ! everything to work...

  if (.not. nl%do_accum) return

! Update accumulations of ATM velocity and pressure

  if (allocated(vc_accum)) then

     !$omp parallel do private(iv,k)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)

        ! Timestep for accumulating velocity and pressure, DTSM, is the small
        ! (acoustic) timestep. The frequency that each IV cell is processed in this
        ! loop should be consistent with (inversely proportional to) DTSM

        do k = lpv(iv),mza
           vc_accum(k,iv) = vc_accum(k,iv) + dtsm * real(vc(k,iv),r8)
        enddo

     enddo
     !$omp end parallel do

  endif

  !$omp parallel do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Timestep for accumulating velocity and pressure, DTSM, is the small
     ! (acoustic) timestep. The frequency that each IV cell is processed in this
     ! loop should be consistent with (inversely proportional to) DTSM

     if (allocated(wc_accum)) then
        do k = lpw(iw),mza
           wc_accum(k,iw) = wc_accum(k,iw) + dtsm * real(wc(k,iw),r8)
        enddo
     endif

     if (allocated(press_accum)) then
        do k = lpw(iw),mza
           press_accum(k,iw) = press_accum(k,iw) + dtsm * press(k,iw)
        enddo
     endif

  enddo
  !$omp end parallel do

! Update accumulations of ATM temperature and specific humidity

  if (mrl_endl(istp) > 0) then

     !$omp parallel do private(iw,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Timestep for accumulating temperature and moisture, DTLM, is the long
        ! timestep. The frequency that each IW cell is processed in this
        ! loop should be consistent with (inversely proportional to) DTLM

        if (allocated(tair_accum)) then
           do k = lpw(iw),mza
              tair_accum(k,iw) = tair_accum(k,iw) + dtlm * real(tair(k,iw),r8)
           enddo
        endif

        if (allocated(rr_v_accum)) then
           do k = lpw(iw),mza
              rr_v_accum(k,iw) = rr_v_accum(k,iw) + dtlm * real(rr_v(k,iw),r8)
           enddo
        endif

     enddo
     !$omp end parallel do

  endif

! Update accumulations of ATM SFC and TOP radiative fluxes

  if (mrl_begl(istp) > 0 .and. ilwrtyp + iswrtyp > 0) then

     !$omp parallel do private(iw)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Timestep for accumulating fluxes is DTLM for the given IW cell;
        ! the frequency that each IW cell is processed in this loop should be
        ! consistent with (inversely proportional to) DTLM

        if (iswrtyp > 0) then
           rshort_accum      (iw) = rshort_accum      (iw) + dtlm * real(rshort(iw),r8)
           rshortup_accum    (iw) = rshortup_accum    (iw) + dtlm * real(rshort(iw) * albedt(iw),r8)
           rshort_top_accum  (iw) = rshort_top_accum  (iw) + dtlm * real(rshort_top(iw),r8)
           rshortup_top_accum(iw) = rshortup_top_accum(iw) + dtlm * real(rshortup_top(iw),r8)
        endif

        if (ilwrtyp > 0) then
           rlong_accum      (iw) = rlong_accum      (iw) + dtlm * real(rlong(iw),r8)
           rlongup_accum    (iw) = rlongup_accum    (iw) + dtlm * real(rlongup(iw),r8)
           rlongup_top_accum(iw) = rlongup_top_accum(iw) + dtlm * real(rlongup_top(iw),r8)
        endif

     enddo
     !$omp end parallel do

  endif

  if (mrl_endl(istp) > 0 .and. isfcl > 0) then

     ! Update accumulations of SFC grid cells

     !$omp parallel do
     do iwsfc = 2,mwsfc

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

           vels_accum(iwsfc) =    vels_accum(iwsfc) + dtlm * real(sfcg%vels   (iwsfc),r8)
        airtemp_accum(iwsfc) = airtemp_accum(iwsfc) + dtlm * real(sfcg%airtemp(iwsfc),r8)
         airrrv_accum(iwsfc) =  airrrv_accum(iwsfc) + dtlm * real(sfcg%airrrv (iwsfc),r8)

        cantemp_accum(iwsfc) = cantemp_accum(iwsfc) + dtlm * real(sfcg%cantemp(iwsfc),r8)
         canrrv_accum(iwsfc) =  canrrv_accum(iwsfc) + dtlm * real(sfcg%canrrv (iwsfc),r8)
         sfluxt_accum(iwsfc) =  sfluxt_accum(iwsfc) + dtlm * real(sfcg%sfluxt (iwsfc),r8)
         sfluxr_accum(iwsfc) =  sfluxr_accum(iwsfc) + dtlm * real(sfcg%sfluxr (iwsfc),r8)
            pcp_accum(iwsfc) =     pcp_accum(iwsfc) +        real(sfcg%pcpg   (iwsfc),r8) ! sfcg%pcpg   is amount in timestep, not flux
         runoff_accum(iwsfc) =  runoff_accum(iwsfc) +        real(sfcg%runoff (iwsfc),r8) ! sfcg%runoff is amount in timestep, not flux

        if (airtempk_dmin(iwsfc) > sfcg%airtemp(iwsfc)) airtempk_dmin(iwsfc) = sfcg%airtemp(iwsfc)
        if (airtempk_dmax(iwsfc) < sfcg%airtemp(iwsfc)) airtempk_dmax(iwsfc) = sfcg%airtemp(iwsfc)

        if (cantempk_dmin(iwsfc) > sfcg%cantemp(iwsfc)) cantempk_dmin(iwsfc) = sfcg%cantemp(iwsfc)
        if (cantempk_dmax(iwsfc) < sfcg%cantemp(iwsfc)) cantempk_dmax(iwsfc) = sfcg%cantemp(iwsfc)

     enddo
     !$omp end parallel do

     ! Update accumulations of LAND cells

     !$omp parallel do private(iwsfc, nls, soiltempk, tempk, fracliq, fldval, &
     !$omp                     w_comb, qw_comb, hcapsoil)
     do iland = 2,mland
        iwsfc = iland + omland

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        nls = land%nlev_sfcwater(iland)

        call qwtk(land%soil_energy(nzg,iland),        &
                  land%soil_water(nzg,iland)*1.e3,    &
                  land%specifheat_drysoil(nzg,iland), &
                  soiltempk, fracliq)

        if (nls > 0) then
           call qtk(land%sfcwater_energy(nls,iland),tempk,fracliq)

           fldval = (1. - land%vf(iland)) * tempk &
                        + land%vf(iland)  * land%veg_temp(iland)
        else
           fldval = (1. - land%vf(iland)) * soiltempk &
                        + land%vf(iland)  * land%veg_temp(iland)
        endif

        skintemp_accum(iwsfc) = skintemp_accum(iwsfc) + dtlm * real(fldval,r8)

        if (vegtempk_dmin(iland) > land%veg_temp(iland)) vegtempk_dmin(iland) = land%veg_temp(iland)
        if (vegtempk_dmax(iland) < land%veg_temp(iland)) vegtempk_dmax(iland) = land%veg_temp(iland)

        if (soiltempk_dmin(iland) > soiltempk) soiltempk_dmin(iland) = soiltempk
        if (soiltempk_dmax(iland) < soiltempk) soiltempk_dmax(iland) = soiltempk

        ! sfctemp_accum and fracliq_accum variables over land
        ! (used in FAST_CANOPY) are computed with sfcwater(1) and soil(nzg)
        ! in thermal and phase equilibrium with each other

        ! Combined sfcwater and soil water mass and energy per square meter

        w_comb  = land%sfcwater_mass(1,iland) + land%soil_water(nzg,iland) * 1.e3 * dslz(nzg)
        qw_comb = land%sfcwater_energy(1,iland) * land%sfcwater_mass(1,iland) &
                + land%soil_energy(nzg,iland) * dslz(nzg)

        ! Soil heat capacity per square meter

        hcapsoil = land%specifheat_drysoil(nzg,iland) * dslz(nzg)

        ! Diagnose equilibrium temperature and fractional liquid/ice water phases

        call qwtk(qw_comb,w_comb,hcapsoil,tempk,fracliq)

        sfctemp_accum(iwsfc) = sfctemp_accum(iwsfc) + dtlm * real(tempk,r8)
        fracliq_accum(iwsfc) = fracliq_accum(iwsfc) + dtlm * real(fracliq,r8)

     enddo
     !$omp end parallel do

     ! Update accumulations of LAKE cells

     !$omp parallel do private (iwsfc, tempk, fracliq)
     do ilake = 2,mlake
        iwsfc = ilake + omlake

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        call qtk(lake%lake_energy(ilake),tempk,fracliq)
        skintemp_accum(iwsfc) = skintemp_accum(iwsfc) + dtlm * real(tempk,r8)
        sfctemp_accum (iwsfc) = sfctemp_accum (iwsfc) + dtlm * real(tempk,r8)
        fracliq_accum (iwsfc) = fracliq_accum (iwsfc) + dtlm * real(fracliq,r8)
     enddo
     !$omp end parallel do

     ! Update accumulations of SEA cells

     !$omp parallel do private(iwsfc,nls,fldval)
     do isea = 2,msea
        iwsfc = isea + omsea

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        nls = sea%nlev_seaice(isea)

        if (nls > 0) then
           fldval = (1.0 - sea%seaicec(isea)) * sea%seatc(isea) &
                  +        sea%seaicec(isea)  * sea%seaice_tempk(nls,isea)
        else
           fldval = sea%seatc(isea)
        endif

        skintemp_accum(iwsfc) = skintemp_accum(iwsfc) + dtlm * real(fldval,r8)
        sfctemp_accum (iwsfc) = sfctemp_accum (iwsfc) + dtlm * real(fldval,r8)
        fracliq_accum (iwsfc) = fracliq_accum (iwsfc) + dtlm * real(1.0 - sea%seaicec(isea),r8)

     enddo
     !$omp end parallel do

  endif

  ! Apply daily min/max accumulations

  if (mod(time_istp8p,86400.0_r8) < dtsm .and. isfcl > 0) then

     ! SFC grid cells

     !$omp parallel do
     do iwsfc = 2,mwsfc

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        airtemp_dmin_accum(iwsfc) = airtemp_dmin_accum(iwsfc) + real(airtempk_dmin(iwsfc),r8) * 86400._r8
        airtemp_dmax_accum(iwsfc) = airtemp_dmax_accum(iwsfc) + real(airtempk_dmax(iwsfc),r8) * 86400._r8

        cantemp_dmin_accum(iwsfc) = cantemp_dmin_accum(iwsfc) + real(cantempk_dmin(iwsfc),r8) * 86400._r8
        cantemp_dmax_accum(iwsfc) = cantemp_dmax_accum(iwsfc) + real(cantempk_dmax(iwsfc),r8) * 86400._r8

     enddo
     !$omp end parallel do

     ! LAND cells

     !$omp parallel do
     do iland = 2,mland
        iwsfc = iland + omland

        ! Skip this cell if running in parallel and cell rank is not MYRANK
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

        vegtemp_dmin_accum(iland) = vegtemp_dmin_accum(iland) + real(vegtempk_dmin(iland),r8) * 86400._r8
        vegtemp_dmax_accum(iland) = vegtemp_dmax_accum(iland) + real(vegtempk_dmax(iland),r8) * 86400._r8

        soiltemp_dmin_accum(iland) = soiltemp_dmin_accum(iland) + real(soiltempk_dmin(iland),r8) * 86400._r8
        soiltemp_dmax_accum(iland) = soiltemp_dmax_accum(iland) + real(soiltempk_dmax(iland),r8) * 86400._r8

     enddo
     !$omp end parallel do

  endif

end subroutine flux_accum

!===============================================================================

End Module mem_flux_accum
