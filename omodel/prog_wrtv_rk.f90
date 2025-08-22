module wrtv_rk

  use consts_coms, only: r8

  implicit none

  real, allocatable :: b20(:)

  real(r8), parameter :: onethird = 1.0_r8 / 3.0_r8

  real(r8), parameter :: xstg3(3) = [ onethird, 0.5_r8, 1.0_r8 ]
  real(r8), parameter :: xstg2(2) =           [ 0.5_r8, 1.0_r8 ]

  ! MPI communication tags
  integer, parameter :: itag_scpt  = 20
  integer, parameter :: itag_gxyps = 21
  integer, parameter :: itag_monot = 22

  ! Monotonic / positive-definite tolerances
  real, parameter :: eps0 = 1.e-32
  real, parameter :: onep = 1.000001
  real, parameter :: onem = 0.999999

  ! For testing - may want to add to namelist!
  integer, parameter  :: iorderv        =  3
  logical, parameter  :: centered_monot = .false.

  private
  public :: init_wrtv_rk, prog_wrtv_rk

contains

!===============================================================================

subroutine init_wrtv_rk()

  use mem_grid,  only: mza
  use misc_coms, only: dtsm
  use mem_rayf,  only: dorayfw, rayf_cofw, krayfw_bot

  implicit none

  integer :: k

  allocate(b20(mza)) ; b20 = 1.0

  if (dorayfw) then
     do k = krayfw_bot, mza-1
        b20(k) = 1.0 + real(dtsm) * rayf_cofw(k)
     enddo
  endif

end subroutine init_wrtv_rk

!===============================================================================

subroutine prog_wrtv_rk()

! This dynamic core version for hexagonal cells combines the Perot (2002)
! finite-volume method for evaluating momentum advective and diffusive flux
! convergences in T cells, the Miura (2007) piecewise-linear advection algorithm,
! and the Walko and Avissar (2008a,b) time differencing method.

! The above description needs amending: The time differencing method is RK3.

  use mem_ijtabs,   only: jtab_v, jtab_w, jtw_wadj, jtm_vadj, jtv_prog, istp, &
                          jtv_wadj, jtv_lbcp, jtw_prog, jtw_lbcp, itab_w, itab_v
  use mem_basic,    only: rho, thil, wc, press, vmc, vc, wmc, vxe, vye, vze, &
                          vmasc
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, volt, xev, nve2_max, lve2, arv
  use mem_tend,     only: thilt, vmxet, vmyet, vmzet, vmt
  use misc_coms,    only: iparallel, dtsm, dtlm, nrk_wrtv, dn01d, th01d, &
                          initial, deltax, nxp, mdomain
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
  use obnd,         only: lbcopy_m, lbcopy_w, lbcopy_v
  use oname_coms,   only: nl
  use vel_t3d,      only: vel_t3d_hex
  use consts_coms,  only: p00i, rocp
  use mem_turb,     only: akmodx, akhodx, khtopv, kmtopv
  use pbl_drivers,  only: solve_eddy_diff_heat, solve_eddy_diff_vxe
  use mem_rayf,     only: dorayf, krayf_bot, rayf_cof, dorayfmix, rayf_mix_top_vxe, &
                          vc03d
  use mem_nudge,    only: vmanud, nudflag

  implicit none

  integer  :: j, iv, iw, k, ksw, iw1, iw2, iwn, jv, istage
  integer  :: i2dv, i2dt
  real     :: dts, v4, rs
  real(r8) :: dt8

! automatic arrays

  real :: thil_upv(mza,mva) ! Upstreamed THIL at each V interface
  real :: vxe_upv (mza,mva) ! Upstreamed VXE  at each V interface
  real :: vye_upv (mza,mva) ! Upstreamed VYE  at each V interface
  real :: vze_upv (mza,mva) ! Upstreamed VZE  at each V interface

  real :: thil_upw(mza,mwa) ! Upstreamed THIL at each W level
  real :: vxe_upw (mza,mwa) ! Upstreamed VXE  at each W level
  real :: vye_upw (mza,mwa) ! Upstreamed VYE  at each W level
  real :: vze_upw (mza,mwa) ! Upstreamed VZE  at each W level

  real :: vmxet_rk(mza,mwa)
  real :: vmyet_rk(mza,mwa)
  real :: vmzet_rk(mza,mwa)
  real :: thilt_rk(mza,mwa)

  real(r8) :: rho0(mza,mwa)
  real(r8) :: rth0(mza,mwa)
  real     :: wmc0(mza,mwa)
  real     :: vmc0(mza,mva)
  real     :: thl0(mza,mva)

! real :: gzps_th (mza), gzps_vx (mza), gzps_vy (mza), gzps_vz (mza)
! real :: gzzps_th(mza), gzzps_vx(mza), gzzps_vy(mza), gzzps_vz(mza)

  real :: thilt_short(mza,mwa)
  real :: vmxet_short(mza,mwa)
  real :: vmyet_short(mza,mwa)
  real :: vmzet_short(mza,mwa)
  real :: vmt_short  (mza,mva)

  ! For prognostic shaved-cell method
  real :: vmxe0(nve2_max,mwa)
  real :: vmye0(nve2_max,mwa)
  real :: vmze0(nve2_max,mwa)

  real :: unit_dist, fracx, rayfx
  real :: vmca(mza,mva)

  i2dt = 0
  if (nl%thil_horiz_adv_order == 2) i2dt = 2
  if (nl%thil_horiz_adv_order == 3) i2dt = 5

  i2dv = 0
  if (nl%wind_horiz_adv_order == 2) i2dv = 2
  if (nl%wind_horiz_adv_order == 3) i2dv = 5

  ! Divergence/vorticity damping if computed each long timestep

  if ( (.not. nl%divh_damp_short) .and. (istp == 1) ) then
     call divh_damp(vmt, dtlm)
  endif

  if ( (.not. nl%vort_damp_short) .and. (istp == 1) ) then
     call vort_damp(vmt, dtlm)
  endif

  !$omp parallel
  !$omp do private(iw,k,v4,ksw,rs,jv,iv,iwn)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     thil_upw(      mza,iw) = 0.
     thil_upw(lpw(iw)-1,iw) = 0.

     do k = lpw(iw), mza
        v4 = real(volt(k,iw))

        ! Store starting values for R-K timestepping
        rth0(k,iw) = thil(k,iw) * rho(k,iw)
        rho0(k,iw) = rho (k,iw)
        wmc0(k,iw) = wmc (k,iw)
        thl0(k,iw) = thil(k,iw)

        ! Include long-timestep tendencies in each short timestep
        ! and weight by volume
        thilt_short(k,iw) = thilt(k,iw) * v4
        vmxet_short(k,iw) = vmxet(k,iw) * v4
        vmyet_short(k,iw) = vmyet(k,iw) * v4
        vmzet_short(k,iw) = vmzet(k,iw) * v4
     enddo

     ! Save initial earth-cartesian momentum for prognostic shaved-cell method

     do ksw = 1, lve2(iw)
        k = ksw + lpw(iw) - 1
        rs = real(rho(k,iw))

        vmxe0(ksw,iw) = vxe(k,iw) * rs
        vmye0(ksw,iw) = vye(k,iw) * rs
        vmze0(ksw,iw) = vze(k,iw) * rs
     enddo

     ! Vertical diffusion now performed every short timestep for heat/momentum

     call solve_eddy_diff_heat(iw, thilt_short(:,iw))
     call solve_eddy_diff_vxe (iw, vmxet_short(:,iw), vmyet_short(:,iw), &
                                   vmzet_short(:,iw))

     ! Rayleigh friction on velocity gradient at model top

     if (dorayfmix) then
        call rayf_mix_top_vxe (iw, vmxet_short(:,iw), vmyet_short(:,iw), &
                                   vmzet_short(:,iw))
     endif

     ! Rayleigh friction on thil - only for horizontally homogeneous initialization

     if (dorayf .and. initial == 1) then
        do k = krayf_bot, mza
           v4 = real(volt(k,iw))
           thilt_short(k,iw) = thilt_short(k,iw) &
                + rayf_cof(k) * dn01d(k) * (th01d(k) - thil(k,iw)) * v4
        enddo
     endif

     ! Horizontal diffusion

     do jv = 1, itab_w(iw)%npoly
        iv  = itab_w(iw)%iv(jv)
        iwn = itab_w(iw)%iw(jv)

        ! Loop over T levels
        do k = lpv(iv), khtopv(iv)
           thilt_short(k,iw) = thilt_short(k,iw) + akhodx(k,iv) * (thil(k,iwn) - thil(k,iw))
        enddo
        do k = lpv(iv), kmtopv(iv)
           vmxet_short(k,iw) = vmxet_short(k,iw) + akmodx(k,iv) * (vxe(k,iwn) - vxe(k,iw))
           vmyet_short(k,iw) = vmyet_short(k,iw) + akmodx(k,iv) * (vye(k,iwn) - vye(k,iw))
           vmzet_short(k,iw) = vmzet_short(k,iw) + akmodx(k,iv) * (vze(k,iwn) - vze(k,iw))
        enddo
     enddo

  enddo
  !$omp end do nowait

  !$omp do private(iv)
  do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)

     vmc0     (:,iv) = vmc(:,iv)
     vmt_short(:,iv) = vmt(:,iv)

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! Rayleigh friction on vmc

  if (dorayf) then

     unit_dist = sqrt(sqrt(4./3.)) * deltax  ! approx = 1.07457 * deltax; for mdomain = 4 only

     !$omp parallel do private(iv,iw1,iw2,fracx,rayfx,k)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        if (mdomain /= 4) then

           ! Vertical loop over V points
           do k = krayf_bot, mza
              vmt_short(k,iv) = vmt_short(k,iv) + rayf_cof(k) &
                   * 0.5 * (rho(k,iw1) + rho(k,iw2)) * (vc03d(k,iv) - vc(k,iv))
           enddo

        else

           ! Coefficient for extra RAYF damping at ends of channel with cyclic BC's

           fracx = abs(xev(iv)) / (0.5 * real(nxp) * unit_dist)
           rayfx = (-2. + 3. * fracx) * (0.2 * rayf_cof(mza))  ! max of (0.2 * rayf_cof) at top bnd

           ! Vertical loop over V points
           do k = lpv(iv), mza
              vmt_short(k,iv) = vmt_short(k,iv) + max(rayf_cof(k),rayfx) &
                   * 0.5 * (rho(k,iw1) + rho(k,iw2)) * (vc03d(k,iv) - vc(k,iv))
           enddo

        endif

     enddo
     !$omp end parallel do

  endif ! (dorayf)

  ! Divergence/vorticity damping if computed each short timestep

  if (nl%divh_damp_short) call divh_damp(vmt_short, dtsm)
  if (nl%vort_damp_short) call vort_damp(vmt_short, dtsm)

  do istage = 1, nrk_wrtv

     if (nrk_wrtv <= 2) then
        dt8 = xstg2(istage) * dtsm
     else
        dt8 = xstg3(istage) * dtsm
     endif

     dts = real(dt8)

     !$omp parallel
     !$omp do private(iw,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        do k = lpw(iw), mza
           vmxet_rk(k,iw) = vmxet_short(k,iw)
           vmyet_rk(k,iw) = vmyet_short(k,iw)
           vmzet_rk(k,iw) = vmzet_short(k,iw)
           thilt_rk(k,iw) = thilt_short(k,iw)
        enddo

     enddo
     !$omp end do nowait

     !$omp do private(iv,k)
     do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

        if (nudflag > 0 .and. nl%nud_preserve_total_mass) then
           do k = lpv(iv), mza
              vmca(k,iv) = vmc(k,iv) * arv(k,iv) + vmanud(k,iv)
           enddo
        else
           do k = lpv(iv), mza
              vmca(k,iv) = vmc(k,iv) * arv(k,iv)
           enddo
        endif

        ! Save half-forward velocities from the final R-K stage for
        ! scalar transport

        if (istage == nrk_wrtv) then
           do k = lpv(iv), mza
              vmasc(k,iv) = vmasc(k,iv) + vmca(k,iv)
           enddo
        endif

     enddo
     !$omp end do nowait
     !$omp end parallel

     if (istage == nrk_wrtv .and. nl%ithil_monot > 0) then

        call fluxes_thil_monot( thil, thl0, thilt_short, thil_upw, thil_upv, &
                                rho, rho0, vmca, wmc, dts, nl%thil_horiz_adv_order, i2dt )
     else

        call fluxes_wrtv( thil, thil_upw, thil_upv, vmc, wmc, &
                          nl%thil_horiz_adv_order, i2dt )
     endif

     call fluxes_wrtv( vxe,  vxe_upw, vxe_upv, vmc, wmc, &
                       nl%wind_horiz_adv_order, i2dv )

     call fluxes_wrtv( vye,  vye_upw, vye_upv, vmc, wmc, &
                       nl%wind_horiz_adv_order, i2dv )

     call fluxes_wrtv( vze,  vze_upw, vze_upv, vmc, wmc, &
                       nl%wind_horiz_adv_order, i2dv )


     ! MAIN LOOP OVER W COLUMNS FOR UPDATING WM, WC, RHO, THIL, AND PRESS

     !$omp do private(iw,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Prognose vertical velocity, density, thil, and diagnose pressure

        call prog_wrt_begs( iw, istage, dt8,                      &
                            thil_upv, vxe_upv, vye_upv, vze_upv,  &
                            thil_upw, vxe_upw, vye_upw, vze_upw,  &
                            thilt_rk, vmxet_rk, vmyet_rk, vmzet_rk, &
                            rho0, rth0, wmc0, vmca )

     enddo
     !$omp end do nowait
     !$omp end parallel

     ! MPI SEND/RECV and LBC copy of quantities needed for prog_v:
     ! PRESS, RHO, VMXET_SHORT, VMYET_SHORT, and VMZET_SHORT

     if (iparallel == 1) then
        call mpi_send_w(dvara1=press, dvara2=rho, &
                        rvara1=vmxet_rk, rvara2=vmyet_rk, rvara3=vmzet_rk)
        call mpi_recv_w(dvara1=press, dvara2=rho, &
                        rvara1=vmxet_rk, rvara2=vmyet_rk, rvara3=vmzet_rk)
     endif

     call lbcopy_w(a1=vmxet_rk,  a2=vmyet_rk,  a3=vmzet_rk, &
                   d1=rho,       d2=press)

     ! MAIN LOOP OVER V POINTS TO UPDATE VMC AND VC

     !$omp parallel do private(iv)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)

        call prog_v_begs(iv, dts, vmc0, vmxet_rk, vmyet_rk, vmzet_rk, vmt_short)

     enddo
     !$omp end parallel do

     ! MPI SEND/RECV and LBC of VC and VMC

     if (iparallel == 1) then
        call mpi_send_v(rvara1=vmc, rvara2=vc)
        call mpi_recv_v(rvara1=vmc, rvara2=vc)
     endif
     call lbcopy_v(vmc=vmc, vc=vc)

     ! UPDATE EARTH CARTESIAN VELOCITIES

     !$omp parallel do private(iw)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (lve2(iw) > 0) then
           call update_vxe_undrgnd( iw, dts, vmxe0, vmye0, vmze0, &
                                    vmxet_rk, vmyet_rk, vmzet_rk )
        endif

        call vel_t3d_hex(iw)

     enddo
     !$omp end parallel do

     ! Parallel send-recieve of Earth Cartesian velocities, W, and THIL

     if (iparallel == 1) then
        call mpi_send_w(rvara1=vxe, rvara2=vye, rvara3=vze, &
                        rvara4=wmc, rvara5=wc, rvara6=thil)
        call mpi_recv_w(rvara1=vxe, rvara2=vye, rvara3=vze, &
                        rvara4=wmc, rvara5=wc, rvara6=thil)
     endif
     call lbcopy_w(a1=vxe, a2=vye, a3=vze, a4=wmc, a5=wc, a6=thil)

  enddo ! istage loop

! SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
! (Example of how to plot "external" field; one not available in module memory)
!
! if (mod(time8,op%frqplt) < dtlm .and. istp == 900) then
!    allocate (op%extfld(mza,mwa))
!    op%extfld(:,:) = vxe(:,:)
!    op%extfldname = 'VXE'
!    call plot_fields(11)
!    deallocate (op%extfld)
! endif
! END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

end subroutine prog_wrtv_rk


!=========================================================================


subroutine prog_wrt_begs( iw, istage, dt8,                      &
                          thil_upv, vxe_upv, vye_upv, vze_upv,  &
                          thil_upw, vxe_upw, vye_upw, vze_upw,  &
                          thilt_rk, vmxet_rk, vmyet_rk, vmzet_rk, &
                          rho0, rth0, wmc0, vmca )

  use mem_ijtabs,  only: itab_w
  use mem_basic,   only: wmc, rho, thil, wc, press, vxe, vye, wmasc, &
                         alpha_press, pwfac
  use misc_coms,   only: nrk_wrtv, mdomain
  use consts_coms, only: cpocv, rocv, fcoriol, pi1, pio180
  use mem_grid,    only: mza, mva, mwa, lpv, lpw, arw, wnx, wny, wnz, volt, &
                         volti, volwi, glatw, glonw, &
                         zwgt_top8, zwgt_bot8, gdz_wgtp8, gdz_wgtm8, &
                         gdz_wgtp, gdz_wgtm
  use tridiag,     only: tridiffo
  use oname_coms,  only: nl
  use mem_nudge,   only: rhot_nud, nudflag, wmanud

  implicit none

  integer,  intent(in) :: iw, istage
  real(r8), intent(in) :: dt8

  real, intent(in)    :: thil_upv(mza,mva)
  real, intent(in)    :: vxe_upv (mza,mva)
  real, intent(in)    :: vye_upv (mza,mva)
  real, intent(in)    :: vze_upv (mza,mva)

  real, intent(in)    :: thil_upw(mza,mwa)
  real, intent(in)    :: vxe_upw (mza,mwa)
  real, intent(in)    :: vye_upw (mza,mwa)
  real, intent(in)    :: vze_upw (mza,mwa)

  real, intent(inout) :: thilt_rk(mza,mwa)
  real, intent(inout) :: vmxet_rk(mza,mwa)
  real, intent(inout) :: vmyet_rk(mza,mwa)
  real, intent(inout) :: vmzet_rk(mza,mwa)

  real(r8), intent(in) :: rth0(mza,mwa)
  real(r8), intent(in) :: rho0(mza,mwa)
  real,     intent(in) :: wmc0(mza,mwa)
  real,     intent(in) :: vmca(mza,mva)

  integer :: jv, iv
  integer :: k, ka, kbv

  real :: c8, c9, dts
  real :: rad0_swtc, rad_swtc, topo_swtc

  ! Vertical implicit scheme weighting parameters

  real(r8), parameter :: fr  = .55 ! rho
  real(r8), parameter :: fp2 = .55 ! press in W
  real(r8), parameter :: fp3 = .70 ! press in V

  real(r8), parameter :: pc2 = fp2 * cpocv
  real(r8), parameter :: pc3 = fp3 * cpocv

  ! Automatic arrays

! real(r8) :: wmarw    (mza)
!  real     :: hflx_rho (mza)
  real(r8) :: hflx_rho (mza)
  real(r8) :: delex_rho(mza)
  real(r8) :: delex_rhothil(mza)
  real(r8) :: po_swtc(mza)

  real(r8) :: rhothil   (mza)
! real(r8) :: thilt_rk  (mza)
  real(r8) :: vflux_thil(mza)

  real(r8) :: press_ex(mza)
  real(r8) :: rho_ex  (mza)

  real :: delex_wm(mza)
  real :: wmf     (mza)
  real :: mass, wma(mza)

  real :: vflux_vxe(mza)
  real :: vflux_vye(mza)
  real :: vflux_vze(mza)

  real(r8) :: b5(mza), b10(mza)

  real :: b7, b8, b21, b22, b23, b24, b25, b26
  real :: b15(mza), b31(mza), b32(mza), b33(mza)

  dts = real(dt8)

  ka = lpw(iw)

  ! Set bottom & top vertical advective mass, momentum, and heat fluxes to zero

! wma(ka-1) = 0.
! wma(mza)  = 0.

  vflux_thil(ka-1) = 0._r8
  vflux_thil(mza)  = 0._r8

  vflux_vxe(ka-1) = 0.
  vflux_vxe(mza)  = 0.

  vflux_vye(ka-1) = 0.
  vflux_vye(mza)  = 0.

  vflux_vze(ka-1) = 0.
  vflux_vze(mza)  = 0.

  ! Nudging tendency terms

  if (nudflag == 1 .and. nl%nud_preserve_total_mass) then

     do k = ka, mza-1
!       wmarw     (k) = wmanud(k,iw)
        vflux_thil(k) = wmanud(k,iw) * thil_upw(k,iw)
        vflux_vxe (k) = wmanud(k,iw) * vxe_upw(k,iw)
        vflux_vye (k) = wmanud(k,iw) * vye_upw(k,iw)
        vflux_vze (k) = wmanud(k,iw) * vze_upw(k,iw)
     enddo

     do k = lpw(iw), mza
        hflx_rho(k)    =                     wmanud (k-1,iw) - wmanud (k,iw)
        thilt_rk(k,iw) = thilt_rk(k,iw) + vflux_thil(k-1) - vflux_thil(k)
        vmxet_rk(k,iw) = vmxet_rk(k,iw) + vflux_vxe (k-1) - vflux_vxe (k)
        vmyet_rk(k,iw) = vmyet_rk(k,iw) + vflux_vye (k-1) - vflux_vye (k)
        vmzet_rk(k,iw) = vmzet_rk(k,iw) + vflux_vze (k-1) - vflux_vze (k)
     enddo

  else

     do k = ka, mza
!       thilt_rk(k) = thilt_short(k,iw)
!       hflx_rho(k) = 0.
        hflx_rho(k) = 0._r8
     enddo

  endif

  ! Sum advective tendencies over V neighbors of this cell

  do jv = 1, itab_w(iw)%npoly
     iv  = itab_w(iw)%iv(jv)
     kbv = lpv(iv)

     if ( itab_w(iw)%dirv(jv) > 0) then

        do k = kbv, mza
           hflx_rho(k)    = hflx_rho(k)    + vmca(k,iv)
           thilt_rk(k,iw) = thilt_rk(k,iw) + vmca(k,iv) * thil_upv(k,iv)
           vmxet_rk(k,iw) = vmxet_rk(k,iw) + vmca(k,iv) * vxe_upv (k,iv)
           vmyet_rk(k,iw) = vmyet_rk(k,iw) + vmca(k,iv) * vye_upv (k,iv)
           vmzet_rk(k,iw) = vmzet_rk(k,iw) + vmca(k,iv) * vze_upv (k,iv)
        enddo

     else

        do k = kbv, mza
           hflx_rho(k)    = hflx_rho(k)    - vmca(k,iv)
           thilt_rk(k,iw) = thilt_rk(k,iw) - vmca(k,iv) * thil_upv(k,iv)
           vmxet_rk(k,iw) = vmxet_rk(k,iw) - vmca(k,iv) * vxe_upv (k,iv)
           vmyet_rk(k,iw) = vmyet_rk(k,iw) - vmca(k,iv) * vye_upv (k,iv)
           vmzet_rk(k,iw) = vmzet_rk(k,iw) - vmca(k,iv) * vze_upv (k,iv)
        enddo

     endif

  enddo

  if (mdomain > 1) then

     ! Coriolis and large-scale PGF tendencies for limited-area run

     do k = ka, mza
        mass = real( rho(k,iw) * volt(k,iw) )
        vmxet_rk(k,iw) = vmxet_rk(k,iw) + mass * fcoriol * (vye(k,iw) - nl%v_geostrophic)
        vmyet_rk(k,iw) = vmyet_rk(k,iw) - mass * fcoriol * (vxe(k,iw) - nl%u_geostrophic)
     enddo

  else

     ! Coriolis tendencies for global run

     do k = ka, mza
        mass = real( rho(k,iw) * volt(k,iw) )
        vmxet_rk(k,iw) = vmxet_rk(k,iw) + mass * fcoriol * vye(k,iw)
        vmyet_rk(k,iw) = vmyet_rk(k,iw) - mass * fcoriol * vxe(k,iw)
     enddo

  endif

  ! Explicit density and theta tendency terms

  do k = ka, mza
     rhothil(k) = rho(k,iw) * thil(k,iw)

     b10(k) = dt8 * volti(k,iw)
     b5 (k) = alpha_press(k,iw) * real(rhothil(k)) ** rocv
     b15(k) = b10(k) * b5(k)

     ! Explicit density tendency
     delex_rho(k) = b10(k) * hflx_rho(k)

     ! Explicit density-thil tendency
     delex_rhothil(k) = b10(k) * thilt_rk(k,iw)
  enddo

  ! Include nudging terms in density tendency

  if (nudflag > 0 .and. .not. nl%nud_preserve_total_mass) then
     do k = ka, mza
        delex_rho(k) = delex_rho(k) + dt8 * rhot_nud(k,iw)
     enddo
  endif

  ! R-K time offset

  if (istage > 1) then
     do k = ka, mza
        delex_rho    (k) = delex_rho    (k) + rho0(k,iw) - rho(k,iw)
        delex_rhothil(k) = delex_rhothil(k) + rth0(k,iw) - rhothil(k)
     enddo
  endif

  ! Updated pressure and density from explicit tendency terms

  do k = ka, mza
     press_ex(k) = b5(k) * (rhothil(k) + pc2 * delex_rhothil(k))
     rho_ex  (k) =           rho(k,iw) + fr  * delex_rho    (k)
  enddo

  ! Explicit W momentum tendency

  do k = ka, mza-1

     delex_wm(k) = wmc0(k,iw) + dts * ( &

          ! Pressure gradient at W level

          + real( press_ex(k) - press_ex(k+1) ) * pwfac(k,iw) &

          ! Z-weighted buoyancy forcing in vertical between T levels

          - real( gdz_wgtm8(k) * rho_ex(k) + gdz_wgtp8(k) * rho_ex(k+1) ) &

          ! Volume-weighted momentum forcings in vertical between T levels

          + ( wnx(iw) * (vmxet_rk(k,iw) + vmxet_rk(k+1,iw)) &
            + wny(iw) * (vmyet_rk(k,iw) + vmyet_rk(k+1,iw)) &
            + wnz(iw) * (vmzet_rk(k,iw) + vmzet_rk(k+1,iw)) ) * volwi(k,iw) )

  enddo

  ! Fill matrix coefficients for implicit update of WM

  c8  =        dts * real(pc2)
  c9  = -dts * dts * real(fr)

  ! Loop over W pts
  do k = ka, mza-1
     b7  = dts * volwi(k,iw)
     b8  = c8  * pwfac(k,iw)

     b21 = b7 * wc(k-1,iw)
     b22 = b7 * wc(k+1,iw)

     b23 = b8 * b15(k)
     b24 = b8 * b15(k+1)

     b25 = c9 * volti(k  ,iw) * gdz_wgtm(k)
     b26 = c9 * volti(k+1,iw) * gdz_wgtp(k)

     b31(k) =        - arw(k-1,iw) * (b21 + b25 +  b23        * thil_upw(k-1,iw))
     b32(k) = b20(k) + arw(k  ,iw) * (b25 - b26 + (b23 + b24) * thil_upw(k  ,iw))
     b33(k) =          arw(k+1,iw) * (b22 + b26 -        b24  * thil_upw(k+1,iw))
  enddo

  ! Solve implicit tri-diagonal matrix equation for delta WM (del_wm)

  if (ka <= mza-1) then
     call tridiffo(mza,ka,mza-1,b31,b32,b33,delex_wm,wmf)
  endif

  ! Fluxes from updated vertical momentum

  wma(ka-1) = 0.
  wma(mza)  = 0.

  if (nudflag > 0 .and. nl%nud_preserve_total_mass) then
     do k = ka, mza-1
        wma(k) = wmf(k) * arw(k,iw) + wmanud(k,iw)
     enddo
  else
     do k = ka, mza-1
        wma(k) = wmf(k) * arw(k,iw)
     enddo
  endif

  do k = ka, mza-1
     vflux_thil(k) = wma(k) * thil_upw(k,iw)
     vflux_vxe (k) = wma(k) * vxe_upw (k,iw)
     vflux_vye (k) = wma(k) * vye_upw (k,iw)
     vflux_vze (k) = wma(k) * vze_upw (k,iw)
  enddo

  ! Add vertical momentum at (t+fw) to array for long-timestep scalar transport

  if (istage == nrk_wrtv) then
     do k = ka, mza-1
        wmasc(k,iw) = wmasc(k,iw) + wma(k)
     enddo
  endif

  ! For shallow water test cases 2 & 5, rho & press are
  ! interpreted as water depth & height

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     do k = ka, mza
        po_swtc(k) = rho(k,iw) + fp3 * delex_rho(k)
     enddo
  endif

  ! Vertical loop over T points

  do k = ka, mza

     ! Change of rho from (t) to (t+1)

!    rho(k,iw) = rho(k,iw) + delex_rho(k) + b10(k) * (wma(k-1) - wma(k))
     rho(k,iw) = rho0(k,iw) + dt8 * volti(k,iw) * (hflx_rho(k) + wma(k-1) - wma(k))

     ! Change of rhothil from (t) to (t+1)

     delex_rhothil(k) = delex_rhothil(k) &
                      + b10(k) * (vflux_thil(k-1) - vflux_thil(k))

     ! Update pressure from (t) to (t+tp)

     press(k,iw) = b5(k) * (rhothil(k) + pc3 * delex_rhothil(k))

     ! Update thil from (t) to (t+1)

     thil(k,iw) = ( rhothil(k) + delex_rhothil(k) ) / rho(k,iw)

     ! Update volume-weighted momentum tendencies at T due to implicit flux correction

     vmxet_rk(k,iw) = vmxet_rk(k,iw) + vflux_vxe(k-1) - vflux_vxe(k)
     vmyet_rk(k,iw) = vmyet_rk(k,iw) + vflux_vye(k-1) - vflux_vye(k)
     vmzet_rk(k,iw) = vmzet_rk(k,iw) + vflux_vze(k-1) - vflux_vze(k)

  enddo

  ! For shallow water test cases 2 & 5, rho & press are
  ! interpreted as water depth & height

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     topo_swtc = 0.

     if (nl%test_case == 5) then
        rad0_swtc = pi1 / 9.

        rad_swtc = sqrt((glonw(iw) * pio180 + 0.5 * pi1)**2 &

                 + (glatw(iw) * pio180 - pi1 / 6.) ** 2)

        topo_swtc = max(0., 2000. * (1. - rad_swtc / rad0_swtc))
     endif

     do k = ka,mza
        press(k,iw) = po_swtc(k) + topo_swtc
     enddo
  endif

  ! Vertical loop over W points

  do k = ka, mza-1
     wmc(k,iw) = wmf(k)
     wc (k,iw) = wmf(k) / real(zwgt_bot8(k) * rho(k,iw) + zwgt_top8(k) * rho(k+1,iw))
  enddo

  ! Set top & bottom values of WC

  wc(ka-1,iw) = wc(ka,iw)

  ! Bottom boundary for THIL

  do k = 1, ka-1
     thil(k,iw) = thil(ka,iw)
  enddo

end subroutine prog_wrt_begs

!============================================================================

subroutine prog_v_begs( iv, dts, vmc0, vmxet_rk, vmyet_rk, vmzet_rk, vmt_short )

  use mem_ijtabs,  only: itab_v
  use mem_basic,   only: vc, press, vmc, rho, pvfac
  use consts_coms, only: gravo2
  use mem_grid,    only: mza, mva, mwa, lpv, vnx, vny, vnz, volvi
  use oname_coms,  only: nl
!  use mem_rayf
  

  implicit none

  integer, intent(in) :: iv
  real,    intent(in) :: dts
  real,    intent(in) :: vmxet_rk(mza,mwa)
  real,    intent(in) :: vmyet_rk(mza,mwa)
  real,    intent(in) :: vmzet_rk(mza,mwa)

  real,    intent(in) :: vmt_short(mza,mva)
  real,    intent(in) :: vmc0     (mza,mva)

  integer :: k, kb
  integer :: iw1, iw2
  real    :: pgf(mza)

  ! Extract neighbor indices and coefficients for this point in the V stencil

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  kb = lpv(iv)

  ! Vertical loop over V points to compute pressure gradient force.

  do k = kb, mza
     pgf(k) = real( press(k,iw1) - press(k,iw2) ) * pvfac(k,iv)
  enddo

  ! For shallow water test cases 2 & 5, rho & press are interpreted
  ! as water depth & height.

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     do k = kb, mza
        pgf(k) = pgf(k) * real( rho(k,iw1) + rho(k,iw2) ) * gravo2
     enddo
  endif

  ! Update VMC

  do k = kb, mza

     vmc(k,iv) = vmc0(k,iv) + dts * ( vmt_short(k,iv) + pgf(k)   &
               + ( vnx(iv) * (vmxet_rk(k,iw1) + vmxet_rk(k,iw2)) &
                 + vny(iv) * (vmyet_rk(k,iw1) + vmyet_rk(k,iw2)) &
                 + vnz(iv) * (vmzet_rk(k,iw1) + vmzet_rk(k,iw2)) ) * volvi(k,iv) )

     vc(k,iv) = 2.0 * vmc(k,iv) / real( rho(k,iw1) + rho(k,iw2) )

  enddo

end subroutine prog_v_begs

!===========================================================================

subroutine fluxes_wrtv( scp, scp_upw, scp_upv, vmc, wmc, iorderh, i2d )

  use mem_grid,     only: mza, mwa, mva
  use mem_ijtabs,   only: jtab_v, jtab_w, jtv_wadj, jtw_prog
  use misc_coms,    only: iparallel
  use grad_lib,     only: grad_t2d, grad_t2d_quadratic
  use olam_mpi_atm, only: mpi_post_direct_recv_w, mpi_finish_direct_recv_w, &
                          mpi_post_direct_send_w, mpi_finish_direct_send_w
  use obnd,         only: lbcopy_w
  import,           only: itag_gxyps, iorderv

  implicit none

  integer, intent(in)  :: iorderh, i2d
  real,    intent(in)  :: scp    (mza,mwa)
  real,    intent(out) :: scp_upw(mza,mwa)
  real,    intent(out) :: scp_upv(mza,mva)
  real,    intent(in)  :: wmc    (mza,mwa)
  real,    intent(in)  :: vmc    (mza,mva)

  integer              :: j, iv, iw
  real                 :: gxyps(mza,mwa,i2d)

  ! COMPUTE HORIZONTAL POLYNOMIAL RECONSTRUCTION COEFFICIENTS AT EACH W CELL

  if (iorderh > 1) then
     if (iparallel == 1) call mpi_post_direct_recv_w(gxyps, itag_gxyps)

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (iorderh == 2) then
           call grad_t2d(iw, scp, gxyps(:,iw,1), gxyps(:,iw,2))
        elseif (iorderh == 3) then
           call grad_t2d_quadratic(iw, scp, gxyps, bounded=.false.)
        endif

     enddo
     !$omp end parallel do

     if (iparallel == 1) call mpi_post_direct_send_w(gxyps, itag_gxyps)
  endif

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH PRIMARY W FACE

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_w_column( iw )

  enddo
  !$omp end parallel do

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH V FACE

  if (iorderh > 1 .and. iparallel == 1) call mpi_finish_direct_recv_w(itag_gxyps)
  call lbcopy_w( aa = gxyps )

  !$omp parallel do private(iv)
  do j = 1, jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

     call scalar_v_column( iv )

  enddo
  !$omp end parallel do

  if (iorderh > 1 .and. iparallel == 1) call mpi_finish_direct_send_w(itag_gxyps)


contains


  subroutine scalar_v_column( iv )

    use mem_grid,   only: lpv, mza
    use mem_ijtabs, only: itab_v
    use mem_adv,    only: xx0_v, yy0_v, xy0_v
    import,         only: scp, scp_upv, vmc, gxyps, iorderh

    implicit none

    integer, intent(in) :: iv
    integer             :: iw1, iw2, k
    real                :: scp1(mza), scp2(mza)

    iw1 = itab_v(iv)%iw(1)
    iw2 = itab_v(iv)%iw(2)

    if (iorderh == 1) then  ! just for testing!

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1)
          scp2(k) = scp(k,iw2)
       enddo

    elseif (iorderh == 2) then

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1) &
                  + itab_v(iv)%dxps(1) * gxyps(k,iw1,1) &
                  + itab_v(iv)%dyps(1) * gxyps(k,iw1,2)

          scp2(k) = scp(k,iw2) &
                  + itab_v(iv)%dxps(2) * gxyps(k,iw2,1) &
                  + itab_v(iv)%dyps(2) * gxyps(k,iw2,2)
       enddo

    elseif (iorderh == 3) then

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1) &
                  + itab_v(iv)%dxps(1) * gxyps(k,iw1,1) &
                  + itab_v(iv)%dyps(1) * gxyps(k,iw1,2) &
                  + xx0_v(1,iv)        * gxyps(k,iw1,3) &
                  + xy0_v(1,iv)        * gxyps(k,iw1,4) &
                  + yy0_v(1,iv)        * gxyps(k,iw1,5)

          scp2(k) = scp(k,iw2)                          &
                  + itab_v(iv)%dxps(2) * gxyps(k,iw2,1) &
                  + itab_v(iv)%dyps(2) * gxyps(k,iw2,2) &
                  + xx0_v(2,iv)        * gxyps(k,iw2,3) &
                  + xy0_v(2,iv)        * gxyps(k,iw2,4) &
                  + yy0_v(2,iv)        * gxyps(k,iw2,5)
       enddo

    endif

    do k = lpv(iv), mza
       scp_upv(k,iv) = scp1(k)
       if (vmc(k,iv) < 0.) scp_upv(k,iv) = scp2(k)
    enddo

  end subroutine scalar_v_column


  subroutine scalar_w_column( iw )

    use mem_grid,   only: lpw, mza
    use grad_lib,   only: grad_z_linear, grad_z_quadratic
    import,         only: scp, wmc, scp_upw, iorderv

    implicit none

    integer, intent(in) :: iw
    integer             :: kb, k
    real                :: scpb(mza), scpt(mza)

    kb = lpw(iw)

    scp_upw(kb-1,iw) = scp(kb ,iw)
    scp_upw(mza ,iw) = scp(mza,iw)

    ! Scalar upwinded values at top and bottom faces (3rd order)

    if (iorderv == 1) then

       do k = kb, mza-1
          scpb(k) = scp(k,iw)
          scpt(k) = scp(k,iw)
       enddo

    elseif (iorderv == 2) then

       call grad_z_linear(iw, scp(:,iw), scpb, scpt)

    elseif (iorderv == 3) then

       call grad_z_quadratic(iw, scp(:,iw), scpb, scpt)

    endif

    do k = kb, mza-1
       scp_upw(k,iw) = scpt(k)
       if (wmc(k,iw) < 0.) scp_upw(k,iw) = scpb(k+1)
    enddo

  end subroutine scalar_w_column


end subroutine fluxes_wrtv

!=========================================================================

subroutine fluxes_thil_monot( scp, scp0, sct, scp_upw, scp_upv, &
                              rho, rho0, vmca, wmc, dtr, iorderh, i2d)

  use consts_coms,  only: r8
  use mem_grid,     only: mza, mwa, mva, lpv, lpw, volti
  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, jtv_wadj, jtw_prog
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_post_direct_recv_w, mpi_finish_direct_recv_w, &
                          mpi_post_direct_send_w, mpi_finish_direct_send_w
  use grad_lib,     only: grad_t2d, grad_t2d_quadratic
  use obnd,         only: lbcopy_w
  import,           only: itag_gxyps, iorderv, centered_monot, eps0, itag_scpt, &
                          itag_gxyps, itag_monot, onep, onem
  implicit none

  integer,  intent(in)  :: iorderh, i2d
  real,     intent(in)  :: scp    (mza,mwa)
  real,     intent(in)  :: scp0   (mza,mwa)
  real,     intent(in)  :: sct    (mza,mwa)
  real,     intent(out) :: scp_upw(mza,mwa)
  real,     intent(out) :: scp_upv(mza,mva)
  real(r8), intent(in)  :: rho    (mza,mwa)
  real(r8), intent(in)  :: rho0   (mza,mwa)
  real,     intent(in)  :: vmca   (mza,mva)
  real,     intent(in)  :: wmc    (mza,mwa)
  real,     intent(in)  :: dtr

  integer :: j, iv, iw, iw1, iw2, k
  real    :: scale

  real :: scpt       (mza,mwa)
  real :: scale_inout(mza,mwa,2)
  real :: gxyps      (mza,mwa,i2d)
  real :: sfluxvh    (mza,mva)
  real :: scp_hiv    (mza,mva)
  real :: sfluxwh    (mza,mwa)
  real :: scp_hiw    (mza,mwa)

  if (iparallel == 1) then
     call      mpi_post_direct_recv_w(scpt,        itag_scpt)
     if (iorderh > 1) &
          call mpi_post_direct_recv_w(gxyps,       itag_gxyps)
     call      mpi_post_direct_recv_w(scale_inout, itag_monot)
  endif

  !$omp parallel do private(iw,k)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza
        scpt(k,iw) = scp0(k,iw) + dtr * volti(k,iw) * sct(k,iw) / real(rho0(k,iw))
     enddo

  enddo
  !$omp end parallel do

  if (iparallel == 1) call mpi_post_direct_send_w(scpt, itag_scpt)

  ! COMPUTE HORIZONTAL POLYNOMIAL RECONSTRUCTION COEFFICIENTS AT EACH W CELL

  if (iorderh > 1) then

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (iorderh == 2) then
           call grad_t2d(iw, scp, gxyps(:,iw,1), gxyps(:,iw,2))
        elseif (iorderh == 3) then
           call grad_t2d_quadratic(iw, scp, gxyps)
        endif

     enddo
     !$omp end parallel do

     if (iparallel == 1) call mpi_post_direct_send_w(gxyps, itag_gxyps)
  endif

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH PRIMARY W FACE

  if (iparallel == 1) call mpi_finish_direct_recv_w(itag_scpt)

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_wflux_monot( iw )

  enddo
  !$omp end parallel do

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH V FACE

  if (iorderh > 1 .and. iparallel == 1) call mpi_finish_direct_recv_w(itag_gxyps)
  call lbcopy_w( aa = gxyps )

  !$omp parallel do private(iv)
  do j = 1, jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

     call scalar_vflux_monot( iv )

  enddo
  !$omp end parallel do

  if (iorderh > 1 .and. iparallel == 1) call mpi_finish_direct_send_w(itag_gxyps)

  ! COMPUTE_FLUX_LIMITER

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_monot_limiter( iw )

  enddo
  !$omp end parallel do

  if (iparallel == 1) call mpi_post_direct_send_w(scale_inout, itag_monot)

  ! APPLY FLUX RENORMALIZATION TO VERTICAL FLUXES

  !$omp parallel do private(iw,k,scale)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza-1
        scale = min( scale_inout(k+1,iw,1), scale_inout(k,iw,2) )
        if (sfluxwh(k,iw) < 0.) scale = min( scale_inout(k,iw,1), scale_inout(k+1,iw,2))
        scp_upw(k,iw) = scp_upw(k,iw) + scp_hiw(k,iw) * scale
     enddo

  enddo
  !$omp end parallel do

  ! APPLY FLUX RENORMALIZATION TO HORIZONTAL FLUXES

  if (iparallel == 1) call mpi_finish_direct_recv_w(itag_monot)
  call lbcopy_w( aa=scale_inout )

  !$omp parallel do private(iv,iw1,iw2,k,scale)
  do j = 1, jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     do k = lpv(iv), mza
        scale = min( scale_inout(k,iw2,1), scale_inout(k,iw1,2) )
        if (sfluxvh(k,iv) < 0.) scale = min( scale_inout(k,iw1,1), scale_inout(k,iw2,2) )
        scp_upv(k,iv) = scp_upv(k,iv) + scp_hiv(k,iv) * scale
     enddo

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call      mpi_finish_direct_send_w(itag_scpt)
     if (iorderh > 1) &
          call mpi_finish_direct_send_w(itag_gxyps)
     call      mpi_finish_direct_send_w(itag_monot)
  endif


contains


  subroutine scalar_vflux_monot( iv )

    use mem_grid,   only: lpv, mza
    use mem_ijtabs, only: itab_v
    use mem_adv,    only: xx0_vu, yy0_vu, xy0_vu
    import,         only: scp, scpt, scp_upv, scp_hiv, vmca, sfluxvh, gxyps, &
                          iorderh, centered_monot
    implicit none

    integer, intent(in) :: iv
    integer             :: iw1, iw2, k
    real                :: scp1(mza), scp2(mza)

    iw1 = itab_v(iv)%iw(1)
    iw2 = itab_v(iv)%iw(2)

     ! Low-order horizontal fluxes

    do k = lpv(iv), mza
       scp_upv(k,iv) = scpt(k,iw1)
       if (vmca(k,iv) < 0.) scp_upv(k,iv) = scpt(k,iw2)
    enddo

    ! High-order horizontal fluxes

    if (iorderh == 1) then

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1)
          scp2(k) = scp(k,iw2)
       enddo

    elseif (iorderh == 2) then

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1) &
                  + itab_v(iv)%dxps(1) * gxyps(k,iw1,1) &
                  + itab_v(iv)%dyps(1) * gxyps(k,iw1,2)

          scp2(k) = scp(k,iw2) &
                  + itab_v(iv)%dxps(2) * gxyps(k,iw2,1) &
                  + itab_v(iv)%dyps(2) * gxyps(k,iw2,2)
       enddo

    elseif (iorderh == 3) then

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1) &
                  + itab_v(iv)%dxps(1) * gxyps(k,iw1,1) &
                  + itab_v(iv)%dyps(1) * gxyps(k,iw1,2) &
                  + xx0_vu(1,iv)       * gxyps(k,iw1,3) &
                  + xy0_vu(1,iv)       * gxyps(k,iw1,4) &
                  + yy0_vu(1,iv)       * gxyps(k,iw1,5)

          scp2(k) = scp(k,iw2)                          &
                  + itab_v(iv)%dxps(2) * gxyps(k,iw2,1) &
                  + itab_v(iv)%dyps(2) * gxyps(k,iw2,2) &
                  + xx0_vu(2,iv)       * gxyps(k,iw2,3) &
                  + xy0_vu(2,iv)       * gxyps(k,iw2,4) &
                  + yy0_vu(2,iv)       * gxyps(k,iw2,5)
       enddo

    endif

    if (centered_monot) then

       do k = lpv(iv), mza
          scp_hiv(k,iv) = 0.5 * (scp1(k) + scp2(k)) - scp_upv(k,iv)
          sfluxvh(k,iv) = vmca(k,iv) * scp_hiv(k,iv)
       enddo

    else

       do k = lpv(iv), mza
          scp_hiv(k,iv) = scp1(k)
          if (vmca(k,iv) < 0.) scp_hiv(k,iv) = scp2(k)

          scp_hiv(k,iv) = scp_hiv(k,iv) - scp_upv(k,iv)
          sfluxvh(k,iv) = vmca(k,iv) * scp_hiv(k,iv)
       enddo

    endif

  end subroutine scalar_vflux_monot



  subroutine scalar_wflux_monot( iw )

    use mem_grid,   only: lpw, lpv, mza, volti, arw
    use grad_lib,   only: grad_z_linear, grad_z_quadratic
    import,         only: scp, scpt, scp_upw, scp_hiw, sfluxwh, wmc, &
                          centered_monot, iorderv
    implicit none

    integer, intent(in) :: iw

    integer             :: kb, k
    real                :: scp1(mza), scp2(mza)

    kb = lpw(iw)

    ! Low order vertical flux

    scp_upw(kb-1,iw) = scp(kb ,iw)
    scp_upw(mza ,iw) = scp(mza,iw)

    do k = kb, mza-1
       scp_upw(k,iw) = scpt(k,iw)
       if (wmc(k,iw) < 0.) scp_upw(k,iw) = scpt(k+1,iw)
    enddo

    ! Higher-order vertical fluxes

    if (iorderv == 1) then

       do k = kb, mza-1
          scp1(k) = scp(k,iw)
          scp2(k) = scp(k,iw)
       enddo

    elseif (iorderv == 2) then

       call grad_z_linear(iw, scp(:,iw), scp1, scp2, itopbc=1, ibotbc=1)

    elseif (iorderv == 3) then

       call grad_z_quadratic(iw, scp(:,iw), scp1, scp2, itopbc=1, ibotbc=1, bounded=.false.)

    endif

    sfluxwh(kb-1,iw) = 0.
    sfluxwh(mza ,iw) = 0.

    if (centered_monot) then

       do k = kb, mza-1
          scp_hiw(k,iw) = 0.5 * (scp2(k) + scp1(k+1)) - scp_upw(k,iw)
          sfluxwh(k,iw) = arw(k,iw) * wmc(k,iw) * scp_hiw(k,iw)
       enddo

    else

       do k = kb, mza-1
          scp_hiw(k,iw) = scp2(k)
          if (wmc(k,iw) < 0.) scp_hiw(k,iw) = scp1(k+1)

          scp_hiw(k,iw) = scp_hiw(k,iw) - scp_upw(k,iw)
          sfluxwh(k,iw) = arw(k,iw) * wmc(k,iw) * scp_hiw(k,iw)
       enddo

    endif

  end subroutine scalar_wflux_monot



  subroutine scalar_monot_limiter( iw )

    use mem_grid,   only: lpw, lpv, mza, volti, arw
    use mem_ijtabs, only: itab_w
    import,         only: scp, scpt, vmca, wmc, scp_upv, scp_upw, scp_hiv, scp_hiw, &
                          sfluxvh, sfluxwh, scale_inout, rho, dtr, eps0, onep, onem
    implicit none

    integer, intent(in) :: iw

    integer             :: kb, k, jv, iv, iwn
    real                :: rscp_low, flux_in, flux_out, sc_in, sc_out
    real                :: scale, dtrp
    real                :: scp1(mza), scp2(mza)
    real                :: smax(mza), smin(mza)
    real                :: fluxdiv_low(mza)
    real                :: fluxdiv_in (mza)
    real                :: fluxdiv_out(mza)
    real                :: wmca(mza), rhos(mza)

    dtrp = onep * dtr

    kb = lpw(iw)

    wmca(kb-1) = 0.
    do k = kb, mza-1
       wmca(k) = wmc(k,iw) * arw(k,iw)
    enddo
    wmca(mza) = 0.

     ! Find upwind-biased tracer max/min for each cell in this column

    do k = kb, mza
       smax(k) = scpt(k,iw)
       smin(k) = scpt(k,iw)
       rhos(k) = rho(k,iw)
    enddo

    do k = kb+1, mza
       if (wmca(k-1) > 0.) then
          smax(k) = max(smax(k), scpt(k-1,iw))
          smin(k) = min(smin(k), scpt(k-1,iw))
       endif
    enddo

    do k = kb, mza-1
       if (wmca(k) < 0.) then
          smax(k) = max(smax(k), scpt(k+1,iw))
          smin(k) = min(smin(k), scpt(k+1,iw))
       endif
    enddo

     ! Scalar horizontal flux divergence and upwind-biased max/min

     fluxdiv_low(:) = 0.
     fluxdiv_in (:) = 0.
     fluxdiv_out(:) = 0.

     do jv = 1, itab_w(iw)%npoly
        iv  = itab_w(iw)%iv(jv)
        iwn = itab_w(iw)%iw(jv)

        do k = lpv(iv), mza

           fluxdiv_low(k) = fluxdiv_low(k) &
                          + itab_w(iw)%dirv(jv) * vmca(k,iv) * (scp_upv(k,iv) - scpt(k,iw))

           fluxdiv_in (k) = fluxdiv_in (k) + max( itab_w(iw)%dirv(jv) * sfluxvh(k,iv), 0.)
           fluxdiv_out(k) = fluxdiv_out(k) + min( itab_w(iw)%dirv(jv) * sfluxvh(k,iv), 0.)

           ! Upwind-biased tracer max/min
           if ( itab_w(iw)%dirv(jv) * vmca(k,iv) > 0.) then
              smax(k) = max(smax(k), scpt(k,iwn))
              smin(k) = min(smin(k), scpt(k,iwn))
           endif

        enddo
     enddo

     ! Scalar vertical flux divergence and flux limiters

     do k = kb, mza
        fluxdiv_low(k) = fluxdiv_low(k) + wmca(k-1) * (scp_upw(k-1,iw) - scpt(k,iw)) &
                                        - wmca(k  ) * (scp_upw(k  ,iw) - scpt(k,iw))

        fluxdiv_in (k) = fluxdiv_in (k) + max(sfluxwh(k-1,iw),0.) - min(sfluxwh(k,iw),0.)
        fluxdiv_out(k) = fluxdiv_out(k) + min(sfluxwh(k-1,iw),0.) - max(sfluxwh(k,iw),0.)

        rscp_low = scpt(k,iw) * rhos(k) + dtr * volti(k,iw) * fluxdiv_low(k)

        flux_in  = max(dtrp * volti(k,iw) * fluxdiv_in (k),  eps0)
        flux_out = min(dtrp * volti(k,iw) * fluxdiv_out(k), -eps0)

        sc_in  = max(rhos(k) * smax(k) - onep * rscp_low - eps0, 0.)
        sc_out = min(rhos(k) * smin(k) - onem * rscp_low + eps0, 0.)

        scale_inout(k,iw,1) = min(1., sc_in  / max(flux_in , 0.01*sc_in ) )
        scale_inout(k,iw,2) = min(1., sc_out / min(flux_out, 0.01*sc_out) )
     enddo

   end subroutine scalar_monot_limiter

 end subroutine fluxes_thil_monot

!=========================================================================

subroutine update_vxe_undrgnd( iw, dts, vmxe0, vmye0, vmze0, &
                               vmxet_rk, vmyet_rk, vmzet_rk )

  use mem_basic,  only: vxe, vye, vze, rho
  use mem_grid,   only: mwa, mza, lpw, lve2, nve2_max, volti

  implicit none

  integer, intent(in) :: iw
  real,    intent(in) :: dts

  real,    intent(in) :: vmxe0(nve2_max,mwa)
  real,    intent(in) :: vmye0(nve2_max,mwa)
  real,    intent(in) :: vmze0(nve2_max,mwa)

  real,    intent(in) :: vmxet_rk(mza,mwa)
  real,    intent(in) :: vmyet_rk(mza,mwa)
  real,    intent(in) :: vmzet_rk(mza,mwa)

  integer :: ka, ks, k
  real    :: rsi

  if (lve2(iw) > 0) then
     ka = lpw(iw)

     do ks = 1, lve2(iw)
        k   = ks + ka - 1
        rsi = 1.0 / real(rho(k,iw))

        vxe(k,iw) = (vmxe0(ks,iw) + dts * ( vmxet_rk(k,iw) * volti(k,iw) ) ) * rsi
        vye(k,iw) = (vmye0(ks,iw) + dts * ( vmyet_rk(k,iw) * volti(k,iw) ) ) * rsi
        vze(k,iw) = (vmze0(ks,iw) + dts * ( vmzet_rk(k,iw) * volti(k,iw) ) ) * rsi
     enddo

  endif

end subroutine update_vxe_undrgnd

!=========================================================================

end module wrtv_rk
