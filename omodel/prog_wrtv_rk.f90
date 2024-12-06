module wrtv_rk

  use consts_coms, only: r8
  implicit none

  real, allocatable :: b20(:)

  real(r8), parameter :: onethird = 1.0_r8 / 3.0_r8

  real(r8), parameter :: xstg3(3) = [ onethird, 0.5_r8, 1.0_r8 ]
  real(r8), parameter :: xstg2(2) =           [ 0.5_r8, 1.0_r8 ]

  private
  public init_wrtv_rk, prog_wrtv_rk

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
  use mem_basic,    only: rho, thil, wc, press, vmc, vc, wmc, &
                          vxe, vye, vze, vmasc
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, dzto2, dztsqo6, volt, xev, &
                          nve2_max, lve2, arv
  use mem_tend,     only: thilt, vmxet, vmyet, vmzet, vmt
  use misc_coms,    only: iparallel, dtsm, dtlm, nrk_wrtv, dn01d, th01d, &
                          initial, nrk_scal, deltax, nxp, mdomain
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v
  use obnd,         only: lbcopy_m, lbcopy_w, lbcopy_v
  use oname_coms,   only: nl
  use mem_adv,      only: gxps_scp, gyps_scp, gxyps_scp, gxxps_scp, gyyps_scp, &
                          xx0_v, yy0_v, xy0_v
  use grad_lib,     only: grad_z_quad, grad_t2d, grad_t2d_quad
  use vel_t3d,      only: vel_t3d_hex
  use consts_coms,  only: p00i, rocp, r8
  use mem_turb,     only: akmodx, akhodx, khtopv, kmtopv
  use pbl_drivers,  only: solve_eddy_diff_heat, solve_eddy_diff_vxe
  use mem_rayf,     only: dorayf, krayf_bot, rayf_cof, dorayfmix, rayf_mix_top_vxe, &
                          vc03d
  use mem_nudge,    only: vmanud, nudflag

  implicit none

  integer  :: j, iv, iw, k, ksw, iw1, iw2, iwn, jv, istage
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

  real(r8) :: rho0(mza,mwa)
  real(r8) :: rth0(mza,mwa)
  real     :: wmc0(mza,mwa)
  real     :: vmc0(mza,mva)

  real :: gzps_th (mza), gzps_vx (mza), gzps_vy (mza), gzps_vz (mza)
  real :: gzzps_th(mza), gzzps_vx(mza), gzzps_vy(mza), gzzps_vz(mza)

  real :: thilt_short(mza,mwa)
  real :: vmxet_short(mza,mwa)
  real :: vmyet_short(mza,mwa)
  real :: vmzet_short(mza,mwa)
  real :: vmt_short  (mza,mva)

  ! For prognostic shaved-cell method
  real :: vmtrk(nve2_max,mva)
  real :: vmxe0(nve2_max,mwa)
  real :: vmye0(nve2_max,mwa)
  real :: vmze0(nve2_max,mwa)

  real :: gxps_vxe(mza,mwa)
  real :: gyps_vxe(mza,mwa)
  real :: gxps_vye(mza,mwa)
  real :: gyps_vye(mza,mwa)
  real :: gxps_vze(mza,mwa)
  real :: gyps_vze(mza,mwa)

  real :: unit_dist, fracx, rayfx

  real(r8) :: vmca(mza,mva)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  real :: gxxps_vxe(mza,mwa)
!!  real :: gyyps_vxe(mza,mwa)
!!  real :: gxyps_vxe(mza,mwa)
!!
!!  real :: gxxps_vye(mza,mwa)
!!  real :: gyyps_vye(ma,mwa)
!!  real :: gxyps_vye(mza,mwa)
!!
!!  real :: gxxps_vze(mza,mwa)
!!  real :: gyyps_vze(mza,mwa)
!!  real :: gxyps_vze(mza,mwa)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  !$omp end do
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

     !$omp parallel do private(iw,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        do k = lpw(iw), mza
           vmxet_rk(k,iw) = vmxet_short(k,iw)
           vmyet_rk(k,iw) = vmyet_short(k,iw)
           vmzet_rk(k,iw) = vmzet_short(k,iw)
        enddo

        if (nl%horiz_adv_order <= 2) then

           call grad_t2d(iw, thil, gxps_scp(:,iw),  gyps_scp(:,iw))

           if (nl%iscal_monot > 0) then
              call limit_h2(iw, thil, gxps_scp(:,iw), gyps_scp(:,iw))
           endif

        else

           call grad_t2d_quad(iw, thil, gxps_scp (:,iw), gyps_scp (:,iw), &
                              gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))

           if (nl%iscal_monot > 0) then
              call limit_h3(iw, thil, gxps_scp(:,iw), gyps_scp(:,iw), &
                            gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))
           endif

        endif

        ! velocity second order for now
        call grad_t2d(iw, vxe, gxps_vxe(:,iw), gyps_vxe(:,iw))
        call grad_t2d(iw, vye, gxps_vye(:,iw), gyps_vye(:,iw))
        call grad_t2d(iw, vze, gxps_vze(:,iw), gyps_vze(:,iw))

!!      call grad_t2d_quad(iw, vxe, gxps_vxe (:,iw), gyps_vxe (:,iw), &
!!                         gxxps_vxe(:,iw), gxyps_vxe(:,iw), gyyps_vxe(:,iw))
!!
!!      call grad_t2d_quad(iw, vye, gxps_vye (:,iw), gyps_vye(:,iw), &
!!                         gxxps_vye(:,iw), gxyps_vye(:,iw), gyyps_vye(:,iw))
!!
!!      call grad_t2d_quad(iw, vze, gxps_vze (:,iw), gyps_vze(:,iw), &
!!                         gxxps_vze(:,iw), gxyps_vze(:,iw), gyyps_vze(:,iw))

!!      call grad_t2d_v(iw, vxe_upv, vye_upv, vze_upv,  &
!!                      gxps_vxe(:,iw), gyps_vxe(:,iw), &
!!                      gxps_vye(:,iw), gyps_vye(:,iw), &
!!                      gxps_vze(:,iw), gyps_vze(:,iw)  )

     enddo
     !$omp end parallel do

     ! MPI send of THIL, VXE, VYE, and VZE gradient components

     if (iparallel == 1) then

        if (nl%horiz_adv_order <= 2) then
           call mpi_send_w(rvara1=gxps_scp, rvara2=gyps_scp, &
                           rvara3=gxps_vxe, rvara4=gyps_vxe, &
                           rvara5=gxps_vye, rvara6=gyps_vye, &
                           rvara7=gxps_vze, rvara8=gyps_vze  )
        else
           call mpi_send_w(rvara1=gxps_scp,   rvara2=gyps_scp, &
                           rvara3=gxxps_scp,  rvara4=gxyps_scp, rvara5=gyyps_scp, &
                           rvara6=gxps_vxe,   rvara7=gyps_vxe, &
                           rvara8=gxps_vye,   rvara9=gyps_vye, &
                           rvara10=gxps_vze,  rvara11=gyps_vze )
        endif

!!      call mpi_send_w(rvara1=gxps_scp, rvara2=gyps_scp, &
!!                      rvara3=gxps_vxe, rvara4=gyps_vxe, &
!!                      rvara5=gxps_vye, rvara6=gyps_vye, &
!!                      rvara7=gxps_vze, rvara8=gyps_vze, &
!!                      rvara9=gxxps_vxe, rvara10=gyyps_vxe, rvara11=gxyps_vxe, &
!!                      rvara12=gxxps_vye, rvara13=gyyps_vye,  rvara14=gxyps_vye,&
!!                      rvara15=gxxps_vze, rvara16=gyyps_vze,  rvara17=gxyps_vze)
     endif

     !$omp parallel private(gzps_th, gzps_vx, gzps_vy, gzps_vz,&
     !$omp                  gzzps_th,gzzps_vx,gzzps_vy,gzzps_vz)
     !$omp do private(iw,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Evaluate vertical gradient of thil, vxe, vye, and vze for begs

        call grad_z_quad(iw, thil(:,iw), gzps_th, gzzps_th)
        call grad_z_quad(iw, vxe (:,iw), gzps_vx, gzzps_vx)
        call grad_z_quad(iw, vye (:,iw), gzps_vy, gzzps_vy)
        call grad_z_quad(iw, vze (:,iw), gzps_vz, gzzps_vz)

        if (nl%iscal_monot > 0) then
           call limit_z(iw, thil(:,iw), gzps_th, gzzps_th)
        endif

        ! Evaluate upwinded thil, vxe, vye, and vze at each w face

        do k = lpw(iw), mza-1

           if (wmc(k,iw) > 0.) then

              thil_upw(k,iw) = thil(k,iw)  &
                             +   dzto2(k) *  gzps_th(k) &
                             + dztsqo6(k) * gzzps_th(k)

              vxe_upw(k,iw)  = vxe(k,iw)                &
                             +   dzto2(k) *  gzps_vx(k) &
                             + dztsqo6(k) * gzzps_vx(k)

              vye_upw(k,iw)  = vye(k,iw)                &
                             +   dzto2(k) *  gzps_vy(k) &
                             + dztsqo6(k) * gzzps_vy(k)

              vze_upw(k,iw)  = vze(k,iw)                &
                             +   dzto2(k) *  gzps_vz(k) &
                             + dztsqo6(k) * gzzps_vz(k)
           else

              thil_upw(k,iw) = thil(k+1,iw)                 &
                             -   dzto2(k+1) *  gzps_th(k+1) &
                             + dztsqo6(k+1) * gzzps_th(k+1)

              vxe_upw(k,iw) = vxe(k+1,iw)                  &
                            -   dzto2(k+1) *  gzps_vx(k+1) &
                            + dztsqo6(k+1) * gzzps_vx(k+1)

              vye_upw(k,iw) = vye(k+1,iw)                  &
                            -   dzto2(k+1) *  gzps_vy(k+1) &
                            + dztsqo6(k+1) * gzzps_vy(k+1)

              vze_upw(k,iw) = vze(k+1,iw)                  &
                            -   dzto2(k+1) *  gzps_vz(k+1) &
                            + dztsqo6(k+1) * gzzps_vz(k+1)
           endif

        enddo

     enddo
     !$omp end do nowait
     !$omp end parallel

! Finish MPI recv of THIL, VXE, VYE, and VZE gradient components

     if (iparallel == 1) then

        if (nl%horiz_adv_order <= 2) then
           call mpi_recv_w(rvara1=gxps_scp, rvara2=gyps_scp, &
                           rvara3=gxps_vxe, rvara4=gyps_vxe, &
                           rvara5=gxps_vye, rvara6=gyps_vye, &
                           rvara7=gxps_vze, rvara8=gyps_vze  )
        else
           call mpi_recv_w(rvara1=gxps_scp,   rvara2=gyps_scp, &
                           rvara3=gxxps_scp,  rvara4=gxyps_scp, rvara5=gyyps_scp, &
                           rvara6=gxps_vxe,   rvara7=gyps_vxe, &
                           rvara8=gxps_vye,   rvara9=gyps_vye, &
                           rvara10=gxps_vze,  rvara11=gyps_vze )
        endif

!!      call mpi_recv_w(rvara1=gxps_scp, rvara2=gyps_scp, &
!!                      rvara3=gxps_vxe, rvara4=gyps_vxe, &
!!                      rvara5=gxps_vye, rvara6=gyps_vye, &
!!                      rvara7=gxps_vze, rvara8=gyps_vze, &
!!                      rvara9=gxxps_vxe, rvara10=gyyps_vxe, rvara11=gxyps_vxe, &
!!                      rvara12=gxxps_vye, rvara13=gyyps_vye,  rvara14=gxyps_vye,&
!!                      rvara15=gxxps_vze, rvara16=gyyps_vze,  rvara17=gxyps_vze)
     endif

     ! Lateral boundary copy of THIL, VXE, VYE, and VZE gradient components

     if (nl%horiz_adv_order <= 2) then
        call lbcopy_w(a1=gxps_scp, a2=gyps_scp, &
                      a3=gxps_vxe, a4=gyps_vxe, &
                      a5=gxps_vye, a6=gyps_vye, &
                      a7=gxps_vze, a8=gyps_vze  )
     else
        call lbcopy_w(a1=gxps_scp,  a2=gyps_scp, &
                      a3=gxxps_scp, a4=gxyps_scp,  a5=gyyps_scp,  &
                      a6=gxps_vxe,  a7=gyps_vxe, &
                      a8=gxps_vye,  a9=gyps_vye, &
                      a10=gxps_vze, a11=gyps_vze )
     endif

     ! EVALUATE UPWINDED THIL, VXE, VYE, AND VZE AT EACH V FACE FOR
     ! HORIZONTAL FLUX COMPUTATION

     !$omp parallel
     !$omp do private(iv,k,iw1,iw2)
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

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        do k = lpv(iv), mza

           if (vmc(k,iv) > 0.) then

              if (nl%horiz_adv_order <= 2) then
                 thil_upv(k,iv) = thil(k,iw1) &
                                + itab_v(iv)%dxps(1) * gxps_scp(k,iw1) &
                                + itab_v(iv)%dyps(1) * gyps_scp(k,iw1)
              else
                 thil_upv(k,iv) = thil(k,iw1) &
                                + itab_v(iv)%dxps(1) *  gxps_scp(k,iw1) &
                                + itab_v(iv)%dyps(1) *  gyps_scp(k,iw1) &
                                + xx0_v(1,iv)        * gxxps_scp(k,iw1) &
                                + yy0_v(1,iv)        * gyyps_scp(k,iw1) &
                                + xy0_v(1,iv)        * gxyps_scp(k,iw1)
              endif

              vxe_upv(k,iv)  = vxe(k,iw1) &
                             + itab_v(iv)%dxps(1) * gxps_vxe(k,iw1) &
                             + itab_v(iv)%dyps(1) * gyps_vxe(k,iw1)

              vye_upv(k,iv)  = vye(k,iw1) &
                             + itab_v(iv)%dxps(1) * gxps_vye(k,iw1) &
                             + itab_v(iv)%dyps(1) * gyps_vye(k,iw1)

              vze_upv(k,iv)  = vze(k,iw1) &
                             + itab_v(iv)%dxps(1) * gxps_vze(k,iw1) &
                             + itab_v(iv)%dyps(1) * gyps_vze(k,iw1)

!!            vxe_upv(k,iv)  = vxe(k,iw1) &
!!                           + itab_v(iv)%dxps(1) * gxps_vxe(k,iw1) &
!!                           + itab_v(iv)%dyps(1) * gyps_vxe(k,iw1) &
!!                           + xx0_v(1,iv)        * gxxps_vxe(k,iw1) &
!!                           + yy0_v(1,iv)        * gyyps_vxe(k,iw1) &
!!                           + xy0_v(1,iv)        * gxyps_vxe(k,iw1)
!!
!!            vye_upv(k,iv)  = vye(k,iw1) &
!!                           + itab_v(iv)%dxps(1) * gxps_vye(k,iw1) &
!!                           + itab_v(iv)%dyps(1) * gyps_vye(k,iw1) &
!!                           + xx0_v(1,iv)        * gxxps_vye(k,iw1) &
!!                           + yy0_v(1,iv)        * gyyps_vye(k,iw1) &
!!                           + xy0_v(1,iv)        * gxyps_vye(k,iw1)
!!
!!            vze_upv(k,iv)  = vze(k,iw1) &
!!                           + itab_v(iv)%dxps(1) * gxps_vze(k,iw1) &
!!                           + itab_v(iv)%dyps(1) * gyps_vze(k,iw1) &
!!                           + xx0_v(1,iv)        * gxxps_vze(k,iw1) &
!!                           + yy0_v(1,iv)        * gyyps_vze(k,iw1) &
!!                           + xy0_v(1,iv)        * gxyps_vze(k,iw1)

           else

              if (nl%horiz_adv_order <= 2) then
                 thil_upv(k,iv) = thil(k,iw2) &
                                + itab_v(iv)%dxps(2) * gxps_scp(k,iw2) &
                                + itab_v(iv)%dyps(2) * gyps_scp(k,iw2)
              else
                 thil_upv(k,iv) = thil(k,iw2) &
                                + itab_v(iv)%dxps(2) *  gxps_scp(k,iw2) &
                                + itab_v(iv)%dyps(2) *  gyps_scp(k,iw2) &
                                + xx0_v(2,iv)        * gxxps_scp(k,iw2) &
                                + yy0_v(2,iv)        * gyyps_scp(k,iw2) &
                                + xy0_v(2,iv)        * gxyps_scp(k,iw2)
              endif

              vxe_upv(k,iv)  = vxe(k,iw2) &
                             + itab_v(iv)%dxps(2) * gxps_vxe(k,iw2) &
                             + itab_v(iv)%dyps(2) * gyps_vxe(k,iw2)

              vye_upv(k,iv)  = vye(k,iw2) &
                             + itab_v(iv)%dxps(2) * gxps_vye(k,iw2) &
                             + itab_v(iv)%dyps(2) * gyps_vye(k,iw2)

              vze_upv(k,iv)  = vze(k,iw2) &
                             + itab_v(iv)%dxps(2) * gxps_vze(k,iw2) &
                             + itab_v(iv)%dyps(2) * gyps_vze(k,iw2)

!!            vxe_upv(k,iv)  = vxe(k,iw2) &
!!                           + itab_v(iv)%dxps(2) * gxps_vxe(k,iw2) &
!!                           + itab_v(iv)%dyps(2) * gyps_vxe(k,iw2) &
!!                           + xx0_v(2,iv)        * gxxps_vxe(k,iw1) &
!!                           + yy0_v(2,iv)        * gyyps_vxe(k,iw1) &
!!                           + xy0_v(2,iv)        * gxyps_vxe(k,iw1)
!!
!!            vye_upv(k,iv)  = vye(k,iw2) &
!!                           + itab_v(iv)%dxps(2) * gxps_vye(k,iw2) &
!!                           + itab_v(iv)%dyps(2) * gyps_vye(k,iw2) &
!!                           + xx0_v(2,iv)        * gxxps_vye(k,iw1) &
!!                           + yy0_v(2,iv)        * gyyps_vye(k,iw1) &
!!                           + xy0_v(2,iv)        * gxyps_vye(k,iw1)
!!
!!            vze_upv(k,iv)  = vze(k,iw2) &
!!                           + itab_v(iv)%dxps(2) * gxps_vze(k,iw2) &
!!                           + itab_v(iv)%dyps(2) * gyps_vze(k,iw2) &
!!                           + xx0_v(2,iv)        * gxxps_vze(k,iw1) &
!!                           + yy0_v(2,iv)        * gyyps_vze(k,iw1) &
!!                           + xy0_v(2,iv)        * gxyps_vze(k,iw1)

           endif

        enddo

     enddo
     !$omp end do

     ! MAIN LOOP OVER W COLUMNS FOR UPDATING WM, WC, RHO, THIL, AND PRESS

     !$omp do private(iw,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        ! Prognose vertical velocity, density, thil, and diagnose pressure

        call prog_wrt_begs( iw, istage, dt8,                      &
                            thil_upv, vxe_upv, vye_upv, vze_upv,  &
                            thil_upw, vxe_upw, vye_upw, vze_upw,  &
                            thilt_short, vmxet_rk, vmyet_rk, vmzet_rk, &
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
           call update_vxe_undrgnd(iw, dts, vmxe0, vmye0, vmze0, vmtrk, &
                                   vmxet_rk, vmyet_rk, vmzet_rk)
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
                          thilt_short, vmxet_rk, vmyet_rk, vmzet_rk, &
                          rho0, rth0, wmc0, vmca )

  use mem_ijtabs,  only: itab_w
  use mem_basic,   only: wmc, rho, thil, wc, press, vxe, vye, wmasc, &
                         alpha_press, pwfac
  use misc_coms,   only: nrk_wrtv, mdomain
  use consts_coms, only: cpocv, rocv, fcoriol, pi1, pio180, r8
  use mem_grid,    only: mza, mva, mwa, lpv, lpw, arw, wnx, wny, wnz, volt, &
                         volti, volwi, glatw, glonw, arv, &
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

  real, intent(in)    :: thilt_short(mza,mwa)
  real, intent(inout) :: vmxet_rk(mza,mwa)
  real, intent(inout) :: vmyet_rk(mza,mwa)
  real, intent(inout) :: vmzet_rk(mza,mwa)

  real(r8), intent(in) :: rth0(mza,mwa)
  real(r8), intent(in) :: rho0(mza,mwa)
  real,     intent(in) :: wmc0(mza,mwa)
  real(r8), intent(in) :: vmca(mza,mva)

  integer :: jv, iv
  integer :: k, ka, kbv, ksw
  integer :: npoly

  real :: c8, c9, dts
  real :: rad0_swtc, rad_swtc, topo_swtc

  ! Vertical implicit scheme weighting parameters

  real(r8), parameter :: fr  = .55 ! rho
  real(r8), parameter :: fp2 = .55 ! press in W
  real(r8), parameter :: fp3 = .70 ! press in V

  real(r8), parameter :: pc2 = fp2 * cpocv
  real(r8), parameter :: pc3 = fp3 * cpocv

  ! Automatic arrays

  real(r8) :: wmarw    (mza)
  real(r8) :: hflx_rho (mza)
  real(r8) :: delex_rho(mza)
  real(r8) :: delex_rhothil(mza)
  real(r8) :: po_swtc(mza)

  real(r8) :: rhothil   (mza)
  real(r8) :: thilt_rk  (mza)
  real(r8) :: vflux_thil(mza)

  real(r8) :: press_ex(mza)
  real(r8) :: rho_ex  (mza)

  real :: delex_wm(mza)
  real :: wmf     (mza)
  real :: mass, wma

  real :: vflux_vxe(mza)
  real :: vflux_vye(mza)
  real :: vflux_vze(mza)

  real(r8) :: b5(mza), b10(mza)

  real :: b7, b8, b21, b22, b23, b24, b25, b26
  real :: b15(mza), b31(mza), b32(mza), b33(mza)

  dts = real(dt8)

  ka = lpw(iw)

  ! Set bottom & top vertical advective mass, momentum, and heat fluxes to zero

  wmarw(ka-1) = 0._r8
  wmarw(mza)  = 0._r8

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
        wmarw     (k) = wmanud(k,iw)
        vflux_thil(k) = wmarw(k) * thil_upw(k,iw)

        vflux_vxe (k) = wmanud(k,iw) * vxe_upw(k,iw)
        vflux_vye (k) = wmanud(k,iw) * vye_upw(k,iw)
        vflux_vze (k) = wmanud(k,iw) * vze_upw(k,iw)
     enddo

     do k = lpw(iw), mza
        hflx_rho(k)    =                     wmarw     (k-1) - wmarw     (k)
        thilt_rk(k)    = thilt_short(k,iw) + vflux_thil(k-1) - vflux_thil(k)
        vmxet_rk(k,iw) = vmxet_rk(k,iw)    + vflux_vxe (k-1) - vflux_vxe (k)
        vmyet_rk(k,iw) = vmyet_rk(k,iw)    + vflux_vye (k-1) - vflux_vye (k)
        vmzet_rk(k,iw) = vmzet_rk(k,iw)    + vflux_vze (k-1) - vflux_vze (k)
     enddo

  else

     do k = ka, mza
        thilt_rk(k) = thilt_short(k,iw)
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
           thilt_rk(k)    = thilt_rk(k)    + vmca(k,iv) * thil_upv(k,iv)
           vmxet_rk(k,iw) = vmxet_rk(k,iw) + vmca(k,iv) * vxe_upv (k,iv)
           vmyet_rk(k,iw) = vmyet_rk(k,iw) + vmca(k,iv) * vye_upv (k,iv)
           vmzet_rk(k,iw) = vmzet_rk(k,iw) + vmca(k,iv) * vze_upv (k,iv)
        enddo

     else

        do k = kbv, mza
           hflx_rho(k)    = hflx_rho(k)    - vmca(k,iv)
           thilt_rk(k)    = thilt_rk(k)    - vmca(k,iv) * thil_upv(k,iv)
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
     delex_rhothil(k) = b10(k) * thilt_rk(k)
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

  do k = ka, mza-1
     wma      = wmf(k) * arw(k,iw)
     wmarw(k) = wma

     vflux_thil(k) = wmarw(k) * thil_upw(k,iw)
     vflux_vxe (k) = wma      * vxe_upw (k,iw)
     vflux_vye (k) = wma      * vye_upw (k,iw)
     vflux_vze (k) = wma      * vze_upw (k,iw)
  enddo

  ! Add vertical momentum at (t+fw) to array for long-timestep scalar transport

  if (istage == nrk_wrtv) then

     if (nudflag > 0 .and. nl%nud_preserve_total_mass) then
        do k = ka, mza-1
           wmasc(k,iw) = wmasc(k,iw) + wmarw(k) + wmanud(k,iw)
        enddo
     else
        do k = ka, mza-1
           wmasc(k,iw) = wmasc(k,iw) + wmarw(k)
        enddo
     endif

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

     rho(k,iw) = rho(k,iw) + delex_rho(k) + b10(k) * (wmarw(k-1) - wmarw(k))

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
  use mem_rayf

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


!=========================================================================


subroutine update_vxe_undrgnd(iw, dts, vmxe0, vmye0, vmze0, pg2, &
                              vmxet_rk, vmyet_rk, vmzet_rk)

  use mem_basic,  only: vxe, vye, vze, rho
  use mem_grid!,   only: mwa, mva, mza, lpw, lve2, lpv, nve2_max, volti
  use mem_ijtabs, only: itab_w
!  use vel_t3d,    only: exm, eym, ezm

  implicit none

  integer, intent(in) :: iw
  real,    intent(in) :: dts

  real,    intent(in) :: vmxe0(nve2_max,mwa)
  real,    intent(in) :: vmye0(nve2_max,mwa)
  real,    intent(in) :: vmze0(nve2_max,mwa)

  real,    intent(in) :: pg2(nve2_max,mva)

  real,    intent(in) :: vmxet_rk(mza,mwa)
  real,    intent(in) :: vmyet_rk(mza,mwa)
  real,    intent(in) :: vmzet_rk(mza,mwa)

  integer :: npoly, ka, jv, iv, kb, k, kv, ks, kt
  real    :: rsi
  real    :: pgfxe(nve2_max), pgfye(nve2_max),pgfze(nve2_max)

  if (lve2(iw) > 0) then

     npoly = itab_w(iw)%npoly
     ka    = lpw(iw)
     kt    = lpw(iw) + lve2(iw) - 1

     pgfxe(1:lve2(iw)) = 0.
     pgfye(1:lve2(iw)) = 0.
     pgfze(1:lve2(iw)) = 0.

!!     do jv = 1, npoly
!!        iv = itab_w(iw)%iv(jv)
!!        kb = lpv(iv)
!!
!!        do k = kb, kt
!!           kv = k - kb + 1
!!           ks = k - ka + 1
!!
!!           pgfxe(ks) = pgfxe(ks) + itab_w(iw)%ecvec_vx(jv) * pg2(kv,iv)
!!           pgfye(ks) = pgfye(ks) + itab_w(iw)%ecvec_vy(jv) * pg2(kv,iv)
!!           pgfze(ks) = pgfze(ks) + itab_w(iw)%ecvec_vz(jv) * pg2(kv,iv)
!!        enddo
!!     enddo

     do ks = 1, lve2(iw)
        k = ks + ka - 1

        rsi  = 1.0 / real(rho(k,iw))

        vxe(k,iw) = (vmxe0(ks,iw) + dts * ( vmxet_rk(k,iw) * volti(k,iw) &
                                          )) * rsi
!                                          + pgfxe(ks) * exm(ks,iw))) * rsi
        vye(k,iw) = (vmye0(ks,iw) + dts * ( vmyet_rk(k,iw) * volti(k,iw) &
                                          )) * rsi
!                                          + pgfye(ks) * eym(ks,iw))) * rsi
        vze(k,iw) = (vmze0(ks,iw) + dts * ( vmzet_rk(k,iw) * volti(k,iw) &
                                          )) * rsi
!                                          + pgfze(ks) * ezm(ks,iw))) * rsi
     enddo

  endif

end subroutine update_vxe_undrgnd




subroutine wind_vec_at_v(iv, vxe_v, vye_v, vze_v)

  use mem_grid,   only: mza, lpv, vnx, vny, vnz
  use mem_ijtabs, only: itab_v
  use mem_basic,  only: vc, vxe, vye, vze

  implicit none

  integer, intent( in) :: iv
  real,    intent(out) :: vxe_v(mza), vye_v(mza), vze_v(mza)
  integer              :: iw1, iw2, k
  real                 :: vx, vy, vz
  real                 :: vcc, vdiff

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  do k = lpv(iv), mza
     vx = 0.5 * (vxe(k,iw1) + vxe(k,iw2))
     vy = 0.5 * (vye(k,iw1) + vye(k,iw2))
     vz = 0.5 * (vze(k,iw1) + vze(k,iw2))

     vcc   = vnx(iv) * vx + vny(iv) * vy + vnz(iv) * vz
     vdiff = vc(k,iv) - vcc

     vxe_v(k) = vx + vdiff * vnx(iv)
     vye_v(k) = vy + vdiff * vny(iv)
     vze_v(k) = vz + vdiff * vnz(iv)
  enddo

end subroutine wind_vec_at_v


subroutine grad_t2d_v(iw, vxe_v, vye_v, vze_v, &
                          gxps_vxe, gyps_vxe, gxps_vye, gyps_vye, &
                          gxps_vze, gyps_vze)

  use mem_basic,  only: vxe, vye, vze
  use mem_grid,   only: mza, mva, lpv, gxps_coef, gyps_coef
  use mem_ijtabs, only: itab_w

  implicit none

  integer, intent( in) :: iw
  real,    intent( in) :: vxe_v(mza,mva), vye_v(mza,mva), vze_v(mza,mva)

  real,    intent(out) :: gxps_vxe(mza), gyps_vxe(mza)
  real,    intent(out) :: gxps_vye(mza), gyps_vye(mza)
  real,    intent(out) :: gxps_vze(mza), gyps_vze(mza)

  integer :: jv, iv, k
  real    :: dvx, dvy, dvz, gx, gy

  gxps_vxe = 0.
  gyps_vxe = 0.

  gxps_vye = 0.
  gyps_vye = 0.

  gxps_vze = 0.
  gyps_vze = 0.

! Loop over V faces of this T cell

  do jv = 1, itab_w(iw)%npoly
     iv = itab_w(iw)%iv(jv)

! Vertical loop over T levels
! Zero-gradient lateral B.C. below lpv(iv)

     gx = 2.0 * gxps_coef(iw,jv)
     gy = 2.0 * gyps_coef(iw,jv)

     do k = lpv(iv), mza
        dvx = vxe_v(k,iv) - vxe(k,iw)
        dvy = vye_v(k,iv) - vye(k,iw)
        dvz = vze_v(k,iv) - vze(k,iw)

        gxps_vxe(k) = gxps_vxe(k) + gx * dvx
        gyps_vxe(k) = gyps_vxe(k) + gy * dvx

        gxps_vye(k) = gxps_vye(k) + gx * dvy
        gyps_vye(k) = gyps_vye(k) + gy * dvy

        gxps_vze(k) = gxps_vze(k) + gx * dvz
        gyps_vze(k) = gyps_vze(k) + gy * dvz

!!        gxps_vxe(k) = gxps_vxe(k) + gxps_coefv(iw,jv) * dvx
!!        gyps_vxe(k) = gyps_vxe(k) + gyps_coefv(iw,jv) * dvx
!!
!!        gxps_vye(k) = gxps_vye(k) + gxps_coefv(iw,jv) * dvy
!!        gyps_vye(k) = gyps_vye(k) + gyps_coefv(iw,jv) * dvy
!!
!!        gxps_vze(k) = gxps_vze(k) + gxps_coefv(iw,jv) * dvz
!!        gyps_vze(k) = gyps_vze(k) + gyps_coefv(iw,jv) * dvz
     enddo

  enddo

end subroutine grad_t2d_v


end module wrtv_rk
