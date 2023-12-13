subroutine scalar_transport_rk(rho_old)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, &
                          jtv_wadj, jtw_prog, jtw_wadj
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, volti, &
                          dztsqo6, dzto2, arw, arv
  use misc_coms,    only: dtlm, iparallel
  use var_tables,   only: num_scalar, scalar_tab
  use mem_turb,     only: akhodx, khtopv
  use oname_coms,   only: nl
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use grad_lib,     only: grad_t2d, grad_t2d_quad, grad_z_quad
  use mem_adv,      only: gxps_scp, gyps_scp, gxyps_scp, gxxps_scp, gyyps_scp, &
                          xx0_v, yy0_v, xy0_v
  use mem_basic,    only: rho, vmsca=>vmsc, wmsca=>wmsc, thil, theta, rr_v
  use misc_coms,    only: nrk_scal
  use tridiag,      only: tridif_prep, tridif_fini
  use mem_nudge,    only: nudflag, rhot_nud
  use consts_coms,  only: r8

  use mem_para,     only: myrank
  use mem_micro,    only: rr_c, rr_r

  implicit none

  real(r8), intent(in) :: rho_old(mza,mwa)

  integer  :: j,iw,iw1,iw2,istage
  integer  :: n,k,kb,iv,iwn,jv
  integer  :: iorder,imonot
  real     :: dirv
  real     :: dtl
  real(r8) :: rhof(mza)

! Automatic arrays:

  real(r8) :: rscp0(mza,mwa)
  real(r8) :: rhoi (mza,mwa,nrk_scal)
  real(r8) :: fluxdiv(mza), vflux(mza)

  real :: scp_upv(mza,mva)
  real :: scp_upw(mza)
  real :: hflux(mza)
  real :: gzps(mza), gzzps(mza)

  real, allocatable :: scpmax(:,:)
  real, allocatable :: scpmin(:,:)

  logical :: do_implic(mwa)
  real    :: ff_implic(mwa)
  integer :: nn_implic, nn, ni

  real, pointer, contiguous :: scp(:,:)
  real, pointer, contiguous :: sct(:,:)

  real, parameter :: twothird = 2. / 3.
  real, parameter :: onethird = 1. / 3.

  real :: c1, c2, scp2
  real :: delex_scp(mza), del_scp(mza)

  integer, allocatable :: iw_implic(:), ii_implic(:)

  real :: wup(mza), wdn(mza)
  real, allocatable :: b1 (:,:), b2 (:,:), b3(:,:)

  real :: b2a(mza), b3a(mza)
  integer :: kwtop(mwa)

  real, parameter :: frk3(3) = (/ onethird, 0.5, 1.0 /)
  real, parameter :: frk2(2) =           (/ 0.5, 1.0 /)

  iorder = nl%horiz_adv_order
  imonot = nl%iscal_monot

  if (imonot > 0) then
     allocate(scpmax(mza,mwa))
     allocate(scpmin(mza,mwa))
  endif

  !$omp parallel
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
     do k = lpw(iw), mza-1
        wmsca(k,iw) = wmsca(k,iw) * arw(k,iw)
     enddo
  enddo
  !$omp end do nowait

  !$omp do private(iv,k)
  do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)
     do k = lpv(iv), mza
        vmsca(k,iv) = vmsca(k,iv) * arv(k,iv)
     enddo
  enddo
  !$omp end do nowait
  !$omp end parallel

  call check_cfls( vmsca, wmsca, rho_old, nrk_scal, ff_implic )

  do_implic = .false.
  nn_implic = 0

  dtl = dtlm

  !$omp parallel private(rhof)
  !$omp do private(iw,kb,k) reduction(+:nn_implic)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     kb = lpw(iw)

     if ( nudflag > 0 .and. nl%nud_preserve_mix_ratio .and. &
          itab_w(iw)%mrlw <= nl%max_nud_mrl ) then

        do k = kb, mza
           rhof(k) = rho(k,iw) - dtl * rhot_nud(k,iw)
        enddo

     else

        do k = kb, mza
           rhof(k) = rho(k,iw)
        enddo

     endif

     if (nrk_scal == 3) then

        do k = kb, mza
           rhoi(k,iw,1) = 3._r8 / (2._r8 * rho_old(k,iw) + rhof(k))
           rhoi(k,iw,2) = 2._r8 / (rho_old(k,iw) + rhof(k))
           rhoi(k,iw,3) = 1._r8 / rhof(k)
        enddo

     else

        do k = kb, mza
           rhoi(k,iw,1) = 2._r8 / (rho_old(k,iw) + rhof(k))
           rhoi(k,iw,2) = 1._r8 / rhof(k)
        enddo

     endif

     if (ff_implic(iw) > 0.025) then
        do_implic(iw) = .true.
        nn_implic = nn_implic + 1
     endif

     kwtop(iw) = maxval( khtopv( itab_w(iw)%iv( 1:itab_w(iw)%npoly ) ) )
  enddo
  !$omp end do nowait
  !$omp end parallel

  if (nn_implic > 0) then

     allocate(iw_implic(nn_implic))
     allocate(ii_implic(mwa))
     allocate(b1(mza,nn_implic))
     allocate(b2(mza,nn_implic))
     allocate(b3(mza,nn_implic))

     nn = 0
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        if (do_implic(iw)) then
           nn = nn + 1
           ii_implic(iw) = nn
           iw_implic(nn) = iw
        endif
     enddo

     !$omp parallel private(b2a,b3a,wup,wdn)
     wup(mza) = 0.0
     wdn(mza) = 0.0

     !$omp do private(iw,kb,k,c1,c2)
     do j = 1, nn_implic; iw = iw_implic(j)

        kb = lpw(iw)

        wup(kb-1) = 0.0
        wdn(kb-1) = 0.0

        do k = kb, mza-1
           wup(k) = max(wmsca(k,iw), 0.0)
           wdn(k) = min(wmsca(k,iw), 0.0)
        enddo

        c1 = ff_implic(iw) * dtl

        do k = kb, mza
           c2 = c1 * volti(k,iw) * real(rhoi(k,iw,nrk_scal))
           b1 (k,j) =     - c2 * wup(k-1)
           b2a(k)   = 1.0 + c2 * (wup(k) - wdn(k-1))
           b3a(k)   =       c2 * wdn(k)
        enddo

        call tridif_prep(mza,kb,mza,b1(:,j),b2a,b3a,b2(:,j),b3(:,j))
     enddo
     !$omp end do nowait
     !$omp end parallel

  endif

! LOOP OVER SCALARS HERE

  do n = 1, num_scalar

     ! Point SCP and SCT to scalar table arrays

     scp => scalar_tab(n)%var_p(:,:)
     sct => scalar_tab(n)%var_t(:,:)

     !$omp parallel private(hflux)
     !$omp do private(iw,kb,k,jv,iv,iwn)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        kb = lpw(iw)

        ! Store initial scalar mass before R-K loop

        do k = kb, mza
           rscp0(k,iw) = scp(k,iw) * rho_old(k,iw)
        enddo

        ! Compute horizontal flux divergence and add to long-timestep tendency

        if (scalar_tab(n)%do_sgsmix .and. kwtop(iw) >= kb) then

           hflux(kb:kwtop(iw)) = 0.

           do jv = 1, itab_w(iw)%npoly
              iv  = itab_w(iw)%iv(jv)
              iwn = itab_w(iw)%iw(jv)

              ! Vertical loop over T levels
              do k = lpv(iv), khtopv(iv)
                 hflux(k) = hflux(k) + akhodx(k,iv) * (scp(k,iwn) - scp(k,iw))
              enddo

           enddo

           ! Vertical loop over T levels
           do k = kb, kwtop(iw)
              sct(k,iw) = sct(k,iw) + hflux(k) * volti(k,iw)
           enddo

        endif

     enddo
     !$omp end do nowait
     !$omp end parallel

     ! Loop over R-K stages

     do istage = 1, nrk_scal

        if (nrk_scal == 3) then
           dtl = frk3(istage) * real( dtlm )
        else
           dtl = frk2(istage) * real( dtlm )
        endif

        !$omp parallel do private(iw,kb,k,jv,iv,iwn)
        do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

           if (iorder <= 2) then
              call grad_t2d(iw, scp, gxps_scp(:,iw), gyps_scp(:,iw))
           else
              call grad_t2d_quad(iw, scp, gxps_scp(:,iw), gyps_scp(:,iw), &
                                 gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))
           endif

           if (imonot > 0) then
              if (istage == nrk_scal) then

                 if (iorder <= 2) then
                    call limit_h2(iw, scp, gxps_scp(:,iw), gyps_scp(:,iw))
                 else
                    call limit_h3(iw, scp, gxps_scp(:,iw), gyps_scp(:,iw),&
                                  gxxps_scp(:,iw), gxyps_scp(:,iw), gyyps_scp(:,iw))
                 endif

              else ! istage /= nrk_scal

                 kb = lpw(iw)

                 scpmax(kb,iw) = max(scp(kb,iw), scp(kb+1,iw))
                 scpmin(kb,iw) = min(scp(kb,iw), scp(kb+1,iw))
                 do k = kb+1, mza-1
                    scpmax(k,iw) = max(scp(k-1,iw), scp(k,iw), scp(k+1,iw))
                    scpmin(k,iw) = min(scp(k-1,iw), scp(k,iw), scp(k+1,iw))
                 enddo
                 scpmax(mza,iw) = max(scp(mza-1,iw), scp(mza,iw))
                 scpmin(mza,iw) = min(scp(mza-1,iw), scp(mza,iw))

                 ! Loop over V neighbors of this W cell
                 do jv = 1, itab_w(iw)%npoly
                    iv  = itab_w(iw)%iv(jv)
                    iwn = itab_w(iw)%iw(jv)
                    do k = lpv(iv), mza
                       scpmax(k,iw) = max(scpmax(k,iw), scp(k,iwn))
                       scpmin(k,iw) = min(scpmin(k,iw), scp(k,iwn))
                    enddo
                 enddo

              endif
           endif

        enddo
        !$omp end parallel do

        ! MPI send/recv of SCP horizontal gradient components

        if (iparallel == 1) then

           if (iorder <= 2) then
              call mpi_send_w(rvara1=gxps_scp, rvara2=gyps_scp)
              call mpi_recv_w(rvara1=gxps_scp, rvara2=gyps_scp)
           else
              call mpi_send_w(rvara1=gxps_scp, rvara2=gyps_scp, &
                              rvara3=gxxps_scp, rvara4=gxyps_scp, rvara5=gyyps_scp)
              call mpi_recv_w(rvara1=gxps_scp, rvara2=gyps_scp, &
                              rvara3=gxxps_scp, rvara4=gxyps_scp, rvara5=gyyps_scp)
           endif

        endif

        if (iorder <= 2) then
           call lbcopy_w(a1=gxps_scp, a2=gyps_scp)
        else
           call lbcopy_w(a1=gxps_scp, a2=gyps_scp, &
                         a3=gxxps_scp, a4=gxyps_scp, a5=gyyps_scp)
        endif

        ! Horizontal loop over all V points surrounding primary W/T columns
        ! to compute upwinded scalar value at the V points

        !$omp parallel private(delex_scp,del_scp,fluxdiv,vflux,scp_upw,gzps,gzzps)
        !$omp do private(iv,k,iw1,iw2)
        do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)
           iw1 = itab_v(iv)%iw(1)
           iw2 = itab_v(iv)%iw(2)

           if (iorder <= 2) then

              do k = lpv(iv), mza

                 if (vmsca(k,iv) > 0.) then
                    scp_upv(k,iv) = vmsca(k,iv) * (           scp(k,iw1) &
                                  + itab_v(iv)%dxps(1) * gxps_scp(k,iw1) &
                                  + itab_v(iv)%dyps(1) * gyps_scp(k,iw1) )
                 else
                    scp_upv(k,iv) = vmsca(k,iv) * (           scp(k,iw2) &
                                  + itab_v(iv)%dxps(2) * gxps_scp(k,iw2) &
                                  + itab_v(iv)%dyps(2) * gyps_scp(k,iw2) )
                 endif

              enddo

           else

              do k = lpv(iv), mza

                 if (vmsca(k,iv) > 0.) then
                    scp_upv(k,iv) = vmsca(k,iv) * (            scp(k,iw1) &
                                  + itab_v(iv)%dxps(1) *  gxps_scp(k,iw1) &
                                  + itab_v(iv)%dyps(1) *  gyps_scp(k,iw1) &
                                  + xx0_v(1,iv)        * gxxps_scp(k,iw1) &
                                  + yy0_v(1,iv)        * gyyps_scp(k,iw1) &
                                  + xy0_v(1,iv)        * gxyps_scp(k,iw1) )

                 else
                    scp_upv(k,iv) = vmsca(k,iv) * (            scp(k,iw2) &
                                  + itab_v(iv)%dxps(2) *  gxps_scp(k,iw2) &
                                  + itab_v(iv)%dyps(2) *  gyps_scp(k,iw2) &
                                  + xx0_v(2,iv)        * gxxps_scp(k,iw2) &
                                  + yy0_v(2,iv)        * gyyps_scp(k,iw2) &
                                  + xy0_v(2,iv)        * gxyps_scp(k,iw2) )
                 endif

              enddo
           endif

        enddo
        !$omp end do

        !$omp do private (iw,kb,jv,iv,dirv,k,ni,iwn,scp2)
        do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

           kb = lpw(iw)

           ! Vertical scalar gradients

           call grad_z_quad(iw, scp(:,iw), gzps, gzzps)

           ! Limit gradients for monotonic advection

           if (imonot > 0 .and. istage == nrk_scal) then
              call limit_z(iw, scp(:,iw), gzps, gzzps)
           endif

           ! Compute scalar values at top and bottom faces

           do k = lpw(iw), mza-1
              if (wmsca(k,iw) > 0.) then
                 scp_upw(k) = wmsca(k,iw) * (scp(k,iw) &
                            +   dzto2(k) *  gzps(k)    &
                            + dztsqo6(k) * gzzps(k))
              else
                 scp_upw(k) = wmsca(k,iw) *   (scp(k+1,iw) &
                            -   dzto2(k+1) *  gzps(k+1)    &
                            + dztsqo6(k+1) * gzzps(k+1))
              endif
           enddo

           vflux(kb-1) = 0._r8
           vflux(mza ) = 0._r8
           do k = kb, mza-1
              vflux(k) = real( scp_upw(k), r8)
           enddo

           do k = kb, mza
              fluxdiv(k) = vflux(k-1) - vflux(k)
           enddo

           ! Loop over neighbor V points of this W cell
           do jv = 1, itab_w(iw)%npoly
              iv   = itab_w(iw)%iv  (jv)
              dirv = itab_w(iw)%dirv(jv)

              ! Horizontal advective fluxes
              do k = lpv(iv), mza
                 fluxdiv(k) = fluxdiv(k) + dirv * scp_upv(k,iv)
              enddo
           enddo

           if (do_implic(iw) .and. istage == nrk_scal) then

              do k = kb, mza
                 scp2 = (rscp0(k,iw) + dtl * (sct(k,iw) + volti(k,iw) * fluxdiv(k))) * rhoi(k,iw,istage)
                 delex_scp(k) = scp2 - scp(k,iw)
              enddo

              ni = ii_implic(iw)
              call tridif_fini(mza,kb,mza,b1(:,ni),b2(:,ni),b3(:,ni),delex_scp,del_scp)

              do k = kb, mza
                 scp(k,iw) = scp(k,iw) + del_scp(k)
              enddo

           elseif (istage < nrk_scal .and. imonot > 0) then

              do k = kb, mza
                 scp2      = (rscp0(k,iw) + dtl * volti(k,iw) * fluxdiv(k)) * rhoi(k,iw,istage)
                 scp(k,iw) = min( max(scp2, scpmin(k,iw)), scpmax(k,iw) ) + dtl * sct(k,iw) * rhoi(k,iw,istage)
              enddo

           else

              do k = kb, mza
                 scp(k,iw) = (rscp0(k,iw) + dtl * (sct(k,iw) + volti(k,iw) * fluxdiv(k))) * rhoi(k,iw,istage)
              enddo

           endif

           if (nl%zero_neg_scalars .and. scalar_tab(n)%pdef) then
              do k = kb, mza
                 scp(k,iw) = max(scp(k,iw), 0.0)
              enddo
           endif

           scp(1:kb-1,iw) = scp(kb,iw)

        enddo
        !$omp end do nowait
        !$omp end parallel

        if (istage < nrk_scal) then
           if (iparallel == 1) then
              call mpi_send_w(rvara1=scp)
              call mpi_recv_w(rvara1=scp)
           endif
           call lbcopy_w(a1=scp)
        endif

     enddo ! istage

  enddo ! n

end subroutine scalar_transport_rk

!===========================================================================

subroutine limit_z(iw,sc,gzps,gzzps)

  use mem_grid, only: mza, lpw, dzto2, dztsqo6, dztsqo12

  implicit none

  integer, intent(in)    :: iw
  real,    intent(inout) :: gzps(mza), gzzps(mza)
  real,    intent(in)    :: sc(mza)

  real    :: scmax(mza), scmin(mza)
  real    :: gx, gxx, g1, g2, gmax, gmin, fact1, fact2, fact
! real    :: zz, g
  integer :: k, ka

  ka = lpw(iw)

  scmax(ka) = max(sc(ka),sc(ka+1)) - sc(ka)
  scmin(ka) = min(sc(ka),sc(ka+1)) - sc(ka)
  do k = ka+1, mza-1
     scmax(k) = max(sc(k-1),sc(k),sc(k+1)) - sc(k)
     scmin(k) = min(sc(k-1),sc(k),sc(k+1)) - sc(k)
  enddo
  scmax(mza) = max(sc(mza-1),sc(mza)) - sc(mza)
  scmin(mza) = min(sc(mza-1),sc(mza)) - sc(mza)

  do k = ka, mza
     gx  =   dzto2(k) *  gzps(k)
     gxx = dztsqo6(k) * gzzps(k)

     g1 = gxx - gx
     g2 = gxx + gx

     gmax = max(g1,g2)
     gmin = min(g1,g2)

!!     if (abs(gzzps(k)) > 1.e-25) then
!!        zz = -0.5 * gzps(k) / gzzps(k)
!!
!!        if (abs(zz) < dzto2(k)) then
!!           g = zz * gzps(k) + (zz * zz - dztsqo12(k)) * gzzps(k)
!!           gmax = max(gmax, g)
!!           gmin = min(gmin, g)
!!        endif
!!     endif

     fact1 = 1.0
     if (gmax > scmax(k)) fact1 = scmax(k) / gmax

     fact2 = 1.0
     if (gmin < scmin(k)) fact2 = scmin(k) / gmin

!    fact = 0.9999998 * min(fact1,fact2)
!    fact = 0.999999 * min(fact1,fact2)
!    fact = 0.9999 * min(fact1,fact2)
     fact = min(fact1,fact2)

     gzps (k) = fact * gzps (k)
     gzzps(k) = fact * gzzps(k)
  enddo

end subroutine limit_z




subroutine limit_h2(iw,scp,gxps,gyps)

  use mem_grid,   only: mza, lpw, mwa
  use mem_adv,    only: dxpsw_v, dypsw_v
  use mem_ijtabs, only: itab_w

  implicit none

  real,    intent(inout) :: gxps(mza), gyps(mza)
  real,    intent(in)    :: scp(mza,mwa)
  integer, intent(in)    :: iw

  real    :: gmax(mza), gmin(mza), g, fact1, fact2, fact
  real    :: scpmax(mza), scpmin(mza), scmax, scmin
  integer :: k, n, iv, iwn

  gmax = 0.0
  gmin = 0.0

  scpmax = scp(:,iw)
  scpmin = scp(:,iw)

!!  gxmax = 0.0
!!  gxmin = 0.0
!!
!!  gymax = 0.0
!!  gymin = 0.0
!!
!!  do n = 1, itab_w(iw)%npoly
!!     do k = lpw(iw), mza
!!
!!        gx = dxpsw_v(n,iw) * gxps(k)
!!        gy = dxpsw_v(n,iw) * gxps(k)
!!
!!        gxmax(k) = max(gxmax(k),gx)
!!        gxmin(k) = min(gxmin(k),gx)
!!
!!        gymax(k) = max(gymax(k),gy)
!!        gymin(k) = min(gymin(k),gy)
!!
!!        endif
!!
!!     enddo
!!  enddo
!!
!!  do k = lpw(iw), mza
!!     fact1 = 1.0
!!     if (gxmax(k) > scmax(k)) fact1 = scmax(k) / gxmax(k)
!!
!!     fact2 = 1.0
!!     if (gxmin(k) < scmin(k)) fact2 = scmin(k) / gxmin(k)
!!
!!     fact = 0.9999998 * min(fact1,fact2)
!!
!!     gxps(k) = fact * gxps(k)
!!
!!     fact1 = 1.0
!!     if (gymax(k) > scmax(k)) fact1 = scmax(k) / gymax(k)
!!
!!     fact2 = 1.0
!!     if (gymin(k) < scmin(k)) fact2 = scmin(k) / gymin(k)
!!
!!     fact = 0.9999998 * min(fact1,fact2)
!!
!!     gyps(k) = fact * gyps(k)
!!  enddo

  ! now combined!!!

  do n = 1, itab_w(iw)%npoly
     iv   = itab_w(iw)%iv(n)
     iwn  = itab_w(iw)%iw(n)

     do k = lpw(iw), mza
        g = dxpsw_v(n,iw) * gxps(k) + dypsw_v(n,iw) * gyps(k)
        gmax(k) = max(gmax(k),g)
        gmin(k) = min(gmin(k),g)

        scpmax(k) = max(scpmax(k),scp(k,iwn))
        scpmin(k) = min(scpmin(k),scp(k,iwn))
     enddo
  enddo

  do k = lpw(iw), mza
     scmax = scpmax(k) - scp(k,iw)
     scmin = scpmin(k) - scp(k,iw)

     fact1 = 1.0
     if (gmax(k) > scmax) fact1 = scmax / gmax(k)

     fact2 = 1.0
     if (gmin(k) < scmin) fact2 = scmin / gmin(k)

!    fact = 0.9999998 * min(fact1,fact2)
!    fact = 0.999999 * min(fact1,fact2)
!    fact = 0.9999 * min(fact1,fact2)
     fact = min(fact1,fact2)

     gxps(k) = fact * gxps(k)
     gyps(k) = fact * gyps(k)
  enddo

end subroutine limit_h2





subroutine limit_h3(iw,scp,gxps,gyps,gxxps,gxyps,gyyps)

  use mem_grid,   only: mza, lpw, mwa
  use mem_adv,    only: dxpsw_v, dypsw_v, dxxpsw_v, dxypsw_v, dyypsw_v
!                       xx0, yy0, xy0, dssq
  use mem_ijtabs, only: itab_w

  implicit none

  integer, intent(in)    :: iw
  real,    intent(in)    :: scp(mza,mwa)
  real,    intent(inout) :: gxps(mza), gyps(mza)
  real,    intent(inout) :: gxxps(mza), gxyps(mza), gyyps(mza)

  real    :: scpmax(mza), scpmin(mza), scmax, scmin
  real    :: gmax(mza), gmin(mza), g, fact1, fact2, fact
  integer :: k, n, iv, iwn

!  real :: dxpsm1, dypsm1, dxxpsm1, dxypsm1, dyypsm1, g1
!  real :: dxpsm2, dypsm2, dxxpsm2, dxypsm2, dyypsm2, g2
!  real :: denom, xx, yy

  gmax = 0.0
  gmin = 0.0

  scpmax = scp(:,iw)
  scpmin = scp(:,iw)

  do n = 1, itab_w(iw)%npoly
     iv   = itab_w(iw)%iv(n)
     iwn  = itab_w(iw)%iw(n)

     do k = lpw(iw), mza
        g =  dxpsw_v(n,iw) *  gxps(k) +  dypsw_v(n,iw) * gyps(k) &
          + dxxpsw_v(n,iw) * gxxps(k) + dxypsw_v(n,iw) * gxyps(k) + dyypsw_v(n,iw) * gyyps(k)

        gmax(k) = max(gmax(k),g)
        gmin(k) = min(gmin(k),g)

        scpmax(k) = max(scpmax(k),scp(k,iwn))
        scpmin(k) = min(scpmin(k),scp(k,iwn))
     enddo
  enddo

  do k = lpw(iw), mza
     scmax = scpmax(k) - scp(k,iw)
     scmin = scpmin(k) - scp(k,iw)

!!     denom = 4. * gxxps(k) * gyyps(k) - gxyps(k)**2
!!
!!     if (abs(denom) > 1.e-25) then
!!        yy = (gxps(k) * gxyps(k) - 2. * gxxps(k) * gyps(k)) / denom
!!        xx = (gyps(k) * gxyps(k) - 2. * gyyps(k) * gxps(k)) / denom
!!
!!        if (xx * xx + yy * yy < dssq(iw)) then
!!
!!           g = xx * gxps(k) + yy * gyps(k)  &
!!             + (xx * xx - xx0(iw)) * gxxps(k) + (yy * yy - yy0(iw)) * gyyps(k) &
!!             + (xx * yy - xy0(iw)) * gxyps(k)
!!
!!           gmax(k) = max(gmax(k),g)
!!           gmin(k) = min(gmin(k),g)
!!
!!        endif
!!     endif

     fact1 = 1.0
     if (gmax(k) > scmax) fact1 = scmax / gmax(k)

     fact2 = 1.0
     if (gmin(k) < scmin) fact2 = scmin / gmin(k)

!    fact = 0.9999998 * min(fact1,fact2)
!    fact = 0.999998 * min(fact1,fact2)
!    fact = 0.9999 * min(fact1,fact2)
     fact = min(fact1,fact2)

     gxps (k) = fact *  gxps(k)
     gyps (k) = fact *  gyps(k)
     gxxps(k) = fact * gxxps(k)
     gxyps(k) = fact * gxyps(k)
     gyyps(k) = fact * gyyps(k)
  enddo

end subroutine limit_h3

!===============================================================================

  subroutine check_cfls(vmsca, wmsca, rho, nrk, ff_implic)

    ! Diagnose outflow CFL number; also used for advection CFL stability check

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, mza, mva, mwa
    use consts_coms, only: r8
    use misc_coms,   only: iparallel, dtlm
    use mem_para,    only: myrank, mgroupsize
    use mem_basic,   only: wc
!$  use omp_lib,     only: omp_get_max_threads, omp_get_thread_num

#ifdef OLAM_MPI
    use mpi
#endif

    implicit none

    integer,           intent(in ) :: nrk
    real,              intent(in ) :: vmsca(mza,mva)
    real,              intent(in ) :: wmsca(mza,mwa)
    real(r8),          intent(in ) :: rho  (mza,mwa)
    real,              intent(out) :: ff_implic(mwa)

    integer :: j, iv, k, iw, n
    integer :: ier, inode, iwnode, ihnode
    integer :: mlocc, mlocw, mloch, nthreads, myid
    real    :: dt, fimp
    real    :: cfl(mza), cflh(mza), dovr, cflw, cfl0
    real    :: ct, ch, wm
    integer :: itk, itw, ihk, ihw, iwk, iww

    integer, allocatable :: imax(:,:), iwmax(:,:), ihmax(:,:)
    real,    allocatable :: cfl_max(:), w_max(:), cflh_max(:)

#ifdef OLAM_MPI
    real :: sbuf(9), rbuf(9,mgroupsize)
#endif

    myid     = 1
    nthreads = 1
 !$ nthreads = omp_get_max_threads()

    allocate(cfl_max(nthreads))
    allocate(imax (2,nthreads))

    allocate(cflh_max(nthreads))
    allocate(ihmax (2,nthreads))

    allocate(w_max  (nthreads))
    allocate(iwmax(2,nthreads))

    if (nrk == 3) then
       cfl0 = 1.50
    else
       cfl0 = 1.25
    endif

    !$omp parallel private(myid,cfl,cflh,ct,itk,itw,ch,ihk,ihw,wm,iwk,iww)
    !$ myid = omp_get_thread_num() + 1

    ct  = 0.0
    itk = 1
    itw = 1

    ch  = 0.0
    ihk = 1
    ihw = 1

    wm  = 0.0
    iwk = 1
    iww = 1

    !$omp do private(iw,dt,k,n,iv,dovr,cflw,fimp)
    do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

       dt = dtlm

       cflh = 0.
       fimp = 0.

       do n = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(n)
          do k = lpv(iv), mza
             cflh(k) = cflh(k) - min(itab_w(iw)%dirv(n) * vmsca(k,iv), 0.0)
          enddo
       enddo

       do k = lpw(iw), mza
          dovr = dt * volti(k,iw) / real(rho(k,iw))
          cflh(k) = cflh(k) * dovr
          cflw    = (max(wmsca(k,iw),0.0) - min(wmsca(k-1,iw),0.0)) * dovr
          cfl (k) = cflh(k) + cflw
          fimp    = max(fimp, min(1.0, max(0.0, cflw - min(0.5, cfl0 - cfl(k)))))
       enddo

       ff_implic(iw) = fimp

       do k = lpw(iw), mza
          if (cfl(k) > ct .or. cfl(k) /= cfl(k)) then
             ct  = cfl(k)
             itk = k
             itw = itab_w(iw)%iwglobe
          endif

          if (cflh(k) > ch .or. cflh(k) /= cflh(k)) then
             ch  = cflh(k)
             ihk = k
             ihw = itab_w(iw)%iwglobe
          endif

          if (abs(wc(k,iw)) > abs(wm) .or. wc(k,iw) /= wc(k,iw)) then
             wm  = wc(k,iw)
             iwk = k
             iww = itab_w(iw)%iwglobe
          endif
       enddo

    enddo
    !$omp end do nowait

    cfl_max(myid) = ct
    imax (:,myid) = (/ itk, itw /)

    cflh_max(myid) = ch
    ihmax (:,myid) = (/ ihk, ihw /)

    w_max  (myid) = wm
    iwmax(:,myid) = (/ iwk, iww /)
    !$omp end parallel

    mlocc = 1
    mloch = 1
    mlocw = 1

    !$ if (nthreads > 1) then
    !$
    !$    mlocw = maxloc(abs(w_max), dim=1)
    !$    mlocc = maxloc(  cfl_max , dim=1)
    !$    mloch = maxloc( cflh_max , dim=1)
    !$
    !$    do n = 1, nthreads
    !$       if ( cfl_max(n) /=  cfl_max(n)) mlocc = n
    !$       if (   w_max(n) /=    w_max(n)) mlocw = n
    !$       if (cflh_max(n) /= cflh_max(n)) mloch = n
    !$    enddo
    !$
    !$ endif

#ifdef OLAM_MPI
    if (iparallel == 1) then
       sbuf(1)   = cfl_max(mlocc)
       sbuf(2:3) = real(imax(:,mlocc))

       sbuf(4)   = w_max(mlocw)
       sbuf(5:6) = real(iwmax(:,mlocw))

       sbuf(7) = cflh_max(mloch)
       sbuf(8:9) = real(ihmax(:,mloch))

       call MPI_Gather(sbuf, 9, MPI_REAL, rbuf, 9, MPI_REAL, &
                       0, MPI_COMM_WORLD, ier)
    endif
#endif

    if (myrank == 0) then
       inode  = 1
       ihnode = 1
       iwnode = 1

#ifdef OLAM_MPI
       if (iparallel == 1) then

          if (any( rbuf(1,:) /= rbuf(1,:) )) then
             do n = 1, mgroupsize
                if (rbuf(1,n) /= rbuf(1,n)) then
                   inode          = n
                   cfl_max(mlocc) = rbuf(1,inode)
                   imax (:,mlocc) = nint(rbuf(2:3,inode))
                   exit
                endif
             enddo
          else
             inode          = maxloc(rbuf(1,:), dim=1)
             cfl_max(mlocc) = rbuf(1,inode)
             imax (:,mlocc) = nint(rbuf(2:3,inode))
          endif

          if (any( rbuf(4,:) /= rbuf(4,:) )) then
             do n = 1, mgroupsize
                if (rbuf(4,n) /= rbuf(4,n)) then
                   iwnode         = n
                   w_max  (mlocw) = rbuf(4,iwnode)
                   iwmax(:,mlocw) = nint(rbuf(5:6,iwnode))
                   exit
                endif
             enddo
          else
             iwnode         = maxloc(abs(rbuf(4,:)), dim=1)
             w_max  (mlocw) = rbuf(4,iwnode)
             iwmax(:,mlocw) = nint(rbuf(5:6,iwnode))
          endif

          if (any( rbuf(7,:) /= rbuf(7,:) )) then
             do n = 1, mgroupsize
                if (rbuf(7,n) /= rbuf(7,n)) then
                   ihnode         = n
                   cflh_max(mloch) = rbuf(7,ihnode)
                   ihmax (:,mloch) = nint(rbuf(8:9,ihnode))
                   exit
                endif
             enddo
          else
             ihnode          = maxloc(rbuf(7,:), dim=1)
             cflh_max(mloch) = rbuf(7,ihnode)
             ihmax (:,mloch) = nint(rbuf(8:9,ihnode))
          endif

       endif
#endif

       write(*,'(5x,A,f0.3,3(A,I0))') "Max total CFL = ", cfl_max(mlocc),  &
            " at node ", inode-1, ", iwglobe=", imax(2,mlocc), ", k=", imax(1,mlocc)

       write(*,'(5x,A,f0.3,3(A,I0))') "Max horiz CFL = ", cflh_max(mloch),  &
            " at node ", ihnode-1, ", iwglobe=", ihmax(2,mloch), ", k=", ihmax(1,mloch)

       write(*,'(5x,A,f0.3,3(A,I0))') "Max  W (m/s)  = ", w_max(mlocw),  &
            " at node ", iwnode-1, ", iwglobe=", iwmax(2,mlocw), ", k=", iwmax(1,mlocw)
    endif

  end subroutine check_cfls
