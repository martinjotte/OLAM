subroutine scalar_transport_orig(rho_old)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, &
                          jtv_wadj, jtw_prog, jtw_wadj
  use mem_grid,     only: mza, mva, mwa, lpv, lpw, arv, arw, volti, &
                          zwgt_bot8, zwgt_top8
  use misc_coms,    only: dtlm, iparallel
  use var_tables,   only: num_scalar, scalar_tab
  use mem_turb,     only: akhodx, khtopv
  use oname_coms,   only: nl
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use obnd,         only: lbcopy_w
  use consts_coms,  only: r8
  use grad_lib,     only: grad_t2d, grad_z
  use mem_basic,    only: vmsc, wmsc, vxesc, vyesc, vzesc, rho
  use mem_nudge,    only: nudflag, rhot_nud

  implicit none

  real(r8), intent(in) :: rho_old(mza,mwa)

  integer  :: j,iw,iw1,iw2,iwd
  integer  :: n,k,kb,kd,iv,iwn,jv
  real     :: dirv

! Automatic arrays:

  real :: vsc(mza)
  real :: wsc(mza)

  real :: vmsca(mza,mva)
  real :: wmsca(mza,mwa)
  real :: scp_prev(mza,mwa)

!  integer :: iwdepv(mza,mva)

!  integer :: kdepw(mza,mwa)

  real :: scp_upv(mza,mva)
  real :: scp_upw(mza,mwa)

  real(r8) :: hflux(mza), vflux(mza)
  real(r8) :: scp_tend(mza)
  real(r8) :: rhoi(mza,mwa)
  real(r8) :: dt
  real     :: dt4

  real :: gxps_scp(mza,mwa), gyps_scp(mza,mwa), gzps_scp(mza,mwa)
  real :: dxps_w(mza,mwa), dyps_w(mza,mwa), dzps_w(mza,mwa)
  real :: dxps_v(mza,mva), dyps_v(mza,mva), dzps_v(mza,mva)

  real, pointer, contiguous :: scp(:,:)
  real, pointer, contiguous :: sct(:,:)

! Horizontal loop over all primary W columns to diagnose
! face-normal vertical velocity at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp parallel private(wsc,vsc)
  !$omp do private(iw,kb,dt4,k)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     kb  = lpw(iw)
     dt4 = dtlm

     do k = kb, mza-1
        wsc  (k)    = wmsc(k,iw) &
                    / real( zwgt_bot8(k) * rho_old(k,iw) + zwgt_top8(k) * rho_old(k+1,iw) )
        wmsca(k,iw) = wmsc(k,iw) * arw(k,iw)
     enddo

     wsc    (kb-1)    = 0.0
     wmsca  (kb-1,iw) = 0.0
     scp_upw(kb-1,iw) = 0.0

     wsc    (mza)    = 0.0
     wmsca  (mza,iw) = 0.0
     scp_upw(mza,iw) = 0.0

     call donorpointw(iw, dt4, wsc, vxesc(:,iw), vyesc(:,iw), vzesc(:,iw), &
                      dxps_w(:,iw), dyps_w(:,iw), dzps_w(:,iw))

     if ( nudflag > 0 .and. nl%nud_preserve_mix_ratio .and. &
          itab_w(iw)%mrlw <= nl%max_nud_mrl ) then

        do k = kb, mza
           rhoi(k,iw) = 1._r8 / ( rho(k,iw) - real( dt4 * rhot_nud(k,iw), r8 ) )
        enddo

     else

        do k = kb, mza
           rhoi(k,iw) = 1._r8 / rho(k,iw)
        enddo

     endif

  enddo
  !$omp end do nowait

! Loop over V columns (the immediate neighbors of all primary W points) to
! diagnose face-normal velocity components at (t + 1/2) from mass fluxes
! (OK to use density at time t)

  !$omp do private(iv,iw1,iw2,kb,dt4,k)
  do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)
     kb  = lpv(iv)
     dt4 = dtlm

     do k = kb, mza
        vsc  (k)    = 2.0 * vmsc(k,iv) / real( rho_old(k,iw1) + rho_old(k,iw2) )
        vmsca(k,iv) = vmsc(k,iv) * arv(k,iv)
     enddo

     call donorpointv(iv, dt4, vsc, vxesc, vyesc, vzesc, &
                      dxps_v(:,iv), dyps_v(:,iv), dzps_v(:,iv))

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! Compute outflow CFL numbers for long timestep stability check;
  ! also needed by Thuburn monotonic scheme

  call comp_cfls_long( vmsca, wmsca, rho_old )

! LOOP OVER SCALARS HERE
  do n = 1, num_scalar

     ! Point SCP and SCT to scalar table arrays
     scp => scalar_tab(n)%var_p(:,:)
     sct => scalar_tab(n)%var_t(:,:)

     ! Save previous SCP for diffusion
     scp_prev = scp

! Evaluate horizontal gradients of scalar field on a local polar
! stereographic projection

     !$omp parallel do private(iw)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        call grad_t2d(iw, scp, gxps_scp(:,iw), gyps_scp(:,iw))
        call grad_z  (iw, scp(:,iw), gzps_scp(:,iw))

     enddo
     !$omp end parallel do

! MPI send of SCP horizontal gradient components

     if (iparallel == 1) then
        call mpi_send_w(rvara1=gxps_scp, rvara2=gyps_scp, rvara3=gzps_scp)
     endif

! Horizontal loop over all primary W columns to compute the
! upwinded scalar value at each W interface. This can be
! computed before scalar gradients are received at the borders

     !$omp parallel do private(iw,k,kd)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        do k = lpw(iw), mza-1

           if (wmsc(k,iw) > 0.) then
              kd = k
           else
              kd = k+1
           endif

           scp_upw(k,iw) =  wmsca (k,iw) * (     scp(kd,iw) &
                         +  dxps_w(k,iw) *  gxps_scp(kd,iw) &
                         +  dyps_w(k,iw) *  gyps_scp(kd,iw) &
                         +  dzps_w(k,iw) *  gzps_scp(kd,iw) )
        enddo
     enddo
     !$omp end parallel do

! MPI recv of SCP horizontal gradient components

     if (iparallel == 1) then
        call mpi_recv_w(rvara1=gxps_scp, rvara2=gyps_scp, rvara3=gzps_scp)
     endif

     call lbcopy_w(a1=gxps_scp, a2=gyps_scp, a3=gzps_scp)

! Horizontal loop over all V points surrounding primary W/T columns
! to compute upwinded scalar value at the V points

     !$omp parallel do private(iv,iw1,iw2,k,iwd)
     do j = 1,jtab_v(jtv_wadj)%jend; iv = jtab_v(jtv_wadj)%iv(j)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        do k = lpv(iv), mza

           if (vmsca(k,iv) > 0) then
              iwd = iw1
           else
              iwd = iw2
           endif

           scp_upv(k,iv) =  vmsca(k,iv) * (    scp(k,iwd) &
                         + dxps_v(k,iv) * gxps_scp(k,iwd) &
                         + dyps_v(k,iv) * gyps_scp(k,iwd) &
                         + dzps_v(k,iv) * gzps_scp(k,iwd) )
        enddo

     enddo
     !$omp end parallel do

! Horizontal loop over W/T points

     !$omp parallel private(hflux,vflux,scp_tend)
     !$omp do private (iw,kb,dt,jv,iv,iwn,dirv,k)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        kb = lpw(iw)
        dt = dtlm

        vflux(kb-1)     = 0.0
        vflux(kb:mza-1) = real( scp_upw(kb:mza,iw), r8)
        vflux(mza)      = 0.0

! Loop over neighbor V points of this W cell

        do k = kb, mza
           hflux(k) = vflux(k-1) - vflux(k)
        enddo

        ! Loop over V faces
        do jv = 1, itab_w(iw)%npoly

           iv   = itab_w(iw)%iv(jv)
           iwn  = itab_w(iw)%iw(jv)
           dirv = itab_w(iw)%dirv(jv)

           ! Vertical loop over T levels to compute horizontal advective fluxes
           do k = lpv(iv), mza
              hflux(k) = hflux(k) + dirv * scp_upv(k,iv)
           enddo

           ! Vertical loop over T levels to compute horizontal difusive fluxes
           do k = lpv(iv), khtopv(iv)
              hflux(k) = hflux(k) + akhodx(k,iv) * (scp_prev(k,iwn) - scp_prev(k,iw))
           enddo

        enddo

        ! Vertical loop over T levels to include tendencies from horizontal
        ! diffusion and horiztonal and vertical advection
        do k = kb, mza
           scp_tend(k) = sct(k,iw) + real(volti(k,iw),r8) * hflux(k)
        enddo

        if (nl%zero_neg_scalars .and. scalar_tab(n)%pdef) then

           ! Vertical loop over T levels to update scalar mixing ratios,
           ! resetting negative values to zero
           do k = kb, mza
              scp(k,iw) = max( (scp(k,iw) * rho_old(k,iw) + dt * scp_tend(k)) * rhoi(k,iw), 0._r8)
           enddo

        else

           ! Vertical loop over T levels to update scalar mixing ratios
           do k = kb, mza
              scp(k,iw) = (scp(k,iw) * rho_old(k,iw) + dt * scp_tend(k)) * rhoi(k,iw)
           enddo

        endif

     enddo
     !$omp end do nowait
     !$omp end parallel

  enddo ! n

end subroutine scalar_transport_orig


subroutine donorpointv(iv, dt, vs, vxe, vye, vze, dxps, dyps, dzps)

  use mem_ijtabs,  only: itab_v, jtv_wadj
  use mem_grid,    only: mza, mwa, lpv, unx, uny, unz, xev, yev, zev, dzto4
  use misc_coms,   only: mdomain
  use consts_coms, only: eradi

  implicit none

  integer, intent(in)  :: iv
  real,    intent(in)  :: dt
  real,    intent(in)  :: vs(mza)
  real,    intent(in)  :: vxe(mza,mwa), vye(mza,mwa), vze(mza,mwa)
  real,    intent(out) :: dxps(mza), dyps(mza), dzps(mza)

  real :: dto2
  real :: dxps1, dyps1, dxps2, dyps2
  real :: cosv1, sinv1, cosv2, sinv2
  real :: wnx_v, wny_v, wnz_v
  real :: vxeface, vyeface, vzeface

  real :: ufacev(mza), wfacev

  integer :: kb, k, iw1, iw2

  iw1 = itab_v(iv)%iw(1)
  iw2 = itab_v(iv)%iw(2)

  kb = lpv(iv)

  dto2 = .5 * dt

  dxps1 = itab_v(iv)%dxps(1)
  dyps1 = itab_v(iv)%dyps(1)

  dxps2 = itab_v(iv)%dxps(2)
  dyps2 = itab_v(iv)%dyps(2)

  cosv1 = itab_v(iv)%cosv(1)
  sinv1 = itab_v(iv)%sinv(1)

  cosv2 = itab_v(iv)%cosv(2)
  sinv2 = itab_v(iv)%sinv(2)

  if (mdomain <= 1) then
     wnx_v = xev(iv) * eradi
     wny_v = yev(iv) * eradi
     wnz_v = zev(iv) * eradi
  else
     wnx_v = 0.
     wny_v = 0.
     wnz_v = 1.
  endif

! Vertical loop over T/V levels

  do k = kb, mza

! Average 3 earth velocity components from T points to V face

     vxeface = .5 * (vxe(k,iw1) + vxe(k,iw2))
     vyeface = .5 * (vye(k,iw1) + vye(k,iw2))
     vzeface = .5 * (vze(k,iw1) + vze(k,iw2))

! Project earth velocity components at V face onto U and W directions

     ufacev(k) = unx(iv) * vxeface + uny(iv) * vyeface + unz(iv) * vzeface
     wfacev    = wnx_v   * vxeface + wny_v   * vyeface + wnz_v   * vzeface

! Compute z displacement component for V face relative to T point

     dzps(k)  = min( max( -dto2 * wfacev, -dzto4(k) ), dzto4(k) )

! Compute X, Y displacement components for V face relative to T point

     if (vs(k) > 0.0) then

        dxps(k) = -dto2 * (vs(k) * cosv1 - ufacev(k) * sinv1) + dxps1
        dyps(k) = -dto2 * (vs(k) * sinv1 + ufacev(k) * cosv1) + dyps1

     else

        dxps(k) = -dto2 * (vs(k) * cosv2 - ufacev(k) * sinv2) + dxps2
        dyps(k) = -dto2 * (vs(k) * sinv2 + ufacev(k) * cosv2) + dyps2

     endif

  enddo

end subroutine donorpointv


!===============================================================================


subroutine donorpointw(iw, dt, ws, vxe, vye, vze, dxps, dyps, dzps)

  use mem_ijtabs, only: itab_w, jtw_prog
  use mem_grid,   only: mza, lpw, dzto2

  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: dt
  real,    intent(in)  :: ws(mza)
  real,    intent(in)  :: vxe(mza), vye(mza), vze(mza)
  real,    intent(out) :: dxps(mza), dyps(mza), dzps(mza)

  real :: dto2
  real :: unx_w, uny_w, vnx_w, vny_w, vnz_w
  real :: vxeface, vyeface, vzeface
  real :: uface,vface

  real, parameter :: onethird = 1.0 / 3.0

  integer :: kb, k

  kb = lpw(iw)

  dto2 = 0.5 * dt

  unx_w = itab_w(iw)%unx_w
  uny_w = itab_w(iw)%uny_w

  vnx_w = itab_w(iw)%vnx_w
  vny_w = itab_w(iw)%vny_w
  vnz_w = itab_w(iw)%vnz_w

! Vertical loop over W levels

  do k = kb, mza-1

! Average 3 velocity components from T points to W point

     vxeface = .5 * (vxe(k) + vxe(k+1))
     vyeface = .5 * (vye(k) + vye(k+1))
     vzeface = .5 * (vze(k) + vze(k+1))

! Project earth velocity components at W face onto U and V directions

     uface = unx_w * vxeface + uny_w * vyeface
     vface = vnx_w * vxeface + vny_w * vyeface + vnz_w * vzeface

! Compute x and y displacement components for W face relative to T point

     dxps(k) = -dto2 * uface
     dyps(k) = -dto2 * vface

! Compute z displacement components for W face relative to T point

     if (ws(k) > 0.0) then
        dzps(k) = max( -dto2 * ws(k) + dzto2(k  ), 0.0)
     else
        dzps(k) = min( -dto2 * ws(k) - dzto2(k+1), 0.0)
     endif

  enddo

end subroutine donorpointw

!===============================================================================

subroutine comp_cfls_long(vmsca, wmsca, rho_old)

  ! Diagnose outflow CFL number; also used for advection CFL stability check

  use mem_basic,   only: wc, vc
  use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
  use mem_grid,    only: lpv, lpw, volti, glatw, glonw, mza, mva, mwa
  use consts_coms, only: r8
  use misc_coms,   only: iparallel, dtlm
  use mem_para,    only: myrank, mgroupsize
  !$ use omp_lib,  only: omp_get_max_threads, omp_get_thread_num

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  real,     intent(in) :: vmsca  (mza,mva)
  real,     intent(in) :: wmsca  (mza,mwa)
  real(r8), intent(in) :: rho_old(mza,mwa)

  integer :: j, iv, k, iw, n
  integer :: ier, inode, iwnode
  integer :: mlocc, mlocw, nthreads, myid
  real    :: dt
  real    :: cfl(mza)

  integer, allocatable :: imax(:,:), iwmax(:,:)
  real,    allocatable :: cfl_max(:), w_max(:)

  real    :: cflmx, wmx
  integer :: jcmax(2), jwmax(2)

#ifdef OLAM_MPI
  real :: sbuf(6), rbuf(6,mgroupsize)
#endif

  myid     = 1
  nthreads = 1
  !$ nthreads = omp_get_max_threads()

  allocate(cfl_max(nthreads))
  allocate(imax (2,nthreads))

  allocate(w_max  (nthreads))
  allocate(iwmax(2,nthreads))

  !$omp parallel private(myid,cfl,cflmx,wmx,jcmax,jwmax) &
  !$omp          shared(cfl_max,w_max,imax,iwmax)

  !$ myid = omp_get_thread_num() + 1
  cflmx = -1.0
  wmx   =  0.0

  !$omp do private(iw,dt,k,n,iv)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     dt = dtlm

     do k = lpw(iw), mza
        cfl(k) = max(wmsca(k,iw),0.0) - min(wmsca(k-1,iw),0.0)
     enddo

     do n = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(n)
        do k = lpv(iv), mza
           cfl(k) = cfl(k) - min(itab_w(iw)%dirv(n) * vmsca(k,iv), 0.0)
        enddo
     enddo

     do k = lpw(iw), mza
        cfl(k) = cfl(k) * dt * volti(k,iw) / real(rho_old(k,iw))
     enddo

     do k = lpw(iw), mza

        if (cfl(k) > cflmx .or. cfl(k) /= cfl(k)) then
           cflmx      = cfl(k)
           jcmax(1:2) = [ k, itab_w(iw)%iwglobe ]
        endif

        if (abs(wc(k,iw)) > abs(wmx) .or. wc(k,iw) /= wc(k,iw)) then
           wmx        = wc(k,iw)
           jwmax(1:2) = [ k, itab_w(iw)%iwglobe ]
        endif

        if ( cfl(k) > 1.0 .or. cfl(k) /= cfl(k) ) then
           n = itab_w(iw)%npoly
           write(*,*)
           write(*,'(4(A,I0),2(A,f0.3),/,2(A,f0.3),A,7(f0.3,1x))')           &
                "!!! CFL VIOLATION at node ", myrank,                        &
                ", iw=", iw, ", iwglobe=", itab_w(iw)%iwglobe, ", k=", k,    &
                " lat=", glatw(iw), " lon=", glonw(iw),                      &
                "!!! CFL = ", cfl(k),                             &
                ", W = ", wc(k,iw), ", VC = ", vc(k,itab_w(iw)%iv(1:n))
        endif

     enddo

  enddo
  !$omp end do nowait

  cfl_max(myid) = cflmx
  imax (1:2,myid) = jcmax(1:2)

  w_max  (myid) = wmx
  iwmax(1:2,myid) = jwmax(1:2)
  !$omp end parallel

  mlocc = 1
  mlocw = 1

  !$ if (nthreads > 1) then
  !$
  !$    mlocw = maxloc(abs(w_max), dim=1)
  !$    mlocc = maxloc(cfl_max, dim=1)
  !$
  !$    do n = 1, nthreads
  !$       if (cfl_max(n) /= cfl_max(n)) mlocc = n
  !$       if (  w_max(n) /=   w_max(n)) mlocw = n
  !$    enddo
  !$ endif

#ifdef OLAM_MPI
  if (iparallel == 1) then
     sbuf(1)   = cfl_max(mlocc)
     sbuf(2:3) = real(imax(:,mlocc))
     sbuf(4)   = w_max(mlocw)
     sbuf(5:6) = real(iwmax(:,mlocw))

     call MPI_Gather(sbuf, 6, MPI_REAL, rbuf, 6, MPI_REAL, &
                     0, MPI_COMM_WORLD, ier)
  endif
#endif

  if (myrank == 0) then
     inode  = 1
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
     endif
#endif

     write(*,'(5x,A,f0.3,3(A,I0))') "Max CFL = ", cfl_max(mlocc),  &
          " at node ", inode-1, ", iwglobe=", imax(2,mlocc), ", k=", imax(1,mlocc)

     write(*,'(5x,A,f0.3,3(A,I0))') "Max  W  = ", w_max(mlocw),  &
          " at node ", iwnode-1, ", iwglobe=", iwmax(2,mlocw), ", k=", iwmax(1,mlocw)
  endif

end subroutine comp_cfls_long
