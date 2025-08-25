module scalar_transport

  implicit none

  ! MPI communication tags
  integer, parameter :: itag_gxyps = 31
  integer, parameter :: itag_sct   = 32
  integer, parameter :: itag_scp   = 33
  integer, parameter :: itag_monot = 34
  integer, parameter :: itag_pd    = 35
  integer, parameter :: itag_cfl   = 36

  ! Monotonic / positive-definite tolerances
  real, parameter :: eps0 = 1.e-32
  real, parameter :: eps1 = 3.e-07
  real, parameter :: onep = 1.0 + eps1
  real, parameter :: onem = 1.0 - eps1

  ! For testing - may want to add to namelist!
  integer, parameter  :: iorderv        =  3
  logical, parameter  :: centered_monot = .false.

  private
  public :: scalar_transport_rk

contains

!===========================================================================

subroutine scalar_transport_rk(rho_old)

  use mem_ijtabs,   only: jtab_w, itab_w, jtv_wadj, jtw_prog, jtw_wadj
  use mem_grid,     only: mza, mwa, lpv, lpw, volti
  use misc_coms,    only: dtlm, iparallel, nrk_scal
  use var_tables,   only: num_scalar, scalar_tab
  use mem_turb,     only: akhodx, khtopv
  use oname_coms,   only: nl
  use olam_mpi_atm, only: mpi_post_direct_recv_w, mpi_post_direct_send_w, &
                          mpi_finish_direct_recv_w, mpi_finish_direct_send_w
  use obnd,         only: lbcopy_w
  use mem_basic,    only: rho, vmasc, wmasc
  use tridiag,      only: tridif_prep, tridif_fini
  use mem_nudge,    only: nudflag, rhot_nud
  use mem_para,     only: myrank

#ifdef OLAM_MPI
  use mpi_f08
#endif

  implicit none

  real, intent(in) :: rho_old(mza,mwa)

  integer  :: j,iw,istage
  integer  :: n,k,kb,iv,iwn,jv
  integer  :: iorder,imonot,i2d
  real     :: dtl, dtr
  real     :: dtlr(nrk_scal)
  real     :: rhof(mza)

! Automatic arrays:

  real :: scp0(mza,mwa)
  real :: rhos(mza,mwa,nrk_scal)
  real :: fluxdiv(mza)

  real    :: ff_implic(nrk_scal,mwa)
  integer :: nn_implic(nrk_scal)
  integer :: nn

  real, pointer, contiguous :: scp(:,:)
  real, pointer, contiguous :: sct(:,:)

  real :: c1, c2
  real :: delex_scp(mza), del_scp(mza), dscp0(mza)

  real :: wup(mza), wdn(mza)

  type implic_vars
     integer, allocatable :: iw_implic(:)
     integer, allocatable :: ii_implic(:)
     real,    allocatable :: b1(:,:)
     real,    allocatable :: b2(:,:)
     real,    allocatable :: b3(:,:)
  end type implic_vars

  type(implic_vars) :: vv_implic(nrk_scal)

  real :: b2a(mza), b3a(mza)
  integer :: kwtop(mwa)

  real, parameter :: one3    = 1. / 3.
  real, parameter :: frk3(3) = [ one3, 0.5, 1.0 ]
  real, parameter :: frk2(2) =       [ 0.5, 1.0 ]

  iorder = nl%scal_horiz_adv_order
  imonot = nl%iscal_monot
  dtl    = dtlm

  i2d = 0
  if (iorder == 2) i2d = 2
  if (iorder == 3) i2d = 5

  if (nrk_scal == 3) then
     do n = 1, nrk_scal
        dtlr(n) = dtl * frk3(n)
     enddo
  else
     do n = 1, nrk_scal
        dtlr(n) = dtl * frk2(n)
     enddo
  endif

  call check_cfls( vmasc, wmasc, rho_old, nrk_scal, ff_implic)

  nn_implic(:) = 0

  !$omp parallel do private(iw,kb,k,rhof,n) reduction(+:nn_implic)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

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
           rhos(k,iw,1) = 2./3. * rho_old(k,iw) + 1./3. * rhof(k)
           rhos(k,iw,2) = 0.5 * (rho_old(k,iw) + rhof(k))
           rhos(k,iw,3) = rhof(k)
        enddo

     else

        do k = kb, mza
           rhos(k,iw,1) = 0.5 * (rho_old(k,iw) + rhof(k))
           rhos(k,iw,2) = rhof(k)
        enddo

     endif

     do n = 1, nrk_scal
        if (ff_implic(n,iw) > 0.001) then
           nn_implic(n) = nn_implic(n) + 1
        endif
     enddo

     kwtop(iw) = maxval( khtopv( itab_w(iw)%iv( 1:itab_w(iw)%npoly ) ) )
  enddo
  !$omp end parallel do

  do n = 1, nrk_scal

     if (nn_implic(n) > 0) then

        allocate(vv_implic(n)%ii_implic(mwa))
        allocate(vv_implic(n)%iw_implic(nn_implic(n)))
        allocate(vv_implic(n)%b1   (mza,nn_implic(n)))
        allocate(vv_implic(n)%b2   (mza,nn_implic(n)))
        allocate(vv_implic(n)%b3   (mza,nn_implic(n)))

        nn = 0
        do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
           if (ff_implic(n,iw) > 0.001) then
              nn = nn + 1
              vv_implic(n)%ii_implic(iw) = nn
              vv_implic(n)%iw_implic(nn) = iw
           endif
        enddo

     endif
  enddo

  if (any(nn_implic(:) > 0)) then
     wup(mza) = 0.0
     wdn(mza) = 0.0

     do n = 1, nrk_scal
        do j = 1, nn_implic(n); iw = vv_implic(n)%iw_implic(j)

           kb = lpw(iw)

           wup(kb-1) = 0.0
           wdn(kb-1) = 0.0

           do k = kb, mza-1
              wup(k) = max(wmasc(k,iw), 0.)
              wdn(k) = min(wmasc(k,iw), 0.)
           enddo

           c1 = ff_implic(n,iw) * dtlr(n)

           do k = kb, mza
              c2 = c1 * volti(k,iw) / rhos(k,iw,n)

              vv_implic(n)%b1(k,j) =     - c2 * wup(k-1)
              b2a            (k)   = 1.0 + c2 * (wup(k) - wdn(k-1))
              b3a            (k)   =       c2 * wdn(k)
           enddo

           call tridif_prep(mza,kb,mza,vv_implic(n)%b1(:,j),b2a,b3a, &
                            vv_implic(n)%b2(:,j),vv_implic(n)%b3(:,j))
        enddo
     enddo
  endif

  !********************************************************
  ! MAIN LOOP OVER ALL SCALAR SPECIES
  !********************************************************

  do n = 1, num_scalar

     ! POINT SCP AND SCT TO SCALAR TABLE ARRAYS

     scp => scalar_tab(n)%var_p(:,:)
     sct => scalar_tab(n)%var_t(:,:)

     ! STORE INITIAL CONCENTRATIONS FOR R-K TIMESTEPPING. THIS NEEDS TO INCLUDE
     ! BORDER CELLS FOR THE MONOTONIC AND POSITIVE DEFINITE OPTIONS.

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_wadj)%jend; iw = jtab_w(jtw_wadj)%iw(j)
        scp0(:,iw) = scp(:,iw)
     enddo
     !$omp end parallel do

     ! UPDATE TRACER TENDENCIES TO INCLUDE HORIZONTAL DIFFUSION. ALL OF THE
     ! TRACER PHYSICS TENDENCIES IN OLAM ARE SPLIT FROM THE R-K TIMESTEPPING
     ! AND DONE FORWARD-IN-TIME. DIFFUSION COULD BE EASILY INTEGRATED WITH
     ! EACH R-K STEP TO MATCH ADVECTION IF NECESSARY.

     if ( scalar_tab(n)%do_sgsmix ) then   ! Flag to skip diffusion for this species

        !$omp parallel do private(iw,fluxdiv,jv,iv,iwn,k)
        do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

           if (kwtop(iw) >= lpw(iw)) then   ! Diffusion is active in this column

              fluxdiv(lpw(iw):kwtop(iw)) = 0.

              ! Loop over neighbor cells
              do jv = 1, itab_w(iw)%npoly
                 iv  = itab_w(iw)%iv(jv)
                 iwn = itab_w(iw)%iw(jv)

                 ! Vertical loop over T levels
                 do k = lpv(iv), khtopv(iv)
                    fluxdiv(k) = fluxdiv(k) + akhodx(k,iv) * (scp(k,iwn) - scp(k,iw))
                 enddo
              enddo

              ! Vertical loop over T levels
              do k = lpw(iw), kwtop(iw)
                 sct(k,iw) = sct(k,iw) + fluxdiv(k) * volti(k,iw)
              enddo

           endif

        enddo
        !$omp end parallel do

     endif

     ! TENDENCY ARRAY SCT NEEDS TO BE COMMUNICATED TO BORDER CELLS FOR
     ! THE MONOTONIC AND POSITIVE DEFINITE OPTIONS. THIS DOES NOT NEED
     ! TO BE COMPLETED UNTIL THE FINAL RK STEP.

     if ( iparallel == 1 .and. imonot > 0 .and. scalar_tab(n)%pdef ) then
        call mpi_post_direct_recv_w(sct, itag_sct)
        call mpi_post_direct_send_w(sct, itag_sct)
     endif

     !*****************************************************
     ! MAIN LOOP OVER ALL R-K STAGES
     !*****************************************************

     do istage = 1, nrk_scal

        dtr = dtlr(istage)

        ! FOR THE FINAL R-K STAGE WITH MONOTONIC OR POSITIVE DEFINITE OPTIONS,
        ! FINISH MPI COMMUNICATION OF TENDENCY ARRAY AND UPDATE THE TRACER FOR
        ! THE LOW-ORDER FLUXES FOLLOWING WANG, SKAMAROCK, AND FEINGOLD, 2009.

        if (istage == nrk_scal .and. imonot > 0 .and. scalar_tab(n)%pdef) then

           if ( iparallel == 1 ) then
              call mpi_finish_direct_recv_w(itag=itag_sct)
              call mpi_finish_direct_send_w(itag=itag_sct)
           endif
           call lbcopy_w( a1=sct )

           !$omp parallel do private(iw,k)
           do j = 1, jtab_w(jtw_wadj)%jend; iw = jtab_w(jtw_wadj)%iw(j)
              do k = lpw(iw), mza
                 scp0(k,iw) = scp0(k,iw) + dtr * sct(k,iw) / rho_old(k,iw)
              enddo
              scp0(2:lpw(iw)-1,iw) = scp0(lpw(iw),iw)
           enddo
           !$omp end parallel do

        endif

        ! FINISH COMMUNICATION OF TRACER BORDER CELLS FROM PREVIOUS R-K STAGE

        if (iparallel == 1 .and. istage > 1) then
           call mpi_finish_direct_recv_w(itag_scp)
           call mpi_finish_direct_send_w(itag_scp)
        endif
        call lbcopy_w( a1=scp )

        ! INTEGRATE SCALARS FORWARD IN TIME USING EXPLICIT ADVECTION AND PHYSICS TENDENCIES

        if (istage == nrk_scal .and. imonot == 1) then

           call update_scalar_monot( scp, scp0, rhos(:,:,nrk_scal), vmasc, wmasc, &
                                     dtr, iorder, i2d, istage, nrk_scal, ff_implic )

        elseif (istage == nrk_scal .and. imonot == 2 .and. scalar_tab(n)%pdef) then

           call update_scalar_pd( scp, scp0, rhos(:,:,nrk_scal), vmasc, wmasc, &
                                  dtr, iorder, i2d, istage, nrk_scal, ff_implic )
        else

           call update_scalar( scp, scp0, sct, rhos(:,:,nrk_scal), vmasc, wmasc, &
                               dtr, iorder, i2d, istage, nrk_scal, ff_implic )
        endif

        ! POST RECEIVE OF TRACER BORDER CELLS

        if (istage < nrk_scal) call mpi_post_direct_recv_w(scp, itag_scp)

        ! IMPLICIT VERTICAL TRANSPORT

        !$omp parallel do private(iw,kb,dscp0,k,c1,delex_scp,del_scp)
        do j = 1, nn_implic(istage); iw = vv_implic(istage)%iw_implic(j)

           kb = lpw(iw)

           dscp0(mza ) = 0.
           dscp0(kb-1) = 0.

           do k = kb, mza-1
              dscp0(k) = scp0(k+1,iw) - scp0(k,iw)
           enddo

           do k = kb, mza
              c1 = ff_implic(istage,iw) * dtr * volti(k,iw) / rhos(k,iw,istage)

              delex_scp(k) = scp(k,iw) - scp0(k,iw) - c1 * ( min(wmasc(k  ,iw),0.) * dscp0(k  ) &
                                                           + max(wmasc(k-1,iw),0.) * dscp0(k-1) )
           enddo

           call tridif_fini(mza, kb, mza, vv_implic(istage)%b1(:,j), vv_implic(istage)%b2(:,j), &
                                          vv_implic(istage)%b3(:,j), delex_scp, del_scp )
           do k = kb, mza
              scp(k,iw) = scp0(k,iw) + del_scp(k)
           enddo
           scp(2:kb-1,iw) = scp(kb,iw)


        enddo
        !$omp end parallel do

        ! COMMUNICATE UPDATED TRACER BORDER CONCENTRATIONS FOR NEXT R-K STAGE

        if (istage < nrk_scal) call mpi_post_direct_send_w(scp, itag_scp)

     enddo ! istage

  enddo ! n

end subroutine scalar_transport_rk

!===========================================================================

subroutine update_scalar( scp, scp0, sct, rhos, vmasc, wmasc, &
                          dtr, iorderh, i2d, istage, nrk, ff_implic )

  use mem_grid,     only: mza, mwa, mva
  use mem_ijtabs,   only: jtab_v, jtab_w, itab_w, jtv_wadj, jtw_prog, &
                          iip, mloops, jtw_lbcp
  use misc_coms,    only: iparallel
  use grad_lib,     only: grad_t2d, grad_t2d_quadratic
  use olam_mpi_atm, only: mpi_post_direct_recv_w, mpi_post_direct_send_w, &
                          mpi_finish_direct_recv_w, mpi_finish_direct_send_w

  implicit none

  integer, intent(in)    :: iorderh, i2d, istage, nrk
  real,    intent(inout) :: scp  (mza,mwa)
  real,    intent(in)    :: scp0 (mza,mwa)
  real,    intent(inout) :: sct  (mza,mwa)
  real,    intent(in)    :: rhos (mza,mwa)
  real,    intent(in)    :: wmasc(mza,mwa)
  real,    intent(in)    :: vmasc(mza,mva)
  real,    intent(in)    :: dtr, ff_implic(nrk,mwa)

  integer :: j, iv, iw, i, ii, iwp

  real :: scp_upv(mza,mva)
  real :: gxyps  (mza,mwa,i2d)

  if (iparallel == 1) then
     if (iorderh > 1) call mpi_post_direct_recv_w(gxyps, itag_gxyps)
  endif

  !$omp parallel private(ii)

  ! COMPUTE HORIZONTAL POLYNOMIAL RECONSTRUCTION COEFFICIENTS AT EACH W CELL

  if (iorderh > 1) then

     !$omp do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (iorderh == 2) then
           call grad_t2d(iw, scp, gxyps(:,iw,1), gxyps(:,iw,2))
        elseif (iorderh == 3) then
           call grad_t2d_quadratic(iw, scp, gxyps)
        endif

     enddo
     !$omp end do

     if (iparallel == 1) then
        !$omp single
        call mpi_post_direct_send_w(gxyps, itag_gxyps)
        !$omp end single nowait
     endif

     if (jtab_w(jtw_lbcp)%jend > 1) then

        !$omp do private(iw,iwp,i)
        do j = 1, jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
           iwp = itab_w(iw)%iwp
           do i = 1, i2d
              gxyps(:,iw,i) = gxyps(:,iwp,i)
           enddo
        enddo
        !$omp end do

     endif

  endif

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH V FACE

  ! Split V loop into an "interior" loop that can overlap with communication
  ! and a "border" loop that requires MPI communication be completed.
  do ii = 1, iip

     if (ii == 2 .and. iorderh > 1) then
        !$omp single
        call mpi_finish_direct_recv_w(itag_gxyps)
        !$omp end single
     endif

     !$omp do private(iv)
     do j = 1, jtab_v(mloops+ii)%jend; iv = jtab_v(mloops+ii)%iv(j)

        call scalar_v_column( iv )

     enddo
     !$omp end do

  enddo  ! interior/border loop

  ! MAIN W LOOP TO COMPUTE VERTICAL FLUXES AND MARCH TRACER TO NEXT R-K STAGE

  !$omp do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_w_column( iw, ff_implic(istage,iw) )

  enddo
  !$omp end do nowait
  !$omp end parallel

  if (iparallel == 1) then
     if (iorderh > 1) call mpi_finish_direct_send_w(itag_gxyps)
  endif


contains


  subroutine scalar_v_column( iv )

    use mem_grid,   only: lpv, mza
    use mem_ijtabs, only: itab_v
    use mem_adv,    only: xx0_v, yy0_v, xy0_v
!   import,         only: scp, scp_upv, vmasc, gxyps, iorderh

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
       if (vmasc(k,iv) < 0.) scp_upv(k,iv) = scp2(k)
    enddo

  end subroutine scalar_v_column


  subroutine scalar_w_column( iw, fimpl )

    use mem_grid,   only: lpw, lpv, mza, volti
    use mem_ijtabs, only: itab_w
    use grad_lib,   only: grad_z_linear, grad_z_quadratic
!   import,         only: scp, wmasc, vmasc, scp_upv, scp0, sct, rhos, dtr, iorderv

    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: fimpl
    integer             :: kb, k, jv, iv
    real                :: scpb(mza), scpt(mza)
    real                :: scp_upw(mza)
    real                :: fluxdiv(mza)
    real                :: fexpl

    kb = lpw(iw)

    fexpl = 1.0 - fimpl

    scp_upw(kb-1) = scp(kb ,iw)
    scp_upw(mza)  = scp(mza,iw)

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
       scp_upw(k) = scpt(k)
       if (wmasc(k,iw) < 0.) scp_upw(k) = scpb(k+1)
    enddo

    ! Scalar horizontal flux divergence

    fluxdiv(:) = 0.

    ! Loop over neighbor V points of this W cell
    do jv = 1, itab_w(iw)%npoly
       iv = itab_w(iw)%iv(jv)

       do k = lpv(iv), mza
          fluxdiv(k) = fluxdiv(k) &
                     + itab_w(iw)%dirv(jv) * vmasc(k,iv) * (scp_upv(k,iv) - scp0(k,iw))
       enddo
    enddo

    do k = kb, mza

       ! Scalar vertical flux divergence
       fluxdiv(k) = fluxdiv(k) + fexpl * ( wmasc(k-1,iw) * (scp_upw(k-1) - scp0(k,iw)) &
                                         - wmasc(k  ,iw) * (scp_upw(k  ) - scp0(k,iw)) )

       ! Integrate forward in time
       scp(k,iw) = scp0(k,iw) + dtr * ( sct(k,iw) + volti(k,iw) * fluxdiv(k) ) / rhos(k,iw)

    enddo

    scp(1:kb-1,iw) = scp(kb,iw)

  end subroutine scalar_w_column


end subroutine update_scalar

!===========================================================================

subroutine update_scalar_monot( scp, scp0, rhos, vmasc, wmasc, &
                                dtr, iorderh, i2d, istage, nrk, ff_implic )

  use mem_grid,     only: mza, mwa, mva, lpv
  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, jtw_prog, &
                          iip, mloops, jtw_lbcp
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_post_direct_recv_w, mpi_finish_direct_recv_w, &
                          mpi_post_direct_send_w, mpi_finish_direct_send_w
  use grad_lib,     only: grad_t2d, grad_t2d_quadratic

  implicit none

  integer,          intent(in)    :: iorderh, i2d, istage, nrk
  real,    intent(inout) :: scp     (mza,mwa)
  real,    intent(in)    :: scp0    (mza,mwa)
  real,    intent(in)    :: rhos    (mza,mwa)
  real,    intent(in)    :: wmasc   (mza,mwa)
  real,    intent(in)    :: vmasc   (mza,mva)
  real,    intent(in)    :: dtr, ff_implic(nrk,mwa)

  integer :: j, iv, iw, iw1, iw2, k, ii, i, iwp
  real    :: scale

  real :: sfluxvh    (mza,mva)
  real :: scale_inout(mza,mwa,2)
  real :: gxyps      (mza,mwa,i2d)
  real :: scp_upv    (mza,mva)
  real :: scp_hiv    (mza,mva)
  real :: scp_upw    (mza,mwa)

  if (iparallel == 1) then
     if (iorderh > 1) call mpi_post_direct_recv_w(gxyps,       itag_gxyps)
     call                  mpi_post_direct_recv_w(scale_inout, itag_monot)
  endif

  !$omp parallel private(ii)

  ! COMPUTE HORIZONTAL POLYNOMIAL RECONSTRUCTION COEFFICIENTS AT EACH W CELL

  if (iorderh > 1) then

     !$omp do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (iorderh == 2) then
           call grad_t2d(iw, scp, gxyps(:,iw,1), gxyps(:,iw,2))
        elseif (iorderh == 3) then
           call grad_t2d_quadratic(iw, scp, gxyps, bounded=.false.)
        endif

     enddo
     !$omp end do

     if (iparallel == 1) then
        !$omp single
        call mpi_post_direct_send_w(gxyps, itag_gxyps)
        !$omp end single nowait
     endif

     if (jtab_w(jtw_lbcp)%jend > 1) then

        !$omp do private(iw,iwp,i)
        do j = 1, jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
           iwp = itab_w(iw)%iwp
           do i = 1, i2d
              gxyps(:,iw,i) = gxyps(:,iwp,i)
           enddo
        enddo
        !$omp end do

     endif

  endif

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH V FACE

  ! Split V loop into an "interior" loop that can overlap with communication
  ! and a "border" loop that requires MPI communication be completed.
  do ii = 1, iip

     if (ii == 2 .and. iorderh > 1) then
        !$omp single
        call mpi_finish_direct_recv_w(itag_gxyps)
        !$omp end single
     endif

     !$omp do private(iv)
     do j = 1, jtab_v(mloops+ii)%jend; iv = jtab_v(mloops+ii)%iv(j)

        call scalar_vflux_monot( iv )

     enddo
     !$omp end do

  enddo  ! interior/border loop

  ! LOOP TO COMPUTE VERTICAL FLUXES AND SCALE FACTOR

  !$omp do private (iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_wflux_and_scales_monot( iw, ff_implic(istage,iw) )

  enddo
  !$omp end do

  ! Post MPI send of flux scale factors

  if (iparallel == 1) then
     !$omp single
     call mpi_post_direct_send_w(scale_inout, itag_monot)
     !$omp end single nowait
  endif

  ! Lateral boundary copy of scale factor

  if (jtab_w(jtw_lbcp)%jend > 1) then

     !$omp do private(iw,iwp)
     do j = 1, jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
        iwp = itab_w(iw)%iwp
        scale_inout(:,iw,1) = scale_inout(:,iwp,1)
        scale_inout(:,iw,2) = scale_inout(:,iwp,2)
     enddo
     !$omp end do

  endif

  ! APPLY FLUX RENORMALIZATION TO HORIZONTAL FLUXES

  ! Split V loop into an "interior" loop that can overlap with communication
  ! and a "border" loop that requires MPI communication be completed.
  do ii = 1, iip

     if (ii == 2) then
        !$omp single
        call mpi_finish_direct_recv_w(itag_monot)
        !$omp end single
     endif

     !$omp do private(iv,iw1,iw2,k,scale)
     do j = 1, jtab_v(mloops+ii)%jend; iv = jtab_v(mloops+ii)%iv(j)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        do k = lpv(iv), mza
           scale = min( scale_inout(k,iw2,1), scale_inout(k,iw1,2) )
           if (sfluxvh(k,iv) < 0.) scale = min( scale_inout(k,iw1,1), scale_inout(k,iw2,2) )
           scp_upv(k,iv) = scp_upv(k,iv) + scp_hiv(k,iv) * scale
        enddo
     enddo
     !$omp end do

  enddo  ! interior/border loop

  ! COMPUTE FLUX DIVERGENCE AND FORWARD INTEGRATE SCALAR CONCENTRATION

  !$omp do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_tend_monot( iw, ff_implic(istage,iw) )

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! Make sure MPI sendss are completed before we leave this routine and
  ! automatic memory is released

  if (iparallel == 1) then
     if (iorderh > 1) call mpi_finish_direct_send_w(itag_gxyps)
     call                  mpi_finish_direct_send_w(itag_monot)
  endif


contains


  subroutine scalar_vflux_monot( iv )

    use mem_grid,   only: lpv, mza
    use mem_ijtabs, only: itab_v
    use mem_adv,    only: xx0_vu, yy0_vu, xy0_vu
!   import,         only: scp, scp0, scp_upv, scp_hiv, vmasc, sfluxvh, gxyps, &
!                         iorderh, centered_monot
    implicit none

    integer, intent(in) :: iv
    integer             :: iw1, iw2, k
    real                :: scp1(mza), scp2(mza)

    iw1 = itab_v(iv)%iw(1)
    iw2 = itab_v(iv)%iw(2)

     ! Low-order horizontal fluxes

    do k = lpv(iv), mza
       scp_upv(k,iv) = scp0(k,iw1)
       if (vmasc(k,iv) < 0.) scp_upv(k,iv) = scp0(k,iw2)
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
          sfluxvh(k,iv) = vmasc(k,iv) * scp_hiv(k,iv)
       enddo

    else

       do k = lpv(iv), mza
          scp_hiv(k,iv) = scp1(k)
          if (vmasc(k,iv) < 0.) scp_hiv(k,iv) = scp2(k)

          scp_hiv(k,iv) = scp_hiv(k,iv) - scp_upv(k,iv)
          sfluxvh(k,iv) = vmasc(k,iv) * scp_hiv(k,iv)
       enddo

    endif

  end subroutine scalar_vflux_monot


  subroutine scalar_wflux_and_scales_monot( iw, fimpl )

    use mem_grid,   only: lpw, lpv, mza, volti
    use mem_ijtabs, only: itab_w
    use grad_lib,   only: grad_z_linear, grad_z_quadratic
    use tridiag,    only: tridif_fini
!   import,         only: scp, scp0, wmasc, vmasc, scp_upv, scp_upw, sfluxvh, onep, &
!                         scale_inout, rhos, dtr, iorderv, centered_monot, eps0, eps1
    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: fimpl

    integer             :: kb, k, jv, iv, iwn
    real                :: rscp_low, flux_in, flux_out, sc_in, sc_out
    real                :: scale, dtrp, fexpl, epsr
    real                :: scpb(mza), scpt(mza), scp_hiw(mza)
    real                :: sfluxwh(mza), smax(mza), smin(mza)
    real                :: fluxdiv_low(mza)
    real                :: fluxdiv_in (mza)
    real                :: fluxdiv_out(mza)
    real                :: imp_correct, dwi

    dtrp = onep * dtr

    kb = lpw(iw)

    fexpl = 1. - fimpl

    ! Low order vertical flux

    scp_upw(kb-1,iw) = scp0(kb ,iw)
    scp_upw(mza ,iw) = scp0(mza,iw)

    do k = kb, mza-1
       scp_upw(k,iw) = scp0(k,iw)
       if (wmasc(k,iw) < 0.) scp_upw(k,iw) = scp0(k+1,iw)
    enddo

    ! Higher-order vertical fluxes

    if (iorderv == 1) then

       do k = kb, mza
          scpb(k) = scp(k,iw)
          scpt(k) = scp(k,iw)
       enddo

    elseif (iorderv == 2) then

       call grad_z_linear(iw, scp(:,iw), scpb, scpt, itopbc=1, ibotbc=1)

    elseif (iorderv == 3) then

       call grad_z_quadratic(iw, scp(:,iw), scpb, scpt, itopbc=1, ibotbc=1, bounded=.false.)

    endif

    sfluxwh(kb-1) = 0.
    sfluxwh(mza)  = 0.

    if (centered_monot) then

       do k = kb, mza-1
          scp_hiw(k) = 0.5 * (scpt(k) + scpb(k+1)) - scp_upw(k,iw)
          sfluxwh(k) = fexpl * wmasc(k,iw) * scp_hiw(k)
       enddo

    else

       do k = kb, mza-1
          scp_hiw(k) = scpt(k)
          if (wmasc(k,iw) < 0.) scp_hiw(k) = scpb(k+1)

          scp_hiw(k) = scp_hiw(k) - scp_upw(k,iw)
          sfluxwh(k) = fexpl * wmasc(k,iw) * scp_hiw(k)
       enddo

    endif

     ! Find upwind-biased tracer max/min for each cell in this column

     do k = kb, mza
        smax(k) = scp0(k,iw)
        smin(k) = scp0(k,iw)
     enddo

     do k = kb+1, mza
        if (wmasc(k-1,iw) > 0.) then
           smax(k) = max(smax(k), scp0(k-1,iw))
           smin(k) = min(smin(k), scp0(k-1,iw))
        endif
     enddo

     do k = kb, mza-1
        if (wmasc(k,iw) < 0.) then
           smax(k) = max(smax(k), scp0(k+1,iw))
           smin(k) = min(smin(k), scp0(k+1,iw))
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
                          + itab_w(iw)%dirv(jv) * vmasc(k,iv) * (scp_upv(k,iv) - scp0(k,iw))

           fluxdiv_in (k) = fluxdiv_in (k) + max( itab_w(iw)%dirv(jv) * sfluxvh(k,iv), 0.)
           fluxdiv_out(k) = fluxdiv_out(k) + min( itab_w(iw)%dirv(jv) * sfluxvh(k,iv), 0.)

           ! Upwind-biased tracer max/min
           if ( itab_w(iw)%dirv(jv) * vmasc(k,iv) > 0.) then
              smax(k) = max(smax(k), scp0(k,iwn))
              smin(k) = min(smin(k), scp0(k,iwn))
           endif

        enddo
     enddo

     ! Scalar vertical flux divergence and flux renormalization scales

     if (fimpl > .001) then

        do k = kb, mza
           dwi = fimpl * (wmasc(k,iw) - wmasc(k-1,iw))

           fluxdiv_low(k) = fluxdiv_low(k) + fexpl * ( wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw))   &
                                                     - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw)) ) &
                                           + scp0(k,iw) * dwi

           fluxdiv_in (k) = fluxdiv_in (k) + max(sfluxwh(k-1),0.) - min(sfluxwh(k),0.)
           fluxdiv_out(k) = fluxdiv_out(k) + min(sfluxwh(k-1),0.) - max(sfluxwh(k),0.)

           imp_correct = rhos(k,iw) + dtr * volti(k,iw) * dwi
           smax(k) = imp_correct * smax(k)
           smin(k) = imp_correct * smin(k)
        enddo

     else

        do k = kb, mza
           fluxdiv_low(k) = fluxdiv_low(k) + wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw)) &
                                           - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw))

           fluxdiv_in (k) = fluxdiv_in (k) + max(sfluxwh(k-1),0.) - min(sfluxwh(k),0.)
           fluxdiv_out(k) = fluxdiv_out(k) + min(sfluxwh(k-1),0.) - max(sfluxwh(k),0.)

           smax(k) = rhos(k,iw) * smax(k)
           smin(k) = rhos(k,iw) * smin(k)
        enddo

     endif

    ! Flux limiter / renormalization

    do k = lpw(iw), mza
       rscp_low = scp0(k,iw) * rhos(k,iw) + dtr * volti(k,iw) * fluxdiv_low(k)

       flux_in  = max(dtrp * volti(k,iw) * fluxdiv_in (k),  eps0)
       flux_out = min(dtrp * volti(k,iw) * fluxdiv_out(k), -eps0)

       epsr = eps1 * abs(rscp_low)

       sc_in  = max(smax(k) - (rscp_low + epsr) - eps0, 0.)
       sc_out = min(smin(k) - (rscp_low - epsr) + eps0, 0.)

       scale_inout(k,iw,1) = min(1., sc_in  / max(flux_in , 0.01*sc_in ) )
       scale_inout(k,iw,2) = min(1., sc_out / min(flux_out, 0.01*sc_out) )
    enddo

     ! Apply flux renormalization to vertical fluxes

     do k = lpw(iw), mza-1
        scale = min( scale_inout(k+1,iw,1), scale_inout(k,iw,2) )
        if (sfluxwh(k) < 0.) scale = min( scale_inout(k,iw,1), scale_inout(k+1,iw,2))
        scp_upw(k,iw) = scp_upw(k,iw) + scp_hiw(k) * scale
     enddo

  end subroutine scalar_wflux_and_scales_monot


  subroutine scalar_tend_monot( iw, fimpl )

    use mem_grid,   only: mza, lpw, lpv, volti
    use mem_ijtabs, only: itab_w
!   import,         only: scp, scp0, wmasc, vmasc, scp_upv, scp_upw, rhos, dtr

    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: fimpl

    integer             :: kb, jv, iv, k
    real                :: fluxdiv(mza), fexpl

    kb = lpw(iw)
    fluxdiv(:) = 0.0

    ! Scalar horizontal flux divergence

    do jv = 1, itab_w(iw)%npoly
       iv = itab_w(iw)%iv(jv)

       do k = lpv(iv), mza
          fluxdiv(k) = fluxdiv(k) &
                     + itab_w(iw)%dirv(jv) * vmasc(k,iv) * (scp_upv(k,iv) - scp0(k,iw))
       enddo
    enddo

    if (fimpl > .001) then
       fexpl = 1. - fimpl

       do k = kb, mza
          ! Scalar vertical flux divergence
          fluxdiv(k) = fluxdiv(k) + fexpl * ( wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw)) &
                                            - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw)) )
          ! Integrate forward in time
          scp(k,iw) = scp0(k,iw) + dtr * volti(k,iw) * fluxdiv(k) / rhos(k,iw)
       enddo

    else

       do k = kb, mza
          ! Scalar vertical flux divergence
          fluxdiv(k) = fluxdiv(k) + wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw)) &
                                  - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw))

          ! Integrate forward in time
          scp(k,iw) = scp0(k,iw) + dtr * volti(k,iw) * fluxdiv(k) / rhos(k,iw)
       enddo

       scp(2:kb-1,iw) = scp(kb,iw)
    endif

  end subroutine scalar_tend_monot


end subroutine update_scalar_monot

!===========================================================================

subroutine update_scalar_pd( scp, scp0, rhos, vmasc, wmasc, &
                             dtr, iorderh, i2d, istage, nrk, ff_implic)

  use mem_grid,     only: mza, mwa, mva, lpv
  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, jtw_prog, &
                          iip, mloops, jtw_lbcp
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_post_direct_recv_w, mpi_finish_direct_recv_w, &
                          mpi_post_direct_send_w, mpi_finish_direct_send_w
  use grad_lib,     only: grad_t2d, grad_t2d_quadratic

  implicit none

  integer, intent(in)    :: iorderh, i2d, istage, nrk
  real,    intent(inout) :: scp  (mza,mwa)
  real,    intent(in)    :: scp0 (mza,mwa)
  real,    intent(in)    :: rhos (mza,mwa)
  real,    intent(in)    :: vmasc(mza,mva)
  real,    intent(in)    :: wmasc(mza,mwa)
  real,    intent(in)    :: dtr, ff_implic(nrk,mwa)

  integer :: j, iv, iw, iw1, iw2, k, ii, i, iwp
  real    :: scale

  real :: sfluxvh  (mza,mva)
  real :: scale_out(mza,mwa)
  real :: gxyps    (mza,mwa,i2d)
  real :: scp_upv  (mza,mva)
  real :: scp_hiv  (mza,mva)
  real :: scp_upw  (mza,mwa)

  if (iparallel == 1) then
     if (iorderh > 1) call mpi_post_direct_recv_w(gxyps,     itag_gxyps)
     call                  mpi_post_direct_recv_w(scale_out, itag_pd)
  endif

  !$omp parallel private(ii)

  ! COMPUTE HORIZONTAL POLYNOMIAL RECONSTRUCTION COEFFICIENTS AT EACH W CELL

  if (iorderh > 1) then

     !$omp do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        if (iorderh == 2) then
           call grad_t2d(iw, scp, gxyps(:,iw,1), gxyps(:,iw,2))
        elseif (iorderh == 3) then
           call grad_t2d_quadratic(iw, scp, gxyps)
        endif

     enddo
     !$omp end do

     if (iparallel == 1) then
        !$omp single
        call mpi_post_direct_send_w(gxyps, itag_gxyps)
        !$omp end single nowait
     endif

     if (jtab_w(jtw_lbcp)%jend > 1) then

        !$omp do private(iw,iwp,i)
        do j = 1, jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
           iwp = itab_w(iw)%iwp
           do i = 1, i2d
              gxyps(:,iw,i) = gxyps(:,iwp,i)
           enddo
        enddo
        !$omp end do

     endif

  endif

  ! COMPUTE UPWINDED TRACER CONCENTRATIONS AT EACH V FACE

  ! Split V loop into an "interior" loop that can overlap with communication
  ! and a "border" loop that requires MPI communication be completed.
  do ii = 1, iip

     if (ii == 2 .and. iorderh > 1) then
        !$omp single
        call mpi_finish_direct_recv_w(itag_gxyps)
        !$omp end single
     endif

     !$omp do private(iv)
     do j = 1, jtab_v(mloops+ii)%jend; iv = jtab_v(mloops+ii)%iv(j)

        call scalar_vflux_pd( iv )

     enddo
     !$omp end do

  enddo  ! interior/border loop

  ! LOOP TO COMPUTE VERTICAL FLUXES AND SCALE FACTOR

  !$omp do private (iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_wflux_and_scale_pd( iw, ff_implic(istage,iw) )

  enddo
  !$omp end do

  ! Post MPI send of flux scale factors

  if (iparallel == 1) then
     !$omp single
     call mpi_post_direct_send_w(scale_out, itag_pd)
     !$omp end single nowait
  endif

  ! Lateral boundary copy of scale factor

  if (jtab_w(jtw_lbcp)%jend > 1) then

     !$omp do private(iw,iwp)
     do j = 1, jtab_w(jtw_lbcp)%jend; iw = jtab_w(jtw_lbcp)%iw(j)
        iwp = itab_w(iw)%iwp
        scale_out(:,iw) = scale_out(:,iwp)
     enddo
     !$omp end do

  endif

  ! APPLY FLUX RENORMALIZATION TO HORIZONTAL FLUXES

  ! Split V loop into an "interior" loop that can overlap with communication
  ! and a "border" loop that requires MPI communication be completed.
  do ii = 1, iip

     if (ii == 2) then
        !$omp single
        call mpi_finish_direct_recv_w(itag_pd)
        !$omp end single
     endif

     !$omp do private(iv,iw1,iw2,k,scale)
     do j = 1, jtab_v(mloops+ii)%jend; iv = jtab_v(mloops+ii)%iv(j)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        do k = lpv(iv), mza
           scale = scale_out(k,iw1)
           if (sfluxvh(k,iv) < 0.) scale = scale_out(k,iw2)
           scp_upv(k,iv) = scp_upv(k,iv) + scp_hiv(k,iv) * scale
        enddo
     enddo
     !$omp end do

  enddo  ! interior/border loop

  ! COMPUTE FLUX DIVERGENCE AND FORWARD INTEGRATE SCALAR CONCENTRATION

  !$omp do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     call scalar_tend_pd( iw, ff_implic(istage,iw) )

  enddo
  !$omp end do nowait
  !$omp end parallel

  ! Make sure MPI sendss are completed before we leave this routine and
  ! automatic memory is released

  if (iparallel == 1) then
     if (iorderh > 1) call mpi_finish_direct_send_w(itag_gxyps)
     call                  mpi_finish_direct_send_w(itag_pd)
  endif


contains


  subroutine scalar_vflux_pd( iv )

    use mem_grid,   only: lpv, mza
    use mem_ijtabs, only: itab_v
    use mem_adv,    only: xx0_v, yy0_v, xy0_v
!   import,         only: scp, scp0, scp_upv, scp_hiv, vmasc, sfluxvh, gxyps, iorderh

    implicit none

    integer, intent(in) :: iv
    integer             :: iw1, iw2, k
    real                :: scp1(mza), scp2(mza)

    iw1 = itab_v(iv)%iw(1)
    iw2 = itab_v(iv)%iw(2)

    ! High-order horizontal fluxes

    if (iorderh == 1) then

       do k = lpv(iv), mza
          scp1(k) = scp(k,iw1)
          scp2(k) = scp(k,iw1)
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

    ! Low-order horizontal fluxes and high-order corrections

    do k = lpv(iv), mza
       scp_hiv(k,iv) = scp1(k)
       if (vmasc(k,iv) < 0.) scp_hiv(k,iv) = scp2(k)

       scp_upv(k,iv) = scp0(k,iw1)
       if (vmasc(k,iv) < 0.) scp_upv(k,iv) = scp0(k,iw2)

       scp_hiv(k,iv) = scp_hiv(k,iv) - scp_upv(k,iv)
       sfluxvh(k,iv) = vmasc(k,iv) * scp_hiv(k,iv)
    enddo

  end subroutine scalar_vflux_pd


  subroutine scalar_wflux_and_scale_pd( iw, fimpl )

    use mem_grid,   only: lpw, lpv, mza, volti
    use mem_ijtabs, only: itab_w
    use grad_lib,   only: grad_z_linear, grad_z_quadratic
!   import,         only: scp, scp0, wmasc, vmasc, scp_upv, scp_upw, eps0, &
!                         sfluxvh, scale_out, rhos, dtr, iorderv, onep, onem
    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: fimpl

    integer             :: kb, k, jv, iv
    real                :: rscp_low, flux_out, sc_out
    real                :: scale, mdtrp, fexpl
    real                :: scpb(mza), scpt(mza), scp_hiw(mza)
    real                :: sfluxwh(mza)
    real                :: fluxdiv_low(mza)
    real                :: fluxdiv_out(mza)

    mdtrp = -onep * dtr

    kb = lpw(iw)

    fexpl = 1. - fimpl

    ! Low order vertical flux

    scp_upw(kb-1,iw) = scp(kb ,iw)
    scp_upw(mza ,iw) = scp(mza,iw)

    do k = kb, mza-1
       scp_upw(k,iw) = scp0(k,iw)
       if (wmasc(k,iw) < 0.) scp_upw(k,iw) = scp0(k+1,iw)
    enddo

    ! Higher-order vertical fluxes

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

    sfluxwh(kb-1) = 0.
    sfluxwh(mza)  = 0.

    do k = kb, mza-1
       scp_hiw(k) = scpt(k)
       if (wmasc(k,iw) < 0.) scp_hiw(k) = scpb(k+1)

       scp_hiw(k) = scp_hiw(k) - scp_upw(k,iw)
       sfluxwh(k) = fexpl * wmasc(k,iw) * scp_hiw(k)
    enddo

    ! Scalar horizontal flux divergence

    fluxdiv_low(:) = 0.
    fluxdiv_out(:) = 0.

    do jv = 1, itab_w(iw)%npoly
       iv = itab_w(iw)%iv(jv)

       do k = lpv(iv), mza
          fluxdiv_low(k) = fluxdiv_low(k) &
                         + itab_w(iw)%dirv(jv) * vmasc(k,iv) * (scp_upv(k,iv) - scp0(k,iw))

          fluxdiv_out(k) = fluxdiv_out(k) + min( itab_w(iw)%dirv(jv) * sfluxvh(k,iv), 0.)
       enddo
    enddo

    ! Scalar vertical flux divergence

    if (fimpl > .001) then

       do k = kb, mza
          fluxdiv_low(k) = fluxdiv_low(k) + fexpl * ( wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw))   &
                                                    - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw)) ) &

                                          + fimpl * scp0(k,iw) * (wmasc(k,iw) - wmasc(k-1,iw))

          fluxdiv_out(k) = fluxdiv_out(k) + min(sfluxwh(k-1),0.) - max(sfluxwh(k),0.)
       enddo

    else

       do k = kb, mza
          fluxdiv_low(k) = fluxdiv_low(k) + wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw)) &
                                          - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw))

          fluxdiv_out(k) = fluxdiv_out(k) + min(sfluxwh(k-1),0.) - max(sfluxwh(k),0.)
       enddo

    endif

    ! Flux limiter / renormalization

    do k = kb, mza
       rscp_low = scp0(k,iw) * rhos(k,iw) + dtr * volti(k,iw) * fluxdiv_low(k)

       flux_out = max( mdtrp * volti(k,iw) * fluxdiv_out(k), eps0 )

       sc_out = max(onem * rscp_low - eps0, 0.)

       scale_out(k,iw) = min(1., sc_out / max(flux_out, 0.01*sc_out) )
    enddo

    ! Apply flux renormalization to vertical fluxes

    do k = kb, mza-1
       scale = scale_out(k,iw)
       if (sfluxwh(k) < 0.) scale = scale_out(k+1,iw)
       scp_upw(k,iw) = scp_upw(k,iw) + scp_hiw(k) * scale
    enddo

  end subroutine scalar_wflux_and_scale_pd


  subroutine scalar_tend_pd( iw, fimpl )

    use mem_grid,   only: mza, lpw, lpv, volti
    use mem_ijtabs, only: itab_w
!   import,         only: scp, scp0, wmasc, vmasc, scp_upv, scp_upw, rhos, dtr

    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: fimpl
    integer             :: kb, jv, iv, k
    real                :: fluxdiv(mza), fexpl

    kb = lpw(iw)
    fluxdiv(:) = 0.0

    ! Scalar horizontal flux divergence

    do jv = 1, itab_w(iw)%npoly
       iv = itab_w(iw)%iv(jv)

       do k = lpv(iv), mza
          fluxdiv(k) = fluxdiv(k) &
                     + itab_w(iw)%dirv(jv) * vmasc(k,iv) * (scp_upv(k,iv) - scp0(k,iw))
       enddo
    enddo

    if (fimpl > .001) then
       fexpl = 1. - fimpl

       do k = kb, mza
          ! Scalar vertical flux divergence
          fluxdiv(k) = fluxdiv(k) + fexpl * ( wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw)) &
                                            - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw)) )
          ! Integrate forward in time
          scp(k,iw) = scp0(k,iw) + dtr * volti(k,iw) * fluxdiv(k) / rhos(k,iw)
       enddo

    else

       do k = kb, mza
          ! Scalar vertical flux divergence
          fluxdiv(k) = fluxdiv(k) + wmasc(k-1,iw) * (scp_upw(k-1,iw) - scp0(k,iw)) &
                                  - wmasc(k  ,iw) * (scp_upw(k  ,iw) - scp0(k,iw))

          ! Integrate forward in time
          scp(k,iw) = scp0(k,iw) + dtr * volti(k,iw) * fluxdiv(k) / rhos(k,iw)
       enddo

       scp(2:kb-1,iw) = scp(kb,iw)
    endif

  end subroutine scalar_tend_pd


end subroutine update_scalar_pd

!===============================================================================

  subroutine check_cfls(vmasc, wmasc, rho, nrk, ff_implic)

    ! Diagnose outflow CFL number; also used for advection CFL stability check

    use mem_ijtabs,  only: jtab_w, jtw_prog, itab_w
    use mem_grid,    only: lpv, lpw, volti, mza, mva, mwa
    use misc_coms,   only: iparallel, dtlm
    use mem_para,    only: myrank, mgroupsize
    use mem_basic,   only: wc
!   use olam_mpi_atm,only: mpi_post_direct_recv_w, mpi_post_direct_send_w, &
!                          mpi_finish_direct_recv_w, mpi_finish_direct_send_w
!$  use omp_lib,     only: omp_get_max_threads, omp_get_thread_num

#ifdef OLAM_MPI
    use mpi_f08
#endif

    implicit none

    integer,        intent(in ) :: nrk
    real,           intent(in ) :: vmasc(mza,mva)
    real,           intent(in ) :: wmasc(mza,mwa)
    real,           intent(in ) :: rho  (mza,mwa)
    real,           intent(out) :: ff_implic(nrk,mwa)

    integer :: j, iv, k, iw, n
    integer :: ier, inode, iwnode, ihnode
    integer :: mlocc, mlocw, mloch, nthreads, myid
    real    :: dt, fimp(nrk)
    real    :: cfl(mza), cflh(mza), dovr, cflw
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

    dt = dtlm

    !$omp parallel private(myid,cfl,cflh,ct,itk,itw,ch,ihk,ihw,wm,iwk,iww,fimp)
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

    !$omp do private(iw,k,n,iv,dovr,cflw)
    do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

       cflh = 0.
       fimp = 0.

       do n = 1, itab_w(iw)%npoly
          iv = itab_w(iw)%iv(n)
          do k = lpv(iv), mza
             cflh(k) = cflh(k) - min( itab_w(iw)%dirv(n) * vmasc(k,iv), 0.0)
          enddo
       enddo

       do k = lpw(iw), mza
          dovr = dt * volti(k,iw) / rho(k,iw)
          cflh(k) = cflh(k) * dovr
          cflw    = (max(wmasc(k,iw),0.) - min(wmasc(k-1,iw),0.)) * dovr
          cfl (k) = cflh(k) + cflw

          cflw = max(cflw, .1)

          fimp(nrk)   = max(fimp(nrk),  (cfl(k) - 0.98) / cflw)
          fimp(nrk-1) = max(fimp(nrk-1),(cfl(k) - 1.98) / cflw)
          if (nrk == 3) &
              fimp(1) = max(fimp(1),    (cfl(k) - 2.98) / cflw)
       enddo

       do n = 1, nrk
          ff_implic(n,iw) = min(1.0, fimp(n))
          if (ff_implic(n,iw) < .02) ff_implic(n,iw) = 0.
       enddo

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

end module scalar_transport
