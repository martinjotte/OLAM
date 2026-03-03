subroutine divh_damp(vmt, dtm)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, &
                          jtv_prog, jtw_prog
  use mem_grid,     only: lpv, mza, mva, mwa, arv, volt, volti, &
                          lpw, dzit, dniv, zm, zt, zfacit, arw0
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
  use mem_basic,    only: vmc
  use obnd,         only: lbcopy_w, lbcopy_m
  use misc_coms,    only: iparallel
  use oname_coms,   only: nl
  use mem_rayf,     only: dorayfdiv, krayfdiv_bot, rayfdiv_zmin, rayfdiv_expon
  use grad_lib,     only: laplacian2d
  use consts_coms,  only: r8
! use oplot_comsm   only: op

  implicit none

  real,     intent(inout) :: vmt(mza,mva)
  real(r8), intent(inout) :: dtm

  integer                 :: iv, iw, iw1, iw2, k, j, jv, nn, n, ka
  integer                 :: il, ih, ilev
  real                    :: fact2, factr
  real(r8)                :: facth, factb

  real, allocatable, save :: fdivh(:,:)   ! 4th or 6th-order damping factor
  real, allocatable, save :: fdivl(:,:)   ! 2nd-order damping factor
  real, allocatable       :: zf   (:)     ! Rayleigh-friction profile
  real, allocatable, save :: div2d(:,:,:) ! divergence, laplacian of divergence

  logical,           save :: firstime = .true.
  logical                 :: dodivdamp

  dodivdamp = (nl%divh_damp_fact > 1.e-7) .and. (nl%divh_damp_level > 0)

  if ((.not. dorayfdiv) .and. (.not. dodivdamp)) return

  if (firstime) then
     firstime = .false.

     ka = 1
     if (.not. dodivdamp) ka = krayfdiv_bot

     ilev = 1
     if (dodivdamp .and. (nl%divh_damp_level > 1)) ilev = 2

     allocate( div2d(ka:mza,mwa,ilev) )

     !$omp parallel do private(n)
     do iw = 1, mwa
        do n = 1, ilev
           div2d(:,iw,n) = 0.
        enddo
     enddo
     !$omp end parallel do

     ! Rayleigh friction divergence damping profile

     allocate( zf(mza) ) ; zf(:) = 0.0

     if ( dorayfdiv ) then
        do k = krayfdiv_bot, mza
           zf(k) = ((zt(k) - rayfdiv_zmin) / (zm(mza) - rayfdiv_zmin)) ** rayfdiv_expon
        enddo
     endif

     ! Factors for vorticity damping

     if (dorayfdiv .or. (dodivdamp .and. (nl%divh_damp_level == 1))) then
        ka = krayfdiv_bot
        if ( dodivdamp .and. (nl%divh_damp_level == 1) ) ka = 1
        allocate( fdivl(ka:mza,mva) )
     endif

     if (dodivdamp .and. (nl%divh_damp_level > 1)) then
        facth = -(-nl%divh_damp_fact)**nl%divh_damp_level / dtm
        allocate( fdivh(mza,mva) )
     endif

     fact2 = 0.
     if (dodivdamp .and. (nl%divh_damp_level == 1)) fact2 = nl%divh_damp_fact / real(dtm)

     factr = 0.
     if (dorayfdiv .and. dodivdamp) then
        factr = max(nl%divh_damp_fact,nl%rayfdiv_fact) / real(dtm)
     elseif (dorayfdiv) then
        factr = nl%rayfdiv_fact / real(dtm)
     endif

     !$omp parallel do private(iv,iw1,iw2,ka,k,factb)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        if (dorayfdiv .or. (dodivdamp .and. (nl%divh_damp_level == 1))) then
           ka = krayfdiv_bot
           if (dodivdamp .and. (nl%divh_damp_level == 1)) ka = lpv(iv)

           do k = ka, mza
              fdivl(k,iv) = (fact2 * (1. - zf(k)) + factr * zf(k)) &
                          * dniv(iv) * zfacit(k) * dzit(k) * min(volt(k,iw1),volt(k,iw2))
           enddo
        endif

        if (dodivdamp .and. (nl%divh_damp_level > 1)) then
           factb = facth * min(arw0(iw1),arw0(iw2))**(nl%divh_damp_level-1)

           do k = lpv(iv), mza
              fdivh(k,iv) = factb * dniv(iv) * zfacit(k) * dzit(k) * min(volt(k,iw1),volt(k,iw2)) * (1. - zf(k))
           enddo
        endif

     enddo
     !$omp end parallel do

     deallocate(zf)

  endif  ! firstime

  ! Compute horizontal divergence

  !$omp parallel do private(iw,jv,iv,ka,k)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     do jv = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(jv)

        ka = lpv(iv)
        if (.not. dodivdamp) ka = max(ka,krayfdiv_bot)

        do k = ka, mza
           div2d(k,iw,1) = div2d(k,iw,1) - arv(k,iv) * itab_w(iw)%dirv(jv) * vmc(k,iv)
        enddo
     enddo

     ka = lpw(iw)
     if (.not. dodivdamp) ka = max(ka,krayfdiv_bot)

     do k = ka, mza
        div2d(k,iw,1) = div2d(k,iw,1) * volti(k,iw)
     enddo

  enddo
  !$omp end parallel do

  ! MPI send/recv of div2d

  if (iparallel == 1) then
     call mpi_send_w(svara1=div2d(:,:,1))
     call mpi_recv_w(svara1=div2d(:,:,1))
  endif
  call lbcopy_w(s1=div2d(:,:,1))

!!SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - - - -
!!(Example of how to plot "external" field; one not available in module memory)
!
! if (mod(time8+dtlm,op%frqplt) < dtlm .and. istp == nstp) then
!    if (.not. allocated(op%extfld)) allocate (op%extfld(mza,mwa))
!    op%extfld(:,:) = div2d(:,:,1)
!    op%extfldname = 'DIVERG'
!    write(io6,*) "Calling plot_fields..."
!    call plot_fields(1)
!    deallocate (op%extfld)
! endif
!!END SPECIAL PLOT SECTION - - - - - - - - - - - - - - - - - -

  if ( dorayfdiv .or. (dodivdamp .and. nl%divh_damp_level == 1) ) then

     !$omp parallel do private(iv,iw1,iw2,ka,k)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        ka = krayfdiv_bot
        if ( dodivdamp .and. (nl%divh_damp_level == 1) ) ka = lpv(iv)

        ! Damp (remove) horizontal divergence

        do k = ka, mza
           vmt(k,iv) = vmt(k,iv) + fdivl(k,iv) * (div2d(k,iw2,1) - div2d(k,iw1,1))
        enddo

     enddo
     !$omp end parallel do

  endif

  ! We are finished here if not doing higher-order divergence damping

  if ( .not. (dodivdamp .and. (nl%divh_damp_level > 1) ) ) return

  ! Loop to compute laplacian of horizontal divergence

  do nn = 1, nl%divh_damp_level - 1
     il = mod(nn-1,2)+1
     ih = mod(nn  ,2)+1

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

        call laplacian2d(iw, div2d(:,:,il), div2d(:,iw,ih))

     enddo
     !$end parallel do

     if (iparallel == 1) then
        call mpi_send_w(svara1=div2d(:,:,ih))
        call mpi_recv_w(svara1=div2d(:,:,ih))
     endif
     call lbcopy_w(s1=div2d(:,:,ih))

  enddo

  ! Divergence damping momentum tendencies

  !$omp parallel do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! Damp (remove) horizontal divergence
     do k = lpv(iv), mza
        vmt(k,iv) = vmt(k,iv) + fdivh(k,iv) * (div2d(k,iw2,ih) - div2d(k,iw1,ih))
     enddo

  enddo
  !$omp end parallel do

end subroutine divh_damp



subroutine vort_damp(vmt, dtm)

  use mem_ijtabs,   only: itab_m, jtab_m, jtm_prog, itab_v, jtab_v, jtv_prog, itabg_v
  use mem_grid,     only: mma, mva, mza, lpm, lpw, lpv, arm0, dnv, dniu, &
                          vnx, vny, vnz, vnxo2, vnyo2, vnzo2, arw0i, zfacit, volt, dzit, arw0
  use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
  use mem_basic,    only: vc, rho, vxe, vye, vze
  use obnd,         only: lbcopy_m
  use misc_coms,    only: iparallel
  use oname_coms,   only: nl
  use consts_coms,  only: r8
  use mem_para,     only: myrank
  use grad_lib,     only: laplacian_m

  use oplot_coms

  implicit none

  real,        intent(inout) :: vmt(mza,mva)
  real(r8),    intent(inout) :: dtm

  integer                    :: iv, iw, iw1, iw2, iw3, i, k, j, jv, im, ka, nn
  integer                    :: im1, im2
  real                       :: vcm
  real                       :: vxbar(mza), vybar(mza), vzbar(mza)
  real                       :: fn(mza)
  integer                    :: ilev, il, ih

  logical,              save :: firstime = .true.
  real,    allocatable, save :: fvortl(:), fvorth(:)
  integer, allocatable, save :: lpmmin(:)
  real,    allocatable, save :: vfact (:,:)
  real,    allocatable, save :: vort2d(:,:,:)


  if (nl%vort_damp_fact < 1.e-7 .or. nl%vort_damp_level <= 0) return

  if (firstime) then
     firstime = .false.

     allocate( fvorth(mva) )
     if (nl%vort_damp_level > 2) allocate( fvortl(mva) )

     allocate( vfact(3,mma) )
     allocate( lpmmin (mma) )

     if (nl%vort_damp_level == 1) then
        allocate( vort2d(mza,mma,1) )
     else
        allocate( vort2d(mza,mma,2) )
     endif

     !$omp parallel
     !$omp do
     do im = 1, mma
        vort2d(:,im,1) = 0.
        if (nl%vort_damp_level > 1) vort2d(:,im,2) = 0.
     enddo
     !$omp end do nowait

     !$omp do private(im,jv,iv)
     do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

        lpmmin(im) = minval( lpv( itab_m(im)%iv(1:3) ) )

        do jv = 1, 3
           iv = itab_m(im)%iv(jv)

           if (itab_v(iv)%im(2) == im) then
              vfact(jv,im) =  dnv(iv) / arm0(im)
           else
              vfact(jv,im) = -dnv(iv) / arm0(im)
           endif
        enddo

     enddo
     !$omp end do nowait

     !$omp do private(iv,im1,im2)
     do j = 1, jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)
        im1 = itab_v(iv)%im(1)
        im2 = itab_v(iv)%im(2)

        fvorth(iv) = 0.50 * dniu(iv) / real(dtm) &
                   * ( -nl%vort_damp_fact * min(arm0(im1),arm0(im2)) )**nl%vort_damp_level

        if (nl%vort_damp_level > 2) then
           fvortl(iv) = 0.50 * dniu(iv) / real(dtm) &
                      * ( -nl%vort_damp_fact * min(arm0(im1),arm0(im2)) )**(nl%vort_damp_level-1)
        endif

     enddo
     !$omp end do nowait
     !$omp end parallel

  endif  ! firstime

  !$omp parallel private(vxbar,vybar,vzbar,fn)
  !$omp do private(im,jv,iw,k,iv,vcm)
  do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

     if (lpmmin(im) < lpm(im)) then

        vxbar(lpmmin(im):lpm(im)-1) = 0.
        vybar(lpmmin(im):lpm(im)-1) = 0.
        vzbar(lpmmin(im):lpm(im)-1) = 0.
        fn   (lpmmin(im):lpm(im)-1) = 0.

        do jv = 1, 3
           iw = itab_m(im)%iw(jv)
           do k = max(lpmmin(im),lpw(iw)), lpm(im)-1
              fn   (k) = fn   (k) + 1.125
              vxbar(k) = vxbar(k) + vxe(k,iw)
              vybar(k) = vybar(k) + vye(k,iw)
              vzbar(k) = vzbar(k) + vze(k,iw)
           enddo
        enddo

     endif

     vort2d(:,im,1) = 0.

     do jv = 1, 3
        iv = itab_m(im)%iv(jv)

        do k = lpmmin(im), lpv(iv)-1
           vcm = ( vnx(iv) * vxbar(k) &
                 + vny(iv) * vybar(k) &
                 + vnz(iv) * vzbar(k) ) / fn(k)

           vort2d(k,im,1) = vort2d(k,im,1) + vcm * vfact(jv,im)
        enddo

        do k = lpv(iv), mza
           vort2d(k,im,1) = vort2d(k,im,1) + vc(k,iv) * vfact(jv,im)
        enddo

     enddo

  enddo
  !$omp end do nowait
  !$omp end parallel

  if (iparallel == 1) then
     call mpi_send_m( rvara1=vort2d(:,:,1) )
     call mpi_recv_m( rvara1=vort2d(:,:,1) )
  endif
  call lbcopy_m( a1=vort2d(:,:,1) )

  if ( nl%vort_damp_level == 1 ) then

     !$omp parallel do private(iv,iw1,iw2,im1,im2,k)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)

        iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
        im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)

        do k = lpv(iv), mza
           vmt(k,iv) = vmt(k,iv) + fvorth(iv) * (vort2d(k,im2,1) - vort2d(k,im1,1)) * real(rho(k,iw1)+rho(k,iw2))
        enddo

     enddo
     !$omp end parallel do

  else

     do nn = 1, nl%vort_damp_level - 1
        il = mod(nn-1,2)+1
        ih = mod(nn  ,2)+1

        !$omp parallel do private(im)
        do j = 1, jtab_m(jtm_prog)%jend; im = jtab_m(jtm_prog)%im(j)

           call laplacian_m(im, vort2d(:,:,il), vort2d(:,im,ih))

        enddo
        !$omp end parallel do

        if (iparallel == 1) then
           call mpi_send_m( rvara1=vort2d(:,:,ih) )
           call mpi_recv_m( rvara1=vort2d(:,:,ih) )
        endif
        call lbcopy_m( a1=vort2d(:,:,ih) )

     enddo

     !$omp parallel do private(iv,iw1,iw2,im1,im2,ka,k)
     do j = 1,jtab_v(jtv_prog)%jend; iv = jtab_v(jtv_prog)%iv(j)

        iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
        im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)

        if (nl%vort_damp_level == 3) then
           ka = max(lpm(im1),lpm(im2))
        else
           ka = lpv(iv)
        endif

        do k = lpv(iv), ka-1
           vmt(k,iv) = vmt(k,iv) + fvortl(iv) * (vort2d(k,im2,il) - vort2d(k,im1,il)) * real(rho(k,iw1)+rho(k,iw2))
        enddo

        do k = ka, mza
           vmt(k,iv) = vmt(k,iv) + fvorth(iv) * (vort2d(k,im2,ih) - vort2d(k,im1,ih)) * real(rho(k,iw1)+rho(k,iw2))
        enddo

     enddo
     !$omp end parallel do

  endif

end subroutine vort_damp
