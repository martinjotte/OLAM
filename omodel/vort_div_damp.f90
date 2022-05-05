subroutine divh_damp(mrl, vmt_short)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, &
                          jtv_prog, jtw_prog
  use mem_grid!,     only: lpv, mza, mva, mwa, arv, volt, volti, &
              !            lpw, dniu, dzit, dniv, zm, zt, zfacit, arw0
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
  use mem_basic!,    only: vmc, vmp, wc
!  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_w, lbcopy_m
  use misc_coms,    only: iparallel, dtsm
  use oname_coms,   only: nl
  use mem_rayf,     only: dorayfdiv, krayfdiv_bot, rayfdiv_zmin, rayfdiv_expon
  use grad_lib,     only: laplacian2d

  implicit none

  integer, intent(in)     :: mrl
  real,    intent(inout)  :: vmt_short(mza,mva)

  integer                 :: iv, iw, iw1, iw2, k, j, jv, iwn, iw3, iw4
  real                    :: fact1, fact2, fact3
  real                    :: div2d(mza,mwa)
  real                    :: del2d(mza,mwa)

  real, allocatable       :: zf   (:)
  real, allocatable, save :: fdiv1(:,:)
  real, allocatable, save :: fdiv2(:,:)

  integer, allocatable       :: kmax1(:), kmax2(:)
  integer, allocatable, save :: kmaxv(:)


  logical,           save :: firstime = .true.
  logical                 :: dodivdamp

  dodivdamp = nl%divh_damp_fact > 1.e-7

  if ((.not. dorayfdiv) .and. (.not. dodivdamp)) return

  if (firstime) then
     firstime = .false.

     allocate(kmax1(mwa))
     allocate(kmax2(mwa))
     allocate(kmaxv(mva))

     if (dodivdamp) then
        allocate(fdiv2(mza,mva)) ; fdiv2 = 0.
        fact2 = nl%divh_damp_fact**2 / real(dtsm(1))
     endif

     if (dodivdamp .or. dorayfdiv) then
        allocate(fdiv1(mza,mva)) ; fdiv1 = 0.
        fact3 = nl%divh_damp_fact / real(dtsm(1))
     endif

     if (dorayfdiv) then
        allocate(zf(krayfdiv_bot:mza))
        do k = krayfdiv_bot, mza
           zf(k) = ((zt(k) - rayfdiv_zmin) / (zm(mza) - rayfdiv_zmin)) ** rayfdiv_expon
        enddo
        fact1 = nl%rayfdiv_fact / real(dtsm(1))
     endif

     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
        kmax1(iw) = lpw(iw) + max(lsw(iw)-1,lve2(iw))

        do jv = 1, itab_w(iw)%npoly
           iwn = itab_w(iw)%iw(jv)
           kmax1(iw) = max( kmax1(iw), lpw(iwn) + max(lsw(iwn)-1,lve2(iwn)) )
        enddo
     enddo

     if (iparallel == 1) then
        call mpi_send_w(mrl, i1dvara1=kmax1)
        call mpi_recv_w(mrl, i1dvara1=kmax1)
     endif
     call lbcopy_w(mrl, iv1=kmax1)

     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
        kmax2(iw) = kmax1(iw)

        do jv = 1, itab_w(iw)%npoly
           iwn = itab_w(iw)%iw(jv)
           kmax2(iw) = max(kmax2(iw), kmax1(iw))
        enddo
     enddo

     if (iparallel == 1) then
        call mpi_send_w(mrl, i1dvara1=kmax2)
        call mpi_recv_w(mrl, i1dvara1=kmax2)
     endif
     call lbcopy_w(mrl, iv1=kmax2)

     !$omp parallel do private(iv,iw1,iw2,k)
     do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)
!       iw3 = itab_v(iv)%iw(3)
!       iw4 = itab_v(iv)%iw(4)

        kmaxv(iv) = max(kmax2(iw1), kmax2(iw2)) !, kmax2(iw3), kmax2(iw4))

        do k = lpv(iv), kmaxv(iv)
           fdiv1(k,iv) = fact3 * min(volt(k,iw1),volt(k,iw2)) * dzit(k) * dniv(iv) * zfacit(k)
        enddo

        if (dodivdamp) then
           do k = lpv(iv), mza
              fdiv2(k,iv) = -fact2 * min(volt(k,iw1),volt(k,iw2)) * dzit(k) * zfacit(k) &
                          * dniv(iv) * min(arw0(iw1),arw0(iw2))
           enddo
        endif

        if (dorayfdiv) then
           do k = max(krayfdiv_bot,kmaxv(iv)+1), mza
              if (dodivdamp) fdiv2(k,iv) = fdiv2(k,iv) * (1.0 - zf(k))

              ! assumes no blockage by terrain in Rayleigh friction layer
              fdiv1(k,iv) = fact1 * min(arw0(iw1),arw0(iw2)) * dniv(iv) * zf(k)
           enddo
        endif

     enddo
     !$omp end parallel do

     deallocate(kmax1,kmax2)

  endif

  !$omp parallel do private(iw,jv,iv,k)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     div2d(:,iw) = 0.0

     do jv = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(jv)

        do k = lpv(iv), mza
           div2d(k,iw) = div2d(k,iw) - arv(k,iv) * itab_w(iw)%dirv(jv) * vmc(k,iv)
        enddo
     enddo

     do k = lpw(iw), mza
        div2d(k,iw) = div2d(k,iw) * volti(k,iw)
     enddo

  enddo
  !$omp end parallel do

  ! MPI send/recv of div2d
  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=div2d)
     call mpi_recv_w(mrl, rvara1=div2d)
  endif
  call lbcopy_w(mrl, a1=div2d)

  if (dodivdamp) then

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
        call laplacian2d(iw, div2d, del2d(:,iw))
     enddo
     !$end parallel do

     ! MPI send/recv of del2d
     if (iparallel == 1) then
        call mpi_send_w(mrl, rvara1=del2d)
        call mpi_recv_w(mrl, rvara1=del2d)
     endif
     call lbcopy_w(mrl, a1=del2d)

  endif

  !$omp parallel do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     if (dodivdamp) then
        do k = kmaxv(iv)+1, mza
           vmt_short(k,iv) = vmt_short(k,iv) + fdiv2(k,iv) * (del2d(k,iw2) - del2d(k,iw1))
        enddo
     endif

     do k = lpv(iv), kmaxv(iv)
        vmt_short(k,iv) = vmt_short(k,iv) + fdiv1(k,iv) * (div2d(k,iw2) - div2d(k,iw1))
     enddo

     if (dorayfdiv) then
        do k = max(krayfdiv_bot, kmaxv(iv)+1), mza
           vmt_short(k,iv) = vmt_short(k,iv) + fdiv1(k,iv) * (div2d(k,iw2) - div2d(k,iw1))
        enddo
     endif

  enddo
  !$omp end parallel do

end subroutine divh_damp



subroutine vort_damp(mrl)

  use mem_ijtabs,   only: jtab_v, jtab_m, itab_v, itab_m, jtv_prog, jtm_prog
  use mem_grid,     only: lpv, lpm, mza, mma, dnu, dnv, arm0, zfacit
  use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
  use mem_basic,    only: vc, rho
  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_m
  use misc_coms,    only: iparallel
  use oname_coms,   only: nl

  implicit none

  integer, intent(in)     :: mrl
  integer                 :: iv, iw1, iw2, k, j, jv, im, kb
  integer                 :: iv1, iv2, iv3
  integer                 :: im1, im2, im3
  real                    :: vort(mza,mma), del_vort(mza,mma)
  real                    :: arm0i, cv, vbar, dnui
  real, allocatable, save :: cm(:)
  logical,           save :: firstime = .true.

  real, parameter         :: onethird = 1./3.
  real, parameter         :: twothird = 2./3.

  if (nl%akmin_vort < 1.e-7) return

  if (firstime) then
     firstime = .false.

     allocate(cm(mma))

     cv = nl%akmin_vort  * 0.075

     !$omp parallel do private(im,iv1,iv2,iv3)
     do j = 1,jtab_m(jtm_prog)%jend(1); im = jtab_m(jtm_prog)%im(j)

        iv1  = itab_m(im)%iv(1)
        iv2  = itab_m(im)%iv(2)
        iv3  = itab_m(im)%iv(3)

        cm(im) = cv * arm0(im)**twothird * dnu(iv1) * dnu(iv2) * dnu(iv3) / &
                      ( dnu(iv1) * dnu(iv2) * dnv(iv3) &
                      + dnu(iv1) * dnv(iv2) * dnu(iv3) &
                      + dnv(iv1) * dnu(iv2) * dnu(iv3) )

     enddo
     !$omp end parallel do

  endif  ! firstime

  !$omp parallel do private(im,kb,jv,iv,k,arm0i)
  do j = 1,jtab_m(jtm_prog)%jend(mrl); im = jtab_m(jtm_prog)%im(j)

     kb = lpm(im)

     vort(:,im) = 0.

     ! Loop over V neighbors to evaluate circulation around M (at time T)
     do jv = 1, 3
        iv = itab_m(im)%iv(jv)

        if (itab_v(iv)%im(2) == im) then

           do k = kb,mza
              vort(k,im) = vort(k,im) + vc(k,iv) * dnv(iv)
           enddo

        else

           do k = kb,mza
              vort(k,im) = vort(k,im) - vc(k,iv) * dnv(iv)
           enddo

        endif
     enddo

     ! Convert circulation to relative vertical vorticity at M
     ! (DNV lacks the zfact factor and ARM0 lacks the zfact**2 factor, so we
     ! divide their quotient by zfact)

     arm0i = 1. / arm0(im)

     do k = kb,mza
        vort(k,im) = vort(k,im) * arm0i * zfacit(k)
     enddo

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_m(mrl, rvara1=vort)
     call mpi_recv_m(mrl, rvara1=vort)
  endif
  call lbcopy_m(mrl, a1=vort)

  !$omp parallel do private(im,im1,im2,im3,k,vbar)
  do j = 1,jtab_m(jtm_prog)%jend(mrl); im = jtab_m(jtm_prog)%im(j)

     im1  = itab_m(im)%im(1)
     im2  = itab_m(im)%im(2)
     im3  = itab_m(im)%im(3)

     del_vort(2:lpm(im)-1,im) = 0.0

     do k = lpm(im), mza
        vbar = onethird * (vort(k,im1) + vort(k,im2) + vort(k,im3))
        del_vort(k,im) = cm(im) * (vbar - vort(k,im))
     enddo

  enddo
  !$omp end parallel do

  if (iparallel == 1) then
     call mpi_send_m(mrl, rvara1=del_vort)
     call mpi_recv_m(mrl, rvara1=del_vort)
  endif
  call lbcopy_m(mrl, a1=del_vort)


  !$omp parallel do private(iv,iw1,iw2,im1,im2,dnui,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
     im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2)

     dnui = 1.0 / dnu(iv)

     ! Horizontal filter for vertical vorticity

     do k = lpv(iv), mza
        vmt(k,iv) = vmt(k,iv) + real( rho(k,iw1) + rho(k,iw2) ) &
                              * (del_vort(k,im2) - del_vort(k,im1)) * dnui
     enddo

  enddo
  !$omp end parallel do

end subroutine vort_damp
