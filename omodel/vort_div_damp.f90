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

     allocate(fdiv1(mza,mva)) ; fdiv1 = 0.
     fact3 = max(nl%divh_damp_fact,.05) / real(dtsm(1))

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

     !$omp parallel do private(iv,iw1,iw2,iw3,iw4,k)
     do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)
        iw3 = itab_v(iv)%iw(3)
        iw4 = itab_v(iv)%iw(4)

        kmaxv(iv) = max(kmax2(iw1), kmax2(iw2), kmax2(iw3), kmax2(iw4))

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

              ! assumes no blockage by terrain in top Rayleigh friction layer
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



!!subroutine vort_damp(mrl)
!!
!!  use mem_ijtabs,   only: jtab_v, itab_v, jtv_prog, jtab_m, itab_m, jtm_vadj
!!  use mem_grid,     only: lpv, mza, vnxo2, vnyo2, vnzo2, dniv, dnv, mma, lpm, &
!!                          lpw, lve2
!!  use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
!!  use mem_basic,    only: vc, vxe, vye, vze, rho
!!  use mem_tend,     only: vmt
!!  use obnd,         only: lbcopy_m
!!  use misc_coms,    only: iparallel, dtlm
!!  use oname_coms,   only: nl
!!
!!  implicit none
!!
!!  integer, intent(in) :: mrl
!!  integer             :: iv, iw1, iw2, k, j, jv, kb, im, im1, im2, k2
!!  real                :: facv, facs, vce, dti, dv
!!  real                :: vortp_ex(mza,mma)
!!
!!  if (nl%vort_damp_fact < 1.e-8) return
!!
!!  !$omp parallel do private(im,kb,jv,iv,iw1,iw2,dv,k)
!!  do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)
!!
!!     kb = lpm(im)
!!
!!     vortp_ex(:,im) = 0.
!!
!!     ! Loop over V neighbors to evaluate circulation around M (at time T)
!!     do jv = 1, 3
!!        iv = itab_m(im)%iv(jv)
!!
!!        iw1 = itab_v(iv)%iw(1)
!!        iw2 = itab_v(iv)%iw(2)
!!
!!        if (itab_v(iv)%im(2) == im) then
!!           dv = +dnv(iv)
!!        else
!!           dv = -dnv(iv)
!!        endif
!!
!!        do k = kb, mza
!!           vortp_ex(k,im) = vortp_ex(k,im) + dv * ( vc(k,iv)      &
!!                          - vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
!!                          - vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
!!                          - vnzo2(iv) * (vze(k,iw1) + vze(k,iw2)) )
!!        enddo
!!
!!     enddo
!!  enddo
!!  !$omp end parallel do
!!
!!  ! MPI send/recv of vortp
!!  if (iparallel == 1) then
!!     call mpi_send_m(mrl, rvara1=vortp_ex)
!!     call mpi_recv_m(mrl, rvara1=vortp_ex)
!!  endif
!!  call lbcopy_m(mrl, a1=vortp_ex)
!!
!!  !$omp parallel do private(iv,iw1,iw2,im1,im2,k2,dti,facv,facs,k,vce)
!!  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
!!
!!     iw1 = itab_v(iv)%iw(1)
!!     iw2 = itab_v(iv)%iw(2)
!!     im1 = itab_v(iv)%im(1)
!!     im2 = itab_v(iv)%im(2)
!!
!!     dti = 1.0 / real(dtlm(itab_v(iv)%mrlv))
!!
!!     facv = 0.5 * nl%vort_damp_fact * dniv(iv) * dti
!!     do k = lpv(iv), mza
!!        vmt(k,iv) = vmt(k,iv) + facv * real(rho(k,iw1) + rho(k,iw2)) &
!!                              * (vortp_ex(k,im1) - vortp_ex(k,im2))
!!     enddo
!!
!!  enddo
!!  !$omp end parallel do
!!
!!end subroutine vort_damp
!!
!!
!!
!!
!!subroutine rayf_vels(mrl)
!!
!!  use mem_basic,  only: vxe, vye, vze, vc, rho
!!  use mem_tend,   only: vmt
!!  use mem_rayf,   only: dorayf, rayf_cof, krayf_bot
!!  use mem_ijtabs, only: jtab_v, itab_v, jtv_prog
!!  use mem_grid,   only: mza, vcn_ew, vcn_ns, vnx, vny, vnz
!!  use misc_coms,  only: initial, u01d, v01d
!!
!!  implicit none
!!
!!  integer, intent(in) :: mrl
!!  integer :: j, iv, iw1, iw2, iw3, iw4, k
!!  real    :: vcbar(mza)
!!
!!  if (.not. dorayf) return
!!
!!  !$omp parallel private(vcbar)
!!  !$omp do private(iv,iw1,iw2,iw3,iw4,k)
!!  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
!!
!!     iw1 = itab_v(iv)%iw(1)
!!     iw2 = itab_v(iv)%iw(2)
!!     iw3 = itab_v(iv)%iw(3)
!!     iw4 = itab_v(iv)%iw(4)
!!
!!     ! For initial = 2 or 3, just smooth V to a local average (similiar to
!!     ! damping the Laplacian). Todo: for latitudinally homogeneous
!!     ! initialization, smooth to the zonally-averaged fields valid at that
!!     ! date and time
!!
!!     if (initial == 1) then
!!
!!        do k = krayf_bot, mza
!!           vcbar(k) = u01d(k) * vcn_ew(k) + v01d(k) * vcn_ns(k)
!!        enddo
!!
!!     else
!!
!!        do k = krayf_bot, mza
!!           vcbar(k) = 0.25 * ( vnx(iv) * (vxe(k,iw1) + vxe(k,iw2) + vxe(k,iw3) + vxe(k,iw4)) &
!!                             + vny(iv) * (vye(k,iw1) + vye(k,iw2) + vye(k,iw3) + vye(k,iw4)) &
!!                             + vnz(iv) * (vze(k,iw1) + vze(k,iw2) + vze(k,iw3) + vze(k,iw4)) )
!!        enddo
!!
!!     endif
!!
!!     do k = krayf_bot, mza
!!        vmt(k,iv) = vmt(k,iv) + 0.5 * rayf_cof(k) * real(rho(k,iw1) + rho(k,iw2)) &
!!                              * (vcbar(k) - vc(k,iv))
!!     enddo
!!
!!     ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions
!!     !
!!     ! fracx = abs(xev(iv)) / (real(nxp-1) * .866 * deltax)
!!     ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
!!     ! rayfx = 0.   ! Default: no extra RAYF
!!     ! do k = krayf_bot, mza
!!     !    vmt(k,iv) = vmt(k,iv) &
!!     !              + max(rayf_cof(k),rayfx) * dn03d(k) *  (vc03d(k,iv) - vc(k,iv))
!!     ! enddo
!!
!!  enddo
!!  !$omp end do
!!  !$omp end parallel
!!
!!end subroutine rayf_vels
