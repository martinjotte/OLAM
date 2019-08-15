subroutine divh_damp(mrl)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, &
                          jtv_prog, jtw_prog
  use mem_grid,     only: lpv, mza, mva, mwa, dnu, arv, volt, volti, &
                          arw0, lpw, dniu, dzit
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
  use mem_basic,    only: vmc, vmp
  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_w, lbcopy_m
  use misc_coms,    only: iparallel, dtlm
  use oname_coms,   only: nl
  use mem_rayf,     only: dorayfdiv, rayf_cofdiv, krayfdiv_bot
  use grad_lib,     only: laplacian2d

  implicit none

  integer, intent(in)     :: mrl
  integer                 :: iv, iw, iw1, iw2, k, j, jv, ktop
  real                    :: facd
  real                    :: fact_div(mza)
  real                    :: div2d(mza,mwa)
  real                    :: del2d(mza,mwa)
  logical                 :: firstime = .true.
  real, allocatable       :: asymw(:)
  real, allocatable, save :: asymv(:)
  real, parameter         :: c1 = -1.0 / sqrt(3.)


  if (firstime) then
     firstime = .false.

     allocate(asymw(mwa))
     allocate(asymv(mva))

     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
        asymw(iw) = minval(dnu(itab_w(iw)%iv(1:itab_w(iw)%npoly))) &
                  / maxval(dnu(itab_w(iw)%iv(1:itab_w(iw)%npoly)))
     enddo

     if (iparallel == 1) then
        call mpi_send_w(mrl, r1dvara1=asymw)
        call mpi_recv_w(mrl, r1dvara1=asymw)
     endif
     call lbcopy_w(mrl, v1=asymw)

     do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)
        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)
        asymv(iv) = min(asymw(iw1),asymw(iw2))
     enddo

     deallocate(asymw)
  endif

  !$omp parallel do private(iw,jv,iv,k)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     div2d(:,iw) = 0.0

     do jv = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(jv)

        do k = lpv(iv), mza
           div2d(k,iw) = div2d(k,iw) - arv(k,iv) * itab_w(iw)%dirv(jv) &
                                     * (1.5 * vmc(k,iv) - 0.5 * vmp(k,iv))
        enddo
     enddo

     do k = lpw(iw), mza
        div2d(k,iw) = div2d(k,iw) * real(volti(k,iw))
     enddo

  enddo
  !$omp end parallel do

  ! MPI send/recv of div2d
  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=div2d)
     call mpi_recv_w(mrl, rvara1=div2d)
  endif
  call lbcopy_w(mrl, a1=div2d)

  !$omp parallel do private(iw)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
     call laplacian2d(iw, div2d, del2d(:,iw))
  enddo
  !$end parallel do

  ! MPI send/recv of div2d
  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=del2d)
     call mpi_recv_w(mrl, rvara1=del2d)
  endif
  call lbcopy_w(mrl, a1=del2d)

  ktop = mza
  if (dorayfdiv) ktop = krayfdiv_bot - 1

  !$omp parallel private(fact_div)
  !$omp do private(iv,iw1,iw2,facd,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     facd = c1 * min(arw0(iw1),arw0(iw2)) * dniu(iv) * asymv(iv) &
          / real(dtlm(itab_v(iv)%mrlv))

     do k = lpv(iv), ktop
        fact_div(k) = facd * nl%divh_damp_fact**2
     enddo

     if (dorayfdiv) then
        do k = ktop+1, mza
           fact_div(k) = facd * max(nl%divh_damp_fact,rayf_cofdiv(k))**2
        enddo
     endif

     do k = lpv(iv), mza
        vmt(k,iv) = vmt(k,iv) + fact_div(k) * real(min(volt(k,iw1),volt(k,iw2))) &
                              * dzit(k) * (del2d(k,iw2) - del2d(k,iw1))
     enddo

  enddo
  !$omp end do
  !$omp end parallel

end subroutine divh_damp



subroutine vort_damp(mrl)

  use mem_ijtabs,   only: jtab_v, itab_v, jtv_prog, jtab_m, itab_m, jtm_vadj
  use mem_grid,     only: lpv, mza, vnxo2, vnyo2, vnzo2, dniv, dnv, mma, lpm
  use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
  use mem_basic,    only: vc, vxe, vye, vze, rho
  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_m
  use misc_coms,    only: iparallel, dtlm
  use oname_coms,   only: nl

  implicit none

  integer, intent(in) :: mrl
  integer             :: iv, iw1, iw2, k, j, jv, kb, im, im1, im2
  real                :: facv
  real                :: vortp_ex(mza,mma)

  !$omp parallel do private(im,kb,jv,iv,k,iw1,iw2)
  do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)

     kb = lpm(im)

     vortp_ex(:,im) = 0.

     ! Loop over V neighbors to evaluate circulation around M (at time T)
     do jv = 1, 3
        iv = itab_m(im)%iv(jv)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        if (itab_v(iv)%im(2) == im) then

           do k = kb,mza
              vortp_ex(k,im) = vortp_ex(k,im) + dnv(iv) * ( vc(k,iv) &
                             - vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                             - vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
                             - vnzo2(iv) * (vze(k,iw1) + vze(k,iw2)) )
           enddo

        else

           do k = kb,mza
              vortp_ex(k,im) = vortp_ex(k,im) - dnv(iv) * ( vc(k,iv) &
                             - vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                             - vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
                             - vnzo2(iv) * (vze(k,iw1) + vze(k,iw2)) )
           enddo

        endif
     enddo
  enddo
  !$omp end parallel do

  ! MPI send/recv of vortp
  if (iparallel == 1) then
     call mpi_send_m(mrl, rvara1=vortp_ex)
     call mpi_recv_m(mrl, rvara1=vortp_ex)
  endif
  call lbcopy_m(mrl, a1=vortp_ex)

  !$omp parallel do private(iv,iw1,iw2,im1,im2,facv,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)
     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)

     facv = 0.5 * nl%vort_damp_fact * dniv(iv) / real(dtlm(itab_v(iv)%mrlv))

     do k = lpv(iv), mza
        vmt(k,iv) = vmt(k,iv) + facv * real(rho(k,iw1) + rho(k,iw2)) &
                                     * (vortp_ex(k,im1) - vortp_ex(k,im2))
     enddo

  enddo
  !$omp end parallel do

end subroutine vort_damp




subroutine rayf_vels(mrl)

  use mem_basic,  only: vxe, vye, vze, vc, rho
  use mem_tend,   only: vmt
  use mem_rayf,   only: dorayf, rayf_cof, krayf_bot
  use mem_ijtabs, only: jtab_v, itab_v, jtv_prog
  use mem_grid,   only: mza, vcn_ew, vcn_ns, vnx, vny, vnz, lpv
  use misc_coms,  only: initial, u01d, v01d

  implicit none

  integer, intent(in) :: mrl
  integer :: j, iv, iw1, iw2, iw3, iw4, k
  real    :: vcbar(mza)

  if (.not. dorayf) return

  !$omp parallel private(vcbar)
  !$omp do private(iv,iw1,iw2,iw3,iw4,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)
     iw3 = itab_v(iv)%iw(3)
     iw4 = itab_v(iv)%iw(4)

     ! For initial = 2 or 3, just smooth V to a local average (similiar to
     ! damping the Laplacian). Todo: for latitudinally homogeneous
     ! initialization, smooth to the zonally-averaged fields valid at that
     ! date and time

     if (initial == 1) then

        do k = krayf_bot, mza
           vcbar(k) = u01d(k) * vcn_ew(k) + v01d(k) * vcn_ns(k)
        enddo

     else

        do k = krayf_bot, mza
           vcbar(k) = 0.25 * ( vnx(iv) * (vxe(k,iw1) + vxe(k,iw2) + vxe(k,iw3) + vxe(k,iw4)) &
                             + vny(iv) * (vye(k,iw1) + vye(k,iw2) + vye(k,iw3) + vye(k,iw4)) &
                             + vnz(iv) * (vze(k,iw1) + vze(k,iw2) + vze(k,iw3) + vze(k,iw4)) )
        enddo

     endif

     do k = krayf_bot, mza
        vmt(k,iv) = vmt(k,iv) + 0.5 * rayf_cof(k) * real(rho(k,iw1) + rho(k,iw2)) &
                              * (vcbar(k) - vc(k,iv))
     enddo

     ! Alternate form: Extra RAYF at open ends of channel with cyclic end boundary conditions
     !
     ! fracx = abs(xev(iv)) / (real(nxp-1) * .866 * deltax)
     ! rayfx = .2 * (-2. + 3. * fracx) * rayf_cof(mza)
     ! rayfx = 0.   ! Default: no extra RAYF
     ! do k = lpv(iv), mza
     !    vmt(k,iv) = vmt(k,iv) &
     !              + max(rayf_cof(k),rayfx) * dn03d(k) *  (vc03d(k,iv) - vc(k,iv))
     ! enddo

  enddo
  !$omp end do
  !$omp end parallel

end subroutine rayf_vels
