subroutine vort_div_damp(mrl)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, &
                          jtv_prog, jtw_prog, jtab_m, itab_m, jtm_vadj
  use mem_grid,     only: lpv, mza, mva, mwa, vnxo2, vnyo2, vnzo2, dnu, &
                          arw0, arw0i, dniv, dnv, mma, lpm, arm0, dniu
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, mpi_send_m, mpi_recv_m
  use mem_basic,    only: vc, vxe, vye, vze, rho, vmc, vmp
  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_w, lbcopy_m
  use misc_coms,    only: iparallel, rinit, dtlm, io6
  use oname_coms,   only: nl

  implicit none

  integer, intent(in) :: mrl
  integer             :: iv, iw, iw1, iw2, k, j, jv, ktop, kb, im, im1, im2
  real                :: dnufac, facd, facv
  real                :: vc_ec   (mza,mva)
  real                :: div2d_ex(mza,mwa)
  real                :: vortp_ex(mza,mma)
  real                :: cd, dnusum, f1, f2, dtli

  !$omp parallel 
  !$omp do private(iw1,iw2,k) schedule(guided)
  do iv = 2, mva

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     if (iw1 < 2 .or. iw2 < 2) cycle

     do k = lpv(iv), mza
        vc_ec(k,iv) = vc(k,iv) - vnxo2(iv) * (vxe(k,iw1) + vxe(k,iw2)) &
                               - vnyo2(iv) * (vye(k,iw1) + vye(k,iw2)) &
                               - vnzo2(iv) * (vze(k,iw1) + vze(k,iw2))
     enddo
  enddo
  !$omp end do

  !$omp do private(iw,jv,iv,k,dnufac)
  do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

     div2d_ex(:,iw) = 0.0

     do jv = 1, itab_w(iw)%npoly
        iv = itab_w(iw)%iv(jv)

         dnufac = dnu(iv) * itab_w(iw)%dirv(jv)

        do k = lpv(iv), mza
           div2d_ex(k,iw) = div2d_ex(k,iw) - dnufac * vc_ec(k,iv)
        enddo
     enddo

  enddo
  !$omp end do
  !$omp end parallel

  ! MPI send/recv of div2d
  if (iparallel == 1) then
     call mpi_send_w(mrl, rvara1=div2d_ex)
     call mpi_recv_w(mrl, rvara1=div2d_ex)
  endif
  call lbcopy_w(mrl, a1=div2d_ex)

  !$omp parallel do private(im,kb,jv,iv,k)
  do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)

     kb = lpm(im)

     vortp_ex(:,im) = 0.

     ! Loop over V neighbors to evaluate circulation around M (at time T)
     do jv = 1, 3
        iv = itab_m(im)%iv(jv)

        if (itab_v(iv)%im(2) == im) then

           do k = kb,mza
              vortp_ex(k,im) = vortp_ex(k,im) + dnv(iv) * vc_ec(k,iv)
           enddo

        else

           do k = kb,mza
              vortp_ex(k,im) = vortp_ex(k,im) - dnv(iv) * vc_ec(k,iv)
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

  !$omp parallel do private(iv,iw1,iw2,im1,im2,dtli,facd,facv,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     im1 = itab_v(iv)%im(1)
     im2 = itab_v(iv)%im(2)

     dtli = 0.25 / real(dtlm(itab_v(iv)%mrlv))

     facd = nl%divh_damp_fact * dniv(iv) * dtli
     facv = nl%vort_damp_fact * dniv(iv) * dtli

     do k = lpv(iv), mza
        vmt(k,iv) = vmt(k,iv) + real(rho(k,iw1) + rho(k,iw2)) * ( facd * (div2d_ex(k,iw2) - div2d_ex(k,iw1)) &
                                                                + facv * (vortp_ex(k,im1) - vortp_ex(k,im2)) )
     enddo

  enddo
  !$omp end parallel do

end subroutine vort_div_damp
