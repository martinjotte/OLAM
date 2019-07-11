subroutine divh_damp(mrl)

  use mem_ijtabs,   only: jtab_v, jtab_w, itab_v, itab_w, jtv_wadj, &
                          jtv_prog, jtw_prog
  use mem_grid,     only: lpv, mza, mva, mwa, vnxo2, vnyo2, vnzo2, dnu, &
                          arw0, arw0i, dniv
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_basic,    only: vc, vxe, vye, vze, rho
  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_w
  use misc_coms,    only: iparallel, rinit
  use oname_coms,   only: nl
  use mem_rayf,     only: dorayfdiv, krayfdiv_bot, rayf_cofdiv

  implicit none

  integer, intent(in)     :: mrl
  integer                 :: iv, iw, iw1, iw2, k, j, jv, ktop
  real                    :: dnufac
  real                    :: vc_ec   (mza,mva)
  real                    :: div2d_ex(mza,mwa)
  logical,           save :: firstime = .true.
  real, allocatable, save :: aoc(:), c1(:), c2(:)
  real                    :: cd, dnusum, f1, f2
  real, parameter         :: onethird = 1./3.

  if (nl%akmin_div2d < 1.e-7 .and. .not. dorayfdiv) return

  if (firstime) then
     firstime = .false.

     allocate(aoc(mwa)) ; aoc = rinit
     allocate(c1(mva))
     allocate(c2(mva))

     !$omp parallel do private(iw,dnusum,jv,iv)
     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
        dnusum = 0.

        do jv = 1, itab_w(iw)%npoly
           iv = itab_w(iw)%iv(jv)
           dnusum = dnusum + dnu(iv)
        enddo

        ! Area-over-circumference factors to convert divergence
        ! to mean velocity of divergence

        aoc(iw) = arw0(iw) / dnusum

     enddo
     !$omp end parallel do

     if (iparallel == 1) then
        call mpi_send_w(1, r1dvara1=aoc)
        call mpi_recv_w(1, r1dvara1=aoc)
     endif
     call lbcopy_w(1, v1=aoc)

     cd = nl%akmin_div2d * 0.075

     !$omp parallel do private(iv,iw1,iw2)
     do j = 1,jtab_v(jtv_prog)%jend(1); iv = jtab_v(jtv_prog)%iv(j)

        iw1 = itab_v(iv)%iw(1)
        iw2 = itab_v(iv)%iw(2)

        c1(iv) =  cd * (arw0(iw1)**onethird)**2 * dniv(iv)
        c2(iv) =  cd * (arw0(iw2)**onethird)**2 * dniv(iv)
     enddo
     !$omp end parallel do

  endif  ! firstime

  !$omp parallel
  !$omp do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

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

        dnufac = dnu(iv) * arw0i(iw) * itab_w(iw)%dirv(jv)

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

  !$omp parallel do private(iv,iw1,iw2,k)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1)
     iw2 = itab_v(iv)%iw(2)

     ! Upper atmosphere layers (rayfdiv)

     if (dorayfdiv) then

        ktop = krayfdiv_bot - 1

        do k = krayfdiv_bot, mza

           f1 = max( rayf_cofdiv(k) * aoc(iw1), c1(iv) )
           f2 = max( rayf_cofdiv(k) * aoc(iw2), c2(iv) )

           vmt(k,iv) = vmt(k,iv) + real(rho(k,iw2)) * f2 * div2d_ex(k,iw2) &
                                 - real(rho(k,iw1)) * f1 * div2d_ex(k,iw1)
        enddo

     else

        ktop = mza

     endif

     ! All atmosphere layers (new)

     do k = lpv(iv), ktop
        vmt(k,iv) = vmt(k,iv) + real(rho(k,iw2)) * c2(iv) * div2d_ex(k,iw2) &
                              - real(rho(k,iw1)) * c1(iv) * div2d_ex(k,iw1)
     enddo

  enddo
  !$omp end parallel do

end subroutine divh_damp




subroutine vort_damp(mrl)

  use mem_ijtabs,   only: jtab_v, jtab_m, itab_v, itab_m, jtv_prog, jtm_vadj
  use mem_grid,     only: mva, lpv, lpm, mza, mma, dnu, dnv, arm0, zfacit
  use olam_mpi_atm, only: mpi_send_m, mpi_recv_m
  use mem_basic,    only: vc, rho
  use mem_tend,     only: vmt
  use obnd,         only: lbcopy_m
  use misc_coms,    only: iparallel
  use oname_coms,   only: nl

  implicit none

  integer, intent(in)     :: mrl
  integer                 :: iv, iw1, iw2, k, j, jv, im, kb
  integer                 :: iv1, iv2, iv3, iv4
  integer                 :: im1, im2, im3, im4, im5, im6
  real                    :: vortp(mza,mma)
  real                    :: arm0i, vort_big1, vort_big2, cv
  real, allocatable, save :: c1(:), c2(:)
  logical,           save :: firstime = .true.
  real, parameter         :: onethird = 1./3.

  if (nl%akmin_vort < 1.e-7) return

  if (firstime) then
     firstime = .false.

     allocate(c1(mva))
     allocate(c2(mva))

     cv = nl%akmin_vort  * 0.075

     !$omp parallel do private(iv,iv1,iv2,iv3,iv4,im1,im2,iw1,iw2)
     do j = 1,jtab_v(jtv_prog)%jend(1); iv = jtab_v(jtv_prog)%iv(j)

        iv1  = itab_v(iv)%iv(1)
        iv2  = itab_v(iv)%iv(2)
        iv3  = itab_v(iv)%iv(3)
        iv4  = itab_v(iv)%iv(4)

        im1  = itab_v(iv)%im(1)
        im2  = itab_v(iv)%im(2)

        iw1  = itab_v(iv)%iw(1)
        iw2  = itab_v(iv)%iw(2)

        c1(iv) = -cv * (arm0(im1)**onethird)**2 * dnu(iv1) * dnu(iv2) / &
                  ( dnu(iv1) * dnu(iv2) * dnv(iv)  &
                  + dnu(iv)  * dnu(iv2) * dnv(iv1) &
                  + dnu(iv)  * dnu(iv1) * dnv(iv2) )

        c2(iv) =  cv * (arm0(im2)**onethird)**2 * dnu(iv3) * dnu(iv4) / &
                  ( dnu(iv3) * dnu(iv4) * dnv(iv)  &
                  + dnu(iv)  * dnu(iv4) * dnv(iv3) &
                  + dnu(iv)  * dnu(iv3) * dnv(iv4) )

     enddo
     !$omp end parallel do

  endif  ! firstime

  !$omp parallel do private(im,kb,jv,iv,k,arm0i)
  do j = 1,jtab_m(jtm_vadj)%jend(mrl); im = jtab_m(jtm_vadj)%im(j)

     kb = lpm(im)

     vortp(:,im) = 0.

     ! Loop over V neighbors to evaluate circulation around M (at time T)
     do jv = 1, 3
        iv = itab_m(im)%iv(jv)

        if (itab_v(iv)%im(2) == im) then

           do k = kb,mza
              vortp(k,im) = vortp(k,im) + vc(k,iv) * dnv(iv)
           enddo

        else

           do k = kb,mza
              vortp(k,im) = vortp(k,im) - vc(k,iv) * dnv(iv)
           enddo

        endif
     enddo

     ! Convert circulation to relative vertical vorticity at M 
     ! (DNV lacks the zfact factor and ARM0 lacks the zfact**2 factor, so we
     ! divide their quotient by zfact)

     arm0i = 1. / arm0(im)

     do k = kb,mza
        vortp(k,im) = vortp(k,im) * arm0i * zfacit(k)
     enddo

  enddo
  !$omp end parallel do 

  if (iparallel == 1) then
     call mpi_send_m(mrl, rvara1=vortp)
     call mpi_recv_m(mrl, rvara1=vortp)
  endif
  call lbcopy_m(mrl, a1=vortp)

  !$omp parallel do private(iv,iw1,iw2,im1,im2,im3,im4,im5,im6,k,vort_big1,vort_big2)
  do j = 1,jtab_v(jtv_prog)%jend(mrl); iv = jtab_v(jtv_prog)%iv(j)

     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
     im1 = itab_v(iv)%im(1); im2 = itab_v(iv)%im(2); im3 = itab_v(iv)%im(3)
     im4 = itab_v(iv)%im(4); im5 = itab_v(iv)%im(5); im6 = itab_v(iv)%im(6)

     do k = lpv(iv), mza

        vort_big1 = onethird * (vortp(k,im2) + vortp(k,im3) + vortp(k,im4))
        vort_big2 = onethird * (vortp(k,im1) + vortp(k,im5) + vortp(k,im6))

        ! Horizontal filter for vertical vorticity

        vmt(k,iv) = vmt(k,iv) + real( rho(k,iw1) + rho(k,iw2) ) &
                  * ( (vort_big1 - vortp(k,im1)) * c1(iv) &
                    + (vort_big2 - vortp(k,im2)) * c2(iv) )
     enddo

  enddo
  !$omp end parallel do

end subroutine vort_damp
