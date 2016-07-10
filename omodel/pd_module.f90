module pdtrans
  implicit none

  real, allocatable :: beta(:,:)
  real, allocatable :: dtom(:,:)


contains


  subroutine alloc_pdtrans(imonot, mza, mwa)
    use misc_coms, only: rinit
    implicit none
    
    integer, intent(in) :: imonot, mza, mwa

    if (imonot == 2) then
       allocate(beta(mza,mwa))
       allocate(dtom(mza,mwa))
    endif

  end subroutine alloc_pdtrans



  subroutine comp_and_apply_pd_lims(scp, scp_upv, scp_upw, iwdepv, kdepw, vmsca, wmsca, mrl)
    use mem_grid,     only: mwa, mva, mza, lpv, lpw
    use mem_ijtabs,   only: jtab_v, jtv_wadj, jtab_w, jtw_prog, itab_w
    use misc_coms,    only: iparallel
    use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
    use obnd,         only: lbcopy_w

    implicit none

    real,    intent(in)    :: scp    (mza,mwa)
    real,    intent(inout) :: scp_upv(mza,mva)
    real,    intent(inout) :: scp_upw(mza,mwa)
    integer, intent(in)    :: iwdepv (mza,mva)
    integer, intent(in)    :: kdepw  (mza,mwa)
    real,    intent(in)    :: vmsca  (mza,mva)
    real,    intent(in)    :: wmsca  (mza,mwa)
    integer, intent(in)    :: mrl
    integer                :: iw, j, iv, jv, k, iwd, kd, km, kp

    !$omp parallel

    !$omp do private(iw,k,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       beta(1:lpw(iw)-1,iw) = 0.0

       do k = lpw(iw), mza
          kp = min(k+1,mza)
          beta(k,iw) = max(wmsca(k,  iw),0.0) * scp_upw(k,  iw)  &
                     - min(wmsca(k-1,iw),0.0) * scp_upw(k-1,iw)
       enddo

       do jv = 1, itab_w(iw)%npoly
          iv  = itab_w(iw)%iv(jv)

          do k = lpv(iv), mza
             kp = min(k+1,mza)
             beta(k,iw) = beta(k,iw) - min(itab_w(iw)%dirv(jv) * vmsca(k,iv), 0.0) * scp_upv(k,iv)
          enddo

       enddo

       do k = lpw(iw), mza
          beta(k,iw) = 0.999998 * scp(k,iw) / (dtom(k,iw) * max(beta(k,iw),1.e-15))
          beta(k,iw) = min(1.0, max(beta(k,iw), 0.0))
       enddo

    enddo
    !$omp end do

    !$omp single
    if (iparallel == 1) call mpi_send_w(mrl, rvara1=beta)
    !$omp end single

    !$omp do private(iw,k,kd)
    do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! Loop over W/M level
       do k = lpw(iw), mza-1
          kd = kdepw(k,iw)
          scp_upw(k,iw) = scp_upw(k,iw) * beta(kd,iw)
       enddo

    enddo
    !$omp end do

    !$omp single
    if (iparallel == 1) call mpi_recv_w(mrl, rvara1=beta)
    call lbcopy_w(mrl, a1=beta)
    !$omp end single

    !$omp do private(iv,k,iwd)
    do j = 1, jtab_v(jtv_wadj)%jend(mrl); iv = jtab_v(jtv_wadj)%iv(j)

       ! Loop over T/V level
       do k = lpv(iv), mza
          iwd = iwdepv(k,iv)
          scp_upv(k,iv) = scp_upv(k,iv) * beta(k,iwd)
       enddo

    enddo
    !$omp end do nowait

    !$omp end parallel
  end subroutine comp_and_apply_pd_lims

end module pdtrans
