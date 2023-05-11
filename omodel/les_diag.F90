subroutine les_diag()

  use mem_basic,  only: thil,rr_v,wc,press,ue,ve
  use mem_grid,   only: mza,lpw,zm,arw0
  use mem_ijtabs, only: jtab_w, jtw_prog
  use misc_coms,  only: iparallel, naddsc
  use mem_para,   only: myrank
  use consts_coms,only: r8
  use mem_addsc

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer :: j, iw, k, i, ierror
  real :: arwo2

  real(8), dimension(:), pointer, contiguous :: w_avg, u_avg, v_avg, t_avg, q_avg, p_avg, a_avg
  real(8), dimension(:), pointer, contiguous :: ww_avg, uu_avg, vv_avg, tt_avg, qq_avg, pp_avg, aa_avg
  real(8), dimension(:), pointer, contiguous :: uw_avg, vw_avg, tw_avg, qw_avg, pw_avg, aw_avg

  real, dimension(mza) :: wprime,uprime,vprime,tprime,qprime,pprime,aprime

  real(8), target :: buffera(mza,7)
  real(8), target :: bufferb(mza,13)

  real(8), allocatable, save :: arwm_tot(:)
  logical, save :: firstime = .true.

  if (firstime) then
     allocate(arwm_tot(mza))

     arwm_tot = 0.0
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        do k = lpw(iw), mza
           arwm_tot(k) = arwm_tot(k) + real(arw0(iw),8)
        enddo
     enddo

#ifdef OLAM_MPI
     if (iparallel == 1) call MPI_Allreduce(MPI_IN_PLACE, arwm_tot, mza, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
#endif

     arwm_tot(2:) = 1.d0 / arwm_tot(2:)

     firstime = .false.
  endif

  buffera = 0.0_r8
  bufferb = 0.0_r8

  w_avg => buffera(:,1)
  u_avg => buffera(:,2)
  v_avg => buffera(:,3)
  t_avg => buffera(:,4)
  q_avg => buffera(:,5)
  p_avg => buffera(:,6)
  a_avg => buffera(:,7)

  ww_avg => bufferb(:,1)
  uu_avg => bufferb(:,2)
  vv_avg => bufferb(:,3)
  tt_avg => bufferb(:,4)
  qq_avg => bufferb(:,5)
  pp_avg => bufferb(:,6)
  aa_avg => bufferb(:,7)

  uw_avg => bufferb(:, 8)
  vw_avg => bufferb(:, 9)
  tw_avg => bufferb(:,10)
  qw_avg => bufferb(:,11)
  pw_avg => bufferb(:,12)
  aw_avg => bufferb(:,13)

  !$omp parallel do private(iw,k) reduction(+:buffera)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
     do k = lpw(iw), mza
        buffera(k,1) = buffera(k,1) + arw0(iw) * wc(k,iw)
        buffera(k,2) = buffera(k,2) + arw0(iw) * ue(k,iw)
        buffera(k,3) = buffera(k,3) + arw0(iw) * ve(k,iw)
        buffera(k,4) = buffera(k,4) + arw0(iw) * thil(k,iw)
        buffera(k,5) = buffera(k,5) + arw0(iw) * rr_v(k,iw)
        buffera(k,6) = buffera(k,6) + arw0(iw) * press(k,iw)
        if (naddsc >= 1) &
             buffera(k,7) = buffera(k,7) + arw0(iw) * addsc(1)%sclp(k,iw)
     enddo
  end do
  !$omp end parallel do

  do i = 1, 7
     buffera(2:,i) = arwm_tot(2:) * buffera(2:,i)
  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_AllReduce(MPI_IN_PLACE, buffera, mza*7, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
  endif
#endif

  !$omp parallel private(wprime,uprime,vprime,tprime,qprime,pprime,aprime)
  !$omp do private(iw,k,arwo2) reduction(+:bufferb)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza
        wprime(k) = wc(k,iw) - w_avg(k)
        uprime(k) = ue(k,iw) - u_avg(k)
        vprime(k) = ve(k,iw) - v_avg(k)
        tprime(k) = thil(k,iw) - t_avg(k)
        qprime(k) = rr_v(k,iw) - q_avg(k)
        pprime(k) = press(k,iw) - p_avg(k)
        if (naddsc >= 1) &
             aprime(k) = addsc(1)%sclp(k,iw) - a_avg(k)

        bufferb(k,1) = bufferb(k,1) + arw0(iw) * wprime(k) * wprime(k)
        bufferb(k,2) = bufferb(k,2) + arw0(iw) * uprime(k) * uprime(k)
        bufferb(k,3) = bufferb(k,3) + arw0(iw) * vprime(k) * vprime(k)
        bufferb(k,4) = bufferb(k,4) + arw0(iw) * tprime(k) * tprime(k)
        bufferb(k,5) = bufferb(k,5) + arw0(iw) * qprime(k) * qprime(k)
        bufferb(k,6) = bufferb(k,6) + arw0(iw) * pprime(k) * pprime(k)
        if (naddsc >= 1) &
             bufferb(k,7) = bufferb(k,7) + arw0(iw) * aprime(k) * aprime(k)
     enddo

     arwo2 = 0.5 * arw0(iw)

     do k = lpw(iw), mza-1
        bufferb(k, 8) = bufferb(k, 8) + arwo2 * wprime(k) * (uprime(k+1) + uprime(k))
        bufferb(k, 9) = bufferb(k, 9) + arwo2 * wprime(k) * (vprime(k+1) + vprime(k))
        bufferb(k,10) = bufferb(k,10) + arwo2 * wprime(k) * (tprime(k+1) + tprime(k))
        bufferb(k,11) = bufferb(k,11) + arwo2 * wprime(k) * (qprime(k+1) + qprime(k))
        bufferb(k,12) = bufferb(k,12) + arwo2 * wprime(k) * (pprime(k+1) + pprime(k))
        if (naddsc >= 1) &
             bufferb(k,13) = bufferb(k,13) + arwo2 * wprime(k) * (aprime(k+1) + aprime(k))
     enddo

  enddo
  !$omp end do
  !$omp end parallel

  do i = 1, 13
     bufferb(2:,i) = arwm_tot(2:) * bufferb(2:,i)
  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (myrank == 0) then
        call MPI_Reduce(MPI_IN_PLACE, bufferb, mza*13, MPI_REAL8, MPI_SUM, 0, &
                        MPI_COMM_WORLD, ierror)
     else
        call MPI_Reduce(bufferb, bufferb, mza*13, MPI_REAL8, MPI_SUM, 0, &
                        MPI_COMM_WORLD, ierror)
     endif
  endif
#endif

  if (myrank == 0) then
     write(*,*)
     write(*,'(A8,4A9,4A11)') "Z  ", "<U> ", "<V> ", "<THil>", "<Qv>  ", "<w'w'>  ", &
                              "<u'u'>  ", "<t'w'>  ", "<p'p'>  "

  do k = mza, 2, -1
  write(*,'(F8.2, 2F9.4, F9.3, F9.6, 4G11.4)') zm(k), u_avg(k), v_avg(k), &
       t_avg(k), q_avg(k), ww_avg(k), uu_avg(k), tw_avg(k), pp_avg(k)
  enddo
  write(*,*)
  endif

end subroutine les_diag
