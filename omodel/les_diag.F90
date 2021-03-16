subroutine les_diag()

  use mem_basic!,  only: rr_w,rho,thil,rr_v,wc,wmc,press,vmc,vc,vxe,vye,vze
  use mem_micro!,  only: rr_c,rr_r
  use mem_grid!,   only: mza,mwa,lpw,volti,mva,lpv,glatw,glonw,glatv,glonv
!  use mem_tend,   only: thilt,rr_wt,vmxet,vmyet,vmzet
  use misc_coms!,  only: io6, iparallel
  use mem_ijtabs!, only: itab_v, itab_w
  use mem_para,   only: myrank
  use mem_cuparm
  use mem_land
  use mem_sea
  use mem_turb
  use mem_tend
  use mem_addsc
  use mem_para,  only: mgroupsize, myrank

#ifdef IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer :: i,k,istop,ierror,j
!  integer, intent(in) :: icall

  integer :: jwl, iwl, kw, np, n, ks
  real :: r4, vss(7), vsx(7)

  integer :: ivmax, ikmax, iw1, iw2, iv, iw
  real    :: dvmax, dv

  integer :: iwmax, ikwmax
  real :: vmax, wmax

  real(8), dimension(:), pointer, contiguous :: w_avg, u_avg, v_avg, t_avg, q_avg, p_avg, a_avg
  real(8), dimension(:), pointer, contiguous :: ww_avg, uu_avg, vv_avg, tt_avg, qq_avg, pp_avg, aa_avg
  real(8), dimension(:), pointer, contiguous :: uw_avg, vw_avg, tw_avg, qw_avg, pw_avg, aw_avg

  real, dimension(mza) :: wprime,uprime,vprime,tprime,qprime,pprime,aprime

  real(8), target :: buffera(mza,7)
  real(8), target :: bufferb(mza,13)

  real(8), allocatable, save :: arwm_tot(:)
  logical, save :: firstime = .true.

! real(8), allocatable, save :: afrac(:,:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!
!  return
  !!!!!!!!!!!!!!!!!!!!!!


  if (firstime) then
     allocate(arwm_tot(mza))

     arwm_tot = 0.0
     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
        do k = lpw(iw), mza
           arwm_tot(k) = arwm_tot(k) + real(arw0(iw),8)
        enddo
     enddo

#ifdef OLAM_MPI
     if (iparallel == 1) call MPI_Allreduce(MPI_IN_PLACE, arwm_tot, mza, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
#endif

     arwm_tot(2:) = 1.d0 / arwm_tot(2:)

!!     allocate(afrac  (mwa))
!!     allocate(afraco2(mwa))
!!     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
!!        do k = lpw(iw), mza
!!        afrac  (k,iw) = real(arw0(iw),8) * arwm_tot
!!        afraco2(iw) = 0.5d0 * afrac(iw)
!!     enddo

     firstime = .false.
  endif

  buffera = 0.0
  bufferb = 0.0

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
  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
     do k = lpw(iw), mza
        buffera(k,1) = buffera(k,1) + arw0(iw) * wc(k,iw)
        buffera(k,2) = buffera(k,2) + arw0(iw) * ue(k,iw)
        buffera(k,3) = buffera(k,3) + arw0(iw) * ve(k,iw)
        buffera(k,4) = buffera(k,4) + arw0(iw) * thil(k,iw)
        buffera(k,5) = buffera(k,5) + arw0(iw) * rr_v(k,iw)
        buffera(k,6) = buffera(k,6) + arw0(iw) * press(k,iw)
!       buffera(k,7) = buffera(k,7) + arw0(iw) * addsc(1)%sclp(k,iw)
     enddo
  end do
  !$omp end parallel do

  do i = 1, 6
     buffera(2:,i) = arwm_tot(2:) * buffera(2:,i)
  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_AllReduce(MPI_IN_PLACE, buffera, mza*7, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierror)
  endif
#endif

  !$omp parallel private(wprime,uprime,vprime,tprime,qprime,pprime,aprime)
  !$omp do private(iw,k) reduction(+:bufferb)
  do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)

     do k = lpw(iw), mza
        wprime(k) = wc(k,iw) - w_avg(k)
        uprime(k) = ue(k,iw) - u_avg(k)
        vprime(k) = ve(k,iw) - v_avg(k)
        tprime(k) = thil(k,iw) - t_avg(k)
        qprime(k) = rr_v(k,iw) - q_avg(k)
        pprime(k) = press(k,iw) - p_avg(k)
!       aprime(k) = addsc(1)%sclp(k,iw) - a_avg(k)

        bufferb(k,1) = bufferb(k,1) + arw0(iw) * wprime(k) * wprime(k)
        bufferb(k,2) = bufferb(k,2) + arw0(iw) * uprime(k) * uprime(k)
        bufferb(k,3) = bufferb(k,3) + arw0(iw) * vprime(k) * vprime(k)
        bufferb(k,4) = bufferb(k,4) + arw0(iw) * tprime(k) * tprime(k)
        bufferb(k,5) = bufferb(k,5) + arw0(iw) * qprime(k) * qprime(k)
        bufferb(k,6) = bufferb(k,6) + arw0(iw) * pprime(k) * pprime(k)
!       bufferb(k,7) = bufferb(k,7) + arw0(iw) * aprime(k) * aprime(k)
     enddo

     do k = lpw(iw), mza-1
        bufferb(k, 8) = bufferb(k, 8) + arw0(iw) * wprime(k) * (uprime(k+1) + uprime(k))
        bufferb(k, 9) = bufferb(k, 9) + arw0(iw) * wprime(k) * (vprime(k+1) + vprime(k))
        bufferb(k,10) = bufferb(k,10) + arw0(iw) * wprime(k) * (tprime(k+1) + tprime(k))
        bufferb(k,11) = bufferb(k,11) + arw0(iw) * wprime(k) * (qprime(k+1) + qprime(k))
        bufferb(k,12) = bufferb(k,12) + arw0(iw) * wprime(k) * (pprime(k+1) + pprime(k))
!       bufferb(k,13) = bufferb(k,13) + arw0(iw) * wprime(k) * (aprime(k+1) + aprime(k))
     enddo

  enddo
  !$omp end do
  !$omp end parallel

  do i = 1, 6
     bufferb(2:,i) = arwm_tot(2:) * bufferb(2:,i)
  enddo

  do i = 8, 12
     bufferb(2:,i) = 0.5 * arwm_tot(2:) * bufferb(2:,i)
  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then
     if (myrank == 0) then
        call MPI_Reduce(MPI_IN_PLACE, bufferb, mza*13, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
     else
        call MPI_Reduce(bufferb, bufferb, mza*13, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
     endif
  endif
#endif

  if (myrank == 0) then
  do k = mza, 2, -1
  write(*,'(I4, 2F9.4, F10.6, F9.3, F9.1, F9.6,F8.4, 20g12.4)')  k, u_avg(k), v_avg(k), w_avg(k), &
       t_avg(k), p_avg(k), q_avg(k), ww_avg(k), uu_avg(k), tw_avg(k), pp_avg(k)
  enddo
  write(*,*)
  endif

end subroutine les_diag
