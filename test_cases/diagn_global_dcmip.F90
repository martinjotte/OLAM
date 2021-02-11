subroutine diagn_global_dcmip()

! This subroutine computes error norms for some NCAR DCMIP test cases
     
use mem_basic
use mem_ijtabs
use misc_coms
use consts_coms
use mem_grid
use oplot_coms, only: op
use mem_addsc,  only: addsc
use oname_coms, only: nl
use mem_para,   only: mgroupsize, myrank

#ifdef OLAM_MPI
use mpi
#endif

implicit none

integer, save :: ncall = 0

integer :: j,iw,kk,k
real :: raxis,u,v,w

real, save, dimension(36000) :: ge1,ge2,ge3,ge4,ge5,ge6,vctr18, &
                                q1l1,q1l2,q1li, &
                                q2l1,q2l2,q2li, &
                                q3l1,q3l2,q3li, &
                                q4l1,q4l2,q4li, &
                                q5l1,q5l2,q5li
character(len=2) :: title

real, save :: aspect = .7
real, save :: scalelab = .012
real, save :: timebeg,timeend,timedif,timeinc
real :: exner, temp
real(r8) :: enk, ent
real(r8) :: enk_sum, ent_sum

real(r8), save :: enk_sum_init, ent_sum_init, tmass_sum_init
real(r8), save :: q1mass_sum_init, q2mass_sum_init
real(r8), save :: q3mass_sum_init, q4mass_sum_init, q5mass_sum_init

real(r8) :: tmass,q1,q2,q3,q4,q5,tmass_sum

real(r8), save, allocatable :: q1_tr(:,:),q2_tr(:,:),q3_tr(:,:),q4_tr(:,:),q5_tr(:,:)

real(r8), save :: sum_absq1tr, sum_q1tr2, absq1tr_max
real(r8), save :: sum_absq2tr, sum_q2tr2, absq2tr_max
real(r8), save :: sum_absq3tr, sum_q3tr2, absq3tr_max
real(r8), save :: sum_absq4tr, sum_q4tr2, absq4tr_max
real(r8), save :: sum_absq5tr, sum_q5tr2, absq5tr_max

real(r8) :: q1mass_sum, q1mas_sum_init, sum_absq1dif, sum_q1dif2, absq1dif_max
real(r8) :: q2mass_sum, q2mas_sum_init, sum_absq2dif, sum_q2dif2, absq2dif_max
real(r8) :: q3mass_sum, q3mas_sum_init, sum_absq3dif, sum_q3dif2, absq3dif_max
real(r8) :: q4mass_sum, q4mas_sum_init, sum_absq4dif, sum_q4dif2, absq4dif_max
real(r8) :: q5mass_sum, q5mas_sum_init, sum_absq5dif, sum_q5dif2, absq5dif_max

real(r8) :: real_mixing, overshooting, range_pres, area_sum
real(r8) :: q1min, q1max, q2min, q2max, aa, bb, eps, c, root
real(r8) :: corr1, corr2, corr3, cf

integer :: nsends, ierror, k4900
real(r8), allocatable :: sendbuf(:), recvbuf(:,:)

k4900 = minloc(abs(zt(:)-4900.),1)

ncall = ncall + 1

! On the first call to this subroutine, compute the plotting time increment

if (ncall == 1) then

   allocate(q1_tr(mza,mwa),q2_tr(mza,mwa),q3_tr(mza,mwa),q4_tr(mza,mwa),q5_tr(mza,mwa))
   
   timebeg = time_istp8   / 86400.
   timeend = timmax8 / 86400.
   timedif = timeend - timebeg

   if (timedif < .03) then
      timeinc = .001
   elseif (timedif < .06) then
      timeinc = .002
   elseif (timedif < .1) then
      timeinc = .004

   elseif (timedif < .3) then
      timeinc = .01
   elseif (timedif < .6) then
      timeinc = .02
   elseif (timedif < 1.) then
      timeinc = .04

   elseif (timedif < 3.) then
      timeinc = .1
   elseif (timedif < 6.) then
      timeinc = .2
   elseif (timedif < 10.) then
      timeinc = .4

   elseif (timedif < 30.) then
      timeinc = 1.
   elseif (timedif < 60.) then
      timeinc = 2.
   elseif (timedif < 100.) then
      timeinc = 4.

   elseif (timedif < 300.) then
      timeinc = 10.
   elseif (timedif < 600.) then
      timeinc = 20.
   elseif (timedif < 1000.) then
      timeinc = 40.
   endif
      
endif

vctr18(ncall) = time_istp8 / 86400.

! Initialize summation and max/min quantities to zero

if (ncall == 1) then

   sum_absq1tr = 0.0_r8
   sum_absq2tr = 0.0_r8
   sum_absq3tr = 0.0_r8
   sum_absq4tr = 0.0_r8
   sum_absq5tr = 0.0_r8

   sum_q1tr2 = 0.0_r8
   sum_q2tr2 = 0.0_r8
   sum_q3tr2 = 0.0_r8
   sum_q4tr2 = 0.0_r8
   sum_q5tr2 = 0.0_r8

   absq1tr_max = 0.0_r8
   absq2tr_max = 0.0_r8
   absq3tr_max = 0.0_r8
   absq4tr_max = 0.0_r8
   absq5tr_max = 0.0_r8
endif

 tmass_sum = 0.0_r8

q1mass_sum = 0.0_r8
q2mass_sum = 0.0_r8
q3mass_sum = 0.0_r8
q4mass_sum = 0.0_r8
q5mass_sum = 0.0_r8

sum_absq1dif = 0.0_r8
sum_absq2dif = 0.0_r8
sum_absq3dif = 0.0_r8
sum_absq4dif = 0.0_r8
sum_absq5dif = 0.0_r8

sum_q1dif2 = 0.0_r8
sum_q2dif2 = 0.0_r8
sum_q3dif2 = 0.0_r8
sum_q4dif2 = 0.0_r8
sum_q5dif2 = 0.0_r8

absq1dif_max = 0.0_r8
absq2dif_max = 0.0_r8
absq3dif_max = 0.0_r8
absq4dif_max = 0.0_r8
absq5dif_max = 0.0_r8

real_mixing  = 0.0_r8
overshooting = 0.0_r8
range_pres   = 0.0_r8
area_sum     = 0.0_r8

q1min=0.0_r8
q1max=1.0_r8
q2min= 0.9_r8 - 0.8_r8*q1min**2
q2max= 0.9_r8 - 0.8_r8*q1max**2

aa = (q2max-q2min)/(q1max-q1min)
bb = q2min - q1min*aa
eps = 1.e-7

! Horizontal loop over all primary IW points

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
!---------------------------------------------------------------------

! Vertical loop over all active T levels

   do k = lpw(iw),mza-1

! Sums and max values for q1

      q1 = 1.0_r8
      q2 = 1.0_r8
      q3 = 1.0_r8
      q4 = 1.0_r8
      q5 = 1.0_r8

      if (naddsc >= 1) q1 = addsc(1)%sclp(k,iw)
      if (naddsc >= 2) q2 = addsc(2)%sclp(k,iw)
      if (naddsc >= 3) q3 = addsc(3)%sclp(k,iw)
      if (naddsc >= 4) q4 = addsc(4)%sclp(k,iw)
      if (naddsc >= 5) q5 = addsc(5)%sclp(k,iw)

! On the first call to this subroutine, save initial values

      if (ncall == 1) then
         q1_tr(k,iw) = q1
         q2_tr(k,iw) = q2
         q3_tr(k,iw) = q3
         q4_tr(k,iw) = q4
         q5_tr(k,iw) = q5

         sum_absq1tr = sum_absq1tr + abs(q1_tr(k,iw)) * rho(k,iw) * volt(k,iw)
         sum_absq2tr = sum_absq2tr + abs(q2_tr(k,iw)) * rho(k,iw) * volt(k,iw)
         sum_absq3tr = sum_absq3tr + abs(q3_tr(k,iw)) * rho(k,iw) * volt(k,iw)
         sum_absq4tr = sum_absq4tr + abs(q4_tr(k,iw)) * rho(k,iw) * volt(k,iw)
         sum_absq5tr = sum_absq5tr + abs(q5_tr(k,iw)) * rho(k,iw) * volt(k,iw)

         sum_q1tr2 = sum_q1tr2 + q1_tr(k,iw)**2 * rho(k,iw) * volt(k,iw)
         sum_q2tr2 = sum_q2tr2 + q2_tr(k,iw)**2 * rho(k,iw) * volt(k,iw)
         sum_q3tr2 = sum_q3tr2 + q3_tr(k,iw)**2 * rho(k,iw) * volt(k,iw)
         sum_q4tr2 = sum_q4tr2 + q4_tr(k,iw)**2 * rho(k,iw) * volt(k,iw)
         sum_q5tr2 = sum_q5tr2 + q5_tr(k,iw)**2 * rho(k,iw) * volt(k,iw)

         absq1tr_max = max(absq1tr_max, abs(q1_tr(k,iw)))
         absq2tr_max = max(absq2tr_max, abs(q2_tr(k,iw)))
         absq3tr_max = max(absq3tr_max, abs(q3_tr(k,iw)))
         absq4tr_max = max(absq4tr_max, abs(q4_tr(k,iw)))
         absq5tr_max = max(absq5tr_max, abs(q5_tr(k,iw)))

      endif

       tmass_sum =  tmass_sum +      rho(k,iw) * volt(k,iw)

      q1mass_sum = q1mass_sum + q1 * rho(k,iw) * volt(k,iw)
      q2mass_sum = q2mass_sum + q2 * rho(k,iw) * volt(k,iw)
      q3mass_sum = q3mass_sum + q3 * rho(k,iw) * volt(k,iw)
      q4mass_sum = q4mass_sum + q4 * rho(k,iw) * volt(k,iw)
      q5mass_sum = q5mass_sum + q5 * rho(k,iw) * volt(k,iw)

      sum_absq1dif = sum_absq1dif + abs(q1 - q1_tr(k,iw)) * rho(k,iw) * volt(k,iw)
      sum_absq2dif = sum_absq2dif + abs(q2 - q2_tr(k,iw)) * rho(k,iw) * volt(k,iw)
      sum_absq3dif = sum_absq3dif + abs(q3 - q3_tr(k,iw)) * rho(k,iw) * volt(k,iw)
      sum_absq4dif = sum_absq4dif + abs(q4 - q4_tr(k,iw)) * rho(k,iw) * volt(k,iw)
      sum_absq5dif = sum_absq5dif + abs(q5 - q5_tr(k,iw)) * rho(k,iw) * volt(k,iw)

      sum_q1dif2 = sum_q1dif2 + (q1 - q1_tr(k,iw))**2 * rho(k,iw) * volt(k,iw)
      sum_q2dif2 = sum_q2dif2 + (q2 - q2_tr(k,iw))**2 * rho(k,iw) * volt(k,iw)
      sum_q3dif2 = sum_q3dif2 + (q3 - q3_tr(k,iw))**2 * rho(k,iw) * volt(k,iw)
      sum_q4dif2 = sum_q4dif2 + (q4 - q4_tr(k,iw))**2 * rho(k,iw) * volt(k,iw)
      sum_q5dif2 = sum_q5dif2 + (q5 - q5_tr(k,iw))**2 * rho(k,iw) * volt(k,iw)

      absq1dif_max = max(absq1dif_max, abs(q1 - q1_tr(k,iw)))
      absq2dif_max = max(absq2dif_max, abs(q2 - q2_tr(k,iw)))
      absq3dif_max = max(absq3dif_max, abs(q3 - q3_tr(k,iw)))
      absq4dif_max = max(absq4dif_max, abs(q4 - q4_tr(k,iw)))
      absq5dif_max = max(absq5dif_max, abs(q5 - q5_tr(k,iw)))

   enddo

   if (nl%test_case == 11) then

      do k = k4900-2, k4900+2

         q1 = 1.0_r8
         q2 = 1.0_r8
         if (naddsc >= 1) q1 = addsc(1)%sclp(k,iw)
         if (naddsc >= 2) q2 = addsc(2)%sclp(k,iw)

         c = (432.0 * q1 + 6.0 * sqrt(750.0 * (2.0*q2 - 1.0)**3 + 5184.0 * q1**2))**(1./3.) / 12.0

         root = c+((5.0/24.0)-(5.0/12.0)*q2)/c
         root = min( max(root, 0.0), 1.0)
         
         corr1 = 0.9 - 0.8*q1**2
         corr2 = aa * q1 + bb
         cf    = 0.9 - 0.8*root**2

        corr3 = sqrt( (root-q1)*(root-q1)/(0.9**2) + (cf-q2)*(cf-q2)/(0.8**2) )
!         corr3 = sqrt( (q1_tr(k,iw)-q1)**2/(0.9**2) + (q2_tr(k,iw)-q2)**2/(0.8**2) )

         ! Add to the correct mixing region

!         write(*,*) q1, root, q2, cf

         if (q2 <= corr1+eps .and. q2 >= corr2-eps) then
            real_mixing = real_mixing + corr3 * arw0(iw)
         elseif (q1 <= q1max+eps .and. q1 >= q1min-eps .and. q2 <= q2min+eps .and. q2 >= q2max-eps) then
            range_pres = range_pres + corr3 * arw0(iw)
         else
            overshooting = overshooting + corr3 * arw0(iw)
         endif
         
         area_sum = area_sum + arw0(iw)
      
      enddo
   endif
   
enddo

#ifdef OLAM_MPI

if (iparallel == 1) then

   if (ncall == 1) then
      nsends = 40
   else
      nsends = 25
   endif

   allocate( sendbuf( nsends) )

   if (myrank == 0) then
      allocate( recvbuf( nsends, mgroupsize) )
   else
      allocate( recvbuf( 1, 1) )
   endif

   sendbuf( 1) =  tmass_sum

   sendbuf( 2) = q1mass_sum
   sendbuf( 3) = q2mass_sum
   sendbuf( 4) = q3mass_sum
   sendbuf( 5) = q4mass_sum
   sendbuf( 6) = q5mass_sum

   sendbuf( 7) = sum_absq1dif
   sendbuf( 8) = sum_absq2dif
   sendbuf( 9) = sum_absq3dif
   sendbuf(10) = sum_absq4dif
   sendbuf(11) = sum_absq5dif

   sendbuf(12) = sum_q1dif2
   sendbuf(13) = sum_q2dif2
   sendbuf(14) = sum_q3dif2
   sendbuf(15) = sum_q4dif2
   sendbuf(16) = sum_q5dif2

   sendbuf(17) = absq1dif_max
   sendbuf(18) = absq2dif_max
   sendbuf(19) = absq3dif_max
   sendbuf(20) = absq4dif_max
   sendbuf(21) = absq5dif_max

   sendbuf(22) = real_mixing
   sendbuf(23) = range_pres
   sendbuf(24) = overshooting
   sendbuf(25) = area_sum

   if (ncall == 1) then

      sendbuf(26) = sum_absq1tr
      sendbuf(27) = sum_absq2tr
      sendbuf(28) = sum_absq3tr
      sendbuf(29) = sum_absq4tr
      sendbuf(30) = sum_absq5tr

      sendbuf(31) = sum_q1tr2
      sendbuf(32) = sum_q2tr2
      sendbuf(33) = sum_q3tr2
      sendbuf(34) = sum_q4tr2
      sendbuf(35) = sum_q5tr2

      sendbuf(36) = absq1tr_max
      sendbuf(37) = absq2tr_max
      sendbuf(38) = absq3tr_max
      sendbuf(39) = absq4tr_max
      sendbuf(40) = absq5tr_max

   endif

   call mpi_gather( sendbuf, nsends, mpi_real8, recvbuf, nsends, mpi_real8, 0, MPI_COMM_WORLD, ierror)

   if (myrank == 0) then

      tmass_sum  = sum( recvbuf(1,1:mgroupsize) )

      q1mass_sum = sum( recvbuf(2,1:mgroupsize) )
      q2mass_sum = sum( recvbuf(3,1:mgroupsize) )
      q3mass_sum = sum( recvbuf(4,1:mgroupsize) )
      q4mass_sum = sum( recvbuf(5,1:mgroupsize) )
      q5mass_sum = sum( recvbuf(6,1:mgroupsize) )

      sum_absq1dif = sum( recvbuf( 7,1:mgroupsize) )
      sum_absq2dif = sum( recvbuf( 8,1:mgroupsize) )
      sum_absq3dif = sum( recvbuf( 9,1:mgroupsize) )
      sum_absq4dif = sum( recvbuf(10,1:mgroupsize) )
      sum_absq5dif = sum( recvbuf(11,1:mgroupsize) )

      sum_q1dif2 = sum( recvbuf(12,1:mgroupsize) )
      sum_q2dif2 = sum( recvbuf(13,1:mgroupsize) )
      sum_q3dif2 = sum( recvbuf(14,1:mgroupsize) )
      sum_q4dif2 = sum( recvbuf(15,1:mgroupsize) )
      sum_q5dif2 = sum( recvbuf(16,1:mgroupsize) )

      absq1dif_max = maxval( recvbuf(17,1:mgroupsize) )
      absq2dif_max = maxval( recvbuf(18,1:mgroupsize) )
      absq3dif_max = maxval( recvbuf(19,1:mgroupsize) )
      absq4dif_max = maxval( recvbuf(20,1:mgroupsize) )
      absq5dif_max = maxval( recvbuf(21,1:mgroupsize) )

      real_mixing  = sum( recvbuf(22,1:mgroupsize) )
      range_pres   = sum( recvbuf(23,1:mgroupsize) )
      overshooting = sum( recvbuf(24,1:mgroupsize) )
      area_sum     = sum( recvbuf(25,1:mgroupsize) )

      if (ncall == 1) then

         sum_absq1tr = sum( recvbuf(26,1:mgroupsize) )
         sum_absq2tr = sum( recvbuf(27,1:mgroupsize) )
         sum_absq3tr = sum( recvbuf(28,1:mgroupsize) )
         sum_absq4tr = sum( recvbuf(29,1:mgroupsize) )
         sum_absq5tr = sum( recvbuf(30,1:mgroupsize) )

         sum_q1tr2 = sum( recvbuf(31,1:mgroupsize) )
         sum_q2tr2 = sum( recvbuf(32,1:mgroupsize) )
         sum_q3tr2 = sum( recvbuf(33,1:mgroupsize) )
         sum_q4tr2 = sum( recvbuf(34,1:mgroupsize) )
         sum_q5tr2 = sum( recvbuf(35,1:mgroupsize) )

         absq1tr_max = maxval( recvbuf(36,1:mgroupsize) )
         absq2tr_max = maxval( recvbuf(37,1:mgroupsize) )
         absq3tr_max = maxval( recvbuf(38,1:mgroupsize) )
         absq4tr_max = maxval( recvbuf(39,1:mgroupsize) )
         absq5tr_max = maxval( recvbuf(40,1:mgroupsize) )

      endif  ! ncall == 1

      deallocate( recvbuf )
         
   endif  ! myrank == 0

   deallocate( sendbuf )
   
endif  ! iparallel == 1

#endif

! Only need to compute/plot on one node

if (myrank /= 0) return

real_mixing = real_mixing / area_sum
range_pres  = range_pres / area_sum
overshooting= overshooting / area_sum

write(*,*) " Mixing Diagnostics "
write(*,*) " Real Mixing        ", real_mixing
write(*,*) " Range Pres Unmixing", range_pres
write(*,*) " Overshooting       ", overshooting
write(*,*)

! first time in this routine, save initial values

if (ncall == 1) then
    tmass_sum_init =  tmass_sum

   q1mass_sum_init = q1mass_sum
   q2mass_sum_init = q2mass_sum
   q3mass_sum_init = q3mass_sum
   q4mass_sum_init = q4mass_sum
   q5mass_sum_init = q5mass_sum
endif

! Time series of normalized global integrals

ge1(ncall) = ( tmass_sum -  tmass_sum_init) /  tmass_sum_init

ge2(ncall) = (q1mass_sum - q1mass_sum_init) / q1mass_sum_init
ge3(ncall) = (q2mass_sum - q2mass_sum_init) / q2mass_sum_init
ge4(ncall) = (q3mass_sum - q3mass_sum_init) / q3mass_sum_init
ge5(ncall) = (q4mass_sum - q4mass_sum_init) / q4mass_sum_init
ge6(ncall) = (q5mass_sum - q5mass_sum_init) / q5mass_sum_init

! First 3 normalized global errors

q1l1(ncall) = sum_absq1dif / sum_absq1tr
q2l1(ncall) = sum_absq2dif / sum_absq2tr
q3l1(ncall) = sum_absq3dif / sum_absq3tr
q4l1(ncall) = sum_absq4dif / sum_absq4tr
q5l1(ncall) = sum_absq5dif / sum_absq5tr

q1l2(ncall) = sqrt(sum_q1dif2) / sqrt(sum_q1tr2)
q2l2(ncall) = sqrt(sum_q2dif2) / sqrt(sum_q2tr2)
q3l2(ncall) = sqrt(sum_q3dif2) / sqrt(sum_q3tr2)
q4l2(ncall) = sqrt(sum_q4dif2) / sqrt(sum_q4tr2)
q5l2(ncall) = sqrt(sum_q5dif2) / sqrt(sum_q5tr2)

q1li(ncall) = absq1dif_max / absq1tr_max
q2li(ncall) = absq2dif_max / absq2tr_max
q3li(ncall) = absq3dif_max / absq3tr_max
q4li(ncall) = absq4dif_max / absq4tr_max
q5li(ncall) = absq5dif_max / absq5tr_max

!write(6,110) &
!sum_absq1dif,     sum_absq1tr, &
!sqrt(sum_q1dif2), sqrt(sum_q1tr2), &
!absq1dif_max,     absq1tr_max

!110 format(20e10.2)

!write(6,111) &
! ge1(ncall), ge2(ncall), ge3(ncall), ge4(ncall), ge5(ncall), &
!q1l1(ncall),q2l1(ncall),q3l1(ncall),q4l1(ncall), &
!q1l2(ncall),q2l2(ncall),q3l2(ncall),q4l2(ncall), &
!q1li(ncall),q2li(ncall),q3li(ncall),q4li(ncall)

!111 format(20e10.2)

if (time_istp8 + 0.5 * dtlong > timmax8) then

! Reopen the current graphics output workstation if it is closed

   call o_reopnwk()

! Plot time series

!-------------------------------------------------------------------
   call plotback()
   call oplot_xy2('0','N',aspect,scalelab,10,0           &
                 ,ncall,  vctr18,ge1                     &
                 ,'time(days)','total mass deviation'    &
                 ,timebeg,timeend,timeinc,5  ,-1.e-12,1.e-12,0.1e-12,10  )
   call o_frame()
!-------------------------------------------------------------------

! ADDSC_1

   if (nl%naddsc >= 1) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab,10,0        &
                    ,ncall,  vctr18,ge2                  &
                    ,'time(days)','q1 mass deviation'    &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_2

   if (nl%naddsc >= 2) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab,10,0        &
                    ,ncall,  vctr18,ge3                  &
                    ,'time(days)','q2 mass deviation'    &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_3

   if (nl%naddsc >= 3) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab,10,0        &
                    ,ncall,  vctr18,ge4                  &
                    ,'time(days)','q3 mass deviation'    &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_4

   if (nl%naddsc >= 4) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab,10,0        &
                    ,ncall,  vctr18,ge5                  &
                    ,'time(days)','q4 mass deviation'    &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_5

   if (nl%naddsc >= 5) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab,10,0        &
                    ,ncall,  vctr18,ge6                  &
                    ,'time(days)','q5 mass deviation'    &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! Tracer mass norms

! ADDSC_1

   if (nl%naddsc >= 1) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0    &
                    ,ncall,  vctr18,q1l1                  &
                       ,'time(days)','q1 mass l1 norm'    &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q1l2               &
                    ,'time(days)','q1 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0    &
                    ,ncall,  vctr18,q1li                  &
                    ,'time(days)','q1 mass l_inf norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1     )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_2

   if (nl%naddsc >= 2) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q2l1               &
                    ,'time(days)','q2 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q2l2               &
                    ,'time(days)','q2 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0   &
                    ,ncall,  vctr18,q2li                 &
                    ,'time(days)','q2 mass linf norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1    )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_3

   if (nl%naddsc >= 3) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q3l1               &
                    ,'time(days)','q3 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q3l2               &
                    ,'time(days)','q3 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0   &
                    ,ncall,  vctr18,q3li                 &
                    ,'time(days)','q3 mass linf norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1    )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_4

   if (nl%naddsc >= 4) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q4l1               &
                    ,'time(days)','q4 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q4l2               &
                    ,'time(days)','q4 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0   &
                    ,ncall,  vctr18,q4li                 &
                    ,'time(days)','q4 mass linf norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1    )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_5

   if (nl%naddsc >= 5) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q5l1               &
                    ,'time(days)','q5 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0 &
                    ,ncall,  vctr18,q5l2               &
                    ,'time(days)','q5 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab,10,0   &
                    ,ncall,  vctr18,q5li                 &
                    ,'time(days)','q5 mass linf norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1    )
      call o_frame()
!-------------------------------------------------------------------
   endif

! Close the current workstation if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

   call o_clswk()

endif

return
end subroutine diagn_global_dcmip
