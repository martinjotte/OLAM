subroutine diagn_global()

use mem_basic
use mem_ijtabs
use misc_coms
use consts_coms
use mem_grid
use oplot_coms,  only: op
use mem_addsc,   only: addsc
use oname_coms, only: nl

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
real(kind=8) :: enk, ent
real(kind=8) :: enk_sum, ent_sum

real(kind=8), save :: enk_sum_init, ent_sum_init, tmass_sum_init
real(kind=8), save :: q1mass_sum_init, q2mass_sum_init
real(kind=8), save :: q3mass_sum_init, q4mass_sum_init, q5mass_sum_init

real(kind=8) :: tmass,q1,q2,q3,q4,q5,tmass_sum

real(kind=8), save, allocatable :: q1_tr(:,:),q2_tr(:,:),q3_tr(:,:),q4_tr(:,:),q5_tr(:,:)

real(8), save :: sum_absq1tr, sum_q1tr2, absq1tr_max
real(8), save :: sum_absq2tr, sum_q2tr2, absq2tr_max
real(8), save :: sum_absq3tr, sum_q3tr2, absq3tr_max
real(8), save :: sum_absq4tr, sum_q4tr2, absq4tr_max
real(8), save :: sum_absq5tr, sum_q5tr2, absq5tr_max

real(8) :: q1mass_sum, q1mas_sum_init, sum_absq1dif, sum_q1dif2, absq1dif_max
real(8) :: q2mass_sum, q2mas_sum_init, sum_absq2dif, sum_q2dif2, absq2dif_max
real(8) :: q3mass_sum, q3mas_sum_init, sum_absq3dif, sum_q3dif2, absq3dif_max
real(8) :: q4mass_sum, q4mas_sum_init, sum_absq4dif, sum_q4dif2, absq4dif_max
real(8) :: q5mass_sum, q5mas_sum_init, sum_absq5dif, sum_q5dif2, absq5dif_max

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

   sum_absq1tr = 0.
   sum_absq2tr = 0.
   sum_absq3tr = 0.
   sum_absq4tr = 0.
   sum_absq5tr = 0.

   sum_q1tr2 = 0.
   sum_q2tr2 = 0.
   sum_q3tr2 = 0.
   sum_q4tr2 = 0.
   sum_q5tr2 = 0.

   absq1tr_max = 0.
   absq2tr_max = 0.
   absq3tr_max = 0.
   absq4tr_max = 0.
   absq5tr_max = 0.
endif

 tmass_sum = 0.

q1mass_sum = 0.
q2mass_sum = 0.
q3mass_sum = 0.
q4mass_sum = 0.
q5mass_sum = 0.

sum_absq1dif = 0.
sum_absq2dif = 0.
sum_absq3dif = 0.
sum_absq4dif = 0.
sum_absq5dif = 0.

sum_q1dif2 = 0.
sum_q2dif2 = 0.
sum_q3dif2 = 0.
sum_q4dif2 = 0.
sum_q5dif2 = 0.

absq1dif_max = 0.
absq2dif_max = 0.
absq3dif_max = 0.
absq4dif_max = 0.
absq5dif_max = 0.

! Horizontal loop over all IW points

!----------------------------------------------------------------------
do iw = 2,mwa
!---------------------------------------------------------------------

! Vertical loop over all active T levels

   do k = lpw(iw),mza-1

! Sums and max values for q1

      q1 = 1.0d0
      q2 = 1.0d0
      q3 = 1.0d0
      q4 = 1.0d0
      q5 = 1.0d0

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
   
enddo

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
   call oplot_xy2('0','N',aspect,scalelab        &
                 ,ncall,  vctr18,ge1                 &
                 ,'time(days)','total mass deviation' &
                 ,timebeg,timeend,timeinc,5  ,-1.e-12,1.e-12,0.1e-12,10  )
   call o_frame()
!-------------------------------------------------------------------

! ADDSC_1

   if (nl%naddsc >= 1) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab        &
                    ,ncall,  vctr18,ge2                 &
                    ,'time(days)','q1 mass deviation' &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_2

   if (nl%naddsc >= 2) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab        &
                    ,ncall,  vctr18,ge3                 &
                    ,'time(days)','q2 mass deviation' &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_3

   if (nl%naddsc >= 3) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab        &
                    ,ncall,  vctr18,ge4                 &
                    ,'time(days)','q3 mass deviation' &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_4

   if (nl%naddsc >= 4) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab        &
                    ,ncall,  vctr18,ge5                 &
                    ,'time(days)','q4 mass deviation' &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_5

   if (nl%naddsc >= 5) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2('0','N',aspect,scalelab        &
                    ,ncall,  vctr18,ge6                 &
                    ,'time(days)','q5 mass deviation' &
                    ,timebeg,timeend,timeinc,5  ,-1.e-7,1.e-7,0.1e-7,10  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! Tracer mass norms

! ADDSC_1

   if (nl%naddsc >= 1) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q1l1               &
                       ,'time(days)','q1 mass l1 norm'    &
                 ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q1l2               &
                    ,'time(days)','q1 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q1li               &
                    ,'time(days)','q1 mass l_inf norm' &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_2

   if (nl%naddsc >= 2) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q2l1               &
                    ,'time(days)','q2 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q2l2               &
                    ,'time(days)','q2 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q2li               &
                    ,'time(days)','q2 mass linf norm'  &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_3

   if (nl%naddsc >= 3) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q3l1               &
                    ,'time(days)','q3 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q3l2               &
                    ,'time(days)','q3 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q3li               &
                    ,'time(days)','q3 mass linf norm'  &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_4

   if (nl%naddsc >= 4) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q4l1               &
                    ,'time(days)','q4 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q4l2               &
                    ,'time(days)','q4 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q4li               &
                    ,'time(days)','q4 mass linf norm'  &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
   endif

! ADDSC_5

   if (nl%naddsc >= 5) then
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q5l1               &
                    ,'time(days)','q5 mass l1 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q5l2               &
                    ,'time(days)','q5 mass l2 norm'    &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
      call o_frame()
!-------------------------------------------------------------------
      call plotback()
      call oplot_xy2log10('0','N',aspect,scalelab      &
                    ,ncall,  vctr18,q5li               &
                    ,'time(days)','q5 mass linf norm'  &
                    ,timebeg,timeend,timeinc,5  ,-4,1  )
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
end subroutine diagn_global

