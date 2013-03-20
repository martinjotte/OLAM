subroutine diagn_global_swtc()

! This subroutine computes error norms for shallow water test cases 1, 2, and 5
     
use mem_basic
use mem_micro
use micro_coms
use mem_ijtabs
use misc_coms
use consts_coms
use mem_grid
use oplot_coms,  only: op
use oname_coms,  only: nl
use mem_para,    only: mgroupsize, myrank

use mem_swtc5_refsoln_cubic

#ifdef OLAM_MPI
use mpi
#endif

implicit none

integer, save :: ncall = 0, ncall_tot

integer :: iw

real :: height, vol

real :: sum_abshdif, sum_abshtr, sum_hdif2, sum_htr2, abshdif_max, abshtr_max

real, save, allocatable :: ge1(:), ge2(:), ge3(:), vctr18(:)
real, save, allocatable :: height_init(:)

real, save :: aspect = .7
real, save :: scalelab = .014
real, save :: timebeg,timeend,timedif,timeinc

real :: zanal0_swtc5,zanal_swtc5

integer :: nsends, ierror
real, allocatable :: sendbuf(:), recvbuf(:,:)

ncall = ncall + 1

! On the first call to this subroutine, compute the plotting time increment

if (ncall == 1) then

   ncall_tot = int(timmax8 / dtlm(1)) + 10

   allocate (ge1(ncall_tot), ge2(ncall_tot), ge3(ncall_tot), vctr18(ncall_tot))
   allocate (height_init(nwa))

   timebeg = time8   / 86400.
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

! Initialize summation and max/min quantities to zero

sum_abshdif = 0.
sum_abshtr  = 0.

sum_hdif2   = 0.
sum_htr2    = 0.

abshdif_max = 0.
abshtr_max  = 0.

! Horizontal loop over all IW points

!----------------------------------------------------------------------
do iw = 2,mwa
!---------------------------------------------------------------------

   if (nl%test_case == 5) then

! Due to problem with reference solution (or its interpolation), 
! omit high latitude points

      if (glatw(iw) > 85. .or. glatw(iw) < -85.) cycle

! Get SWTC5 reference solution (for days 0 and 15) for this point

      call npr_bicubics(zanal00_swtc5,glatw(iw),glonw(iw),zanal0_swtc5)
      call npr_bicubics(zanal15_swtc5,glatw(iw),glonw(iw),zanal_swtc5)

   else

! Case for SWTC2

      zanal0_swtc5 = 0.
      zanal_swtc5 = 0.
   
   endif

! Get height and volume for current IW point

   height = rho(2,iw)
   vol = volt(2,iw)

! first time in this routine, save initial values: uh(iw) = u; vh(iw) = v

   if (ncall == 1) then
      height_init(iw) = height
   endif

! Find sums and max values for height

   sum_abshdif = sum_abshdif + vol *  &
                 abs(height - height_init(iw) - zanal_swtc5 + zanal0_swtc5)

   sum_hdif2 = sum_hdif2 + vol *  &
               (height - height_init(iw) - zanal_swtc5 + zanal0_swtc5)**2

   abshdif_max = max(abshdif_max,  &
                 abs(height - height_init(iw) - zanal_swtc5 + zanal0_swtc5))

   if (nl%test_case == 2) then
      sum_abshtr = sum_abshtr + vol * abs(height_init(iw))
      sum_htr2   = sum_htr2   + vol * height_init(iw)**2
      abshtr_max = max(abshtr_max,    height_init(iw))
   else
      sum_abshtr = sum_abshtr + vol * abs(zanal_swtc5)
      sum_htr2   = sum_htr2   + vol * zanal_swtc5**2
      abshtr_max = max(abshtr_max,    zanal_swtc5)
   endif

enddo

#ifdef OLAM_MPI
if (iparallel == 1) then

   nsends = 6

   allocate( sendbuf( nsends) )
   allocate( recvbuf( nsends, mgroupsize) )

   sendbuf(1) = sum_abshdif
   sendbuf(2) = sum_hdif2
   sendbuf(3) = abshdif_max
   sendbuf(4) = sum_abshtr
   sendbuf(5) = sum_htr2
   sendbuf(6) = abshtr_max

   call mpi_allgather( sendbuf, nsends, mpi_real, &
                       recvbuf, nsends, mpi_real, MPI_COMM_WORLD, ierror )

   sum_abshdif = sum   ( recvbuf(1,1:mgroupsize) )
   sum_hdif2   = sum   ( recvbuf(2,1:mgroupsize) )
   abshdif_max = maxval( recvbuf(3,1:mgroupsize) )
   sum_abshtr  = sum   ( recvbuf(4,1:mgroupsize) )
   sum_htr2    = sum   ( recvbuf(5,1:mgroupsize) )
   abshtr_max  = maxval( recvbuf(6,1:mgroupsize) )

   deallocate( sendbuf )
   deallocate( recvbuf )

endif
#endif

! First 3 normalized global errors (Williamson et al 1992 Eqs. 82-84)

ge1(ncall) = sum_abshdif / sum_abshtr
ge2(ncall) = sqrt(sum_hdif2) / sqrt(sum_htr2)
ge3(ncall) = abshdif_max / abshtr_max

vctr18(ncall) = time8 / 86400.

if (nl%test_case == 2 .and. time8 + .5 * dtlong <  432000.        ) return
if (nl%test_case == 5 .and. time8 + .5 * dtlong < 1296000. - 7200.) return

! Reopen the current graphics output workstation if it is closed

call o_reopnwk()

!----------------------------------------------------------------------
! Plot 3 height norms on single frame
!----------------------------------------------------------------------

 call plotback()
 call oplot_xy2log10('0','N',aspect,scalelab &
               ,ncall,  vctr18,ge1          &
               ,'time(days)',' '             &
               ,timebeg,timeend,timeinc,5  ,-6,-1  )
 call oplot_xy2log10('0','N',aspect,scalelab &
               ,ncall,  vctr18,ge2          &
               ,' ',' '                      &
               ,timebeg,timeend,timeinc,5  ,-6,-1  )
 call oplot_xy2log10('0','N',aspect,scalelab &
               ,ncall,  vctr18,ge3          &
               ,' ',' '                      &
               ,timebeg,timeend,timeinc,5  ,-6,-1  )
 call o_frame()

! Close the current workstation if output
! is to a NCAR graphics meta file. This allows viewing the complete
! meta file (including the last frame) during a run and in case the
! simulation crashes.

if ((trim(runtype) /= 'PLOTONLY') .and. (op%plttype == 0)) then
   call o_clswk()
endif

return
end subroutine diagn_global_swtc
