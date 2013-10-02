!===============================================================================
! OLAM version 4.0

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine olam_run(name_name)

#ifdef IEEE_ARITHMETIC
use, intrinsic :: ieee_arithmetic
#endif

use misc_coms,   only: io6, time8, iflag, runtype, hfilin, time_istp8, nzp,    &
                       expnme, mdomain, ngrids, initial, iswrtyp, ilwrtyp,     &
                       meshtype, timmax8, alloc_misc, iparallel, ipar_out,     &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim, &
                       time_prevhist, rinit, rinit8, debug_fp, init_nans

use leaf_coms,   only: nzg, nzs, isfcl, nwl, mwl
use mem_leaf,    only: fill_jland
use sea_coms,    only: nws, mws
use mem_sea,     only: fill_jsea

use mem_ijtabs,  only: istp, mrls, fill_jtabs, itab_u, itab_v, itab_w
use oplot_coms,  only: op
use mem_grid,    only: nma, nua, nva, nwa, mma, mua, mva, mwa, mza, zm, zt
use mem_basic,   only: alloc_basic, wc, vxe, vye, vze
use micro_coms,  only: gnu
use mem_nudge,   only: nudflag, nudnxp, fill_jnudge
use mem_rayf,    only: rayf_init
use mem_sflux,   only: init_fluxcells, fill_jflux, mseaflux, mlandflux
use mem_para,    only: myrank
use consts_coms, only: r8
use oname_coms,  only: nl
use olam_mpi_atm,only: olam_alloc_mpi, mpi_send_w, mpi_recv_w, alloc_mpi_sndrcv_bufs
use ed_misc_coms,only: ed2_active, ed2_namelist
use hcane_rz,    only: init_hurr_step, hurricane_init
use obnd,        only: trsets, lbcopy_w
use var_tables,  only: nvar_par, vtab_r, nptonv
use mem_swtc5_refsoln_cubic

use mem_average_vars, only: reset_mavg_vars, reset_davg_vars

implicit none

character(len=*), intent(in) :: name_name

integer :: i,ifm,nndtflg,ifileok,ierr,iplt_file,mrl, imavg_file, idavg_file
integer :: mwa_prog, mua_prog, mva_prog
real :: w1,w2,t1,t2,wtime_start
real, external :: walltime
character(len=128) :: mavgfile, davgfile

wtime_start = walltime(0.)
w1 = walltime(wtime_start)
call cpu_time(t1)

istp  = 1
time8 = 0.0_r8
iflag = 0

! Read, check, and copy namelist variables

write(io6,'(/,a)') 'olam_run calling namelist read'
call read_nl(name_name)

write(io6,'(/,a)') 'olam_run checking namelist values'
call oname_check()

write(io6,'(/,a)') 'olam_run calling namelist copy'
call copy_nl('ALL_CASES')

if (runtype == 'HISTORY') then
   write(io6,'(/,a)') 'olam_run reading common values from history file'
   call history_start('COMMIO')
elseif ((runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE')) then
   write(io6,'(/,a)') 'olam_run reading common values from plot file'
   hfilin = op%plt_files(1)
   call history_start('COMMIO')
else
   call copy_nl('NOT_HISTORY')
endif

! Set constants for initializing model arrays to zeros or, if supported, NANs
! depending on namelist variable init_nans

rinit  = 0.0
rinit8 = 0.0_r8

#ifdef IEEE_ARITHMETIC
if (init_nans) then
   if (ieee_support_nan(1.0)) then
      rinit = ieee_value(1.0, ieee_signaling_nan)
   endif
endif
#endif

#ifdef IEEE_ARITHMETIC
if (init_nans) then
   if (ieee_support_nan(1.0_r8)) then
      rinit8 = ieee_value(1.0_r8, ieee_signaling_nan)
   endif
endif
#endif

! Should we enable parallel output?

ipar_out = 0

#if defined(OLAM_MPI) && defined(OLAM_PARALLEL_HDF5)
if (iparallel == 1 .and. nl%ipar_out == 1) then
   ipar_out = 1
endif
#endif

! If debugging, halt on illegal floating operations if supported

#ifdef IEEE_ARITHMETIC
if (ieee_support_halting(ieee_invalid) .and. debug_fp) then
   call ieee_set_halting_mode(ieee_usual, .true.)
endif
#endif

! Get abs seconds of simulation start and current simulation time
! Note that itime1 has format of hhmm, and date_abs_secs2 needs hhmmss

time_istp8 = time8

call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,s1900_init)
s1900_sim = s1900_init + time8

! Print initial banner

write(io6,'(a1)')         ' '
write(io6,'(a1,78a1)')    ' ',('*',i=1,78)
write(io6,'(2a1,a42)')    ' ','*','    OLAM version 4.0'
write(io6,'(2a1)')        ' ','*'
write(io6,'(2a1,a3,a64)') ' ','*','   ',EXPNME
write(io6,'(a1,78a1)')    ' ',('*',i=1,78)

! MAKESFC/MAKEGRID/PLOTONLY/PARCOMBINE runs should be single-processor

if ((runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE') .or.  &
    (runtype == 'MAKESFC' ) .or. (runtype == 'MAKEGRID') ) then
   if (iparallel == 1) then
      write(io6,*) trim(runtype)//' will only be done on a single process.'
      iparallel = 0
      ipar_out  = 0
      if (myrank > 0) go to 1000
   endif
endif

! Initialize various oplot parameters

call oplot_init()

! Generate or read from files the full-domain ATMOS, LAND, SEA, and FLUX grids

call gridinit()

! If RUNTYPE = 'MAKESFC' or 'MAKEGRID', run is finished; EXIT

if (runtype == 'MAKESFC' .or. runtype == 'MAKEGRID') then
   call gridset_print()
   write(io6,*)
   write(io6,*) trim(runtype) // ' run complete'
   go to 1000
endif

! Initial allocation of communication buffer arrays based on # of nodes

call alloc_mpi_sndrcv_bufs()

! reading only the needed to para_decomp from gridfile

call gridfile_read_pd()

! If run is parallel, assign each grid cell (ATM and, if ISFCL = 1, LAND and SEA)
! to one of multiple subdomains. If not, the grid remains unchanged

write(io6,'(/,a)') 'olam_run calling para_decomp'

call para_decomp()

write(io6,'(/,a,2i7)') 'olam_run after para_decomp',nwl,nws

! Set up itab data types and grid coordinate arrays for current node, and 
! reallocate memory for current node

call para_init() 

write(io6,'(/,a)') 'olam_run after para_init'

call gridset_print()

mwa_prog = 0
do i=1,mwa
   if (itab_w(i)%irank == myrank) mwa_prog = mwa_prog + 1
enddo
write(io6,*)
write(io6,'(a,i8)') ' # of prognostic W points on this node = ', mwa_prog
   
if (meshtype == 1) then
   
   mua_prog = 0
   do i=1,mua
      if (itab_u(i)%irank == myrank) mua_prog = mua_prog + 1
   enddo
   write(io6,'(a,i8)') ' # of prognostic U points on this node = ', mua_prog
   write(io6,*)

else

   mva_prog = 0
   do i=1,mva
      if (itab_v(i)%irank == myrank) mva_prog = mva_prog + 1
   enddo
   write(io6,'(a,i8)') ' # of prognostic V points on this node = ', mva_prog
      
endif

write(io6,'(/,a)' ) 'Local model indices:'
write(io6,'(a,i8)') ' mma = ',mma
write(io6,'(a,i8)') ' mua = ',mua
write(io6,'(a,i8)') ' mva = ',mva
write(io6,'(a,i8)') ' mwa = ',mwa

if (isfcl == 1) then
   write(io6,'(a,i8)')   ' mwl = ',mwl
   write(io6,'(a,i8)')   ' mws = ',mws
   write(io6,'(a,i8)')   ' mlandflux = ',mlandflux
   write(io6,'(a,i8)')   ' mseaflux  = ',mseaflux
endif

! Initialize dtlm, dtsm, ndtrat, and nacoust, 
! and compute the timestep schedule for all grid operations.

write(io6,'(/,a)') 'olam_run calling modsched'

call modsched()

write(io6,'(/,a)') 'olam_run calling fill_jtabs'

call fill_jtabs(mma,mua,mva,mwa)

if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) call fill_jnudge()

! Allocate and fill jsea, jland, and jflux data structures

if (isfcl == 1) then
   write(io6,'(/,a)') 'olam_run calling fill_jsea, fill_jland, fill_jflux'

   call fill_jsea()
   call fill_jland()
   call fill_jflux()
endif

! Allocate column initial state arrays

call alloc_misc(mza)

write(io6,'(/,a)') 'olam_run calling jnmbinit'

call jnmbinit()

!------------------------------------------------------------------
! If we got here, we are doing an actual simulation or PLOTONLY run
!------------------------------------------------------------------

! Allocate remainder of main model memory (for 'INITIAL', 'HISTORY', or
! 'PLOTONLY' run).  Allocate variable tables and fill variable tables for
! history files and parallel communication.

write(io6,'(/,a)') 'olam_run calling olam_mem_alloc'

call olam_mem_alloc()

if (iparallel == 1) then
   write(io6,'(/,a)') 'olam_run calling olam_alloc_mpi'

   call olam_alloc_mpi(mza,mrls)

   if (isfcl == 1) then
      write(io6,'(/,a)') 'olam_run calling olam_alloc_mpi_land'

      call olam_alloc_mpi_land(mrls)

      write(io6,'(/,a)') 'olam_run calling olam_alloc_mpi_sea'

      call olam_alloc_mpi_sea(mrls)

      write(io6,'(/,a)') 'olam_run after olam_alloc_mpi_sea'
   endif
endif

! Initialize primary atmospheric fields

if (initial == 1) then
   write(io6,'(/,a)') 'olam_run calling inithh'
   call inithh()            ! Horizontally-homogeneous initialization
elseif (initial == 2) then
   write(io6,'(/,a)') 'olam_run calling isan driver(0)'
   call isan_driver(0) 
! (If in future, initialization is not automatically done for history restart,
!  isan_driver(0) will still need to be called when nudging is to be done in 
!  order to set current value of IFGFILE.)

elseif (initial == 3) then
   write(io6,'(/,a)') 'olam_run calling fldslhi'
   call fldslhi()           ! Longitudinally-homogeneous initialization
endif

!------------------------------------------------------------
! Call fill_swtc5 to read in and initialize reference
! solution for time = 0 and time = 15d.

if (nl%test_case == 2 .or. nl%test_case == 5) then
   call swtc_init()
endif
!------------------------------------------------------------

! Diagnose earth-coordinate velocities

if (runtype == 'INITIAL') then
   mrl = 1
   call diagvel_t3d(mrl)

   if (iparallel == 1) then
      call mpi_send_w('V', vxe=vxe, vye=vye, vze=vze)
      call mpi_recv_w('V', vxe=vxe, vye=vye, vze=vze)
   endif

   call lbcopy_w(mrl, a1=vxe, a2=vye, a3=vze)
endif

! A good place to initialize added scalars

! Initialize 3d microphysics fields (if level = 3) and other microphysics
! quantities

write(io6,'(/,a)') 'olam_run calling micinit'

call micinit()

!-------------------------------------------------------------------------------
if (runtype == "INITIAL") then
   if (init_hurr_step > 0) call hurricane_init()
endif
!-------------------------------------------------------------------------------

call trsets()  

! For parallel run, send and receive initialized scalars

mrl = 1

if (iparallel == 1) then
   call mpi_send_w('S', domrl=mrl)  ! Send scalars
   call mpi_recv_w('S', domrl=mrl)  ! Recv scalars
endif

do i = 1, nvar_par
   call lbcopy_w(mrl, a1=vtab_r(nptonv(i))%rvar2_p)
enddo

! Start up radiation scheme

if (iswrtyp > 0 .or. ilwrtyp > 0) then
   write(io6,'(/,a)') 'olam_run calling radinit'
   call radinit()
endif

! Check if LEAF will be used

if (isfcl == 1) then

! Start up LEAF

   write(io6,'(/,a)') 'olam_run calling leaf4_startup'
   call leaf4_startup()

! Start up SEA

   write(io6,'(/,a)') 'olam_run calling sea_startup'
   call sea_startup()

! Start up ED2

#ifdef USE_ED2
   if (ed2_active == 1) then
      call ed_1st_master(0,1,0,0,trim(ed2_namelist))
      call ed_driver(1)
   endif
#endif

! Initialize leaf fields that depend on atmosphere

   write(io6,'(/,a)') 'olam_run calling leaf4_init_atm'
   call leaf4_init_atm()

#ifdef USE_ED2
   if (ed2_active == 1) then
      write(io6,'(/,a)') 'olam_run calling ed_driver 2'
      call ed_driver(2)
   endif
#endif

   write(io6,'(/,a)') 'olam_run calling sea_init_atm'
   call sea_init_atm()

endif

! If using variable initialization and polygon nudging, read most recent
! and next observational analyses, and fill nudging polygons for both

if (initial == 2 .and. nudflag == 1) then
   write(io6,'(/,a)') 'olam_run calling isan_driver(1)'
   call isan_driver(1)
endif

! Initialize Rayleigh friction profile

write(io6,'(/,a)') 'olam_run calling rayf_init'
call rayf_init(mua,mva,mwa,mza)

! For shallow water test case 5, read in and initialize reference
! solution for time = 0 and time = 15d.

if (nl%test_case == 5) then
   call fill_swtc5()
endif

! If this is 'PLOTONLY' run, loop through input history files, plot 
! specified fields, and exit

if ((runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE')) then

   write(io6,'(/,a)') 'beginning '//trim(runtype)//' loop'

   do iplt_file = 1,op%nplt_files

      hfilin = op%plt_files(iplt_file)
      call history_start('COMMIO')
      call history_start('HISTREAD')

      ! Earth-coordinate winds not saved in history file

      mrl = 1
      call diagvel_t3d(mrl)

      if (iparallel == 1) then
         call mpi_send_w('V', vxe=vxe, vye=vye, vze=vze)
         call mpi_recv_w('V', vxe=vxe, vye=vye, vze=vze)
      endif

      call lbcopy_w(mrl, a1=vxe, a2=vye, a3=vze)

      if (runtype == 'PLOTONLY') then

! If month-average or day-average plots are specified, do them when
! iplt_file = 1 and then exit iplt_file do loop.  IT IS ASSUMED THAT 
! MONTH-AVERAGE AND DAY-AVERAGE FILES ARE SINGLE NON-PARALLEL FILES.

         if (nl%ioutput_mavg == 1 .and. nl%nmavg_files > 0) then
            do imavg_file = 1,nl%nmavg_files
               mavgfile = nl%mavg_files(imavg_file)
               call read_mavg_vars(mavgfile)
               call plot_fields(0)
            enddo
            exit
         endif

         if (nl%ioutput_davg == 1 .and. nl%ndavg_files > 0) then
            do idavg_file = 1,nl%ndavg_files
               davgfile = nl%davg_files(idavg_file)
               call read_davg_vars(davgfile)
               call plot_fields(0)
            enddo
            exit
         endif

! If month-average and day-average plots are NOT specified, do regular plots
! from history files

         call plot_fields(0)

! For shallow water test cases, compute error norms

         if (nl%test_case == 2 .or. nl%test_case == 5) then
            call diagn_global_swtc()
         endif

      else
         call history_write('STATE')
      endif

   enddo

   write(io6,'(/,a)') trim(runtype) //' run complete'
   go to 1000
endif

! REPLACE INITIAL FIELDS WITH HISTORY READ

if (runtype == 'HISTORY') then
   write(io6,*) 'olam_run calling history_start'
   call history_start('HISTREAD')
   write(io6,*) 'olam_run finished history_start'

   ! Earth cartesian velocities not saved in history file

   mrl = 1
   call diagvel_t3d(mrl)

   if (iparallel == 1) then
      call mpi_send_w('V', vxe=vxe, vye=vye, vze=vze)
      call mpi_recv_w('V', vxe=vxe, vye=vye, vze=vze)
   endif

   call lbcopy_w(mrl, a1=vxe, a2=vye, a3=vze)
endif

! Initialize PBL quantities

call pbl_init()

write(io6,'(/,a)') 'olam_run calling plot_fields'

call plot_fields(0)

if (runtype /= 'HISTORY') then
   write(io6,'(/,a)') 'olam_run calling history_write'
   call history_write('STATE')
   write(io6,'(/,a)') 'olam_run finished history_write'
endif

time_prevhist = time8
w2 = walltime(wtime_start)
call cpu_time(t2)
write(io6,'(a,2f13.3)') '++++++++++ CPU - wall time: olam initialization: ',t2-t1,w2-w1

write(io6,'(/,a)') 'olam_run completed initialization'

! Exit if doing a zero time run
if (time8 >= timmax8) go to 1000

! Initialize field average arrays

if (nl%ioutput_mavg == 1) call reset_mavg_vars()
if (nl%ioutput_davg == 1) call reset_davg_vars()

! Call the model time integration driver

write(io6,'(/,a)') 'olam_run calling model'
call model()
write(io6,'(/,a)') 'subroutine model returned to olam_run'

! OLAM finished, clean up some last things...

1000 continue

write(io6,'(/,a)') 'olam_run calling o_clsgks'

call o_clsgks()

write(io6,'(/,a)') 'olam_run returning to olammain'

return
end subroutine olam_run

!=====================================================================

subroutine model()

use misc_coms, only: io6, time8, timmax8, dtlm, time_istp8, simtime,  &
                     current_time, s1900_init, s1900_sim
use consts_coms, only: r8
use oname_coms,  only: nl

implicit none

!   +------------------------------------------------------------------
!   ! This routine drives the entire time integration process
!   !   for a non-parallel run.
!   +------------------------------------------------------------------
!

integer :: mstp
real :: wtime_start,t1,wtime1,wtime2,t2,wtime_tot
real, external :: walltime
character(len=40) :: stepc1,stepc2,stepc3,stepc4,stepc5
type(simtime) :: begtime
real(r8) :: time8p, bias

 write(io6,*) 'starting subroutine MODEL'

wtime_start = walltime(0.)

! Start the timesteps

mstp   = 0
bias   = 1.e-7_r8 * dtlm(1) ! A number small compared to the timestep
time8p = time8 + bias       ! Slightly forward biased time

do while (time8p < timmax8)

   begtime = current_time

! CPU timing information

   call cpu_time(t1)
   wtime1 = walltime(wtime_start)

! Start the timestep schedule to loop through all grids and advance them
! in time an increment equal to dtlm(1).

   call timestep()

   mstp = mstp + 1
   time8 = time8 + dtlm(1)
   time_istp8 = time8
   time8p = time8 + bias  ! Slightly forward biased time
   s1900_sim = s1900_init + time8

   call update_model_time(current_time, dtlm(1))

   wtime2 = walltime(wtime_start)
   call cpu_time(t2)

   write (stepc1,'(i8)'   ) mstp
   write (stepc2,'(f13.1)') time8
   write (stepc3,'(f10.3)') time8/86400.
   write (stepc4,'(f9.3)' ) t2-t1
   write (stepc5,'(f9.3)' ) wtime2-wtime1

   stepc1 = ' [nstep = '//trim(adjustl(stepc1))//']'
   stepc2 = '   [simtime = '//trim(adjustl(stepc2))//' sec'
   stepc3 = ' = '//trim(adjustl(stepc3))//' days]'
   stepc4 = '   [cpu,wall(sec) = '//trim(adjustl(stepc4))
   stepc5 = ' , '//trim(adjustl(stepc5))//']'
   
   write(io6,'(a)')  &
      trim(stepc1)//trim(stepc2)//trim(stepc3)//trim(stepc4)//trim(stepc5)

! Add current contribution to time-averaged variables

   if (nl%ioutput_mavg == 1) call inc_mavg_vars()
   if (nl%ioutput_davg == 1) call inc_davg_vars()

! Check schedule for I/O operations and perform those that are due 

   call olam_output()
   
enddo

wtime_tot = walltime(wtime_start)
write(io6, '(//,a,f13.3)') ' -----Total elapsed time: ',wtime_tot

return
end subroutine model

!===========================================================================

subroutine olam_output()

use misc_coms,   only: io6, time8, dtlm, iflag, frqstate, timmax8, initial, &
                       s1900_sim, time_prevhist, iyear1, imonth1, idate1,   &
                       itime1
use leaf_coms,   only: isfcl, iupdndvi, indvifile, s1900_ndvi
use sea_coms,    only: iupdsst, iupdseaice, isstfile, iseaicefile,  &
                       s1900_sst, s1900_seaice, isstflg, iseaiceflg
use oplot_coms,  only: op
use mem_nudge,   only: nudflag
use isan_coms,   only: ifgfile, s1900_fg
use consts_coms, only: r8
use oname_coms,  only: nl
use hcane_rz,    only: init_hurr_step, hurricane_track

use mem_average_vars, only: reset_mavg_vars, reset_davg_vars

implicit none

integer  :: ierr, ifm, ifileok
real(r8) :: time8p, bias
integer  :: outyear, outmonth, outdate, outhour

bias    = 1.e-7_r8 * dtlm(1) ! A number small compared to the timestep
time8p  = time8 + bias       ! Slightly forward biased time

!-------------------- SPECIAL - HURRICANE TRACKING ------------------
if (init_hurr_step == 1 .or. init_hurr_step == 2) then
    if (time8p > timmax8) then
       call hurricane_track(3)
    elseif (mod( time8p, 3600.0_r8) < dtlm(1)) then
       call hurricane_track(2)
    else
       call hurricane_track(1)
    endif
endif
!-------------------------------------------------------------------

if (mod(time8p,op%frqplt) < dtlm(1) .or. iflag == 1) then
   call plot_fields(0)
endif

if (mod(time8p,frqstate) < dtlm(1)  .or.  &
   time8 >= timmax8 - bias .or. iflag == 1) then
   call history_write('INST')
   time_prevhist = time8
endif

! Output of time-averaged quantities

if (mod(time8p,86400.0_r8) < dtlm(1)) then
   call date_add_to8(iyear1,imonth1,idate1,itime1,time8p,'s',outyear,  &
        outmonth,outdate,outhour)

   if (nl%ioutput_davg == 1) then
      call norm_davg_vars()
      call write_davg_vars(outyear,outmonth,outdate)
      call reset_davg_vars()
   endif

   if (nl%ioutput_mavg == 1 .and. outdate == 1) then
      call norm_mavg_vars()
      call write_mavg_vars(outyear,outmonth)
      call reset_mavg_vars()
   endif

endif

! For idealized test cases, compute error norms

if (nl%test_case == 2 .or. nl%test_case == 5) then
   call diagn_global_swtc()
endif

if (isfcl == 1 .and. iupdsst == 1) then
   if (isstflg == 1) then
      if (s1900_sim >= s1900_sst(isstfile) .and. time8p < timmax8) then
         call sst_database_read(1)
      endif
   elseif (isstflg == 2) then
      if (s1900_sim >= s1900_fg(isstfile) .and. time8p < timmax8) then
         call read_sst_analysis(1)
      endif
   endif
endif

if (isfcl == 1 .and. iupdseaice == 1) then
   if (iseaiceflg == 1) then
      if (s1900_sim >= s1900_seaice(iseaicefile) .and. time8p < timmax8) then
         call seaice_database_read(1)
      endif
   elseif (iseaiceflg == 2) then
      if (s1900_sim >= s1900_fg(iseaicefile) .and. time8p < timmax8) then
         call read_seaice_analysis(1)
      endif
   endif
endif

if (isfcl == 1 .and. iupdndvi == 1) then
   if (s1900_sim >= s1900_ndvi(indvifile) .and. time8p < timmax8) then
      call ndvi_database_read(1)
   endif
endif

if (initial == 2 .and. nudflag == 1) then
   if (s1900_sim >= s1900_fg(ifgfile) .and. time8p < timmax8) then

      write(io6,*) ' '
      write(io6,*) 'olam_output ',time8,s1900_sim,s1900_fg(ifgfile),timmax8

      call isan_driver(1)

   endif
endif

if (iflag == 1) stop 'IFLAG'

return
end subroutine olam_output
