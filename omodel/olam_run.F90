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

use misc_coms,   only: io6, time8, iflag, runtype, hfilin, time_istp8, &
                       expnme, mdomain, ngrids, initial, iswrtyp, ilwrtyp, &
                       meshtype, nzp, timmax8, alloc_misc, iparallel, &
                       iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim, &
                       time_prevhist, rinit, rinit8, debug_fp, init_nans

use leaf_coms,   only: nzg, nzs, isfcl, nwl, mwl
use mem_leaf,    only: fill_jland
use sea_coms,    only: nws, mws
use mem_sea,     only: fill_jsea

use mem_ijtabs,  only: istp, mrls, fill_jtabs, itab_u, itab_v, itab_w
use oplot_coms,  only: op
use mem_grid,    only: nma, nua, nva, nwa, mma, mua, mva, mwa, mza, zm, zt
use mem_basic,   only: alloc_basic, wc
use micro_coms,  only: gnu
use ed_options,  only: ied_offline
use mem_nudge,   only: nudflag
use mem_rayf,    only: rayf_init
use mem_sflux,   only: init_fluxcells, fill_jflux, mseaflux, mlandflux
use mem_para,    only: myrank
use consts_coms, only: r8
use olam_mpi_atm, only: olam_alloc_mpi, mpi_send_w, mpi_recv_w, alloc_mpi_sndrcv_bufs

implicit none

character(len=*), intent(in) :: name_name

integer :: i,ifm,nndtflg,ifileok,ierr,iplt_file
integer :: mwa_prog, mua_prog, mva_prog
real :: w1,w2,t1,t2,wtime_start
real, external :: walltime

wtime_start = walltime(0.)
w1 = walltime(wtime_start)
call cpu_time(t1)

istp = 1
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
   hfilin = trim(op%plt_files(1))
   call history_start('COMMIO')
else
   call copy_nl('NOT_HISTORY')
endif

! Set constants for initializing model arrays to zeros or, if supported, NANs
! depending on namelist variable init_nans

rinit = 0.0
rinit8 = 0.0_r8

if (init_nans) then
#ifdef IEEE_ARITHMETIC
   if (ieee_support_nan(1.0)) then
      rinit = ieee_value(1.0, ieee_signaling_nan)
   endif
#endif
endif

if (init_nans) then
#ifdef IEEE_ARITHMETIC
   if (ieee_support_nan(1.0_r8)) then
      rinit8 = ieee_value(1.0_r8, ieee_signaling_nan)
   endif
#endif
endif

! If debugging, halt on illegal floating operations if supported

#ifdef IEEE_ARITHMETIC
if (ieee_support_halting(ieee_invalid) .and. debug_fp) then
   call ieee_set_halting_mode(ieee_usual, .true.)
endif
#endif

! Get abs seconds of simulation start and current simulation time

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
      if (myrank > 0) go to 1000
   endif
endif

! Initialize various oplot parameters

call oplot_init()

! Generate or read from files the full-domain ATMOS, LAND, SEA, and FLUX grids

call gridinit()

! If RUNTYPE = 'MAKESFC' or 'MAKEGRID', run is finished; EXIT

if (trim(runtype) == 'MAKESFC' .or. trim(runtype) == 'MAKEGRID') then
   write(io6,*) trim(runtype)//' run complete'
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
write(io6,'(a,i8)')   ' mma = ',mma
write(io6,'(a,i8)')   ' mua = ',mua
write(io6,'(a,i8)')   ' mva = ',mva
write(io6,'(a,i8)')   ' mwa = ',mwa
write(io6,*)

mwa_prog = 0
do i=1,mwa
   if (itab_w(i)%irank == myrank) mwa_prog = mwa_prog + 1
enddo
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
   write(io6,*)
      
endif

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

write(io6,'(/,a)') 'olam_run calling olam_alloc_mpi'

call olam_alloc_mpi(mza,mrls)

if (isfcl == 1) then
   write(io6,'(/,a)') 'olam_run calling olam_alloc_mpi_land'

   call olam_alloc_mpi_land(mrls)

   write(io6,'(/,a)') 'olam_run calling olam_alloc_mpi_sea'

   call olam_alloc_mpi_sea(mrls)

   write(io6,'(/,a)') 'olam_run after olam_alloc_mpi_sea'

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

!-------------------- SPECIAL - HURRICANE TRACKING ------------------
! call hurricane_track()
!-------------------------------------------------------------------

! A good place to initialize added scalars

! Initialize 3d microphysics fields (if level = 3) and other microphysics
! quantities

write(io6,'(/,a)') 'olam_run calling micinit'

call micinit()

! For parallel run, send and receive initialized scalars

if (iparallel == 1) then
   call mpi_send_w('S')  ! Send scalars
   call mpi_recv_w('S')  ! Recv scalars
endif

if (iswrtyp == 3 .or. ilwrtyp == 3) then
   write(io6,'(/,a)') 'olam_run calling harr_radinit'
   call harr_radinit()
endif

! Check if LEAF3 will be used

if (isfcl == 1) then

! Start up LEAF3

   write(io6,'(/,a)') 'olam_run calling leaf3_startup'
   call leaf3_startup()

! Start up SEA

   write(io6,'(/,a)') 'olam_run calling sea_startup'
   call sea_startup()

! Initialize meteorological drivers for an offline ED run.

   if(ied_offline == 1)then
      call init_offline_met()
      call read_offline_met_init()
   endif

! Initialize leaf fields that depend on atmosphere

   write(io6,'(/,a)') 'olam_run calling leaf3_init_atm'
   call leaf3_init_atm()

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
call rayf_init(mza,zm,zt)

!!x    call tkeinit(nza,nwa)

! If this is 'PLOTONLY' run, loop through input history files, plot 
! specified fields, and exit

if ((runtype == 'PLOTONLY') .or. (runtype == 'PARCOMBINE')) then

   write(io6,'(/,a)') 'beginning '//trim(runtype)//' loop'

   do iplt_file = 1,op%nplt_files

      hfilin = trim(op%plt_files(iplt_file))
      call history_start('COMMIO')
      call history_start('HISTREAD')

      if (runtype == 'PLOTONLY') then
         call plot_fields(0)
      else
         call history_write('STATE')
      endif

   enddo

   write(io6,'(/,a)') trim(runtype)//' run complete'
   go to 1000
endif

! REPLACE INITIAL FIELDS WITH HISTORY READ

if (trim(runtype) == 'HISTORY') then
   write(io6,*) 'olam_run calling history_start'
   call history_start('HISTREAD')
   write(io6,*) 'olam_run finished history_start'
endif

write(io6,'(/,a)') 'olam_run calling plot_fields'

call plot_fields(0)

if (trim(runtype) /= 'HISTORY') then
   write(io6,'(/,a)') 'olam_run calling history_write'
   call history_write('STATE')
   write(io6,'(/,a)') 'olam_run finished history_write'
endif

time_prevhist = time8
w2 = walltime(wtime_start)
call cpu_time(t2)
write(io6,'(a,2f12.3)') '++++++++++ CPU - wall time: olam initialization: ',t2-t1,w2-w1

write(io6,'(/,a)') 'olam_run completed initialization'

! Exit if doing a zero time run
if (time8 >= timmax8) go to 1000

! Call the model time integration driver

if(ied_offline == 0)then
   write(io6,'(/,a)') 'olam_run calling model'
   call model()
   write(io6,'(/,a)') 'subroutine model returned to olam_run'
else
   write(io6,'(/,a)') 'olam_run calling ed_offline_model'
   call ed_offline_model()
   write(io6,'(/,a)') 'subroutine ed_offline_model returned to olam_run'
endif

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

 write(io6,*) 'starting subroutine MODEL'

wtime_start = walltime(0.)

! Start the timesteps

mstp = 0

do while (time8 < timmax8)

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
   s1900_sim = s1900_init + time8

   call update_model_time(current_time, dtlm(1))

   wtime2 = walltime(wtime_start)
   call cpu_time(t2)

   write (stepc1,'(i8)'   ) mstp
   write (stepc2,'(f13.1)') time8
   write (stepc3,'(f9.2)' ) time8/86400.
   write (stepc4,'(f8.2)' ) t2-t1
   write (stepc5,'(f8.2)' ) wtime2-wtime1

   stepc1 = ' [nstep = '//trim(adjustl(stepc1))//']'
   stepc2 = '   [simtime = '//trim(adjustl(stepc2))//' sec'
   stepc3 = ' = '//trim(adjustl(stepc3))//' days]'
   stepc4 = '   [cpu,wall(sec) = '//trim(adjustl(stepc4))
   stepc5 = ' , '//trim(adjustl(stepc5))//']'
   
   write(io6,'(a)')  &
      trim(stepc1)//trim(stepc2)//trim(stepc3)//trim(stepc4)//trim(stepc5)

   call olam_output
   
enddo

wtime_tot = walltime(wtime_start)
write(io6, '(//,a,f10.0)') ' -----Total elapsed time: ',wtime_tot

return
end subroutine model

!===========================================================================

subroutine ed_offline_model()

  use misc_coms,   only: io6, time8, timmax8, dtlm, simtime, current_time, &
                         frqstate, s1900_init, s1900_sim
  use ed_options,  only: frq_phenology
  use consts_coms, only: r8
  implicit none

  real :: wtime_start
  real, external :: walltime
  integer :: mstp
  type(simtime) :: begtime
  real :: t1
  real :: wtime1
  real :: wtime2
  real :: t2
  character(len=40) :: stepc1
  character(len=40) :: stepc2
  character(len=40) :: stepc3
  character(len=40) :: stepc4
  character(len=40) :: stepc5
  real :: wtime_tot

  write(io6,*) 'starting subroutine ED_OFFLINE_MODEL'

  wtime_start = walltime(0.)

! Start the timesteps
  mstp = 0
  if (time8 < epsilon(time8)) call write_ed_output()

  do while (time8 < timmax8)

     begtime = current_time

     if(current_time%time < dtlm(1))  &
     
       write(io6,*)'Simulating:',  &
          current_time%month,'/',current_time%date,'/',current_time%year

! CPU timing information

     call cpu_time(t1)
     wtime1 = walltime(wtime_start)

     ! Get radiative fluxes
     call radiate_offline()

     ! Get sensible and latent heat fluxes, and ustar.
     call sfluxes_offline()

     ! Update the land surface
     call leaf3()

     mstp = mstp + 1
     time8 = time8 + dtlm(1)
     s1900_sim = s1900_init + time8

     call update_model_time(current_time, dtlm(1))

     if(mod(time8,real(frq_phenology,r8)) < dtlm(1))call ed_vegetation_dynamics()

     if(current_time%date == 1 .and. current_time%time < dtlm(1))then
        ! Read new met driver files only if this is the first timestep 
        ! on the first day of a month.
        call read_offline_met()
     endif

     ! Update the meterological drivers.
     call update_offline_met()

!     wtime2 = walltime(wtime_start)
!     call cpu_time(t2)

!     write (stepc1,'(i8)'   ) mstp
!     write (stepc2,'(F13.1)') time8
!     write (stepc3,'(F9.2)' ) time8/86400.
!     write (stepc4,'(F8.2)' ) t2-t1
!     write (stepc5,'(F8.2)' ) wtime2-wtime1

!     stepc1 = ' [nstep = '//trim(adjustl(stepc1))//']'
!     stepc2 = '   [simtime = '
!     stepc3 = ' = '//trim(adjustl(stepc3))//' days]'
!     stepc4 = '   [cpu,wall(sec) = '//trim(adjustl(stepc4))
!     stepc5 = ' , '//trim(adjustl(stepc5))//']'
   
!     write(io6, ('(a)'))  &
!          trim(stepc1)//trim(stepc2)//trim(stepc3)//trim(stepc4)//trim(stepc5)

     if (mod(current_time%time,frqstate) < dtlm(1))then
        call write_ed_output()
     endif

  enddo

  wtime_tot = walltime(wtime_start)
  write(io6, '(//,a,f10.0)') ' -----Total elapsed time: ',wtime_tot


  return
end subroutine ed_offline_model

!===========================================================================

subroutine olam_output()

use misc_coms,   only: io6, time8, dtlm, iflag, frqstate, timmax8, initial, &
                       s1900_sim, time_prevhist
use leaf_coms,   only: isfcl, iupdndvi, indvifile, s1900_ndvi
use sea_coms,    only: iupdsst, iupdseaice, isstfile, iseaicefile,  &
                       s1900_sst, s1900_seaice
use oplot_coms,  only: op
use mem_nudge,   only: nudflag
use isan_coms,   only: ifgfile, s1900_fg
use consts_coms, only: r8
implicit none

integer :: ierr,ifm,ifileok
real(kind=r8) :: frqplt8

frqplt8 = op%frqplt

if (mod(time8,frqplt8) < dtlm(1) .or. iflag == 1) then
   call plot_fields(0)
endif

if (mod(time8,real(frqstate,r8)) < dtlm(1)  .or.  &
   time8  >=  timmax8 - .01*dtlm(1) .or. iflag == 1) then
   call history_write('INST')
   time_prevhist = time8
endif

if (isfcl == 1 .and. iupdsst == 1) then
   if (s1900_sim >= s1900_sst(isstfile) .and. time8+.001 < timmax8) then
      call sst_database_read(1)
   endif
endif

if (isfcl == 1 .and. iupdseaice == 1) then
   if (s1900_sim >= s1900_seaice(iseaicefile) .and. time8+.001 < timmax8) then
      call seaice_database_read(1)
   endif
endif

if (isfcl == 1 .and. iupdndvi == 1) then
   if (s1900_sim >= s1900_ndvi(indvifile) .and. time8+.001 < timmax8) then
      call ndvi_database_read(1)
   endif
endif

if (initial == 2 .and. nudflag == 1) then
   if (s1900_sim >= s1900_fg(ifgfile) .and. time8+.001 < timmax8) then

      write(io6,*) ' '
      write(io6,*) 'olam_output ',time8,s1900_sim,s1900_fg(ifgfile),timmax8

      call isan_driver(1)

   endif
endif

!-------------------- SPECIAL - HURRICANE TRACKING ------------------
!if (mod(time8,real(1800.,r8)) < dtlm(1)) then
!   call hurricane_track()
!endif
!-------------------------------------------------------------------

if (iflag == 1) stop 'IFLAG'

return
end subroutine olam_output
