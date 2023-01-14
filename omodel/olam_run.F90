subroutine olam_run(name_name)

  use, intrinsic :: ieee_arithmetic

  use misc_coms,   only: io6, mstp, time8, time8p, time_istp8, time_istp8p,      &
                         iflag, expnme, mdomain, initial, iswrtyp, ilwrtyp,      &
                         runtype, hfilin, timmax8, alloc_misc, iparallel,        &
                         iyear1, imonth1, idate1, itime1, s1900_init, s1900_sim, &
                         time_prevhist, rinit, rinit8, debug_fp, init_nans,      &
                         do_chem, time_bias, dtlong, ioutput, nrk_wrtv

  use olam_mpi_atm,only: olam_mpi_atm_start, olam_mpi_atm_stop, &
                         mpi_send_w, mpi_recv_w, mpi_send_v, mpi_recv_v

  use olam_mpi_sfc,only: olam_mpi_sfc_start, olam_mpi_sfc_stop, &
                         mpi_send_wsfc, mpi_recv_wsfc

  use olam_mpi_nud,only: olam_mpi_nud_start, olam_mpi_nud_stop

  use mem_para,    only: myrank, compute_pario_points

  use leaf_coms,   only: isfcl, iupdndvi
  use sea_coms,    only: iupdsst, iupdseaice
  use mem_ijtabs,  only: istp, fill_jtabs, itab_v, itab_w
  use mem_sfcg,    only: mwsfc, alloc_sfcgrid2, filltab_sfcg, fill_jtab_sfcg
  use sea_swm,     only: swm_init, swm_diagvel
  use oplot_coms,  only: op, iplt_file
  use mem_grid,    only: mma, mva, mwa, mza, alloc_gridz_other, volvi, lpvmax
  use mem_nudge,   only: nudflag, nudnxp, fill_jnudge, o3nudflag
  use mem_rayf,    only: rayf_init
  use consts_coms, only: r8, init_consts
  use oname_coms,  only: nl
  use hcane_rz,    only: ncycle_hurrinit, icycle_hurrinit, hurricane_init, &
                         timmax_hurrinit, vortex_center_diagnose, vortex_azim_avg, &
                         vortex_reloc3d, vortex_relocated, htc0, &
                         hlat, hlon, hlat0, hlon0, hlat_hist, hlon_hist
  use obnd,        only: set_scalars_bottom, set_scalars_lbc, lbcopy_v
  use mem_plot,    only: alloc_plot, copy_plot
  use lite_vars,   only: prepare_lite, lite_write, lite_read
  use mem_addgrid, only: init_addgrid
  use mem_land,    only: land, mland
  use mem_sea,     only: sea, msea
  use vel_t3d,     only: diagvel_t3d, diag_uzonal_umerid
  use mem_adv,     only: alloc_adv
  use mem_co2,     only: co2init
  use wrtv_rk,     only: init_wrtv_rk
  use wrtv_orig,   only: init_wrtv_orig

  use cgrid_spcs,  only: cgrid_spcs_init
  use aero_data,   only: map_aero
  use emis_defn,   only: emis_init
  use depv_defn,   only: depv_init
  use mem_megan,   only: megan_init
  use cgrid_conv,  only: conv_cgrid
  use soa_defn,    only: map_soa
  use aq_data,     only: aq_data_init
  use pbl_drivers, only: pbl_init

  use precursor_data, only: map_precursor

  use mem_sfcnud,  only: sfcnud_write, alloc_sfcnud, sfcnud_read_init, read_gw_spinup

  use mem_average_vars, only: reset_davg_vars
  use mem_swtc5_refsoln_cubic

  implicit none

  character(len=*), intent(in) :: name_name

  integer :: i, idavg_file, ilite_file
  integer :: mwa_prog, mva_prog
  real :: w1,w2,t1,t2,wtime_start
  real, external :: walltime
  character(len=128) :: davgfile, litefile
  logical :: result

  wtime_start = walltime(0.)
  w1 = walltime(wtime_start)
  call cpu_time(t1)

  iflag = 0
  istp  = 1
  time8 = 0.0_r8

  ! Set remaining physical constants

  call init_consts(mdomain, nl%test_case, nl%rlat0)

  ! Check and copy namelist variables

  write(io6,'(/,a)') 'olam_run checking namelist values'
  call oname_check()

  write(io6,'(/,a)') 'olam_run calling namelist copy'
  call copy_nl()

  if (runtype == 'HISTORY' .or. runtype == 'HISTADDGRID') then
     write(io6,'(/,a)') 'olam_run reading some namelist values from history file'
     call history_start('COMMIO')
  elseif (runtype == 'PLOTONLY') then
     write(io6,'(/,a)') 'olam_run reading common values from plot file'
     hfilin = op%plt_files(1)
     call history_start('COMMIO')
  endif

  ! Initialize time variables

  time_bias   = 0.01_r8 * dtlong    ! Time bias less than the smallest acoustic timestep
  time8p      = time8 + time_bias   ! Slightly forward biased time
  time_istp8  = time8
  time_istp8p = time8p              ! Slightly forward biased time

  ! Set constants for initializing model arrays to zeros or, if supported, NANs
  ! depending on namelist variable init_nans

  rinit  = 0.0
  rinit8 = 0.0_r8

  if (init_nans) then
     if (ieee_support_nan(1.0)) then
        rinit = ieee_value(1.0, ieee_signaling_nan)
     endif
  endif

  if (init_nans) then
     if (ieee_support_nan(1.0_r8)) then
        rinit8 = ieee_value(1.0_r8, ieee_signaling_nan)
     endif
  endif

  ! If debugging, halt on illegal floating operations if supported

  if (ieee_support_halting(ieee_invalid) .and. debug_fp) then
     call ieee_set_halting_mode(ieee_usual, .true.)
  endif

  ! Get abs seconds of simulation start and current simulation time
  ! Note that itime1 has format of hhmm, and date_abs_secs2 needs hhmmss

  call date_abs_secs2(iyear1,imonth1,idate1,itime1*100,s1900_init)
  s1900_sim = s1900_init + time8

  ! Print initial banner

  write(io6,'(a1)')         ' '
  write(io6,'(a1,78a1)')    ' ',('*',i=1,78)
  write(io6,'(2a1,a42)')    ' ','*','    OLAM version 6.0'
  write(io6,'(2a1)')        ' ','*'
  write(io6,'(2a1,a3,a64)') ' ','*','   ',EXPNME
  write(io6,'(a1,78a1)')    ' ',('*',i=1,78)

  ! MAKEGRID runs must be single-processor

  if ( runtype == 'MAKEGRID' .or. runtype == 'MAKEADDGRID' .or. runtype == 'MAKEGRID_PLOT' ) then
     if (iparallel == 1) then
        write(io6,*) trim(runtype)//' will only be done on a single process.'
        iparallel  = 0
        if (myrank > 0) go to 1000
     endif
  endif

  ! Initialize various oplot parameters

  call oplot_init()

  ! If RUNTYPE = 'MAKEGRID', generate full-domain ATM and SFC grids
  ! If RUNTYPE = 'MAKEADDGRID', generate full-domain ATM grid and use OLD SFC grid

  if (runtype == 'MAKEGRID' .or. runtype == 'MAKEADDGRID' .or. runtype == 'MAKEGRID_PLOT') then
     call gridinit()

     call gridset_print()

     if (runtype /= 'MAKEGRID_PLOT') call plot_fields(0)

     write(io6,*)
     write(io6,*) trim(runtype) // ' run complete'
     go to 1000
  endif

  ! If land/sea models are active, read SFCGFILE iw points

  if (mdomain <= 1) then
     write(io6,*)
     write(io6,*) 'olam_run calling sfcgfile_read_pd'
     call sfcgfile_read_pd()
  endif

  ! Read from GRIDFILE the fields that are needed by para_decomp

  write(io6,*)
  write(io6,*) 'olam_run calling gridfile_read_pd'

  call gridfile_read_pd()

  ! If run is parallel, assign each grid cell (ATM and SFCG)
  ! to one of multiple subdomains. If not, the grid remains unchanged

  write(io6,*)
  write(io6,*) 'olam_run calling para_decomp'

  call para_decomp()

  ! Set up itab data types and grid coordinate arrays for current node, and
  ! reallocate memory for current node, and finish gridfile reads

  write(io6,*)
  write(io6,*) 'olam_run calling init_atmgrid'

  call init_atmgrid()

  if (mdomain <= 1) then
     write(io6,*)
     write(io6,*) 'olam_run calling init_sfcgrid'
     call init_sfcgrid()
  endif

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
     write(io6,*)
     write(io6,*) 'olam_run calling init_nudgrid'
     call init_nudgrid()
  endif

  ! Select points that are writen to disk for parallel I/O (each
  ! point can only be written once even if defined on multiple nodes)

  call compute_pario_points()

  call grav_init()

  write(io6,*)
  write(io6,*) 'olam_run calling alloc_gridz_other'

  call alloc_gridz_other()

  call gridset_print()

  if (mdomain <= 1) call landgrid_print()

  mwa_prog = 0
  do i=1,mwa
     if (itab_w(i)%irank == myrank) mwa_prog = mwa_prog + 1
  enddo
  write(io6,*)
  write(io6,'(a,i8)') ' # of prognostic W points on this node = ', mwa_prog

  mva_prog = 0
  do i=1,mva
     if (itab_v(i)%irank == myrank) mva_prog = mva_prog + 1
  enddo
  write(io6,'(a,i8)') ' # of prognostic V points on this node = ', mva_prog

  write(io6,'(/,a)' ) 'Local model indices:'
  write(io6,'(a,i8)') ' mma = ',mma
  write(io6,'(a,i8)') ' mva = ',mva
  write(io6,'(a,i8)') ' mwa = ',mwa

  if (mdomain <= 1) then
     write(io6,'(a,i8)') ' mwsfc = ', mwsfc
     write(io6,'(a,i8)') ' mland = ', mland
     write(io6,'(a,i8)') ' msea  = ', msea
  endif

  ! Initialize long and short timesteps, and compute the timestep schedule for all operations

  write(io6,'(/,a)') 'olam_run calling modsched'

  call modsched()

  write(io6,'(/,a)') 'olam_run calling fill_jtabs'

  call fill_jtabs(mma,mva,mwa)

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) call fill_jnudge()

  ! Allocate and fill jtab_sfcg data structure

!!  if (isfcl == 1) then
     write(io6,'(/,a)') 'olam_run calling fill_jtab_sfcg'
     call fill_jtab_sfcg(mwa)
!!  endif

  ! Allocate column initial state arrays

  call alloc_misc(mza)

  write(io6,'(/,a)') 'olam_run calling jnmbinit and micinit_tabs'

  call jnmbinit()
  call micinit_tabs()

  ! Setup CMAQ chemical species

  if (do_chem == 1) then
     result = cgrid_spcs_init()
     call hrinit()
     call map_aero()
     call map_soa()
     call map_precursor()
     call aq_data_init()
  endif

  ! Allocate remainder of main model memory (for 'INITIAL', 'HISTORY',
  ! or 'PLOTONLY' run).
  ! Allocate variable tables and fill variable tables for history files
  ! and parallel communication.

  write(io6,'(/,a)') 'olam_run calling olam_mem_alloc'

  call olam_mem_alloc()

  if (iparallel == 1) then

     write(io6,'(/,a)') 'olam_run calling olam_mpi_atm_start'
     call olam_mpi_atm_start()

     write(io6,'(/,a)') 'olam_run calling olam_mpi_sfc_start'
     call olam_mpi_sfc_start()

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        write(io6,'(/,a)') 'olam_run calling olam_mpi_sfc_start'
        call olam_mpi_nud_start()
     endif

  endif

  if (iparallel == 1) then
     call mpi_send_v(rvara1=volvi, i1dvara1=lpvmax)
     call mpi_recv_v(rvara1=volvi, i1dvara1=lpvmax)
  endif
  call lbcopy_v(vmc=volvi, iv1=lpvmax)

  ! Allocate memory for advection and pre-compute arrays for 3rd order advection
  ! (needs MPI to have been initialized)

  call alloc_adv()

  ! Initialize 3d microphysics fields (if miclevel = 3) and other microphysics
  ! quantities

  write(io6,'(/,a)') 'olam_run calling micinit_fields'

  call micinit_fields()

  ! Initialize spectral nudging

  if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) call init_spec_nudge()

  ! Initialize primary atmospheric fields

  if (initial == 1) then
     if (runtype == 'INITIAL') then
        write(io6,'(/,a)') 'olam_run calling inithh'
        call inithh()            ! Horizontally-homogeneous initialization
     endif
  elseif (initial == 2) then
     if ((runtype == 'INITIAL') .or. (nudflag == 1 .or. o3nudflag == 1)) then
        write(io6,'(/,a)') 'olam_run calling isan driver(0)'
        call isan_driver(0)
     endif
  elseif (initial == 3) then
     if (runtype == 'INITIAL') then
        write(io6,'(/,a)') 'olam_run calling fldslhi'
        call fldslhi()           ! Longitudinally-homogeneous initialization
     endif
  endif

  !------------------------------------------------------------
  ! Call fill_swtc5 to read in and initialize reference
  ! solution for time = 0 and time = 15d.

  if (nl%test_case == 2 .or. nl%test_case == 5) then
     call swtc_init()
  endif
  !------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Call NCAR interface with zero (i.e., initial) time.  This initializes
  ! both momentum and scalars

  if (nl%test_case > 7 .and. nl%test_case < 500) call olam_dcmip_init()
  !----------------------------------------------------------------------

  ! Initial diagnosis of vxe,vye,vze

  if (runtype == 'INITIAL') then
     call diagvel_t3d()
  endif

  ! Initialize cloud fraction

  if (runtype == 'INITIAL') then
     call calc_3d_cloud_fraction()
  endif

  !----------------------------------------------------------------------
  ! For NCAR DCMIP supercell test case, save initial velocity and thil fields

  if (nl%test_case == 131) call dcmip_save_initfields()
  !----------------------------------------------------------------------

  ! Initialize CMAQ chemical species

  if (do_chem == 1 .and. runtype == 'INITIAL') then
     write(io6,'(/,1x,a)') 'Initializing chemical concentrations'
     call init_cgrid()
  endif

  ! A good place to initialize added scalars

  ! Initialize 3d CO2 field (if co2flg = 1)

  write(io6,'(/,a)') 'olam_run calling co2init'

  call co2init()

  ! Set dummy values for scalars below lpw

  call set_scalars_bottom()

  ! For parallel run, send and receive initialized scalars

  if (iparallel == 1) then
     call mpi_send_w(scalars='S')  ! Send scalars
     call mpi_recv_w(scalars='S')  ! Recv scalars
  endif

  ! Lateral boundary copy of scalars for limited-area run

  call set_scalars_lbc()

  ! Start up radiation scheme

  if (iswrtyp > 0 .or. ilwrtyp > 0) then
     write(io6,'(/,a)') 'olam_run calling radinit'
     call radinit()
  endif

  ! Check if LEAF will be used

  if (isfcl == 1) then

     ! Allocate time-dependent SFC grid arrays

     write(io6,'(/,a)') 'olam_run calling alloc_sfcgrid2'
     call alloc_sfcgrid2(mwsfc)
     call filltab_sfcg()

     ! Start up LAND model

     write(io6,'(/,a)') 'olam_run calling land_startup'
     call land_startup()

     ! Start up LAKE model

     write(io6,'(/,a)') 'olam_run calling lake_startup'
     call lake_startup()

     ! Start up SEA model

     write(io6,'(/,a)') 'olam_run calling sea_startup'
     call sea_startup()

     ! Start up POM1D model

     write(io6,'(/,a)') 'olam_run calling pom_startup'
     call pom_startup()

     write(io6,'(/,a)') 'olam_run calling swm_init'
     call swm_init()

     ! Average atmospheric fields to SFC grid cells

     write(io6,'(/,a)') 'olam_run calling sfcg_avgatm'
     call sfcg_avgatm()

     ! Initialize leaf fields

     write(io6,'(/,a)') 'olam_run calling leaf4_init_atm'
     call leaf4_init_atm()

     ! Initialize lake fields

     write(io6,'(/,a)') 'olam_run calling lake_init_atm'
     call lake_init_atm()

     ! Initialize ocean fields

     write(io6,'(/,a)') 'olam_run calling sea_init_atm'
     call sea_init_atm()

     ! Initialize POM1D fields

     write(io6,'(/,a)') 'olam_run calling pom_init'
     call pom_init()

     if (isfcl == 1) call swm_diagvel()

     if (nl%igw_spinup == 1) then

        ! If making or using surface nudging files for SOIL model spin-up, allocate memory

        call alloc_sfcnud()

        ! If using surface nudging files, initialize file info and read current file

        if (runtype == 'INITIAL' .or. runtype == 'HISTORY' .or. runtype == 'HISTADDGRID') then
           call sfcnud_read_init()
        endif

     elseif (nl%igw_spinup == 2) then

        ! If initializing soil water and energy from spin-up simulation results,
        ! read results and apply to initial condition

        call read_gw_spinup()

     endif

     ! MPI send/recv of time-dependent SFCG, LAND, LAKE, SEA quantities

     if (iparallel == 1) then
        call mpi_send_wsfc(set='sfc_driv_end')
        call mpi_recv_wsfc(set='sfc_driv_end')
        ! Also, send/recv vsfc?
     endif

  endif

  ! If using variable initialization and polygon nudging, read most recent
  ! and next observational analyses, and fill nudging polygons for both

  if (initial == 2 .and. (nudflag == 1 .or. o3nudflag == 1)) then
     write(io6,'(/,a)') 'olam_run calling isan_driver(1)'
     call isan_driver(1)
  endif

  ! Initialize Rayleigh friction profile

  write(io6,'(/,a)') 'olam_run calling rayf_init'
  call rayf_init()

  ! For shallow water test case 5, read in and initialize reference
  ! solution for time = 0 and time = 15d.

  if (nl%test_case == 5) then
     call fill_swtc5()
  endif

  ! Initialize PBL quantities

  call pbl_init()

  ! Extra initializations for small-timestep solver

  if (nrk_wrtv == 1) then
     call init_wrtv_orig()
  else
     call init_wrtv_rk()
  endif

  ! Initialize emissions/deposition if doing chemistry

  if (do_chem == 1) then
     write(io6,'(/,1x,a)') 'Initializing chemical emissions/deposition'
     call emis_init()
     call depv_init()

     if (isfcl > 0) then
        call megan_init()
     endif
  endif

  ! Memory for storing past values for plotting

  call alloc_plot()

  ! If this is 'PLOTONLY' run, loop through input history files, plot
  ! specified fields, and exit

  if (runtype == 'PLOTONLY') then

     write(io6,'(/,a)') 'beginning '//trim(runtype)//' loop'

     do iplt_file = 1,op%nplt_files

        hfilin = op%plt_files(iplt_file)
        call history_start('COMMIO')
        call history_start('HISTREAD')

        ! If day-average plots are specified, do them when
        ! iplt_file = 1 and then exit iplt_file do loop.  IT IS ASSUMED THAT
        ! DAY-AVERAGE FILES ARE SINGLE NON-PARALLEL FILES.

        if (nl%ioutput_davg == 1 .and. nl%ndavg_files > 0) then
           do idavg_file = 1,nl%ndavg_files
              davgfile = nl%davg_files(idavg_file)
              call read_davg_vars(davgfile)
              call plot_fields(0)
           enddo
           exit
        endif

        if (nl%ioutput_lite == 1 .and. nl%nlite_files > 0) then
           do ilite_file = 1,nl%nlite_files
              litefile = nl%lite_files(ilite_file)
              call lite_read(litefile)
            !  call plot_fields(0)
              call timeseries_plots('L')
           enddo
           exit
        endif

        ! If day-average plots are NOT specified, do regular
        ! plots from history files.  First, save a copy of some fields; used
        ! later to plot difference fields

        call copy_plot(iplt_file)

        ! The following commented-out IF block is a template for cases where
        ! fields need not be plotted for every call to copy_plot

        ! if (iplt_file == 2) then
        ! if (iplt_file == 20 .or. &
        !     iplt_file == 30 .or. &
        !     iplt_file == 40 .or. &
        !     iplt_file == 50) then

       ! if (mod(iplt_file,240) == 0) then

           call plot_fields(0)

           if (nl%ioutput_latlon == 1 .or. nl%latlonplot == 1) then
              call diag_uzonal_umerid()
              call fields3_ll()
           endif

        ! endif

        ! Write sfcnud file to be used for nudged spin-up simulation

        if (nl%igw_spinup == 1 .and. isfcl == 1) then
           ! if (mod(iplt_file,5) == 0) then    ! Template: modify for # of years in climate simulation
              call sfcnud_write()
           ! endif
        endif

        ! For shallow water test cases, compute error norms

        if (nl%test_case == 2 .or. nl%test_case == 5) then
           call diagn_global_swtc()
        endif

        ! call timeseries_plots('H')

     enddo

     write(io6,'(/,a)') trim(runtype) //' run complete'
     go to 1000

  endif

  ! If adding new grids (new mesh refinements), read information from old
  ! gridfile and sfcgfile, and initialize procedure for mapping
  ! model fields from OLD history file to NEW model grid

  if (runtype == 'HISTADDGRID') then
     call gridfile_read_oldgrid()

     print*, 'calling init_addgrid'
     call init_addgrid()
  endif

  icycle_hurrinit = 0

  if (ncycle_hurrinit > 0) call hurricane_init()

  70 continue ! HURRICANE INITIALIZATION CYCLE BEGINS HERE

  if (ncycle_hurrinit > 0) then
     icycle_hurrinit = icycle_hurrinit + 1
  endif

  ! For HISTORY start or for second and later cycles of hurricane dynamic
  ! initialization, replace initial fields with HISTORY read.

  ! (Trying to history restart while hurricane initialization cycles are in
  ! progress does not work correctly.)

  if (runtype == 'HISTORY' .or. runtype == 'HISTADDGRID' .or. &
     icycle_hurrinit > 1) then

     if (runtype == 'INITIAL' .and. icycle_hurrinit > 1) then
        hfilin = trim(htc0)
     else
        hfilin = nl%hfilin
     endif

     write(io6,*) 'olam_run calling history_start'
     call history_start('COMMIO')
     call history_start('HISTREAD')
     write(io6,*) 'olam_run finished history_start'

     ! Reset time variables

     time8p      = time8 + time_bias   ! Slightly forward biased time
     time_istp8  = time8
     time_istp8p = time8p              ! Slightly forward biased time

     mstp = 0

     ! reset convective variables if we have turned off convection

     call reset_cuparm()

     ! If not updating SST/SEAICE/NDVI, copy the curent values to
     ! the past and future arrays

     if (iupdsst /= 1) then
        sea%seatp(:) = sea%seatc(:)
        sea%seatf(:) = sea%seatc(:)
     endif

     if (iupdseaice /= 1) then
        sea%seaicep(:) = sea%seaicec(:)
        sea%seaicef(:) = sea%seaicec(:)
     endif

     if (iupdndvi /= 1) then
        land%veg_ndvip(:) = land%veg_ndvic(:)
        land%veg_ndvif(:) = land%veg_ndvic(:)
     endif

  endif

  ! If this is not a history start AND if it is not the second or later
  ! cycle of hurricane initialization, write initial history file

  if (runtype /= 'HISTORY' .and. icycle_hurrinit < 2) then
     if (icycle_hurrinit == 1 .and. ncycle_hurrinit > 1) then
        write(io6,'(/,a)') 'olam_run calling history_write with HTC0 vtype'
        call history_write('HTC0')
        write(io6,'(/,a)') 'olam_run finished history_write'
     else if (ioutput /= 0) then
        write(io6,'(/,a)') 'olam_run calling history_write with STATE vtype'
        call history_write('STATE')
        write(io6,'(/,a)') 'olam_run finished history_write'
     endif
  endif

  ! Initialize cloud fraction in case there is any initial saturation

  if (runtype == 'INITIAL') then
     call calc_3d_cloud_fraction()
  endif

  time_prevhist = time8
  w2 = walltime(wtime_start)
  call cpu_time(t2)
  write(io6,'(a,2f13.3)') '++++++++++ CPU - wall time: olam initialization: ',t2-t1,w2-w1

  write(io6,'(/,a)') 'olam_run completed initialization'

  if (icycle_hurrinit == ncycle_hurrinit) then ! Includes case where ncycle_hurrinit = 0

     ! Setup lite variable output

     if (nl%ioutput_lite == 1) then
        call prepare_lite()
        if (runtype /= 'HISTORY' .and. runtype /= 'HISTADDGRID') call lite_write()
     endif

     ! Initialize field average arrays

     call reset_davg_vars() ! call unconditionally because used in accum

  endif

  ! If hurricane dynamic initialization or tracking is to be done...

  if (ncycle_hurrinit > 0) then

     ! Initialize hurricane location

     if (runtype == 'INITIAL' .and. icycle_hurrinit == 1) then
        hlat = hlat0      ! Value specified in namelist
        hlon = hlon0      ! Value specified in namelist
        call vortex_center_diagnose()
     elseif (runtype == 'HISTORY' .or. runtype == 'HISTADDGRID') then
        hlat = hlat_hist  ! Current value in history file
        hlon = hlon_hist  ! Current value in history file
     endif

     ! In 'INITIAL' run, on second or later cycle, relocate hurricane fields
     ! (that were prognosed on previous cycle) back to initial location

     if (runtype == 'INITIAL' .and. icycle_hurrinit > 1) then
        print*, 'olam_run calling vortex_relocated'

        call vortex_relocated()

        print*, 'returned from vortex_relocated'

     endif

     ! If INITIAL runtype, write secondary history file with 'HTC' in name.

     if (runtype == 'INITIAL' .and. ncycle_hurrinit > 1) then
        write(io6,'(/,a)') 'olam_run calling history_write with HTC1 vtype'
        call history_write('HTC1')
        write(io6,'(/,a)') 'olam_run finished history_write with HTC1 vtype'
     endif

  endif

  ! Plot initial fields

  write(io6,'(/,a)') 'olam_run calling plot_fields'

  call copy_plot(0)
  call plot_fields(0)
  call fields3_ll()
  if (mdomain == 5 .and. nl%les_diag_freq > 0._r8) call les_diag()

  ! Compute azimuthal averages of dynamic, thermodynamic, and moisture
  ! fields, and plot averages (in radial-height cross sections)

  if (ncycle_hurrinit > 0) then
     call vortex_azim_avg('plot')
  endif

  ! Exit if doing a zero time run
  if (time8 >= timmax8) go to 1000

  ! Call the model time integration driver

  write(io6,'(/,a)') 'olam_run calling model'
  call model()
  write(io6,'(/,a)') 'subroutine model returned to olam_run'

  ! PREPARE FOR NEXT HURRICANE DYNAMIC INITIALIZATION CYCLE (if
  ! any more are to be done)

  if (icycle_hurrinit > 0 .and. icycle_hurrinit < ncycle_hurrinit) then

     ! Remap 3D hurricane fields from grid cells at present location to
     ! grid cells at initial location; store in arrays

     call vortex_reloc3d()

     ! Return to starting point of next initialization cycle

     goto 70

  endif

  ! OLAM finished, clean up some last things...

  1000 continue

  write(io6,'(/,a)') 'olam_run calling o_clsgks'

  call o_clsgks()

  if (iparallel == 1) then
     write(io6,'(/,a)') 'olam_run closing MPI communication buffers'

     call olam_mpi_atm_stop()

     call olam_mpi_sfc_stop()

     if (mdomain == 0 .and. nudflag > 0 .and. nudnxp > 0) then
        call olam_mpi_nud_stop()
     endif

  endif

  write(io6,'(/,a)') 'olam_run returning to olammain'

end subroutine olam_run

!=====================================================================

subroutine model()

  use misc_coms, only: io6, mstp, time8, time8p, time_istp8, time_istp8p, &
                       time_bias, timmax8, dtlm, simtime, current_time, &
                       s1900_init, s1900_sim
  use consts_coms, only: r8
  use oname_coms,  only: nl
  use hcane_rz,    only: ncycle_hurrinit, icycle_hurrinit, timmax_hurrinit

  implicit none

  !   +------------------------------------------------------------------
  !   ! This routine drives the entire time integration process
  !   !   for a non-parallel run.
  !   +------------------------------------------------------------------
  !

  real :: wtime_start, t1, wtime1, wtime2, t2, wtime_tot
  real, external :: walltime
  real(r8) :: timmax8_model
  character(len=40) :: stepc1,stepc2,stepc3,stepc4,stepc5
  type(simtime) :: begtime

  write(io6,*) 'starting subroutine MODEL'

  wtime_start = walltime(0.)

  if (ncycle_hurrinit > 0 .and. icycle_hurrinit < ncycle_hurrinit) then
     timmax8_model = timmax_hurrinit
  else
     timmax8_model = timmax8
  endif

  ! Start the timesteps

  mstp = 0

  do while (time8p < timmax8_model)

     begtime = current_time

     ! CPU timing information

     call cpu_time(t1)
     wtime1 = walltime(wtime_start)

     ! Start the timestep schedule to loop through all grids and advance them
     ! in time an increment equal to dtlm.

     call timestep()

     mstp = mstp + 1
     time8       = time8 + dtlm
     time8p      = time8 + time_bias   ! Slightly forward biased time
     time_istp8  = time8
     time_istp8p = time8p              ! Slightly forward biased time

     s1900_sim = s1900_init + time8

     call update_model_time(current_time, dtlm)

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

     if (nl%ioutput_davg == 1) call inc_davg_vars()

     ! Check schedule for I/O operations and perform those that are due

     call olam_output()

  enddo

  wtime_tot = walltime(wtime_start)
  write(io6, '(//,a,f13.3)') ' -----Total elapsed time: ',wtime_tot

end subroutine model

!===========================================================================

subroutine olam_output()

  use misc_coms,   only: io6, time8, time8p, dtlm, iflag, frqstate, timmax8, &
                         initial, s1900_sim, time_prevhist, mdomain, &
                         iyear1, imonth1, idate1, itime1, do_chem
  use leaf_coms,   only: isfcl, iupdndvi, indvifile, s1900_ndvi
  use sea_coms,    only: iupdsst, iupdseaice, isstfile, iseaicefile, &
                         s1900_sst, s1900_seaice, isstflg, iseaiceflg
  use oplot_coms,  only: op
  use mem_nudge,   only: nudflag, o3nudflag
  use isan_coms,   only: ifgfile, s1900_fg
  use consts_coms, only: r8
  use oname_coms,  only: nl
  use mem_plot,    only: copy_plot
  use hcane_rz,    only: ncycle_hurrinit, icycle_hurrinit, vortex_center_diagnose, &
                         vortex_azim_avg
  use lite_vars,   only: lite_write
  use mem_megan,   only: megan_store_lai

  use mem_average_vars, only: reset_davg_vars
  use mem_sfcnud,  only: s1900_sfcnud, isfcnudfile, sfcnud_read

  implicit none

  integer :: outyear, outmonth, outdate, outhour

  ! Track location of hurricane every timestep

  if (icycle_hurrinit > 0) then
     call vortex_center_diagnose()
  endif

  ! Save a copy of some fields; used later to compute & plot difference fields,
  ! and plot fields

  if (mod(time8p,op%frqplt) < dtlm .or. time8p >= timmax8 .or. iflag == 1) then
     call copy_plot(0)
     call plot_fields(0)

     if (ncycle_hurrinit > 0) then
        call vortex_azim_avg('plot')
     endif
  endif

  call date_add_to8(iyear1,imonth1,idate1,itime1,time8p,'s',outyear,  &
                    outmonth,outdate,outhour)

  ! Output full history restart file

  if (mod(time8p,frqstate) < dtlm   .or. &
     (outdate == 1 .and. mod(time8p,86400.0_r8) < dtlm .and. nl%igw_spinup /= 1) .or. &
     time8p >= timmax8 .or. iflag == 1) then
     call history_write('INST')
     time_prevhist = time8
  endif

  ! Ouput lat/lon interpolated quantities

  if (nl%ioutput_latlon == 1 .and. mod(time8p,nl%frqlatlon) < dtlm) then
     call fields3_ll()
  endif

  ! Output of "lite" quantities

  if (nl%ioutput_lite == 1 .and. mod(time8p,nl%frqlite) < dtlm) then
     call lite_write()
  endif

  ! Output of time-averaged quantities

  if (mod(time8p,86400.0_r8) < dtlm) then
     call date_add_to8(iyear1,imonth1,idate1,itime1,time8p,'s',outyear,  &
          outmonth,outdate,outhour)

     if (nl%ioutput_davg == 1) then
        call norm_davg_vars()
        call write_davg_vars(outyear,outmonth,outdate)
     endif
     call reset_davg_vars() ! Always call because used in subroutine flux_accum

  endif

  ! Print LES area-averaged statistics

  if (mdomain == 5 .and. nl%les_diag_freq > 1.e-5_r8) then
     if (mod(time8p,nl%les_diag_freq) < dtlm) then
        call les_diag()
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

  ! If using surface nudging files for SOIL model spin-up, read next file at proper time

  if (isfcl == 1 .and. nl%igw_spinup == 1) then
     if (s1900_sim >= s1900_sfcnud(isfcnudfile) .and. time8p < timmax8) then
        call sfcnud_read(1)
     endif
  endif

  if (isfcl == 1 .and. iupdndvi == 1) then
     if (s1900_sim >= s1900_ndvi(indvifile) .and. time8p < timmax8) then
        call ndvi_database_read(1)
        if (do_chem == 1) call megan_store_lai()
     endif
  endif

  if (initial == 2 .and. (nudflag == 1 .or. o3nudflag == 1)) then
     if (s1900_sim >= s1900_fg(ifgfile) .and. time8p < timmax8) then

        write(io6,*) ' '
        write(io6,*) 'Calling isan_driver for next analysis file.'
        write(io6,*)  time8, s1900_sim, s1900_fg(ifgfile), timmax8

        call isan_driver(1)

     endif
  endif

  if (iflag == 1) stop 'IFLAG'

end subroutine olam_output
