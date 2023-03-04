subroutine plot_trajecfile(iplt)

  ! This subroutine plots trajectories from multiple hurricane forecasts that
  ! are read from a data file originally downloaded from the National Hurricane
  ! Center.  The raw file is usually amended to include best track information
  ! and one or more OLAM simulated trajectories.  The name of the data file, the
  ! number of model simulations it contains, and a duration in hours over which
  ! to plot the trajectories must be specified below.  Additional customization
  ! of this subroutine may be needed for a large number of model simulations.

  use oplot_coms,  only: op
  use consts_coms, only: pio180, erad

  implicit none

  integer, intent(in) :: iplt

  real :: zeh, reh, xeh, yeh, bsize_hour, bsize_modelid
  character(len=3) :: title

  integer :: imodel, imodlab, ihour, icolor
  character(2) :: basin
  character(4) :: pastid, modelid
  integer :: cyclone_num, yymmddhh, technum, tau, idlat, idlon, vmax, mslp

  integer, parameter :: nmodels = 40, nhours = 50
  logical, save :: newcall = .true.

  integer, save :: kmodel

  integer,      allocatable, save :: mhour(:,:)
  integer,      allocatable, save :: numhours(:)
  real,         allocatable, save :: mlat(:,:),  mlon(:,:)
  real,         allocatable, save :: xsims(:,:), ysims(:,:)
  character(4), allocatable, save :: model_id(:,:)

  real :: xpt, ypt

  ! if (iplt > 1) return
  ! return

  if (newcall) then
     newcall = .false.

     allocate( mhour   (nmodels,nhours) )
     allocate( numhours(nmodels)        )
     allocate( mlat    (nmodels,nhours) )
     allocate( mlon    (nmodels,nhours) )
     allocate( xsims   (nmodels,nhours) )
     allocate( ysims   (nmodels,nhours) )
     allocate( model_id(nmodels,nhours) )

     ! Read locations from multiple forecast models at OLAM initialization time

     imodel = 0
     pastid = '0000'

     open(32,file='aal_olam_2022sep27_12z.dati4',status='old',form='formatted')

     do 

!STD    read(32,'(a2,2x,i2,2x,i10,2x,i2,2x,a4,2x,i3,2x,i3,3x,i4,3x,i3,2x,i4)') &  ! STD read
        read(32,'(a2,2x,i2,2x,i10,2x,i2,2x,a4,2x,i3,2x,i4,3x,i5,3x,i3,2x,i4)') &  ! lat,lon x 10
           basin, cyclone_num, yymmddhh, technum, modelid, tau, idlat, idlon, vmax, mslp

        ! Check model id

        if (basin == 'XX') exit  ! XX is our code for last line of file

        if (modelid /= pastid) then
           pastid = modelid

           imodel = imodel + 1
           ihour = 0
        endif

        ihour = ihour + 1

        numhours(imodel) = ihour

        model_id(imodel,ihour) = modelid
        mhour   (imodel,ihour) = tau
!STD    mlat    (imodel,ihour) =  0.1 * real(idlat)
!STD    mlon    (imodel,ihour) = -0.1 * real(idlon)  ! + 0.01 * real(imodel - 17) ! offset for visual separation
        mlat    (imodel,ihour) =  0.01 * real(idlat)
        mlon    (imodel,ihour) = -0.01 * real(idlon) ! + 0.01 * real(imodel - 17) ! offset for visual separation
                                                                                ! (OL01 is model 17)
        ! Find "earth" coordinates of sims hurricane center

        zeh = erad * sin(mlat(imodel,ihour) * pio180)
        reh = erad * cos(mlat(imodel,ihour) * pio180)  ! distance from earth center
        xeh = reh  * cos(mlon(imodel,ihour) * pio180)
        yeh = reh  * sin(mlon(imodel,ihour) * pio180)

        ! Transform hurricane earth coords to whatever projection is in use

        call oplot_transform(iplt,xeh,yeh,zeh,xsims(imodel,ihour),ysims(imodel,ihour))

        write(6,'(3i5,2x,a,i6,2f6.1,2f12.1)') &
           imodel, ihour, numhours(imodel), &
           model_id(imodel,ihour), mhour(imodel,ihour), &
           mlat(imodel,ihour),mlon(imodel,ihour), &
           xsims(imodel,ihour),ysims(imodel,ihour)

     enddo
     close(32)

     kmodel = imodel

  endif

  ! Set character font, line width, and plotted size

  call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
  call o_pcsetr('CL',2.)  ! line width for characters
  bsize_hour    = .010 * (op%hp2 - op%hp1)
  bsize_modelid = .012 * (op%hp2 - op%hp1)
  call o_gslwsc(2.) ! line width

  ! Plot hurricane tracks from other models (now all in one group; iplt is 1 only)

  imodlab = 0
  do imodel = 1, kmodel

     ! Set plot line color

     if     (model_id(imodel,1) == 'AVNO') then; cycle; icolor = 142 ! GFS Model Forecast
     elseif (model_id(imodel,1) == 'AC00') then; cycle; icolor = 118 ! GFS Ensemble Control Forecast
     elseif (model_id(imodel,1) == 'AEMN') then; cycle; icolor = 108 ! GFS Ensemble Mean Forecast
     elseif (model_id(imodel,1) == ' CMC') then; cycle; icolor = 137 ! Canadian Global Forecast Model
     elseif (model_id(imodel,1) == 'CEMN') then; cycle; icolor =  91 ! Canadian Ensemble Mean Forecast
     elseif (model_id(imodel,1) == 'COTC') then; cycle; icolor = 123 ! U.S. Navy COAMPS-TC Model Forecast
     elseif (model_id(imodel,1) == 'CTCX') then; cycle; icolor = 111 ! U.S. Navy Experimental COAMPS-TC Forecast
     elseif (model_id(imodel,1) == 'NVGM') then; cycle; icolor = 126 ! U.S. Navy NAVGEM Model Forecast
     elseif (model_id(imodel,1) == ' NGX') then; cycle; icolor = 113 ! U.S. Navy NOGAPS Model Forecast (deprecated)
     elseif (model_id(imodel,1) == 'EGRR') then; cycle; icolor = 121 ! UKMET/UKX Model Forecast
     elseif (model_id(imodel,1) == ' UKM') then; cycle; icolor = 134 ! UKMET Model
     elseif (model_id(imodel,1) == 'UKM2') then; cycle; icolor = 105 ! UKMET Model (interpolated 12 hours)
     elseif (model_id(imodel,1) == 'HMON') then; cycle; icolor = 120 ! Hurricanes in a Multiscale Ocean-Coupled Non-Hydrostatic...
     elseif (model_id(imodel,1) == 'HWRF') then; cycle; icolor = 143 ! Hurricane WRF
     elseif (model_id(imodel,1) == 'OFCL') then; cycle; icolor =  10 ! NHC Official Forecast 
     elseif (model_id(imodel,1) == 'BEST') then; cycle;        icolor =  97 ! BEST TRACK
     elseif (model_id(imodel,1) == 'OLAM') then;        icolor =  92 ! OLAM
     elseif (model_id(imodel,1) == 'OL01') then; cycle;        icolor =  88 ! OLAM-1
     elseif (model_id(imodel,1) == 'OL02') then; cycle;        icolor =  78 ! OLAM-2
     elseif (model_id(imodel,1) == 'OL03') then; cycle;        icolor =  73 ! OLAM-3
     elseif (model_id(imodel,1) == 'OL04') then;         icolor =  68 ! OLAM-4
     elseif (model_id(imodel,1) == 'OL05') then; cycle;        icolor =  68 ! OLAM-5
     endif

     call o_sflush()

     call o_gsplci(icolor)
     call o_gstxci(icolor)
     call o_gsfaci(icolor)

     ! Plot color code

     imodlab = imodlab + 1
     xpt = op%xmin + (op%xmax - op%xmin) * 0.15 * real((imodlab-1)/4)         ! Uses geographic coordinates in meters
     ypt = op%ymin - (op%ymax - op%ymin) * (0.133 + 0.033 * mod(imodlab-1,4)) ! that are set for a given plot frame

     call o_plchhq(xpt,ypt,model_id(imodel,1),bsize_modelid,0.,-1.)

     ! Plot results

     do ihour = 1,numhours(imodel)

!if (mod(mhour(imodel,ihour),6) /= 0) cycle
 if (mhour(imodel,ihour) < 18) cycle

        if (mhour(imodel,ihour) > 42.1)  exit    ! Plot only up to a selected time
  !     if (mhour(imodel,ihour) > 72.1)  exit    ! Plot only up to a selected time
  !     if (mhour(imodel,ihour) > 120.1) exit    ! Plot only up to a selected time

!        if (ihour == 1) then
        if (mhour(imodel,ihour) == 18) then
           call o_frstpt(xsims(imodel,ihour),ysims(imodel,ihour))
        else
           call o_vector(xsims(imodel,ihour),ysims(imodel,ihour))
        endif

        write(title,'(i3)') mhour(imodel,ihour)
        if (ihour > 1 .or. model_id(imodel,1) == 'BEST') then   ! plot 0 HRS only for BEST TRACK
           call o_plchhq (xsims(imodel,ihour),ysims(imodel,ihour),trim(adjustl(title)),bsize_hour,0.,0.)
        endif

     enddo

  enddo

end subroutine plot_trajecfile

