subroutine plot_harvey(iplt)

  use mem_basic,   only: press
  use mem_grid,    only: mwa, lpw, glatw, glonw
  use misc_coms,   only: time8
  use oplot_coms,  only: op
  use consts_coms, only: pio180, erad

  implicit none

  integer :: iplt

  integer :: idat,num,ifhour,nfhour,iwlp,iw
  logical, save :: newcall = .true.
  real, save, dimension(0:300) :: flat,flon
  real :: zeh,reh,xeh,yeh,press_lowest,bsize
  character(len=2) :: title

  integer :: iline, imodel, jmodel, ihour, im1, im2, icolor
  character(2) :: basin
  character(4) :: pastid, modelid
  integer :: cyclone_num, yymmddhh, technum, tau, idlat, idlon, vmax, mslp

  integer, parameter :: nmodels = 200, nhours = 50, nolam = 3

  integer, save :: mhour(nmodels,nhours)
  integer, save :: numhours(nmodels)
  real, save ::  mlat(nmodels,nhours),  mlon(nmodels,nhours)
  real,save :: xsims(nmodels,nhours), ysims(nmodels,nhours)
  character(4), save :: model_id(nmodels,nhours)

  real :: xpt1, xpt2, xpt3, ypt

  real, save, dimension(100) :: lonsim1,latsim1,lonsim2,latsim2
  real, save, dimension(100) :: xobs,yobs,xsim1,ysim1,xsim2,ysim2

    RETURN  ! Do this return for main 2-day run or normal plotonly run,
            ! and comment out this return when doing spaghetti plots

  ! if (iplt > 1) return

  if (newcall) then
     newcall = .false.

     ! Read Harvey location

     open(32,file='harvey_location_24aug2017_12UTC',status='old',form='formatted')
  !  open(32,file='harvey_location_27aug2017_12UTC',status='old',form='formatted')
  !  open(32,file='harvey_location_28aug2017_00UTC',status='old',form='formatted')
     do idat = 0,17
        read(32,'(i4,2f8.2)') num,flat(idat),flon(idat)
     enddo
     close(32)

     ! Read locations from multiple forecast models at OLAM initialization time

     imodel = 0
     pastid = '0000'

     open(32,file='aug24_12z.dat2_olamA',status='old',form='formatted')
     do iline = 1,462  ! Number of lines in file

        read(32,'(a2,2x,i2,2x,i10,2x,i2,2x,a4,2x,i3,2x,i3,3x,i4,3x,i3,2x,i4)') &
           basin, cyclone_num, yymmddhh, technum, modelid, tau, idlat, idlon, vmax, mslp

        ! Check model id

        if (modelid /= pastid) then
           pastid = modelid

           imodel = imodel + 1
           ihour = 0
        endif

        ihour = ihour + 1

        numhours(imodel) = ihour

        model_id(imodel,ihour) = modelid
        mhour   (imodel,ihour) = tau
        mlat    (imodel,ihour) =  0.1 * real(idlat)
        mlon    (imodel,ihour) = -0.1 * real(idlon)

        ! Find "earth" coordinates of sims hurricane center

        zeh = erad * sin(mlat(imodel,ihour) * pio180)
        reh = erad * cos(mlat(imodel,ihour) * pio180)  ! distance from earth center
        xeh = reh  * cos(mlon(imodel,ihour) * pio180)
        yeh = reh  * sin(mlon(imodel,ihour) * pio180)

        ! Transform hurricane earth coords to whatever projection is in use

        call oplot_transform(iplt,xeh,yeh,zeh,mlon(imodel,ihour),mlat(imodel,ihour), &
                             xsims(imodel,ihour),ysims(imodel,ihour))

        write(6,'(3i5,2x,a,i6,2f6.1,2f12.1)') &
           imodel, ihour, numhours(imodel), &
           model_id(imodel,ihour), mhour(imodel,ihour), &
           mlat(imodel,ihour),mlon(imodel,ihour), &
           xsims(imodel,ihour),ysims(imodel,ihour)

     enddo
     close(32)

! BASIN, CYCLONE NUMBER, YYYYMMDDHH, TECHNUM/MIN, TECH, TAU (HOUR), LAT, LON, VMAX, MSLP, TY, ...

     ! Assign hurricane locations for sim 1 (build_harvey5_534)

     lonsim1( 0) =  -93.10 ; latsim1( 0) =  23.70 ! observed; use this
  !  lonsim1( 0) =  -92.96 ; latsim1( 0) =  23.93 ! based on CFSR; do not use
     lonsim1( 1) =  -93.80 ; latsim1( 1) =  24.27
     lonsim1( 2) =  -94.44 ; latsim1( 2) =  24.79
     lonsim1( 3) =  -95.38 ; latsim1( 3) =  25.65
     lonsim1( 4) =  -96.20 ; latsim1( 4) =  26.25
     lonsim1( 5) =  -96.94 ; latsim1( 5) =  27.04
     lonsim1( 6) =  -97.66 ; latsim1( 6) =  27.54
     lonsim1( 7) =  -98.26 ; latsim1( 7) =  27.82
     lonsim1( 8) =  -98.74 ; latsim1( 8) =  28.18
     lonsim1( 9) =  -99.42 ; latsim1( 9) =  28.40
     lonsim1(10) = -100.23 ; latsim1(10) =  28.50
     lonsim1(11) = -100.76 ; latsim1(11) =  28.34
     lonsim1(12) = -101.55 ; latsim1(12) =  28.33
     lonsim1(13) = -102.27 ; latsim1(13) =  28.08
     lonsim1(14) =  0.; latsim1(14) =   0.
     lonsim1(15) =  0.; latsim1(15) =   0.
     lonsim1(16) =  0.; latsim1(16) =   0.
     lonsim1(17) =  0.; latsim1(17) =   0.
     lonsim1(18) =  0.; latsim1(18) =   0.
     lonsim1(19) =  0.; latsim1(19) =   0.
     lonsim1(20) =  0.; latsim1(20) =   0.
     lonsim1(21) =  0.; latsim1(21) =   0.

     ! Assign hurricane locations for sim 2 (build_harvey5B_534)

     lonsim2( 0) =  -93.10 ; latsim1( 0) =  23.70 ! observed; use this
  !  lonsim2( 0) =  -92.96 ; latsim2( 0) =  23.93 ! based on CFSR; do not use
     lonsim2( 1) =  -93.76 ; latsim2( 1) =  24.28
     lonsim2( 2) =  -94.42 ; latsim2( 2) =  24.81
     lonsim2( 3) =  -95.48 ; latsim2( 3) =  25.72
     lonsim2( 4) =  -96.29 ; latsim2( 4) =  26.26
     lonsim2( 5) =  -96.91 ; latsim2( 5) =  27.03
     lonsim2( 6) =  -97.57 ; latsim2( 6) =  27.47
     lonsim2( 7) =  -98.24 ; latsim2( 7) =  27.85
     lonsim2( 8) =  -98.86 ; latsim2( 8) =  28.13
     lonsim2( 9) =  -99.34 ; latsim2( 9) =  28.30
     lonsim2(10) =  0.; latsim2(10) =   0.
     lonsim2(11) =  0.; latsim2(11) =   0.
     lonsim2(12) =  0.; latsim2(12) =   0.
     lonsim2(13) =  0.; latsim2(13) =   0.
     lonsim2(14) =  0.; latsim2(14) =   0.
     lonsim2(15) =  0.; latsim2(15) =   0.
     lonsim2(16) =  0.; latsim2(16) =   0.
     lonsim2(17) =  0.; latsim2(17) =   0.
     lonsim2(18) =  0.; latsim2(18) =   0.
     lonsim2(19) =  0.; latsim2(19) =   0.
     lonsim2(20) =  0.; latsim2(20) =   0.
     lonsim2(21) =  0.; latsim2(21) =   0.

     do idat = 0,12

        ! Find "earth" coordinates of observed hurricane center

        zeh = erad * sin(flat(idat) * pio180)
        reh = erad * cos(flat(idat) * pio180)  ! distance from earth center
        xeh = reh  * cos(flon(idat) * pio180)
        yeh = reh  * sin(flon(idat) * pio180)

        ! The following version uses offsets to shift from the 'best track'
        ! initial location to the effective GFS initial location

        !n zeh = erad * sin((flat(idat)-.095) * pio180)
        !n reh = erad * cos((flat(idat)-.095) * pio180)  ! distance from earth center
        !n xeh = reh  * cos((flon(idat)-.212) * pio180)
        !n yeh = reh  * sin((flon(idat)-.212) * pio180)

        ! Transform hurricane earth coords to whatever projection is in use

        call oplot_transform(iplt,xeh,yeh,zeh,flon(idat),flat(idat),xobs(idat),yobs(idat))

        ! Find "earth" coordinates of sim1 hurricane center

        zeh = erad * sin(latsim1(idat) * pio180)
        reh = erad * cos(latsim1(idat) * pio180)  ! distance from earth center
        xeh = reh  * cos(lonsim1(idat) * pio180)
        yeh = reh  * sin(lonsim1(idat) * pio180)

        ! Transform hurricane earth coords to whatever projection is in use

        call oplot_transform(iplt,xeh,yeh,zeh,lonsim1(idat),latsim1(idat),xsim1(idat),ysim1(idat))

        ! Find "earth" coordinates of sim2 hurricane center

        zeh = erad * sin(latsim2(idat) * pio180)
        reh = erad * cos(latsim2(idat) * pio180)  ! distance from earth center
        xeh = reh  * cos(lonsim2(idat) * pio180)
        yeh = reh  * sin(lonsim2(idat) * pio180)

        ! Transform hurricane earth coords to whatever projection is in use

        call oplot_transform(iplt,xeh,yeh,zeh,lonsim2(idat),latsim2(idat),xsim2(idat),ysim2(idat))

     enddo

  endif

  ! Set character font, line width, and plotted size

  call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
  call o_pcsetr('CL',2.)
  bsize = .016 * (op%hp2 - op%hp1) ! * 0.3

  call gslwsc(2.) ! line width

  ! Plot hurricane tracks from other models (now all in one group; iplt is 1 only)

  im1 = 1; im2 = 18

  do imodel = im1, im2
     jmodel = imodel

     ! Set plot line color

     call o_sflush()

     if (jmodel ==  1) icolor =  98
     if (jmodel ==  2) icolor = 101
     if (jmodel ==  3) icolor = 104
     if (jmodel ==  4) icolor = 106
     if (jmodel ==  5) icolor = 110
     if (jmodel ==  6) icolor = 113
     if (jmodel ==  7) icolor = 115
     if (jmodel ==  8) icolor = 117
     if (jmodel ==  9) icolor = 118
     if (jmodel == 10) icolor =   2
     if (jmodel == 11) icolor = 122
     if (jmodel == 12) icolor = 126
     if (jmodel == 13) icolor = 127
     if (jmodel == 14) icolor = 135
     if (jmodel == 15) icolor = 136
     if (jmodel == 16) icolor = 141
     if (jmodel == 17) icolor = 143
     if (jmodel == 18) icolor = 145

     call o_gsplci(icolor)
     call o_gstxci(icolor)
     call o_gsfaci(icolor)

     ! Plot color code

     xpt1 = 2.*235.e3

     ypt = 2.*215.e3 - 2.*18.e3 * real(imodel - im1)

     call o_plchhq(xpt1,ypt,model_id(jmodel,1),1.1*bsize,0.,-1.)

     do ihour = 1,numhours(jmodel)

        ! Plot only up to a selected time

        if (mhour(jmodel,ihour) > 72.1) exit

        if (ihour == 1) then
           call o_frstpt(xsims(jmodel,ihour),ysims(jmodel,ihour))
        else
           call o_vector(xsims(jmodel,ihour),ysims(jmodel,ihour))
        endif

print*, 'prt ',jmodel,icolor,mhour(jmodel,ihour),xsims(jmodel,ihour),ysims(jmodel,ihour)

     enddo

  enddo

  nfhour = nint(time8/21600.)

  ! Set color to black

  call o_gsplci(10)
  call o_gstxci(10)
  call o_gsfaci(10)
  call o_sflush()

  ! Plot observed hurricane eye location with black color

  do ifhour = 0,nfhour
     write(title,'(i2)') ifhour
     print*, 'plot_harvey_4 ',ifhour,xobs(ifhour),yobs(ifhour),bsize,trim(adjustl(title))
     call o_plchhq (xobs(ifhour),yobs(ifhour),trim(adjustl(title)),bsize,0.,0.)
  enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Return to not plot modeled location(s)

RETURN

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  ! Set color to red

  call o_gsplci(1)
  call o_gstxci(1)
  call o_gsfaci(1)
  call o_sflush()

  ! Plot initial sim1 hurricane eye location with red color

  print*, 'plot_harvey_5 ',xsim1(0),ysim1(0),bsize,trim(adjustl(title))
  write(title,'(a1)') '+'
  call o_plchhq (xsim1(0),ysim1(0),trim(adjustl(title)),bsize,0.,0.)

  ! Plot sim1 hurricane eye location with red color, excluding initial time

  do ifhour = 1,nfhour
     write(title,'(i2)') ifhour
     print*, 'plot_harvey_6 ',ifhour,xsim1(ifhour),ysim1(ifhour),bsize,trim(adjustl(title))
     call o_plchhq (xsim1(ifhour),ysim1(ifhour),trim(adjustl(title)),bsize,0.,0.)
  enddo

  ! Set color to green

  call o_gsplci(3)
  call o_gstxci(3)
  call o_gsfaci(3)
  call o_sflush()

  ! Plot initial sim2 hurricane eye location with red color

  print*, 'plot_harvey_7 ',xsim2(0),ysim2(0),bsize,trim(adjustl(title))
  write(title,'(a1)') '+'
  call o_plchhq (xsim2(0),ysim2(0),trim(adjustl(title)),bsize,0.,0.)

  ! Plot sim1 hurricane eye location with red color, excluding initial time

  do ifhour = 1,nfhour
     write(title,'(i2)') ifhour
     print*, 'plot_harvey_8 ',ifhour,xsim2(ifhour),ysim2(ifhour),bsize,trim(adjustl(title))
     call o_plchhq (xsim2(ifhour),ysim2(ifhour),trim(adjustl(title)),bsize,0.,0.)
  enddo

end subroutine plot_harvey

