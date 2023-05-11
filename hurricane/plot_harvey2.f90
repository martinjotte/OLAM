subroutine plot_harvey(iplt)

  use mem_basic,   only: press
  use mem_grid,    only: mwa, lpw, glatw, glonw
  use misc_coms,   only: time8, mstp
  use oplot_coms,  only: op
  use consts_coms, only: pio180, erad
  use hcane_rz,    only: hlata, hlona, htima

  implicit none

  integer, intent(in) :: iplt

  logical, save :: newcall = .true.
  real :: zeh, reh, xeh, yeh, press_lowest, bsize, rhour, xs, ys
  character(len=2) :: title

  integer :: iline, imodel, jmodel, ihour, im1, im2, icolor, kstp, khour, lhour
  character(2) :: basin
  character(4) :: pastid, modelid
  integer :: cyclone_num, yymmddhh, technum, tau, idlat, idlon, vmax, mslp

  integer, parameter :: nmodels = 40, nhours = 50

  integer, save :: mhour(nmodels,nhours)
  integer, save :: numhours(nmodels)
  real, save ::  mlat(nmodels,nhours),  mlon(nmodels,nhours)
  real,save :: xsims(nmodels,nhours), ysims(nmodels,nhours)
  character(4), save :: model_id(nmodels,nhours)

  real :: xpt1, xpt2, xpt3, ypt

 !   RETURN  ! Do this return for main 2-day run or normal plotonly run,
            ! and comment out this return when doing spaghetti plots

  ! if (iplt > 1) return

  if (newcall) then
     newcall = .false.

     ! Read locations from multiple forecast models at OLAM initialization time

     imodel = 0
     pastid = '0000'

     open(32,file='aug24_12z.dat2_olamA',status='old',form='formatted')
     do

        read(32,'(a2,2x,i2,2x,i10,2x,i2,2x,a4,2x,i3,2x,i3,3x,i4,3x,i3,2x,i4)') &
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
        mlat    (imodel,ihour) =  0.1 * real(idlat)
        mlon    (imodel,ihour) = -0.1 * real(idlon)      + 0.01 * real(imodel - 17) ! offset for visual separation
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

  endif

  ! Set character font, line width, and plotted size

  call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
  call o_pcsetr('CL',2.)
  bsize = .016 * (op%hp2 - op%hp1) ! * 0.3

  call gslwsc(2.) ! line width

  ! Plot hurricane tracks from other models (now all in one group; iplt is 1 only)

  im1 = 16; im2 = 19

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
     if (jmodel == 16) icolor =  10
     if (jmodel == 17) icolor = 143
     if (jmodel == 18) icolor = 144
     if (jmodel == 19) icolor = 145

     call o_gsplci(icolor)
     call o_gstxci(icolor)
     call o_gsfaci(icolor)

     ! Plot color code

     xpt1 = 235.e3
     ypt  = 215.e3 - 18.e3 * real(imodel - im1)

     ! call o_plchhq(xpt1,ypt,model_id(jmodel,1),1.1*bsize,0.,-1.)

     ! Plot results

     do ihour = 1,numhours(jmodel)

        if (mhour(jmodel,ihour) > 72.1) exit    ! Plot only up to a selected time

        if (ihour == 1) then
           call o_frstpt(xsims(jmodel,ihour),ysims(jmodel,ihour))
        else
           call o_vector(xsims(jmodel,ihour),ysims(jmodel,ihour))
        endif

        write(title,'(i2)') mhour(jmodel,ihour)
        if (ihour > 1 .or. jmodel == 16) then   ! BEST TRACK is model 16
           call o_plchhq (xsims(jmodel,ihour),ysims(jmodel,ihour),trim(adjustl(title)),bsize,0.,0.)
        endif

     enddo

  enddo

  ! Process time series of hurricane locations for present simulation

  call o_sflush()

  icolor = 127

  call o_gsplci(icolor)
  call o_gstxci(icolor)
  call o_gsfaci(icolor)

  lhour = 0

  do kstp = 0, mstp

     ! Find "earth" coordinates of sims hurricane center

     zeh = erad * sin(hlata(kstp) * pio180)
     reh = erad * cos(hlata(kstp) * pio180)  ! distance from earth center
     xeh = reh  * cos(hlona(kstp) * pio180)
     yeh = reh  * sin(hlona(kstp) * pio180)

     ! Transform hurricane earth coords to whatever projection is in use

     call oplot_transform(iplt,xeh,yeh,zeh,xs,ys)

     if (kstp == 0) then
        call o_frstpt(xs,ys)
     else
        call o_vector(xs,ys)
     endif

     rhour = htima(kstp) / 3600.
     khour = int(rhour)

     if (khour > 0 .and. khour > lhour) then
        write(title,'(i2)') khour
        call o_plchhq (xs,ys,trim(adjustl(title)),bsize,0.,0.)
        lhour = khour
     endif

  enddo

end subroutine plot_harvey

