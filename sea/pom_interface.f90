subroutine pom_startup()

  use pom2k1d,     only: nzpom, dtl2, cbc, cbcmin, cbcmax, y, yy, &
                         vonk, z0b, pom, alloc_pom, filltab_pom

  use consts_coms, only: omega2, pio180

  use mem_sfcg,    only: sfcg

  use mem_sea,     only: msea, omsea

  use leaf_coms,   only: dt_leaf

  implicit none

  integer :: isea, iwsfc, k

  ! Later, get these from OLAM simulation parameters

  dtl2 = dt_leaf * 2.0

  call alloc_pom(msea)
  call filltab_pom()

  do isea = 2,msea
     iwsfc = isea + omsea

     if (sfcg%pom_active(iwsfc)) then
        pom%cor(isea) = omega2 * sin(sfcg%glatw(iwsfc) * pio180)
     endif
  enddo

  do k = 1,nzpom-1
     cbc(k) = (vonk / log((yy(k) - y(k+1)) / z0b))**2
     cbc(k) = max(cbcmin, cbc(k))

     ! If the following is invoked, then it is probable that the wrong
     ! choice of z0b or vertical spacing has been made:

     cbc(k) = min(cbcmax, cbc(k))
  enddo
  cbc(nzpom) = cbc(nzpom-1)

end subroutine pom_startup

!===============================================================================

subroutine pom_init()

  use mem_sfcg, only: sfcg
  use mem_sea,  only: sea, msea, omsea
  use pom2k1d,  only: nzpom, pom, y, yy, small

  implicit none

  integer :: k, isea, iwsfc
  real    :: tlen  ! turbulent length scale [m]
  real    :: seatempc

  do isea = 2,msea
     iwsfc = isea + omsea

     ! Fill surface temperature of all POM cells, even if outside pom_active regions

     seatempc = sea%seatc(isea) - 273.15

     pom%potmp(1,isea) = seatempc

     if (sfcg%pom_active(iwsfc)) then

        tlen = 0.1 * y(pom%kba(isea)) 

        pom%wubot(isea) = 0.
        pom%wvbot(isea) = 0.

        ! Default initialization: Use SST at surface and vary linearly with depth
        ! to 1 deg C at 3000 m depth; use constant 1 deg C below 3000 m.

        do k = 1,nzpom
           pom%potmp(k,isea) = seatempc + (1.0 - seatempc) &
                             * min(1.0, (yy(k) - yy(1)) / -3000.)
           pom%salin(k,isea) = 35.
           pom%q2   (k,isea) = small
           pom%q2l  (k,isea) = pom%q2(k,isea) * tlen
           pom%u    (k,isea) = 0.
           pom%v    (k,isea) = 0.

           pom%potmpb(k,isea) = pom%potmp(k,isea)
           pom%salinb(k,isea) = pom%salin(k,isea)
           pom%q2b   (k,isea) = pom%q2   (k,isea)
           pom%q2lb  (k,isea) = pom%q2l  (k,isea)
           pom%ub    (k,isea) = pom%u    (k,isea)
           pom%vb    (k,isea) = pom%v    (k,isea)

           pom%kh    (k,isea) = tlen * sqrt(pom%q2b(k,isea))
           pom%km    (k,isea) = pom%kh(k,isea)
           pom%kq    (k,isea) = pom%kh(k,isea)
        enddo

     endif
  enddo
 
end subroutine pom_init

!===============================================================================

subroutine plot_pom()

  use mem_sea, only: omsea
  use pom2k1d, only: nzpom, pom, yy

  implicit none

  real, parameter :: aspect = .7
  real, parameter :: scalelab = .014
  real, save :: valmin, valmax, valinc

  real, allocatable :: vctr18(:), val(:,:)

  integer :: labincx, labincy, iwsfc, isea, iv, k, jpom, kb

  integer, parameter :: icolor(6) = [16, 109, 12, 11, 9, 8]
! orange, blue green, dark red, purple, dark green, dark blue

  real :: xmin, xmax, xinc
  real :: ymin, ymax, yinc

  integer, parameter :: ijpom(4) = [162230, 167916, 159390, 164544] ! selected POM1D columns to plot
  character*1, parameter :: ip(4) = ['1','2','3','4']

  labincx = 5
  labincy = 5

  xmin = -11.
  xmax = 41.
  xinc = 2.
  
  ymin =  1.05 * yy(nzpom)
  ymax = -0.05 * yy(nzpom)
  yinc = 100.

  call o_reopnwk()

  call plotback()

  do jpom = 1,4
     iwsfc = ijpom(jpom)
     isea = iwsfc - omsea

     kb = pom%kba(isea)

     allocate(vctr18(kb-1))
     allocate(val(kb-1,6))

     do k = 1,kb-1
        vctr18(k) = yy(k)

        val(k,1) = pom%potmp(k,isea)
        val(k,2) = pom%salin(k,isea)
        val(k,3) = pom%q2(k,isea) * 100.
        val(k,4) = pom%q2l(k,isea) * 100.
        val(k,5) = pom%u(k,isea) * 10.
        val(k,6) = pom%v(k,isea) * 10.
     enddo

     do iv = 1,6

        call oplot_xy2l(ip(jpom),'N','a','N',aspect,scalelab,icolor(iv),0, &
                       kb-1,  val(:,iv), &
                       vctr18, &
                       'vals','Z', &
                       xmin, xmax, xinc, labincx, &
                       ymin, ymax, yinc, labincy  )

     enddo

     deallocate(vctr18,val)

  enddo

  call o_frame()

  call o_clswk()

end subroutine plot_pom

!===============================================================================

subroutine oplot_xy2l(panel,frameoff,pltborder,colorbar0,aspect,scalelab,linecolor,ndashes, &
   n,xval,yval,xlab,ylab, &
   xmin,xmax,xinc,labincx  ,ymin,ymax,yinc,labincy)
 
use oplot_coms, only: op
use misc_coms,  only: io6

! This routine is a substitute for NCAR Graphics routine ezxy to allow 
! control over fonts, labels, axis labels, line width, scaling, etc.
! Pass in a value of 1 for n to not plot (only draw frame and ticks)

implicit none

character(len=1), intent(in) :: panel,frameoff,pltborder,colorbar0
integer, intent(in) :: n,labincx,labincy,linecolor,ndashes
real, intent(in) :: aspect,scalelab,xmin,xmax,xinc,ymin,ymax,yinc
real, intent(in) :: xval(n),yval(n)
character(len=*), intent(in) :: xlab,ylab

integer :: i,logy,itickvalq
real :: xl,yl,dx,dy,tickval,xmargin,ymargin,sizelab,xlabx,ylaby,tickvalq
character(len=20)  :: numbr,numbr2

integer :: icyc
real :: dashlen, dist, remain, step, xfac, yfac, asp2, x, y, eps

! Set plot color (black)

 call o_sflush()
 call o_gsplci(10)
 call o_gsfaci(10)
 call o_gstxci(10)
 call o_gslwsc(1.) ! line width

! Scale local working window (0,1,0,1) 
! to plotter coordinates (op%hp1,op%hp2,op%vp1,op%vp2)

 call oplot_panel(panel,frameoff,pltborder,colorbar0,aspect,'N')
 call o_set(op%hp1,op%hp2,op%vp1,op%vp2,0.,1.,0.,1.,1)

! Draw frame

 call o_frstpt(op%fx1,op%fy1)
 call o_vector(op%fx2,op%fy1)
 call o_vector(op%fx2,op%fy2)
 call o_vector(op%fx1,op%fy2)
 call o_vector(op%fx1,op%fy1)

! Specify font # and scale font size to designated plotter coordinates

 call o_sflush()
 call o_pcseti ('FN',4)  ! set font number to 4 (font 2 is similar but wider spacing)
 call o_pcsetr('CL',1.)  ! set character line width to 1

sizelab = scalelab * (op%hp2 - op%hp1)

! Write x axis label

if (pltborder == 'a' .or. &   
    panel     == 'N' .or. &
    panel     == '1' .or. &
    panel     == '2' .or. &
    panel     == '9') then

   xlabx = .5 * (op%fx1 + op%fx2)
   call o_plchhq(xlabx,op%xlaby,trim(xlab),sizelab, 0.,0.)
endif

! Write y axis label

if (pltborder == 'a' .or. &   
     panel    == 'N' .or. &
     panel    == '1' .or. &
     panel    == '3' .or. &
     panel    == '5') then

   ylaby = .5 * (op%fy1 + op%fy2)
   call o_plchhq(op%ylabx,ylaby,trim(ylab),sizelab,90.,0.)
endif

! Scale local working window (xmin,xmax,0.,1.) 
! to plotter coordinates (op%h1,op%h2,op%vp1,op%vp2)

 call o_set(op%h1,op%h2,op%vp1,op%vp2,xmin,xmax,0.,1.,1)

! Plot and label X-axis ticks

tickval = nint(xmin/xinc) * xinc
if (tickval < xmin - .001 * xinc) tickval = tickval + xinc
   
do while (tickval < xmax + .001 * xinc)
   
   if (mod(nint(tickval/xinc),labincx) == 0) then  ! Only for long ticks

      dy = .014

      ! Encode and plot current X tick label

      if (pltborder == 'a' .or. &
          panel     == 'N' .or. &
          panel     == '1' .or. &
          panel     == '2' .or. &
          panel     == '9') then

         if (xinc * labincx >= .999) then
            write (numbr,'(i6)') nint(tickval)
         elseif (xinc * labincx >= .0999) then
            write (numbr,'(f5.1)') tickval
         elseif (xinc * labincx >= .00999) then
            write (numbr,'(f5.2)') tickval
         elseif (xinc * labincx >= .000999) then
            write (numbr,'(f6.3)') tickval
         else
            write (numbr,'(f7.4)') tickval
         endif
   
         call o_plchhq(tickval,op%xtlaby,trim(adjustl(numbr)),sizelab,0.,0.)
      endif

   else                                          ! Only for short ticks
      dy = .007
   endif
   
! Plot current X tick
   
   call o_frstpt(tickval,op%fy1)
   call o_vector(tickval,op%fy1 + dy)
   call o_frstpt(tickval,op%fy2)
   call o_vector(tickval,op%fy2 - dy)
   
   tickval = tickval + xinc
enddo

! Scale local working window (0.,1.,ymin,ymax) 
! to plotter coordinates (op%hp1,op%hp2,op%v1,op%v2)

 call o_set(op%hp1,op%hp2,op%v1,op%v2,0.,1.,ymin,ymax,1)

! Plot and label Y-axis ticks

tickval = nint(ymin/yinc) * yinc

if ((tickval - ymax) / (ymin - ymax) > 1.001) tickval = tickval + yinc
   
do while ((tickval - ymin) / (ymax - ymin) < 1.001)
   
   if (mod(nint(tickval/yinc),labincy) == 0) then  ! Only for long ticks

      dx = .014
 
      ! Encode and plot current Y tick label

      if (pltborder == 'a' .or. &
          panel     == 'N' .or. &
          panel     == '1' .or. &
          panel     == '3' .or. &
          panel     == '5') then

         ! Encode current Y tick label

         if (abs(yinc * labincy) >= .999) then
            write (numbr,'(i6)') nint(tickval)
         elseif (yinc * labincy >= .0999) then
            write (numbr,'(f5.1)') tickval
         elseif (yinc * labincy >= .00999) then
            write (numbr,'(f5.2)') tickval
         elseif (yinc * labincy >= .000999) then
            write (numbr,'(f6.3)') tickval
         else
       
            logy = int(log10(yinc * labincy)) - 2
            tickvalq = tickval * 10. ** (-logy)         
            itickvalq = nint(tickvalq)

            ! If significand is at or above 10, reduce

            do while (abs(itickvalq) >= 10)
               logy = logy + 1
               tickvalq = tickvalq * .1
               itickvalq = nint(tickvalq)
            enddo

            ! Use only one of the following 2 lines

            ! write (numbr,'(f4.1)') tickvalq   ! If real significand is required
            write (numbr,'(i2)') itickvalq    ! If integer significand is ok

            write (numbr2,'(i3)') logy
            numbr = trim(adjustl(numbr))
            numbr2 = trim(adjustl(numbr2))

            ! Determine whether significand, power of 10, or both are to be plotted

            if (itickvalq == 1) then 
               numbr = '10:S3:'//trim(numbr2)//'        '
            elseif (itickvalq == -1) then 
               numbr = '-10:S3:'//trim(numbr2)//'        '
            elseif (itickvalq /= 0) then
               numbr = trim(adjustl(numbr))//'x'//'10:S3:'//trim(adjustl(numbr2))//'        '
            endif

         endif

         ! Plot Y tick label

         ! call o_plchhq(op%ytlabx,tickval,numbr(1:len_trim(numbr)),sizelab,0.,1.)
         call o_plchhq(op%ytlabx,tickval,trim(adjustl(numbr)),sizelab,0.,1.)

      endif ! frameoff/panel

   else                                          ! Only for short ticks
      dx = .007
   endif

! Plot current Y tick

   call o_frstpt(op%fx1     ,tickval)
   call o_vector(op%fx1 + dx,tickval)
   call o_frstpt(op%fx2     ,tickval)
   call o_vector(op%fx2 - dx,tickval)
   
   tickval = tickval + yinc
enddo

! Scale local working window (xmin,xmax,ymin,ymax)
!  to plotter coordinates (op%h1,op%h2,op%v1,op%v2)

 call o_set(op%h1,op%h2,op%v1,op%v2,xmin,xmax,ymin,ymax,1)

! Plot values

! Set plot color (linecolor)

 call o_sflush()
 call o_gsplci(linecolor)
 call o_gsfaci(linecolor)
 call o_gstxci(linecolor)
 call o_gslwsc(1.) ! line width

 if (panel == '3')then
    if     (linecolor == 16) then
       call o_plchhq(xmin,ymax+(ymax-ymin)*0.70,'Temperature (deg C)',1.5*sizelab,0.,-1.)
    elseif (linecolor == 109) then
       call o_plchhq(xmin,ymax+(ymax-ymin)*0.60,'Salinity (g/kg) ',1.5*sizelab,0.,-1.)
    elseif (linecolor == 12) then
       call o_plchhq(xmin,ymax+(ymax-ymin)*0.50,'TKE*2 (m^2/s^2 x 100) ',1.5*sizelab,0.,-1.)
    elseif (linecolor == 11) then
       call o_plchhq(xmin,ymax+(ymax-ymin)*0.40,'TKE*L*2 (m^3/s^2 x 100) ',1.5*sizelab,0.,-1.)
    elseif (linecolor == 9) then
       call o_plchhq(xmin,ymax+(ymax-ymin)*0.30,'U (m/s x 10) ',1.5*sizelab,0.,-1.)
    elseif (linecolor == 8) then
       call o_plchhq(xmin,ymax+(ymax-ymin)*0.20,'V (m/s x 10) ',1.5*sizelab,0.,-1.)
    endif
 endif

if (ndashes <= 0) then

   call o_frstpt(xval(1),yval(1))
   do i = 2,n
      call o_vector(xval(i),yval(i))
   enddo

else

   eps = 1.e-6 * (xmax - xmin)

   xfac = (op%h2 - op%h1) / (xmax - xmin)
   yfac = (op%v2 - op%v1) / (ymax - ymin)

   asp2 = (yfac / xfac)**2

   dashlen = (xmax - xmin) / real(2 * ndashes)
   remain = dashlen

   x = xval(1)
   y = yval(1)

   call o_frstpt(x,y)

   icyc = 1
   i = 2

   do while (i < n)

      dist = sqrt((xval(i) - x)**2 + asp2 * (yval(i) - y)**2) 

      if (remain > dist + eps) then

         step = dist

         x = xval(i)
         y = yval(i)

         if (icyc > 0) then
            call o_vector(x,y)
         else
            call o_frstpt(x,y)
         endif

         remain = remain - step
         i = i + 1

      else

         step = remain

         x = x + (xval(i) - x) * step / dist
         y = y + (yval(i) - y) * step / dist

         if (icyc > 0) then
            call o_vector(x,y)
         else
            call o_frstpt(x,y)
         endif

         remain = dashlen
         icyc = -icyc

      endif

   enddo

endif

end subroutine oplot_xy2l

