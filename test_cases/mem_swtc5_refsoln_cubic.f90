Module mem_swtc5_refsoln_cubic

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

real, allocatable :: zanal00_swtc5(:,:),zanal15_swtc5(:,:)

integer :: nlon, nlat

real(dp) :: zglon
real(dp) :: pi

real(dp), allocatable :: zglat(:)
real(dp), allocatable :: wt(:,:,:)

Contains

   subroutine fill_swtc5()

   use misc_coms, only: io6

   implicit none

   integer :: ifile,ilat,ilon

   character(80) :: fname_swtc5
   logical :: fexists

   pi = 2._dp*asin(1._dp)

   nlon = 1280
   nlat = 640

   allocate(zanal00_swtc5(nlon,nlat),zanal15_swtc5(nlon,nlat))

   allocate (zglat(-1:nlat+2),wt(4, 2, 0:nlat))

   do ifile = 1,2

      if (ifile == 1) then
         fname_swtc5 = '../../swtc5_refsoln/hmodc5_0000_E000570.dat'
      elseif (ifile == 2) then
         fname_swtc5 = '../../swtc5_refsoln/hmodc5_0360_E000570.dat'
      endif

      inquire(file=fname_swtc5, exist=fexists)
      if (.not. fexists) then
         write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(io6,*) '!!!   Trying to open file:'
         write(io6,*) '!!!   '//trim(fname_swtc5)
         write(io6,*) '!!!   but it does not exist.  The run is ended.    '
         write(io6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         stop 'in fill_swtc5'
      endif

      open(29,file=trim(fname_swtc5),form='FORMATTED',status='OLD',action='READ')

! Read data to fill interior of arrays

      do ilat = 1, nlat
         do ilon = 1, nlon
            if (ifile == 1) then
               read(29,'(3e30.16e3)') &
                  zglon,zglat(nlat-ilat+1),zanal00_swtc5(ilon,nlat-ilat+1)
            elseif (ifile == 2) then
               read(29,'(3e30.16e3)') &
                  zglon,zglat(nlat-ilat+1),zanal15_swtc5(ilon,nlat-ilat+1)
            endif

            zglat(nlat-ilat+1) = zglat(nlat-ilat+1) * pi / 180._dp
         enddo
      enddo

      close (29)
   enddo

   zglat(-1)     = -pi - zglat(2)
   zglat( 0)     = -pi - zglat(1)
   zglat(nlat+1) =  pi - zglat(nlat)
   zglat(nlat+2) =  pi - zglat(nlat-1)

  ! Compute the weights for the horizontal interpolation from the Gaussian grid
  ! of the IFS to the triangular grid; weights for bicubic interpolation

   do ilat = 0,nlat
      call npr_lcbas(zglat(ilat-1), wt(1,1,ilat), wt(1,2,ilat))
   enddo

   return
   end subroutine fill_swtc5

!===============================================================================

   subroutine npr_lcbas(grd,bas1,bas2)

   implicit none

!------------------------------arguments--------------------------------

! input arguments

   real(kind=dp) :: grd(4) ! grid stencil

! output arguments

   real(kind=dp) :: bas1(4) ! grid values on stencil
   real(kind=dp) :: bas2(4) ! lagrangian basis functions

! grid value differences used in weights

   real(kind=dp) :: x0mx1, x0mx2, x0mx3, x1mx2, x1mx3, x2mx3

   x0mx1 = grd(1) - grd(2)
   x0mx2 = grd(1) - grd(3)
   x0mx3 = grd(1) - grd(4)
   x1mx2 = grd(2) - grd(3)
   x1mx3 = grd(2) - grd(4)
   x2mx3 = grd(3) - grd(4)

   bas1(1) = grd(1)
   bas1(2) = grd(2)
   bas1(3) = grd(3)
   bas1(4) = grd(4)

   bas2(1) =  1._dp / (x0mx1 * x0mx2 * x0mx3)
   bas2(2) = -1._dp / (x0mx1 * x1mx2 * x1mx3)
   bas2(3) =  1._dp / (x0mx2 * x1mx2 * x2mx3)
   bas2(4) = -1._dp / (x0mx3 * x1mx3 * x2mx3)

   end subroutine npr_lcbas

!===============================================================================

   subroutine npr_bicubics(sg,glatw,glonw,s)

   implicit none

   integer indexs(2)
   real :: sg(nlon,nlat), s

   real(kind=dp) :: dlon, dloni, sixth, s1, s2, s3, s4
   real(kind=dp) :: xlat, xlon, x1, x2, x3, x4, y1, y2, y3, y4
   real(kind=dp) :: pi
   integer ilat, ilon, ilonm1
   integer ilonp1, ilonp2, jlon, jlonm1, jlonp1, jlonp2, mlat

   real :: glatw,glonw

   pi    = 2._dp*asin(1._dp)
   dlon  = 2._dp*pi/nlon
   dloni = 1._dp/dlon
   sixth = 1._dp/6._dp
   mlat  = nlat + nlat - 1

   xlat = glatw*pi/180._dp
   xlon = glonw*pi/180._dp
   if (xlon < 0) xlon = xlon + 2._dp*pi

   call npr_getindex(xlat, xlon, indexs)
   ilat = indexs(1)
   ilon = indexs(2)

   ilonp1 = mod(ilon, nlon) + 1
   ilonp2 = mod(ilonp1, nlon) + 1
   ilonm1 = mod(ilon+nlon-2, nlon) + 1
   jlon   = mod(ilon+nlon/2-1, nlon) + 1
   jlonp1 = mod(jlon, nlon) + 1
   jlonp2 = mod(jlonp1, nlon) + 1
   jlonm1 = mod(jlon+nlon-2, nlon) + 1

! Compute the x-interpolants for the four latitudes.

   x2 = (xlon - dlon*(ilon - 1))*dloni
   x1 = -x2 - 1._dp
   x3 = -x2 + 1._dp
   x4 =  x2 - 2._dp

   if (ilat <= 1) then
      s1 = sg(jlonp1,2-ilat)*x2 &
         + sg(jlon  ,2-ilat)*x3
   else
      s1 = sg(ilonp1,ilat-1)*x2 &
         + sg(ilon  ,ilat-1)*x3
   endif

   if (ilat == 0) then
      s2 = sg(jlonm1,1)*x2*x3*x4*sixth  &
         + sg(jlon  ,1)*x1*x3*x4*0.5_dp &
         + sg(jlonp1,1)*x1*x2*x4*0.5_dp &
         + sg(jlonp2,1)*x1*x2*x3*sixth
   else
      s2 = sg(ilonm1,ilat)*x2*x3*x4*sixth  &
         + sg(ilon  ,ilat)*x1*x3*x4*0.5_dp &
         + sg(ilonp1,ilat)*x1*x2*x4*0.5_dp &
         + sg(ilonp2,ilat)*x1*x2*x3*sixth
   endif

   if (ilat == nlat) then
      s3 = sg(jlonm1,nlat)*x2*x3*x4*sixth  &
         + sg(jlon  ,nlat)*x1*x3*x4*0.5_dp &
         + sg(jlonp1,nlat)*x1*x2*x4*0.5_dp &
         + sg(jlonp2,nlat)*x1*x2*x3*sixth
   else
      s3 = sg(ilonm1,ilat+1)*x2*x3*x4*sixth  &
         + sg(ilon  ,ilat+1)*x1*x3*x4*0.5_dp &
         + sg(ilonp1,ilat+1)*x1*x2*x4*0.5_dp &
         + sg(ilonp2,ilat+1)*x1*x2*x3*sixth
   endif

   if (ilat > nlat-1) then
      s4 = sg(jlonp1,mlat-ilat)*x2 &
         + sg(jlon  ,mlat-ilat)*x3
   else
      s4 = sg(ilonp1,ilat+2)*x2 &
         + sg(ilon  ,ilat+2)*x3
   endif

! Compute the y-interpolant.

   y1 = xlat - wt(1,1,ilat)
   y2 = xlat - wt(2,1,ilat)
   y3 = xlat - wt(3,1,ilat)
   y4 = xlat - wt(4,1,ilat)

   s = s1*y2*y3*y4*wt(1,2,ilat) &
     + s2*y1*y3*y4*wt(2,2,ilat) &
     + s3*y1*y2*y4*wt(3,2,ilat) &
     + s4*y1*y2*y3*wt(4,2,ilat)

!write(6,'(a,2i6,15e12.3)') 'bc1 ',ilat,ilon,glatw,glonw,xlat,xlon,s, &
!    y1,y2,y3,y4,wt(1,2,ilat),wt(2,2,ilat),wt(3,2,ilat),wt(4,2,ilat)

   return
   end subroutine npr_bicubics

!==============================================================================

   subroutine npr_getindex(gilat, gilon, indexs)

   implicit none

   real(kind=dp) :: gilon, gilat
   integer :: indexs(2)
   real(kind=dp) :: pi, dlati,  dlon, dloni,  xlat,  xlon
   integer :: ilat ,ilon

   pi    = 2._dp*asin(1._dp)
   dlati = 1._dp/(zglat(nlat/2+1) - zglat(nlat/2))
   dlon  = 2._dp*pi/nlon
   dloni = 1._dp/dlon

   xlat = gilat
   ilat = (xlat - zglat(1))*dlati + 1._dp
   if (xlat >= zglat(ilat+1)) ilat = ilat + 1

   indexs(1) = ilat

   xlon = gilon
   if (xlon < 0) xlon = xlon + 2._dp*pi
   ilon = xlon*dloni + 1._dp
   if (ilon == nlon+1) ilon = nlon

   indexs(2) = ilon

   return
   end subroutine npr_getindex

End Module mem_swtc5_refsoln_cubic


