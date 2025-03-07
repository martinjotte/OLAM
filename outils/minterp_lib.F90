module minterp_lib

  use ll_bins, only: binset_vars
  implicit none

  type minterp_vars   ! data structure to store values at each M point
                      ! for quickly computing weighting coefficients
     real :: w1factx, w1facty
     real :: w2factx, w2facty
     real :: xg3
     real :: yg3

     real :: xemin, xemax
     real :: yemin, yemax
     real :: zemin, zemax

     real :: cosmlat, sinmlat
     real :: cosmlon, sinmlon
  end type minterp_vars

  type(minterp_vars), allocatable :: minterp(:)

  Type minterp_wghts  ! data structure to store weighting coefs
     integer :: imglobe ! M "triangle" that contains point (global index)
     real    :: wt(3)   ! 3 weighting coefs for each itab_m(im)%w(1:3)
  End Type minterp_wghts

  type(binset_vars) :: bsetm

  interface get_weights_lonlat
     module procedure               &
          get_weights_lonlat_point, &
          get_weights_lonlat_1d,    &
          get_weights_lonlat_2d,    &
          get_weights_lonlat_2db,   &
          get_weights_lonlat_grid,  &
          get_weights_lonlat_gridb
  end interface get_weights_lonlat

  interface interp_level
     module procedure     &
          interp_level_1, &
          interp_level_2
  end interface interp_level

  interface interp_column
     module procedure      &
          interp_column_1, &
          interp_column_2
  end interface interp_column

  interface interp_column_tolev
     module procedure            &
          interp_column_tolev_1, &
          interp_column_tolev_2
  end interface interp_column_tolev

  integer, parameter :: ihuge = huge(1)

  private :: binset_vars

contains

!=========================================================================

subroutine init_minterp()

  use mem_ijtabs,  only: itab_m, mrls
  use mem_grid,    only: mma, xew, yew, zew, xem, yem, zem, glonw, glatw
  use consts_coms, only: eradi, r8
  use mem_para,    only: myrank
  use ll_bins,     only: itab_grid_vars, latlon_bins
  use misc_coms,   only: nxp

  implicit none

  integer :: im, iw, n
  real    :: dxe(3), dye(3), dze(3), xw(3), yw(3), denom, denomi
  real    :: raxis, raxisi, res

  type(itab_grid_vars), allocatable :: itab_m0(:)

  allocate( minterp( mma ) )
  allocate( itab_m0( mma ) )

  !$omp parallel do private(raxis, raxisi, denom, denomi, &
  !$omp                     dxe, dye, dze, xw, yw, n, iw)
  do im = 2, mma

     if ( itab_m(im)%irank /= myrank .or. itab_m(im)%npoly /= 3 ) then
        itab_m0(im)%np = 0  ! flag to indicate to skip this point
        cycle
     endif

     itab_m0(im)%np = 3

     do n = 1, 3
        iw = itab_m(im)%iw(n)
        itab_m0(im)%glats(n) = glatw(iw)
        itab_m0(im)%glons(n) = glonw(iw)
     enddo

     raxis = sqrt( xem(im)**2 + yem(im)**2 )
     minterp(im)%sinmlat = zem(im) * eradi
     minterp(im)%cosmlat = raxis   * eradi

     ! For points less than 100 m from Earth's polar axis, make arbitrary
     ! assumption that longitude = 0 deg.  This is just to settle on a PS
     ! planar coordinate system in which to do the algebra.

     if (raxis >= 1.e2) then
        raxisi = 1.0 / raxis
        minterp(im)%sinmlon = yem(im) * raxisi
        minterp(im)%cosmlon = xem(im) * raxisi
     else
        minterp(im)%sinmlon = 0.
        minterp(im)%cosmlon = 1.
     endif

     dxe = xew( itab_m(im)%iw(1:3) ) - xem(im)
     dye = yew( itab_m(im)%iw(1:3) ) - yem(im)
     dze = zew( itab_m(im)%iw(1:3) ) - zem(im)

     minterp(im)%xemin = 1.01 * minval( dxe ) + xem(im)
     minterp(im)%xemax = 1.01 * maxval( dxe ) + xem(im)

     minterp(im)%yemin = 1.01 * minval( dye ) + yem(im)
     minterp(im)%yemax = 1.01 * maxval( dye ) + yem(im)

     minterp(im)%zemin = 1.01 * minval( dze ) + zem(im)
     minterp(im)%zemax = 1.01 * maxval( dze ) + zem(im)

     call de_gn_mult( 3, dxe, dye, dze, &
                      minterp(im)%cosmlat, minterp(im)%sinmlat, &
                      minterp(im)%cosmlon, minterp(im)%sinmlon, xw, yw )

     denom = (xw(1)-xw(3)) * (yw(2) - yw(3)) &
           - (xw(2)-xw(3)) * (yw(1) - yw(3))

     denomi = 1.0 / denom

     minterp(im)%w1factx = (yw(2) - yw(3)) * denomi
     minterp(im)%w1facty = (xw(3) - xw(2)) * denomi

     minterp(im)%w2factx = (yw(3) - yw(1)) * denomi
     minterp(im)%w2facty = (xw(1) - xw(3)) * denomi

     minterp(im)%xg3 = xw(3)
     minterp(im)%yg3 = yw(3)

  enddo
  !$omp end parallel do

  res = 7150.e3 / real( nxp * 2**(mrls-1) ) * 7.

  call latlon_bins(mma, res, itab_m0, bsetm)

end subroutine init_minterp

!=========================================================================

subroutine get_weights_lonlat_point(glon, glat, wgts, xe, ye, ze)

  use  ll_bins,    only: gridcells_from_latlon_bins
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use  misc_coms,  only: iparallel
  use  mpi
#endif

  implicit none

  real,                intent(in)  :: glon
  real,                intent(in)  :: glat
  type(minterp_wghts), intent(out) :: wgts
  real,      optional, intent(in)  :: xe, ye, ze

  integer, pointer                 :: nm
  integer, pointer, contiguous     :: impts(:)
  integer                          :: j
  real                             :: rads, raxis, xeg, yeg, zeg
  real                             :: cosglat, singlat
  real                             :: cosglon, singlon

#ifdef OLAM_MPI
  integer                          :: ig, ierr
#endif

  wgts%imglobe = 1

  call gridcells_from_latlon_bins(glat, glon, bsetm, nm, impts)

  if (present(xe) .and. present(ye) .and. present(ze)) then

     xeg = xe
     yeg = ye
     zeg = ze

  else

     rads    = glat * pio180
     cosglat = cos(rads)
     singlat = sin(rads)

     rads    = glon * pio180
     cosglon = cos(rads)
     singlon = sin(rads)

     raxis = erad  * cosglat
     zeg   = erad  * singlat
     xeg   = raxis * cosglon
     yeg   = raxis * singlon

  endif

  do j = 1, nm
     call mweights_xyze(impts(j), xeg, yeg, zeg, wgts)
     if (wgts%imglobe > 1) exit
  enddo

#ifdef OLAM_MPI
  if (iparallel == 1) then

     ig = ihuge
     if (wgts%imglobe > 1) ig = wgts%imglobe
     call MPI_Allreduce(MPI_IN_PLACE, ig, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)
     if (ig < ihuge) wgts%imglobe = ig

  endif
#endif

end subroutine get_weights_lonlat_point

!=========================================================================

subroutine get_weights_lonlat_1d(glon, glat, wgts, n, xe, ye, ze, iskip)

  use  ll_bins,    only: gridcells_from_latlon_bins
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use  misc_coms,  only: iparallel
  use  mpi
#endif

  implicit none

  integer,             intent(in)  :: n
  real,                intent(in)  :: glon(n)
  real,                intent(in)  :: glat(n)
  type(minterp_wghts), intent(out) :: wgts(n)
  real,      optional, intent(in)  :: xe(n), ye(n), ze(n)
  logical,   optional, intent(in)  :: iskip(n)

  integer, pointer                 :: nm
  integer, pointer, contiguous     :: impts(:)
  integer                          :: i, j
  real                             :: rads, raxis, xeg, yeg, zeg
  real                             :: cosglat, singlat
  real                             :: cosglon, singlon

#ifdef OLAM_MPI
  integer, allocatable             :: ig(:)
  integer                          :: ierr
#endif

  !$omp parallel do private(i,rads,cosglat,cosglon,singlat,singlon, &
  !$                        raxis,xeg,yeg,zeg,nm,impts,j)
  do i = 1, n
     wgts(i)%imglobe = 1

     if ( present(iskip) ) then
        if (iskip(i)) cycle
     endif

     call gridcells_from_latlon_bins(glat(i), glon(i), bsetm, nm, impts)

     if (present(xe) .and. present(ye) .and. present(ze)) then

        xeg = xe(i)
        yeg = ye(i)
        zeg = ze(i)

     else

        rads    = glat(i) * pio180
        cosglat = cos(rads)
        singlat = sin(rads)

        rads    = glon(i) * pio180
        cosglon = cos(rads)
        singlon = sin(rads)

        raxis = erad  * cosglat
        zeg   = erad  * singlat
        xeg   = raxis * cosglon
        yeg   = raxis * singlon

     endif

     do j = 1, nm
        call mweights_xyze(impts(j), xeg, yeg, zeg, wgts(i))
        if (wgts(i)%imglobe > 1) exit
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then

     allocate(ig(n)) ; ig = ihuge
     do i = 1, n
        if (wgts(i)%imglobe > 1) ig(i) = wgts(i)%imglobe
     enddo

     call MPI_Allreduce(MPI_IN_PLACE, ig, n, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)

     do i = 1, n
        if (ig(i) < ihuge) wgts(i)%imglobe = ig(i)
     enddo

  endif
#endif

end subroutine get_weights_lonlat_1d

!=========================================================================

subroutine get_weights_lonlat_2d(glon, glat, wgts, n1, n2, xe, ye, ze, iskip)

  use  ll_bins,    only: gridcells_from_latlon_bins
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use  misc_coms,  only: iparallel
  use  mpi
#endif

  implicit none

  integer,             intent(in)  :: n1, n2
  real,                intent(in)  :: glon(n1,n2)
  real,                intent(in)  :: glat(n1,n2)
  type(minterp_wghts), intent(out) :: wgts(n1,n2)
  real,      optional, intent(in)  :: xe(n1,n2), ye(n1,n2), ze(n1,n2)
  logical,   optional, intent(in)  :: iskip(n1,n2)

  integer, pointer                 :: nm
  integer, pointer, contiguous     :: impts(:)
  integer                          :: in1, in2, j
  real                             :: rads, raxis, xeg, yeg, zeg
  real                             :: cosglat, singlat
  real                             :: cosglon, singlon

#ifdef OLAM_MPI
  integer, allocatable             :: ig (:,:)
  integer                          :: ierr
#endif

  !$omp parallel do collapse(2) private(in2,in1,nm,impts,xeg,yeg,zeg,rads, &
  !$                                    raxis,cosglat,singlat,cosglon,singlon,j)
  do in2 = 1, n2
     do in1 = 1, n1
        wgts(in1,in2)%imglobe = 1

        if ( present(iskip) ) then
           if (iskip(in1,in2)) cycle
        endif

        call gridcells_from_latlon_bins(glat(in1,in2), glon(in1,in2), bsetm, nm, impts)

        if (present(xe) .and. present(ye) .and. present(ze)) then

           xeg = xe(in1,in2)
           yeg = ye(in1,in2)
           zeg = ze(in1,in2)

        else

           rads    = glat(in1,in2) * pio180
           cosglat = cos(rads)
           singlat = sin(rads)

           rads    = glon(in1,in2) * pio180
           cosglon = cos(rads)
           singlon = sin(rads)

           raxis = erad  * cosglat
           zeg   = erad  * singlat
           xeg   = raxis * cosglon
           yeg   = raxis * singlon

        endif

        do j = 1, nm
           call mweights_xyze(impts(j), xeg, yeg, zeg, wgts(in1,in2))
           if (wgts(in1,in2)%imglobe > 1) exit
        enddo

     enddo
  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then

     allocate(ig(n1,n2)) ; ig = ihuge
     do in2 = 1, n2
        do in1 = 1, n1
           if (wgts(in1,in2)%imglobe > 1) ig(in1,in2) = wgts(in1,in2)%imglobe
        enddo
     enddo

     call MPI_Allreduce(MPI_IN_PLACE, ig, n1*n2, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)

     do in2 = 1, n2
        do in1 = 1, n1
           if (ig(in1,in2) < ihuge) wgts(in1,in2)%imglobe = ig(in1,in2)
        enddo
     enddo

  endif
#endif

end subroutine get_weights_lonlat_2d

!=========================================================================

subroutine get_weights_lonlat_2db(glon, glat, wgts, n1, n2, xe, ye, ze, iskip)

  use  ll_bins,    only: gridcells_from_latlon_bins
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use  misc_coms,  only: iparallel
  use  mpi
#endif

  implicit none

  integer,             intent(in)  :: n1, n2
  real,                intent(in)  :: glon(n1,n2)
  real,                intent(in)  :: glat(n1,n2)
  type(minterp_wghts), intent(out) :: wgts(n1*n2)
  real,      optional, intent(in)  :: xe(n1,n2), ye(n1,n2), ze(n1,n2)
  logical,   optional, intent(in)  :: iskip(n1,n2)

  integer, pointer                 :: nm
  integer, pointer, contiguous     :: impts(:)
  integer                          :: in1, in2, ii, j
  real                             :: rads, raxis, xeg, yeg, zeg
  real                             :: cosglat, singlat
  real                             :: cosglon, singlon

#ifdef OLAM_MPI
  integer, allocatable             :: ig (:)
  integer                          :: ierr
#endif

  !$omp parallel do collapse(2) private(in2,in1,ii,nm,impts,xeg,yeg,zeg,rads, &
  !$                                    raxis,cosglat,singlat,cosglon,singlon,j)
  do in2 = 1, n2
     do in1 = 1, n1
        ii = (in2-1)*n1+in1

        wgts(ii)%imglobe = 1

        if ( present(iskip) ) then
           if (iskip(in1,in2)) cycle
        endif

        call gridcells_from_latlon_bins(glat(in1,in2), glon(in1,in2), bsetm, nm, impts)

        if (present(xe) .and. present(ye) .and. present(ze)) then

           xeg = xe(in1,in2)
           yeg = ye(in1,in2)
           zeg = ze(in1,in2)

        else

           rads    = glat(in1,in2) * pio180
           cosglat = cos(rads)
           singlat = sin(rads)

           rads    = glon(in1,in2) * pio180
           cosglon = cos(rads)
           singlon = sin(rads)

           raxis = erad  * cosglat
           zeg   = erad  * singlat
           xeg   = raxis * cosglon
           yeg   = raxis * singlon

        endif

        do j = 1, nm
           call mweights_xyze(impts(j), xeg, yeg, zeg, wgts(ii))
           if (wgts(ii)%imglobe > 1) exit
        enddo

     enddo
  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then

     allocate(ig(n1*n2)) ; ig = ihuge
     do ii = 1, n1*n2
        if (wgts(ii)%imglobe > 1) ig(ii) = wgts(ii)%imglobe
     enddo

     call MPI_Allreduce(MPI_IN_PLACE, ig, n1*n2, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)

     do ii = 1, n1*n2
        if (ig(ii) < ihuge) wgts(ii)%imglobe = ig(ii)
     enddo

  endif
#endif

end subroutine get_weights_lonlat_2db

!=========================================================================

subroutine get_weights_lonlat_grid(glon, glat, wgts, nlon, nlat, xe, ye, ze)

  use  ll_bins,    only: gridcells_from_latlon_bins
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use  misc_coms,  only: iparallel
  use  mpi
#endif

  implicit none

  integer,             intent(in)  :: nlat, nlon
  real,                intent(in)  :: glon(nlon)
  real,                intent(in)  :: glat(nlat)
  type(minterp_wghts), intent(out) :: wgts(nlon,nlat)
  real,      optional, intent(in)  :: xe(nlon,nlat), ye(nlon,nlat), ze(nlon,nlat)

  integer, pointer                 :: nm
  integer, pointer, contiguous     :: impts(:)
  integer                          :: ilon, ilat, j
  real                             :: rads, raxis, xeg, yeg, zeg
  real                             :: cosglat(nlat), singlat(nlat)
  real                             :: cosglon(nlon), singlon(nlon)

#ifdef OLAM_MPI
  integer, allocatable             :: ig(:,:)
  integer                          :: ierr
#endif

  if (.not. (present(xe) .and. present(ye) .and. present(ze))) then

     do ilat = 1, nlat
        rads          = glat(ilat) * pio180
        cosglat(ilat) = cos(rads)
        singlat(ilat) = sin(rads)
     enddo

     do ilon = 1, nlon
        rads          = glon(ilon) * pio180
        cosglon(ilon) = cos(rads)
        singlon(ilon) = sin(rads)
     enddo

  endif

  !$omp parallel do private(zeg,raxis,ilon,xeg,yeg,nm,impts,j)
  do ilat = 1, nlat

     if (.not. (present(xe) .and. present(ye) .and. present(ze))) then
        zeg   = erad * singlat(ilat)
        raxis = erad * cosglat(ilat)
     endif

     do ilon = 1, nlon

        if (present(xe) .and. present(ye) .and. present(ze)) then
           xeg = xe(ilon,ilat)
           yeg = ye(ilon,ilat)
           zeg = ze(ilon,ilat)
        else
           xeg = raxis * cosglon(ilon)
           yeg = raxis * singlon(ilon)
        endif

        wgts(ilon,ilat)%imglobe = 1

        call gridcells_from_latlon_bins(glat(ilat), glon(ilon), bsetm, nm, impts)

        do j = 1, nm
           call mweights_xyze(impts(j), xeg, yeg, zeg, wgts(ilon,ilat))
           if (wgts(ilon,ilat)%imglobe > 1) exit
        enddo
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then

     allocate(ig(nlon,nlat)) ; ig = ihuge
     do ilat = 1, nlat
        do ilon = 1, nlon
           if (wgts(ilon,ilat)%imglobe > 1) ig(ilon,ilat) = wgts(ilon,ilat)%imglobe
        enddo
     enddo

     call MPI_Allreduce(MPI_IN_PLACE, ig, nlon*nlat, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)

     do ilat = 1, nlat
        do ilon = 1, nlon
           if (ig(ilon,ilat) < ihuge) wgts(ilon,ilat)%imglobe = ig(ilon,ilat)
        enddo
     enddo

  endif
#endif

end subroutine get_weights_lonlat_grid

!=========================================================================

subroutine get_weights_lonlat_gridb(glon, glat, wgts, nlon, nlat, xe, ye, ze)

  use  ll_bins,    only: gridcells_from_latlon_bins
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use  misc_coms,  only: iparallel
  use  mpi
#endif

  implicit none

  integer,             intent(in)  :: nlat, nlon
  real,                intent(in)  :: glon(nlon)
  real,                intent(in)  :: glat(nlat)
  type(minterp_wghts), intent(out) :: wgts(nlon*nlat)
  real,      optional, intent(in)  :: xe(nlon,nlat), ye(nlon,nlat), ze(nlon,nlat)

  integer, pointer                 :: nm
  integer, pointer, contiguous     :: impts(:)
  integer                          :: ilon, ilat, ii, j
  real                             :: rads, raxis, xeg, yeg, zeg
  real                             :: cosglat(nlat), singlat(nlat)
  real                             :: cosglon(nlon), singlon(nlon)

#ifdef OLAM_MPI
  integer, allocatable             :: ig(:)
  integer                          :: ierr
#endif

  if (.not. (present(xe) .and. present(ye) .and. present(ze))) then

     do ilat = 1, nlat
        rads          = glat(ilat) * pio180
        cosglat(ilat) = cos(rads)
        singlat(ilat) = sin(rads)
     enddo

     do ilon = 1, nlon
        rads          = glon(ilon) * pio180
        cosglon(ilon) = cos(rads)
        singlon(ilon) = sin(rads)
     enddo

  endif

  !$omp parallel do private(zeg,raxis,ilon,xeg,yeg,nm,impts,j)
  do ilat = 1, nlat

     if (.not. (present(xe) .and. present(ye) .and. present(ze))) then
        zeg   = erad * singlat(ilat)
        raxis = erad * cosglat(ilat)
     endif

     do ilon = 1, nlon

        ii = (ilat-1)*nlon + ilon

        if (present(xe) .and. present(ye) .and. present(ze)) then
           xeg = xe(ilon,ilat)
           yeg = ye(ilon,ilat)
           zeg = ze(ilon,ilat)
        else
           xeg = raxis * cosglon(ilon)
           yeg = raxis * singlon(ilon)
        endif

        wgts(ii)%imglobe = 1

        call gridcells_from_latlon_bins(glat(ilat), glon(ilon), bsetm, nm, impts)

        do j = 1, nm
           call mweights_xyze(impts(j), xeg, yeg, zeg, wgts(ii))
           if (wgts(ii)%imglobe > 1) exit
        enddo
     enddo

  enddo
  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then

     allocate(ig(nlon*nlat)) ; ig = ihuge
     do ii = 1, nlon*nlat
        if (wgts(ii)%imglobe > 1) ig(ii) = wgts(ii)%imglobe
     enddo

     call MPI_Allreduce(MPI_IN_PLACE, ig, nlon*nlat, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)

     do ii = 1, nlon*nlat
        if (ig(ii) < ihuge) wgts(ii)%imglobe = ig(ii)
     enddo

  endif
#endif

end subroutine get_weights_lonlat_gridb

!=========================================================================

subroutine mweights_xyze(im, xe, ye, ze, w)

  use mem_grid,   only: xem, yem, zem
  use mem_ijtabs, only: itab_m

  implicit none

  integer,             intent(in ) :: im
  real,                intent(in ) :: xe, ye, ze
  type(minterp_wghts), intent(out) :: w

  real :: dxe, dye, dze
  real :: xp, yp, dxp3, dyp3
  real :: wgts(3)

  if (xe < minterp(im)%xemin .or. xe > minterp(im)%xemax) return
  if (ye < minterp(im)%yemin .or. ye > minterp(im)%yemax) return
  if (ze < minterp(im)%zemin .or. ze > minterp(im)%zemax) return

  dxe = xe - xem(im)
  dye = ye - yem(im)
  dze = ze - zem(im)

  call de_gn( dxe, dye, dze, &
              minterp(im)%cosmlat, minterp(im)%sinmlat, &
              minterp(im)%cosmlon, minterp(im)%sinmlon, xp, yp )

  dxp3 = xp - minterp(im)%xg3
  dyp3 = yp - minterp(im)%yg3

  wgts(1) = dxp3 * minterp(im)%w1factx + dyp3 * minterp(im)%w1facty
  if (wgts(1) < -5.e-6) return

  wgts(2) = dxp3 * minterp(im)%w2factx + dyp3 * minterp(im)%w2facty
  if (wgts(2) < -5.e-6) return

  wgts(3) = 1.0 - wgts(1) - wgts(2)
  if (wgts(3) < -5.e-6) return

  ! if we arrived here, point is in/on the current triangle

  w%imglobe = itab_m(im)%imglobe
  w%wt      = wgts

end subroutine mweights_xyze

!=========================================================================

subroutine mweights_lonlat(im, qlon, qlat, w)
  implicit none

  integer,             intent(in ) :: im
  real,                intent(in ) :: qlon, qlat
  type(minterp_wghts), intent(out) :: w

  real :: xe, ye, ze

  call ec_e(qlon, qlat, xe, ye, ze)

  call mweights_xyze(im, xe, ye, ze, w)

end subroutine mweights_lonlat

!=========================================================================

subroutine interp_level_1(npts, wts, fieldin, fieldout, docomm)

  use mem_grid,     only: mwa
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: itab_m, itabg_m
  use mem_para,     only: myrank

  implicit none

  integer,             intent(in)    :: npts
  type(minterp_wghts), intent(in)    :: wts     (npts)
  real,                intent(inout) :: fieldin (mwa)
  real,                intent(inout) :: fieldout(npts)
  logical,   optional, intent(in)    :: docomm

  integer :: ipt, im
  logical :: dompi

  if (iparallel == 1) then
     dompi = .false.
     if (present(docomm)) dompi = docomm

     if (dompi) then
        call mpi_send_w(r1dvara1=fieldin)
        call mpi_recv_w(r1dvara1=fieldin)
     endif
  endif

  !$omp parallel do private(ipt,im)
  do ipt = 1, npts
     if (itabg_m(wts(ipt)%imglobe)%irank == myrank) then

        im = itabg_m(wts(ipt)%imglobe)%im_myrank
        fieldout(ipt) = wts(ipt)%wt(1) * fieldin( itab_m(im)%iw(1) ) &
                      + wts(ipt)%wt(2) * fieldin( itab_m(im)%iw(2) ) &
                      + wts(ipt)%wt(3) * fieldin( itab_m(im)%iw(3) )

     endif
  enddo
  !$omp end parallel do

end subroutine interp_level_1

!=========================================================================

subroutine interp_level_2(n1, n2, wts, fieldin, fieldout, docomm)

  use mem_grid,     only: mwa
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: itab_m, itabg_m
  use mem_para,     only: myrank

  implicit none

  integer,             intent(in)    :: n1, n2
  type(minterp_wghts), intent(in)    :: wts     (n1,n2)
  real,                intent(inout) :: fieldin (mwa)
  real,                intent(inout) :: fieldout(n1,n2)
  logical,   optional, intent(in)    :: docomm

  integer :: j, i, im
  logical :: dompi

  if (iparallel == 1) then
     dompi = .false.
     if (present(docomm)) dompi = docomm

     if (dompi) then
        call mpi_send_w(r1dvara1=fieldin)
        call mpi_recv_w(r1dvara1=fieldin)
     endif
  endif

  !$omp parallel do collapse(2) private(j,i,im)
  do j = 1, n2
     do i = 1, n1
        if (itabg_m(wts(i,j)%imglobe)%irank == myrank) then

           im = itabg_m(wts(i,j)%imglobe)%im_myrank
           fieldout(i,j) = wts(i,j)%wt(1) * fieldin( itab_m(im)%iw(1) ) &
                         + wts(i,j)%wt(2) * fieldin( itab_m(im)%iw(2) ) &
                         + wts(i,j)%wt(3) * fieldin( itab_m(im)%iw(3) )

        endif
     enddo
  enddo
  !$omp end parallel do

end subroutine interp_level_2

!================================================================================

subroutine interp_column_1(npts, wts, nlevin, nlevout, fieldin, fieldout, docomm, set_bottom)

  use mem_grid,     only: mwa, mza, lpw
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: jtw_prog, itab_m, itabg_m, itab_w
  use mem_para,     only: myrank

  implicit none

  integer,             intent(in)    :: npts, nlevin, nlevout
  type(minterp_wghts), intent(in)    :: wts     (npts)
  real,                intent(inout) :: fieldin (nlevin,mwa)
  real,                intent(inout) :: fieldout(npts,nlevout)
  logical,   optional, intent(in)    :: docomm
  logical,   optional, intent(in)    :: set_bottom

  integer :: kin, kout, ipt, k, iw, im
  logical :: dompi, dobot

  dompi = .false.
  if (iparallel == 1) then
     if (present(docomm)) dompi = docomm
  endif

  dobot = .false.
  if (nlevin == mza) then
     if (present(set_bottom)) dobot = set_bottom
  endif

  if (dobot) then
     !$omp parallel do private(iw,k)
     do iw = 2, mwa
        if (dompi .and. itab_w(iw)%irank /= myrank) cycle
        do k = lpw(iw)-1, 1, -1
           fieldin(k,iw) = fieldin(lpw(iw),iw)
        enddo
     enddo
     !$omp end parallel do
  endif

  if (dompi) then
     call mpi_send_w(svara1=fieldin)
     call mpi_recv_w(svara1=fieldin)
  endif

  !$omp parallel do private(ipt,im,kout,kin)
  do ipt = 1, npts
     if (itabg_m(wts(ipt)%imglobe)%irank == myrank) then

        im = itabg_m(wts(ipt)%imglobe)%im_myrank
        do kout = 1, nlevout
           kin = kout + nlevin - nlevout
           fieldout(ipt,kout) = wts(ipt)%wt(1) * fieldin( kin, itab_m(im)%iw(1) ) &
                              + wts(ipt)%wt(2) * fieldin( kin, itab_m(im)%iw(2) ) &
                              + wts(ipt)%wt(3) * fieldin( kin, itab_m(im)%iw(3) )
        enddo

     endif
  enddo
  !$omp end parallel do

end subroutine interp_column_1

!================================================================================

subroutine interp_column_2(n1, n2, wts, nlevin, nlevout, fieldin, fieldout, docomm, set_bottom)

  use mem_grid,     only: mwa, mza, lpw
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: jtw_prog, itab_m, itabg_m, itab_w
  use mem_para,     only: myrank

  implicit none

  integer,             intent(in)    :: n1, n2, nlevin, nlevout
  type(minterp_wghts), intent(in)    :: wts     (n1,n2)
  real,                intent(inout) :: fieldin (nlevin,mwa)
  real,                intent(inout) :: fieldout(n1,n2,nlevout)
  logical,   optional, intent(in)    :: docomm
  logical,   optional, intent(in)    :: set_bottom

  integer :: kin, kout, i, j, k, iw, im
  logical :: dompi, dobot

  dompi = .false.
  if (iparallel == 1) then
     if (present(docomm)) dompi = docomm
  endif

  dobot = .false.
  if (nlevin == mza) then
     if (present(set_bottom)) dobot = set_bottom
  endif

  if (dobot) then
     !$omp parallel do private(iw,k)
     do iw = 2, mwa
        if (dompi .and. itab_w(iw)%irank /= myrank) cycle
        do k = lpw(iw)-1, 1, -1
           fieldin(k,iw) = fieldin(lpw(iw),iw)
        enddo
     enddo
     !$omp end parallel do
  endif

  if (dompi) then
     call mpi_send_w(svara1=fieldin)
     call mpi_recv_w(svara1=fieldin)
  endif

  !$omp parallel do collapse(2) private(j,i,im,kout,kin)
  do j = 1, n2
     do i = 1, n1
        if (itabg_m(wts(i,j)%imglobe)%irank == myrank) then

           im = itabg_m(wts(i,j)%imglobe)%im_myrank
           do kout = 1, nlevout
              kin = kout + nlevin - nlevout
              fieldout(i,j,kout) = wts(i,j)%wt(1) * fieldin( kin, itab_m(im)%iw(1) ) &
                                 + wts(i,j)%wt(2) * fieldin( kin, itab_m(im)%iw(2) ) &
                                 + wts(i,j)%wt(3) * fieldin( kin, itab_m(im)%iw(3) )
           enddo

        endif
     enddo
  enddo
  !$omp end parallel do

end subroutine interp_column_2

!================================================================================

subroutine interp_column_tolev_1(npts, wts, nlevin, fieldin, fieldout, lev, docomm, set_bottom)

  use mem_grid,     only: mwa, mza, lpw
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: jtw_prog, itab_m, itabg_m, itab_w
  use mem_para,     only: myrank

  implicit none

  integer,             intent(in)    :: npts, nlevin, lev
  type(minterp_wghts), intent(in)    :: wts     (npts)
  real,                intent(inout) :: fieldin (nlevin,mwa)
  real,                intent(inout) :: fieldout(npts)
  logical,   optional, intent(in)    :: docomm
  logical,   optional, intent(in)    :: set_bottom

  integer :: ipt, k, iw, im
  logical :: dompi, dobot
  real    :: field2d(mwa)

  dompi = .false.
  if (iparallel == 1) then
     if (present(docomm)) dompi = docomm
  endif

  dobot = .false.
  if (nlevin == mza) then
     if (present(set_bottom)) dobot = set_bottom
  endif

  !$omp parallel do private(iw,k)
  do iw = 2, mwa
     if (dompi .and. itab_w(iw)%irank /= myrank) cycle
     k = lev
     if (dobot) k = max(k,lpw(iw))
     field2d(iw) = fieldin(k,iw)
  enddo
  !$omp end parallel do

  if (dompi) then
     call mpi_send_w(r1dvara1=field2d)
     call mpi_recv_w(r1dvara1=field2d)
  endif

  !$omp parallel do private(ipt,im)
  do ipt = 1, npts
     if (itabg_m(wts(ipt)%imglobe)%irank == myrank) then

        im = itabg_m(wts(ipt)%imglobe)%im_myrank
        fieldout(ipt) = wts(ipt)%wt(1) * field2d( itab_m(im)%iw(1) ) &
                      + wts(ipt)%wt(2) * field2d( itab_m(im)%iw(2) ) &
                      + wts(ipt)%wt(3) * field2d( itab_m(im)%iw(3) )
     endif
  enddo
  !$omp end parallel do

end subroutine interp_column_tolev_1

!================================================================================

subroutine interp_column_tolev_2(n1, n2, wts, nlevin, fieldin, fieldout, lev, docomm, set_bottom)

  use mem_grid,     only: mwa, mza, lpw
  use misc_coms,    only: iparallel
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w
  use mem_ijtabs,   only: jtw_prog, itab_m, itabg_m, itab_w
  use mem_para,     only: myrank

  implicit none

  integer,             intent(in)    :: n1, n2, nlevin, lev
  type(minterp_wghts), intent(in)    :: wts     (n1,n2)
  real,                intent(inout) :: fieldin (nlevin,mwa)
  real,                intent(inout) :: fieldout(n1,n2)
  logical,   optional, intent(in)    :: docomm
  logical,   optional, intent(in)    :: set_bottom

  integer :: j, i, k, iw, im
  logical :: dompi, dobot
  real    :: field2d(mwa)

  dompi = .false.
  if (iparallel == 1) then
     if (present(docomm)) dompi = docomm
  endif

  dobot = .false.
  if (nlevin == mza) then
     if (present(set_bottom)) dobot = set_bottom
  endif

  !$omp parallel do private(iw,k)
  do iw = 2, mwa
     if (dompi .and. itab_w(iw)%irank /= myrank) cycle
     k = lev
     if (dobot) k = max(k,lpw(iw))
     field2d(iw) = fieldin(k,iw)
  enddo
  !$omp end parallel do

  if (dompi) then
     call mpi_send_w(r1dvara1=field2d)
     call mpi_recv_w(r1dvara1=field2d)
  endif

  !$omp parallel do collapse(2) private(j,i,im)
  do j = 1, n2
     do i = 1, n2
        if (itabg_m(wts(i,j)%imglobe)%irank == myrank) then

           im = itabg_m(wts(i,j)%imglobe)%im_myrank
           fieldout(i,j) = wts(i,j)%wt(1) * field2d( itab_m(im)%iw(1) ) &
                         + wts(i,j)%wt(2) * field2d( itab_m(im)%iw(2) ) &
                         + wts(i,j)%wt(3) * field2d( itab_m(im)%iw(3) )
        endif
     enddo
  enddo
  !$omp end parallel do

end subroutine interp_column_tolev_2

!================================================================================

end module minterp_lib
