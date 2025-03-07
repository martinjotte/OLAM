Module ll_bins

  implicit none

!!  integer, parameter :: nset = 4
!!
!!  integer, parameter :: nlat (nset) = [ 18,   90,  450, 1800 ] ! Number of bins spanning pole-to-pole
!!                                                               ! distance for each bin set
!!
!!! integer, parameter :: maxb (nset) = [ 200,  200,  200, huge(1) ] ! Max allowed population in a bin
!!                                                                   ! for each bin set
!!
!!  real,    parameter :: delat(nset) = real(nlat(:)) / 180.  ! # bins per deg of lat for each bin set
!!
!!  real,    parameter :: bres (nset) = 0.5 / delat(:)

  Type bin_vars
     integer              :: nw
     integer, allocatable :: iw(:)
  End type bin_vars

  Type blat_vars
     integer                     :: nlon    ! # bins per lon circle for each bin set and lat band
     real                        :: delon   ! # bins per deg of lon for each bin set and lat band
     type(bin_vars), allocatable :: bins(:)
  End type blat_vars

  Type binset_vars
     integer                      :: nlat    ! # bins per lat
     real                         :: delat   ! # bins per deg of lat
     type(blat_vars), allocatable :: blat(:)
  End type binset_vars

  Type itab_grid_vars
     integer :: np        ! number of bounding points listed
     real    :: glats(8)  ! lats of points that describe the bounds of this cell
     real    :: glons(8)  ! lons of points that describe the bounds of this cell
  End Type itab_grid_vars

Contains

!=========================================================================

subroutine latlon_bins(nwa, res, itab_w0, bset)

  use consts_coms, only: erad, pi1
  implicit none

  integer,              intent(in   ) :: nwa
  real,                 intent(in   ) :: res
  type(itab_grid_vars), intent(inout) :: itab_w0(nwa)
  type(binset_vars),    intent(  out) :: bset

  integer :: jsup, i, j, iw, nlato3, nlon, nw
  integer :: ix, iy
  real    :: dx, dy

  integer :: nn
  integer, pointer :: iipts(:)

  integer :: isetm, j0, j1, jj, i0, i1, ii, np
  real    :: ymin, ymax, xmin, xmax
  real    :: epsx, epsy

  real, parameter :: fuzz = 0.05

  bset%nlat  = nint( erad * pi1 / res )
  bset%delat = real(bset%nlat) / 180.

  allocate(bset%blat(bset%nlat))
  nlato3 = bset%nlat / 3

  !$omp parallel
  !$omp do private(jsup,nlon)
  do j = 1, bset%nlat
     jsup = bset%nlat + 1 - j

     nlon = 6 * min(j, jsup, nlato3)

     bset%blat(j)%nlon  = nlon
     bset%blat(j)%delon = real(nlon) / 360.

     allocate(bset%blat(j)%bins(nlon)); bset%blat(j)%bins(:)%nw = 0
  enddo
  !$opm end do

  !$omp do private(np,xmin,xmax,j,ymin,ymax,j0,j1,epsx,epsy,i0,i1,ii,i)
  do iw = 2, nwa
     if (itab_w0(iw)%np < 1) cycle  ! flag to skip this cell

     np = itab_w0(iw)%np

     if ( any( abs(itab_w0(iw)%glats(1:np)) > 89.75 ) ) then

        xmin = -180.0
        xmax =  180.0
        epsx = 0.0

     else

        ! unwrap lons if they span -180,+180

        do j = 2, np
           if (itab_w0(iw)%glons(1) < -90. .and. itab_w0(iw)%glons(j) > 90.) then
              itab_w0(iw)%glons(j) = itab_w0(iw)%glons(j) - 360.
           elseif (itab_w0(iw)%glons(1) > 90. .and. itab_w0(iw)%glons(j) < -90.) then
              itab_w0(iw)%glons(j) = itab_w0(iw)%glons(j) + 360.
           endif
        enddo

        ! Get range of lats/lons this cell spans

        xmin = minval( itab_w0(iw)%glons(1:np) )
        xmax = maxval( itab_w0(iw)%glons(1:np) )
        epsx = fuzz * (xmax - xmin)

     endif

     ymin = minval( itab_w0(iw)%glats(1:np) )
     ymax = maxval( itab_w0(iw)%glats(1:np) )
     epsy = fuzz * (ymax - ymin)

     ! Loop over all bins this cell spans

     j0 = max(1,         floor(bset%delat * (ymin - epsy + 90.)) + 1)
     j1 = min(bset%nlat, floor(bset%delat * (ymax + epsy + 90.)) + 1)

     do j = j0, j1
        i0 = floor(bset%blat(j)%delon * (xmin - epsx + 180.)) + 1
        i1 = floor(bset%blat(j)%delon * (xmax + epsx + 180.)) + 1

        i1 = min(i1,i0+bset%blat(j)%nlon-1)
        do ii = i0, i1
           i = modulo(ii-1, bset%blat(j)%nlon) + 1
           !$omp atomic
           bset%blat(j)%bins(i)%nw = bset%blat(j)%bins(i)%nw + 1
        enddo
     enddo

  enddo
  !$omp end do

!  ! Loop through all bins of all sets and allocate IW(:) array

  !$omp do private(i)
  do j = 1, bset%nlat
     do i = 1, bset%blat(j)%nlon
        allocate( bset%blat(j)%bins(i)%iw( bset%blat(j)%bins(i)%nw ) )
        bset%blat(j)%bins(i)%nw = 0
     enddo
  enddo
  !$omp end do

  ! Loop through all ATM iw cells again and sort them into the global lat-lon bins

  !$omp do private(np,xmin,xmax,j,ymin,ymax,j0,j1,epsx,epsy,i0,i1,ii,i,nw)
  do iw = 2, nwa
     if (itab_w0(iw)%np < 1) cycle  ! flag to skip this cell

     np = itab_w0(iw)%np

     if ( any( abs(itab_w0(iw)%glats(1:np)) > 89.75 ) ) then

        xmin = -180.0
        xmax =  180.0
        epsx = 0.0

     else

        ! Lons have already been unwrapped,
        ! get range of lats/lons this cell spans

        xmin = minval( itab_w0(iw)%glons(1:np) )
        xmax = maxval( itab_w0(iw)%glons(1:np) )
        epsx = fuzz * (xmax - xmin)

     endif

     ymin = minval( itab_w0(iw)%glats(1:np) )
     ymax = maxval( itab_w0(iw)%glats(1:np) )
     epsy = fuzz * (ymax - ymin)

     ! Loop over all bins this cell spans

     j0 = max(1,         floor(bset%delat * (ymin - epsy + 90.)) + 1)
     j1 = min(bset%nlat, floor(bset%delat * (ymax + epsy + 90.)) + 1)

     do j = j0, j1
        i0 = floor(bset%blat(j)%delon * (xmin - epsx + 180.)) + 1
        i1 = floor(bset%blat(j)%delon * (xmax + epsx + 180.)) + 1
        i1 = min(i1,i0+bset%blat(j)%nlon-1)

        do ii = i0, i1
           i = modulo(ii-1, bset%blat(j)%nlon) + 1

           !$omp atomic capture
           bset%blat(j)%bins(i)%nw = bset%blat(j)%bins(i)%nw + 1
           nw                      = bset%blat(j)%bins(i)%nw
           !$omp end atomic

           bset%blat(j)%bins(i)%iw( nw ) = iw
        enddo
     enddo

  enddo
  !$omp end do
  !$omp end parallel


end subroutine latlon_bins

!=========================================================================

subroutine gridcells_from_latlon_bins( glat, glon, bset, nn, iipts )

  implicit none

  real,                         intent(in)  :: glat
  real,                         intent(in)  :: glon
  type(binset_vars), target,    intent(in)  :: bset
  integer, pointer,             intent(out) :: nn
  integer, pointer, contiguous, intent(out) :: iipts(:)

  integer :: iset, j, i

  j = floor(bset%delat * (glat + 90.)) + 1
  j = max(1, min(bset%nlat, j))

  i = floor(bset%blat(j)%delon * (glon + 180.)) + 1
  i = modulo(i-1, bset%blat(j)%nlon) + 1

  nn => bset%blat(j)%bins(i)%nw

  if (nn > 0) then
     iipts => bset%blat(j)%bins(i)%iw
  else
     iipts => null()
  endif

end subroutine gridcells_from_latlon_bins

end module ll_bins
