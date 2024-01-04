module area_remap_ll

  ! This module contains routines for conservatively remapping a global
  ! lat/lon field to the olam mesh using area overlaps. The global integral
  ! of the input field should be conserved. These routines are useful for
  ! interpolating emissions, flux, and surface datasets to the olam grid.

  type overlap_vars
     integer              :: ncells = 0
     integer, allocatable :: i(:)
     integer, allocatable :: j(:)
     real,    allocatable :: area(:)
  end type overlap_vars

contains

!============================================================================

subroutine get_ll_overlaps_atm( nlon, nlat, dlon, dlat, overlaps, tolerance, &
                                  swlon, swlat, lat_beg, lat_end, irev_ns )

  use mem_ijtabs, only: jtab_w, itab_w, jtw_prog
  use mem_grid,   only: mwa, glatw, glonw, glatm, glonm, arw0, &
                        xem, yem, zem, xew, yew, zew

  integer,            intent(in)  :: nlon, nlat
  real,               intent(in)  :: dlon, dlat
  type(overlap_vars), intent(out) :: overlaps(mwa)
  real,     optional, intent(in)  :: tolerance
  real,     optional, intent(in)  :: swlon
  real,     optional, intent(in)  :: swlat
  real,     optional, intent(in)  :: lat_beg
  real,     optional, intent(in)  :: lat_end
  logical,  optional, intent(in)  :: irev_ns

  integer, allocatable :: ii(:) ! temporary storage array
  integer, allocatable :: jj(:) ! temporary storage array
  real,    allocatable :: aa(:) ! temporary storage array

  integer :: j, iw, np, n, im
  real    :: alat_beg, alat_end

  integer, parameter :: maxvert = 7

  real :: lats(maxvert), lons(maxvert)
  real :: xes(maxvert), yes(maxvert), zes(maxvert)

  ! avoid cells close to poles

  alat_beg = -84.0
  if (present(lat_beg)) alat_beg = lat_beg

  alat_end = 84.0
  if (present(lat_end)) alat_end = lat_end

  !$omp parallel private(ii,jj,aa)
  !$omp do private(iw,np,n,im,xes,yes,zes,lons,lats)
  do j = 1, jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Skip cells near poles; small emissions and many lat-lon overlaps
     ! Just find lat/lon point closes to olam cell center

     if (glatw(iw) < alat_beg .or. glatw(iw) > alat_end ) then

        call get_ll_nearest( glonw(iw), glatw(iw), arw0(iw), overlaps(iw), &
                             dlon, dlat, nlon, nlat, swlon=swlon, swlat=swlat, irev_ns=irev_ns )
     else

        np = itab_w(iw)%npoly

        do n = 1, np
           im = itab_w(iw)%im(n)

           xes(n) = xem(im)
           yes(n) = yem(im)
           zes(n) = zem(im)

           lons(n) = glonm(im)
           lats(n) = glatm(im)
        enddo

        call get_ll_overlaps( xew(iw), yew(iw), zew(iw), arw0(iw), &
                              np, lons, lats, xes, yes, zes, &
                              dlon, dlat, nlon, nlat, &
                              ii, jj, aa, overlaps(iw), &
                              swlon=swlon, swlat=swlat, irev_ns=irev_ns, tolerance=tolerance )
     endif

  enddo
  !$omp end do
  if (allocated(ii)) deallocate(ii)
  if (allocated(jj)) deallocate(jj)
  if (allocated(aa)) deallocate(aa)
  !$omp end parallel

end subroutine get_ll_overlaps_atm

!============================================================================

subroutine get_ll_overlaps_sfc( nlon, nlat, dlon, dlat, overlaps, tolerance, &
                                swlon, swlat, lat_beg, lat_end, irev_ns, primary_only, &
                                iland, isea, ilake )

  use mem_sfcg,  only: mwsfc, sfcg, itab_wsfc
  use mem_land,  only: mland, omland
  use mem_sea,   only: msea, omsea
  use mem_lake,  only: mlake, omlake
  use misc_coms, only: iparallel
  use mem_para,  only: myrank

  integer,            intent(in)  :: nlon, nlat
  real,               intent(in)  :: dlon, dlat
  type(overlap_vars), intent(out) :: overlaps(:)
  real,     optional, intent(in)  :: tolerance
  real,     optional, intent(in)  :: swlon
  real,     optional, intent(in)  :: swlat
  real,     optional, intent(in)  :: lat_beg
  real,     optional, intent(in)  :: lat_end
  logical,  optional, intent(in)  :: irev_ns
  logical,  optional, intent(in)  :: primary_only
  logical,  optional, intent(in)  :: iland
  logical,  optional, intent(in)  :: isea
  logical,  optional, intent(in)  :: ilake

  integer, allocatable :: ii(:) ! temporary storage array
  integer, allocatable :: jj(:) ! temporary storage array
  real,    allocatable :: aa(:) ! temporary storage array

  integer :: is, iwsfc, np, n, im, iend
  real    :: alat_beg, alat_end
  logical :: prim_only

  integer, parameter :: maxvert = 7

  real :: lats(maxvert), lons(maxvert)
  real :: xes(maxvert), yes(maxvert), zes(maxvert)

  ! avoid cells close to poles

  alat_beg = -84.0
  if (present(lat_beg)) alat_beg = lat_beg

  alat_end = 84.0
  if (present(lat_end)) alat_end = lat_end

  prim_only = .true.
  if (present(primary_only)) prim_only = primary_only

  if ( count( [present(iland), present(isea), present(ilake)] ) > 1 ) then
     write(*,*) "Error in get_ll_overlaps_sfc:"
     write(*,*) "only one of iland, isea, or ilake can be specified"
     stop
  endif

  if (present(iland)) then
     iend = mland
     if (size(overlaps) /= mland) then
        write(*,*) "Error in get_ll_overlaps_sfc:"
        write(*,*) "overlaps must be dimensioned to mland if iland is specified."
        stop
     endif
  elseif (present(ilake)) then
     iend = mlake
     if (size(overlaps) /= mlake) then
        write(*,*) "Error in get_ll_overlaps_sfc:"
        write(*,*) "overlaps must be dimensioned to mlake if ilake is specified."
        stop
     endif
  elseif (present(isea)) then
     iend = msea
     if (size(overlaps) /= msea) then
        write(*,*) "Error in get_ll_overlaps_sfc:"
        write(*,*) "overlaps must be dimensioned to msea if isea is specified."
        stop
     endif
  else
     iend = mwsfc
     if (size(overlaps) /= mwsfc) then
        write(*,*) "Error in get_ll_overlaps_sfc:"
        write(*,*) "overlaps must be dimensioned to mwsfc if iland, ilake, or isea is not specified."
        stop
     endif
  endif

  !$omp parallel private(ii,jj,aa)
  !$omp do private(iwsfc,np,n,im,xes,yes,zes,lons,lats)
  do is = 2, iend

     iwsfc = is
     if (present(iland)) iwsfc = is + omland
     if (present(ilake)) iwsfc = is + omlake
     if (present(isea )) iwsfc = is + omsea

     ! Skip this cell if running in parallel and cell rank is not MYRANK

     if (prim_only) then
        if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle
     endif

     ! Skip cells near poles; small emissions and many lat-lon overlaps
     ! Just find lat/lon point closes to olam cell center

     if ( sfcg%glatw(iwsfc) < alat_beg .or. sfcg%glatw(iwsfc) > alat_end ) then

        call get_ll_nearest( sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), sfcg%area(iwsfc), overlaps(is), &
                             dlon, dlat, nlon, nlat, swlon=swlon, swlat=swlat, irev_ns=irev_ns )
     else

        np = itab_wsfc(iwsfc)%npoly

        do n = 1, np
           im = itab_wsfc(iwsfc)%imn(n)

           xes(n) = sfcg%xem(im)
           yes(n) = sfcg%yem(im)
           zes(n) = sfcg%zem(im)

           lons(n) = sfcg%glonm(im)
           lats(n) = sfcg%glatm(im)
        enddo

        call get_ll_overlaps( sfcg%xew(iwsfc), sfcg%yew(iwsfc), sfcg%zew(iwsfc), sfcg%area(iwsfc), &
                              np, lons, lats, xes, yes, zes, &
                              dlon, dlat, nlon, nlat, &
                              ii, jj, aa, overlaps(is), &
                              swlon=swlon, swlat=swlat, irev_ns=irev_ns, tolerance=tolerance )
     endif

  enddo
  !$omp end do
  if (allocated(ii)) deallocate(ii)
  if (allocated(jj)) deallocate(jj)
  if (allocated(aa)) deallocate(aa)
  !$omp end parallel

end subroutine get_ll_overlaps_sfc

!============================================================================

subroutine get_ll_overlaps( xew, yew, zew, arw0,              &
                            np, glonm, glatm, xem, yem, zem,  &
                            dlon, dlat, nlon, nlat,           &
                            ipoints, jpoints, areafrc,        &
                            overlaps, swlon, swlat, irev_ns, tolerance )

  use, intrinsic :: iso_fortran_env, only: r8=>real64
  use               consts_coms,     only: eradsq, pio180
  use               polygon_lib,     only: polygon_overlap

  implicit none

  real,                 intent(in)    :: xew, yew, zew              ! Earth cart coord of OLAM cell center
  real,                 intent(in)    :: arw0
  integer,              intent(in)    :: np                         ! # of OLAM cell vertices
  real,                 intent(in)    :: glatm(np), glonm(np)       ! lat/lon of OLAM cell vertices
  real,                 intent(in)    :: xem(np), yem(np), zem(np)  ! Earth cart coord of OLAM cell vertices

  real,                 intent(in)    :: dlon, dlat    ! lon/lat spacing
  integer,              intent(in)    :: nlon, nlat    ! # of lons/lats

  integer, allocatable, intent(inout) :: ipoints(:) ! temporary storage array
  integer, allocatable, intent(inout) :: jpoints(:) ! temporary storage array
  real,    allocatable, intent(inout) :: areafrc(:) ! temporary storage array

  type (overlap_vars),  intent(out)   :: overlaps

  real,       optional, intent(in)    :: swlon, swlat ! coordinates of the center of the lowest left cell
  real,       optional, intent(in)    :: tolerance
  logical,    optional, intent(in)    :: irev_ns

  real :: coswlon, sinwlon, coswlat, sinwlat
  real :: glons(np)
  real :: x(max(4,np)), y(max(4,np))
  real :: dxe(np), dye(np), dze(np)

  real(r8) :: xf(np), yf(np)
  real(r8) :: xg(4), yg(4)

  real :: lons(4), lats(4)
  real :: area, areaij, atol, tol
  real :: lon0, lat0

  integer :: ng, ijsize
  integer :: i, j, n, jj
  integer :: is, ie, js, je

  integer, allocatable :: itmp(:)
  integer, allocatable :: jtmp(:)
  real,    allocatable :: atmp(:)

  lon0 = 0.
  if (present(swlon)) lon0 = swlon - 0.5*dlon

  lat0 = -90.
  if (present(swlat)) lat0 = swlat - 0.5*dlat

  if (.not. allocated(ipoints)) allocate(ipoints(1000))
  if (.not. allocated(jpoints)) allocate(jpoints(1000))
  if (.not. allocated(areafrc)) allocate(areafrc(1000))

  glons = glonm

  ! Area overlap tolerance (to avoid tiny overlaps)

  tol = 1.e-4
  if (present(tolerance)) tol = max(tolerance, 0.)

  ! "unwrap" an OLAM cell that straddles -180/180

  if (any(glons < -90.) .and. any(glons > 90.)) then
     do j = 1, np
        if (glons(j) < -90.) glons(j) = glons(j) + 360.
     enddo
  endif

  ! get starting/ending indices of lat/lon points that may overlap with this cell

  is = 1 + floor( (minval(glons(1:np)) - lon0) / dlon )
  ie = 1 + floor( (maxval(glons(1:np)) - lon0) / dlon )

  js = 1 + floor( (minval(glatm(1:np)) - lat0) / dlat )
  je = 1 + floor( (maxval(glatm(1:np)) - lat0) / dlat )

  ! Calculate the iw cell vertices on a gnomonic tangent plane
  ! centered at Earth cartesian coordinates (xew,yew,zew)

  call get_sincos_latlon(coswlon,sinwlon,coswlat,sinwlat,xew,yew,zew)

  do n = 1, np
     dxe(n) = xem(n) - xew
     dye(n) = yem(n) - yew
     dze(n) = zem(n) - zew
  enddo

  call de_gn_mult(np, dxe, dye, dze, coswlat, sinwlat, coswlon, sinwlon, x, y)

  xf = x(1:np)
  yf = y(1:np)

  ! Loop over all possible overlapping lat/lon cells to determine area overlaps

  ng = 0

  do j = js, je

     jj = j
     if (present(irev_ns)) then
        if (irev_ns) jj = nlat - j + 1
     endif

     do i = is, ie

        lons = [ real(i-1)*dlon, real(i-1)*dlon, real(i)*dlon, real(i  )*dlon ] + lon0
        lats = [ real(j-1)*dlat, real(j  )*dlat, real(j)*dlat, real(j-1)*dlat ] + lat0

        ! x,y coordinates of emissions cell on a gnomonic tangent plane
        ! centered at Earth cartesian coordinates (xew,yew,zew)

        call ll_gn_mult(4, lats, lons, coswlat, sinwlat, coswlon, sinwlon, &
                        xew, yew, zew, x, y)

        xg = x(1:4)
        yg = y(1:4)

        ! area of lat/lon cell

        areaij =  pio180 * eradsq * dlon * abs( sin(lats(2)*pio180) - sin(lats(1)*pio180) )

        ! compute any overlaps

        call polygon_overlap(np, 4, xf, yf, xg, yg, arw0, areaij, area)

        ! minimum overlap area

        atol = tol * min(arw0, areaij)

        if (area > atol) then

           ng = ng + 1
           ijsize = size(ipoints)

           if (ng > ijsize) then
              allocate(itmp(ijsize+1000))
              itmp(1:ijsize) = ipoints(1:ijsize)
              call move_alloc(itmp, ipoints)

              allocate(jtmp(ijsize+1000))
              jtmp(1:ijsize) = jpoints(1:ijsize)
              call move_alloc(jtmp, jpoints)

              allocate(atmp(ijsize+1000))
              atmp(1:ijsize) = areafrc(1:ijsize)
              call move_alloc(atmp, areafrc)
           endif

           ipoints(ng) = modulo(i-1,nlon) + 1
           jpoints(ng) = jj
           areafrc(ng) = area

        endif

     enddo
  enddo

  overlaps%ncells = ng

  if (ng > 0) then
     allocate( overlaps%i   (ng) )
     allocate( overlaps%j   (ng) )
     allocate( overlaps%area(ng) )

     overlaps%i    = ipoints(1:ng)
     overlaps%j    = jpoints(1:ng)
     overlaps%area = areafrc(1:ng)
  endif

end subroutine get_ll_overlaps

!============================================================================

subroutine get_ll_nearest( glonw, glatw, arw0, overlaps, dlon, dlat, &
                           nlon, nlat, swlon, swlat, irev_ns )

  implicit none

  real,                 intent(in)  :: glatw, glonw
  real,                 intent(in)  :: arw0

  type (overlap_vars),  intent(out) :: overlaps

  real,                 intent(in)  :: dlon, dlat
  integer,              intent(in)  :: nlon, nlat
  real,       optional, intent(in)  :: swlon, swlat ! coordinates of the center of the lowest left cell
  logical,    optional, intent(in)  :: irev_ns

  real    :: lon0, lat0
  real    :: stax, stay
  integer :: i, j

  lon0 = 0.5 * dlon
  if (present(swlon)) lon0 = swlon

  lat0 = 0.5 * dlat - 90.
  if (present(swlat)) lat0 = swlat

  stax = (glonw - lon0) / dlon
  i    = modulo(nint(stax),nlon) + 1

  stay = (glatw - lat0) / dlat
  j    = max(1, min(nlat, nint(stay) + 1))

  if (present(irev_ns)) then
     if (irev_ns) j = nlat - j + 1
  endif

  allocate( overlaps%i   (1) )
  allocate( overlaps%j   (1) )
  allocate( overlaps%area(1) )

  overlaps%ncells = 1
  overlaps%i      = i
  overlaps%j      = j
  overlaps%area   = arw0

end subroutine get_ll_nearest

!============================================================================

end module area_remap_ll
