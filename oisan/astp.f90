subroutine pressure_stage()

  use isan_coms,   only: nprx, npry, nprz, gdatdy, gdatdx, pnpr, inproj, &
                         plats, xswlat, xswlon, irev_ns, o_press, o_theta, &
                         o_rho, o_rrw, o_uzonal, o_umerid, o_ozone, pbc, z_pbc
  use hdf5_utils,  only: shdf5_info, shdf5_irec, shdf5_close
  use misc_coms,   only: i_o3, io6, runtype
  use mem_ijtabs,  only: jtab_w, jtw_init, itab_w
  use mem_zonavg,  only: zonz, zont, zonr, zonu, zono, zonp_vect
  use consts_coms, only: p00, rocp, eps_vap, cvocp, p00kord, eps_vapi, t00
  use therm_lib,   only: eslf, esif
  use mem_grid,    only: mza, mwa, zt, glatw, glonw
  use analysis_lib,only: gdtost_ll
  use micro_coms,  only: miclevel

  implicit none

  real, parameter :: mwair  = 28.9628             ! molecular weight of air
  real, parameter :: mwo3   = 48.0                ! molecular weight of ozone
  real, parameter :: cnvto3 = mwair / mwo3 * 1.e6 ! ozone mixing ratio to ppmV

  logical, parameter :: log_interp_p = .true.  ! F: linearly interpolate pressure vertically
                                               ! T: logarithmically interpolate pressure vertically

  logical, parameter :: log_interp_q = .false. ! F: linearly interpolate mixing ratio vertically
                                               ! T: logarithmically interpolate mixing ratio vertically

  real :: dprat1, dprat2, dg, vapor_press, airtemp, frac

  integer :: k, j, iw, iz, kbc
  integer :: lzon_bot, kzonoff, npd
  integer :: levp, nprz_rh, nbot_o3
  integer :: ndims, idims(3), chunkdims(3)
  logical :: isrh, rh_is_percent, sh_is_gkg

  real, allocatable :: pcol_p     (:)
  real, allocatable :: pcol_exneri(:)
  real, allocatable :: pvect      (:)
  real, allocatable :: plog       (:)

  character(10) :: varname  = ' '
  character(10) :: vnams(5) = ' '

  ! Arrays to store input lat/lon/pressure data

  real, allocatable :: a2d   (:,:)
  real, allocatable :: pcol_z(:,:)
  real, allocatable :: field (:,:)

  ! Find analysis pressure level close to 500 mb

  do kbc = 1, nprz-1
     if (pnpr(kbc) < 50100.) exit
  enddo
  pbc = pnpr(kbc)

  ! Determine index of lowest ZONAVG pressure level that is at least 1/2
  ! ZONAVG pressure level higher than highest input pressure data level
  ! (i.e., maximum zonp_vect value that is less than 82.5% of pnpr(nprz),
  ! which is in hPa)

  lzon_bot = min(23, nint(31. - 6. *  log10( pnpr(nprz) )) + 1)
  npd      = nprz + 25 - lzon_bot
  kzonoff  = npd - 22

  allocate( pcol_p     (npd) )
  allocate( pcol_exneri(npd) )
  allocate( pvect      (npd) )
  allocate( plog       (npd) )

  ! Fill column array of pressure level values in MKS

  pcol_p(1) = 120000.      ! Phony underground level
  pcol_p(2) = 110000.      ! Phony underground level

  do k = 3, nprz+2
     pcol_p(k) = pnpr(k-2) ! Analysis levels
  enddo

  ! If necessary, fill upper pressure levels from ZONAVG data

  do levp = lzon_bot, 22
     k = levp + kzonoff
     pcol_p(k) = zonp_vect(levp)
  enddo

  ! Inverse exner function at analysis levels (defined without cp factor)

  do k = 1, npd
     pcol_exneri(k) = ( p00 / pcol_p(k) )**rocp
     plog       (k) = log( pcol_p(k) )
  enddo

  ! Allocate space for analysis arrays

  allocate( a2d (nprx,npry) )
  allocate( pcol_z(npd,mwa) )
  allocate( field(nprz,mwa) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read geopotential height
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  varname = 'GEO'
  call shdf5_info(varname, ndims, idims, chunk_dims=chunkdims)

  if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
     write(*,*) "Geopotential height (GEO) not found in analysis file."
     stop
  endif

  write(io6,*) "Reading geopotential height " // trim(varname)

  if (chunkdims(ndims) > 1) then
     write(io6,*)
     write(io6,*) "!! NOTE: HDF5 data is stored in chunks that extend over multiple levels."
     write(io6,*) "!! OLAM reads analysis data along a horizontal slice to save memory."
     write(io6,*) "!! Reading is slow with chunks that extend over multiple vetical levels."
     write(io6,*) "!! You can reprocess your HDF5 data with the command:"
     write(io6,*)
     write(io6,*) "h5repack -l CHUNK=1xNYxNX -f SHUF -f GZIP=4 infile outfile"
     write(io6,*)
     write(io6,*) "!! Where NY and NX are the number of longitudes and latitudes in the file."
     write(io6,*) "!! You can also reprocess the data with a more recent version of grib2olam"
     write(io6,*) "!! which stores compressed data approprately for OLAM."
     write(io6,*)
  endif

  do iz = 1, nprz

     call shdf5_irec( ndims, idims, varname, rvar2=a2d, &
                      start=[1,1,iz], counts=[nprx,npry,1] )

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        call gdtost_ll( nprx, npry, glonw(iw), glatw(iw), pcol_z(iz+2,iw), &
                        xswlon, xswlat, gdatdx, gdatdy, inproj, r2d=a2d, &
                        plats=plats, irev_ns=irev_ns )

     enddo
     !$omp end parallel do

  enddo

  ! Special for extrapolated levels at bottom

  dprat2 = 10000. / ((pcol_p(3) - pcol_p(4)) * 1.05)
  dprat1 = 10000. / ((pcol_p(3) - pcol_p(4)) * 1.13)

  !$omp parallel do private(iw,dg,k)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Set phony underground height field
     dg           = pcol_z(4,iw) - pcol_z(3,iw)
     pcol_z(2,iw) = pcol_z(3,iw) - dg * dprat2
     pcol_z(1,iw) = pcol_z(2,iw) - dg * dprat1

     ! Set any extra levels above analysis from climatological zonavg arrays
     if (npd > nprz+2) then
        call hinterp_zonavg(glatw(iw), zonz, npd, pcol_z(:,iw), lzon_bot, kzonoff)
     endif

     ! Interpolate pressure to model levels
     if (.not. log_interp_p) then
        call hintrp_cc(npd, pcol_p, pcol_z(:,iw), mza, o_press(:,iw), zt)
     else
        call hintrp_cc(npd, plog, pcol_z(:,iw), mza, o_press(:,iw), zt)
        do k = 1, mza
           o_press(k,iw) = exp( o_press(k,iw) )
        enddo
     endif

     ! Store analysis 500 mb height (level kbc)
     z_pbc(iw) = pcol_z(kbc+2,iw)

  enddo
  !$omp end parallel do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  varname = 'TEMP'
  call shdf5_info(varname, ndims, idims)

  if (ndims /= 3 .or. any( idims /= [nprx,npry,nprz] ) ) then
     write(*,*) "Air temperature (TEMP) not found in analysis file."
     stop
  endif

  write(io6,*) "Reading air temperature " // trim(varname)

  do iz = 1, nprz

     call shdf5_irec( ndims, idims, varname, rvar2=a2d, &
                      start=[1,1,iz], counts=[nprx,npry,1] )

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        call gdtost_ll( nprx, npry, glonw(iw), glatw(iw), field(iz,iw), &
                        xswlon, xswlat, gdatdx, gdatdy, inproj, r2d=a2d, &
                        plats=plats, irev_ns=irev_ns )
     enddo
     !$omp end parallel do

  enddo

  !$omp parallel private(pvect)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Copy temperature to expanded column array
     do k = 1, nprz
        pvect(k+2) = field(k,iw)
     enddo

     ! Set any extra levels above analysis from climatological zonavg arrays
     if (npd > nprz+2) then
        call hinterp_zonavg(glatw(iw), zont, npd, pvect, lzon_bot, kzonoff)
     endif

     ! Convert column to potential temperature
     do k = 3, npd
        pvect(k) = pvect(k) * pcol_exneri(k)
     enddo

     ! Set theta of phony underground levels
     pvect(1:2) = pvect(3)

     ! Vertically interpolate potential temperature to model levels
     call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_theta(:,iw), zt)

  enddo
  !$omp end do
  !$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read water vapor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! First check if specific humidity is in analysis file
  ndims      = 3
  idims(1:3) = [nprx,npry,nprz]
  vnams(1:3) = ['SHV   ', 'RH    ', 'RELHUM']
  call check_names(3, vnams, ndims, idims, varname)

  if (len_trim(varname) == 0) then
     write(*,*) "Water vapor not found in analysis file."
     stop
  endif

  isrh = .false.
  if (varname(1:1) == 'R') isrh = .true.

  write(io6,*) "Reading water vapor " // trim(varname)

  rh_is_percent = .false.
  sh_is_gkg     = .false.
  nprz_rh       = nprz

  do iz = 1, nprz

     call shdf5_irec( ndims, idims, varname, rvar2=a2d, &
                      start=[1,1,iz], counts=[nprx,npry,1] )

     if (iz == 1) then

        if (isrh) then

           ! RH should be stored as a ratio (0-1). If there are larger values
           ! at the lowest level assume RH is in percent and convert to decimal
           ! RH. This is probably not necessary with recent versions of
           ! grib2olam but we will keep this check anyway.

           if ( any( a2d > 2.0 ) ) then
              rh_is_percent = .true.
              write(io6,*) '    Converting relative humidity ( % ) to mixing ratio'
           else
              write(io6,*) '    Converting relative humidity (0-1) to mixing ratio'
           endif

        else

           ! If any specific humidities at lowest level are greater than 1,
           ! assume it is g/kg and convert to kg/kg. This is probably not
           ! necessary with recent versions of grib2olam but we will keep
           ! this check anyway.

           if ( any( a2d > 1.0 ) ) then
              sh_is_gkg = .true.
              write(io6,*) '    Converting g/kg specific humidity to kg/kg mixing ratio'
           else
              write(io6,*) '    Converting specific humidity to mixing ratio'
           endif

        endif

     endif

     ! Some older reanalyses do not report humidity up to the top of the model.
     ! Check for the highest level that reports humidity (missing values
     ! assumed to be negative here)

     if ( a2d(1,1) < -1.0 ) then
        nprz_rh = iz - 1
        exit
     endif

     if (rh_is_percent) a2d = 0.01 * a2d        ! rh % to fraction
     if (sh_is_gkg)     a2d = 0.001 * a2d       ! g/kg to kg/kg sh
     if (.not. isrh)    a2d = a2d / (1.0 - a2d) ! sh to mixing ratio

     !$omp parallel do private(iw,airtemp,vapor_press)
     do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        ! "field" contains air temperature, which may be needed to convert
        ! RH to mixing ratio. Store in temporary variable if needed.

        if (isrh) airtemp = field(iz,iw)

        call gdtost_ll( nprx, npry, glonw(iw), glatw(iw), field(iz,iw), &
                        xswlon, xswlat, gdatdx, gdatdy, inproj, r2d=a2d, &
                        plats=plats, irev_ns=irev_ns )

        if (isrh) then

           ! Compute ambient vapor pressure based on R.H.
           ! and saturation vapor pressure

           if (airtemp >= t00) then
              vapor_press = field(iz,iw) * eslf( airtemp-t00 )
           else if (airtemp > t00 - 20.) then
              frac = (airtemp - t00 + 20.) / 20.
              vapor_press = field(iz,iw) * ( eslf( airtemp-t00 ) *       frac  &
                                           + esif( airtemp-t00 ) * (1. - frac) )
           else
              vapor_press = field(iz,iw) * esif( airtemp-t00 )
           endif

           ! Do not allow vapor pressure to exceed ambient pressure

           vapor_press = min( 0.9*pnpr(iz), vapor_press )

           ! Compute mixing ratio from vapor press and ambient press

           field(iz,iw) = eps_vap * vapor_press / ( pnpr(iz) - vapor_press )

        endif

     enddo
     !$omp end parallel do

  enddo

  if (nprz_rh < nprz) then
     write(io6,'(A,I0)'  ) '  Humidity is only reported up to level ', nprz_rh
     write(io6,'(A,I0,A)') '  out of ', nprz, ' in the analysis file. Missing'
     write(io6,'(A)'     ) '  humidity will be set from the Mclatchy soundings.'
  endif

  !$omp parallel private(pvect)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Copy humidity to expanded column array
     do k = 1, nprz_rh
        pvect(k+2) = field(k,iw)
     enddo

     ! Fill any missing humidity levels from the Mclatchy soundings
     if (nprz_rh < nprz) then
        call fill_wvap_plevs_mclat(glatw(iw), nprz_rh+3, nprz+2, npd, pcol_p, pvect)
     endif

     ! Set any extra levels above analysis from climatological zonavg arrays
     if (npd > nprz+2) then
        call hinterp_zonavg(glatw(iw), zonr, npd, pvect, lzon_bot, kzonoff)
     endif

     ! Phony underground levels
     pvect(1:2) = pvect(3)

     ! Set lower limit on water vapor mixing ratio
     pvect = max(pvect, 1.e-8)

     ! If logarithmically interpolating mixing ratio
     if (log_interp_q) then
        do k = 1, npd
           pvect(k) = log(pvect(k))
        enddo
     endif

     ! Vertically interpolate humidity from pressure levels to model levels
     call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_rrw(:,iw), zt)

     ! If logarithmically interpolating mixing ratio
     if (log_interp_q) then
        do k = 1, mza
           o_rrw(k,iw) = exp( o_rrw(k,iw) )
        enddo
     endif

     ! Now that we have water vapor, pressure, and theta, compute density
     if (miclevel == 0) then
        do k = 1, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / o_theta(k,iw)
        enddo
     else
        do k = 1, mza
           o_rho(k,iw) = o_press(k,iw)**cvocp * p00kord / &
                         ( o_theta(k,iw) * (1.0 + eps_vapi * o_rrw(k,iw)) )
        enddo
     endif

  enddo
  !$iomp end do
  !$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read zonal wind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! East-west velocity may be called UP, UE, or U
  ndims      = 3
  idims(1:3) = [nprx,npry,nprz]
  vnams(1:3) = ['U ', 'UP', 'UE']
  call check_names(3, vnams, ndims, idims, varname)

  if (len_trim(varname) == 0) then
     write(*,*) "Zonal wind (U) not found in analysis file."
     stop
  endif

  write(io6,*) "Reading zonal wind " // trim(varname)

  do iz = 1, nprz

     call shdf5_irec( ndims, idims, varname, rvar2=a2d, &
                      start=[1,1,iz], counts=[nprx,npry,1] )

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        call gdtost_ll( nprx, npry, glonw(iw), glatw(iw), field(iz,iw), &
                        xswlon, xswlat, gdatdx, gdatdy, inproj, r2d=a2d, &
                        plats=plats, irev_ns=irev_ns )
     enddo
     !$omp end parallel do

  enddo

  !$omp parallel private(pvect)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Copy zonal wind to expanded column array
     do k = 1, nprz
        pvect(k+2) = field(k,iw)
     enddo

     ! Set any extra levels above analysis from climatological zonavg arrays
     if (npd > nprz+2) then
        call hinterp_zonavg(glatw(iw), zonu, npd, pvect, lzon_bot, kzonoff)
     endif

     ! Phony underground levels
     pvect(1:2) = pvect(3)

     ! Vertically interpolate wind from pressure levels to model levels
     call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_uzonal(:,iw), zt)

  enddo
  !$omp end do
  !$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read meridional wind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! North-south velocity may be called VP, VE, or V
  ndims      = 3
  idims(1:3) = [nprx,npry,nprz]
  vnams(1:3) = ['V ', 'VP', 'VE']
  call check_names(3, vnams, ndims, idims, varname)

  if (len_trim(varname) == 0) then
     write(*,*) "Meridional wind (V) not found in analysis file."
     stop
  endif

  write(io6,*) "Reading meridional wind " // trim(varname)

  do iz = 1, nprz

     call shdf5_irec( ndims, idims, varname, rvar2=a2d, &
                      start=[1,1,iz], counts=[nprx,npry,1] )

     !$omp parallel do private(iw)
     do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        call gdtost_ll( nprx, npry, glonw(iw), glatw(iw), field(iz,iw), &
                        xswlon, xswlat, gdatdx, gdatdy, inproj, r2d=a2d, &
                        plats=plats, irev_ns=irev_ns )
     enddo
     !$omp end parallel do

  enddo

  !$omp parallel private(pvect)
  !$omp do private(iw,k)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Copy meridional wind to expanded column array
     do k = 1, nprz
        pvect(k+2) = field(k,iw)
     enddo

     ! Set any extra levels above analysis
     if (npd > nprz+2) pvect(nprz+3:npd) = 0.0

     ! Phony underground levels
     pvect(1:2) = pvect(3)

     ! Vertically interpolate wind from pressure levels to model levels
     call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_umerid(:,iw), zt)

  enddo
  !$omp end do
  !$omp end parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read ozone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (i_o3 > 0) then

     nbot_o3 = nprz + 1

     ! Check for ozone variable in file

     ndims      = 3
     idims(1:3) = [nprx,npry,nprz]
     vnams(1:2) = ['O3MR ', 'OZONE']
     call check_names(2, vnams, ndims, idims, varname)

     if (len_trim(varname) > 0) then

        write(io6,*) "Reading ozone mixing ratio " // trim(varname)

        do iz = 1, nprz

           call shdf5_irec( ndims, idims, varname, rvar2=a2d, &
                            start=[1,1,iz], counts=[nprx,npry,1] )

           ! Special for OZONE:
           ! Some analysis only report ozone ABOVE a certain level. Check for
           ! the lowest level at which ozone is reported in the analysis file:

           if (a2d(1,1) < -1.0) then

              nbot_o3 = iz + 1

           else

              !$omp parallel do private(iw)
              do j = 1, jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

                 call gdtost_ll( nprx, npry, glonw(iw), glatw(iw), field(iz,iw), &
                                 xswlon, xswlat, gdatdx, gdatdy, inproj, r2d=a2d, &
                                 plats=plats, irev_ns=irev_ns )
              enddo
              !$omp end parallel do

           endif

        enddo

     endif

     !$omp parallel private(pvect)
     !$omp do private(iw,k)
     do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        ! Copy ozone to expanded column array
        do k = nbot_o3, nprz
           pvect(k+2) = field(k,iw)
        enddo

        ! Fill any missing ozone levels from Mclatchy soundings
        if ( nbot_o3 > 1 ) then
           call fill_ozone_plevs_mclat(glatw(iw), 3, nbot_o3+1, npd, pcol_p, pvect)
        endif

        ! Set any extra levels above analysis from climatological zonavg arrays
        if (npd > nprz+2) then
           call hinterp_zonavg(glatw(iw), zono, npd, pvect, lzon_bot, kzonoff)
        endif

        ! Phony underground levels
        pvect(1:2) = pvect(3)

        ! Vertically interpolate ozone from pressure levels to model levels
        call hintrp_cc(npd, pvect, pcol_z(:,iw), mza, o_ozone(:,iw), zt)

        ! Convert to ppmV
        do k = 1, mza
           o_ozone(k,iw) = max( 1.e-30, cnvto3 * o_ozone(k,iw) )
        enddo

     enddo
     !$omp end do
     !$omp end parallel

  endif

end subroutine pressure_stage



subroutine check_names(nv, vnames, ndims, idims, varname)

  use hdf5_utils, only: shdf5_info
  implicit none

  integer,      intent(in)  :: nv
  character(*), intent(in)  :: vnames(nv)
  integer,      intent(in)  :: ndims
  integer,      intent(in)  :: idims(ndims)
  character(*), intent(out) :: varname

  integer :: i, ndims2, idims2(4)

  varname = ' '

  if (ndims==0) then
     write(*,*) "Error in check_names:"
     write(*,*) "     Set ndims to expected array size in routine check_names"
     return
  endif

  do i = 1, nv
     call shdf5_info(vnames(i), ndims2, idims2)

     if (ndims2 == ndims .and. all(idims2(1:ndims) == idims(1:ndims))) then
        varname = vnames(i)
        exit
     endif
  enddo

end subroutine check_names




subroutine hinterp_zonavg(glat, zona, npd, field, lzon_bot, kzonoff)

  use mem_zonavg, only: nlata, nplev
  implicit none

  integer, intent(in)    :: npd, lzon_bot, kzonoff
  real,    intent(in)    :: glat, zona(nlata,nplev)
  real,    intent(inout) :: field(npd)

  integer :: ilat, k, levp
  real    :: rlat, wt2

  rlat = .4 * (glat + 93.75)
  ilat = int(rlat)
  wt2  = rlat - real(ilat)

  do levp = lzon_bot, 22
     k = levp + kzonoff
     field(k) = (1. - wt2) * zona(ilat,levp) + wt2 * zona(ilat+1,levp)
  enddo

end subroutine hinterp_zonavg




subroutine fill_wvap_plevs_mclat(glat, nbot, ntop, nlevs, pcol, qvect)

  use mem_mclat, only: sslat, mclat, ypp_mclat, mclat_spline
  implicit none

  real,    intent(in)  :: glat
  integer, intent(in)  :: nbot, ntop, nlevs
  real,    intent(in)  :: pcol (nlevs)
  real,    intent(out) :: qvect(nlevs)

  real    :: mcol(33,6)
  integer :: lv, nmax

  ! Fills pressure levels from nbot to ntop with water vapor values
  ! from the Mclatchy soundings

  nmax = min(ntop,nlevs)

  if (nbot > nmax) return

  ! This assumes that mclat_spline has already been called
  ! with the current date

  call spline2_vec(13, 33*6, sslat, mclat, ypp_mclat, glat, mcol)

  ! Compute water vapor mixing ratios by dividing by dry density

  do lv = 1,33
     mcol(lv,4) = mcol(lv,4) / (mcol(lv,6) - mcol(lv,4))
  enddo

  ! Vertically interpolate Mclatchy water vapor BY PRESSURE to analysis levels

  lv = nmax - nbot + 1

  call pintrp_ee(33, mcol(:,4), mcol(:,2), lv, qvect(nbot:nmax), pcol(nbot:nmax))

end subroutine fill_wvap_plevs_mclat




subroutine fill_ozone_plevs_mclat(glat, nbot, ntop, nlevs, pcol, ovect)

  use mem_mclat, only: sslat, mclat, ypp_mclat, mclat_spline
  implicit none

  real,    intent(in)  :: glat
  integer, intent(in)  :: nbot, ntop, nlevs
  real,    intent(in)  :: pcol (nlevs)
  real,    intent(out) :: ovect(nlevs)

  real    :: mcol(33,6)
  integer :: lv, nmax

  ! Fills pressure levels from nbot to ntop with ozone values
  ! from the Mclatchy soundings

  nmax = min(ntop,nlevs)

  if (nmax > nmax) return

  ! This assumes that mclat_spline has already been called
  ! with the current date

  call spline2_vec(13, 33*6, sslat, mclat, ypp_mclat, glat, mcol)

  ! Compute ozone mixing ratios by dividing by dry density

  do lv = 1,33
     mcol(lv,5) = mcol(lv,5) / (mcol(lv,6) - mcol(lv,4))
  enddo

  ! Vertically interpolate Mclatchy ozone BY PRESSURE to analysis levels

  lv = nmax - nbot + 1

  call pintrp_ee(33, mcol(:,5), mcol(:,2), lv, ovect(nbot:nmax), pcol(nbot:nmax))

end subroutine fill_ozone_plevs_mclat
