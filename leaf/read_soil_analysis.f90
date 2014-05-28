subroutine read_soil_analysis(soil_tempc)

  use misc_coms,  only: io6, s1900_sim, s1900_init, isubdomain
  use leaf_coms,  only: nzg, mwl, slzt, soilcp, slmsts, slcpd, soilstate_db
  use mem_leaf,   only: land, itab_wl
  use consts_coms,only: pio180, piu180, erad, cliq1000, alli1000, cice,   &
                         cice1000, r8
  use max_dims,   only: pathlen
  use isan_coms,  only: nfgfiles, s1900_fg, fnames_fg
  use hdf5_utils, only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close
  use mem_para,   only: myrank

  implicit none

  real, intent(inout) :: soil_tempc(nzg,mwl)

  integer            :: ifgfile, nf, ipoffset
  character(pathlen) :: fname
  character(16)      :: ext
  logical            :: exists, isnow, itemp, iswat
  integer            :: ndims, idims(3)
  integer            :: nx, ny, ngnd, inproj
  real               :: xswlat, xswlon, gdatdx, gdatdy
  real               :: xoffpix, yoffpix, xperdeg, yperdeg
  real               :: glat, glon, rio, rjo, dss
  real               :: snowdens
  integer            :: iwl, io, jo, i, j, k, kk, ntext
  logical            :: has_snow, has_soilt, has_soilw, soil_volw

  real, allocatable  :: snow (:,:)   ! snow mass  [kg/m2]
  real, allocatable  :: soilt(:,:,:) ! soil temp  [K]
  real, allocatable  :: soilw(:,:,:) ! soil water [frac of sat.]

  real, allocatable  :: a2d(:,:), a3d(:,:,:), tcol(:), wcol(:)
  real               :: wgt0(-1:1,-1:1), wgts(-1:1,-1:1)
  real, allocatable  :: snowmask(:,:), tempmask(:,:), soilwmask(:,:)
  real, allocatable  :: zcol(:), ztmp(:)

  has_snow  = .false.
  has_soilt = .false.
  has_soilw = .false.
  soil_volw = .false. ! Is volumetric soil water input or fraction of saturation
  ngnd      =  0

! Loop over analysis files and search for the one that corresponds
! to the model start time

  ifgfile = 0
  do nf = 1, nfgfiles
     if (s1900_fg(nf) <= s1900_sim) then
        ifgfile = nf
     endif
  enddo
  
  if (ifgfile < 1) then
     write(io6,*) ' '
     write(io6,*) 'Unable to find analysis file for soil initialization'
     write(io6,*) 'Using default initialization'
     return
  elseif (s1900_fg(ifgfile) < s1900_init - 1800.0_r8) then
     write(io6,*) ' '
     write(io6,*) 'Available analysis file is more than 1800 seconds'
     write(io6,*) 'prior to simulation initialization time.'
     write(io6,*) 'Ising default soil initialization'
     return
  endif

  fname = fnames_fg(ifgfile)

! Make sure we have a HDF5 analysis file

  ext = trim( fname( index(fname,'.',back=.true.)+1:) )

  if ( .not. ( &
       ext == 'h5' .or. ext == 'hdf5' .or. ext == 'hdf' .or. &
       ext == 'H5' .or. ext == 'HDF5' .or. ext == 'HDF' )) then
     write(io6,*) "read_soil: Analysis file is not hdf5 format"
     write(io6,*) "Using default soil initialization instead."
     return
  endif

! Process selected analysis file

  inquire(file=fname, exist=exists)
  if (.not. exists) then
     write(io6,*) "read_soil: Error opening analysis file " // trim(fname)
     write(io6,*) "Using default soil initialization instead."
     return
  endif

  call shdf5_open (fname, 'R')

  ndims    = 1
  idims(1) = 1
  idims(2) = 1
  idims(3) = 1

  call shdf5_irec(ndims, idims, 'nx'   , ivars=nx)
  call shdf5_irec(ndims, idims, 'ny'   , ivars=ny)
  call shdf5_irec(ndims, idims, 'iproj', ivars=inproj)
  call shdf5_irec(ndims, idims, 'swlat', rvars=xswlat)
  call shdf5_irec(ndims, idims, 'swlon', rvars=xswlon)
  call shdf5_irec(ndims, idims, 'dx'   , rvars=gdatdx)
  call shdf5_irec(ndims, idims, 'dy'   , rvars=gdatdy)

  ! Check if ngnd, the # of soil levels, is in the analysis file and read it

  call shdf5_info('ngnd', ndims, idims)
  if (ndims > 0) call shdf5_irec(ndims, idims, 'ngnd' , ivars=ngnd)

  ! Check if sdepths, the soil depth array, is in the analysis file and read it
  
  call shdf5_info('sdepths', ndims, idims)
  if (ndims > 0) then

     if (ngnd == 0) ngnd = idims(1)
     allocate( ztmp(ngnd), zcol(ngnd) )
     call shdf5_irec(ndims, idims, 'sdepths', rvara=ztmp)

     ! olam stores the soil arrays from bottom to top, so we need to reverse
     ! the input soil depth array, and convert to m

     do k = 1, ngnd
        kk = ngnd - k + 1
        zcol(kk) = ztmp(k) * 0.01
     enddo

     deallocate(ztmp)
  else
     write(io6,*) "Gridded analysis dataset does not contain soil depth information."
     write(io6,*) "Soil initialization from analysis file will be skipped."
     write(io6,*)
  endif

  ! Check data domain size and location

  if (inproj /= 1) then
     write(io6,*) 'You must input a lat-lon grid for soil initialization.'
     write(io6,*) "Using default soil initialization instead."
     call shdf5_close()
     return
  endif

  ! We make the requirement that a full global domain of data be
  ! read in.  Check this here.  Following the convention for the NCEP/DOE
  ! Reanalysis2 data, assume that data exists at both latitudinal boundaries
  ! (-90. and 90. degrees) but that the longitudinal boundary is not repeated.
  ! If either is not the case, this check will stop execution. 

  if (abs(nx      * gdatdx - 360.) > .1 .or. &
      abs ((ny-1) * gdatdy - 180.) > .1) then
     write(io6,*) 'Gridded soil data does not have global coverage.'
     write(io6,*) "Using default soil initialization instead."
     call shdf5_close()
     return
  endif

  ! Compute longitudinal offset index, which is the index in the expanded
  ! arrays where the first input data point (at xswlon) is located. We
  ! follow the GRIB standard and keep the SW corner at (0,-90)

  ipoffset = xswlon / gdatdx + 1

  ! Check if snow mass is in the analysis file, and read it

  call shdf5_info('SNOWMASS', ndims, idims)

  if (ndims > 0) then

     allocate(a2d(nx,ny))
     allocate(snow(nx+3,ny+2))

     call shdf5_irec(ndims, idims, 'SNOWMASS', rvara = a2d)
     call prfill(nx, ny, ipoffset, a2d, snow)

     has_snow = .true.
     deallocate(a2d)
  endif

  ! Check if soil temperature and/or soil moisture are in the analysis file, 
  ! and read them

  if (ngnd > 0 .and. allocated(zcol)) then
     allocate(a3d(nx,ny,ngnd))
     
     call shdf5_info('SOILT', ndims, idims)

     if (ndims > 0) then
        allocate(soilt(nx+3,ny+2,ngnd))

        call shdf5_irec(ndims, idims, 'SOILT', rvara = a3d)
        call prfill3(nx, ny, ngnd, ipoffset, a3d, soilt)
        has_soilt = .true.
     endif

     call shdf5_info('SOILW', ndims, idims)

     if (ndims > 0) then
        allocate(soilw(nx+3,ny+2,ngnd))

        call shdf5_irec(ndims, idims, 'SOILW', rvara = a3d)
        call prfill3(nx, ny, ngnd, ipoffset, a3d, soilw)
        has_soilw = .true.

     else

        call shdf5_info('SOILVOLW', ndims, idims)
        
        if (ndims > 0) then
           allocate(soilw(nx+3,ny+2,ngnd))

           call shdf5_irec(ndims, idims, 'SOILVOLW', rvara = a3d)
           call prfill3(nx, ny, ngnd, ipoffset, a3d, soilw)
           has_soilw = .true.
           soil_volw = .true.
        endif
     endif

     deallocate(a3d)
  endif

  call shdf5_close()

  if (.not. has_snow) then
     write(io6,*) "read_soil: Analysis file does not contain snow mass."
     write(io6,*) "Skipping snow depth initialization."
     write(io6,*)
  endif

  if (.not. has_soilt) then
     write(io6,*) "read_soil: Analysis file does not contain soil temperature."
     write(io6,*) "Using default soil temperatures instead."
     write(io6,*)
  endif
  
  if (.not. has_soilw) then
     write(io6,*) "read_soil: Analysis file does not contain soil moisture."
     write(io6,*) "Using default soil water instead."
     write(io6,*)
  endif

  ! Make sure that the input fields don't contain all missing data. This
  ! can happen if we specified to save the soil fields during degribbing,
  ! but the grib file didn't have any soil or snow fields

  if (has_snow) then
     if ( count( snow(:,:) < -1.0 .or. snow(:,:) > 1.e30) >= nx*ny ) then
        write(io6,*) "read_soil: Snow mass from analysis contains mostly missing data."
        write(io6,*) "Skipping snow depth initialization."
        has_snow = .false.
     endif
  endif
  
  if (has_soilt) then
     if ( count( soilt(:,:,:) < -1.0  .or. soilt(:,:,:) > 1.e30) >= nx*ny*ngnd ) then
        write(io6,*) "read_soil: Soil temperature from analysis contains mostly missing data."
        write(io6,*) "Skipping soil temperature initialization."
        has_soilt = .false.
     endif
  endif
  
  if (has_soilw) then
     if ( count( soilw(:,:,:) < -1.0 .or. soilw(:,:,:) > 1.e30) >= nx*ny*ngnd ) then
        write(io6,*) "read_soil: Soil moisture from analysis contains mostly missing data."
        write(io6,*) "Skipping soil moisture initialization."
        has_soilw = .false.
     endif
  endif

  ! If all data missing, return
  
  if ((.not. has_snow) .and. (.not. has_soilt) .and. (.not. has_soilw)) return

  ! Set masks to indicate missing data

  if (has_snow) then
     allocate(snowmask(nx+3,ny+2))
     
     do i = 1, nx+3
        do j = 1, ny+2
           if (snow(i,j) < -1.0 .or. snow(i,j) > 1.e30) then
              snowmask(i,j) = 0.0
           else
              snowmask(i,j) = 1.0
              snow    (i,j) = max( snow(i,j), 0.0 )
           endif
        enddo
     enddo
  endif

  if (has_soilt) then
     allocate(tempmask(nx+3,ny+2))

     do i = 1, nx+3
        do j = 1, ny+2
           if ( any( soilt(i,j,:) < 100.0 .or. soilt(i,j,:) > 1.e30 )) then
              tempmask(i,j) = 0.0
           else
              tempmask(i,j)   = 1.0
              soilt   (i,j,:) = max(soilt(i,j,:), 200.0)
           endif
        enddo
     enddo
  endif

  if (has_soilw) then
     allocate(soilwmask(nx+3,ny+2))

     do i = 1, nx+3
        do j = 1, ny+2
           if ( any( soilw(i,j,:) < -1.0 .or. soilw(i,j,:) > 1.e30 )) then
              soilwmask(i,j) = 0.0
           else
              soilwmask(i,j)   = 1.0
              soilw    (i,j,:) = max(soilw(i,j,:), 0.0)
           endif
        enddo
     enddo
  endif

  allocate(tcol(ngnd))
  allocate(wcol(ngnd))

  ! GRIB data is unstaggered, so there should be no offsets, but we 
  ! offset the grid by 1 row and column so we take that into account here

  xoffpix = 1.0
  yoffpix = 1.0

  xperdeg = 1.0 / gdatdx
  yperdeg = 1.0 / gdatdy

! Fill seaice array

do iwl = 2, mwl

   ! Skip this cell if running in parallel and primary rank of IWL /= MYRANK

   if (isubdomain == 1 .and. itab_wl(iwl)%irank /= myrank) cycle

   glat = land%glatw(iwl)
   glon = land%glonw(iwl)

   ! Convert OLAM's -180,180 longitude range to GRIB's 0,360 

   if (glon < 0.0) glon = 360.0 + glon 
   glon = max(0.001,min(359.999,glon))

   ! Find the nearest analysis point to this land cell

   rio = 1. + (glon      ) * xperdeg + xoffpix
   rjo = 1. + (glat + 90.) * yperdeg + yoffpix

   io = nint(rio)
   jo = nint(rjo)

   ! Weight the 9 nearest analysis points by their inverse distance squared

   do i = -1, 1
      do j = -1, 1
         dss = (rio - real(io+i))**2 + (rjo - real(jo+j))**2
         wgt0(i,j) = 1.0 / max(dss,0.05)
      enddo
   enddo

   ! Interpolate snow depth to this land cell if any of the 9 closest
   ! analysis points have non-missing snow data

   if (has_snow .and. any( snowmask(io-1:io+1,jo-1:jo+1) > 0.5 )) then
      
      wgts = wgt0 * snowmask(io-1:io+1,jo-1:jo+1)
      wgts = wgts / sum(wgts)

      land%sfcwater_mass  (1,iwl) = sum(wgts * max(snow(io-1:io+1,jo-1:jo+1), 0.))

      land%sfcwater_energy(1,iwl) = min(0., (land%can_temp(iwl) - 273.15) * cice)

      ! snow density calculation comes from CLM3.0 documentation 
      ! which is based on Anderson 1975 NWS Technical Doc # 19 

      snowdens = 50.0
      if (land%can_temp(iwl) > 258.15) snowdens =   &
           50.0 + 1.5 * (land%can_temp(iwl) - 258.15)**1.5

      land%sfcwater_depth(1,iwl) = land%sfcwater_mass(1,iwl) / snowdens

   endif

   ! Interpolate soil temperature to this land cell if any of the 9 closest
   ! analysis points have non-missing temperature data

   if (has_soilt .and. any( tempmask(io-1:io+1,jo-1:jo+1) > 0.5 )) then
      
      wgts = wgt0 * tempmask(io-1:io+1,jo-1:jo+1)
      wgts = wgts / sum(wgts)

      ! olam has soil depths starting from the bottom

      do k = 1, ngnd
         kk = ngnd - k + 1
         tcol(kk) = sum(wgts * soilt(io-1:io+1,jo-1:jo+1,k)) - 273.15
      enddo
         
      if (ngnd == 1) then
         soil_tempc(1:nzg,iwl) = tcol(1)
      else
         call hintrp_cc( ngnd, tcol, zcol, nzg, soil_tempc(:,iwl), slzt )
      endif

   endif

   ! Interpolate soil moisture to this land cell if any of the 9 closest
   ! analysis points have non-missing soil moisture data

   if (has_soilw .and. any( soilwmask(io-1:io+1,jo-1:jo+1) > 0.5 )) then
      
      wgts = wgt0 * soilwmask(io-1:io+1,jo-1:jo+1)
      wgts = wgts / sum(wgts)

      ! olam has soil depths starting from the bottom

      do k = 1, ngnd
         kk = ngnd - k + 1
         wcol(kk) = sum(wgts * soilw(io-1:io+1,jo-1:jo+1,k))
      enddo
         
      if (ngnd == 1) then
         land%soil_water(1:nzg,iwl) = wcol(1)
      else
         call hintrp_cc( ngnd, wcol, zcol, nzg, land%soil_water(:,iwl), slzt )
      endif

      if (soil_volw) then

         ! Bound soil moisture between capcity soilcp and porosity slmsts

         do k = 1, nzg
            ntext = land%ntext_soil(k,iwl)
            land%soil_water(k,iwl) = max( soilcp(ntext), min( land%soil_water(k,iwl), slmsts(ntext) ))
         enddo

      else

         ! Convert from fraction of saturation to soil moisture [m^3_wat/m^3_tot]
         do k = 1, nzg
            ntext = land%ntext_soil(k,iwl)
            land%soil_water(k,iwl) = max( soilcp(ntext), land%soil_water(k,iwl)*slmsts(ntext) )
         enddo

      endif

   endif

enddo
   
end subroutine read_soil_analysis
