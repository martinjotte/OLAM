module hcane_rz

  use consts_coms, only: r8

  implicit none

! The purpose of module hcane_rz is to restore the intensity and eyewall
! diameter of a hurricane that is under-resolved in the initial conditions of
! an OLAM simulation.  This situation commonly arises when the only gridded
! initialization dataset is of low resolution (for a hurricane), such as the
! CFSR global gridded reanalysis.  The procedure consists of one or more
! dynamic initialization cycles in which (1) OLAM is first initialized by
! interpolation from the CFS (or other) gridded dataset, (2) OLAM is integrated
! forward in time one or more hours with a user-specified heat source to allow
! the inner region of the hurricane to intensify and for the hurricane as a whole
! to develop convective cells and mature condensate fields, and (3) the simulated
! hurricane is remapped from its final simulated location back to the grid cells
! at its initial location in preparation for the next cycle.
!
! The number of cycles performed, the temporal length of each cycle, and the
! location and intensity of the heat input are all user-specified parameters.
! Plots of the vortex are made at the beginning and end of each cycle to allow
! inspection of the progress and quality of the initialization.

  ! radius_ax is an array of radial distances of the radial-height grid on which
  ! axisymmetric hurricane fields are diagnosed and plotted.  The radial spacing
  ! between these values needs to be at least twice the grid spacing of the OLAM
  ! hexagonal grid where the hurricane is located.

  integer, parameter :: nr = 28 ! Number of radial points in vortex axial averages

  !real :: radius_ax(nr) = (/ &
  !   0.e3,  10.e3,  20.e3,  30.e3,  40.e3,  50.e3,  60.e3,  70.e3,  80.e3,  90.e3, &
  ! 100.e3, 120.e3, 140.e3, 160.e3, 180.e3, 200.e3, 250.e3, 300.e3, 350.e3, 400.e3 /)

  real :: radius_ax(nr) = (/ &
     0.e3,   5.e3,  10.e3,  15.e3,  20.e3,  25.e3,  30.e3,  35.e3,  40.e3,  45.e3, &
    50.e3,  55.e3,  60.e3,  70.e3,  80.e3,  90.e3, 100.e3, 120.e3, 140.e3, 160.e3, &
   180.e3, 200.e3, 240.e3, 280.e3, 320.e3, 360.e3, 400.e3, 440.e3 /)

  integer :: ncycle_hurrinit ! Number of initialization cycles to perform
  integer :: icycle_hurrinit ! Ending integration time of each cycle (s)

  real(r8) :: timmax_hurrinit

  real :: hlat0, hlon0   ! Initial obs hurricane lat/lon (deg)

  real :: rad1_blend     ! Inner radius for relocation blending weights (m)
  real :: rad2_blend     ! Outer radius for relocation blending weights (m)

  real :: zcent_thpert   ! Center height of toroidal heating region (m)
  real :: zhwid_thpert   ! Vertical half-width of toroidal heating region (m)

  real :: rcent_thpert   ! Center radius of toroidal heating region (m)
  real :: rhwid_thpert   ! Radial half-width of toroidal heating region (m)

  real :: maxrate_thpert ! Maximum heating rate in toroidal heating region (K/s)

  real :: vtan_targ      ! Target max tangential wind speed (m/s)
  real :: vtan_max       ! Current max tangential wind speed (m/s)

  real :: pmsl_targ      ! Target hurricane central MSL pressure (Pa)
  real :: pmsl_min       ! Current min hurricane MSL pressure (Pa)

  real :: hlat_hist = 0., hlon_hist = 0. ! Stored on history file
  real :: hlat_reloc, hlon_reloc         ! Relocation point
  real :: hlat, hlon                     ! Input to and updated by vortex_center_diagnose;
                                         ! used by vortex_azim_avg and vortex_add_thetapert

  real :: xeh_reloc, yeh_reloc, zeh_reloc
  real :: cos_hlat_reloc, sin_hlat_reloc
  real :: cos_hlon_reloc, sin_hlon_reloc

  ! Time series of hlat, hlon, model time

  integer :: nhtim  ! Number of model timesteps

  real, allocatable :: hlata(:,:)
  real, allocatable :: hlona(:,:)
  real, allocatable :: htima(:,:)

  character(128) :: htc0 = 'n'  ! pathlen = 128 in max_dims.f90

  ! rad1_blend and rad2_blend control how a relocated vortex is blended in with
  ! model initial fields at the beginning of the second and subsequent dynamic
  ! initialization cycles.  At locations inside radius rad1, 100% of the relocated
  ! fields are used (implying 0% of the model initial fields).  At locations
  ! outside radius rad2, none of the relocated vortex is used, and the model
  ! initial fields remain unchanged.  Between these limits, interpolation weights
  ! vary linearly with radius.

  ! Relocation fields for I/O

  integer, parameter :: nfld = 31    ! number of fields relocated

  integer :: nout_dim ! Horiz # of W points filled by interpolation (array dim)
  integer :: nout     ! Horiz # of W points filled by interpolation (actual count)

  integer, parameter :: nips = 200, njps = 200 ! size of PS arrays

  integer :: nwps(nips,njps)    ! # of iw indices in each ips,jps cell
  integer :: iwps(nips,njps,60) ! iw indices in each ips,jps cell (max of 60)

  integer, allocatable :: iwout(:)

  real, allocatable :: reloc_field(:,:,:)

  integer :: kthpert_min = 2
  integer :: kthpert_max = 1
  integer :: l8k

  real, allocatable :: wt_vert (:)
  real, allocatable :: dradi_ax(:)

  real :: rmax, rmaxsq

Contains

!===============================================================================

  subroutine hurricane_init()

  use misc_coms,  only: timmax8, dtlm
  use mem_grid,   only: mza, zt

  implicit none

  integer         :: nhtim, nhcyc, k
  real            :: delz
  real, parameter :: zexpon_thpert = 2.0    ! []

  nhtim = max(timmax8, timmax_hurrinit) / dtlm
  nhcyc = max(1,ncycle_hurrinit)

  rmax   = 0.999999 * radius_ax(nr)
  rmaxsq = rmax * rmax

  allocate (hlata(0:nhtim,nhcyc))
  allocate (hlona(0:nhtim,nhcyc))
  allocate (htima(0:nhtim,nhcyc)); htima(:,:) = 0.

  do k = 2, mza
     if (zt(k) > zcent_thpert - zhwid_thpert) exit
  enddo
  kthpert_min = k

  do k = kthpert_min, mza
     if (zt(k) > zcent_thpert + zhwid_thpert) exit
  enddo
  kthpert_max = k-1

  allocate(wt_vert(mza)) ; wt_vert = 0.

  do k = kthpert_min, kthpert_max
     delz       = abs( zt(k) - zcent_thpert )
     wt_vert(k) = 1. - (delz / zhwid_thpert) ** zexpon_thpert
  enddo

  do k = 2, mza
     if (zt(k) > 8.e3) exit
  enddo
  l8k = k - 1

  allocate(dradi_ax(nr))
  do k = 1, nr-1
     dradi_ax(k) = 1.0 / (radius_ax(k+1) - radius_ax(k))
  enddo
  dradi_ax(nr) = dradi_ax(nr-1)

  end subroutine hurricane_init

!==================================================================================

  subroutine vortex_center_diagnose()

  ! This subroutine diagnoses the latitude & longitude of the MSL pressure minimum
  ! inside a tropical cyclone.  It requires a nearby starting location (hlat,hlon)
  ! to be provided, either from a user specification of the known cyclone location
  ! or from a recent location previously diagnosed from this subroutine.

  ! The algorithm of subroutine vortex_center_diagnose is not to simply search for
  ! the lowest MSL pressure in the vicinity, but rather to take a weighted average
  ! over multiple grid cells having pressure close to that lowest value.  This yields
  ! a location that is more steady in time, less subject to turbulent fluctuations,
  ! and more tied to the circulation at the scale of the cyclone core region.

  ! This works now for parallel run

  use mem_ijtabs,  only: jtab_w, jtw_prog
  use mem_basic,   only: press, tair
  use misc_coms,   only: mstp, time8, iparallel, io6
  use mem_grid,    only: mwa, zt, lpw, xew, yew, zew, arw0, glatw, glonw
  use consts_coms, only: pio180, erad

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer :: j, iw, ka
  logical :: iskip(mwa)

  real :: reh, xeh, yeh, zeh
  real :: distsq, weight
  real :: rlon, rlat, dlat, dlon

  real :: area_tot, pmsl_avg, pmsl_thresh
  real :: xew_avg, yew_avg, zew_avg, weight_sum
  real :: pmsl(mwa)

#ifdef OLAM_MPI
  real    :: buffer(4)
  integer :: ierr
#endif

  write(io6,'(a,2f10.3)') 'vortex_center_diagnose BEGIN ',hlat,hlon

  ! Find "earth" coordinates of previous vortex center location (hlat,hlon)

  reh = erad * cos(hlat * pio180)  ! distance from earth center
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)
  zeh = erad * sin(hlat * pio180)

  ! Initialize quantities

  area_tot = 0.
  pmsl_min = 2.e5
  pmsl_avg = 0.

  xew_avg = 0.
  yew_avg = 0.
  zew_avg = 0.

  weight_sum = 0.

  iskip = .true.

  ! Horizontal loop over W points
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Skip current W point if its location is far from first guess position
     ! (This is rough check that eliminates most points)

     dlat = abs(glatw(iw) - hlat)
     dlon = abs(glonw(iw) - hlon)

     if (dlat > 5. .and. dlat < 355.) cycle
     if (dlon > 5. .and. dlon < 355.) cycle

     ! Distance of current W point to first guess position

     distsq = (xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2

     ! Skip current W point if its location is far from first guess position
     ! (This is finer check)

     if (distsq > 1.e10) cycle

     iskip(iw) = .false.

     ka = lpw(iw)

     pmsl(iw) = press(ka,iw) &
              * (1. - .0065 * zt(ka) / (tair(ka,iw) + .0065 * zt(ka)))**(-5.257)

     ! Determine lowest value of pmsl within search area

     if (pmsl_min > pmsl(iw)) then
         pmsl_min = pmsl(iw)
     endif

     ! Sum grid cell area and (area * pressure) product at each model level within search area

     area_tot = area_tot + arw0(iw)
     pmsl_avg = pmsl_avg + arw0(iw) * pmsl(iw)

  enddo  ! iw

#ifdef OLAM_MPI
  if (iparallel == 1) then
     buffer(1:2) = [area_tot, pmsl_avg]
     call MPI_Allreduce( MPI_IN_PLACE, buffer, 2, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
     area_tot = buffer(1)
     pmsl_avg = buffer(2)
     call MPI_Allreduce( MPI_IN_PLACE, pmsl_min, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr )
  endif
#endif

  ! Compute average pressure for any k level with cells above ground
  ! Compute threshold pressure at 80% of the range from avg to min

  pmsl_avg = pmsl_avg / area_tot
  pmsl_thresh = pmsl_avg + .80 * (pmsl_min - pmsl_avg)

  ! Horizontal loop over W points
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Skip current W point if its location is far from first guess position
     if (iskip(iw)) cycle

     ! If pmsl < threshold value, sum area-pressure weighted grid cell location

     if (pmsl(iw) < pmsl_thresh) then

        weight = arw0(iw) * (pmsl_thresh - pmsl(iw)) ! ** exponent  (if desired)

        xew_avg = xew_avg + weight * xew(iw)
        yew_avg = yew_avg + weight * yew(iw)
        zew_avg = zew_avg + weight * zew(iw)

        weight_sum = weight_sum + weight

     endif
  enddo  ! iw

#ifdef OLAM_MPI
  if (iparallel == 1) then
     buffer(1:4) = [xew_avg, yew_avg, zew_avg, weight_sum]
     call MPI_Allreduce( MPI_IN_PLACE, buffer, 4, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
     xew_avg = buffer(1)
     yew_avg = buffer(2)
     zew_avg = buffer(3)
     weight_sum = buffer(4)
  endif
#endif

  ! Compute mean location

  xew_avg = xew_avg / weight_sum
  yew_avg = yew_avg / weight_sum
  zew_avg = zew_avg / weight_sum

  ! Transform mean location to lat/lon coordinates

  call e_ec(xew_avg,yew_avg,zew_avg,rlon,rlat)

  hlat = rlat
  hlon = rlon

  hlat_hist = rlat
  hlon_hist = rlon

  ! This subroutine is called with mstp = 0 only when icycle_hurrinit = 1.
  ! Therefore, the following hlata and hlona values with mstp = 0 and
  ! icycle_hurrinit > 1 are filled separately (in subroutine vortex_relocated).

  hlata(mstp,icycle_hurrinit) = rlat
  hlona(mstp,icycle_hurrinit) = rlon
  htima(mstp,icycle_hurrinit) = real(time8)

  write(io6,'(a,i8,2f10.3,f12.2)') &
     'vortex_center_diagnose END: mstp, hlat, hlon, pmsl_min   ', mstp, hlat, hlon, pmsl_min

  end subroutine vortex_center_diagnose

!==================================================================================

  subroutine vortex_reloc3d()

  use mem_basic,   only: wc, wmc, thil, theta, tair, rr_w, rr_v, rho, press, &
                         vxe, vye, vze
  use mem_micro,   only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, &
                         q2, q6, q7, &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h
  use mem_ijtabs,  only: itab_m, jtab_w, jtw_init
  use mem_grid,    only: mza, mma, mwa, lpw, xem, yem, zem, xew, yew, zew, dzm
  use consts_coms, only: erad, pio180, eps_virt, rocp, cpor, gocp, cvocp, &
                         p00i, p00k, p00kord, t00
  use max_dims,    only: pathlen
  use therm_lib,   only: rhovsl, rhovsl_inv

  implicit none

  real :: reh, zeh, xeh, yeh

  ! PS index transfer

  integer :: iwiflag(mwa)

  integer :: im,iw,j,k,kr,jnext,iwnext,ips,jps,npoly,ipt,jpt,lpt
  integer :: ifld, iwi, n

  real :: rad, radsq, rad0, rpolyi

  real :: x(3),y(3),z(3)
  real :: xw(7),yw(7)

  real :: a(mza,nfld),b(mza,nfld),c(mza,nfld)
  real :: field(mza,7,nfld),field_avg(mza,nfld)

  real :: v0x,v0y,v1x,v1y,v2x,v2y,dot00,dot01,dot02,dot11,dot12,denomi,u,v
  real :: xwi,ywi,rchksq,thv
  real :: dxe, dye, dze

  logical, save :: first_call = .TRUE.

  if (first_call) then
     first_call = .FALSE.

     nout_dim = 0
     nwps(:,:) = 0

     ! Select relocation point and compute its "earth" coordinates

!    hlat_reloc = hlata(0,1)
!    hlon_reloc = hlona(0,1)
!    hlat_reloc = hlata(5,1)
!    hlon_reloc = hlona(5,1)
     hlat_reloc = hlat0
     hlon_reloc = hlon0

     cos_hlat_reloc = cos(hlat_reloc * pio180)
     sin_hlat_reloc = sin(hlat_reloc * pio180)
     cos_hlon_reloc = cos(hlon_reloc * pio180)
     sin_hlon_reloc = sin(hlon_reloc * pio180)

     reh       = erad * cos_hlat_reloc  ! distance from earth axis
     zeh_reloc = erad * sin_hlat_reloc
     xeh_reloc = reh  * cos_hlon_reloc
     yeh_reloc = reh  * sin_hlon_reloc

     rchksq = (rad2_blend + 50.e3)**2

     ! Loop over all W points for recording iw indices in iwps array

     do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

        ! Distance of this IW point from relocation point

!       rad0 = sqrt((xew(iw)-xeh_reloc)**2 + (yew(iw)-yeh_reloc)**2 + (zew(iw)-zeh_reloc)**2)
        radsq = (xew(iw)-xeh_reloc)**2 + (yew(iw)-yeh_reloc)**2 + (zew(iw)-zeh_reloc)**2

        ! Skip hurricane assimilation for all points outside specified radius

        if (radsq >= rchksq) cycle
!       if (rad0 >= rad2_blend + 50.e3) cycle

        ! Transform current W point to PS coordinates tangent at hurricane center

        dxe = xew(iw) - xeh_reloc
        dye = yew(iw) - yeh_reloc
        dze = yew(iw) - zeh_reloc
        call de_ps( dxe,dye,dze,cos_hlat_reloc,sin_hlat_reloc,&
                    cos_hlon_reloc,sin_hlon_reloc,xw(1),yw(1) )

        ! Get ips,jps indices on PS grid (assumed to be 5 km mesh that is 1000 km wide)

        ips = nint((xw(1) + 500.e3) / 5.e3)
        jps = nint((yw(1) + 500.e3) / 5.e3)

        nout_dim = nout_dim + 1
        nwps(ips,jps) = nwps(ips,jps) + 1

        if (nwps(ips,jps) > 60) then
           print*, 'nwps at ',ips,',',jps,' exceeds 60'
           stop 'stop nwps '
        endif

        ! Store iw index in ips,jps element of iwps array to mark its location

        iwps(ips,jps,nwps(ips,jps)) = iw

     enddo

     write(6,'(/,a,3i8)') 'vortex_reloc3d: nout_dim ',nout_dim,mza,nfld

     allocate (iwout(nout_dim))               ; iwout (:) = 0
     allocate (reloc_field(mza,nout_dim,nfld)); reloc_field(:,:,:) = 0.

  endif

  nout = 0

  ! Find "earth" coordinates of hurricane center current location

  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  zeh = erad * sin(hlat * pio180)
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

  write(6,'(/,a,2f10.3)') 'vortex_reloc3d: hlat, hlon ',hlat, hlon
  write(6,'(/,a,2f10.3)') '    hlat_reloc, hlon_reloc ',hlat_reloc, hlon_reloc

  iwiflag(:) = 0

  ! Loop over M points for interpolating W points

  do im = 2,mma

     ! Distance of this IM point from current eye center

     rad = sqrt((xem(im)-xeh)**2 + (yem(im)-yeh)**2 + (zem(im)-zeh)**2)

     ! Skip hurricane assimilation for all points outside specified radius

     if (rad >= rad2_blend + 30.e3) cycle

     ! Transform M point to PS coordinates tangent at current hurricane center

     call e_ps(xem(im),yem(im),zem(im),hlat,hlon,x(1),y(1))

     npoly = itab_m(im)%npoly
     rpolyi = 1. / real(npoly)

     ! Initialize field average

     field_avg(1:mza,1:nfld) = 0.

     ! Loop over all W points that surround current M point

     do j = 1,npoly

        ! Current W point index

        iw = itab_m(im)%iw(j)

        ! Transform W point to PS coordinates tangent at current hurricane center

        call e_ps(xew(iw),yew(iw),zew(iw),hlat,hlon,xw(j),yw(j))

        ! Vertical loop over T levels

        do k = 2,mza
                                 field(k,j, 1) =   vxe(k,iw)
                                 field(k,j, 2) =   vye(k,iw)
                                 field(k,j, 3) =   vze(k,iw)
                                 field(k,j, 4) =   wmc(k,iw)
                                 field(k,j, 5) =    wc(k,iw)
                                 field(k,j, 6) =   rho(k,iw)
                                 field(k,j, 7) = press(k,iw)
                                 field(k,j, 8) =  thil(k,iw)
                                 field(k,j, 9) = theta(k,iw)
                                 field(k,j,10) =  tair(k,iw)
                                 field(k,j,11) =  rr_w(k,iw)
                                 field(k,j,12) =  rr_v(k,iw)
           if (allocated(rr_c))  field(k,j,13) =  rr_c(k,iw)
           if (allocated(rr_d))  field(k,j,14) =  rr_d(k,iw)
           if (allocated(rr_r))  field(k,j,15) =  rr_r(k,iw)
           if (allocated(rr_p))  field(k,j,16) =  rr_p(k,iw)
           if (allocated(rr_s))  field(k,j,17) =  rr_s(k,iw)
           if (allocated(rr_a))  field(k,j,18) =  rr_a(k,iw)
           if (allocated(rr_g))  field(k,j,19) =  rr_g(k,iw)
           if (allocated(rr_h))  field(k,j,20) =  rr_h(k,iw)
           if (allocated(con_c)) field(k,j,21) = con_c(k,iw)
           if (allocated(con_d)) field(k,j,22) = con_d(k,iw)
           if (allocated(con_r)) field(k,j,23) = con_r(k,iw)
           if (allocated(con_p)) field(k,j,24) = con_p(k,iw)
           if (allocated(con_s)) field(k,j,25) = con_s(k,iw)
           if (allocated(con_a)) field(k,j,26) = con_a(k,iw)
           if (allocated(con_g)) field(k,j,27) = con_g(k,iw)
           if (allocated(con_h)) field(k,j,28) = con_h(k,iw)
           if (allocated(q2))    field(k,j,29) =    q2(k,iw)
           if (allocated(q6))    field(k,j,30) =    q6(k,iw)
           if (allocated(q7))    field(k,j,31) =    q7(k,iw)

        enddo

        ! Special for levels below lpw:

        if (lpw(iw) > 2) then

           kr   = lpw(iw)
           thv  = theta(kr,iw) * ( 1.0 + eps_virt * rr_v(kr,iw) )

           do k = kr-1, 2, -1
              ! ignore condensate
              field(k,j,7) = (field(k+1,j,7)**rocp + gocp * dzm(k) * p00k / thv )**cpor
           enddo

           do k = 2, kr-1
              field(k,j, 4) = wmc(kr,iw)
              field(k,j, 5) = wc (kr,iw)
              field(k,j, 6) = field(k,j,7)**cvocp * p00kord / thv
              field(k,j, 8) = theta(kr,iw)  ! ignore condensate
              field(k,j, 9) = theta(kr,iw)
              field(k,j,10) = theta(kr,iw) * (field(k,j,7)*p00i)**rocp
              field(k,j,11) = rr_v (kr,iw)  ! ignore condensate
              field(k,j,12) = rr_v (kr,iw)
              field(k,j,13:31) = 0.0        ! ignore condensate
           enddo

        endif

        do n = 1, nfld
           do k = 2, mza
              field_avg(k,n) = field_avg(k,n) + field(k,j,n) * rpolyi
           enddo
        enddo

     enddo

     ! Loop over all W points that surround current M point and fill field values

     do j = 1,npoly
        jnext = j + 1
        if (j == npoly) jnext = 1

        iw     = itab_m(im)%iw(j)
        iwnext = itab_m(im)%iw(jnext)

        x(2) = xw(j)
        y(2) = yw(j)

        x(3) = xw(jnext)
        y(3) = yw(jnext)

        ! Loop over vertical levels

        do k = 2,mza

           ! Loop over fields

           do ifld = 1,nfld

              z(1) = field_avg(k,ifld)
              z(2) = field(k,j,ifld)
              z(3) = field(k,jnext,ifld)

              ! Evaluate interpolation coefficients for current trio of points

              call matrix_3x3(1.  , x(1), y(1), &
                              1.  , x(2), y(2), &
                              1.  , x(3), y(3), &
                              z(1), z(2), z(3), &
                              a(k,ifld), b(k,ifld), c(k,ifld))

           enddo

        enddo

        ! Set up some triangle-check coefficients

        v0x = x(2) - x(1)
        v0y = y(2) - y(1)

        v1x = x(3) - x(1)
        v1y = y(3) - y(1)

        dot00 = v0x * v0x + v0y * v0y
        dot01 = v0x * v1x + v0y * v1y
        dot11 = v1x * v1x + v1y * v1y

        denomi = 1. / (dot00 * dot11 - dot01 * dot01)

        ! Get ips,jps indices on PS grid (assumed to be 5 km mesh that is 1000 km wide)
        ! for current M point relative to current hurricane location

        ips = nint((x(1) + 500.e3) / 5.e3)
        jps = nint((y(1) + 500.e3) / 5.e3)

        ! From the ips,jps indices, look up nearby W points in iwps table that
        ! have similar PS coordinates relative to relocation point

        do jpt = jps-1,jps+1
           do ipt = ips-1,ips+1
              do lpt = 1,nwps(ipt,jpt)

                 iwi = iwps(ipt,jpt,lpt)

                 if (iwi < 2 .or. iwi > mwa) &
                    write(6,'(a,6i8)') 'iwi out of bounds ',iwi,mwa,ipt,jpt,lpt,ips

                 ! Distance of this IWI point from relocation point

                 rad0 = sqrt((xew(iwi)-xeh_reloc)**2 + (yew(iwi)-yeh_reloc)**2 + (zew(iwi)-zeh_reloc)**2)

                 ! Skip interpolation for all points outside specified radius

                 if (rad0 >= rad2_blend + 20.e3) cycle

                 ! If this point has already been interpolated, cycle

                 if (iwiflag(iwi) > 0) cycle

                 ! Transform IWI point to PS coordinates tangent at relocation point

                 call e_ps(xew(iwi),yew(iwi),zew(iwi),hlat_reloc,hlon_reloc,xwi,ywi)

                 ! Set up remaining triangle_check coefficients

                 v2x = xwi - x(1)
                 v2y = ywi - y(1)

                 dot02 = v0x * v2x + v0y * v2y
                 dot12 = v1x * v2x + v1y * v2y

                 u = (dot11 * dot02 - dot01 * dot12) * denomi
                 v = (dot00 * dot12 - dot01 * dot02) * denomi

                 ! Check if current qx,qy point is inside or very near current triangle

                 if (u > -.01 .and. v > -.01 .and. u + v < 1.01) then

                    ! Point is inside or very near triangle; loop over vertical levels

                    iwiflag(iwi) = 1
                    nout = nout + 1

                    if (nout > nout_dim) then
                       print*, 'In subroutine vortex_reloc3D, nout > nout_dim. '
                       print*, 'nout, nout_dim = ',nout,nout_dim
                       stop 'stop: nout '
                    endif

                    iwout(nout) = iwi

                    do k = 2,mza

                       ! Interpolate to output field point.

                       ! NEW IN MARCH 2022:  (vxe, vye, vze) are now interpolated in place of
                       ! vtan and vrad (wc is still interpolated as before).  For now, we
                       ! neglect rotation of the 3D velocity vector to correct for the small
                       ! change in the unit vertical vector that accompanies the relocation.
                       ! The only impact is on the horizontal velocity since wc is interpolated
                       ! on its own.

                       reloc_field(k,nout,1:nfld) = a(k,1:nfld)       &
                                                  + b(k,1:nfld) * xwi &
                                                  + c(k,1:nfld) * ywi

                    enddo  ! k

                 endif  ! q point inside triangle

              enddo  ! lpt

           enddo   ! ipt

        enddo  ! jpt

     enddo   ! j

   9 continue

  enddo   ! im

  end subroutine vortex_reloc3d

!=======================================================================================

  subroutine vortex_relocated()

  ! This subroutine replaces a model initial state that contains a poorly-resolved
  ! tropical cyclone from a GFS or other analysis with a cyclone that was
  ! simulated in a previous model run and was relocated to the position
  ! inferred in the initial GFS fields.

  use mem_basic,   only: vc, vmc, wc, wmc, thil, theta, tair, &
                         rr_w, rr_v, rho, press, vxe, vye, vze
  use mem_micro,   only: rr_c, rr_d, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, &
                         q2, q6, q7, &
                         con_c, con_d, con_r, con_p, con_s, con_a, con_g, con_h
  use therm_lib,   only: rhovsl
  use mem_ijtabs,  only: itab_v, jtab_v, jtv_init
  use misc_coms,   only: iparallel
  use mem_grid,    only: mza, mwa, xev, yev, zev, xew, yew, zew, &
                         vnxo2, vnyo2, vnzo2, lpw, lve2
  use consts_coms, only: pio180, grav, rvap, rdry, alvlocp, &
                         cvocp, rocp, p00k, p00i, r8, eps_vapi, p00kord
  use max_dims,    only: pathlen
  use olam_mpi_atm, only: mpi_send_w, mpi_recv_w, &
                          mpi_send_v, mpi_recv_v
  use obnd,         only: lbcopy_v, lbcopy_w
  use vel_t3d,      only: diagvel_t3d

  implicit none

  integer :: iout
  integer :: iw,j,k,iv,iw1,iw2
  real    :: wt1, wt1c
  real    :: rad0, vreloc

  real :: vxer(mza,mwa), vyer(mza,mwa), vzer(mza,mwa)

  hlat = hlat_reloc
  hlon = hlon_reloc

  hlata(0,icycle_hurrinit) = hlat_reloc
  hlona(0,icycle_hurrinit) = hlon_reloc

  vxer = vxe
  vyer = vye
  vzer = vze

  ! Horizontal loop over all active W points in file data

  do iout = 1,nout

     iw = iwout(iout)

     ! Distance of this IW point from relocation point

     rad0 = sqrt((xew(iw)-xeh_reloc)**2 + (yew(iw)-yeh_reloc)**2 + (zew(iw)-zeh_reloc)**2)

     ! Define radial weight coefficients, skipping hurricane assimilation for
     ! all points outside specified radius

     if (rad0 > rad2_blend + 10.e3) cycle

     if (rad0 > rad2_blend) then
        wt1 = 0.
     elseif (rad0 < rad1_blend) then
        wt1 = 1.
     else
        wt1 = (rad2_blend - rad0) / (rad2_blend - rad1_blend)
     endif

     wt1c = 1. - wt1

     ! Vertical loop over T levels

     do k = 2,mza

        ! Transfer relocation array data to OLAM arrays, with weighting

                               vxer(k,iw) =                      reloc_field(k,iout, 1)
                               vyer(k,iw) =                      reloc_field(k,iout, 2)
                               vzer(k,iw) =                      reloc_field(k,iout, 3)
                                wmc(k,iw) =   wmc(k,iw) * wt1c + reloc_field(k,iout, 4) * wt1 ! applies at W level
                                 wc(k,iw) =    wc(k,iw) * wt1c + reloc_field(k,iout, 5) * wt1 ! applies at W level
                                rho(k,iw) =   rho(k,iw) * wt1c + reloc_field(k,iout, 6) * wt1
                              press(k,iw) = press(k,iw) * wt1c + reloc_field(k,iout, 7) * wt1
                               thil(k,iw) =  thil(k,iw) * wt1c + reloc_field(k,iout, 8) * wt1
                              theta(k,iw) = theta(k,iw) * wt1c + reloc_field(k,iout, 9) * wt1
                               tair(k,iw) =  tair(k,iw) * wt1c + reloc_field(k,iout,10) * wt1
                               rr_w(k,iw) =  rr_w(k,iw) * wt1c + reloc_field(k,iout,11) * wt1
                               rr_v(k,iw) =  rr_v(k,iw) * wt1c + reloc_field(k,iout,12) * wt1
        if (allocated(rr_c))   rr_c(k,iw) =  rr_c(k,iw) * wt1c + reloc_field(k,iout,13) * wt1
        if (allocated(rr_d))   rr_d(k,iw) =  rr_d(k,iw) * wt1c + reloc_field(k,iout,14) * wt1
        if (allocated(rr_r))   rr_r(k,iw) =  rr_r(k,iw) * wt1c + reloc_field(k,iout,15) * wt1
        if (allocated(rr_p))   rr_p(k,iw) =  rr_p(k,iw) * wt1c + reloc_field(k,iout,16) * wt1
        if (allocated(rr_s))   rr_s(k,iw) =  rr_s(k,iw) * wt1c + reloc_field(k,iout,17) * wt1
        if (allocated(rr_a))   rr_a(k,iw) =  rr_a(k,iw) * wt1c + reloc_field(k,iout,18) * wt1
        if (allocated(rr_g))   rr_g(k,iw) =  rr_g(k,iw) * wt1c + reloc_field(k,iout,19) * wt1
        if (allocated(rr_h))   rr_h(k,iw) =  rr_h(k,iw) * wt1c + reloc_field(k,iout,20) * wt1
        if (allocated(con_c)) con_c(k,iw) = con_c(k,iw) * wt1c + reloc_field(k,iout,21) * wt1
        if (allocated(con_d)) con_d(k,iw) = con_d(k,iw) * wt1c + reloc_field(k,iout,22) * wt1
        if (allocated(con_r)) con_r(k,iw) = con_r(k,iw) * wt1c + reloc_field(k,iout,23) * wt1
        if (allocated(con_p)) con_p(k,iw) = con_p(k,iw) * wt1c + reloc_field(k,iout,24) * wt1
        if (allocated(con_s)) con_s(k,iw) = con_s(k,iw) * wt1c + reloc_field(k,iout,25) * wt1
        if (allocated(con_a)) con_a(k,iw) = con_a(k,iw) * wt1c + reloc_field(k,iout,26) * wt1
        if (allocated(con_g)) con_g(k,iw) = con_g(k,iw) * wt1c + reloc_field(k,iout,27) * wt1
        if (allocated(con_h)) con_h(k,iw) = con_h(k,iw) * wt1c + reloc_field(k,iout,28) * wt1
        if (allocated(q2))       q2(k,iw) =    q2(k,iw) * wt1c + reloc_field(k,iout,29) * wt1
        if (allocated(q6))       q6(k,iw) =    q6(k,iw) * wt1c + reloc_field(k,iout,30) * wt1
        if (allocated(q7))       q7(k,iw) =    q7(k,iw) * wt1c + reloc_field(k,iout,31) * wt1

     enddo

     ! Special for shaved cells with prognostic method:

     if (lve2(iw) > 0) then
        do k = lpw(iw), lpw(iw) + lve2(iw) - 1
           vxe(k,iw) = vxe(k,iw) * wt1c + reloc_field(k,iout, 1) * wt1
           vye(k,iw) = vye(k,iw) * wt1c + reloc_field(k,iout, 2) * wt1
           vze(k,iw) = vze(k,iw) * wt1c + reloc_field(k,iout, 3) * wt1
        enddo
     endif

  enddo

  ! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_w(dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil)

     call mpi_recv_w(dvara1=press, dvara2=rho, &
                     rvara1=wc, rvara2=wmc, rvara3=thil)
  endif

  ! LBC copy (THETA and TAIR will be copied later with the scalars)

  call lbcopy_w(a1=wc, a2=wmc, a3=thil, d1=press, d2=rho)

  ! Initialize VC field

  do j = 1,jtab_v(jtv_init)%jend; iv = jtab_v(jtv_init)%iv(j)
     iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)

     ! Distance of this point from relocation point

     rad0 = sqrt((xev(iv)-xeh_reloc)**2 + (yev(iv)-yeh_reloc)**2 + (zev(iv)-zeh_reloc)**2)

     ! Skip hurricane assimilation for all points outside specified radius

     if (rad0 > rad2_blend) cycle

     ! Define radial weight coefficients

     if (rad0 < rad1_blend) then
        wt1 = 1.
     elseif (rad0 > rad2_blend) then
        wt1 = 0.
     else
        wt1 = (rad2_blend - rad0) / (rad2_blend - rad1_blend)
     endif

     wt1c = 1. - wt1

     ! Vertical loop over T levels

     do k = 2,mza
        vreloc = vnxo2(iv) * (vxer(k,iw1) + vxer(k,iw2)) &
               + vnyo2(iv) * (vyer(k,iw1) + vyer(k,iw2)) &
               + vnzo2(iv) * (vzer(k,iw1) + vzer(k,iw2))

        vc(k,iv) = vc(k,iv) * wt1c + vreloc * wt1
        vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))
     enddo

  enddo

  ! If using MPI, perform parallel send/recv

  if (iparallel == 1) then
     call mpi_send_v(rvara1=vmc, rvara2=vc)
     call mpi_recv_v(rvara1=vmc, rvara2=vc)
  endif

  ! LBC copy of VMC, VC

  call lbcopy_v(vmc=vmc, vc=vc)

  ! Re-diagnose earth-relative velocities

  call diagvel_t3d()

  end subroutine vortex_relocated

!==============================================================================

  subroutine vortex_add_thetapert()

  ! This subroutine adds an incremental theta perturbation to the core of a
  ! tropical cyclone to gradually restore its intensity in cases where it is
  ! poorly resolved in the initial fields (e.g., the GFS analysis).

  ! This subroutine is called only during dynamic initialization cycles.  It
  ! is only called during the last forward cycle if that is the only cycle.

  use mem_ijtabs,   only: jtw_prog, jtab_w
  use mem_tend,     only: thilt, rr_wt
  use mem_basic,    only: rho, rr_v, tair, theta, press
  use misc_coms,    only: mstp, dtlong, time_istp8p
  use mem_grid,     only: lpw, xew, yew, zew
  use consts_coms,  only: erad, pio180, t00, r8
  use therm_lib,    only: rhovsl, rhovsl_inv
  use mem_para,     only: myrank

  implicit none

  real, parameter :: rexpon_thpert = 2.0    ! []

  integer :: iw, j, k

  real :: reh, zeh, xeh, yeh
  real :: rad, delr, wt_horiz
  real :: vterm, pterm, thtend, thtend0, tfact

  ! Call vortex_azim_avg periodically to update vtan_max, which is used to
  ! modulate the imposed heating rate as the vortex intensity approaches the
  ! target value

  if (mod(mstp,4) == 0) then
     call vortex_get_vtanmax()
  endif

  ! Gradually ramp down heating at end of spinup cycle

  tfact = 1.0

  if ( icycle_hurrinit == ncycle_hurrinit .and. &
       time_istp8p > 0.5_r8 * timmax_hurrinit ) then

     tfact = sqrt( max(0., 2.0 - 2.0 * real(time_istp8p) / real(timmax_hurrinit)) )
  endif

  ! Find "earth" coordinates of current hurricane center location (hlat, hlon)

  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  zeh = erad * sin(hlat * pio180)
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

  vterm = max(0., min(1., (vtan_targ - vtan_max ) / 15.0 ) )
  pterm = max(0., min(1., (pmsl_min  - pmsl_targ) / 10.e2) )

  ! geometrically average press and velocity weights
  thtend0 = maxrate_thpert * tfact * (sqrt( (1.0+vterm) * (1.0+pterm) ) - 1.0)

  if (myrank == 0) then
     write(*,*) "Theta Pert: ", thtend0, tfact, vtan_max, vterm, pterm
  endif

  if (thtend0 < 1.e-6) return

  ! Add perturbation to thermodynamic fields

  !$omp parallel do private(iw,rad,delr,wt_horiz,k,thtend) schedule(guided)
  do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)

     ! Distance of this IW point from eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

     delr = abs(rad - rcent_thpert)

     ! Skip for all points outside limiting perturbation radius

     if (delr >= rhwid_thpert) cycle

     wt_horiz = 1. - (delr / rhwid_thpert) ** rexpon_thpert

     do k = max(lpw(iw),kthpert_min), kthpert_max
        thtend      = thtend0 * wt_horiz * wt_vert(k)
        thilt(k,iw) = thilt(k,iw) + thtend * real(rho(k,iw))
     enddo

  enddo
  !$omp end parallel do

  end subroutine vortex_add_thetapert

!==================================================================================

  subroutine vortex_get_vtanmax()

  ! This subroutine computes axisymmetric azimuthal averages of cyclone radial
  ! and azimuthal wind as a function of radius and height. Azimuthal averages
  ! are computed on each model level up to a specified maximum height and over
  ! intervals between radial distance values defined in the radius_ax array.

  ! This works in parallel for huricane tracking

  use mem_ijtabs,  only: jtw_init, jtab_w
  use mem_basic,   only: vxe, vye, vze
  use mem_grid,    only: mza, lpw, xew, yew, zew
  use consts_coms, only: erad, eradi, pio180
  use misc_coms,   only: iparallel

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  integer :: iw, j, k, irad, kmin

  real :: reh, zeh, xeh, yeh
  real :: wnxh,wnyh,wnzh
  real :: wnxrad,wnyrad,wnzrad,wnxtan,wnytan,wnztan
  real :: rad,wrad1,wrad2,radm,radsq
  real :: vtan_ax0

  ! Axisymmmetric vortex profile arrays

  real, target  :: buffer(mza,nr,2)
  real, pointer :: vtan_ax (:,:)  ! tangential wind (m/s)
  real, pointer :: weight_t(:,:)  ! weight array for T points

!!  real :: vtan_ax0 (mza,nr)
!!  real :: weight_t0(mza,nr)

#ifdef OLAM_MPI
  integer :: ierr
#endif

  ! Find "earth" coordinates of hurricane center

  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  zeh = erad * sin(hlat * pio180)
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

  ! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh * eradi
  wnyh = yeh * eradi
  wnzh = zeh * eradi

  ! Initialize axixymmetric arrays to zero prior to summation

  buffer = 0.0

  vtan_ax  => buffer(:,:,1)
  weight_t => buffer(:,:,2)

  ! Loop over all W points

! TODO: get openmp working here?
!!  !$omp parallel private(vtan_ax0,weight_t0)
!!  vtan_ax0  = 0.
!!  weight_t0 = 0.
!!  !$omp do private(iw,radsq,rad,radm,irad,wrad2,wrad1,wnxrad,wnyrad,wnzrad,&
!!  !$omp            &wnxtan,wnytan,wnztan,k,vtan_ax0) schedule(guided)
  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Distance of this IW point from eye center

!    rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)
     radsq = (xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2

     ! Skip axisymmetric averaging for all points outside specified radius

!    if (rad >= radius_ax(nr) - 1.) cycle
     if (radsq >= rmaxsq) cycle

     rad  = sqrt(radsq)
     radm = 1.0 / rad

     ! Determine interpolation point in radial dimension

     irad = 1
     do while (rad > radius_ax(irad+1))
        irad = irad + 1
     enddo

     ! Make sure irad is in array bounds

     if (irad < 1 .or. irad+1 > nr) cycle

     wrad2 = (rad - radius_ax(irad)) * dradi_ax(irad)
     wrad1 = 1. - wrad2

     ! Unit normal vector components from hurricane center to current IW point

     wnxrad = (xew(iw) - xeh) * radm
     wnyrad = (yew(iw) - yeh) * radm
     wnzrad = (zew(iw) - zeh) * radm

     ! Unit vector components in direction of tangential vortex wind

     wnxtan = wnyh * wnzrad - wnzh * wnyrad
     wnytan = wnzh * wnxrad - wnxh * wnzrad
     wnztan = wnxh * wnyrad - wnyh * wnxrad

     ! Vertical loop over T levels

     do k = lpw(iw), l8k
        weight_t(k,irad)   = weight_t(k,irad)   + wrad1
        weight_t(k,irad+1) = weight_t(k,irad+1) + wrad2

        ! Diagnose axisymmetric component of model's own vortex

        vtan_ax0 = vxe(k,iw) * wnxtan + vye(k,iw) * wnytan + vze(k,iw) * wnztan

        vtan_ax(k,irad)   = vtan_ax(k,irad)   + wrad1 * vtan_ax0
        vtan_ax(k,irad+1) = vtan_ax(k,irad+1) + wrad2 * vtan_ax0
     enddo
  enddo
!!  !$omp end parallel do

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Allreduce( MPI_IN_PLACE, buffer, mza*nr*2, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
  endif
#endif

  ! Convert sums to averages

  vtan_max = 0.

  do irad = nr, 2, -1

     do k = 2, l8k
        if (weight_t(k,irad) > 1.e-6) exit
     enddo
     kmin = k

     do k = kmin, l8k
        vtan_ax(k,irad) = vtan_ax(k,irad) / weight_t(k,irad)
     enddo

     if (kmin <= l8k) then
        vtan_max = max(vtan_max, maxval(vtan_ax(kmin:l8k,irad)))
     endif

  enddo

end subroutine vortex_get_vtanmax

!==================================================================================

  subroutine vortex_azim_avg(plt)

  ! This subroutine computes axisymmetric azimuthal averages, centered around the
  ! latitude & longitude of the minimum PMSL pressure in a tropical cyclone, of
  ! cyclone radial, azimuthal, and vertical wind components and a few scalar
  ! fields as a function of radius and height.  Azimuthal averages are computed
  ! on each model level up to a specified maximum height and over intervals
  ! between radial distance values defined in the radius_ax array.

  use mem_ijtabs,  only: jtw_init, jtab_w
  use mem_basic,   only: thil, theta, tair, rho, rr_w, rr_v, wc, vxe, vye, vze
  use mem_grid,    only: mza, zt, xew, yew, zew, lpw
  use consts_coms, only: erad, pio180, alvlocp
  use therm_lib,   only: rhovsl
  use misc_coms,   only: iparallel
  use mem_para,    only: myrank
  use plotcolors,  only: make_colortable

  implicit none

  character(*), intent(in) :: plt

  integer :: iw, j, k, irad, nzz, kmin
  real :: reh, zeh, xeh, yeh
  real :: wnxh,wnyh,wnzh
  real :: wnxrad,wnyrad,wnzrad,wnxtan,wnytan,wnztan

  real :: rad,wrad1,wrad2
  real :: vtan_ax0, vrad_ax0, ss_liq, ss_thetadif

  ! Axisymmmetric vortex profile arrays

  real ::  thil_ax(mza,nr) ! ice-liquid potential temperature (K)
  real :: theta_ax(mza,nr) ! potential temperature (K)
  real ::  tair_ax(mza,nr) ! temperature (K)
  real ::   rrw_ax(mza,nr) ! total water mixing ratio (kg/kg_dryair)
  real ::   rrv_ax(mza,nr) ! vapor mixing ratio (kg/kg_dryair)
  real ::  vtan_ax(mza,nr) ! tangential wind (m/s)
  real ::  vrad_ax(mza,nr) ! radial wind (m/s)
  real ::     w_ax(mza,nr) ! vertical wind (m/s)
  real :: ssliq_ax(mza,nr) ! supsat wrt liquid (%)
  real :: thdif_ax(mza,nr) ! theta difference from supsat (K)
  real ::  cond_ax(mza,nr) ! condensation mixing ratio (kg/kg_dryair)

  real :: weight_t(mza,nr) ! weight array for T points

  ! Needs work in parallel!!!!!!!!!!!
  if (iparallel == 1) return
  if (myrank /= 0) return

  ! Find "earth" coordinates of hurricane center

  reh = erad * cos(hlat * pio180)  ! distance from earth axis
  zeh = erad * sin(hlat * pio180)
  xeh = reh  * cos(hlon * pio180)
  yeh = reh  * sin(hlon * pio180)

  ! Components of unit vector outward normal to earth surface at hurricane center

  wnxh = xeh / erad
  wnyh = yeh / erad
  wnzh = zeh / erad

  ! Initialize axixymmetric arrays to zero prior to summation

   thil_ax(:,:) = 0.
  theta_ax(:,:) = 0.
   tair_ax(:,:) = 0.
    rrw_ax(:,:) = 0.
    rrv_ax(:,:) = 0.
   vtan_ax(:,:) = 0.
   vrad_ax(:,:) = 0.
      w_ax(:,:) = 0.
  ssliq_ax(:,:) = 0.
  thdif_ax(:,:) = 0.

  weight_t(:,:) = 0.

  ! Loop over all W points

  do j = 1,jtab_w(jtw_init)%jend; iw = jtab_w(jtw_init)%iw(j)

     ! Distance of this IW point from eye center

     rad = sqrt((xew(iw)-xeh)**2 + (yew(iw)-yeh)**2 + (zew(iw)-zeh)**2)

     ! Skip axisymmetric averaging for all points outside specified radius

     if (rad >= radius_ax(nr) - 1.0) cycle

     ! Determine interpolation point in radial dimension

     irad = 1
     do while (rad > radius_ax(irad+1))
        irad = irad + 1
     enddo

     if (irad < 1 .or. irad+1 > nr) then
        print*, 'irad out of bounds1 ', irad, nr
        cycle
     endif

     wrad2 = (rad - radius_ax(irad)) * dradi_ax(irad)
     wrad1 = 1. - wrad2

     ! Unit normal vector components from hurricane center to current IW point

     wnxrad = (xew(iw) - xeh) / rad
     wnyrad = (yew(iw) - yeh) / rad
     wnzrad = (zew(iw) - zeh) / rad

     ! Unit vector components in direction of tangential vortex wind

     wnxtan = wnyh * wnzrad - wnzh * wnyrad
     wnytan = wnzh * wnxrad - wnxh * wnzrad
     wnztan = wnxh * wnyrad - wnyh * wnxrad

     ! Vertical loop over T levels

     do k = lpw(iw),mza

        weight_t(k,irad)   = weight_t(k,irad)   + wrad1
        weight_t(k,irad+1) = weight_t(k,irad+1) + wrad2

        ! Diagnose axisymmetric component of model's own vortex

         thil_ax(k,irad)   =  thil_ax(k,irad)   + wrad1 *  thil(k,iw)
         thil_ax(k,irad+1) =  thil_ax(k,irad+1) + wrad2 *  thil(k,iw)

        theta_ax(k,irad)   = theta_ax(k,irad)   + wrad1 * theta(k,iw)
        theta_ax(k,irad+1) = theta_ax(k,irad+1) + wrad2 * theta(k,iw)

         tair_ax(k,irad)   =  tair_ax(k,irad)   + wrad1 *  tair(k,iw)
         tair_ax(k,irad+1) =  tair_ax(k,irad+1) + wrad2 *  tair(k,iw)

          rrw_ax(k,irad)   =   rrw_ax(k,irad)   + wrad1 *  rr_w(k,iw)
          rrw_ax(k,irad+1) =   rrw_ax(k,irad+1) + wrad2 *  rr_w(k,iw)

          rrv_ax(k,irad)   =   rrv_ax(k,irad)   + wrad1 *  rr_v(k,iw)
          rrv_ax(k,irad+1) =   rrv_ax(k,irad+1) + wrad2 *  rr_v(k,iw)

            w_ax(k,irad)   =     w_ax(k,irad)   + wrad1 *    wc(k,iw)
            w_ax(k,irad+1) =     w_ax(k,irad+1) + wrad2 *    wc(k,iw)

        ss_liq = max(0.,(rr_v(k,iw) * real(rho(k,iw)) / rhovsl(tair(k,iw)-273.15) - 1.0) * 1.e2)

        ssliq_ax(k,irad)   = ssliq_ax(k,irad)   + wrad1 * ss_liq
        ssliq_ax(k,irad+1) = ssliq_ax(k,irad+1) + wrad2 * ss_liq

        ss_thetadif = min(0.,(rhovsl(tair(k,iw)-273.15) - rr_v(k,iw) * real(rho(k,iw)))) &
                    * alvlocp / real(rho(k,iw))

        thdif_ax(k,irad)   = thdif_ax(k,irad)   + wrad1 * ss_thetadif
        thdif_ax(k,irad+1) = thdif_ax(k,irad+1) + wrad2 * ss_thetadif

        ! Diagnose axisymmetric component of model's own vortex

        vtan_ax0 = vxe(k,iw) * wnxtan + vye(k,iw) * wnytan + vze(k,iw) * wnztan
        vrad_ax0 = vxe(k,iw) * wnxrad + vye(k,iw) * wnyrad + vze(k,iw) * wnzrad

        vtan_ax(k,irad)   = vtan_ax(k,irad)   + wrad1 * vtan_ax0
        vtan_ax(k,irad+1) = vtan_ax(k,irad+1) + wrad2 * vtan_ax0

        vrad_ax(k,irad)   = vrad_ax(k,irad)   + wrad1 * vrad_ax0
        vrad_ax(k,irad+1) = vrad_ax(k,irad+1) + wrad2 * vrad_ax0

     enddo

  enddo

  ! Convert sums to averages

!! vtan_max     = 0.
  vtan_ax(:,1) = 0.
  vrad_ax(:,1) = 0.

  do irad = nr,1,-1

     do k = 2, mza
        if (weight_t(k,irad) > 1.e-6) exit
     enddo
     kmin = k

     do k = kmin, mza
         thil_ax(k,irad) =  thil_ax(k,irad) / weight_t(k,irad)
        theta_ax(k,irad) = theta_ax(k,irad) / weight_t(k,irad)
         tair_ax(k,irad) =  tair_ax(k,irad) / weight_t(k,irad)
          rrw_ax(k,irad) =   rrw_ax(k,irad) / weight_t(k,irad)
          rrv_ax(k,irad) =   rrv_ax(k,irad) / weight_t(k,irad)
            w_ax(k,irad) =     w_ax(k,irad) / weight_t(k,irad)
        ssliq_ax(k,irad) = ssliq_ax(k,irad) / weight_t(k,irad)
        thdif_ax(k,irad) = thdif_ax(k,irad) / weight_t(k,irad)

        if (irad > 1) then
           vtan_ax(k,irad) = vtan_ax(k,irad) / weight_t(k,irad)
           vrad_ax(k,irad) = vrad_ax(k,irad) / weight_t(k,irad)
        endif
     enddo

     if (kmin > 2) then

        if (irad == nr) then

           if (kmin > mza) stop 'stop irad T'

           do k = 2, kmin-1
               thil_ax(k,irad) =  thil_ax(kmin,irad)
              theta_ax(k,irad) = theta_ax(kmin,irad)
               tair_ax(k,irad) =  tair_ax(kmin,irad)
                rrw_ax(k,irad) =   rrw_ax(kmin,irad)
                rrv_ax(k,irad) =   rrv_ax(kmin,irad)
                  w_ax(k,irad) =     w_ax(kmin,irad)
              ssliq_ax(k,irad) = ssliq_ax(kmin,irad)
              thdif_ax(k,irad) = thdif_ax(kmin,irad)
               vtan_ax(k,irad) =  vtan_ax(kmin,irad)
               vrad_ax(k,irad) =  vrad_ax(kmin,irad)
           enddo

        else

           do k = 2, kmin-1
               thil_ax(k,irad) =  thil_ax(k,irad+1)
              theta_ax(k,irad) = theta_ax(k,irad+1)
               tair_ax(k,irad) =  tair_ax(k,irad+1)
                rrw_ax(k,irad) =   rrw_ax(k,irad+1)
                rrv_ax(k,irad) =   rrv_ax(k,irad+1)
                  w_ax(k,irad) =     w_ax(k,irad+1)
              ssliq_ax(k,irad) = ssliq_ax(k,irad+1)
              thdif_ax(k,irad) = thdif_ax(k,irad+1)

              if (irad > 1) then
                 vtan_ax(k,irad) = vtan_ax(k,irad+1) &
                                 * radius_ax(irad) / radius_ax(irad+1)
                 vrad_ax(k,irad) = vrad_ax(k,irad+1) &
                                 * radius_ax(irad) / radius_ax(irad+1)
              endif
           enddo

        endif
     endif

!     if (kmin <= mza .and. irad /= 1) then
!        vtan_max = max(vtan_max, maxval(vtan_ax(kmin:mza,irad)))
!     endif

  enddo

  ! print*, 'vtan_max ', vtan_max

  if (trim(plt) /= 'plot') return
  if (myrank /= 0) return

  cond_ax(:,:) = rrw_ax(:,:) - rrv_ax(:,:)

  ! Plot up to selected height

  do k = 2,mza
     nzz = k
     if (zt(k) > 20000.) exit
  enddo

  call make_colortable(500,'lin',280.0,500.0,5.0,1.)
  call make_colortable(501,'lin',  0.0, 80.0,5.0,1.)
  call make_colortable(502,'lin',-20.0, 50.0,2.0,1.)
  call make_colortable(503,'lin', -2.0,  8.0,0.2,1.)
  call make_colortable(504,'lin',250.0,310.0,2.0,1.)
  call make_colortable(505,'lin',  0.0, 30.0,1.0,.1)

  call o_reopnwk()

  !                  panel label    units      field   factor lbc ctab

  call plotback()
  call vortex_rzplot(nzz,'3','theta1 '  ,' (K)'   , theta_ax ,1.   ,1 ,500)
  call vortex_rzplot(nzz,'4','vtan1 '   ,' (m/s)' ,  vtan_ax ,1.   ,0 ,501)
  call vortex_rzplot(nzz,'1','vrad1 '   ,' (m/s)' ,  vrad_ax ,1.   ,0 ,502)
  call vortex_rzplot(nzz,'2','w1 '      ,' (m/s)' ,     w_ax ,1.   ,0 ,503)
  call o_frame()

  call plotback()
  call vortex_rzplot(nzz,'3','tair1 '   ,' (K)'   ,  tair_ax ,1.   ,1 ,504)
  call vortex_rzplot(nzz,'4','rrw1 '    ,' (g/kg)',   rrw_ax ,1.e3 ,1 ,505)
  call vortex_rzplot(nzz,'1','rrv1 '    ,' (g/kg)',   rrv_ax ,1.e3 ,1 ,505)
  call vortex_rzplot(nzz,'2','cond1 '   ,' (g/kg)',  cond_ax ,1.e3 ,1 ,505)
  call o_frame()

!  call vortex_rzplot(nzz,'0','ssliq1 '  ,' (%)'   , ssliq_ax ,1.   ,0 ,443)
!  call vortex_rzplot(nzz,'0','ss_thdif ',' (K)'   , thdif_ax ,1.   ,0 ,464)

  call o_clswk()

  end subroutine vortex_azim_avg

!==================================================================================

  subroutine vortex_rzplot(nzz, panel, label, units, fieldin, factor, lbc, colortab)

  ! This subroutine is a wrapper for plotting radius-height arrays of
  ! dimension (nzz,nr) that contain azimuthal averages of the TC vortex.
  ! It calls oplot_zxy2 to carry out the actual plot.

  use misc_coms
  use mem_grid,   only: mza, zm, zt
  use consts_coms
  use max_dims,   only: pathlen
  use mem_para,   only: myrank

  implicit none

  integer,      intent(in) :: nzz, lbc, colortab
  character(*), intent(in) :: panel, label, units
  real,         intent(in) :: fieldin(mza,nr)
  real,         intent(in) :: factor

  real :: radius(nr), height(nzz), field(nzz,nr)

  integer :: ifill = 1

  real, save :: aspect = 1.0
  real, save :: scalelab = .014

  if (myrank /= 0) return

  radius(1:nr) = radius_ax(1:nr) * 1.e-3

  height(1:nzz) = zt(1:nzz) * 1.e-3
  height(1) = zm(1) * 1.e-3

  field(2:nzz,1:nr) = fieldin(2:nzz,1:nr) * factor
  field(1,1:nr) = field(2,1:nr)
  if (lbc == 0) field(1,1:nr) = 0.

  call oplot_zxy2(trim(panel),'N','a','c',aspect,scalelab,                 &
                  trim(label),trim(units), nr, nzz, radius, height, &
                  'radius (km)','height (km)', &
                  field,    & ! 11N
                  colortab,ifill,0.,radius(nr),20.,5,               &
                  0.0,height(nzz),1.0,5  )

  end subroutine vortex_rzplot

!==================================================================================

  subroutine vortex_trajec_plot(iplt)

  ! This subroutine plots a time series of hurricane locations for present simulation

  use misc_coms,   only: mstp, dtlong
  use oplot_coms,  only: op
  use consts_coms, only: pio180, erad
  use mem_para,    only: myrank

  implicit none

  integer, intent(in) :: iplt

  integer :: icolor, lhour, khour, kstp, kstp_max, jcyc
  real :: reh, zeh, xeh, yeh, rhour, bsize, xs, ys
  character(len=2) :: title

  if (myrank /= 0) return

  bsize = .016 * (op%hp2 - op%hp1) ! * 0.3

  do jcyc = 1,icycle_hurrinit

     if     (ncycle_hurrinit - jcyc == 0) then
        icolor = 127
     elseif (ncycle_hurrinit - jcyc == 1) then
        icolor = 135
     elseif (ncycle_hurrinit - jcyc == 2) then
        icolor = 139
     elseif (ncycle_hurrinit - jcyc == 3) then
        icolor = 102
     else
        icolor = 114
     endif

     call o_sflush()

     call o_gsplci(icolor)
     call o_gstxci(icolor)
     call o_gsfaci(icolor)

     if (jcyc < icycle_hurrinit) then
        kstp_max = int(timmax_hurrinit / dtlong)
     else
        kstp_max = mstp
     endif

     lhour = 0

     do kstp = 0, kstp_max

        ! Find "earth" coordinates of sims hurricane center

        reh = erad * cos(hlata(kstp,jcyc) * pio180)  ! distance from earth center
        zeh = erad * sin(hlata(kstp,jcyc) * pio180)
        xeh = reh  * cos(hlona(kstp,jcyc) * pio180)
        yeh = reh  * sin(hlona(kstp,jcyc) * pio180)

        ! Transform hurricane earth coords to whatever projection is in use

        call oplot_transform(iplt,xeh,yeh,zeh,xs,ys)

        if (kstp == 0) then
           call o_frstpt(xs,ys)
        else
           call o_vector(xs,ys)
        endif

        rhour = htima(kstp,jcyc) / 3600.
        khour = int(rhour)

        if (khour > 0 .and. khour > lhour) then
           write(title,'(i2)') khour
           call o_plchhq (xs,ys,trim(adjustl(title)),bsize,0.,0.)
           lhour = khour
        endif

     enddo
  enddo

  end subroutine vortex_trajec_plot

end module hcane_rz
