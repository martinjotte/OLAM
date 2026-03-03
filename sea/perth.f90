module perth

! Module for the perth5 tide prediction code

! Purpose: to compute the ocean tide height at a given time and location
! from grids of harmonic constants.

! The main user-callable subroutine is perth5.
!
!          SUBROUTINE PERTH5( PLAT, PLON, TIMEMJD, TIDE, ISDATA )
! where:
!   PLAT,PLON - latitude and longitude in degrees for desired location.
!   TIMEMJD   - desired UT time, in decimal Modified Julian Day.
!   TIDE      - the computed tidal height, in same units as input amplitude grids.
!   ISDATA    - logical denoting whether tide data exist at the desired location
!
! Before the PERTH5 subroutine is called, an initialization subroutine must
! be called.  It has some optional arguments, although the defaults should be
! taken for most purposes.  All arguments are input, none are output.
!
!     SUBROUTINE PERTH5_INITIALIZE( DATASETNAME )
! where:
!    DATASETNAME - an array of file names, each file holding one constituent.
!
! R Ray, Dec 2017 - pieced together from a number of older routines.
!
! R Walko, Nov 2024 - extracted a subset of the perth5 software package and hardwired
!                     it for the specific purpose of computing tides from combined
!                     versions 5.5 and 5.6 gridded tide component data (20 files in
!                     all).  This trimmed-down perth5 package is used as a component
!                     of the Ocean-Land-Atmosphere Model (OLAM) to provide tides for
!                     its shallow-water coastal storm surge model.  OLAM reads the
!                     tide component files with hdf5 software.

  implicit none

  private
  public  :: perth5_initialize, perth5_interp, perth5_endinit, perth5_temporal, &
             perth5_tide, ntides

  integer, parameter :: ntides = 20            ! Number of tide constituents (GOT5.5 + GOT5.6)
  integer, parameter :: nx = 2881, ny = 1441   ! dimensions of each constituent in database files
  real, parameter    :: dx = 0.125, dy = 0.125 ! lon/lat increment in constituents
  real, parameter    :: undef = 9999.          ! value of undefined data in database files
  real, parameter    :: undefm1 = 9998.        ! one less than undef 
  integer            :: doodno(7,ntides)       ! 7-vector Doodson number (integers, without 5's)
  character(len=6)   :: constituent(ntides)
  character(len=100) :: datasetname(ntides)

  real, parameter :: pi=3.1415926, pio180=pi/180., piu180 = 180./pi

  real, allocatable :: amp(:,:), phs(:,:)           ! amplitude and phase read from tide database files
  real, allocatable :: acospd(:,:,:), asinpd(:,:,:) ! amp*cos(phs) and amp*sin(phs) of database files 

  real(8) :: f(ntides)    ! (depends on time, not location)
  real(8) :: u(ntides)    ! (depends on time, not location)
  real(8) :: args(ntides) ! (depends on time, not location)
  real(8) :: doodbeta(7)  ! Doodson's 6(+1) tidal arguments (depends on time, not location)

CONTAINS

  subroutine perth5_initialize( tide_database )

     use hdf5_utils, only: shdf5_open, shdf5_exists, shdf5_irec, shdf5_close

     ! Initialize the perth5 codes; must be called once before calling perth5.

     implicit none
     character(len=*), intent(in) :: tide_database

     integer :: i, j, k
     integer :: ndims, idims(2)
     logical :: exists

     constituent(1)  = 'Q1_f'   ; doodno(:,1)  = [ 1,-2, 0, 1, 0, 0, 3 ]
     constituent(2)  = 'O1_f'   ; doodno(:,2)  = [ 1,-1, 0, 0, 0, 0, 3 ]
     constituent(3)  = 'P1_f'   ; doodno(:,3)  = [ 1, 1,-2, 0, 0, 0, 3 ]
     constituent(4)  = 'S1_f'   ; doodno(:,4)  = [ 1, 1,-1, 0, 0, 0, 2 ]
     constituent(5)  = 'K1_f'   ; doodno(:,5)  = [ 1, 1, 0, 0, 0, 0, 1 ]
     constituent(6)  = 'N2_f'   ; doodno(:,6)  = [ 2,-1, 0, 1, 0, 0, 0 ]
     constituent(7)  = 'M2_f'   ; doodno(:,7)  = [ 2, 0, 0, 0, 0, 0, 0 ]
     constituent(8)  = 'S2_f'   ; doodno(:,8)  = [ 2, 2,-2, 0, 0, 0, 0 ]
     constituent(9)  = 'K2_f'   ; doodno(:,9)  = [ 2, 2, 0, 0, 0, 0, 0 ]
     constituent(10) = 'M4_f'   ; doodno(:,10) = [ 4, 0, 0, 0, 0, 0, 0 ] ! new in 5.5
     constituent(11) = 'MS4_f'  ; doodno(:,11) = [ 4, 2,-2, 0, 0, 0, 0 ] ! new in 5.5
     constituent(12) = '2N2_f'  ; doodno(:,12) = [ 2,-2, 0, 2, 0, 0, 0 ] ! new in 5.5
     constituent(13) = 'mu2_f'  ; doodno(:,13) = [ 2,-2, 2, 0, 0, 0, 0 ] ! new in 5.5
     constituent(14) = 'J1_f'   ; doodno(:,14) = [ 1, 2, 0,-1, 0, 0, 1 ] ! new in 5.5
     constituent(15) = 'sig1_f' ; doodno(:,15) = [ 1,-3, 2, 0, 0, 0, 3 ] ! new in 5.5
     constituent(16) = 'OO1_f'  ; doodno(:,16) = [ 1, 3, 0, 0, 0, 0, 1 ] ! new in 5.5
     constituent(17) = 'L2p_f'  ; doodno(:,17) = [ 2, 1, 0, 0, 0, 0, 3 ] ! new in 5.6; deg-3
     constituent(18) = 'M1p_f'  ; doodno(:,18) = [ 1, 0, 0, 0, 0, 0, 2 ] ! new in 5.6; deg-3
     constituent(19) = 'M3_f'   ; doodno(:,19) = [ 3, 0, 0, 0, 0, 0, 2 ] ! new in 5.6; deg-3
     constituent(20) = 'N2p_f'  ; doodno(:,20) = [ 2,-1, 0, 0, 0, 0, 1 ] ! new in 5.6; deg-3

     write(6,'(a)') 'Initializing module PERTH for model-based tide prediction'

     write(6,30) nx, ny, ntides
     30 format(' Tide atlas grid size: nx =',i6,'  x  ny =',i6,i8,' constituents')

     allocate( amp(nx,ny), phs(nx,ny), acospd(nx,ny,ntides), asinpd(nx,ny,ntides) )
     acospd(:,:,:) = undef
     asinpd(:,:,:) = undef

     ! Read and store all tide grids

     ndims = 2
     idims(1) = nx
     idims(2) = ny

     do i = 1,ntides
        datasetname(i) = trim(tide_database) // trim(constituent(i)) // '.h5'
        call shdf5_exists(datasetname(i), exists)

        if (exists) then
           write(6,'(a)') 'Reading tide dataset '//trim(datasetname(i))
           call shdf5_open(datasetname(i),'R',serial=.true.)
           call shdf5_irec(ndims,idims,'amplitude',rvar2=amp)
           call shdf5_irec(ndims,idims,'phase'    ,rvar2=phs)
           call shdf5_close()

           do j = 1,ny
              do k = 1,nx
                 if (amp(k,j) < undefm1) then
                    acospd(k,j,i) = amp(k,j)*cos(phs(k,j)*pio180)
                    asinpd(k,j,i) = amp(k,j)*sin(phs(k,j)*pio180)
                 endif
              enddo
           enddo
        else
           print*, 'Trying to open tide constituent file:'
           print*, trim(datasetname(i))
           print*, 'but it does not exist. The run is ended.'
           stop
        endif
     end do

     f = 1.d0;  u = 0.d0

     write(6,*) 'Subroutine perth5_initialize completed'

  end subroutine perth5_initialize

!==================================================================================

  subroutine perth5_interp( plat, plon, acosp, asinp, isdata )

     implicit none

     real,    intent(in)  :: plat, plon                   ! lat/lon coordinates to interpolate to
     real,    intent(out) :: acosp(ntides), asinp(ntides) ! interpolated values of acospd and asinpd
     logical, intent(out) :: isdata                       ! true if interpolation is successful

     real            :: qlat, qlon, rlat, rlon
     real, parameter :: lonmin = 0., latmin = -90.
     real, parameter :: dx = 0.125, dy = 0.125

     real :: rio, rjo, wx1, wx2, wy1, wy2, w, wsum
     integer :: i, io1, io2, jo1, jo2

     isdata = .true.

     ! Compute indices for desired position

     qlon = plon
     if (qlon >= 360.) qlon = qlon - 360.
     if (qlon <    0.) qlon = qlon + 360.
     qlon = max(  0.0001,min(359.9999,qlon))
     qlat = max(-89.9999,min( 89.9999,plat))

     rjo = (qlat - latmin) / dy + 1.
     jo1 = int(rjo)
     jo2 = jo1 + 1

     rio = (qlon - lonmin) / dx + 1.
     io1 = int(rio)
     io2 = io1 + 1

     wx2 = rio - real(io1)
     wx1 = 1. - wx2
     wy2 = rjo - real(jo1)
     wy1 = 1. - wy2

     wsum = 0.0

     acosp(1:ntides) = 0.
     asinp(1:ntides) = 0.

     if (abs(acospd(io1,jo1,1)) < undefm1) then
        w    = wx1 * wy1
        wsum = w
        do i = 1,ntides
           acosp(i) = w * acospd(io1,jo1,i)
           asinp(i) = w * asinpd(io1,jo1,i)
        enddo
     endif

     if (abs(acospd(io2,jo1,1)) < undefm1) then
        w    = wx2 * wy1
        wsum = wsum + w
        do i = 1,ntides
           acosp(i) = acosp(i) + w * acospd(io2,jo1,i)
           asinp(i) = asinp(i) + w * asinpd(io2,jo1,i)
        enddo
     endif

     if (abs(acospd(io1,jo2,1)) < undefm1) then
        w    = wx1 * wy2
        wsum = wsum + w
        do i = 1,ntides
           acosp(i) = acosp(i) + w * acospd(io1,jo2,i)
           asinp(i) = asinp(i) + w * asinpd(io1,jo2,i)
        enddo
     endif

     if (abs(acospd(io2,jo2,1)) < undefm1) then
        w    = wx2 * wy2
        wsum = wsum + w
        do i = 1,ntides
           acosp(i) = acosp(i) + w * acospd(io2,jo2,i)
           asinp(i) = asinp(i) + w * asinpd(io2,jo2,i)
        enddo
     endif

!     if (wsum > 0.5) then ! Allow lower threshold?
     if (wsum > 1.e-9) then
        acosp(1:ntides) = acosp(1:ntides) / wsum
        asinp(1:ntides) = asinp(1:ntides) / wsum
     else
        isdata = .false.
     endif

   end subroutine perth5_interp

!==================================================================================

  subroutine perth5_endinit()

     implicit none

     if (allocated(amp))    deallocate(amp)
     if (allocated(phs))    deallocate(phs)
     if (allocated(acospd)) deallocate(acospd)
     if (allocated(asinpd)) deallocate(asinpd)

  end subroutine perth5_endinit

!==================================================================================

  subroutine perth5_temporal(timemjd)

     ! Updates variables (f, u, args) that depend on time but not geographical location

     implicit none

     real(8), intent(in) :: timemjd

     integer :: i
     real(8) :: perigee, omega

     call doodson_vector( timemjd, doodbeta(1:6) )
     doodbeta(7) = 90.d0

     perigee =  doodbeta(4)
     omega   = -doodbeta(5)

     call fu_nodal6( omega, perigee, constituent, ntides, f, u )

     ! Compute Doodson's tidal arguments

     do i = 1,ntides
        args(i) = DOT_PRODUCT( doodbeta, dble(doodno(:,i)) )
        args(i) = modulo( args(i), 360.d0 )
     end do

  end subroutine perth5_temporal

!==================================================================================

  subroutine doodson_vector( time_UT, doodbeta )

     ! Routine to evaluate Doodson's 6 astronomical variables.
     ! Computed angles are returned in units of degrees, between [0,360)
     !    doodbeta = [ tau, s, h, p, N', p_s ]
     !
     !    time_UT - Universal Time in decimal Modified Julian Days
     !
     ! This routine computes required astronomical mean longitudes by
     ! calling the IERS's FUNDARG routine.  By this means, any updates by
     ! the IERS to FUNDARG can be readily incorporated here.

     implicit none
     real(8), intent(in)  ::  time_UT
     real(8), intent(out) ::  doodbeta(6)

     real(8), parameter :: pi=3.14159265358979d0, rad=pi/180.d0, twopi=2.d0*pi
     real(8)  L, LP, F, D, omega
     real(8)  tx, tjd, tau, s, h, p, fn, ps, time_TT, tsolar

     real(8), parameter :: tt_ut = 75.d0 ! Terrestrial dynamical Time minus Universal Time [sec]
                                         ! (increasing by year; accurate to within 5s until 2034)

     ! Terrestrial Time:
     time_TT = time_UT + tt_ut/86400.d0

     tjd = time_TT + 2400000.5d0
     tx = ( tjd - 2451545.d0 )/36525.d0  !  TT since J2000, in centuries

     call fundarg( tx, L, LP, F, D, omega )

     s = F + omega                      !  mean longitude of moon
     h = F + omega - D                  !  mean longitude of sun
     p = F + omega - L                  !  longitude of lunar perigee
     fn = -omega                        !  negative longitude of lunar node
     ps = -LP + F - D + omega           !  longitude of solar perigee

     tsolar = (time_UT - int(time_UT))*twopi
     tau = tsolar - s + h

     doodbeta = [ tau, s, h, p, fn, ps ] / rad
     doodbeta = modulo( doodbeta, 360.d0 )

  end subroutine doodson_vector

!==================================================================================

  recursive subroutine fu_nodal6( omega, p, constituent, n, f, u )

     !  Computes tidal nodal (& perigee) corrections "f" and "u" 
     !  for n tidal constituents, given the mean longitudes p and omega.

     !-------------------------------------------------------
     !  Input:
     !     omega - mean longitude of lunar node, in degrees.
     !     p     - mean longitude of lunar perigee.
     !     constituent - array of constituent names.
     !     n  - length of arrays {f,u,constituent}.
     !  Output:
     !     f  - modulation factor for given tide(s).
     !     u  - phase correction for given tide, in degrees.
     !-------------------------------------------------------
     !
     !  Note: omega = -N' (where N' is the 5th Doodson variable); 
     !       omega is decreasing in time whereas N' is increasing.
     !
     !  R Ray  21 August 2003
     !  Revised 9/15/08.
     !  Revised 3/28/11 - fixed Mt, MSt.
     !  Revised 6/19/20 - added degree-3 diurnal/semidl tides.
     !  Revised - from time to time, an overlooked constituent is added.

     implicit none

     integer,          intent(in)  ::  n
     double precision, intent(in)  ::  omega, p
     character(len=*), intent(in)  ::  constituent(n)
     double precision, intent(out) ::  f(n), u(n)

     double precision :: s,h
     double precision :: term1,term2,ftmp(4),utmp(4)
     character(len=8) :: ctmp(4)
     integer          :: i
     real             :: cosn,cos2n,sinn,sin2n,sin3n
     real             :: sinp,cosp,sin2p,cos2p
     double precision, parameter :: two=2.d0, three=3.d0
     double precision, parameter :: pi=3.141592654d0, rad=pi/180.D0

     ! Various trig functions of astronomical longitudes

     sinn = sin(omega*rad)
     cosn = cos(omega*rad)
     sin2n = sin(two*omega*rad)
     cos2n = cos(two*omega*rad)
     sin3n = sin(three*omega*rad)
     sinp  = sin(p*rad)
     cosp  = cos(p*rad)
     sin2p = sin(two*p*rad)
     cos2p = cos(two*p*rad)

     ! Compute standard nodal corrections f and u 

     do i=1,n

        select case ( constituent(i) )

           case ( 'O1' )
              term1 = 0.1886*sinn - 0.0058*sin2n - 0.0065*sin2p
              term2 = 1.0 + 0.1886*cosn - 0.0058*cos2n - 0.0065*cos2p
           case ( 'Q1','sig1' )
              term1 = 0.1886*sinn 
              term2 = 1.0 + 0.1886*cosn 
           case ( 'P1' )
              term1 = -0.0112*sinn 
              term2 = 1.0 - 0.0112*cosn 
           case ( 'K1' )
              term1 = -0.1554*sinn + 0.0031*sin2n
              term2 = 1.0 + 0.1158*cosn - 0.0028*cos2n
           case ( 'J1' )
              term1 = -0.227*sinn
              term2 = 1.0 + 0.169*cosn
           case ( 'OO1' )
              term1 = -0.640*sinn - 0.134*sin2n - 0.150*sin2p
              term2 = 1.0 + 0.640*cosn + 0.134*cos2n + 0.150*cos2p
           case ( 'M2','2N2','mu2','N2','MS4' )
              term1 = -0.03731*sinn + 0.00052*sin2n
              term2 = 1.0 - 0.03731*cosn + 0.00052*cos2n
           case ( 'K2' )
              term1 = -0.3108*sinn - 0.0324*sin2n
              term2 = 1.0 + 0.2853*cosn + 0.0324*cos2n
           case ( 'S2' )
              term1 = 0.00225*sinn
              term2 = 1.0 + 0.00225*cosn
           case ( 'M1p' )          ! Linear 3rd-deg terms
              term1 = -0.01815*sinn
              term2 = 1.0 - 0.27837*cosn
           case ( 'N2p' )
              term1 = 0.1705*sinn - 0.0035*sin2n - 0.0176*sin2p
              term2 = 1.0 + 0.1705*cosn - 0.0035*cos2n - 0.0176*cos2p
           case ( 'L2p' )
              term1 = -0.2495*sinn
              term2 = 1.0 + 0.1315*cosn
           case ( 'M3' )   ! Linear 3rd-deg terms
              term1 = -0.05644*sinn
              term2 = 1.0 - 0.05644*cosn
           case default
              term1 = 0.
              term2 = 1.

        end select

        f(i) = sqrt(term1**2 + term2**2)
        u(i) = atan( term1/term2 )/rad
        if (term1.ne.0.d0) cycle

        ! Following tide is compound & uses recursion

        select case ( constituent(i) )

           case ( 'M4' )
              ctmp(1)='M2'
              call fu_nodal6( omega, p, ctmp, 1, ftmp, utmp )
              f(i) = ftmp(1)**2
              u(i) = 2.0*utmp(1)

        end select

     end do

  end subroutine fu_nodal6

!==================================================================================

  subroutine perth5_tide( acosp, asinp, tide )

     ! Predict the tidal height for given values of acosp, asinp, and timemjd.
     ! timemjd is UT, in decimal MJD. --  e.g., 1 Jan 1990 noon = 47892.5
     ! From OLAM, timemjd = (s1900sim / 86400.) + 15020.

     implicit none
     real,    intent(in)  :: acosp(ntides), asinp(ntides)
     real,    intent(out) :: tide

     real(8) :: sum, x, haddx
     integer :: i

     sum = 0.d0
     do i = 1,ntides
        x = (args(i) + u(i))*pio180
        haddx = acosp(i)*f(i)*cos(x) + asinp(i)*f(i)*sin(x)
        sum = sum + haddx
     end do
     tide = 0.01 * sum  ! Converts from sum [cm] to tide [m]

  end subroutine perth5_tide

end module perth
