module mem_megan

  use CONSTS_MEGAN
  use EMIS_MAPS_MEGAN
  use consts_coms, only: r8, i1
  use leaf_coms,   only: nvtyp
  use rxns_data,   only: mechname

  implicit none

  private :: r8, i1

  real, allocatable :: pfts    (:,:)
  real, allocatable :: avg_temp(:)
  real, allocatable :: avg_ppfd(:)
  real, allocatable :: avg_ppfd_dif(:)
  real, allocatable :: lai_prev(:)
  real, allocatable :: lai_next(:)
  real, allocatable :: megan_emis(:,:)

  ! Length of the time between LAI updates (days)
  real, PARAMETER :: TSTLEN = 30.0

  integer :: ino
  integer :: n_scon_spc
  integer :: n_mech_spc

  real,          pointer, contiguous :: conv_fac(:) => null()
  integer,       pointer, contiguous :: spmh_map(:) => null()
  integer,       pointer, contiguous :: mech_map(:) => null()
  character(16), pointer, contiguous :: mech_spc(:) => null()
  real,          pointer, contiguous :: mech_mwt(:) => null()

  integer, allocatable :: mgn_2_gc_map(:)

  real :: norm_fact(nvtyp,16)

  logical :: has_pft_dataset = .false.

  real :: ugphr_molpsec(n_spca_spc)

contains


  subroutine alloc_megan(mland, mwa)
    use utilio_defn, only: index1
    use cgrid_spcs,  only: n_gc_emis

    implicit none

    integer, intent(in) :: mland
    integer, intent(in) :: mwa
    integer             :: iland

    allocate(avg_temp    (mland)) ; avg_temp = 0.0
    allocate(lai_prev    (mland)) ; lai_prev = 0.0
    allocate(lai_next    (mland)) ; lai_next = 0.0
    allocate(avg_ppfd    (mland)) ; avg_ppfd     = 0.0
    allocate(avg_ppfd_dif(mland)) ; avg_ppfd_dif = 0.0

    ! Store index of NO
    INO = INDEX1('NO', N_MGN_SPC, MGN_SPC)

    ! Map output mechanism

    if (index(mechname,'CB6') > 0) then

       n_scon_spc = n_cb6
       n_mech_spc = n_cb6_spc

       mech_spc => mech_spc_cb6
       spmh_map => spmh_map_cb6
       mech_map => mech_map_cb6
       conv_fac => conv_fac_cb6
       mech_mwt => mech_mwt_cb6

    elseif (index(mechname,'CB05') > 0) then

       n_scon_spc = n_cb05
       n_mech_spc = n_cb05_spc

       mech_spc => mech_spc_cb05
       spmh_map => spmh_map_cb05
       mech_map => mech_map_cb05
       conv_fac => conv_fac_cb05
       mech_mwt => mech_mwt_cb05

    elseif (index(mechname,'SAPRC') > 0) then

       n_scon_spc = n_saprc99
       n_mech_spc = n_saprc99_spc

       spmh_map => spmh_map_saprc99
       mech_map => mech_map_saprc99
       conv_fac => conv_fac_saprc99

    else
       write(*,*) "Invalid chemical mechanism specified in MEGAN."
       stop
    endif

    allocate(mgn_2_gc_map(n_gc_emis))

    allocate(megan_emis(n_mech_spc,mland))
    allocate(pfts            (0:16,mland))

    !$omp parallel do
    do iland = 1, mland
       megan_emis(:,iland) = 0.0
       pfts      (:,iland) = 0.0
    enddo
    !$omp end parallel do

  end subroutine alloc_megan


  subroutine filltab_megan()
    use var_tables, only: increment_vtable
    implicit none

    if (allocated(avg_temp)) then
       call increment_vtable('MEGAN_AVGTEMP', 'LW', rvar1=avg_temp)
    endif

    if (allocated(avg_ppfd)) then
       call increment_vtable('MEGAN_AVGPPFD', 'LW', rvar1=avg_ppfd)
    endif

    if (allocated(avg_ppfd_dif)) then
       call increment_vtable('MEGAN_AVGPPFD_DIF', 'LW', rvar1=avg_ppfd_dif)
    endif

  end subroutine filltab_megan


  integer function olam2pft( lclass, glat )
    implicit none

    integer, intent(in) :: lclass
    real,    intent(in) :: glat

    ! Converts the OLAM land surface types to one of the 16 CAM
    ! plant funcional types needed by the MEGAN emissions model:
    !
    !  1 Needleaf evergreen temperate tree
    !  2 Needleaf deciduous boreal tree
    !  3 Needleaf evergreen boreal tree
    !  4 Broadleaf evergreen tropical tree
    !  5 Broadleaf evergreen temperate tree
    !  6 Broadleaf deciduous tropical tree
    !  7 Broadleaf deciduous temperate tree
    !  8 Broadleaf deciduous boreal tree
    !  9 Broadleaf evergreen temperate shrub
    ! 10 Broadleaf deciduous temperate shrub
    ! 11 Broadleaf deciduous boreal shrub
    ! 12 Cold C3 grass
    ! 13 Cool C3 grass
    ! 14 Warm C3 grass
    ! 15 Corn
    ! 16 Other crops

    select case( lclass )

    case(0)  ! ocean -> unmapped

       olam2pft = 0

    case(1)  ! lakes -> unmapped

       olam2pft = 0

    case(2)  ! glacier -> unmapped

       olam2pft = 0

    case(3)  ! desert -> unmapped

       olam2pft = 0

    case(4)  ! Evergreen needleleaf tree

       if (abs(glat) > 60.0) then
          olam2pft = 3  ! Needleaf evergreen boreal tree
       else
          olam2pft = 1  ! Needleaf evergreen temperate tree
       endif

    case(5)  ! Deciduous needleleaf tree

       olam2pft = 2  ! Needleaf deciduous boreal tree

    case(6)  ! Deciduous broadleaf tree

       if (abs(glat) > 60.0) then
          olam2pft = 8  ! Broadleaf deciduous boreal tree
       elseif (abs(glat) > 30.0) then
          olam2pft = 7  ! Broadleaf deciduous temperate tree
       else
          olam2pft = 6  ! Broadleaf deciduous tropical tree
       endif

    case(7)  ! Evergreen broadleaf tree

       if (abs(glat) > 30.0) then
          olam2pft = 5  ! Broadleaf evergreen temperate tree
       else
          olam2pft = 4  ! Broadleaf evergreen tropical tree
       endif

    case(8:9)  ! Short and tall grass

       if (abs(glat) > 60.0) then
          olam2pft = 12  ! Cold C3 grass
       elseif (abs(glat) > 30.0) then
          olam2pft = 13  ! Cool C3 grass
       else
          olam2pft = 14  ! Warm C3 grass
       endif

    case (10)  ! semi-desert

       if (abs(glat) > 60.0) then
          olam2pft = 11 ! Broadleaf deciduous boreal shrub
       else
          olam2pft = 10 ! Broadleaf deciduous temperate shrub
       endif

    case(11)  ! tundra

       olam2pft = 12  ! Cold C3 grass

    case(12)  ! Evergreen shrub

       olam2pft = 9 ! Broadleaf evergreen temperate shrub

    case(13)  ! Deciduous shrub

       if (abs(glat) > 60.0) then
          olam2pft = 11  ! Broadleaf deciduous boreal shrub
       else
          olam2pft = 10  ! Broadleaf deciduous temperate shrub
       endif

    case(14)  ! Mixed woodland

       if (abs(glat) > 60.0) then
          olam2pft = 8  ! Broadleaf deciduous boreal tree
       elseif (abs(glat) > 30.0) then
          olam2pft = 7  ! Broadleaf deciduous temperate tree
       else
          olam2pft = 6  ! Broadleaf deciduous tropical tree
       endif

    case(15)  ! Crop

       olam2pft = 16 ! Other crops

    case(16)  ! Irrigated crop

       olam2pft = 16 ! Other crops

    case(17)  ! Bog or marsh

       if (abs(glat) > 60.0) then
          olam2pft = 12  ! Cold C3 grass
       elseif (abs(glat) > 30.0) then
          olam2pft = 13  ! Cool C3 grass
       else
          olam2pft = 14  ! Warm C3 grass
       endif

    case(18)  ! Wooded grassland

       if (abs(glat) > 60.0) then
          olam2pft = 8  ! Broadleaf deciduous boreal tree
       elseif (abs(glat) > 30.0) then
          olam2pft = 7  ! Broadleaf deciduous temperate tree
       else
          olam2pft = 6  ! Broadleaf deciduous tropical tree
       endif

    case(19)  ! Urban

       if (abs(glat) > 60.0) then
          olam2pft = 8  ! Broadleaf deciduous boreal tree
       elseif (abs(glat) > 30.0) then
          olam2pft = 7  ! Broadleaf deciduous temperate tree
       else
          olam2pft = 6  ! Broadleaf deciduous tropical tree
       endif

    case(20)   ! Wetland evergreen broadleaf

       if (abs(glat) > 30.0) then
          olam2pft = 5  ! Broadleaf evergreen temperate tree
       else
          olam2pft = 4  ! Broadleaf evergreen tropical tree
       endif

    case(21)   ! Very urban

       olam2pft = 0

    case default

       olam2pft = 0

    end select

  end function olam2pft



  subroutine read_pfts()

    use mem_land,    only: mland, omland
    use mem_sfcg,    only: sfcg
    use misc_coms,   only: io6
    use max_dims,    only: pathlen
    use oname_coms,  only: nl
    use analysis_lib,only: gdtost_ll
    use hdf5_utils,  only: shdf5_open, shdf5_irec, shdf5_info, shdf5_close

    implicit none

    integer            :: ndims, idims(4)
    character(pathlen) :: filename
    integer, parameter :: nx_pft = 7200
    integer, parameter :: ny_pft = 3600

    integer, parameter :: nio = nx_pft + 4
    integer, parameter :: njo = ny_pft + 4

    integer(i1), allocatable :: landmask(:,:)
    integer(i1), allocatable :: pfts_ll (:,:)

    real    :: dlon, dlat, swlat, swlon
    integer :: iproj

    integer :: ips, ip, iland, iwsfc
    logical :: exists
    real    :: grx, gry

    has_pft_dataset = .false.

    filename = nl%megan_pfts_file

    inquire(file=filename, exist=exists)
    if (.not. exists) then
       write(io6,*) "megan_init:  Cannot find plant functional types dataset."
    else
       has_pft_dataset = .true.
    endif

    if (has_pft_dataset) then
       write(io6,*) "Reading " // trim(filename)
       call shdf5_open(filename, 'R', trypario=.true.)

       call shdf5_info('LANDMASK', ndims, idims)
       if ( ndims /= 2 .or. all( idims(1:2) /= [nx_pft, ny_pft] ) ) then
          write(*,*) "Cannot find plant functional types dataset."
          has_pft_dataset = .false.
       endif

       call shdf5_info('PCT_PFT', ndims, idims)
       if ( ndims /= 3 .or. all( idims(1:3) /= [nx_pft, ny_pft, 17] ) ) then
          write(*,*) "Cannot find plant functional types dataset."
          has_pft_dataset = .false.
       endif

       if (.not. has_pft_dataset) call shdf5_close()
    endif

    if (has_pft_dataset) then

       dlon  =    0.05
       dlat  =    0.05
       swlat =  -89.975
       swlon = -179.975
       iproj =    1

       allocate(landmask(nx_pft,ny_pft))
       allocate(pfts_ll (nx_pft,ny_pft))

       ! First read land/sea mask
       call shdf5_irec(2, [nx_pft, ny_pft], 'LANDMASK', bvar2=landmask)

       ! Loop over all pfts in file. Skip 17 (other crop) which is 0 in file
       do ips = 1, 16

          ! pft 16 in file is crops; map to 'other crop' rather than corn
          if (ips==16) then
             ip = ips
          else
             ip = ips - 1
          endif

          ! Read PFT from file
          call shdf5_irec(ndims, idims, 'PCT_PFT', bvar2=pfts_ll, &
                          start=[1,1,ips], counts=[nx_pft,ny_pft,1] )

          where (landmask == 0)
             pfts_ll = -127_i1
          end where

          !$omp parallel do private(iwsfc)
          do iland = 2, mland
             iwsfc = iland + omland

             if (sfcg%glatw(iwsfc) < -70.) then

                ! Set Antarctica to all barren for now
                if (ip == 0) then
                   pfts(ip,iland) = 1.0
                else
                   pfts(ip,iland) = 0.0
                endif

             else

                call gdtost_ll(nx_pft, ny_pft, sfcg%glonw(iwsfc), sfcg%glatw(iwsfc), &
                               pfts(ip,iland), swlon, swlat, dlon, dlat, iproj, b2d=pfts_ll)

                pfts(ip,iland) = pfts(ip,iland) * 0.01

             endif
          enddo

       enddo

       call shdf5_close()

    endif

  end subroutine read_pfts




  subroutine megan_init()

    use mem_land,     only: land, mland, omland
    use mem_sfcg,     only: sfcg
    use leaf4_canopy, only: vegndvi
    use consts_coms,  only: pio180, pi2
    use misc_coms,    only: current_time, runtype
    use cgrid_spcs,   only: n_gc_emis, gc_emis
    use utilio_defn,  only: index1
    use oname_coms,   only: nl
    use max_dims,     only: pathlen
    use leaf_coms,    only: veg_frac

    implicit none

    real    :: ta, tb, tc, td, te
    integer :: iland, iwsfc, n, ip
    real    :: tot

    integer :: day
    real    :: hour, beta, sinbeta

    integer, external :: julday

    real, parameter :: ug2g   = 1.e-6       ! convert microgram to metric gram
    real, parameter :: hr2sec = 1./3600.    ! convert 1/hr to 1/second

    real, parameter :: solcon = 1370.0
    real, parameter :: tdiff0 = 5.0
    real, parameter :: hrmax  = 14.0

    real, parameter :: par_ratio = 0.4  ! typical PAR / SOLAR
    real, parameter :: ppfd_ratio = 4.25  ! typical PPFD / PAR

    day = julday( current_time%month, current_time%date, current_time%year )

    if (nl%megan_use_pfts_file .and. len_trim(nl%megan_pfts_file) > 0) then
       call read_pfts( )
    endif

    !$omp parallel do private (iwsfc,hour,tot,beta,sinbeta,ta,tb,tc,td,te)
    do iland = 2, mland
       iwsfc = iland + omland

! THIS ILAND/IWSFC LOOP IS OVER ALL MLAND POINTS, EVEN THOSE THAT ARE NOT
! PRIMARY ON THIS SUBDOMAIN.  IS THIS OK FOR THIS OPERATION?  IF NOT, AND
! IF LOOP MUST BE LIMITED TO PRIMARY IWSFC POINTS ONLY, IS MPI_SEND/MPI_RECV
! REQUIRED FOR ANY QUANTITIES AT A TIME OTHER THAN THE END OF EACH TIMESTEP
! (WHEN A ROUTINE SEND/RECV IS DONE)?

       hour = current_time%time / 3600.0_r8 + sfcg%glonw(iwsfc) / 15.0

       if ( hour < 0.0 ) then
          hour = hour + 24.0
       elseif ( hour > 24.0 ) then
          hour = hour - 24.0
       endif

       if (has_pft_dataset) then
          tot = sum(pfts(0:16,iland))
       else
          tot = 1.0
       endif

       if ( (.not. has_pft_dataset) .or. (tot < 0.2) ) then

          pfts(   0,iland) = 1.0
          pfts(1:16,iland) = 0.0

          ! Map OLAM land surface classes to what MEGAN expects
          ! Note: This need to be improved, should use original OLSON ecosystem classes
          ip = olam2pft( sfcg%leaf_class(iwsfc), sfcg%glatw(iwsfc) )

          if (ip > 0) then
             pfts(ip,iland) =       veg_frac( sfcg%leaf_class(iwsfc) )
             pfts( 0,iland) = 1.0 - veg_frac( sfcg%leaf_class(iwsfc) )
          endif

       else

          pfts(0:16,iland) = pfts(0:16,iland) / tot

       endif

       ! initialize average temperature based on current temperature modified
       ! by a simple diurnal variation, and the average radiation to be a
       ! fraction of midday clear-sky top-of-atmosphere solar flux

       if (runtype == 'INITIAL') then
          beta    = calcbeta( Day, sfcg%glatw(iwsfc), 12.0 )
          sinbeta = max( sin(beta * pio180), 0.0 )

          avg_temp(iland) = sfcg%cantemp(iwsfc) &
               - tdiff0 * cos( (hour - hrmax) * pi2 / 24.0 ) * sinbeta

          avg_ppfd    (iland) = solcon * sinbeta * 0.25 * par_ratio * ppfd_ratio
          avg_ppfd_dif(iland) = solcon * sinbeta * 0.05 * par_ratio * ppfd_ratio
       endif

       ! set current and next months LAI

       call vegndvi( iland, iwsfc, 0.0, sfcg%leaf_class(iwsfc), 1.0, te,   &
                     land%veg_ndvip(iland), land%veg_ndvif(iland), ta, tb, &
                     lai_prev(iland), tc, td                               )

       call vegndvi( iland, iwsfc, 1.0, sfcg%leaf_class(iwsfc), 1.0, te,   &
                     land%veg_ndvip(iland), land%veg_ndvif(iland), ta, tb, &
                     lai_next(iland), tc, td                               )

    enddo
    !$omp end parallel do

    do n = 1, n_gc_emis
       mgn_2_gc_map(n) = INDEX1( GC_EMIS(n), n_mech_spc, mech_spc )
    enddo

    norm_fact( :3,:) = 0.0
    norm_fact(21:,:) = 0.0

    norm_fact( 4,:) = (/ 0.442, 0.442, 0.442, 0.551, 0.484, 0.551, 0.484, 0.484, &
                         0.442, 0.442, 0.442, 0.405, 0.405, 0.405, 0.362, 0.362 /)
    norm_fact( 5,:) = (/ 0.445, 0.445, 0.445, 0.554, 0.487, 0.554, 0.487, 0.487, &
                         0.445, 0.445, 0.445, 0.408, 0.408, 0.408, 0.365, 0.365 /)
    norm_fact( 6,:) = (/ 0.454, 0.454, 0.454, 0.566, 0.498, 0.566, 0.498, 0.498, &
                         0.454, 0.454, 0.454, 0.417, 0.417, 0.417, 0.373, 0.373 /)
    norm_fact( 7,:) = (/ 0.463, 0.463, 0.463, 0.577, 0.508, 0.577, 0.508, 0.508, &
                         0.463, 0.463, 0.463, 0.425, 0.425, 0.425, 0.380, 0.380 /)
    norm_fact( 8,:) = (/ 0.501, 0.501, 0.501, 0.624, 0.549, 0.624, 0.549, 0.549, &
                         0.501, 0.501, 0.501, 0.459, 0.459, 0.459, 0.411, 0.411 /)
    norm_fact( 9,:) = (/ 0.547, 0.547, 0.547, 0.682, 0.600, 0.682, 0.600, 0.600, &
                         0.547, 0.547, 0.547, 0.502, 0.502, 0.502, 0.449, 0.449 /)
    norm_fact(10,:) = (/ 0.469, 0.469, 0.469, 0.584, 0.514, 0.584, 0.514, 0.514, &
                         0.469, 0.469, 0.469, 0.430, 0.430, 0.430, 0.385, 0.385 /)
    norm_fact(11,:) = (/ 0.545, 0.545, 0.545, 0.679, 0.597, 0.679, 0.597, 0.597, &
                         0.545, 0.545, 0.545, 0.499, 0.499, 0.499, 0.447, 0.447 /)
    norm_fact(12,:) = (/ 0.346, 0.346, 0.346, 0.431, 0.379, 0.431, 0.379, 0.379, &
                         0.346, 0.346, 0.346, 0.317, 0.317, 0.317, 0.284, 0.284 /)
    norm_fact(13,:) = (/ 0.358, 0.358, 0.358, 0.446, 0.392, 0.446, 0.392, 0.392, &
                         0.358, 0.358, 0.358, 0.328, 0.328, 0.328, 0.294, 0.294 /)
    norm_fact(14,:) = (/ 0.448, 0.448, 0.448, 0.559, 0.491, 0.559, 0.491, 0.491, &
                         0.448, 0.448, 0.448, 0.411, 0.411, 0.411, 0.368, 0.368 /)
    norm_fact(15,:) = (/ 0.529, 0.529, 0.529, 0.660, 0.580, 0.660, 0.580, 0.580, &
                         0.529, 0.529, 0.529, 0.485, 0.485, 0.485, 0.434, 0.434 /)
    norm_fact(16,:) = (/ 0.353, 0.353, 0.353, 0.441, 0.387, 0.441, 0.387, 0.387, &
                         0.353, 0.353, 0.353, 0.324, 0.324, 0.324, 0.290, 0.290 /)
    norm_fact(17,:) = (/ 0.358, 0.358, 0.358, 0.446, 0.392, 0.446, 0.392, 0.392, &
                         0.358, 0.358, 0.358, 0.328, 0.328, 0.328, 0.294, 0.294 /)
    norm_fact(18,:) = (/ 0.595, 0.595, 0.595, 0.741, 0.652, 0.741, 0.652, 0.652, &
                         0.595, 0.595, 0.595, 0.545, 0.545, 0.545, 0.488, 0.488 /)
    norm_fact(19,:) = (/ 0.413, 0.413, 0.413, 0.515, 0.453, 0.515, 0.453, 0.453, &
                         0.413, 0.413, 0.413, 0.379, 0.379, 0.379, 0.339, 0.339 /)
    norm_fact(20,:) = (/ 0.456, 0.456, 0.456, 0.568, 0.499, 0.568, 0.499, 0.499, &
                         0.456, 0.456, 0.456, 0.418, 0.418, 0.418, 0.374, 0.374 /)

    do n = 1, n_spca_spc
       ugphr_molpsec(n) = ug2g / spca_mwt(n) * hr2sec
    enddo

  end subroutine megan_init


  subroutine megan_store_lai()

    use mem_land,     only: land, mland, omland
    use mem_sfcg,      only: sfcg
    use leaf4_canopy, only: vegndvi

    implicit none

    real    :: ta, tb, tc, td, te
    integer :: iland, iwsfc

    !$omp parallel do private (iwsfc, ta, tb, tc, td, te)
    do iland = 2, mland
       iwsfc = iland + omland

! THIS ILAND/IWSFC LOOP IS OVER ALL MLAND POINTS, EVEN THOSE THAT ARE NOT
! PRIMARY ON THIS SUBDOMAIN.  IS THIS OK FOR THIS OPERATION?  IF NOT, AND
! IF LOOP MUST BE LIMITED TO PRIMARY IWSFC POINTS ONLY, IS MPI_SEND/MPI_RECV
! REQUIRED FOR ANY QUANTITIES AT A TIME OTHER THAN THE END OF EACH TIMESTEP
! (WHEN A ROUTINE SEND/RECV IS DONE)?

       lai_prev(iland) = lai_next(iland)

       call vegndvi( iland, iwsfc, 1.0, sfcg%leaf_class(iwsfc), 1.0,   &
                     te, land%veg_ndvip(iland), land%veg_ndvif(iland), &
                     ta, tb, lai_next(iland), tc, td                   )

       call vegndvi( iland, iwsfc, 1.0, sfcg%leaf_class(iwsfc), land%veg_height(iland), &
                     te, land%veg_ndvip(iland), land%veg_ndvif(iland), &
                     ta, tb, lai_next(iland), tc, td                   )
    enddo
    !$omp end parallel do

  end subroutine megan_store_lai


  subroutine megan_avg_temp()

    use leaf_coms,   only: dt_leaf
    use mem_sfcg,    only: sfcg
    use mem_land,    only: land, mland, omland
    use consts_coms, only: r8

    implicit none

    real, parameter :: tscali = 2.0 / 86400.0

    real    :: si, so
    integer :: iland, iwsfc

    si = dt_leaf * tscali
    so = 1.0 - si

    !$omp parallel do private (iwsfc)
    do iland = 2, mland
       iwsfc = iland + omland

! THIS ILAND/IWSFC LOOP IS OVER ALL MLAND POINTS, EVEN THOSE THAT ARE NOT
! PRIMARY ON THIS SUBDOMAIN.  IS THIS OK FOR THIS OPERATION?  IF NOT, AND
! IF LOOP MUST BE LIMITED TO PRIMARY IWSFC POINTS ONLY, IS MPI_SEND/MPI_RECV
! REQUIRED FOR ANY QUANTITIES AT A TIME OTHER THAN THE END OF EACH TIMESTEP
! (WHEN A ROUTINE SEND/RECV IS DONE)?

       avg_temp    (iland) = avg_temp    (iland) * so + sfcg%cantemp     (iwsfc) * si
       avg_ppfd    (iland) = avg_ppfd    (iland) * so + land%ppfd        (iland) * si
       avg_ppfd_dif(iland) = avg_ppfd_dif(iland) * so + land%ppfd_diffuse(iland) * si
    enddo
    !$omp end parallel do

  end subroutine megan_avg_temp


  subroutine megan_driver()

    use mem_land,    only: mland

    implicit none

    integer             :: iland

    !$omp parallel do
    do iland = 2, mland

! THIS ILAND/IWSFC LOOP IS OVER ALL MLAND POINTS, EVEN THOSE THAT ARE NOT
! PRIMARY ON THIS SUBDOMAIN.  IS THIS OK FOR THIS OPERATION?  IF NOT, AND
! IF LOOP MUST BE LIMITED TO PRIMARY IWSFC POINTS ONLY, IS MPI_SEND/MPI_RECV
! REQUIRED FOR ANY QUANTITIES AT A TIME OTHER THAN THE END OF EACH TIMESTEP
! (WHEN A ROUTINE SEND/RECV IS DONE)?

       call megan_driver1(iland)
    enddo
    !$omp end parallel do

  end subroutine megan_driver


  subroutine megan_driver1(iland)

    use leaf_coms,   only: nvtyp, veg_frac, tai_max
    use mem_land,    only: omland, land
    use misc_coms,   only: current_time
    use mem_sfcg,    only: sfcg
    use consts_coms, only: pio180

    implicit none

    integer, intent(in) :: iland

    real :: gam_lht             ! LAI activity factor
    real :: gam_age(n_mgn_spc)  ! leaf age activity factor
    real :: gam_tld(n_mgn_spc)  ! light-dependent activity factor
    real :: gam_tli(n_mgn_spc)  ! light-independent activity factor
    real :: gam_smt(n_mgn_spc)  ! soil moisture activity factor

    real :: gam_tli_ip(n_mgn_spc)
    real :: gam_tld_ip(n_mgn_spc)

    real :: em(n_mgn_spc)
    real :: tmper(n_spca_spc)
    real :: outer(n_mech_spc)

    integer :: s, iwsfc, ip, lc

    real :: glat, glon, hour
    real :: LAIc, slai
    real :: beta, sinbeta

    integer :: nmpmg, nmpsp, nmpmc
    integer :: day, jday, jdate

    integer, external :: julday

    real :: ppfd0, ppfd0_dif, ppfd24, ppfd24_dif, temp24, vegtk, ptot, laimax

    jday  = julday( current_time%month, current_time%date, current_time%year )
    jdate = current_time%year*1000 + jday

    iwsfc = iland + omland

    slai = land%veg_lai(iland) * (1.0 - land%snowfac(iland))
    ptot = 1.0 - pfts(0,iland)

    if (slai > 0.01 .and. land%veg_fracarea(iland) > 0.01 .and. ptot > 0.01) then

       glat = sfcg%glatw(iwsfc)
       glon = sfcg%glonw(iwsfc)
       lc   = sfcg%leaf_class(iwsfc)

       ppfd0      = land%ppfd(iland)
       ppfd0_dif  = land%ppfd_diffuse(iland)

       ppfd24     = avg_ppfd(iland)
       ppfd24_dif = avg_ppfd_dif(iland)

       temp24 = avg_temp(iland)
       vegtk  = land%veg_temp(iland)

       laimax = tai_max(lc) * land%veg_lai(iland) / ( land%veg_tai(iland) * veg_frac(lc) )
       LAIc   = min( land%veg_lai(iland) / max(ptot, 0.1), laimax ) * (1.0 - land%snowfac(iland))

       ! compute local solar day/hour

       day  = jday
       hour = current_time%time / 3600.0_r8 + glon / 15.0

       if ( hour < 0.0 ) then
          hour = hour + 24.0
       !  day  = day  -  1
       elseif ( hour > 24.0 ) then
          hour = hour - 24.0
       !  day  = day  +  1
       endif

       ! solar zenith angle
       ! TODO: take into account surface slope?

       ! Sinbeta = cosz(iw)
       ! Sinbeta = land%cosz(iland)
       Beta    = Calcbeta( Day, glat, Hour )
       Sinbeta = SIN(Beta * pio180)

       call gamma_lai( gam_lht, laic )

       call gamma_sm( gam_smt,                  &
                      sfcg%leaf_class(  iwsfc), &
                      land%soil_water(:,iland), &
                      land%wresid_vg (:,iland), &
                      land%wsat_vg   (:,iland), &
                      land%alpha_vg  (:,iland), &
                      land%en_vg     (:,iland), &
                      land%wfrac_low (:,iland)  )

       call gamma_a( gam_age, lai_prev(iland), lai_next(iland), temp24 )

       gam_tld(:) = 0.
       gam_tli(:) = 0.

       do ip = 1, 16
          if (pfts(ip,iland) < 0.01) cycle

          call gamma_ce( gam_tld_ip, gam_tli_ip, sinbeta, temp24, vegtk, ppfd0, &
               ppfd0_dif, ppfd24, ppfd24_dif, ip, laic, norm_fact(sfcg%leaf_class(iwsfc),ip) )

          do s = 1, n_mgn_spc
             gam_tld(s) = gam_tld(s) + pfts(ip,iland) * gam_tld_ip(s)
             gam_tli(s) = gam_tli(s) + pfts(ip,iland) * gam_tli_ip(s)
          enddo
       enddo

       do s = 1, n_mgn_spc
          gam_tld(s) = gam_tld(s) / ptot
          gam_tli(s) = gam_tli(s) / ptot

          em(s) = GAM_AGE(s) * GAM_SMT(s) * &
                  ( (1.0-LDF_FCT(S)) * GAM_TLI(s) * GAM_LHT + LDF_FCT(S) * GAM_TLD(s) )
       enddo

       ! We do soil NOx elsewhere, so turn off NOx here

       em(INO) = 0.0

       ! Conversion from MGN 20 to speciated 150

       do s = 1, n_smap_spc
          nmpmg    = mg20_map(s)
          tmper(s) = em(nmpmg) &
                   * sum(ef_all(1:16,nmpmg) * effs_all(1:16,s) * pfts(1:16,iland))
       enddo

       ! Convert from ug/m^2/hr to mol/m^2/sec using their MW

       do s = 1, n_spca_spc
          tmper(s) = tmper(s) * ugphr_molpsec(s)
       enddo

       ! Conversion from speciated species to MECHANISM species,
       ! with units of mol/sec

       outer = 0.0

       do s = 1, n_scon_spc
          nmpsp = spmh_map(s)         ! Mapping value for SPCA
          nmpmc = mech_map(s)         ! Mapping value for MECHANISM
          outer(nmpmc) = outer(nmpmc) + tmper(nmpsp) * conv_fac(s)
       enddo

       ! Convert to mol/sec
       do s = 1, n_mech_spc
          megan_emis(s,iland) = outer(s) * sfcg%area(iwsfc)
       enddo

    else

       do s = 1, n_mech_spc
          megan_emis(s,iland) = 0.0
       enddo

    endif

  end subroutine megan_driver1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Functions to calculate GAMMA_P, GAMMA_T, GAMMA_L, GAMMA_A for BVOCs
!
!     Scientific algorithm
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensino)
!                 Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.EXT)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_age] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1 (11/27/06)
!
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = PPFD emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!        Derivation:
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }
!           where LDF = light dependent function (non-dimension)
!
!     For ISOPRENE
!                 Assumption: LDF = 1 for isoprene            (11/27/06)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!
!     For NON-ISOPRENE
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
!                        { (1-LDF) + [LDF][GAMMA_P] }
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------------------------------------------------------------------
!.....1) Calculate GAM_L (GAMMA_LAI)
!-----------------------------------------------------------------------
!                            0.49[LAI]
!             GAMMA_LAI = ----------------    (non-dimension)
!                         (1+0.2LAI^2)^0.5
!
!     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
!-----------------------------------------------------------------------

  SUBROUTINE GAMMA_LAI( gam_l, lai )
    IMPLICIT NONE

    REAL, INTENT(IN)  :: lai
    REAL, INTENT(OUT) :: gam_l

    gam_l = 0.49 * lai / sqrt(1.0 + 0.2 * lai**2)

  END SUBROUTINE GAMMA_LAI

!-----------------------------------------------------------------------
!.....5) Calculate GAM_A (GAMMA_age)
!-----------------------------------------------------------------------
!
!             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
!           where Fnew = new foliage fraction
!                 Fgro = growing foliage fraction
!                 Fmat = mature foliage fraction
!                 Fold = old foliage fraction
!                 Anew = relative emission activity for new foliage
!                 Agro = relative emission activity for growing foliage
!                 Amat = relative emission activity for mature foliage
!                 Aold = relative emission activity for old foliage
!
!
!             For foliage fraction
!             Case 1) LAIc = LAIp
!             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
!
!             Case 2) LAIp > LAIc
!             Fnew = 0.0  , Fgro = 0.0
!             Fmat = 1-Fold
!             Fold = (LAIp-LAIc)/LAIp
!
!             Case 3) LAIp < LAIc
!             Fnew = 1-(LAIp/LAIc)                       t <= ti
!                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
!
!             Fmat = LAIp/LAIc                           t <= tm
!                  = (LAIp/LAIc) +
!                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
!
!             Fgro = 1 - Fnew - Fmat
!             Fold = 0.0
!
!           where
!             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
!                = 2.9                                  Tt >  303
!             tm = 2.3*ti
!
!             t  = length of the time step (days)
!             ti = number of days between budbreak and the induction of
!                  emission
!             tm = number of days between budbreak and the initiation of
!                  peak emissions rates
!             Tt = average temperature (K) near top of the canopy during
!                  current time period (daily ave temp for this case)
!
!
!             For relative emission activity
!             Case 1) Constant
!             Anew = 1.0  , Agro = 1.0  , Amat = 1.0  , Aold = 1.0
!
!             Case 2) Monoterpenes
!             Anew = 2.0  , Agro = 1.8  , Amat = 0.95 , Aold = 1.0
!
!             Case 3) Sesquiterpenes
!             Anew = 0.4  , Agro = 0.6  , Amat = 1.075, Aold = 1.0
!
!             Case 4) Methanol
!             Anew = 3.0  , Agro = 2.6  , Amat = 0.85 , Aold = 1.0
!
!             Case 5) Isoprene
!             Anew = 0.05 , Agro = 0.6  , Amat = 1.125, Aold = 1.0
!
!     SUBROUTINE GAMMA_A returns GAMMA_A
!-----------------------------------------------------------------------
  SUBROUTINE GAMMA_A( GAM_A, LAIp, LAIc, Tt )

    IMPLICIT NONE

    ! input
    real,    intent(in)  :: Tt           ! daily average temperature (K)
    real,    intent(in)  :: LAIp         ! previous LAI
    real,    intent(in)  :: LAIc         ! current (next) LAI

    ! output
    REAL,    INTENT(OUT) :: GAM_A(n_mgn_spc)

    ! Local variables

    INTEGER :: AINDX        ! relative emission activity index
    integer :: s            ! species loop index
    REAL    :: Fnew, Fgro, Fmat, Fold
    REAL    :: ti, tm       ! number of days between budbreak
                            ! and induction of emission,
                            ! initiation of peak emissions rates

!... Calculate foliage fraction

    if (LAIp < LAIc) then

!      Calculate ti and tm
       if (Tt .LE. 303.0) then
          ti = 5.0 + 0.7*(300.0-Tt)
       else
          ti = 2.9
       endif
       tm = 2.3*ti

!      Calculate Fnew and Fmat, then Fgro and Fold

       if (ti .GE. tstlen) then
          Fnew = 1.0 - (LAIp/LAIc)
       ELSE
          Fnew = (ti/tstlen) * ( 1.0-(LAIp/LAIc) )
       endif

       if (tm .ge. tstlen) then
          fmat = laip/laic
       else
          fmat = (laip/laic) + ( (tstlen-tm)/tstlen ) * ( 1.0-(laip/laic) )
       endif

       Fgro = 1.0 - Fnew - Fmat
       Fold = 0.0

    ELSEIF (LAIp == LAIc) then

       Fnew = 0.0
       Fgro = 0.1
       Fmat = 0.8
       Fold = 0.1

    ELSE ! (LAIp > LAIc)

       Fnew = 0.0
       Fgro = 0.0
       Fold = ( LAIp-LAIc ) / LAIp
       Fmat = 1.0 - Fold

    ENDIF

    do s = 1, n_mgn_spc

!...  Choose relative emission activity
!--------code by Xuemei Wang 11/04/2007----------------
       AINDX = REA_INDEX(S)
!---------------------------------------------------

!...Calculate GAMMA_A
       GAM_A(s) = Fnew * Anew(AINDX) + Fgro * Agro(AINDX) &
                + Fmat * Amat(AINDX) + Fold * Aold(AINDX)

    enddo

  END SUBROUTINE GAMMA_A

!-----------------------------------------------------------------------
!.....6) Calculate GAM_SMT (GAMMA_SM)
!-----------------------------------------------------------------------
!
!             GAMMA_SM =     1.0   (non-dimension)
!
!
!     SUBROUTINE GAMMA_S returns the GAMMA_SM values
!-----------------------------------------------------------------------
  SUBROUTINE GAMMA_SM( GAM_S, leaf_class, soil_water, wresid_vg, wsat_vg, &
                       alpha_vg, en_vg, wfrac_low )

    use leaf_coms,  only: kroot
    use mem_land,   only: nzg
    use leaf4_soil, only: soil_pot2wat

    IMPLICIT NONE

    REAL,    intent(out) :: GAM_S(n_mgn_spc)
    integer, intent(in)  :: leaf_class
    real,    intent(in)  :: soil_water(nzg)
    real,    intent(in)  :: wresid_vg (nzg)
    real,    intent(in)  :: wsat_vg   (nzg)
    real,    intent(in)  :: alpha_vg  (nzg)
    real,    intent(in)  :: en_vg     (nzg)
    real,    intent(in)  :: wfrac_low (nzg)

    real, parameter :: dsw = 0.06
    real            :: smk, soilwilt
    integer         :: k

    gam_s = 1.0

    do k = kroot(leaf_class), nzg

 ! Now, wilting point needs to be calculated for each grid cell
 ! since clearly defined soil textural class is no longer used.
 ! Temporary new code follows here, and orginal code is commented
 ! out afterward.

! TEMPORARY NEW CODE:

     ! Diagnose wilting point using psi of -150 m as definition
     ! threshold value

       call soil_pot2wat(-150., wresid_vg(k), wsat_vg(k), &
                         alpha_vg(k), en_vg(k), wfrac_low(k), soilwilt)

       if (soil_water(k) <= soilwilt) then
          smk = 0.0
       elseif (soil_water(k) >= soilwilt + dsw) then
          smk = 1.0
       else
          smk = (soil_water(k) - soilwilt) / dsw
       endif


! ORIGINAL CODE:

  !     nts = ntext_soil(k)

  !     if (soil_water(k) <= soilwilt(nts)) then
  !        smk = 0.0
  !     elseif (soil_water(k) >= soilwilt(nts) + dsw) then
  !        smk = 1.0
  !     else
  !        smk = (soil_water(k) - soilwilt(nts)) / dsw
  !     endif

! END ORIGINAL CODE

       gam_s(1) = max(1.0, smk)
    enddo

  END SUBROUTINE GAMMA_SM


  subroutine gamma_ce( gamma_tld, gamma_tli, sinbeta, temp24, veg_temp, &
                       ppfd, ppfd_dif, ppfd24, ppfd24_dif, cantype, lai, cnorm )

    implicit none

    ! input
    integer, intent(in) :: cantype
    real,    intent(in) :: sinbeta, lai
    real,    intent(in) :: temp24, veg_temp
    real,    intent(in) :: ppfd, ppfd_dif, ppfd24, ppfd24_dif
    real,    intent(in) :: cnorm

    ! output
    real, intent(out) :: gamma_tld (n_mgn_spc)
    real, intent(out) :: gamma_tli (n_mgn_spc)

    real :: kb, cluster, laidepth
    real :: sunfrac, sunppfd, shadeppfd
    real :: gamma_p_sun, gamma_p_shade
    real :: topt, eopt
    real :: gamma_t

    integer :: s

    if (ppfd24 > 0.1 .and. ppfd > 0.1) then

       Cluster = Canopychar(9,Cantype)

       ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
       Kb = Cluster * 0.5 / max(sinbeta, 1.e-4)

       ! LAI depth at this layer
       LAIdepth   = 0.5 * LAI

       ! fraction of leaves that are sunlit
       sunfrac = EXP(-Kb * LAIdepth)

       sunppfd   = ppfd
       shadeppfd = ppfd_dif

       gamma_p_sun   = Ea1p99(sunPPFD,   PPFD24,     pstd_sun)   *        sunfrac
       gamma_p_shade = Ea1p99(shadePPFD, PPFD24_dif, pstd_shade) * (1.0 - sunfrac)

       topt = 312.5 + 0.6 * (temp24 - 297.0)
       eopt = EXP(0.1 * (temp24 - 297.))

       do s = 1, n_mgn_spc
          gamma_t      = Ea1t99(veg_temp, topt, eopt, S)
          gamma_tld(s) = ( gamma_p_sun + gamma_p_shade ) * gamma_t * cnorm * lai
       enddo

    else

       do s = 1, n_mgn_spc
          gamma_tld(s) = 0.0
       enddo

    endif

    do s = 1, n_mgn_spc
       gamma_tli(s) = Ealti99(S, veg_temp)
    enddo

  end subroutine gamma_ce


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Ea1t99
!
!   Temperature dependence activity factor for emission type 1
!          (e.g. isoprene, MBO)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  REAL FUNCTION Ea1t99(T1, TOPT, EOPT, SPCNUM)

    IMPLICIT NONE

    INTEGER, intent(in) :: SPCNUM
    real,    intent(in) :: T1, TOPT, EOPT

    REAL, PARAMETER :: Ctm2 = 230.0
    REAL            :: X

    IF (T1 < 260.) THEN

       Ea1t99 = 0.0

    ELSE

       X = (T1 - Topt) / (Topt * T1 * 0.00831)

       Ea1t99 = CLeo(SPCNUM) * Eopt * Ctm2 * Exp(Ctm1(SPCNUM) * X)  &
              / (Ctm2 - Ctm1(SPCNUM) * (1. - EXP(Ctm2 * X)))

    ENDIF

  END FUNCTION  Ea1t99


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Ea1pp
!
! pstd = 200 for sun leaves and 50 for shade leaves
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  REAL FUNCTION Ea1p99(PPFD1, PPFD24, PSTD)

    IMPLICIT NONE

    real, intent(in) :: PPFD1, PPFD24, PSTD
    REAL             :: Alpha, C1

    IF (PPFD24 < 0.01) THEN
       Ea1p99 = 0.
    ELSE
       Alpha  = 0.004 - 0.0005 * LOG(PPFD24)
       C1     = 0.0468 * EXP(0.0005 * (PPFD24 - PSTD)) * (PPFD24 ** 0.6)
       Ea1p99 = (Alpha * C1 * PPFD1) /  &
                sqrt(1. + Alpha**2 * PPFD1**2)
    ENDIF

  END FUNCTION  Ea1p99


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Ealti99
!
! calculate light indepent algorithms
! coded by Xuemei Wang 05 Nov. 2007
!--   GAMMA_TLI =  exp[BETA*(T-Ts)]
!           where BETA   = temperature dependent parameter
!                 Ts     = standard temperature (normally 303K, 30C)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  REAL FUNCTION Ealti99(spcnum, temp)

    IMPLICIT NONE

    INTEGER, intent(in) :: spcnum
    real,    intent(in) :: temp

    REAL, PARAMETER :: Ts = 303.15

    Ealti99 = exp( tdf_prm(spcnum)*(temp-Ts) )

  END FUNCTION Ealti99


!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Calcbeta
!   Calculates the solar zenith angle
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  REAL FUNCTION Calcbeta(Day, Lat, Hour)

    IMPLICIT NONE

    integer, intent(in) :: Day
    real,    intent(in) :: lat, hour

    real :: SinDelta,  CosDelta, A, B, Sinbeta
    real, parameter :: Pi = 3.141592654, Rpi180 = 57.29578
    real, parameter :: twopi = 2.0 * pi
    real, parameter :: Rpi180i = 1.0 / Rpi180

    SinDelta = -SIN(0.40907) * COS(twopi * (Day + 10.0) / 365.0)
    CosDelta = sqrt(1.0 - SinDelta**2)

    A = SIN(Lat * Rpi180i) * SinDelta
    B = COS(Lat * Rpi180i) * Cosdelta
    Sinbeta = A + B * COS(twopi * (Hour - 12.0) / 24.0)
    Calcbeta = ASIN(Sinbeta) * Rpi180

  END FUNCTION Calcbeta


end module mem_megan
