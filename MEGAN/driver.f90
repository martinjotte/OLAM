module mem_megan

  use consts_coms, only: r8
  implicit none

  include 'EACO.EXT'
  include 'SPC_NOCONVER.EXT'
  include 'EFS_PFT.EXT'

  include 'SPC_CB05.EXT'
  include 'MAP_CV2CB05.EXT'

 include 'SPC_SAPRC99.EXT'
 include 'MAP_CV2SAPRC99.EXT'

  private :: r8

  integer, allocatable         :: ipft    (:)
  real,    allocatable, target :: avg_temp(:)
  real,    allocatable, target :: avg_srad(:)
  real,    allocatable         :: lai_prev(:)
  real,    allocatable         :: lai_next(:)

  real,    allocatable         :: megan_emis(:,:)

  ! Length of the time between LAI updates (days)
  real, PARAMETER :: TSTLEN = 30.0

  ! hardwire machanism for now
  character(16), parameter :: mechanism = 'CB05'

  integer, save :: ino
  integer, save :: n_scon_spc
  integer, save :: n_mech_spc

  real,          pointer, save :: conv_fac(:) => null()
  integer,       pointer, save :: spmh_map(:) => null()
  integer,       pointer, save :: mech_map(:) => null()
  character(16), pointer, save :: mech_spc(:) => null()

  integer, allocatable :: mgn_2_gc_map(:)
! integer, allocatable :: mgn_2_nr_map(:)

!  is this necessary?
!  real(r8), parameter :: meganfrq = 3600.0_r8
  
contains

  subroutine alloc_megan(mwl, mwa)
    use utilio_defn, only: index1
    use cgrid_spcs,  only: n_gc_emis, n_nr_emis
    implicit none

    integer, intent(in) :: mwl
    integer, intent(in) :: mwa

    allocate(avg_temp(mwl)) ; avg_temp = 0.0
    allocate(ipft    (mwl)) ; ipft     = 0
    allocate(lai_prev(mwl)) ; lai_prev = 0.0
    allocate(lai_next(mwl)) ; lai_next = 0.0

    allocate(avg_srad(mwa)) ; avg_srad = 0.0

    ! Store index of NO
    INO   = INDEX1('NO',  N_MGN_SPC,MGN_SPC)

    ! Map output mechanism

    if (mechanism == 'CB05') then

       n_scon_spc = n_cb05
       n_mech_spc = n_cb05_spc

       mech_spc => mech_spc_cb05
       spmh_map => spmh_map_cb05
       mech_map => mech_map_cb05
       conv_fac => conv_fac_cb05

    elseif (mechanism == 'SAPRC99') then

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
!   allocate(mgn_2_nr_map(n_nr_emis))

    allocate(megan_emis(mwl,n_mech_spc))

  end subroutine alloc_megan


  subroutine filltab_megan()
    use var_tables, only: vtab_r, num_var, increment_vtable
    implicit none

    if (allocated(avg_temp)) then
       call increment_vtable('MEGAN_AVGTEMP', 'LW')
       vtab_r(num_var)%rvar1_p => avg_temp
    endif

    if (allocated(avg_srad)) then
       call increment_vtable('MEGAN_AVGSRAD', 'AW')
       vtab_r(num_var)%rvar1_p => avg_srad
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


  subroutine megan_init()

    use leaf_coms,    only: mwl
    use mem_leaf,     only: land
    use leaf4_canopy, only: vegndvi
    use consts_coms,  only: pio180, pi2
    use misc_coms,    only: io6, current_time, runtype
    use mem_grid,     only: mwa, glatw
    use cgrid_spcs,   only: n_gc_emis, gc_emis, n_nr_emis, nr_emis
    use utilio_defn,  only: index1
    implicit none

    real    :: ta, tb, tc, td, te
    integer :: iwl, iw, n

    integer :: day
    real    :: hour, beta, sinbeta_max

    real,    external :: calcbeta
    integer, external :: julday

    real, parameter :: solcon = 1370.0
    real, parameter :: tdiff0 = 5.0
    real, parameter :: hrmax  = 14.0

    day = julday( current_time%month, current_time%date, current_time%year )

    !$omp parallel 
    !$omp do private (hour,beta,sinbeta_max,ta,tb,tc,td,te)
    do iwl = 2, mwl

       hour = current_time%time / 3600.0_r8 + land%glonw(iwl) / 15.0

       if ( hour < 0.0 ) then
          hour = hour + 24.0
       elseif ( hour > 24.0 ) then
          hour = hour - 24.0
       endif

       ! map OLAM land surface classes to what MEGAN expects

       ipft(iwl) = olam2pft( land%leaf_class(iwl), land%glatw(iwl) )

       ! initialize average temperature based on current temperature modified
       ! by a simple diurnal variation

       if (runtype == 'INITIAL') then
          beta        = calcbeta( Day, land%glatw(iwl), 12.0 )
          sinbeta_max = max( sin(beta * pio180), 0.0 )
          avg_temp(iwl) = land%cantemp(iwl) &
               - tdiff0 * cos( (hour - hrmax) * pi2 / 24.0 ) * sinbeta_max
       endif

       ! set current and next months LAI

       call vegndvi( iwl, land%leaf_class(iwl), 0.0, 1.0,                  &
                     land%veg_ndvip(iwl), land%veg_ndvif(iwl), ta, tb,     &
                     lai_prev(iwl), tc, td, te                             )

       call vegndvi( iwl, land%leaf_class(iwl), 1.0, 1.0,                  &
                     land%veg_ndvip(iwl), land%veg_ndvif(iwl), ta, tb,     &
                     lai_next(iwl), tc, td, te                             )

    enddo
    !$omp end do

    if (runtype == 'INITIAL') then

       !$omp do private (beta,sinbeta_max)
       do iw = 2, mwa

          ! initialize average radiation to be a fraction of midday clear-sky 
          ! top-of-atmosphere solar flux

          beta         = calcbeta( Day, glatw(iw), 12.0 )
          sinbeta_max  = max( sin(beta * pio180), 0.0 )
          avg_srad(iw) = solcon * sinbeta_max * 0.25
       enddo
       !$omp end do

    endif
    !$omp end parallel

    do n = 1, n_gc_emis
       mgn_2_gc_map(n) = INDEX1( GC_EMIS(n), n_mech_spc, mech_spc )
    enddo

!   do n = 1, n_nr_emis
!      mgn_2_nr_map(n) = INDEX1( NR_EMIS(n), n_mech_spc, mech_spc )
!   enddo

  end subroutine megan_init


  subroutine megan_store_lai()

    use leaf_coms,    only: mwl
    use mem_leaf,     only: land
    use leaf4_canopy, only: vegndvi
    
    implicit none

    real    :: ta, tb, tc, td, te
    integer :: iwl

    !$omp parallel do
    do iwl = 2, mwl

       lai_prev(iwl) = lai_next(iwl)
       
       call vegndvi( iwl, land%leaf_class(iwl), 1.0, land%veg_height(iwl), &
                     land%veg_ndvip(iwl), land%veg_ndvif(iwl), ta, tb,     &
                     lai_next(iwl), tc, td, te                             )

    enddo
    !$omp end parallel do
       
  end subroutine megan_store_lai


  subroutine megan_avg_temp()

    use leaf_coms,   only: mwl, dt_leaf
    use mem_leaf,    only: land
    use misc_coms,   only: time8
    use consts_coms, only: r8
    use mem_radiate, only: rshort
    use mem_ijtabs,  only: jtab_w, jtw_prog

    implicit none

    real, parameter :: tscali = 2.0 / 86400.0

    real    :: si, so
    integer :: iwl, iw, j

    si = dt_leaf * tscali
    so = 1.0 - si

    !$omp parallel 
    !$omp do
    do iwl = 2, mwl
       avg_temp(iwl) = avg_temp(iwl) * so + land%cantemp(iwl) * si
    enddo
    !$omp end do

    !$omp do private(iw)
    do j = 1, jtab_w(jtw_prog)%jend(1)
       iw = jtab_w(jtw_prog)%iw(j)
       avg_srad(iw) = avg_srad(iw) * so + rshort(iw) * si
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine megan_avg_temp
    

  subroutine megan_driver(mrl)

    use leaf_coms,   only: nvtyp
    use mem_leaf,    only: land
    use mem_basic,   only: vxe, vye, vze, press
    use misc_coms,   only: io6, isubdomain, current_time
    use mem_ijtabs,  only: itabg_w
    use mem_leaf,    only: land, itab_wl, itabg_wl
    use leaf_coms,   only: mwl
    use consts_coms, only: pio180
    use mem_radiate, only: rshort
    use leaf_coms,   only: nzg, slmsts, slcpd

    implicit none

    integer, intent(in) :: mrl

    integer, parameter :: nemis = n_mgn_spc 

    real :: gam_lht  ! LAI correction factor

    real, parameter :: gam_smt = 1.0 ! Soil moisture correction factor

    real :: gam_age(n_mgn_spc)  ! leaf age correction factor
    real :: gam_tld(n_mgn_spc)
    real :: gam_tli(n_mgn_spc)

    real :: em(nemis)
    real :: tmper(n_spca_spc)
    real :: outer(n_mech_spc)

    integer :: j, kw, iw, iwl, s

    real :: glat, glon
    real :: d_ppfd, ppfd
    real :: d_temp, temp
    real :: windspd
    real :: pres
    real :: qv
    real :: LAIc
    real :: di
    real :: arf
    real :: hour, beta, sinbeta
    real :: soilm, soilt, xm, precadj

    real :: cfno, cfnog, cf

    integer :: nmpmg, nmpsp, nmpmc, ip
    integer :: day, jday, jdate
    integer :: gday, glen

    real,    external :: calcbeta
    integer, external :: julday

    real, parameter :: ug2g   = 1.e-6       ! convert microgram to metric gram
    real, parameter :: hr2sec = 1./3600.    ! convert 1/hr to 1/second
    real, parameter :: n2no   = 2.142857    ! nitrogen conversion?

    jday  = julday( current_time%month, current_time%date, current_time%year )
    jdate = current_time%year*1000 + jday

    do iwl = 2, mwl

       iw = itab_wl(iwl)%iw
       if (isubdomain == 1) then
          iw = itabg_w(iw)%iw_myrank
       endif
       kw = itab_wl(iwl)%kw

       glat  = land%glatw(iwl)
       glon  = land%glonw(iwl)
       arf   = land%area(iwl)

!      ppfd: srad - short wave from sun (W/m2)
!      assuming 4.5 (umol m-2 s-1) per (W m-2)
!      assume 1/2 of srad is in 400-700nm band
       d_ppfd = avg_srad(iw) * 4.5 * 0.5
       ppfd   = rshort  (iw) * 4.5 * 0.5

       d_temp = avg_temp(iwl)
       temp   = land%cantemp(iwl)
     
       windspd = sqrt( vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2 )
       pres    = press(kw,iw)
       qv      = land%canshv(iwl)
       precadj = 1.0

       LAIc    = land%veg_lai(iwl) * (1.0 - land%snowfac(iwl))
       
       ! drought index (pass in stomatal conductance instead?)
       di      = 0.0

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
       ! Sinbeta = land%cosz(iwl)

       Beta    = Calcbeta( Day, glat, Hour )
       Sinbeta = SIN(Beta * pio180)
       
       if (sinbeta > 0.0001 .and. ppfd > 1.0 .AND. LAIc > 0.01 .and. ipft(iwl) /= 0) then
          
          call gamma_lai( LAIc, gam_lht )
        ! call gamma_s( gam_smt )

          call gamma_a( lai_prev(iwl), lai_next(iwl), d_temp, gam_age )

          call gamme_ce( sinbeta, temp, d_temp, ppfd, d_ppfd, windspd,   &
                         qv, ipft(iwl), laic, pres, di, gam_tld, gam_tli )

          do s = 1, nemis

             em(s) = GAM_AGE(s) * GAM_SMT * RHO_FCT(S)            &
                   * ( (1.0-LDF_FCT(S)) * GAM_TLI(s) * GAM_LHT +  &
                            LDF_FCT(S)  * GAM_TLD(s) )

          enddo

          ! Determine top-layer soil temperature and moisture

          soilm = land%soil_water(nzg,iwl) / slmsts(land%ntext_soil(nzg,iwl))
          call qwtk( land%soil_energy(nzg,iwl),land%soil_water(nzg,iwl)*1.e3, &
                      slcpd(land%ntext_soil(nzg,iwl)), soilt, xm )

          ! determine soil nox correction

!!          CALL SOILNOX(jdate, TEMP, SOILM, SOILT, LAIc, glat, PRECADJ, CFNO, CFNOG)
!!
!!          call growseason( jdate, glat, gday, glen )

          ! Conversion from MGN 20 to speciated 150

          do s = 1, n_smap_spc

             nmpmg = mg20_map(s)
             nmpsp = spca_map(s)

             if (nmpmg /= INO) then

                tmper(nmpsp) = em(nmpmg) * ef_all(ipft(iwl),nmpmg) * effs_all(ipft(iwl),nmpsp)

             else

                tmper(nmpsp) = 0.0

!!                if (gday == 0) then
!!
!!                   ! non growing season
!!                   ! CFNOG for everywhere
!!                   ! Override crop with grass warm = 14
!!
!!                   if (ipft(iwl) == 15 .or. ipft(iwl) == 16) then
!!                      ip = 14
!!                   else
!!                      ip = ipft(iwl)
!!                   endif
!!
!!                   tmper(nmpsp) = em(nmpmg) * CFNOG * ef_all(ip,nmpmg)  &
!!                                * effs_all(ip,nmpsp) * n2no
!!
!!                else
!!
!!                   ! growing season
!!                   ! CFNOG for everywhere except crops
!!                   ! CFNO for crop and corn
!!                
!!                   if (ipft(iwl) == 15 .or. ipft(iwl) == 16) then
!!                      ! crop
!!                      cf = cfno
!!                   else
!!                      cf = cfnog
!!                   endif
!!
!!                   tmper(nmpsp) = em(nmpmg) * cf * ef_all(ipft(iwl),nmpmg)  &
!!                                * effs_all(ipft(iwl),nmpsp) * n2no
!!
!!                endif

             endif

          enddo

          ! Convert from ug/m^2/hr to mol/m^2/hr using their MW

          do s = 1, n_spca_spc
             tmper(s) = tmper(s) / spca_mwt(s) * ug2g
          enddo

          ! Conversion from speciated species to MECHANISM species,
          ! with units of mol/sec

          outer(:) = 0.0
       
          do s = 1, n_scon_spc
             nmpsp = spmh_map(s)         ! Mapping value for SPCA
             nmpmc = mech_map(s)         ! Mapping value for MECHANISM
             outer(nmpmc) = outer(nmpmc) + tmper(nmpsp) * conv_fac(s) * arf * hr2sec
          enddo
       
          megan_emis(iwl,:) = outer(:)

       else

          megan_emis(iwl,:) = 0.0

       endif

    enddo

  end subroutine megan_driver

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
  SUBROUTINE GAMMA_LAI( lai, gam_l )
    IMPLICIT NONE

    REAL, INTENT(IN)  :: lai
    REAL, INTENT(OUT) :: gam_l

    gam_l = (0.49*lai) / ( (1+0.2*(lai**2))**0.5 )

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
  SUBROUTINE GAMMA_A( LAIp, LAIc, Tt, GAM_A )

    IMPLICIT NONE

    ! input
    real,    intent(in)  :: Tt           ! daily average temperature (K)
    real,    intent(in)  :: LAIp         ! previous LAI
    real,    intent(in)  :: LAIc         ! current (next) LAI

    ! output
    REAL,    INTENT(OUT) :: GAM_A(n_mgn_spc)

    ! Local parameters
    
    INTEGER :: AINDX        ! relative emission activity index
    integer :: s            ! species loop index
    REAL    :: Fnew, Fgro, Fmat, Fold
    REAL    :: ti,tm        ! number of days between budbreak
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
       GAM_A(s) = Fnew*Anew(AINDX) + Fgro*Agro(AINDX) &
                + Fmat*Amat(AINDX) + Fold*Aold(AINDX)

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
!!  SUBROUTINE GAMMA_S( GAM_S )
!!    IMPLICIT NONE
!!
!!    REAL, intent(out) :: GAM_S
!!
!!    GAM_S = 1.0
!!
!!  END SUBROUTINE GAMMA_S

      SUBROUTINE SOILNOX( JDATE, &
                          TA, SOILM, SOILT, &
                          LAIc, LAT, &
                          PRECADJ, &
                          CFNO, CFNOG )

!***********************************************************************
!  DESCRIPTION:
!  
!     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
!     to estimate NO emissions 
!     Information needed to estimate NO emissions
!     Julian Day          (integer)    JDATE
!     Surface Temperature (MCIP field) TA    (K)
!     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (LSOIL)
!          (ratio of volume of water per volume of soil)
!     Soil Temperature    (MCIP field) SOILT (K)         (LSOIL)
!     Soil Type           (MCIP field) ISLTYP            (LSOIL)
!
!     saturation values for soil types (constants)       (LSOIL)
!     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
!                and  the precipitation adjustment factor accounts for saturated soils
!     FOR the non-PX version, the basic algorithm remains with a temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse 
!     PULSEDATE julian date for the beginning of an NO pulse 
!     PULSETIME        time for the beginning of an NO pulse
!  
!     The calculation are based on the following paper by J.J. Yienger and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet and dry soils with
!       the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet and dry adjustment is 
!       calculated at each grid cell.  A linear interpolation between the wet and dry adjustment
!       factor is made using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture is determined by
!       taking the MCIP soil moisture field and dividing by the saturation value defined for each
!       soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total. 
!       THe types of Pulses as described in YL95 were used to estimate the NO emission
!       rate.  
!
!    Also see the following paper for more information:
!    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!    by Tom Pierce and Lucille Bender       
!
!    REFERENCES
!
!    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and Nitric Oxide Emissions from Agricultural Processes
!       Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!        Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start,
!     soil type, SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!     
!  REVISION  HISTORY:
!    10/01 : Prototype by GAP
!    10/03 : modified transition to non growing season for jul-oct of the year
!    08/04 : Converted to SMOKE code style by C. Seppanen
!    07/21/11 : Imported form SMOKE-BEIS v3.14 for MEGAN v2.10
! 
!***********************************************************************

!       USE SOILNOX_FX

        IMPLICIT NONE
        
!.........  INCLUDES
!        INCLUDE 'PARMS3.EXT'      ! I/O API constants
!        INCLUDE 'FDESC3.EXT'      ! I/O API file description data structure
!        INCLUDE 'IODECL3.EXT'     ! I/O API function declarations

!.........  ARGUMENTS and their descriptions
        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)

        REAL, INTENT (IN)  ::  TA          !  air temperature (K)
        REAL, INTENT (IN)  ::  SOILM       !  vol. soil moisture relative to saturation [0-1]
        REAL, INTENT (IN)  ::  SOILT       !  soil temperature (K)
        REAL, INTENT (IN)  ::  PRECADJ     !  precip adjustment
        REAL, INTENT (IN)  ::  LAIc        !  leaf area index
        REAL, INTENT (IN)  ::  LAT         !  Latitude
        REAL, INTENT (OUT) ::  CFNO        !  NO correction factor
        REAL, INTENT (OUT) ::  CFNOG       !  NO correction factor for grass
        
!.........  Local ARRAYS
! Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
!       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
! See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!        INTEGER, PARAMETER :: MAXSTYPES = 11
!        REAL SATURATION( MAXSTYPES )
!        DATA SATURATION / 0.395, 0.410, 0.435, 0.485,
!     &                    0.451, 0.420, 0.477, 0.476,
!     &                    0.426, 0.482, 0.482        /

!.........  SCRATCH LOCAL VARIABLES and their descriptions:
        INTEGER         R, C, L      ! counters
        
        REAL            CF           ! NO correction factor
        REAL            TAIR         ! surface temperature
        REAL            TSOI         ! soil temperature
        REAL            CFNOWET, CFNODRY, RATIO


!        CHARACTER(256)  MESG         ! message buffer
!        CHARACTER(16) :: PROGNAME = 'SOILNOX'   !  program name

!***********************************************************************

 
!.....  Loop through cells
!        DO R = 1, NY
!        DO C = 1, NX

          TAIR = TA         ! unit in degree K

!.......  Check max and min bounds for temperature
          IF (TAIR < 200.0) THEN
             cfno  = 0.0
             cfnog = 0.0
             return
          END IF

!.......  CFNOG
          IF( TAIR > 303.00 ) TAIR = 303.00

          IF ( TAIR > 268.8690 ) THEN  
              CFNOG = EXP( 0.04686 * TAIR - 14.30579 ) ! grass (from BEIS2)
          ELSE
              CFNOG = 0.0
          END IF

!.......  CFNO

          TSOI = SOILT
          IF (TSOI <= 273.16) TSOI = 273.16
          IF (TSOI >= 303.16) TSOI = 303.16

          CFNODRY = (1./3.)*(1./30.)*(TSOI-273.16) ! see YL 1995 Equa 9a p. 11452
          IF (TSOI <= 283.16) THEN ! linear cold case
             CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
          ELSE                  ! exponential case
             CFNOWET = EXP(0.103 * (TSOI-273.16)) * EXP(-0.103 * 30.0)
          END IF

          CF = SOILM * CFNOWET + (1.-SOILM) * CFNODRY

          CFNO = CF * FERTLZ_ADJ( JDATE, LAT ) * VEG_ADJ( LAIc ) * PRECADJ

        END SUBROUTINE SOILNOX



!=======================================================================
!     MODULE SOILNOX_FX
!
!     This module contain functions to assist soil NOx calculation.
!     
!
!     CONTAINS: 1)FERTLZ_ADJ
!               2)VEG_ADJ
!               3)GROWSEASON
!
!     Note:
!
!     Requirement:
!
!
!     Imported from SMOKE-BEIS v3.14 and modified
!          by Tan 07/21/11 for MEGAN v2.10
!
!     Function PRECADJ is moved to MET2MGN
!              PULSETYPE is moved to MET2MGN
!              PRECIPFAC is moved to MET2MGN
!
!     History:
!
!=======================================================================
!!
!!      MODULE SOILNOX_FX
!!
!!      IMPLICIT NONE
!!
!!!...  Program I/O parameters
!!
!!!...  External parameters
!!
!!      CONTAINS
!!



!=======================================================================
!=======================================================================
      REAL FUNCTION FERTLZ_ADJ( DATE, LAT )

!***********************************************************************
!  DESCRIPTION:
!    This internal function computes a fertilizer adjustment factor
!    for the given date in yyyyddd format. If it is not growing 
!    season, the adjustment factor is 0; otherwise, it ranges from
!    0.0 to 1.0.
!
!  CALL:
!    GROWSEASON
!
!  HISTORY:
!    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
!***********************************************************************

      IMPLICIT NONE
            
!.... Function arguments
      INTEGER, INTENT(IN) :: DATE
      REAL,    INTENT(IN) :: LAT

!.... Local variables
      INTEGER  GDAY, GLEN

      CHARACTER(256)  MESG         ! message buffer
!-----------------------------------------------------------------------------

      CALL GROWSEASON( DATE, LAT, GDAY, GLEN )

      IF( GDAY == 0 ) THEN
          FERTLZ_ADJ = 0.
      ELSE IF( GDAY >= 1 .AND. GDAY < 30 ) THEN
          ! first month of growing season
          FERTLZ_ADJ = 1.
      ELSE IF( GDAY >= 30 .AND. GDAY <= 366) THEN
          ! later month of growing season
          FERTLZ_ADJ = 1. + 30. / FLOAT(GLEN) - FLOAT(GDAY) / FLOAT(GLEN)
      ELSE
         ! invalid date?
         FERTLZ_ADJ = 0.
      ENDIF

      END FUNCTION FERTLZ_ADJ
!=======================================================================


!=======================================================================
      REAL FUNCTION VEG_ADJ( LAI )

!***********************************************************************
!  DESCRIPTION
!    This internal function computes a vegetation adjustment factor
!    based on LAIv.  See Yienger and Levy 1995
!    VEG_ADJ = (EXP(-0.24*LAIv)+EXP(-0.0525*LAIv))*0.5 
!
!  CALL
!    NONE
!
!  HISTORY:
!***********************************************************************

      IMPLICIT NONE
      
!...  Function arguments
      REAL,    INTENT(IN) :: LAI

      VEG_ADJ = (EXP(-0.24*LAI)+EXP(-0.0525*LAI))*0.5 

      END FUNCTION VEG_ADJ
!=======================================================================



!=======================================================================
      SUBROUTINE GROWSEASON ( DATE, LAT, GDAY, GLEN )

!~***********************************************************************
!~  DESCRIPTION
!~    This internal function computes the day of the growing season
!~    corresponding to the given date in yyyyddd format.
!~
!~  CALL
!~    JULIAN
!~
!~  HISTORY:
!~    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
!~               Variation of growing season depends on latitude
!~               (Guenther)
!~***********************************************************************

      IMPLICIT NONE

!~.......  Function arguments
      INTEGER, INTENT(IN) :: DATE
      REAL,    INTENT(IN) :: LAT

!~.......  External functions
!      INTEGER, EXTERNAL :: JULIAN

!~.......  Local parameters
      INTEGER            :: GSEASON_START
      INTEGER            :: GSEASON_END

!~.......  Local variables
      INTEGER  YEAR, MONTH, DAY
      INTEGER  JDAY, GDAY, GLEN
      INTEGER  GSJULIAN_START
      INTEGER  GSJULIAN_END
      CHARACTER(256)  MESG         ! message buffer
!     INTEGER G2J

!~-----------------------------------------------------------------------------

!     YEAR = INT( FLOAT( DATE ) / 1000. )
      YEAR = DATE / 1000
      JDAY = DATE - YEAR * 1000

      jday = max(jday,1)
      jday = min(jday,366)

!!      IF( JDAY .LT. 1 .OR. JDAY .GT. 366 ) THEN
!!         WRITE( MESG,94010 ) 'Invalid date specified; date = ',  &
!!                             DATE, 'jday = ', JDAY
!!         CALL M3EXIT( 'GROWSEASON', 0, 0, MESG, 2 )
!!      ENDIF

      IF ( LAT .LE. 23.0 .AND. LAT .GE. -23.0 ) THEN
      ! tropical regions, year round

         GSEASON_START = 0101
         GSEASON_END   = 1231

         GSJULIAN_START = G2J(YEAR, GSEASON_START)
         GSJULIAN_END   = G2J(YEAR, GSEASON_END)
         GDAY = JDAY - GSJULIAN_START + 1
         GLEN = GSJULIAN_END - GSJULIAN_START + 1

      ELSE IF ( LAT .LT. -23.0 ) THEN
      ! southern hemisphere
         IF ( LAT .LT. -60.0 ) THEN
         ! antarctic start = 0 end = 0, no growing
            GDAY = 0
            GLEN = 0
         ELSE
         ! southern hemisphere temperate, NOV - MAY
            GSEASON_START = 1101
            GSEASON_END   = 0531

            GSJULIAN_START = G2J(YEAR, GSEASON_START)
            GSJULIAN_END   = G2J(YEAR, GSEASON_END)

            if (jday >= gsjulian_start) then
               gday = jday - gsjulian_start + 1
            elseif (jday <= gsjulian_end) then
               gday = jday + 61
            else
               gday = 0
            endif

            GLEN = 61 + gsjulian_end
            
         ENDIF

      ELSE   ! IF ( LAT .GT. 23.0 ) THEN

      ! northern hemisphere
         IF ( LAT .GT. 65.0 ) THEN
         ! arctic start = 0 end = 0, no growing season
            GDAY = 0
            GLEN = 0
         ELSE
         ! northern hemisphere temperate
         ! start= (lat-23)*4.5            189
         ! end = 365 -((lat-23)*3.3)      226
            GSEASON_START = 0101
            GSEASON_END   = 1231

            GSJULIAN_START = 1
            GSJULIAN_END   = G2J(YEAR, GSEASON_END)

            GSJULIAN_START = GSJULIAN_START + INT( (LAT-23.0)*4.5 )
            GSJULIAN_END   = GSJULIAN_END   - INT( (LAT-23.0)*3.3 )

            IF (JDAY .GE. GSJULIAN_START .AND. JDAY .LE. GSJULIAN_END) THEN
               GDAY = JDAY - GSJULIAN_START + 1
            ELSE
               GDAY = 0
            ENDIF

            GLEN = GSJULIAN_END - GSJULIAN_START + 1

         ENDIF
!      ELSE
!         MESG = 'Invalid LAT'
!         CALL M3EXIT( 'GROWSEASON', 0, 0, MESG, 2 )
      ENDIF

      END SUBROUTINE GROWSEASON

!=======================================================================

      INTEGER FUNCTION G2J( YYYY, MMDD )
      IMPLICIT NONE

!~.......  Function arguments
      INTEGER, INTENT(IN) :: YYYY
      INTEGER, INTENT(IN) :: MMDD

!~.......  External functions
      INTEGER, EXTERNAL :: julday

!~.......  Local parameters
      INTEGER :: MM
      INTEGER :: DD

!     MM = INT( FLOAT( MMDD ) / 100. )
      MM = MMDD / 100
      DD = MMDD - MM * 100
!     G2J = JULIAN( YYYY, MM , DD )
      g2j = julday( MM, DD, YYYY )

      END FUNCTION G2J

end module mem_megan
