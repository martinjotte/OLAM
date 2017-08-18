module mem_megan

  use CONSTS_MEGAN
  use EMIS_MAPS_MEGAN
  use consts_coms, only: r8
  implicit none

  private :: r8

  integer, allocatable         :: ipft    (:)
  real,    allocatable, target :: avg_temp(:)
  real,    allocatable, target :: avg_srad(:)
  real,    allocatable, target :: avg_srad_dif(:)
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

  integer,   allocatable, save :: mgn_2_gc_map(:)


contains


  subroutine alloc_megan(mwl, mwa)
    use utilio_defn, only: index1
    use cgrid_spcs,  only: n_gc_emis

    implicit none

    integer, intent(in) :: mwl
    integer, intent(in) :: mwa

    allocate(avg_temp(mwl)) ; avg_temp = 0.0
    allocate(ipft    (mwl)) ; ipft     = 0
    allocate(lai_prev(mwl)) ; lai_prev = 0.0
    allocate(lai_next(mwl)) ; lai_next = 0.0

    allocate(avg_srad    (mwa)) ; avg_srad     = 0.0
    allocate(avg_srad_dif(mwa)) ; avg_srad_dif = 0.0

    ! Store index of NO
    INO = INDEX1('NO', N_MGN_SPC, MGN_SPC)

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

    allocate(megan_emis(mwl,n_mech_spc))

  end subroutine alloc_megan


  subroutine filltab_megan()
    use var_tables, only: increment_vtable
    implicit none

    if (allocated(avg_temp)) then
       call increment_vtable('MEGAN_AVGTEMP', 'LW', rvar1=avg_temp)
    endif

    if (allocated(avg_srad)) then
       call increment_vtable('MEGAN_AVGSRAD', 'AW', rvar1=avg_srad)
    endif

    if (allocated(avg_srad_dif)) then
       call increment_vtable('MEGAN_AVGSRAD_DIF', 'AW', rvar1=avg_srad_dif)
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
    use misc_coms,    only: current_time, runtype, io6
    use mem_grid,     only: mwa, glatw
    use cgrid_spcs,   only: n_gc_emis, gc_emis
    use utilio_defn,  only: index1
    implicit none

    real    :: ta, tb, tc, td, te
    integer :: iwl, iw, n

    integer :: day
    real    :: hour, beta, sinbeta_max

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

          beta             = calcbeta( Day, glatw(iw), 12.0 )
          sinbeta_max      = max( sin(beta * pio180), 0.0 )
          avg_srad    (iw) = solcon * sinbeta_max * 0.25
          avg_srad_dif(iw) = solcon * sinbeta_max * 0.05
       enddo
       !$omp end do

    endif
    !$omp end parallel

    do n = 1, n_gc_emis
       mgn_2_gc_map(n) = INDEX1( GC_EMIS(n), n_mech_spc, mech_spc )

       !if ( mgn_2_gc_map(n) > 0) then
       !   write(io6,*) "MEGAN: ", trim(gc_emis(n)), " is mapped to ", trim(mech_spc(mgn_2_gc_map(n)))
       !else
       !   write(io6,*) "MEGAN: ", trim(gc_emis(n)), " is not mapped"
       !endif
    enddo

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
    use consts_coms, only: r8
    use mem_radiate, only: ppfd, ppfd_diffuse
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
       avg_srad    (iw) = avg_srad    (iw) * so + ppfd        (iw) * si
       avg_srad_dif(iw) = avg_srad_dif(iw) * so + ppfd_diffuse(iw) * si
    enddo
    !$omp end do
    !$omp end parallel

  end subroutine megan_avg_temp


  subroutine megan_driver(mrl)

    use leaf_coms,   only: mwl
    implicit none

    integer, intent(in) :: mrl
    integer             :: iwl

    !$omp parallel do
    do iwl = 2, mwl
       call megan_driver1(iwl)
    enddo
    !$omp end parallel do

  end subroutine megan_driver


  subroutine megan_driver1(iwl)

    use leaf_coms,   only: nvtyp
    use mem_leaf,    only: land
    use misc_coms,   only: isubdomain, current_time
    use mem_ijtabs,  only: itabg_w
    use mem_leaf,    only: land, itab_wl
    use consts_coms, only: pio180
    use mem_radiate, only: ppfd, ppfd_diffuse
    use leaf_coms,   only: nzg

    implicit none

    integer, intent(in) :: iwl

    real :: gam_lht             ! LAI activity factor
    real :: gam_age(n_mgn_spc)  ! leaf age activity factor
    real :: gam_tld(n_mgn_spc)  ! light-dependent activity factor
    real :: gam_tli(n_mgn_spc)  ! light-independent activity factor
    real :: gam_smt(n_mgn_spc)  ! soil moisture activity factor

    real :: em(n_mgn_spc)
    real :: tmper(n_spca_spc)
    real :: outer(n_mech_spc)

    integer :: kw, iw, s

    real :: glat, glon, hour
    real :: LAIc
    real :: arf
    real :: beta, sinbeta

    integer :: nmpmg, nmpsp, nmpmc
    integer :: day, jday, jdate

    integer, external :: julday

    real, parameter :: ug2g   = 1.e-6       ! convert microgram to metric gram
    real, parameter :: hr2sec = 1./3600.    ! convert 1/hr to 1/second
    real, parameter :: n2no   = 2.142857    ! nitrogen conversion?

    real :: ppfd0, ppfd0_dif, ppfd24, ppfd24_dif, temp0, temp24, vegtk

    jday  = julday( current_time%month, current_time%date, current_time%year )
    jdate = current_time%year*1000 + jday

    iw = itab_wl(iwl)%iw
    if (isubdomain == 1) then
       iw = itabg_w(iw)%iw_myrank
    endif
    kw = itab_wl(iwl)%kw

    LAIc = land%veg_lai(iwl) * (1.0 - land%snowfac(iwl))
    glat = land%glatw(iwl)
    glon = land%glonw(iwl)
    arf  = land%area(iwl)

    if (laic > 0.01 .and. land%vf(iwl) > 0.01 .and. ipft(iwl) /= 0) then

       ppfd0      = ppfd(iw)
       ppfd0_dif  = ppfd_diffuse(iw)
       ppfd24     = avg_srad(iw)
       ppfd24_dif = avg_srad_dif(iw)

       temp24 = avg_temp(iwl)
       temp0  = land%cantemp(iwl)
       vegtk  = land%veg_temp(iwl)

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

       call gamma_lai( gam_lht, laic )

       call gamma_sm( gam_smt, land%leaf_class(iwl), &
            land%soil_water(1:nzg,iwl), land%ntext_soil(1:nzg,iwl) )

       call gamma_a( gam_age, lai_prev(iwl), lai_next(iwl), temp24 )

       call gamma_ce( gam_tld, gam_tli, sinbeta, temp0, temp24, vegtk,      &
            ppfd0, ppfd0_dif, ppfd24, ppfd24_dif, ipft(iwl), laic )

       do s = 1, n_mgn_spc

          em(s) = GAM_AGE(s) * GAM_SMT(s) * &
               ( (1.0-LDF_FCT(S)) * GAM_TLI(s) * GAM_LHT + LDF_FCT(S) * GAM_TLD(s) )
       enddo

       ! We do soil NOx elsewhere, so turn off NOx here

       em(INO) = 0.0

       ! Conversion from MGN 20 to speciated 150

       do s = 1, n_smap_spc

          nmpmg = mg20_map(s)
          nmpsp = spca_map(s)

          tmper(nmpsp) = em(nmpmg) * ef_all(ipft(iwl),nmpmg) * effs_all(ipft(iwl),nmpsp)

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
          outer(nmpmc) = outer(nmpmc) + tmper(nmpsp) * conv_fac(s)
       enddo

       do s = 1, n_mech_spc
          megan_emis(iwl,s) = outer(s) * arf * hr2sec
       enddo

    else

       do s = 1, n_mech_spc
          megan_emis(iwl,s) = 0.0
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
  SUBROUTINE GAMMA_SM( GAM_S, leaf_class, soil_water, ntext_soil )

    use leaf_coms, only: nzg, kroot, soilwilt
    IMPLICIT NONE

    REAL,    intent(out) :: GAM_S(n_mgn_spc)
    integer, intent(in)  :: leaf_class
    real,    intent(in)  :: soil_water(nzg)
    integer, intent(in)  :: ntext_soil(nzg)

    real, parameter :: dsw = 0.06
    real            :: smk
    integer         :: k, nts, s

    gam_s(1) = 0.0

    do k = kroot(leaf_class), nzg
       nts = ntext_soil(k)

       if (soil_water(k) <= soilwilt(nts)) then
          smk = 0.0
       elseif (soil_water(k) >= soilwilt(nts) + dsw) then
          smk = 1.0
       else
          smk = (soil_water(k) - soilwilt(nts)) / dsw
       endif

       gam_s(1) = max(gam_s(1), smk)
    enddo

    do s = 2, n_mgn_spc
       gam_s(s)  = 1.0
    enddo

  END SUBROUTINE GAMMA_SM


  subroutine gamma_ce( gamma_tld, gamma_tli, sinbeta, temp, temp24, veg_temp, &
                       ppfd, ppfd_dif, ppfd24, ppfd24_dif, cantype, lai )
     
    implicit none

    ! input
    integer, intent(in) :: cantype
    real,    intent(in) :: sinbeta, lai
    real,    intent(in) :: temp, temp24, veg_temp
    real,    intent(in) :: ppfd, ppfd_dif, ppfd24, ppfd24_dif

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
          gamma_tld(s) = ( gamma_p_sun + gamma_p_shade ) * gamma_t * cce * lai
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

       X = ((1.0 / Topt) - (1.0 / T1)) / 0.00831

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
