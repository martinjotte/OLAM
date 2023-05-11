module soil_nox

  real, allocatable :: soilnox(:) ! Soil NOx emissions (mol/sec)
  real, allocatable :: pfactor(:) ! Soil NOx pulse factor
  real, allocatable :: drytime(:) ! Time since last precip (hr)

  integer, allocatable      :: lbiome(:)

  integer, parameter        :: nbiom = 23
  character(16), parameter  :: sn_spc = 'NO'
  integer, save             :: sn_map = 0

  ! Steinkamp and Lawrence, 2011 A values, wet biome coefficients
  ! for each of the 24 soil biomes [ng N/m2/s].
  real, parameter :: a_biome(0:nbiom) = (/                         &
       0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.09, 0.09, 0.01, 0.84, &
       0.84, 0.24, 0.42, 0.62, 0.03, 0.36, 0.36, 0.35, 1.66, 0.08, &
       0.44, 0.57, 0.57, 0.57  /)

  ! Canopy wind extinction coefficients
  real,  parameter :: soilexc(0:nbiom) = (/                         &
        0.10, 0.50, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 1.00, &
        1.00, 1.00, 1.00, 2.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, &
        4.00, 2.00, 0.10, 2.00                                      /)

  real, parameter :: sniri(0:nbiom) = (/                            &
       9999., 200.,9999.,9999.,9999.,9999., 200., 200., 200., 200., &
        200., 200., 200., 200., 200., 200., 200., 400., 400., 200., &
        200., 200.,9999., 200.                                     /)

  real, parameter :: snirlu(0:nbiom) = (/                           &
       9999.,9000.,9999.,9999.,9999.,9999.,9000.,9000.,9000.,9000., &
       9000.,9000.,9000.,9000.,9000.,1000.,9000.,9000.,9000.,9000., &
       1000.,9000.,9999.,9000.                                     /)

  real, parameter :: snirac(0:nbiom) = (/                           &
          0., 300.,   0.,   0.,   0.,   0., 100., 100., 100., 100., &
        100., 100., 100., 100.,2000.,2000.,2000.,2000.,2000.,2000., &
       2000., 200., 100., 200.                                     /)

  real, parameter :: snirgss(0:nbiom) = (/                          &
          0.,   0., 100.,1000., 100.,1000., 350., 350., 350., 350., &
        350., 350., 350., 350., 500., 200., 500., 500., 500., 500., &
        200., 150., 400., 150.                                     /)

  real, parameter :: snirgso(0:nbiom) = (/                          &
       2000.,1000.,3500., 400.,3500., 400., 200., 200., 200., 200., &
        200., 200., 200., 200., 200., 200., 200., 200., 200., 200., &
        200., 150., 300., 150.                                     /)

  real, parameter :: snircls(0:nbiom) = (/                          &
       9999.,2500.,9999.,9999.,9999.,9999.,2000.,2000.,2000.,2000., &
       2000.,2000.,2000.,2000.,2000.,9999.,2000.,2000.,2000.,2000., &
       9999.,2000.,9999.,2000.                                     /)

  real, parameter :: snirclo(0:nbiom) = (/                          &
       9999.,1000.,1000.,9999.,1000.,9999.,1000.,1000.,1000.,1000., &
       1000.,1000.,1000.,1000.,1000.,9999.,1000.,1000.,1000.,1000., &
       9999.,1000.,9999.,1000.                                     /)

  real, parameter :: snivsmax(0:nbiom) = (/                         &
         10., 100., 100.,  10., 100.,  10., 100., 100., 100., 100., &
        100., 100., 100., 100., 100., 100., 100., 100., 100., 100., &
        100., 100., 100., 100.                                     /)


contains


  subroutine get_soil_nox()

    use mem_land,    only: mland, omland
    use misc_coms,   only: iparallel
    use mem_sfcg,    only: itab_wsfc
    use mem_para,    only: myrank
    use mem_grid,    only: mwa, mza, lpw
    use mem_ijtabs,  only: jtab_w, jtw_prog
    use mem_radiate, only: cloud_frac
    use olam_mpi_atm,only: mpi_send_w, mpi_recv_w

    implicit none

    integer :: iland, iwsfc, j, iw
    real    :: cfracw(mwa)

    !$omp parallel do private(iw)
    do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
       cfracw(iw) = maxval( cloud_frac(lpw(iw):mza,iw) )
    enddo
    !$omp end parallel do

    if (iparallel == 1) then
       call mpi_send_w(r1dvara1=cfracw)
       call mpi_recv_w(r1dvara1=cfracw)
    endif

    !$omp parallel do
    do iland = 2, mland
       iwsfc = iland + omland

       ! Skip this cell if running in parallel and cell rank is not MYRANK
       if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

       call get_soil_nox_driver(iland, cfracw)

    enddo
    !$omp end parallel do

  end subroutine get_soil_nox



  subroutine get_soil_nox_driver(iland, cfracw)

    use mem_land,    only: land, omland, nzg
    use mem_sfcg,    only: sfcg, itab_wsfc
    use consts_coms, only: t00
    use misc_coms,   only: dtlm
    use therm_lib,   only: qtk
    use mem_grid,    only: mwa

    implicit none

    integer, intent(in) :: iland
    real,    intent(in) :: cfracw(mwa)

    real, parameter :: ng2g    = 1.e-9
    real, parameter :: molwt_N = 14.0
    real, parameter :: ng2molN = ng2g / molwt_N

    logical :: arid
    integer :: j, iw, iwsfc, kb
    real    :: wfps, pcprate, ts_emis
    real    :: lai, tempc, tempk, rhos
    real    :: tsfcwater, liqfrac, cfrac
    real    :: snfac, sndepth, r_canopy, windsqr
    integer :: nlev_sfcwater

    real    :: a_biom, a_fert
    real    :: temp_term, wet_term
    real    :: crf_term

    iwsfc = iland + omland

    snfac = land%snowfac(iland)
    lai   = land%veg_lai(iland) * (1. - snfac)
    tempk = sfcg%cantemp(iwsfc)
    tempc = sfcg%cantemp(iwsfc) - t00
    rhos  = sfcg%rhos   (iwsfc)

    nlev_sfcwater = land%nlev_sfcwater(iland)

    if (nlev_sfcwater == 0) then
       sndepth = 0.0
    else
       sndepth = sum( land%sfcwater_depth(1:nlev_sfcwater,iland) )
    endif

    ! Get surface biome type

    kb = lbiome(iland)

    ! If surface is snow or ice, use snow/ice biome

    if ( nlev_sfcwater > 0 ) then
       call qtk( land%sfcwater_energy(nlev_sfcwater,iland), tsfcwater, liqfrac )
       if (liqfrac < 0.3 .and. sndepth > 0.001) kb = 2
    endif

    ! Get fraction of saturation (water-filled pore space)

    wfps = land%soil_water(nzg,iland) / land%wsat_vg(nzg,iland)

    ! If desert type, mark soil as arid

    if ( sfcg%leaf_class(iwsfc)==3 .or. sfcg%leaf_class(iwsfc)==10 ) then
       arid = .true.
    else
       arid = .false.
    endif

    ! Base emissions per biome. No fertilizer use for now as it is
    ! also partially contained in the edgar42 emissions

    a_biom = a_biome(kb) * ng2molN * sfcg%area(iwsfc)  ! mol N / sec
    a_fert = 0.0

! If the following line gets uncommented, then need to decide whether the land
! cells that contact multiple IW atmosphere cells can use only one of them
! or must use all of them.  In the latter case, it might be best to first average
! edgar/nox values from the multiple atmospheric IW cells.

!   a_fert = edgar42_vars(iw)%nox(9) / arw0(iw) * sfcg%area(iwsfc) * 1000. / 14.

    ! Temperature-dependent term of soil NOx emissions

    temp_term = soiltemp_term( tempc )

    ! Soil moisture scaling of soil NOx emissions

    wet_term  = soilwet_term( wfps, arid )

    ! Cumulative multiplication factor (over baseline emissions)
    ! that accounts for soil pulsing

    pcprate = sfcg%pcpg(iwsfc) / dtlm
    ts_emis = dtlm

    call PULSING( wfps, ts_emis, pcprate, PFACTOR(iland), drytime(iland) )

    ! Compute NOx canopy resistance

    cfrac = 0.0
    do j = 1,itab_wsfc(iwsfc)%nwatm
       iw = itab_wsfc(iwsfc)%iwatm(j)  ! local index
       cfrac = cfrac + itab_wsfc(iwsfc)%arcoarsfc(j) * cfracw(iw)
    enddo

    r_canopy = get_canopy_nox( tempk, rhos, kb, lai,       &
                               sfcg%rshort(iwsfc), land%cosz(iland), cfrac )

    ! Canopy reduction factor

    windsqr  = sfcg%vels(iwsfc)**2
    crf_term = soilcrf( kb, lai, r_canopy, windsqr, land%cosz(iland) )

    ! Soil NOx emissions

    soilnox(iland) = ( A_BIOM + A_FERT )              &
                 * ( TEMP_TERM * WET_TERM * Pfactor(iland) ) &
                 * ( 1.0 - CRF_TERM )

  end subroutine get_soil_nox_driver



  subroutine alloc_soil_nox()

    use mem_land, only: mland

    implicit none

    allocate(soilnox(mland)) ; soilnox = 0.0
    allocate(pfactor(mland)) ; pfactor = 0.0
    allocate(drytime(mland)) ; drytime = 0.0
    allocate(lbiome (mland)) ; lbiome  = 0

  end subroutine alloc_soil_nox


  subroutine filltab_soil_nox()
    use var_tables, only: increment_vtable
    implicit none

!   if (allocated(soilnox)) &
!        call increment_vtable('SOILNOX', 'LW', rvar1=soilnox, hist=.false.)

    if (allocated(pfactor)) &
         call increment_vtable('SNOX_PFACTOR', 'LW', rvar1=pfactor)

    if (allocated(drytime)) &
         call increment_vtable('SNOX_DRYTIME', 'LW', rvar1=drytime)

  end subroutine filltab_soil_nox


  subroutine soil_nox_init()
    use mem_land,    only: land, mland, nzg, omland
    use mem_sfcg,    only: sfcg, itab_wsfc
    use utilio_defn, only: index1, m3exit, xstat1
    use cgrid_spcs,  only: n_gc_emis, gc_emis
    use misc_coms,   only: runtype, iparallel
    use mem_para,    only: myrank

    implicit none

    real          :: wfps
    integer       :: iland, iwsfc
    character(30) :: xmsg

    if (runtype == 'INITIAL') then

       ! Initialize soil nox pulse to 1, and estimate initial drytime
       ! based on top soil layer wetness

       !$omp parallel do private(iland,wfps)
       do iland = 2, mland
          iwsfc = iland + omland

          ! Skip this cell if running in parallel and cell rank is not MYRANK
          if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

          pfactor(iland) = 1.0
          soilnox(iland) = 0.0

          wfps = land%soil_water(nzg,iland) / land%wsat_vg(nzg,iland)

          if (wfps >= 0.5) then
             drytime(iland) = 0.0
          else
             drytime(iland) = 240.0 * (1.0 - wfps / 0.5)
          endif
       enddo
       !$omp end parallel do

    endif

    ! Soil NOx to gas-phase species map

    sn_map = index1( sn_spc, n_gc_emis, gc_emis )
    if ( sn_map == 0 ) then
       xmsg = trim(sn_spc) // ' not found in GC_EMIS table'
       call m3exit( 'SOIL_NOX_INIT', 0, 0, xmsg, xstat1 )
    end if

    ! LEAF to biome map
    do iland = 2, mland
       iwsfc = iland + omland

       ! Skip this cell if running in parallel and cell rank is not MYRANK
       if (iparallel == 1 .and. itab_wsfc(iwsfc)%irank /= myrank) cycle

       lbiome(iland) = olam2sl10( sfcg%leaf_class(iwsfc), sfcg%glatw(iwsfc) )
    enddo

  end subroutine soil_nox_init


  pure real function soiltemp_term ( tc )
    implicit none

    real, intent(in) :: tc   ! Top layer soil temperature [C]

    ! Computes the temperature-dependent term of soil NOx emissions

    ! Based on Ormeci et al., [1999] and Otter et al., [1999]
    ! there exists and entirely exponential relationship between
    ! temperature and soil NOx emissions at constant soil moisture
    ! Therefore we use the following relationship based
    ! on Yienger and Levy et al., [1995] for temperatures 0-30C:
    !
    !      f(T) =  exp( 0.103+/-0.04 * T )
    !        in ng N/m2/s
    !
    !  where T is the temperature in degrees Celsius....Below
    !  0 C, we assume emissions are zero because they are insignificant

    if ( tc <= 0.0 ) then
       soiltemp_term = 0.0
    else
       soiltemp_term = exp( 0.103 * min(tc, 30.0) )
    endif

  end function soiltemp_term



  pure real function soilwet_term ( gwet, arid )
    implicit none

    real,    intent(in) :: gwet  ! water-filled pore space (0 <-> 1)
    logical, intent(in) :: arid  ! is location arid?

    ! Computes the soil moisture scaling of soil NOx emissions

    ! SNOx dependence on soil moisture is best described as a
    ! Poisson function [Parsons et al., 1996; Otter et al., 1999;
    ! Pierce and Aneja, 2000; Kirkman et al., 2001;
    ! van Dijk and Meixner, 2001; van Dijk et al., 2002]:
    !
    !    scaling = a*x*exp(-b*x^2)
    !
    ! where the values of a and b are chosen such that the maximum value
    ! (unity) occurs for WFPS=0.3, which laboratory and field measurements
    ! have found to be the optimal value for emissions in most soils. The
    ! typical range of values are 0.2 (arid) up to 0.45 (floodplain)
    ! [Yang and Meixner, 1997; Ormeci et al., 1999].

    if ( arid ) then
       ! max at wfps = 0.2
       soilwet_term = 8.24 * gwet * exp(-12.5 * gwet * gwet)
    else
       ! max at wfps = 0.3
       soilwet_term = 5.5 * gwet * exp( -5.55 * gwet * gwet)
    endif

  end function soilwet_term



  real function get_canopy_nox( tempk, rho, kk, lai, rad0, cosz, cfrac )
    implicit none

    ! Computes the bulk surface resistance of the canopy to NOx.

    real,    intent(in) :: tempk  ! surface temperature [K]
    real,    intent(in) :: rho    ! surface air density
    integer, intent(in) :: kk     ! biome class
    real,    intent(in) :: lai    ! leaf area index
    real,    intent(in) :: rad0   ! incoming solar [W/m2]
    real,    intent(in) :: cosz
    real,    intent(in) :: cfrac

    REAL, PARAMETER :: XMWH2O = 18.e-3  ! Molecular weight of water [kg]
    REAL, PARAMETER :: HSTAR  = 0.01    ! Henry's law constant
    REAL, PARAMETER :: F0     = 0.1     ! Reactivity factor for biological oxidation
    REAL, PARAMETER :: XMW    = 46.e-3  ! Molecular wt of NO2 (kg)

    REAL            :: DTMP1, DTMP2, DTMP3,  DTMP4, GFACT, GFACI
    REAL            :: RT,    RIX,    RIXX,  RDC,   RLUXX
    REAL            :: RGSX,  RCLX,  TEMPC
    REAL            :: BIO_RESULT

    REAL            :: RI
    REAL            :: RLU
    REAL            :: RAC
    REAL            :: RGSS
    REAL            :: RGSO
    REAL            :: RCLS
    REAL            :: RCLO

    ! Surface temperature [K] and [C]
    TEMPC = tempk - 273.15
    TEMPC = max(TEMPC, -40.0)

    ! Adjust external surface resistances for temperature;
    ! from Wesely [1989], expression given in text on p. 1296.
    RT = 1000.0 * EXP( -TEMPC - 4.0 )

    ! Read the internal resistance RI (minimum stomatal resistance
    ! for water vapor, per unit area of leaf) from the IRI array;
    ! a '9999' value means no deposition to stomata so we impose a
    ! very large value for RI.
    RI = SNIRI(KK)
    IF ( RI > 9998.99 ) RI = 1.e+12

    ! Cuticular resistances IRLU read in from 'drydep.table'
    ! are per unit area of leaf; divide them by the leaf area index
    ! to get a cuticular resistance for the bulk canopy.  If IRLU is
    !'9999' it means there are no cuticular surfaces on which to
    ! deposit so we impose a very large value for RLU.
    IF ( SNIRLU(KK) > 9998.99 .OR. LAI < 1.e-3 ) THEN
       RLU = 1.e+6
    ELSE
       RLU = SNIRLU(KK) / LAI + RT
    ENDIF

    ! The following are the remaining resistances for the Wesely
    ! resistance-in-series model for a surface canopy
    ! (see Atmos. Environ. paper, Fig.1).
    RAC  = MAX( SNIRAC (KK)     , 1.0 )
    RGSS = MAX( SNIRGSS(KK) + RT, 1.0 )
    RGSO = MAX( SNIRGSO(KK) + RT, 1.0 )
    RCLS =      SNIRCLS(KK) + RT
    RCLO =      SNIRCLO(KK) + RT

    IF (  RAC >= 9998.99 ) RAC  = 1.e+12
    IF ( RGSS >= 9998.99 ) RGSS = 1.e+12
    IF ( RGSO >= 9998.99 ) RGSO = 1.e+12
    IF ( RCLS >= 9998.99 ) RCLS = 1.e+12
    IF ( RCLO >= 9998.99 ) RCLO = 1.e+12

    !-------------------------------------------------------------
    ! Adjust stomatal resistances for insolation and temperature:
    !
    ! Temperature adjustment is from Wesely [1989], equation (3).
    !
    ! Light adjustment by the function BIOFIT is described by Wang
    ! [1996].  It combines:
    !
    ! - Local dependence of stomal resistance on the intensity I
    !   of light impinging the leaf; this is expressed as a
    !   multiplicative factor I/(I+b) to the stomatal resistance
    !   where b = 50 W m-2
    !   (equation (7) of Baldocchi et al. [1987])
    ! - Radiative transfer of direct and diffuse radiation in the
    !   canopy using equations (12)-(16) from Guenther et al.
    !   [1995]
    ! - Separate accounting of sunlit and shaded leaves using
    !   equation (12) of Guenther et al. [1995]
    ! - Partitioning of the radiation at the top of the canopy
    !   into direct and diffuse components using a
    !   parameterization to results from an atmospheric radiative
    !   transfer model [Wang, 1996]
    !
    ! The dependent variables of the function BIOFIT are the leaf
    ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and
    ! the fractional cloud cover (CFRAC).  The factor GFACI
    ! integrates the light dependence over the canopy depth; so
    ! be scaled by LAI to yield a bulk canopy value because that's
    ! already done in the GFACI formulation.
    !-------------------------------------------------------------

    ! Internal resistance
    RIX  = RI

    ! Skip the following block if the resistance RIX is high
    IF ( RIX < 9999.0 ) THEN

       IF ( TEMPC > 0. .AND. TEMPC < 39.99999) THEN
          GFACT = 400. / (TEMPC * ( 40.0 - TEMPC ))
       Else
          GFACT = 100.0
       endif

       GFACI = 100.

       IF ( RAD0 > 1.e-3 .AND. LAI > 1.e-3 .and. cosz > 1.e-3 ) THEN

          BIO_RESULT = BIOFIT( LAI, cosz, cfrac )

          GFACI= 1.0 / BIO_RESULT
       ENDIF

       RIX = RIX * GFACT * GFACI
    ENDIF

    ! Compute aerodynamic resistance to lower elements in lower
    ! part of the canopy or structure, assuming level terrain -
    ! equation (5) of Wesely [1989].
    RDC = 100. * (1.0 + 1000.0 / (RAD0 + 10.))

    ! Loop over species; species-dependent corrections to resistances
    ! are from equations (6)-(9) of Wesely [1989].
    !
    ! NOTE: here we only consider NO2 (bmy, 6/22/09)
    RIXX = RIX * DIFFG( TEMPK, rho, XMWH2O ) &
               / DIFFG( TEMPK, rho, XMW    ) &
         + 1.0 / ( HSTAR/3000. + 100.*F0 )

    IF ( RLU < 9999.0 ) THEN
       RLUXX = RLU / ( HSTAR * 1.e-5 + F0 )
    else
       rluxx = 1.e+12
    endif

    ! To prevent virtually zero resistance to species with huge HSTAR,
    ! such as HNO3, a minimum value of RLUXX needs to be set.
    ! The rationality of the existence of such a minimum is
    ! demonstrated by the observed relationship between Vd(NOy-NOx)
    ! and Ustar in Munger et al.[1996]; Vd(HNO3) never exceeds 2 cm/s
    ! in observations. The corresponding minimum resistance is 50 s/m.
    ! was introduced by J.Y. Liang on 7/9/95.
    RGSX = RGSS * RGSO / ( RGSO * HSTAR * 1.e-5 + F0 * RGSS )
    RCLX = RCLS * RCLO / ( RCLO * HSTAR * 1.e-5 + F0 * RCLS )

    ! Get the bulk surface resistance of the canopy
    ! from the network of resistances in parallel and in series
    ! (Fig. 1 of Wesely [1989])
    DTMP1 = 1. / RIXX
    DTMP2 = 1. / RLUXX
    DTMP3 = 1. / ( RAC + RGSX )
    DTMP4 = 1. / ( RDC + RCLX )

    ! Save the within canopy depvel of NOx, used in calculating
    ! the canopy reduction factor for soil emissions [1/s]
    get_canopy_nox = dtmp1 + dtmp2 + dtmp3 + dtmp4

  end function get_canopy_nox



  real function biofit( xlai1, suncos1, cfrac1 )

    real,    intent(in) :: xlai1           ! Leaf area index [cm2/cm2]
    real,    intent(in) :: suncos1         ! Cosine( Solar Zenith Angle )
    real,    intent(in) :: cfrac1          ! Cloud fraction [unitless]

    integer, parameter :: npoly = 20
    real,    parameter :: coeff1(npoly) = (/                                    &
         -3.58E-01, 3.02E+00, 3.85E+00,-9.78E-02,-3.66E+00, 1.20E+01, 2.52E-01, &
         -7.80E+00, 2.26E-01, 2.74E-01, 1.14E+00,-2.19E+00, 2.61E-01,-4.62E+00, &
          6.85E-01,-2.54E-01, 4.37E+00,-2.66E-01,-1.59E-01,-2.06E-01 /)

    ! !REMARKS:
    !  This routine is ancient code from Yuhang Wang.  It was part of the old
    !  Harvard-GISS CTM and was ported into GEOS-Chem.  See this reference for
    !  more information:
    !                                                                             .
    !    Wang, Y., D.J. Jacob, and J.A. Logan, "Global simulation of tropospheric
    !     O3-NOx-hydrocarbon chemistry, 1. Model formulation", J. Geophys. Res.,
    !     103/D9, 10,713-10,726, 1998.

    integer, parameter :: kk = 4
    real               :: term(kk)
    real               :: realterm(npoly)
    integer            :: k, k1, k2, k3

    term(1) = 1.
    term(2) = xlai1
    term(3) = suncos1
    term(4) = cfrac1

    call sunparam_r4(term(2:4))

    k=0
    do k3 = 1, kk
       do k2 = k3, kk
          do k1 = k2, kk
             k = k+1
             realterm(k) = term(k1) * term(k2) * term(k3)
          enddo
       enddo
    enddo

    biofit = 0.0
    do k = 1, npoly
       biofit = biofit + coeff1(k) * realterm(k)
    end do
    if ( biofit .lt. 0.1 ) biofit = 0.1

  end function biofit



  subroutine sunparam_r4( x )
    implicit none

    integer, parameter  :: nn = 3  ! # of variables (LAI, SUNCOS, CLDFRC)

    real, intent(inout) :: x(nn)   ! LAI, SUNCOS, or cloud fraction

    ! !REMARKS:
    !  This routine is ancient code from Yuhang Wang.  It was part of the old
    !  Harvard-GISS CTM and was ported into GEOS-Chem.  See this reference for
    !  more information:
    !                                                                             .
    !    Wang, Y., D.J. Jacob, and J.A. Logan, "Global simulation of tropospheric
    !     O3-NOx-hydrocarbon chemistry, 1. Model formulation", J. Geophys. Res.,
    !     103/D9, 10,713-10,726, 1998.

    real               :: xlow
    integer            :: i

    ! ND = scaling factor for each variable
    integer, parameter :: nd(nn) = (/ 5, 20, 11 /)

    ! X0 = maximum for each variable
    real               :: x0(nn) = (/ 11., 1., 1. /)

    do i = 1, nn
       x(i) = min(x(i),x0(i))

       if (i.ne.3) then
          xlow = x0(i) / real(nd(i))
       else
          xlow = 0.
       endif
       x(i) = max(x(i),xlow)
       x(i) = x(i)/x0(i)
    end do
  end subroutine sunparam_r4



  real function diffg( tk, rho, xm )
    implicit none

    real, intent(in) :: tk      ! Temperature [K]
    real, intent(in) :: rho     ! air density [kg/m3]
    real, intent(in) :: xm      ! Molecular weight of gas [kg]

! Function DiffG calculates the molecular diffusivity [m2/s] in
! air for a gas X of molecular weight XM [kg] at temperature TK [K] and
! density RHO [kg/m3].

    real, parameter  :: xmair  = 28.8e-3
    real, parameter  :: radair = 1.2e-10
    real, parameter  :: pi     = 3.1415926535897932
    real, parameter  :: pii    = 1.0 / pi
    real, parameter  :: radx   = 1.5e-10
    real, parameter  :: rgas   = 8.32
    real, parameter  :: avogad = 6.023e+23
    real, parameter  :: pdiam2 = pi * (radx + radair)**2
    real, parameter  :: const1 = 8.0 * rgas / pi
    real, parameter  :: const2 = 3.0 * pi / 32.0
    real             :: z, frpath, speed

    ! Calculate the mean free path for gas X in air:
    ! eq. 8.5 of Seinfeld [1986];
    z      = xm  / xmair
    frpath = 1.0 / ( pdiam2 * sqrt(1.0 + z) * rho )

    ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
    speed  = sqrt( const1 * tk / xm )

    ! Calculate diffusion coefficient of gas X in air;
    ! eq. 8.9 of Seinfeld [1986]
    diffg = const2 * ( 1.0 + z ) * frpath * speed

  end function diffg



  pure real FUNCTION SoilCrf( K, LAI, CPYNOX, WINDSQR, SUNCOS )
    implicit none

    INTEGER, INTENT(IN) :: K          ! Soil biome type
    REAL,    INTENT(IN) :: LAI        ! Leaf area index [cm2/cm2]
    REAL,    INTENT(IN) :: CPYNOX     ! Bulk sfc resistance to NOx [1/s]
    REAL,    INTENT(IN) :: WINDSQR    ! Square of sfc wind speed [m2/s2]
    REAL,    INTENT(IN) :: SUNCOS     ! Cosine of solar zenith angle

    ! Computes the canopy reduction factor for the soil NOx emissions according
    ! to Jacob \% Bakwin [1991] (and as used in Wang et al [1998]).

    REAL, PARAMETER :: VFDAY   = 1.0e-2
    REAL, PARAMETER :: VFNIGHT = 0.2e-2
    REAL            :: VFNEW

    ! Pick proper ventilation velocity for day or night
    IF ( SUNCOS > 1.e-5 ) THEN
       VFNEW = VFDAY
    ELSE
       VFNEW = VFNIGHT
    ENDIF

    ! If the leaf area index and the bulk surface resistance
    ! of the canopy to NOx deposition are both nonzero ...
    IF ( LAI > 1.e-3 .and. CPYNOX > 1.e-3 ) THEN

       ! Adjust the ventilation velocity.
       ! NOTE: SOILEXC(20) is the canopy wind extinction
       ! coefficient for the tropical rainforest biome.
       VFNEW = VFNEW * SQRT(WINDSQR / 9. * 7. / LAI) * SOILEXC(20) / SOILEXC(K)

       ! Soil canopy reduction factor
       SOILCRF = CPYNOX / ( CPYNOX + VFNEW )

    ELSE

       ! Otherwise set the soil canopy reduction factor to zero
       SOILCRF = 0.0

    ENDIF

  END FUNCTION SoilCrf



  subroutine Pulsing( GWET, TS_EMIS, pcprate, PFACTOR, DRYPERIOD )
    implicit none

    ! Function "pulsing" calculates the increase (or "pulse") of
    ! soil NOx emission that happens after preciptiation falls on dry soil.
    !                                                                             .
    ! According to  Yan et al., [2005] , this pulsing process is thought to
    ! be due to a release of inorganic nitrogen trapped on top of the dry soil
    ! and a subsequent reactivation of water-stressed bacteria, which then
    ! metabolize the excess nitrogen. This can happen in seasonally dry
    ! grasslands and savannahs or over freshly fertilized fields.

    REAL, INTENT(IN)    :: GWET        ! Soil Moisture
    REAL, INTENT(IN)    :: TS_EMIS     ! Emissions timestep [s]
    real, intent(in)    :: pcprate     ! precipitation rate (mm/sec)
    REAL, INTENT(INOUT) :: PFACTOR     ! Pulsing Factor
    REAL, INTENT(INOUT) :: DRYPERIOD   ! Dry period length in hours
    REAL                :: DTSRCE

    ! Soil NOx emissions consist of baseline emissions plus discrete "pulsing"
    ! episodes.  We follow thw Yan et al., [2005] algorithm, where the pulse
    ! (relative to the flux prewetting) is determined by the antecedent dry
    ! period, with a simple logarithmic relationship,
    !
    !     PFACTOR = 13.01 ln ( DRYPERIOD ) -  53.6
    !
    !  where PFACTOR is the magnitude of peak flux relative to prewetting flux,
    !  and DRYPERIOD  is the length of the antecedent dry period in hours.
    !
    !  The pulse decays with
    !
    !     PFACTOR = PFACTOR * EXP( -0.068e+0_hp * DTSRCE )
    !
    !  References:
    !  Yan, X., T. Ohara, and H. Akimoto (2005), Statistical modeling of
    !      global soil NOx emissions, Global Biogeochem. Cycles, 19, GB3019

    ! Emission timestep [s --> hours]
    DTSRCE = TS_EMIS / 3600.0

    ! If soil moisture less than 0.3 and no pulse is taking place
    IF ( GWET < 0.3 .and. PFACTOR < 1.0000001) THEN

       ! If there is precipitation
       IF ( pcprate > 2.e-5 ) THEN

          ! Initialize new pulse factor (dry period hours)
          IF ( DRYPERIOD > 1. ) THEN
             PFACTOR = 13.01 * LOG( DRYPERIOD ) - 53.6
          ELSE
             PFACTOR = -53.6
          ENDIF

          ! If dry period < ~3 days then no pulse
          IF ( PFACTOR < 1.0 ) PFACTOR = 1.0

          ! Reinitialize dry period
          DRYPERIOD = 0.0

       ! If no rain
       ELSE

          ! Add one timestep to dry period
          DRYPERIOD = DRYPERIOD + DTSRCE

       ENDIF

    ! If box is already pulsing , then decay pulse one timestep
    ELSEIF ( PFACTOR >= 1.0000001) THEN

       ! Decay pulse
       PFACTOR   = PFACTOR * EXP( -0.068 * DTSRCE )

       ! Update dry period
       IF ( GWET < 0.3 ) DRYPERIOD = DRYPERIOD + DTSRCE

       ! If end of pulse
       IF ( PFACTOR < 1.0 ) PFACTOR = 1.0

    ENDIF

  END subroutine Pulsing




  pure integer function olam2sl10( lclass, glat )
    implicit none

    integer, intent(in) :: lclass
    real,    intent(in) :: glat

    ! Converts the OLAM land surface types to one of the 24
    ! landcover types adopted by Steinkamp & Lawrence for
    ! simulating soil NOx emissions:
    !
    !  0 Water
    !  1 Permanent wetland
    !  2 Snow and ice
    !  3 Barren (polar)
    !  4 Unclassified
    !  5 Barren (warm desert)
    !  6 Closed shrubland
    !  7 Open shrubland (warm)
    !  8 Open shrubland (polar)
    !  9 Grassland (polar)
    ! 10 Savannah (polar)
    ! 11 Savannah (warm)
    ! 12 Grassland (warm)
    ! 13 Woody savannah
    ! 14 Mixed forest
    ! 15 Evergr. broadl. forest (temperate/polar)
    ! 16 Dec. broadl. forest (temperate/polar)
    ! 17 Dec. needlel. forest
    ! 18 Evergr. needlel. forest
    ! 19 Dec. broadl. forest (tropical)
    ! 20 Evergr. broadl. forest (tropical)
    ! 21 Cropland
    ! 22 Urban and build-up lands
    ! 23 Cropland/nat. veg. mosaic

    select case( lclass )

    case(0)  ! ocean

       olam2sl10 = 0

    case(1)  ! lakes

       olam2sl10 = 0

    case(2)  ! glacier

       olam2sl10 = 2

    case(3)  ! desert

       if (abs(glat) > 60.0) then
          olam2sl10 = 3
       else
          olam2sl10 = 5
       endif

    case(4)  ! Evergreen needleleaf tree

       olam2sl10 = 18

    case(5)  ! Deciduous needleleaf tree

       olam2sl10 = 17

    case(6)  ! Deciduous broadleaf tree

       if (abs(glat) > 23.0) then
          olam2sl10 = 16
       else
          olam2sl10 = 19
       endif

    case(7,20)  ! Evergreen broadleaf tree

       if (abs(glat) > 23.0) then
          olam2sl10 = 15
       else
          olam2sl10 = 20
       endif

    case(8,9)  ! Short and tall grass

       if (abs(glat) > 60.0) then
          olam2sl10 = 9
       else
          olam2sl10 = 12
       endif

    case (10)  ! semi-desert

       if (abs(glat) > 60.0) then
          olam2sl10 = 8
       else
          olam2sl10 = 7
       endif

    case(11)  ! tundra

       olam2sl10 = 8

    case(12,13)  ! Evergreen or Deciduous shrub

       olam2sl10 = 6

    case(14)  ! Mixed woodland

       olam2sl10 = 14

    case(15,16)  ! Crop

       olam2sl10 = 21

    case(17)  ! Bog or marsh

       olam2sl10 = 1

    case(18)  ! Wooded grassland

       olam2sl10 = 13

    case(19,21)  ! Urban or Very Urban

       olam2sl10 = 22

    case default

       olam2sl10 = 0

    end select

  end function olam2sl10

end module soil_nox
