
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!

module ltng_defn

  ! Function: production of NO from lightning

  implicit none

  logical,           save :: ltng_no           ! flag for lightning NO emissions
  integer,           save :: ltng_map          ! map to gas chemistry NO
  real, allocatable, save :: vdemis_lt(:,:)    ! lightning emis
  real, allocatable, save :: column_ltng_no(:) ! column total NO
  character(2), parameter :: ltspc = 'NO'      ! lightning emis species name

  public ltng_no, ltng_map, vdemis_lt, ltng_init, get_ltng, column_ltng_no
  private


contains


  function ltng_init ( ) result ( success )

    use utilio_defn
    use cgrid_spcs, only: n_gc_emis, gc_emis
    use mem_grid,   only: mwa, mza
    use oname_coms, only: nl
    use misc_coms,  only: io6

    implicit none

    logical          :: success
    character( 16 )  :: pname = 'LTNG_INIT'
    character( 120 ) :: xmsg  = ''
    integer          :: status

    success = .true.

    if (nl%ltng_nox == 1) then
       ltng_no = .true.
       write(io6,*) "Enabling lightning NOx production."
    else
       ltng_no = .false.
       write(io6,*) "Disabling lightning NOx production."
    endif

    if ( .not. ltng_no ) return

    allocate( column_ltng_no( mwa ) )
    column_ltng_no = 0.0

    ! Lightning to gas-phase species map

    ltng_map = index1( ltspc, n_gc_emis, gc_emis )
    if ( ltng_map .eq. 0 ) then
       xmsg = trim( ltspc ) // ' not found in gc_emis table'
       call m3exit( pname, 0, 0, xmsg, xstat1 )
    end if

    allocate( vdemis_lt( mza, mwa ), stat = status )
    if ( status .ne. 0 ) then
       xmsg = 'VDEMIS_LT memory allocation failed'
       call m3warn ( pname, 0, 0, xmsg )
       success = .false.
    end if
    vdemis_lt = 0.0

  end function ltng_init


  subroutine get_ltng ( mrl )

    ! Get NO produced from lightning in VDEMIS_LT

    use mem_basic,  only: tair, press
    use mem_grid,   only: mza, lpw, arw0, zm
    use mem_cuparm, only: kcutop, kcubot, conprr
    use mem_ijtabs, only: jtab_w, jtw_prog
    use mem_turb,   only: frac_sea
    use consts_coms,only: t00
    
    implicit none

    integer, intent(in) :: mrl

    real, parameter :: tmin = t00 - 40.0
    real, parameter :: tmax = t00 - 15.0

    ! 147 CG flashes per cm precip in 36 km^2 times moles N per flash
    ! and 1 cm / 10 mm gives moles N per m^2 per mm precip:
    real, parameter :: molesN_per_flash = 350.
    real, parameter :: cg_to_total      = 1. / 4.
    real, parameter :: scale = 147. / 10. / 36000.**2 * molesN_per_flash / cg_to_total
    integer         :: j, iw, k
    real            :: cc_thick, rate, fsea, fland

    if ( .not. ltng_no ) return

    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       vdemis_lt(1:mza,iw) = 0.0
       column_ltng_no (iw) = 0.0
            
       ! First make sure that there was deep convection

       if ( .not. ( kcubot(iw) >= lpw(iw) .and. kcutop(iw) > kcubot(iw) &
                    .and. conprr(iw) > 1.e-9 ) ) cycle

       ! And that the cloud top was sufficiently below freezing to generate
       ! positive and negative charge areas

       if (tair(kcutop(iw),iw) > tmin) cycle

       ! And that the cloud base temperature is warm enough to have
       ! water droplets in the cloud

       if (tair(kcubot(iw),iw) < tmax) cycle

       ! Compute lightning flash rate. Scale is the moles of N produced
       ! per m^2 per mm precip, resulting in moles/sec

       rate = arw0(iw) * conprr(iw) * scale

       ! Reduce lightning over oceans
       
       column_ltng_no(iw) = rate * ( 1.0 - 0.9 * frac_sea(iw) )

       ! Get the vertical distibution of NOx from lightning

       call lightdist( iw, kcutop(iw), column_ltng_no(iw), vdemis_lt(:,iw) )
               
    end do

  end subroutine get_ltng


  pure function get_cg_ratio( ccthick ) result( f_cg )
    implicit none

    real, intent(in) :: ccthick ! Cold cloud thickness [m]
    real             :: cc      ! Cold cloud thickness [km]
    real             :: f_cg    ! ratio of CG to total flashes

    ! COMPUTE RATIO OF CLOUD-TO-GROUND FLASHES TO TOTAL FLASHES
    !
    ! Price & Rind (1993) compute the ratio of Cloud-Ground 
    ! to Total Flashes by the parameterization:
    !
    ! For 5.5 < dz < 14:
    !
    !     f_CG = 1 / (( A*dz^4 + B*dz^3 + C*dz^2 + D*dz + E ) + 1 )
    !
    ! Where:
    !
    ! (1) dz is the depth [km] of the cloud above the freezing 
    !     level.  The cold-cloud thickness (dz) is the depth of 
    !     the layer between the cloud top and the center of the 
    !     highest layer for which the temperature exceeds 273 K. 
    !     The cold-cloud thickness is set to 5.5 km at grid points 
    !     where it is less than 5.5 km.

    ! Convert cold cloud thickness from [m] to [km] (min value: 5.5 km)

    cc = max( ccthick * 1.e-3, 5.5 )
    cc = min( cc, 12.0 )
      
    f_cg = 1.0 / ( 1.0 + 63.09 + cc * ( -36.540 + &
                                 cc * (   7.493 + &
                                 cc * (  -0.648 + &
                                 cc * (   0.021 ) ) ) ) )

  end function get_cg_ratio


  subroutine lightdist( iw, ktop, total, vertprof )

    use misc_coms, only: current_time
    use mem_grid,  only: mza, glatw, zt, zm, lpw
    use mem_turb,  only: frac_sea

    implicit none

    integer, intent(in)  :: iw
    integer, intent(in)  :: ktop
    real,    intent(in)  :: total
    real,    intent(out) :: vertprof(mza)

    real    :: frac(mza)
    real    :: glat, r1, r0, eta, scale
    integer :: mtype, k, in, ka

    ! lightning NOx vertical cummulative source distributions 
    ! from Ott et al [JGR, 2010]; 50 evenly spaced levels from
    ! surface to cloud top for 1) tropical marine, 2) tropical
    ! land, 3) midlatitudes, and 4) subtropics

    real, parameter :: profiles(0:50,4) = reshape( (/         &
         0.0000,  0.0020,  0.0040,  0.0062,  0.0113,  0.0164, &
         0.0220,  0.0319,  0.0418,  0.0524,  0.0670,  0.0816, &
         0.0971,  0.1155,  0.1338,  0.1535,  0.1762,  0.1990, &
         0.2230,  0.2492,  0.2754,  0.3026,  0.3315,  0.3605, &
         0.3911,  0.4237,  0.4564,  0.4901,  0.5248,  0.5595, &
         0.5948,  0.6305,  0.6662,  0.7013,  0.7359,  0.7706, &
         0.8006,  0.8284,  0.8563,  0.8798,  0.9019,  0.9240, &
         0.9406,  0.9559,  0.9712,  0.9796,  0.9871,  0.9946, &
         0.9966,  0.9983,  1.0000,                            &
         0.0000,  0.0007,  0.0014,  0.0021,  0.0038,  0.0055, &
         0.0072,  0.0093,  0.0113,  0.0138,  0.0185,  0.0233, &
         0.0291,  0.0382,  0.0474,  0.0579,  0.0715,  0.0851, &
         0.0998,  0.1168,  0.1338,  0.1525,  0.1736,  0.1947, &
         0.2195,  0.2488,  0.2780,  0.3103,  0.3453,  0.3803, &
         0.4180,  0.4574,  0.4969,  0.5381,  0.5802,  0.6224, &
         0.6653,  0.7085,  0.7517,  0.7941,  0.8363,  0.8784, &
         0.9072,  0.9331,  0.9589,  0.9710,  0.9812,  0.9914, &
         0.9946,  0.9973,  1.0000,                            &
         0.0000,  0.0066,  0.0133,  0.0201,  0.0298,  0.0438, &
         0.0578,  0.0718,  0.0915,  0.1122,  0.1329,  0.1551, &
         0.1812,  0.2072,  0.2333,  0.2619,  0.2916,  0.3213, &
         0.3512,  0.3832,  0.4151,  0.4470,  0.4791,  0.5113, &
         0.5435,  0.5757,  0.6066,  0.6374,  0.6682,  0.6977, &
         0.7254,  0.7532,  0.7809,  0.8048,  0.8280,  0.8513, &
         0.8730,  0.8906,  0.9083,  0.9259,  0.9394,  0.9511, &
         0.9629,  0.9739,  0.9800,  0.9862,  0.9924,  0.9958, &
         0.9972,  0.9986,  1.0000,                            &
         0.0000,  0.0032,  0.0064,  0.0096,  0.0158,  0.0225, &
         0.0292,  0.0402,  0.0527,  0.0652,  0.0814,  0.0999, &
         0.1185,  0.1400,  0.1647,  0.1893,  0.2158,  0.2456, &
         0.2754,  0.3061,  0.3397,  0.3733,  0.4071,  0.4423, &
         0.4775,  0.5127,  0.5479,  0.5831,  0.6184,  0.6519, &
         0.6852,  0.7185,  0.7489,  0.7784,  0.8078,  0.8339, &
         0.8579,  0.8819,  0.9027,  0.9203,  0.9379,  0.9530, &
         0.9639,  0.9748,  0.9842,  0.9890,  0.9938,  0.9981, &
         0.9987,  0.9994,  1.0000 /), (/51,4/)                )
    
    mtype           = 0
    vertprof(1:mza) = 0.0
    glat            = glatw(iw)
    ka              = lpw(iw)
      
    ! Assign profile kind to grid box, following Allen et al. [JGR, 2010]

    select case (current_time%month)

    ! Southern Hemisphere Summer

    case ( 1,2,3,12 )

       if ( abs(glat) .le. 15 ) then
          if (frac_sea(iw) > 0.75) then
             mtype = 1        ! Tropical marine
          else
             mtype = 2        ! Tropical continental
          endif
       else if ( ( glat .gt. 15. ) .and. ( glat .le. 30. ) ) then
          mtype = 4           ! N. Subtropics
       else if ( ( glat .ge. -40. ) .and. ( glat .lt. -15. ) ) then
          mtype = 4           ! S. Subtropics
       else
          mtype = 3           ! Midlatitude
       endif

    ! Equinox months

    case ( 4,5,10,11 )

       if ( abs(glat) .le. 15 ) then
          if (frac_sea(iw) > 0.75) then
             mtype = 1        ! Tropical marine
          else
             mtype = 2        ! Tropical continental
          endif
       else if ( abs(glat) .le. 30 ) then
          mtype = 4           ! Subtropics
       else
          mtype = 3           ! Midlatitude
       endif

       ! Northern Hemisphere Summer

    case ( 6,7,8,9 )

       if ( abs(glat) .le. 15 ) then
          if (frac_sea(iw) > 0.75) then
             mtype = 1        ! Tropical marine
          else
             mtype = 2        ! Tropical continental
          endif
       else if ( ( glat .gt. 15. ) .and. ( glat .le. 40. ) ) then
          mtype = 4           ! N. Subtropics
       else if ( ( glat .ge. -30. ) .and. ( glat .lt. -15. ) ) then
          mtype = 4           ! S. Subtropics
       else
          mtype = 3           ! Midlatitude
       endif
         
    end select

    ! Safety check

    if ( mtype == 0 ) return
      
    ! Use the profile for this type of lightnox event to partition
    ! the total column lightnox onto the vertical grid

    scale = 50.0 / (zm(ktop) - zm(ka-1))

    do k = ka, ktop-1
       eta = (zm(k) - zm(lpw(iw)-1)) * scale  ! 0 <-> 50

       in = max(0, min(49, int(eta)))
       r1 = eta - real(in)
       r0 = 1.0 - r1
  
       frac(k) = r0 * profiles(in,mtype) + r1 * profiles(in+1,mtype)
    enddo

    frac(ka-1) = 0.0
    frac(ktop) = 1.0

    do k = ka, ktop
       vertprof(k) = (frac(k) - frac(k-1)) * total
    enddo

  end subroutine lightdist


  subroutine flashes_cth( fland, fsea, height )
    implicit none

    real, intent(in)  :: height
    real, intent(out) :: fland
    real, intent(out) :: fsea

    ! COMPUTE LIGHTNOX FLASH RATE / SECOND
  
    ! Price & Rind (1992) give the following parameterizations for
    ! lightnox flash rates as a function of convective cloud top
    ! height. LightNOX will therefore occur much more often on land

    ! Flashes/sec over land
    fland = 5.733e-7 * ( height * 1.e-3 )**4.90
      
    ! Flahes/sec over sea
    fsea  = 1.067e-5 * ( height * 1.e-3 )**1.73

  end subroutine flashes_cth

end module ltng_defn
