
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

module sedv_defn

  implicit none

  ! no. of surrogates for aero settling velocities
  integer, parameter :: n_ae_sed_spc = 6 

  ! set up species indices for settling velocity internal array vsed
  integer, parameter :: vgnacc = 1, & ! accumulation mode number
                        vgncor = 2, & ! coarse mode number
                        vgsacc = 3, & ! accumulation mode surface area
                        vgscor = 4, & ! coarse mode surface area
                        vgmacc = 5, & ! accumulation mode mass
                        vgmcor = 6     ! coarse mode mass

  ! follow the Namelist dep vel surrogate name table
  character( 16 ), parameter :: vgae_name( n_ae_sed_spc ) = &
                                (/ 'VNUMACC', &
                                   'VNUMCOR', &
                                   'VSRFACC', &
                                   'VSRFCOR', &
                                   'VMASSJ ', &
                                   'VMASSC '  /)

  integer, allocatable, save :: sedi_sur( : )   ! pointer to surrogate


contains

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine aero_sedi( mrl )

    use mem_grid,    only: mza, lpw, volt, volti, arw
    use mem_ijtabs,  only: jtab_w, jtw_prog
    use misc_coms,   only: io6, dtlm
    use consts_coms, only: r8
    use mem_basic,   only: rho
    use utilio_defn
    use cgrid_conv
    use cgrid_spcs
    use cgrid_defn

    implicit none

    integer, intent(in) :: mrl

    real, parameter :: alpha = 1.1

    integer :: j, iw, s, v, n, off, k
    integer :: astat

    logical, save :: firstime = .true.
    integer, save :: logdev

    character( 120 )           :: xmsg = ' '
    character( 16 ), parameter :: pname = 'aero_sedi'

    real(r8) :: eps, dt8, time
    real     :: dt
    real(r8) :: dtmin( n_ae_sed_spc )
    real     :: vsed_ae( n_ae_sed_spc, mza )
    real     :: conc   ( mza, nspcsd )
    real     :: fplus(mza,n_ae_sed_spc), fminus(mza,n_ae_sed_spc)
    real     :: aplus(mza), aminus(mza)

    if ( firstime ) then

       firstime = .false.
       logdev = init3()

       allocate( sedi_sur( n_ae_spc ), stat = astat )
       if ( astat .ne. 0 ) then
          xmsg = 'Failure allocating sedi_sur'
          call m3exit( pname, 0, 0, xmsg, xstat1 )
       end if

       ! Set the settling vel surrogate pointers according to the depv table

       if ( n_ae_depv .ne. n_ae_spc ) then
          xmsg = 'Aerosol sedimention is assuming n_ae_spc = n_ae_depv'
          call m3exit( pname, 0, 0, xmsg, xstat1 )
       end if

       j = 0
       do v = 1, n_ae_depv   ! assume n_ae_spc = n_ae_depv
          n = index1( ae_depv( v ), n_ae_sed_spc, vgae_name )
          if ( n .ne. 0 ) then
             j = j + 1
             sedi_sur( v ) = n
          else
             ! write( logdev,* ) ' surrogate ', trim( ae_depv( v ) ), &
             !                   ' not used for', v, trim( ae_spc( v ) )
             sedi_sur( v ) = 0
          end if
       end do

       !n = j
       !write( logdev,* ) n, ' Aerosol species with a grav. settling vel'
       !do j = 1, n_ae_spc
       !   n = sedi_sur( j )
       !   if ( n .ne. 0 ) write( logdev,'( i3, 2x, a9, i3, 2x, a )' ) &
       !                   j, ae_spc( j ), n, trim( ae_depv( j ) )
       !end do

    endif  ! firstime

    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       call rev_cgrid_iw( iw, conc )

       call aero_sedv( iw, conc, vsed_ae )

       vsed_ae(:,mza)       = 0.0
       vsed_ae(:,lpw(iw)-1) = 0.0

       do n = 1, n_ae_sed_spc
          dtmin(n) = alpha * minval( volt(lpw(iw):mza-1,iw) &
                                   / arw(lpw(iw):mza-1,iw) / vsed_ae(n,lpw(iw):mza-1) )
       enddo

       do k = lpw(iw), mza-1
          aplus (k) = arw(k  ,iw) * volti(k,iw) * rho(k+1,iw) / rho(k,iw)
          aminus(k) = arw(k-1,iw) * volti(k,iw)
       enddo
       aplus (mza) = 0.0
       aminus(mza) = arw(mza-1,iw) / volt(mza,iw)

       do n = 1, n_ae_sed_spc
          do k = lpw(iw), mza
             fplus (k,n) = aplus (k) * vsed_ae(n,k)
             fminus(k,n) = aminus(k) * vsed_ae(n,k-1)
          enddo
       enddo

       do v = 1, n_ae_spc

          s = n_gc_spc + v
          n = sedi_sur(v)
          if (n > 0) then
             
             time = 0.0_r8
             eps  = 1.e-5_r8 * min( dtmin(n), dtlm(mrl) )

             do while( time + eps < dtlm(mrl) )
             
                dt8 = min( dtmin(n), dtlm(mrl) - time )
                dt  = real(dt8)

                do k = lpw(iw), mza-1
                   cgrid(k,iw,s) = cgrid(k,iw,s) + dt *  &
                                   (fplus(k,n) * cgrid(k+1,iw,s) - fminus(k,n) * cgrid(k,iw,s))
                enddo
                cgrid(mza,iw,s) = cgrid(mza,iw,s) - dt * fminus(mza,n) * cgrid(mza,iw,s)
                
                time = time + dt8

             end do

          end if
       end do

    enddo

  end subroutine aero_sedi

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine aero_sedv ( iw, cgrd, vsed_ae )

    !-----------------------------------------------------------------------
    ! Get accum. and coarse mode grav. settling vel
    ! used Binkowski`s aerosol dry deposition routine as a guide
    ! 08 Feb 13 J.Young: initial
    ! 20 Jun 14 J.Young: restructure
    ! 22 Oct 14 J.Bash:  replaced P0 with STDATMPA from CONST.EXT and shared 
    !                    variables in the asx_data_mod
    !-----------------------------------------------------------------------

    use utilio_defn      
    use aero_data           ! aero variable data
    use soa_defn            ! gas soa data
    use aeromet_data        ! Includes CONST.EXT
    use mem_grid,  only: mza, lpw
    use mem_basic, only: tair, press, rho

    implicit none

    ! Arguments
    integer, intent( in )  :: iw
    real,    intent( in )  :: cgrd( :,: )     ! cgrid subsection (layers,species)
    real,    intent( out ) :: vsed_ae( :,: )  ! settling velocities [ m/s ]

    ! Parameters
    real,    parameter :: t0 = 288.15      ! [ K ] ! starting standard surface temp.
    real,    parameter :: two3 = 2.0 / 3.0

    ! Local variables:

    real :: xxlsgac   ! log of stnd dev
    real :: xxlsgco
    real :: dgacc     ! geometric mean diameter
    real :: dgcor  
    real :: pdensac   ! particle density
    real :: pdensco

    real :: xlm       ! mean free path [ m ]
    real :: amu       ! dynamic viscosity [ kg/m/s ]

    real :: m2_wet, m2_dry
    real :: m3_wet, m3subt, m3_dry

    integer :: l, v, n, j      ! loop counters
    integer :: spc, s          ! species loop counter

    !-----------------------------------------------------------------------

    do l = lpw(iw)+1, mza

       ! Set meteorological data for the grid cell.

       airtemp = tair ( l,iw )
       airpres = press( l,iw )
       airdens = rho  ( l,iw )

       ! extract grid cell concentrations of aero species from CGRID
       ! into aerospc_conc in aero_data module

       call extract_aero( cgrd( l,: ), .true. )  ! set minimum floor

       ! extract soa concentrations from cgrid

       call extract_soa( cgrd( l,: ) )

       ! Calculate aerosol surface area to 2nd moment. 

       do s = 2, n_mode
          n = aerosrf_map( s )
          moment2_conc( s ) = cgrd( l,n ) / pi
       end do

       ! Get the geometric mean diameters and standard deviations of the
       ! "dry" size distribution (also returns MOMENT3_CONC)

       call getpar( .false., .false. ) 
                       !        ! do not fix stnd dev`s to existing value
                       ! exclude H2O and SOA from 3rd moment

       ! add contribution of water and SOA

       ! accum mode

       m3_dry = moment3_conc( 2 )
       m3subt = ( 1.0e-9 * f6pi / aerospc( ah2o_idx )%density ) &
              * aerospc_conc( ah2o_idx,2 )
       do spc = 1, n_vapor
          m3subt = m3subt + ( 1.0e-9 * f6pi / aerospc( apoc_idx )%density ) &
                 * aerospc_conc( soa_aeromap( spc ),2 )
       end do

       m3_wet = m3_dry + m3subt
       m2_dry = moment2_conc( 2 )
       m2_wet = m2_dry * ( m3_wet / m3_dry ) ** two3

       moment3_conc( 2 ) = max( conmin, m3_wet )
       moment2_conc( 2 ) = max( conmin, m2_wet )

       ! coarse mode

       m3_dry = moment3_conc( 3 )
       m3subt = ( 1.0e-9 * f6pi / aerospc( ah2o_idx )%density ) &
              * aerospc_conc( ah2o_idx,3 )
       m3_wet = m3_dry + m3subt
       m2_dry = moment2_conc( 3 )
       m2_wet = m2_dry * ( m3_wet / m3_dry ) ** two3

       moment3_conc( 3 ) = max( conmin, m3_wet )
       moment2_conc( 3 ) = max( conmin, m2_wet )

       ! Get the geometric mean diameters and standard deviations of the
       ! "wet" size distribution

       call getpar( .true., .false. )     
                      !        ! do not fix stnd dev`s to existing value
                      ! include H2O and SOA in 3rd moment

       ! Save getpar values

       xxlsgac = aeromode_sdev( 2 )
       xxlsgco = aeromode_sdev( 3 )

       dgacc   = aeromode_diam( 2 )
       dgcor   = aeromode_diam( 3 )

       pdensac = aeromode_dens( 2 )
       pdensco = aeromode_dens( 3 )
 
       ! Calculate mean free path [ m ]:

       xlm = 6.6328e-8 * STDATMPA * airtemp / ( t0 * airpres )

       ! Calculate dynamic viscosity [ kg/m/s ]:

       amu = 1.458e-6 * airtemp * sqrt( airtemp ) / ( airtemp + 110.4 )

       ! get settling velocities:

       call get_sedv ( xlm, amu, dgacc, dgcor, xxlsgac, xxlsgco, &
                       pdensac, pdensco, vsed_ae(:,l-1)          )

    end do
    
  end subroutine aero_sedv

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine get_sedv ( xlm, amu, dgacc, dgcor, xxlsgac, xxlsgco, &
                        pdensac, pdensco, vsed )

    ! Calculate settling velocity for Aitken, accumulation, and coarse modes.

    use aeromet_data   ! Includes CONST.EXT

    implicit none

    ! *** arguments

    real, intent( in ) :: xlm      ! atmospheric mean free path [m]
    real, intent( in ) :: amu      ! atmospheric dynamic viscosity [kg/(m s)]
    real, intent( in ) :: dgacc    ! accum mode geom mean diameter [m]
    real, intent( in ) :: dgcor    ! coarse mode geom mean diameter [m]
    real, intent( in ) :: xxlsgac  ! accum mode log of modal geom stnd dev`s
    real, intent( in ) :: xxlsgco  ! coarse mode log of modal geom stnd dev`s
    real, intent( in ) :: pdensac  ! avg particle dens in accum mode [kg/m**3]
    real, intent( in ) :: pdensco  ! avg particle dens in coarse mode [kg/m**3]

    real, intent( out ) :: vsed( : )    ! grav settling velocity [ m/s ]

    ! modal Knudsen numbers X bhat

    real :: bknacc   ! accumulation mode 
    real :: bkncor   ! coarse mode

    ! modal sedimentation velocities for 0th (number), 2nd (srf area), and 3rd (mass) moments

    real :: vghat0a, vghat0c
    real :: vghat2a, vghat2c
    real :: vghat3a, vghat3c

    real :: dconst2, dconst3a, dconst3c
    real :: bxlm

    real, parameter :: bhat    = 2.492 ! 2 X Constant from Cunningham slip correction
    real, parameter :: two3    = 2.0 / 3.0

    ! Scalar variables for VARIABLE standard deviations.

    real :: l2sgac, l2sgco   ! log^2( sigmag )

    real :: esac01           ! accumu mode " ** 4
    real :: esco01           ! coarse      "

    real :: esac02           ! accumu mode " ** 8
    real :: esco02           ! coarse      "

    real :: esac04           ! accumu mode " ** 16
    real :: esco04           ! coarse      "

    real :: esac05           ! accumu mode " ** 20
    real :: esco05           ! coarse      "

    real :: esac07           ! accumu mode " ** 28
    real :: esco07           ! coarse      "

    real :: esac12           ! accumu mode " ** 48    
    real :: esco12           ! coarse      "     

    real :: esac16           ! accumu mode " ** 64
    real :: esco16           ! coarse      "

    !-----------------------------------------------------------------------

    ! Calculate Knudsen numbers * bhat

    bxlm = bhat * xlm
    bknacc = bxlm / dgacc
    bkncor = bxlm / dgcor

    ! Calculate functions of variable standard deviation

    l2sgac = xxlsgac * xxlsgac
    l2sgco = xxlsgco * xxlsgco

    esac01  = exp( 0.5 * l2sgac )
    esco01  = exp( 0.5 * l2sgco )

    esac02  = esac01 * esac01
    esco02  = esco01 * esco01

    esac04  = esac02 * esac02
    esco04  = esco02 * esco02

    esac05  = esac04 * esac01
    esco05  = esco04 * esco01

    esac07  = esac05 * esac02
    esco07  = esco05 * esco02

    esac12  = esac07 * esac05
    esco12  = esco07 * esco05

    esac16  = esac12 * esac04
    esco16  = esco12 * esco04

    dconst2  = grav / ( 18.0 * amu )
    dconst3a = dconst2 * pdensac * dgacc * dgacc
    dconst3c = dconst2 * pdensco * dgcor * dgcor

    ! acc mode

    vghat0a  = dconst3a * ( esac04  + bknacc * esac01 )
    vghat2a  = dconst3a * ( esac12  + bknacc * esac05 )
    vghat3a  = dconst3a * ( esac16  + bknacc * esac07 )

    ! coarse mode

    vghat0c  = dconst3c * ( esco04  + bkncor * esco01 )
    vghat2c  = dconst3c * ( esco12  + bkncor * esco05 )
    vghat3c  = dconst3c * ( esco16  + bkncor * esco07 )

    ! settling velocities

    ! vsed of 0th moment for the number 
    vsed( vgnacc ) = vghat0a   ! accum mode
    vsed( vgncor ) = vghat0c   ! coarse mode

    ! vsed of 2nd moment for the surface area 
    vsed( vgsacc ) = vghat2a   ! accum mode
    vsed( vgscor ) = vghat2c   ! coarse mode

    ! vsed of 3rd moment for the mass 
    vsed( vgmacc ) = vghat3a   ! accum mode
    vsed( vgmcor ) = vghat3c   ! coarse mode

  end subroutine get_sedv

end module sedv_defn
