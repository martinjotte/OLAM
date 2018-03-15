
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
                        vgmcor = 6    ! coarse mode mass

  ! follow the Namelist dep vel surrogate name table
  character( 16 ), parameter :: vgae_name( n_ae_sed_spc ) = &
                                (/ 'VNUMACC', &
                                   'VNUMCOR', &
                                   'VSRFACC', &
                                   'VSRFCOR', &
                                   'VMASSJ ', &
                                   'VMASSC '  /)

  integer, allocatable, save :: sedi_sur( : )   ! pointer to surrogate
  real,    allocatable, save :: sedi_sfc( :,: ) ! kg/m^s deposited this timestep

  logical, save :: firstime = .true.
  integer, save :: logdev

  private
  public :: aero_sedi

contains

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine aero_sedi( mrl )

    use mem_ijtabs,  only: jtab_w, jtw_prog
    use mem_grid,    only: mwa
    use utilio_defn, only: index1, xstat1, m3exit, init3
    use cgrid_spcs,  only: n_ae_spc, n_ae_depv, ae_depv

    implicit none

    integer, intent(in) :: mrl

    integer :: j, iw, v, n
    integer :: astat

    character( 120 )           :: xmsg = ' '
    character( 16 ), parameter :: pname = 'aero_sedi'

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

       allocate(sedi_sfc( n_ae_spc, mwa )) ; sedi_sfc( :,: ) = 0.0

    endif  ! firstime

    !$omp parallel do private(iw)
    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       call aero_sedi2( iw )

    enddo
    !$omp end parallel do

  end subroutine aero_sedi

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine aero_sedi2( iw )

    use mem_grid,    only: mza, lpw, volt, volti, arw, &
                           lsw, zm, zfacm2, zfacim2
    use mem_ijtabs,  only: itab_w
    use misc_coms,   only: dtlm
    use cgrid_spcs,  only: n_gc_spc, n_ae_spc
    use cgrid_defn,  only: cgrid
    use mem_turb,    only: frac_sfc

    implicit none

    integer, intent(in) :: iw

    integer :: s, v, n, k, kk, ks
    real    :: dt
    real    :: vsed_ae( mza, n_ae_sed_spc )
    real    :: voa(mza), arp(mza)
    logical :: only1( n_ae_sed_spc )
    real    :: zbotnew(mza, n_ae_sed_spc)
    real    :: aeflux(mza) 
    real    :: fracwkk, areascale

    call aero_sedv( iw, vsed_ae )

    dt = dtlm(itab_w(iw)%mrlw)

    ! Ratio of grid cell volume to top horizontal area arw projected onto W(k-1) level

    do k = lpw(iw) - 1, lpw(iw) + lsw(iw) - 2
       arp(k) = arw(k+1,iw) * zfacim2(k+1) * zfacm2(k)
    enddo

    do k = lpw(iw) + lsw(iw) - 1, mza
       arp(k) = arw(k,iw)
    enddo

    do k = lpw(iw), mza
       voa(k) = volt(k,iw) / arp(k-1)
    enddo

    ! New bottom height of source-cell aerosol after fall for 1 time step

    do n = 1, n_ae_sed_spc
       do k = lpw(iw), mza
          zbotnew(k,n) = zm(k-1) - dt * vsed_ae(k,n)
       enddo
       only1(n) = all( zbotnew(lpw(iw)+1:mza,n) >= zm(lpw(iw)-1:mza-2) )
    enddo

    aeflux(mza) = 0.0

    do v = 1, n_ae_spc
       s = n_gc_spc + v
       n = sedi_sur(v)

       if (n > 0) then

          if (only1(n)) then

             do k = lpw(iw)-1, mza-1

                ! Depth of source grid cell that crosses (k) level
                fracwkk = min(voa(k+1), zm(k) - zbotnew(k+1,n))

                ! Flux across (k) level
                aeflux(k) = cgrid(k+1,iw,s) * fracwkk
             enddo

          else

             aeflux( lpw(iw)-1:mza-1 ) = 0.0

             do k = lpw(iw), mza
             
                do kk = k-1, lpw(iw)-1, -1
                   if (zm(kk) <= zbotnew(k,n)) exit
                
                   ! Horizontal area scale factor for current (kk) level
                   areascale = zfacim2(k-1) * zfacm2(kk)

                   ! Depth of source grid cell that crosses (kk) level
                   fracwkk = min(voa(k), zm(kk) - zbotnew(k,n))

                   ! Add source cell contribution to flux
                   aeflux(kk) = aeflux(kk) + cgrid(k,iw,s) * fracwkk * areascale
                enddo

             enddo

          endif
    
          ! Apply fluxes to obtain new aerosol concentrations

          do k = lpw(iw), mza
             cgrid(k,iw,s) = cgrid(k,iw,s) + volti(k,iw) &
                  * (aeflux(k) * arw(k,iw) - aeflux(k-1) * arp(k-1))
          enddo

          ! Save total deposited to ground (per m^2)

          sedi_sfc(v,iw) = 0.0
          do ks = 1, lsw(iw)
             k  = ks + lpw(iw) - 1
             sedi_sfc(v,iw) = sedi_sfc(v,iw) + aeflux(k-1) * frac_sfc(ks,iw)
          enddo

       endif
    enddo

  end subroutine aero_sedi2

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine aero_sedv ( iw, vsed_ae )

    !-----------------------------------------------------------------------
    ! Get accum. and coarse mode grav. settling vel
    ! used Binkowski`s aerosol dry deposition routine as a guide
    ! 08 Feb 13 J.Young: initial
    ! 20 Jun 14 J.Young: restructure
    ! 22 Oct 14 J.Bash:  replaced P0 with STDATMPA from CONST.EXT and shared 
    !                    variables in the asx_data_mod
    !-----------------------------------------------------------------------

    use aero_data     ! aero variable data
    use const_data, only: stdatmpa, grav
    use mem_grid,   only: mza, lpw
    use mem_basic,  only: tair, press
    use cgrid_defn, only: cgrid
    use getpar_mod, only: getpar

    implicit none

    ! Arguments
    integer, intent( in )  :: iw
    real,    intent( out ) :: vsed_ae( mza,n_ae_sed_spc )  ! settling velocities [ m/s ]

    ! Parameters
    real,    parameter :: t0   = 288.15      ! [ K ] ! starting standard surface temp.
    real,    parameter :: one3 = 1.0 / 3.0
    real,    parameter :: bhat = 2.492       ! 2 X Constant from Cunningham slip correction

    ! Local variables:

    real :: xxlsgac   ! log of stnd dev
    real :: xxlsgco
    real :: dgacc     ! geometric mean diameter
    real :: dgcor  
    real :: pdensac   ! particle density
    real :: pdensco

    real :: airtemp
    real :: airpres

    real :: xlm       ! mean free path [ m ]
    real :: amu       ! dynamic viscosity [ kg/m/s ]

    integer :: l      ! loop counters

    ! modal Knudsen numbers X bhat

    real :: bknacc   ! accumulation mode 
    real :: bkncor   ! coarse mode

    ! modal sedimentation velocities for 0th (number), 2nd (srf area), and 3rd (mass) moments

    real :: vghat0a, vghat0c
    real :: vghat2a, vghat2c
    real :: vghat3a, vghat3c

    real :: dconst2, dconst3a, dconst3c
    real :: bxlm

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

    do l = lpw(iw), mza

       ! Set meteorological data for the grid cell.

       airtemp = tair ( l,iw )
       airpres = press( l,iw )

       ! extract grid cell concentrations of aero species from CGRID
       ! into aerospc_conc in aero_data module
       ! Also converts dry surface area to wet 2nd moment

       call extract_aero( cgrid( l,iw,: ), .true. )  ! set minimum floor

       ! Get the geometric mean diameters and standard deviations of the
       ! "wet" size distribution

       call getpar( .false., noM3=.true. )
                       ! do not fix stnd dev`s to existing value

       ! Save getpar values

       xxlsgac = aeromode_lnsg( 2 )
       xxlsgco = aeromode_lnsg( 3 )

       dgacc   = aeromode_diam( 2 )
       dgcor   = aeromode_diam( 3 )

       pdensac = aeromode_dens( 2 )
       pdensco = aeromode_dens( 3 )

       ! Calculate mean free path [ m ]:

       xlm = 6.6328e-8 * STDATMPA * airtemp / ( t0 * airpres )

       ! Calculate dynamic viscosity [ kg/m/s ]:

       amu = 1.458e-6 * airtemp * sqrt( airtemp ) / ( airtemp + 110.4 )

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
       vsed_ae( l,vgnacc ) = vghat0a   ! accum mode
       vsed_ae( l,vgncor ) = vghat0c   ! coarse mode

       ! vsed of 2nd moment for the surface area 
       vsed_ae( l,vgsacc ) = vghat2a   ! accum mode
       vsed_ae( l,vgscor ) = vghat2c   ! coarse mode

       ! vsed of 3rd moment for the mass 
       vsed_ae( l,vgmacc ) = vghat3a   ! accum mode
       vsed_ae( l,vgmcor ) = vghat3c   ! coarse mode

    end do
    
  end subroutine aero_sedv

end module sedv_defn
