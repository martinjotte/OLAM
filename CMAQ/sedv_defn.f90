
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
    use oname_coms,  only: nl

    implicit none

    integer, intent(in) :: mrl

    integer :: j, iw, v, n
    integer :: astat

    character( 120 )           :: xmsg = ' '
    character( 16 ), parameter :: pname = 'aero_sedi'

    if (nl%do_aesedi == 0) return

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

       allocate(sedi_sfc( n_ae_spc, mwa )) ; sedi_sfc = 0.0

    endif  ! firstime

    !$omp parallel do private(iw)
    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       call aero_sedi2( iw )

    enddo
    !$omp end parallel do

  end subroutine aero_sedi

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  subroutine aero_sedi2( iw )

    use mem_grid,    only: mza, lpw, volt, volti, arw, voa0, &
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

    do k = lpw(iw), lpw(iw) + lsw(iw) - 1
       arp(k-1) = arw(k,iw) * zfacim2(k) * zfacm2(k-1)
       voa(k)   = real(volt(k,iw)) / arp(k-1)
    enddo

    do k = lpw(iw) + lsw(iw), mza
       arp(k-1) = arw (k-1,iw)
       voa(k)   = voa0(k)
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
             cgrid(k,iw,s) = cgrid(k,iw,s) + real(volti(k,iw)) &
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

    use aero_data,  only: aer_str, aer_end, aer_num, aer_m3, mode_map, aer_trac, &
                          n_mode, aeromode, m0_min, m2_min, num_str, srf_str
    use const_data, only: stdatmpa, grav, oneovpi
    use mem_grid,   only: mza, lpw
    use mem_basic,  only: tair, press
    use cgrid_defn, only: cgrid
    use getpar_mod, only: getpar, getdens_conc
    use cgrid_spcs, only: nspcsd, ae_strt, ae_fini

    implicit none

    ! Arguments
    integer, intent( in )  :: iw
    real,    intent( out ) :: vsed_ae( mza,n_ae_sed_spc )  ! settling velocities [ m/s ]

    ! Parameters
    real,    parameter :: t0   = 288.15      ! [ K ] ! starting standard surface temp.
    real,    parameter :: one3 = 1.0 / 3.0
    real,    parameter :: bhat = 2.492       ! 2 X Constant from Cunningham slip correction

    ! Local variables:

    real :: xxlsgac(mza)   ! log of stnd dev
    real :: xxlsgco(mza)
    real :: dgacc  (mza)   ! geometric mean diameter
    real :: dgcor  (mza)
    real :: pdensac(mza)   ! particle density
    real :: pdensco(mza)

    real :: airtemp
    real :: airpres

    real :: xlm
    real :: amu

    integer :: k, i, n, s  ! loop counters

    real :: conc(nspcsd)

    ! modal Knudsen numbers X bhat

    real :: bknacc   ! accumulation mode
    real :: bkncor   ! coarse mode

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

    real :: aeromode_mass(n_mode)
    real :: aeromode_dens(n_mode)
    real :: aeromode_lnsg(n_mode)
    real :: aeromode_diam(n_mode)
    real :: moment0_conc (n_mode)
    real :: moment2_conc (n_mode)
    real :: moment3_conc (n_mode)
    real :: m3           (aer_num)

    !-----------------------------------------------------------------------

    do k = lpw(iw), mza

       ! extract grid cell concentrations of aero species from CGRID

       conc(ae_strt:ae_fini) = cgrid(k,iw,ae_strt:ae_fini)

       ! compute aerosol moments directly from concentration array

       m3 = conc(aer_str:aer_end) * aer_m3

       do i = 2, 3
          n = num_str + i - 1
          s = srf_str + i - 1
          moment3_conc( i ) = Max( sum(m3, mask=(mode_map==i .and. .not. aer_trac) ), &
                                   aeromode( i )%min_m3conc )
          moment0_conc( i ) = Max( conc( n ), m0_min( i ) )
          moment2_conc( i ) = Max( conc( s ) * oneovpi, m2_min( i ) )
       enddo

       ! compute aerosol diameters and geometric standard deviation

       call getpar( .false., moment0_conc, moment2_conc, moment3_conc, &
                     aeromode_lnsg, aeromode_diam, 2, 3 )

       ! Compute mean particle densities

       call getdens_conc(aeromode_mass, aeromode_dens, conc, moment3_conc, 2, 3)

       ! Save getpar values

       xxlsgac(k) = aeromode_lnsg( 2 )
       xxlsgco(k) = aeromode_lnsg( 3 )

       dgacc(k)   = aeromode_diam( 2 )
       dgcor(k)   = aeromode_diam( 3 )

       pdensac(k) = aeromode_dens( 2 )
       pdensco(k) = aeromode_dens( 3 )
    enddo


    do k = lpw(iw), mza

       ! Set meteorological data for the grid cell.

       airtemp = tair ( k,iw )
       airpres = press( k,iw )

       ! Calculate mean free path [ m ]:

       xlm = 6.6328e-8 * STDATMPA * airtemp / ( t0 * airpres )

       ! Calculate dynamic viscosity [ kg/m/s ]:

       amu = 1.458e-6 * airtemp * sqrt( airtemp ) / ( airtemp + 110.4 )

       ! Calculate Knudsen numbers * bhat

       bxlm = bhat * xlm
       bknacc = bxlm / dgacc(k)
       bkncor = bxlm / dgcor(k)

       ! Calculate functions of variable standard deviation

       l2sgac = xxlsgac(k) * xxlsgac(k)
       l2sgco = xxlsgco(k) * xxlsgco(k)

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
       dconst3a = dconst2 * pdensac(k) * dgacc(k) * dgacc(k)
       dconst3c = dconst2 * pdensco(k) * dgcor(k) * dgcor(k)

       ! acc mode settling velocities

       vsed_ae( k,vgnacc ) = dconst3a * ( esac04  + bknacc * esac01 ) ! 0th moment for number
       vsed_ae( k,vgsacc ) = dconst3a * ( esac12  + bknacc * esac05 ) ! 2nd moment for area
       vsed_ae( k,vgmacc ) = dconst3a * ( esac16  + bknacc * esac07 ) ! 3rd moment for mass

       ! coarse mode settling velocities

       vsed_ae( k,vgncor ) = dconst3c * ( esco04  + bkncor * esco01 ) ! 0th moment for number
       vsed_ae( k,vgscor ) = dconst3c * ( esco12  + bkncor * esco05 ) ! 2nd moment for area
       vsed_ae( k,vgmcor ) = dconst3c * ( esco16  + bkncor * esco07 ) ! 3rd moment for mass

    end do

  end subroutine aero_sedv

end module sedv_defn
