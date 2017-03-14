
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

module cgrid_conv

  implicit none

  include 'CONST.EXT'

  logical,              save :: firstime = .true.
  integer,              save :: logdev

  integer,              save :: nqae        ! number of micro-grams/m**3 species
  integer,              save :: nnae        ! number of #/m**3 species
  integer,              save :: nsae        ! number of m**2/m**3 species
  integer, allocatable, save :: qae   ( : ) ! CGRID pointer to micro-grams/m**3 species
  integer, allocatable, save :: nae   ( : ) ! CGRID pointer to #/m**3 species
  integer, allocatable, save :: sae   ( : ) ! CGRID pointer to m**2/m**3 species
  real,    allocatable, save :: molwt ( : ) ! only for "QAE" species
  real,    allocatable, save :: molwti( : ) ! only for "QAE" species

  real, parameter :: gpkg = 1.0e+03        ! g/kg
  real, parameter :: maogpkg = mwair / gpkg
  real, parameter :: gpkgoma = 1.0 / maogpkg
  real, parameter :: maoavo1000 = 1.0e+03 * mwair / avo
  real, parameter :: avooma_001 = 1.0 / maoavo1000

  public :: conv_cgrid, rev_cgrid, conv_cgrid_iw, rev_cgrid_iw, rev_cgrid_sfc_iw
  private
      

contains


  subroutine setup_conv()

    use cgrid_spcs
    use utilio_defn 

    implicit none

    character( 16 ), parameter :: pname = 'setup_conv'
    character( 96 )            :: xmsg  = ' '
    integer                    :: ios, off, s

    logdev = init3()

    if ( n_ae_spc .gt. 0 ) then

       ! create aerosol species pointers to distinguish micro-grams/m**3,
       ! #/m**3 (number density), and m**2/m**3 (surface area) species

       allocate ( qae   ( n_ae_spc ), &
                  nae   ( n_ae_spc ), &
                  sae   ( n_ae_spc ), &
                  molwt ( n_ae_spc ), &
                  molwti( n_ae_spc ), stat = ios )

       if ( ios .ne. 0 ) then
          xmsg = 'Failure allocating QAE, NAE, SAE, or MOLWT'
          call m3exit( pname, 0, 0, xmsg, xstat1 )
       end if

       nqae = 0       ! no. of micro-grams/m**3 species
       nnae = 0       ! no. of  #/m**3 species
       nsae = 0       ! no. of  m**2/m**3 species

       off = ae_strt - 1
       do s = 1, n_ae_spc
          if ( ae_spc( s )( 1:3 ) .eq. 'NUM' ) then
             nnae = nnae + 1
             nae( nnae ) = off + s
          else if ( ae_spc( s )( 1:3 ) .eq. 'SRF' ) then
             nsae = nsae + 1
             sae( nsae ) = off + s
          else
             nqae = nqae + 1
             qae   ( nqae ) = off + s
             molwt ( nqae ) = ae_molwt ( s )
             molwti( nqae ) = ae_molwti( s )
          end if
       end do

    end if

  end subroutine setup_conv


  subroutine conv_cgrid ( mrl )

    !-----------------------------------------------------------------------
    ! Function:
    !   Convert aerosol species to molar units (ppm and m**2/mol)
    !
    ! Revision History:
    !   Written by: J.Young 21 Aug 03
    !   J.Young 31 Jan 05: dyn alloc - establish both horizontal & vertical
    !                      domain specifications in one module
    !   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
    !-----------------------------------------------------------------------

    use cgrid_defn, only: cgrid
    use mem_basic,  only: rho
    use mem_grid,   only: lpw, mza
    use mem_ijtabs, only: jtab_w, jtw_prog

    implicit none

    integer, intent(in) :: mrl

    integer :: iw, k, n, v, j
    real    :: fac
    real    :: rhoi(mza)

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas and non-reactive - no conversions necessary (always ppmV)

    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       do k = lpw(iw), mza
          rhoi(k) = 1.0 / real(rho(k,iw))
       enddo

       ! micro-grams/m**3 aerosol -> ppmv
       ! (Don't divide by MGPG, then multiply by 1.0E+6: 1/MGPG = 10**-6 cancels out
       ! ppm = 10**6)

       if ( nqae .gt. 0 ) then
          do v = 1, nqae
             n = qae( v )
             fac = maogpkg * molwti(v)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * fac * rhoi(k)
             enddo
          enddo
       endif

       ! number/m**3 aerosol -> ppmv
       ! (Don't divide by MGPG, etc. See note above)

       if ( nnae .gt. 0 ) then
          do v = 1, nnae
             n = nae( v )
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * MAOAVO1000 * rhoi(k)
             enddo
          enddo
       endif

       ! aerosol surface area
       ! m**2/m**3 aerosol -> m**2/mol air

       if ( nsae .gt. 0 ) then
          do v = 1, nsae
             n = sae( v )
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * MAOGPKG * rhoi(k)
             enddo
          enddo
       endif
      
    enddo

  end subroutine conv_cgrid


  subroutine rev_cgrid ( mrl )

    !-----------------------------------------------------------------------
    ! Function:
    !   Revert non-molar mixing ratio aerosol species back to densities
    !
    ! Revision History:
    !   Written by: J.Young 21 Aug 03
    !   J.Young 31 Jan 05: dyn alloc - establish both horizontal & vertical
    !                      domain specifications in one module
    !   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
    !-----------------------------------------------------------------------

    use cgrid_defn, only: cgrid
    use mem_basic,  only: rho
    use mem_grid,   only: lpw, mza
    use mem_ijtabs, only: jtab_w, jtw_prog

    implicit none

    integer, intent(in) :: mrl

    integer :: iw, k, n, v, j
    real    :: fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas and non-reactive - no conversions necessary (always ppmV)

    do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

       ! aerosol ppmv -> micro-grams/m**3
       ! (Don't multiply by MGPG, then divide by 1.0E+6: 1/MGPG = 10**-6 cancels out
       ! ppm = 10**6)

       if ( nqae .gt. 0 ) then
          do v = 1, nqae
             n = qae( v )
             fac = gpkgoma * molwt( v )
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * fac * real(rho(k,iw))
             enddo
          enddo
       endif
     
       ! aerosol ppmv -> number/m**3
       ! (Don't multiply by MGPG, etc. See note above)

       if ( nnae .gt. 0 ) then
          do v = 1, nnae
             n = nae( v )
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * AVOOMA_001 * real(rho(k,iw))
             enddo
          enddo
       endif

       ! aerosol surface area
       ! m**2/mol air -> m**2/m**3

       if ( nsae .gt. 0 ) then
          do v = 1, nsae
             n = sae( v )
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * GPKGOMA * real(rho(k,iw))
             enddo
          enddo
       endif

    enddo

  end subroutine rev_cgrid


  subroutine conv_cgrid_iw ( iw, cngrd )

    !-----------------------------------------------------------------------
    ! Function:
    !   Convert aerosol species to molar units (ppm and m**2/mol)
    !
    ! Revision History:
    !   Written by: J.Young 21 Aug 03
    !   J.Young 31 Jan 05: dyn alloc - establish both horizontal & vertical
    !                      domain specifications in one module
    !   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
    !-----------------------------------------------------------------------

    use cgrid_defn, only: cgrid
    use mem_basic,  only: rho
    use mem_grid,   only: lpw, mza
    use cgrid_spcs, only: n_gc_spc, n_nr_spc, nr_strt

    implicit none

    integer, intent(in ) :: iw
    real,    intent(out) :: cngrd(:,:)

    integer :: k, n, v, off
    real    :: fac
    real    :: rhoi(mza)

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    do k = lpw(iw), mza
       rhoi(k) = 1.0 / real(rho(k,iw))
    enddo

    ! Gas - no conversion

    if ( n_gc_spc .gt. 0 ) then
       do v = 1, n_gc_spc
          do k = lpw(iw), mza
             cngrd(k,v) = cgrid(k,iw,v)
          enddo
       enddo
    endif

    ! micro-grams/m**3 aerosol -> ppmv
    ! (Don't divide by MGPG, then multiply by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    if ( nqae .gt. 0 ) then
       do v = 1, nqae
          n = qae( v )
          fac = maogpkg * molwti(v)
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * fac * rhoi(k)
          enddo
       enddo
    endif
     
    ! number/m**3 aerosol -> ppmv
    ! (Don't divide by MGPG, etc. See note above)

    if ( nnae .gt. 0 ) then
       do v = 1, nnae
          n = nae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * MAOAVO1000 * rhoi(k)
          enddo
       enddo
    endif

    ! aerosol surface area
    ! m**2/m**3 -> m**2/mol air

    if ( nsae .gt. 0 ) then
       do v = 1, nsae
          n = sae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * MAOGPKG * rhoi(k)
          enddo
       enddo
    endif

    ! Non-reactives - no conversion

    if ( n_nr_spc .gt. 0 ) then
       off = nr_strt - 1
       do v = 1, n_nr_spc
          do k = lpw(iw), mza
             cngrd(k,off+v) = cgrid(k,iw,off+v)
          enddo
       enddo
    endif

  end subroutine conv_cgrid_iw


  subroutine rev_cgrid_iw ( iw, cngrd )

    !-----------------------------------------------------------------------
    ! Function:
    !   Revert non-molar mixing ratio aerosol species back to densities
    !
    ! Revision History:
    !   Written by: J.Young 21 Aug 03
    !   J.Young 31 Jan 05: dyn alloc - establish both horizontal & vertical
    !                      domain specifications in one module
    !   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
    !-----------------------------------------------------------------------

    use cgrid_defn, only: cgrid
    use mem_basic,  only: rho
    use mem_grid,   only: lpw, mza
    use cgrid_spcs, only: n_gc_spc, n_nr_spc, nr_strt

    implicit none

    integer, intent(in ) :: iw
    real,    intent(out) :: cngrd(:,:)

    integer :: k, n, v, off
    real    :: fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas - no conversion

    if ( n_gc_spc .gt. 0 ) then
       do v = 1, n_gc_spc
          do k = lpw(iw), mza
             cngrd(k,v) = cgrid(k,iw,v)
          enddo
       enddo
    endif

    ! aerosol ppmv -> micro-grams/m**3
    ! (Don't multiply by MGPG, then divide by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    if ( nqae .gt. 0 ) then
       do v = 1, nqae
          n = qae( v )
          fac = gpkgoma * molwt( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * fac * real(rho(k,iw))
          enddo
       enddo
    endif
     
    ! aerosol ppmv -> number/m**3
    ! (Don't multiply by MGPG, etc. See note above)

    if ( nnae .gt. 0 ) then
       do v = 1, nnae
          n = nae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * AVOOMA_001 * real(rho(k,iw))
          enddo
       enddo
    endif

    ! aerosol surface area
    ! m**2/mol air -> m**2/m**3

    if ( nsae .gt. 0 ) then
       do v = 1, nsae
          n = sae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * GPKGOMA * real(rho(k,iw))
          enddo
       enddo
    endif

    ! Non-reactives - no conversion

    if ( n_nr_spc .gt. 0 ) then
       off = nr_strt - 1
       do v = 1, n_nr_spc
          do k = lpw(iw), mza
             cngrd(k,off+v) = cgrid(k,iw,off+v)
          enddo
       enddo
    endif

  end subroutine rev_cgrid_iw


  subroutine rev_cgrid_sfc_iw ( iw, cngrd )

    !-----------------------------------------------------------------------
    ! Function:
    !   Convert aerosol species to molar units (ppm and m**2/mol)
    !
    ! Revision History:
    !   Written by: J.Young 21 Aug 03
    !   J.Young 31 Jan 05: dyn alloc - establish both horizontal & vertical
    !                      domain specifications in one module
    !   16 Feb 11 S.Roselle: replaced I/O API include files with UTILIO_DEFN
    !-----------------------------------------------------------------------

    use cgrid_defn, only: cgrid
    use mem_basic,  only: rho
    use mem_grid,   only: lpw, lsw
    use cgrid_spcs, only: n_gc_spc, n_nr_spc, nr_strt

    implicit none

    integer, intent(in ) :: iw
    real,    intent(out) :: cngrd(:,:)

    integer :: ks, k, n, v, kb, kt, off
    real    :: fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas - no conversion

    if ( n_gc_spc .gt. 0 ) then
       do v = 1, n_gc_spc
          do ks = 1, lsw(iw)
             k  = ks + lpw(iw) - 1
             cngrd(ks,v) = cgrid(k,iw,v)
          enddo
       enddo
    endif

    ! micro-grams/m**3 aerosol -> ppmv
    ! (Don't divide by MGPG, then multiply by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    if ( nqae .gt. 0 ) then
       do v = 1, nqae
          n = qae( v )
          fac = gpkgoma * molwt(v)
          do ks = 1, lsw(iw)
             k  = ks + lpw(iw) - 1
             cngrd(ks,n) = cgrid(k,iw,n) * fac * real(rho(k,iw))
          enddo
       enddo
    endif

    ! number/m**3 aerosol -> ppmv
    ! (Don't divide by MGPG, etc. See note above)

    if ( nnae .gt. 0 ) then
       do v = 1, nnae
          n = nae( v )
          do ks = 1, lsw(iw)
             k  = ks + lpw(iw) - 1
             cngrd(ks,n) = cgrid(k,iw,n) * AVOOMA_001 * real(rho(k,iw))
          enddo
       enddo
    endif

    ! aerosol surface area
    ! m**2/m**3 aerosol -> m**2/mol air

    if ( nsae .gt. 0 ) then
       do v = 1, nsae
          n = sae( v )
          do ks = 1, lsw(iw)
             k  = ks + lpw(iw) - 1
             cngrd(ks,n) = cgrid(k,iw,n) * GPKGOMA * real(rho(k,iw))
          enddo
       enddo
    endif
      
    ! Non-reactives - no conversion

    if ( n_nr_spc .gt. 0 ) then
       off = nr_strt - 1
       do v = 1, n_nr_spc
          do ks = 1, lsw(iw)
             k  = ks + lpw(iw) - 1
             cngrd(ks,off+v) = cgrid(k,iw,off+v)
          enddo
       enddo
    endif

  end subroutine rev_cgrid_sfc_iw


end module cgrid_conv
