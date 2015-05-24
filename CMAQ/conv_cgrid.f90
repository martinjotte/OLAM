
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

  integer,              save :: nqae       ! number of micro-grams/m**3 species
  integer,              save :: nnae       ! number of #/m**3 species
  integer,              save :: nsae       ! number of m**2/m**3 species
  integer, allocatable, save :: qae  ( : ) ! CGRID pointer to micro-grams/m**3 species
  integer, allocatable, save :: nae  ( : ) ! CGRID pointer to #/m**3 species
  integer, allocatable, save :: sae  ( : ) ! CGRID pointer to m**2/m**3 species
  real,    allocatable, save :: molwt( : ) ! only for "QAE" species

  real, parameter :: gpkg = 1.0e+03        ! g/kg
  real, parameter :: maogpkg = mwair / gpkg
  real, parameter :: gpkgoma = 1.0 / maogpkg
  real, parameter :: maoavo1000 = 1.0e+03 * mwair / avo
  real, parameter :: avooma_001 = 1.0 / maoavo1000

  public :: conv_cgrid, rev_cgrid, conv_cgrid_iw, rev_cgrid_iw
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

       allocate ( qae  ( n_ae_spc ), &
                  nae  ( n_ae_spc ), &
                  sae  ( n_ae_spc ), &
                  molwt( n_ae_spc ), stat = ios )

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
             qae( nqae ) = off + s
             molwt( nqae ) = ae_molwt( s )
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
    use mem_grid,   only: mwa, lpw, mza
    use mem_ijtabs, only: jtab_w, jtw_prog

    implicit none

    integer, intent(in) :: mrl

    integer :: nspcs
    integer :: iw, k, n, v, j
    real    :: conv, fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas - no conversion

    ! micro-grams/m**3 aerosol -> ppmv
    ! (Don't divide by MGPG, then multiply by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    nspcs = nqae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = qae( v )
          fac = maogpkg / molwt(v)
          do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * fac / rho(k,iw)
             enddo
          enddo
       enddo

    endif

    ! number/m**3 aerosol -> ppmv
    ! (Don't divide by MGPG, etc. See note above)

    nspcs = nnae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = nae( v )
          do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * MAOAVO1000 / rho(k,iw)
             enddo
          enddo
       enddo

    endif

    ! m**2/m**3 aerosol -> m**2/mol air

    nspcs = nsae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = sae( v )
          do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * MAOGPKG / rho(k,iw)
             enddo
          enddo
       enddo

    endif
      
    ! Non-reactives - no conversion

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
    use mem_grid,   only: mwa, lpw, mza
    use mem_ijtabs, only: jtab_w, jtw_prog

    implicit none

    integer, intent(in) :: mrl

    integer :: nspcs
    integer :: iw, k, n, v, j
    real    :: conv, fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas - no conversion

    ! aerosol ppmv -> micro-grams/m**3
    ! (Don't multiply by MGPG, then divide by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    nspcs = nqae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = qae( v )
          fac = gpkgoma * molwt( v )
          do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * fac * rho(k,iw)
             enddo
          enddo
       enddo
         
    endif
     
    ! aerosol ppmv -> number/m**3
    ! (Don't multiply by MGPG, etc. See note above)

    nspcs = nnae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = nae( v )
          do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * AVOOMA_001 * rho(k,iw)
             enddo
          enddo
       enddo
         
    endif

    nspcs = nsae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = sae( v )
          do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
             do k = lpw(iw), mza
                cgrid(k,iw,n) = cgrid(k,iw,n) * GPKGOMA * rho(k,iw)
             enddo
          enddo
       enddo
         
    endif

    ! Non-reactives - no conversion

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
    use cgrid_spcs, only: n_gc_spc, n_nr_spc, gc_strt, nr_strt

    implicit none

    integer, intent(in ) :: iw
    real,    intent(out) :: cngrd(:,:)

    integer :: nspcs, off
    integer :: k, n, v
    real    :: conv, fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas - no conversion

    nspcs = n_gc_spc
    if ( nspcs .gt. 0 ) then

       off = gc_strt - 1
       do v = 1, nspcs
          do k = lpw(iw), mza
             cngrd(k,off+v) = cgrid(k,iw,off+v)
          enddo
       enddo

    endif

    ! aerosol ppmv -> micro-grams/m**3
    ! (Don't multiply by MGPG, then divide by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    nspcs = nqae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = qae( v )
          fac = maogpkg / molwt(v)
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * fac / rho(k,iw)
          enddo
       enddo
         
    endif
     
    ! aerosol ppmv -> number/m**3
    ! (Don't multiply by MGPG, etc. See note above)

    nspcs = nnae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = nae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * MAOAVO1000 / rho(k,iw)
          enddo
       enddo
         
    endif

    ! m**2/m**3 aerosol -> m**2/mol air

    nspcs = nsae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = sae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * MAOGPKG / rho(k,iw)
          enddo
       enddo

    endif

    ! Non-reactives - no conversion

    nspcs = n_nr_spc
    if ( nspcs .gt. 0 ) then

       off = nr_strt - 1
       do v = 1, nspcs
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
    use cgrid_spcs, only: n_gc_spc, n_nr_spc, gc_strt, nr_strt

    implicit none

    integer, intent(in ) :: iw
    real,    intent(out) :: cngrd(:,:)

    integer :: nspcs, off
    integer :: k, n, v
    real    :: conv, fac

    if ( firstime ) then
       firstime = .false.
       call setup_conv()
    endif

    ! Gas - no conversion

    nspcs = n_gc_spc
    if ( nspcs .gt. 0 ) then

       off = gc_strt - 1
       do v = 1, nspcs
          do k = lpw(iw), mza
             cngrd(k,off+v) = cgrid(k,iw,off+v)
          enddo
       enddo

    endif

    ! aerosol ppmv -> micro-grams/m**3
    ! (Don't multiply by MGPG, then divide by 1.0E+6: 1/MGPG = 10**-6 cancels out
    ! ppm = 10**6)

    nspcs = nqae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = qae( v )
          fac = gpkgoma * molwt( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * fac * rho(k,iw)
          enddo
       enddo
         
    endif
     
    ! aerosol ppmv -> number/m**3
    ! (Don't multiply by MGPG, etc. See note above)

    nspcs = nnae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = nae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * AVOOMA_001 * rho(k,iw)
          enddo
       enddo
         
    endif

    ! m**2/m**3 aerosol -> m**2/mol air

    nspcs = nsae
    if ( nspcs .gt. 0 ) then

       do v = 1, nspcs
          n = sae( v )
          do k = lpw(iw), mza
             cngrd(k,n) = cgrid(k,iw,n) * GPKGOMA * rho(k,iw)
          enddo
       enddo

    endif

    ! Non-reactives - no conversion

    nspcs = n_nr_spc
    if ( nspcs .gt. 0 ) then

       off = nr_strt - 1
       do v = 1, nspcs
          do k = lpw(iw), mza
             cngrd(k,off+v) = cgrid(k,iw,off+v)
          enddo
       enddo

    endif

  end subroutine rev_cgrid_iw

end module cgrid_conv
