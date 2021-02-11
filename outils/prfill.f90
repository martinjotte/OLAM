module prfill_mod

  use consts_coms, only: r8, i1

  private :: r8, i1

  interface prfill
     module procedure prfill_int
     module procedure prfill_int1
     module procedure prfill_real
     module procedure prfill_real8
  end interface prfill

  interface prfill3
     module procedure prfill3_int
     module procedure prfill3_int1
     module procedure prfill3_real
     module procedure prfill3_real8
  end interface prfill3

contains

!===============================================================================

  subroutine prfill_int (nprx, npry, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer, intent(in)    :: nprx, npry
    integer, intent(in)    :: xx(nprx,npry)
    integer, intent(inout) :: dn(nprx+4,npry+4)

    real,    intent(in)    :: gdatdy, xswlat
    integer, intent(in)    :: ipoffset, inproj

    include 'prfill_body.inc'

  end subroutine prfill_int


  subroutine prfill_int1 (nprx, npry, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer,     intent(in)    :: nprx, npry
    integer(i1), intent(in)    :: xx(nprx,npry)
    integer(i1), intent(inout) :: dn(nprx+4,npry+4)

    real,        intent(in)    :: gdatdy, xswlat
    integer,     intent(in)    :: ipoffset, inproj

    include 'prfill_body.inc'

  end subroutine prfill_int1


  subroutine prfill_real (nprx, npry, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer, intent(in)    :: nprx, npry
    real,    intent(in)    :: xx(nprx,npry)
    real,    intent(inout) :: dn(nprx+4,npry+4)

    real,    intent(in)    :: gdatdy, xswlat
    integer, intent(in)    :: ipoffset, inproj

    include 'prfill_body.inc'

  end subroutine prfill_real


  subroutine prfill_real8 (nprx, npry, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer, intent(in)    :: nprx, npry
    real(8), intent(in)    :: xx(nprx,npry)
    real(8), intent(inout) :: dn(nprx+4,npry+4)

    real,    intent(in)    :: gdatdy, xswlat
    integer, intent(in)    :: ipoffset, inproj

    include 'prfill_body.inc'

  end subroutine prfill_real8

!===============================================================================

  subroutine prfill3_int (nprx, npry, nprz, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer, intent(in)    :: nprx, npry, nprz

    integer, intent(in)    :: xx(nprx,npry,nprz)
    integer, intent(inout) :: dn(nprx+4,npry+4,nprz)

    real,    intent(in)    :: gdatdy, xswlat
    integer, intent(in)    :: ipoffset, inproj

    include 'prfill3_body.inc'

  end subroutine prfill3_int


  subroutine prfill3_int1 (nprx, npry, nprz, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer,     intent(in)    :: nprx, npry, nprz

    integer(i1), intent(in)    :: xx(nprx,npry,nprz)
    integer(i1), intent(inout) :: dn(nprx+4,npry+4,nprz)

    real,        intent(in)    :: gdatdy, xswlat
    integer,     intent(in)    :: ipoffset, inproj

    include 'prfill3_body.inc'

  end subroutine prfill3_int1


  subroutine prfill3_real (nprx, npry, nprz, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer, intent(in)    :: nprx, npry, nprz

    real,    intent(in)    :: xx(nprx,npry,nprz)
    real,    intent(inout) :: dn(nprx+4,npry+4,nprz)

    real,    intent(in)    :: gdatdy, xswlat
    integer, intent(in)    :: ipoffset, inproj

    include 'prfill3_body.inc'

  end subroutine prfill3_real


  subroutine prfill3_real8 (nprx, npry, nprz, xx, dn, gdatdy, xswlat, ipoffset, inproj)

    implicit none

    integer, intent(in)    :: nprx, npry, nprz

    real(8), intent(in)    :: xx(nprx,npry,nprz)
    real(8), intent(inout) :: dn(nprx+4,npry+4,nprz)

    real,    intent(in)    :: gdatdy, xswlat
    integer, intent(in)    :: ipoffset, inproj

    include 'prfill3_body.inc'

  end subroutine prfill3_real8

end module prfill_mod
