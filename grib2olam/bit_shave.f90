module bit_shave

contains

!************************************************************************

subroutine scalfield(a,nx,ny,nsd,rmiss)

  implicit none

  integer, intent(in)           :: nx, ny, nsd
  real,    intent(inout)        :: a(nx,ny)
  real,    intent(in), optional :: rmiss
  integer                       :: i, j

  if (present(rmiss)) then

     do j = 1, ny
        do i = 1, nx
           if (a(i,j) /= rmiss) then
              call shave_bits( a(i,j), nsd )
           endif
        enddo
     enddo

  else

     do j = 1, ny
        do i = 1, nx
           call shave_bits( a(i,j), nsd )
        enddo
     enddo

  endif

end subroutine scalfield

!************************************************************************

subroutine shave_bits(x,nsd)

  implicit none

  integer, parameter :: r4 = selected_real_kind(6,37)
  integer, parameter :: i4 = selected_int_kind(8)

  ! Bit shaving algorithm to set floating bits to zero depending on
  ! the desired numerical precision.
  ! This allows the HDF5 shuffle and deflate compression algorithms
  ! to better compress the datasets

  real(r4), intent(inout) :: x
  integer,  intent(in)    :: nsd ! number of significant digits preserved
  integer                 :: j
  integer                 :: bits_masked
  integer(i4)             :: i
  real,     parameter     :: bits_per_decimal_prec = log(10.) / log(2.0)
  integer,  parameter     :: bits_preserved_arr(6) = ceiling(                  &
                                 (/1.,2.,3.,4.,5.,6./) * bits_per_decimal_prec )

  if (nsd < 1) return
  if (nsd > 6) return

  ! skip 0
  if (x == 0._r4) return

  ! skip 2**i (..., 0.25, 0.5, 1.0, 2.0, 4.0, ...)
  if (fraction(x) == 0.5) return

  ! real*4 has 23 precision bits
  bits_masked = 23 - bits_preserved_arr(nsd)

  ! Add half the magnitude of the largest shaved bit to the floating
  ! point value before shaving so that values are rounded properly.
  ! Without this, values will always be rounded down.
  x = x + sign( 2.**(exponent(x) - bits_preserved_arr(nsd) - 2), x )

  ! Fortran bit operations only work on integers
  i = transfer(x,1_i4)

  ! Shave the desired number of bits
  do j = 0, bits_masked-1
     i = ibclr(i,j)
  enddo

  ! Transfer back to real*4 type
  x = transfer(i,1._r4)

end subroutine shave_bits

!************************************************************************

end module bit_shave
