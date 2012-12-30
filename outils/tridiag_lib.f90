module tridiag
  
contains

!===========================================================================

  subroutine tridiffo(m1,ka,kz,cim1,ci,cip1,rhs,soln)

    implicit none

    ! SERIAL TRIDIAGONAL SOLVER

    integer, intent(in) :: m1,ka,kz
    real, intent(in) :: cim1(m1),ci(m1),cip1(m1),rhs(m1)
    real, intent(out) :: soln(m1)
    real :: scr1(m1)  ! automatic array
    real :: scr2(m1)  ! automatic array

    integer :: k
    real    :: cji

    scr1(ka) = cip1(ka) / ci(ka)
    scr2(ka) = rhs(ka) / ci(ka)

    do k = ka+1,kz
       soln(k) = ci(k) - cim1(k) * scr1(k-1)
       cji = 1. / soln(k)
       scr1(k) = cip1(k) * cji
       scr2(k) = (rhs(k) - cim1(k) * scr2(k-1)) * cji
    enddo
    
    soln(kz) = scr2(kz)

    do k = kz-1,ka,-1
       soln(k) = scr2(k) - scr1(k) * soln(k+1)
    enddo

  end subroutine tridiffo

!===========================================================================

  subroutine tridv ( l, d, u, b, x, ka, kz, nz, nsp )

!   Solves tridiagonal system by Thomas algorithm for multiple input vectors
!
!   The associated tri-diagonal system is stored in 3 arrays:
!   D : diagonal
!   L : sub-diagonal
!   U : super-diagonal
!
!   B : right hand side for multiple vectors
!   X : return solution from tridiagonal solver

!     [ D(1) U(1) 0    0    0 ...       0     ]
!     [ L(2) D(2) U(2) 0    0 ...       .     ]
!     [ 0    L(3) D(3) U(3) 0 ...       .     ]
!     [ .       .     .     .           .     ] X(i) = B(i)
!     [ .             .     .     .     0     ]
!     [ .                   .     .     .     ]
!     [ 0                           L(n) D(n) ]

!-----------------------------------------------------------------------

    implicit none
      
! Arguments:
    
    integer, intent(in)  :: ka, kz, nz
    integer, intent(in)  :: nsp

    real,    intent(in)  :: l(nz)      ! subdiagonal
    real,    intent(in)  :: d(nz)      ! diagonal
    real,    intent(in)  :: u(nz)      ! superdiagonal
    real,    intent(in)  :: b(nz,nsp)  ! r.h. side
    real,    intent(out) :: x(nz,nsp)  ! solution

! Local Variables:

    real    ::  gam(kz)
    real    ::  bet
    integer ::  v, k

! Decomposition and forward substitution:

    bet = 1.0 / d( ka )
    do v = 1, nsp
       x( ka,v ) = bet * b(ka,v)
    enddo

    do k = ka+1, kz
       gam(k) = bet * u( k-1 )
       bet = 1.0 / ( d( k ) - l( k ) * gam( k ) )
       do v = 1, nsp
          x( k,v ) = bet * ( b( k,v ) - l( k ) * x( k-1,v ) )
       enddo
    enddo

! Back-substitution:

    do v = 1, nsp
       do k = kz - 1, ka, -1
          x( k,v ) = x( k,v ) - gam( k+1 ) * x( k+1,v )
       enddo
    enddo
     
  end subroutine tridv

!===========================================================================

  subroutine acm_matrix ( a, b, c, d, e, x, nlays, nspcs, m )

!-- Bordered band diagonal matrix solver for ACM2

!-- ACM2 Matrix is in this form:
!   B1  E1
!   C2  B2  E2
!   A31 C3  B3  E3
!   A41 A42 C4  B4 E4
!   A51 A52 A5M C5 B5 E5
!   A61 A62 A5M    C6 B6

!--Upper Matrix is
!  U11 U12
!      U22 U23
!          U33 U34
!              U44 U45
!                  U55 U56
!                      U66

!--Lower Matrix is:
!  1
! L21  1
! L31 L32  1
! L41 L42 L43  1
! L51 L52 L53 L54  1
! L61 L62 L63 L64 L65 1

    implicit none
    
! Arguments:

    real,    intent( in  ) :: a( :,: )   ! left matrix columns
    real,    intent( in  ) :: b( : )     ! diagonal
    real,    intent( in  ) :: c( : )     ! subdiagonal
    real,    intent( in  ) :: e( : )     ! superdiagonal
    real,    intent( in  ) :: d( :,: )   ! R.H.S
    real,    intent( out ) :: x( :,: )   ! returned solution
    
    integer, intent( in  ) :: nlays      ! number of model layers
    integer, intent( in  ) :: nspcs      ! number of species
    integer, intent( in  ) :: m          ! number of A vectors
    
! Locals:

    real :: y  ( nlays, nspcs )
    real :: l  ( nlays, nlays )
    real :: u  ( nlays )
    real :: up1( nlays )
    real :: ru ( nlays )

    real    :: dd, dd1, ysum
    integer :: i, j, v, jj

!-- Define Upper and Lower matrices

    l( 1,1 ) = 1.0
    u( 1 ) = b( 1 )
    ru( 1 ) = 1.0 / b( 1 )

    l( 2,1 ) = c(2) / b( 1 )

    do i = 2, nlays
       l( i,i ) = 1.0
       up1( i-1 ) = e( i-1 )
    end do

    do i = 3, nlays
       l( i,1 ) = a( i,1 ) / b( 1 )
    enddo

    do i = 3, nlays

       jj = min(m, i-2)

       do j = 2, jj
          dd = b( j ) - l( j,j-1 ) * e( j-1 )
          l( i,j ) = ( a( i,j ) - l( i,j-1 ) * e( j-1 ) ) / dd
       end do

       do j = jj + 1, i - 2
          dd = b( j ) - l( j,j-1 ) * e( j-1 )
          l( i,j ) = - l( i,j-1 ) * e( j-1 ) / dd
       end do

       j = i - 1
       dd = b( j ) - l( j,j-1 ) * e( j-1 )
       l( i,j ) = ( c( i ) - l( i,j-1 ) * e( j-1 ) ) / dd
    end do

    do i = 2, nlays
       u( i ) = b( i ) - l( i,i-1 ) * e( i-1 )
       ru( i ) = 1.0 / u( i )
    end do

!-- Forward sub for Ly=d

    do v = 1, nspcs
       y( 1,v ) = d( 1,v )
       do i = 2, nlays
          ysum = d( i,v )
          do j = 1, i-1
             ysum = ysum - l( i,j ) * y( j,v )
          end do
          y( i,v ) = ysum
       end do
    end do
      
! -- Back sub for Ux=y

    do v= 1, nspcs
       x( nlays,v ) = y( nlays,v ) * ru( nlays )
    end do

    do i = nlays - 1, 1, -1
       dd = ru( i )
       dd1 = up1( i )
       do v = 1, nspcs
          x( i,v ) = ( y( i,v ) - dd1 * x( i+1,v ) ) * dd
       end do
    end do

  end subroutine acm_matrix

end module tridiag
