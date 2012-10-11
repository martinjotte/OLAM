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

end module tridiag
