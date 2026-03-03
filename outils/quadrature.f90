module quadrature

  private
  public :: hex_quad

contains

!===============================================================================

  subroutine hex_quad(iw, n, fint, func)

    use mem_grid,    only: xem, yem, zem, glatw, glonw, xew, yew
    use mem_ijtabs,  only: itab_w
    use misc_coms,   only: mdomain
    use consts_coms, only: r8
    use map_proj,    only: ec_ps

    implicit none

    integer,  intent(in)  :: iw, n
    real(r8), intent(out) :: fint

    interface
       real(r8) function func(x,y)
         use consts_coms, only: r8
         real(r8), intent(in) :: x
         real(r8), intent(in) :: y
       end function func
    end interface

    real     :: xm(7), ym(7)
    real(r8) :: u1, v1, u2, v2, u3, v3
    real(r8) :: tri_int
    integer  :: im, jm, npoly, jm2, jm3

    npoly = itab_w(iw)%npoly

    ! Determine local PS coordinates of M points

    do jm = 1,npoly
       im = itab_w(iw)%im(jm)

       if (mdomain <= 1) then
          call ec_ps(xem(im),yem(im),zem(im),glatw(iw),glonw(iw),xm(jm),ym(jm))
       else
          xm(jm) = xem(im) - xew(iw)
          ym(jm) = yem(im) - yew(iw)
       endif
    enddo

    ! Determine integral over each sub triangle of hexagon

    u1   = 0.0_r8
    v1   = 0.0_r8
    fint = 0.0_r8

    jm3 = npoly
    do jm2 = 1,npoly

       u2 = xm(jm2)
       v2 = ym(jm2)

       u3 = xm(jm3)
       v3 = ym(jm3)

       call dtria(n, u1, v1, u2, v2, u3, v3, tri_int, func)
       fint = fint + tri_int

       jm3 = jm2
    enddo

  end subroutine hex_quad

!===============================================================================

  subroutine dtria(n, u1, v1, u2, v2, u3, v3, approx, f)

    use consts_coms, only: r8
    implicit none

    integer,  intent( in) :: n
    real(r8), intent( in) :: u1, v1, u2, v2, u3, v3
    real(r8), intent(out) :: approx

    interface
       real(r8) function f(x,y)
         use consts_coms, only: r8
         real(r8), intent(in) :: x
         real(r8), intent(in) :: y
       end function f
    end interface

    real(r8) :: absc(n), wght(n)
    real(r8) :: u21, u31, v21, v31
    real(r8) :: save1, save2, save3, save4, save5, save6
    real(r8) :: x, y, xx, yy, absjac, temp
    integer  :: nnp1, itest, l, m

    ! abscissas and weights for n = 2

    real(r8) :: x2(2) = [ 0.577350269189626_r8, -0.577350269189626_r8  ]
    real(r8) :: w2(2) = [ 1.000000000000000_r8,  1.000000000000000_r8  ]

    ! abscissas and weights for n = 3

    real(r8) :: x3(3) = [ 0.774596669241483_r8,  0.000000000000000_r8, &
                         -0.774596669241483_r8 ]

    real(r8) :: w3(3) = [ 0.555555555555555_r8,  0.888888888888889_r8, &
                          0.555555555555555_r8 ]

    ! abscissas and weights for n = 4

    real(r8) :: x4(4) = [ 0.861136311594053_r8,  0.339981043584856_r8, &
                         -0.339981043584856_r8, -0.861136311594053_r8  ]

    real(r8) :: w4(4) = [ 0.347854845137454_r8,  0.652145154862546_r8, &
                          0.652145154862546_r8,  0.347854845137454_r8  ]

    ! abscissas and weights for n = 5

    real(r8) :: x5(5) = [ 0.906179845938664_r8,  0.538469310105683_r8, &
                          0.000000000000000_r8, -0.538469310105683_r8, &
                         -0.906179845938664_r8 ]

    real(r8) :: w5(5) = [ 0.236926885056189_r8,  0.478628670499366_r8, &
                          0.568888888888889_r8,  0.478628670499366_r8, &
                          0.236926885056189_r8 ]

    ! abscissas and weights for n = 6

    real(r8) :: x6(6) = [ 0.932469514203152_r8,  0.661209386466264_r8, &
                          0.238619186083197_r8, -0.238619186083197_r8, &
                         -0.661209386466264_r8, -0.932469514203152_r8  ]

    real(r8) :: w6(6) = [ 0.171324492379170_r8,  0.360761573048139_r8, &
                          0.467913934572691_r8,  0.467913934572691_r8, &
                          0.360761573048139_r8,  0.171324492379170_r8  ]

    ! abscissas and weights for n = 7

    real(r8) :: x7(7) = [ 0.949107912342758_r8,  0.741531185599394_r8, &
                          0.405845151377397_r8,  0.000000000000000_r8, &
                         -0.405845151377397_r8, -0.741531185599394_r8, &
                         -0.949107912342758_r8 ]

    real(r8) :: w7(7) = [ 0.129484966168870_r8,  0.279705391489277_r8, &
                          0.381830050505119_r8,  0.417959183673469_r8, &
                          0.381830050505119_r8,  0.279705391489277_r8, &
                          0.129484966168870_r8 ]

    ! abscissas and weights for n = 8

    real(r8) :: x8(8) = [ 0.960289856497536_r8,  0.796666477413627_r8, &
                          0.525532409916329_r8,  0.183434642495650_r8, &
                         -0.183434642495650_r8, -0.525532409916329_r8, &
                         -0.796666477413627_r8, -0.960289856497536_r8  ]

    real(r8) :: w8(8) = [ 0.101228536290376_r8,  0.222381034453374_r8, &
                          0.313706645877887_r8,  0.362683783378362_r8, &
                          0.362683783378362_r8,  0.313706645877887_r8, &
                          0.222381034453374_r8,  0.101228536290376_r8  ]

    nnp1   = n + 1
    itest  = nnp1/2
    approx = 0.0_r8

    ! GET QUADRATURE ABSCISSAS and WEIGHTS

    if     (n == 2) then
       absc = x2
       wght = w2
    elseif (n == 3) then
       absc = x3
       wght = w3
    elseif (n == 4) then
       absc = x4
       wght = w4
    elseif (n == 5) then
       absc = x5
       wght = w5
    elseif (n == 6) then
       absc = x6
       wght = w6
    elseif (n == 7) then
       absc = x7
       wght = w7
    elseif (n == 8) then
       absc = x8
       wght = w8
    else
       call grule(n, absc, wght)
       do l = itest+1, n
          m = nnp1 - l
          absc(l) = -absc(m)
          wght(l) =  wght(m)
       enddo
    endif

    ! COMPUTE COEFFICIENTS FOR AFFINE TRANSFORMATION
    ! AND ABSOLULE VALUE OF JACOBIAN

    u21 = u2 - u1
    v21 = v2 - v1
    u31 = u3 - u1
    v31 = v3 - v1
    absjac = abs(u21*v31 - v21*u31)

    ! INITIALIZE PARAMETERS

    u21 = u21 * 0.50_r8
    v21 = v21 * 0.50_r8
    u31 = u31 * 0.25_r8
    v31 = v31 * 0.25_r8

    ! CACLCULATE CUBATURE SUM

    do l = 1, n

       ! SAVE CONSTANTS TO AVOID REPETITIOUS CALCULATION

       save1 = 1.0_r8 + absc(l)
       save2 = 1.0_r8 - absc(l)
       save3 = u31 * save1
       save4 = v31 * save1
       save5 = save1 * wght(l)
       xx    = u21 * save2 + u1
       yy    = v21 * save2 + v1
       temp  = 0.0_r8

       do m = 1, n

          save6 = 1.0_r8 - absc(m)
          x     = xx + save3 * save6
          y     = yy + save4 * save6
          temp  = temp + wght(m) * f(x,y)

       enddo

       approx = approx + save5 * temp

    enddo

    approx = approx * absjac * 0.125_r8

  end subroutine dtria

!===============================================================================

  subroutine grule (n, x, w)

    use consts_coms, only: r8
    implicit none

    integer,  intent( in) :: n
    real(r8), intent(out) :: x(n)
    real(r8), intent(out) :: w(n)

    ! DETERMINES THE (N+1)/2 NONNEGATIVE POINTS X(I) AND
    ! THE CORRESPONDING WEIGHTS W(I) OF THE N-POINT
    ! GAUSS-LEGENDRE INTEGRATION RULE, NORMALIZED TO THE
    ! INTERVAL \-1,1\. THE X(I) APPEAR IN DESCENDING ORDER.
!
    ! THIS ROUTINE IS FROM 'METHODS OF NUMERICAL INTEGRATION',
    ! P.J. DAVIS AND P. RABINOWITZ, PAGE 369.

    integer  :: m, i, it, k
    real(r8) :: e1, t, x0, pk, pkm1, t1, pkp1, den, d1, n8
    real(r8) :: dpn, d2pn, d3pn, d4pn, u, v, h, p, dp, fx

    real(r8), parameter :: zero = 0.0_r8
    real(r8), parameter :: one  = 1.0_r8
    real(r8), parameter :: two  = 2.0_r8
    real(r8), parameter :: pi   = 4.0_r8 * atan(one)

    m  = (n+1)/2
    n8 = real(n,r8)
    e1 = n8*(n8+one)

    do i = 1, m
       t = dble(4*i-1)*pi/dble(4*n+2)
       x0 = (1.d0-(1.d0-1.d0/n8)/(8.d0*n8*n8))*cos(t)

       ! ITERATE ON THE VALUE  (M.W. JAN. 1982)
       do it = 1, 3
          pkm1=one
          pk=x0

          do k = 2, n
             t1=x0*pk
             pkp1=t1-pkm1-(t1-pkm1)/real(k,r8)+t1
             pkm1=pk
             pk=pkp1
          enddo

          den=one-x0*x0
          d1=n8*(pkm1-x0*pk)
          dpn=d1/den
          d2pn=(2.d0*x0*dpn-e1*pk)/den
          d3pn=(4.d0*x0*d2pn+(2.d0-e1)*dpn)/den
          d4pn=(6.d0*x0*d3pn+(6.d0-e1)*d2pn)/den
          u=pk/dpn
          v=d2pn/dpn
          h=-u*(1.d0+.5d0*u*(v+u*(v*v-u*d3pn/(3.d0*dpn))))
          p=pk+h*(dpn+.5d0*h*(d2pn+h/3.d0*(d3pn+.25d0*h*d4pn)))
          dp=dpn+h*(d2pn+.5d0*h*(d3pn+h*d4pn/3.d0))
          h=h-p/dp
          x0=x0+h
       enddo

       x(i)=x0
       fx=d1-h*e1*(pk+.5d0*h*(dpn+h/3.d0*(d2pn+.25d0*h*(d3pn+.2d0*h*d4pn))))
       w(i)=two*(one-x(i)*x(i))/(fx*fx)
    enddo

    if (m+m > n) x(m) = zero

  end subroutine grule


end module quadrature
