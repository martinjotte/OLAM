Module Lapack_Single

  use Lapack_utils

  ! This module contains a subset of the Lapack/Blas single precision library routines.
  ! Currently none are used in OLAM.

Contains

!==============================================================================
!
!  ILASLC scans a matrix for its last non-zero column.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download ILASLC + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslc.f
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILASLC( M, N, A, LDA )
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA,* )
!       ..
!
!  Purpose:
!  =============
!
! ILASLC scans A for its last non-zero column.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.
!
! param[in] A
!
!          A is REAL array, dimension (LDA,N)
!          The m by n matrix A.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! September 2012
!
!=====================================================================

  INTEGER FUNCTION ILASLC( M, N, A, LDA )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: M, N, LDA

    ! .. Array Arguments ..
    REAL,    INTENT(IN) :: A( LDA,* )

!=====================================================================

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER :: I

    ! .. Executable Statements ..

    ! Quick test for the common case where one corner is non-zero.
    IF ( N == 0 ) THEN
       ILASLC = N
    ELSE IF ( ABS( A(1,N) ) >= EPS .OR. ABS( A(M,N) ) >= EPS ) THEN
       ILASLC = N
    ELSE
       ! Now scan each column from the end, returning with the first non-zero.
       DO ILASLC = N, 1, -1
          DO I = 1, M
             IF ( ABS( A(I,ILASLC) ) >= EPS ) RETURN
          END DO
       END DO
    END IF

  END FUNCTION ILASLC

!=====================================================================
!
!  ILASLR scans a matrix for its last non-zero row.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download ILASLR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslr.f
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILASLR( M, N, A, LDA )
!
!       IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA,* )
!       ..
!
!  Purpose:
!  =============
!
! ILASLR scans A for its last non-zero row.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.
!
! param[in] A
!
!          A is REAL array, dimension (LDA,N)
!          The m by n matrix A.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! September 2012
!
!=====================================================================

  INTEGER FUNCTION ILASLR( M, N, A, LDA )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: M, N, LDA

    ! .. Array Arguments ..
    REAL,    INTENT(IN) :: A( LDA,* )

!=====================================================================

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER :: I, J

    ! .. Executable Statements ..

    ! Quick test for the common case where one corner is non-zero.
    IF ( M == 0 ) THEN
       ILASLR = M
    ELSE IF ( ABS( A(M,1) ) >= EPS .OR. ABS( A(M,N) ) >= EPS ) THEN
       ILASLR = M
    ELSE
       ! Scan up each column tracking the last zero row seen.
       ILASLR = 0
       DO J = 1, N
          I=M
          DO WHILE ( (ABS( A(MAX(I,1),J) ) < EPS) .AND. (I >= 1) )
             I=I-1
          ENDDO
          ILASLR = MAX( ILASLR, I )
       END DO
    END IF

  END FUNCTION ILASLR

!=====================================================================
!
!  SLAMCH
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      REAL             FUNCTION SLAMCH( CMACH )
!
!       IMPLICIT NONE (EXTERNAL, TYPE)
!
!     .. Scalar Arguments ..
!      CHARACTER          CMACH
!     ..
!
!  Purpose:
!  =============
!
! SLAMCH determines single precision machine parameters.
!
!  Arguments:
!  ==========
!
! param[in] CMACH
!
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!          where
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!=====================================================================

  REAL FUNCTION SLAMCH( CMACHS )

    USE IEEE_ARITHMETIC, ONLY: IEEE_QUIET_EQ
    IMPORT,              ONLY: LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: CMACHS

!=====================================================================

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0

    ! .. Local Scalars ..
    REAL      :: RND, EPS, SFMIN, SMALL, RMACH
    CHARACTER :: CMACH

    ! .. Executable Statements ..

    CMACH = CMACHS(1:1)

    ! Assume rounding, not chopping. Always.

    RND = ONE

    IF ( IEEE_QUIET_EQ( ONE, RND ) ) THEN
       EPS = EPSILON(ZERO) * 0.5
    ELSE
       EPS = EPSILON(ZERO)
    END IF

    IF ( LSAME( CMACH, 'E' ) ) THEN
       RMACH = EPS
    ELSE IF ( LSAME( CMACH, 'S' ) ) THEN
       SFMIN = TINY(ZERO)
       SMALL = ONE / HUGE(ZERO)
       IF ( SMALL >= SFMIN ) THEN
          ! Use SMALL plus a bit, to avoid the possibility of rounding
          ! causing overflow when computing  1/sfmin.
          SFMIN = SMALL*( ONE+EPS )
       END IF
       RMACH = SFMIN
    ELSE IF ( LSAME( CMACH, 'B' ) ) THEN
       RMACH = RADIX(ZERO)
    ELSE IF ( LSAME( CMACH, 'P' ) ) THEN
       RMACH = EPS * RADIX(ZERO)
    ELSE IF ( LSAME( CMACH, 'N' ) ) THEN
       RMACH = DIGITS(ZERO)
    ELSE IF ( LSAME( CMACH, 'R' ) ) THEN
       RMACH = RND
    ELSE IF ( LSAME( CMACH, 'M' ) ) THEN
       RMACH = MINEXPONENT(ZERO)
    ELSE IF ( LSAME( CMACH, 'U' ) ) THEN
       RMACH = tiny(zero)
    ELSE IF ( LSAME( CMACH, 'L' ) ) THEN
       RMACH = MAXEXPONENT(ZERO)
    ELSE IF ( LSAME( CMACH, 'O' ) ) THEN
       RMACH = HUGE(ZERO)
    ELSE
       RMACH = ZERO
    END IF

    SLAMCH = RMACH

  END FUNCTION SLAMCH

!=====================================================================
!
!  SCOPY
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*),SY(*)
!       ..
!
!  Purpose:
!  =============
!
!    SCOPY copies a vector, x, to a vector, y.
!    uses unrolled loops for increments equal to 1.
!
!  Authors:
!  ========
!
!  Univ. of Tennessee
!  Univ. of California Berkeley
!  Univ. of Colorado Denver
!  NAG Ltd.
!
!  November 2011
!
!  Further Details:
!  =====================
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!=====================================================================

  SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: INCX, INCY, N

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: SX(*), SY(*)

!=====================================================================

    ! .. Local Scalars ..
    INTEGER :: I,IX,IY,M,MP1

    IF (N <= 0) RETURN

    IF (INCX == 1 .AND. INCY == 1) THEN

       ! code for both increments equal to 1

       ! clean-up loop
       M = MOD(N,7)
       IF (M /= 0) THEN
          DO I = 1,M
             SY(I) = SX(I)
          END DO
          IF (N < 7) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,7
          SY(I) = SX(I)
          SY(I+1) = SX(I+1)
          SY(I+2) = SX(I+2)
          SY(I+3) = SX(I+3)
          SY(I+4) = SX(I+4)
          SY(I+5) = SX(I+5)
          SY(I+6) = SX(I+6)
       END DO

    ELSE

       ! code for unequal increments or equal increments not equal to 1

       IX = 1
       IY = 1
       IF (INCX < 0) IX = (-N+1)*INCX + 1
       IF (INCY < 0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          SY(IY) = SX(IX)
          IX = IX + INCX
          IY = IY + INCY
       END DO

    END IF

  END SUBROUTINE SCOPY

!=====================================================================
!
!  SGELQ2
!
!  Purpose
!  =======
!
!  SGELQ2 computes an LQ factorization of a real m by n matrix A:
!  A = L * Q.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and below the diagonal of the array
!          contain the m by min(m,n) lower trapezoidal matrix L (L is
!          lower triangular if m <= n); the elements above the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) REAL array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) REAL array, dimension (M)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!  and tau in TAU(i).
!
!=====================================================================

  SUBROUTINE SGELQ2( M, N, A, LDA, TAU, WORK, INFO )

    IMPORT, ONLY: XERBLA, SLARF, SLARFG

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ONE = 1.0

    ! .. Local Scalars ..
    INTEGER         :: I, K
    REAL            :: AII

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0
    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -4
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SGELQ2', -INFO )
       RETURN
    END IF

    K = MIN( M, N )

    DO I = 1, K

       ! Generate elementary reflector H(i) to annihilate A(i,i+1:n)

       CALL SLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, TAU( I ) )

       IF ( I < M ) THEN

          ! Apply H(i) to A(i+1:m,i:n) from the right
          AII = A( I, I )
          A( I, I ) = ONE
          CALL SLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), &
                      A( I+1, I ), LDA, WORK )
          A( I, I ) = AII

       END IF
    END DO

  END SUBROUTINE SGELQ2
!
!=====================================================================
!
!  SGELQF
!
!  Purpose
!  =======
!
!  SGELQF computes an LQ factorization of a real M-by-N matrix A:
!  A = L * Q.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and below the diagonal of the array
!          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
!          lower triangular if m <= n); the elements above the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) REAL array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!  and tau in TAU(i).
!
!=====================================================================

  SUBROUTINE SGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

    IMPORT, ONLY: XERBLA, ILAENV, SGELQ2, SLARFB, SLARFT

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Local Scalars ..
    LOGICAL :: LQUERY
    INTEGER :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0
    NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
    LWKOPT = M*NB
    WORK( 1 ) = REAL( LWKOPT )
    LQUERY = ( LWORK == -1 )

    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -4
    ELSE IF ( LWORK < MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
       INFO = -7
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SGELQF', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    K = MIN( M, N )
    IF ( K == 0 ) THEN
       WORK( 1 ) = 1
       RETURN
    END IF

    NBMIN = 2
    NX = 0
    IWS = M

    IF ( NB > 1 .AND. NB < K ) THEN

       ! Determine when to cross over from blocked to unblocked code.

       NX = MAX( 0, ILAENV( 3, 'SGELQF', ' ', M, N, -1, -1 ) )

       IF ( NX < K ) THEN

          ! Determine if workspace is large enough for blocked code.

          LDWORK = M
          IWS = LDWORK*NB

          IF ( LWORK < IWS ) THEN
             ! Not enough workspace to use optimal NB:  reduce NB and
             ! determine the minimum value of NB.
             NB = LWORK / LDWORK
             NBMIN = MAX( 2, ILAENV( 2, 'SGELQF', ' ', M, N, -1, -1 ) )
          END IF

       END IF

    END IF

    IF ( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN

       ! Use blocked code initially

       DO I = 1, K - NX, NB
          IB = MIN( K-I+1, NB )

          ! Compute the LQ factorization of the current block
          ! A(i:i+ib-1,i:n)

          CALL SGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

          IF ( I+IB <= M ) THEN

             ! Form the triangular factor of the block reflector
             ! H = H(i) H(i+1) . . . H(i+ib-1)

             CALL SLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), LDA, &
                          TAU( I ), WORK, LDWORK )

             ! Apply H to A(i+ib:m,i:n) from the right

             CALL SLARFB( 'Right', 'No transpose', 'Forward', 'Rowwise', &
                          M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, &
                          A( I+IB, I ), LDA, WORK( IB+1 ), LDWORK )
          END IF
       END DO

    ELSE

       I = 1

    END IF

    ! Use unblocked code to factor the last or only block.

    IF ( I <= K ) &
         CALL SGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

    WORK( 1 ) = REAL( IWS )

  END SUBROUTINE SGELQF

!=====================================================================
!
!  SGELS
!
!  Purpose
!  =======
!
!  SGELS solves overdetermined or underdetermined real linear systems
!  involving an M-by-N matrix A, or its transpose, using a QR or LQ
!  factorization of A.  It is assumed that A has full rank.
!
!  The following options are provided:
!
!  1. If TRANS = 'N' and m >= n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A*X ||.
!
!  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!     an underdetermined system A * X = B.
!
!  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
!     an undetermined system A**T * X = B.
!
!  4. If TRANS = 'T' and m < n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A**T * X ||.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!  matrix X.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          = 'N': the linear system involves A;
!          = 'T': the linear system involves A**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of
!          columns of the matrices B and X. NRHS >=0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!            if M >= N, A is overwritten by details of its QR
!                       factorization as returned by SGEQRF;
!            if M <  N, A is overwritten by details of its LQ
!                       factorization as returned by SGELQF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the matrix B of right hand side vectors, stored
!          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!          if TRANS = 'T'.
!          On exit, if INFO = 0, B is overwritten by the solution
!          vectors, stored columnwise:
!          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!          squares solution vectors; the residual sum of squares for the
!          solution in each column is given by the sum of squares of
!          elements N+1 to M in that column;
!          if TRANS = 'N' and m < n, rows 1 to N of B contain the
!          minimum norm solution vectors;
!          if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!          minimum norm solution vectors;
!          if TRANS = 'T' and m < n, rows 1 to M of B contain the
!          least squares solution vectors; the residual sum of squares
!          for the solution in each column is given by the sum of
!          squares of elements M+1 to N in that column.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          LWORK >= max( 1, MN + max( MN, NRHS ) ).
!          For optimal performance,
!          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!          where MN = min(M,N) and NB is the optimum block size.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO =  i, the i-th diagonal element of the
!                triangular factor of A is zero, so that A does not have
!                full rank; the least squares solution could not be
!                computed.
!
!=====================================================================

  SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

    IMPORT, ONLY: SGELQF, SGEQRF, SLABAD, SLASCL, SLASET, SORMLQ, &
                  SORMQR, STRTRS, ILAENV, SLAMCH, SLANGE, XERBLA, LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK driver routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER          TRANS
    INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS

    ! .. Array Arguments ..
    REAL               A( LDA,* ), B( LDB,* ), WORK( * )

    ! .. Parameters ..
    REAL               ZERO, ONE
    PARAMETER          ( ZERO = 0.0, ONE = 1.0 )
    REAL, PARAMETER :: EPS = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL            LQUERY, TPSD
    INTEGER            BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
    REAL               ANRM, BIGNUM, BNRM, SMLNUM

    ! .. Local Arrays ..
    REAL               RWORK( 1 )

    ! .. Executable Statements ..

    ! Test the input arguments.

    INFO = 0
    MN = MIN( M, N )
    LQUERY = ( LWORK == -1 )

    IF ( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'T' ) ) ) THEN
       INFO = -1
    ELSE IF ( M < 0 ) THEN
       INFO = -2
    ELSE IF ( N < 0 ) THEN
       INFO = -3
    ELSE IF ( NRHS < 0 ) THEN
       INFO = -4
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -6
    ELSE IF ( LDB < MAX( 1, M, N ) ) THEN
       INFO = -8
    ELSE IF ( LWORK < MAX( 1, MN + MAX( MN, NRHS ) ) .AND. .NOT.LQUERY ) THEN
       INFO = -10
    END IF

    ! Figure out optimal block size

    IF ( INFO == 0 .OR. INFO == -10 ) THEN

       TPSD = .TRUE.
       IF ( LSAME( TRANS, 'N' ) ) TPSD = .FALSE.

       IF ( M >= N ) THEN
          NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
          IF ( TPSD ) THEN
             NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LN', M, NRHS, N, -1 ) )
          ELSE
             NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LT', M, NRHS, N, -1 ) )
          END IF
       ELSE
          NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
          IF ( TPSD ) THEN
             NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LT', N, NRHS, M, -1 ) )
          ELSE
             NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LN', N, NRHS, M, -1 ) )
          END IF
       END IF

       WSIZE = MAX( 1, MN + MAX( MN, NRHS )*NB )
       WORK( 1 ) = REAL( WSIZE )

    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SGELS ', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( MIN( M, N, NRHS ) == 0 ) THEN
       CALL SLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
       RETURN
    END IF

    ! Get machine parameters

    SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
    BIGNUM = ONE / SMLNUM
    CALL SLABAD( SMLNUM, BIGNUM )

    ! Scale A, B if max element outside range [SMLNUM,BIGNUM]

    ANRM = SLANGE( 'M', M, N, A, LDA, RWORK )
    IASCL = 0

    IF ( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN

       ! Scale matrix norm up to SMLNUM
       CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
       IASCL = 1

    ELSE IF ( ANRM > BIGNUM ) THEN

       ! Scale matrix norm down to BIGNUM
       CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
       IASCL = 2

    ELSE IF ( ANRM < EPS ) THEN

       ! Matrix all zero. Return zero solution.
       CALL SLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
       WORK( 1 ) = REAL( WSIZE )
       RETURN

    END IF

    BROW = M
    IF ( TPSD ) BROW = N
    BNRM = SLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
    IBSCL = 0

    IF ( BNRM > ZERO .AND. BNRM < SMLNUM ) THEN

       ! Scale matrix norm up to SMLNUM
       CALL SLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO )
       IBSCL = 1

    ELSE IF ( BNRM > BIGNUM ) THEN

       ! Scale matrix norm down to BIGNUM
       CALL SLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO )
       IBSCL = 2

    END IF

    IF ( M >= N ) THEN

       ! compute QR factorization of A

       CALL SGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

       ! workspace at least N, optimally N*NB

       IF ( .NOT.TPSD ) THEN

          ! Least-Squares Problem min || A * X - B ||

          ! B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

          CALL SORMQR( 'Left', 'Transpose', M, NRHS, N, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          ! B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

          CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS, &
                       A, LDA, B, LDB, INFO )

          IF ( INFO > 0 ) THEN
             RETURN
          END IF

          SCLLEN = N

       ELSE

          ! Overdetermined system of equations A**T * X = B

          ! B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

          CALL STRTRS( 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                       A, LDA, B, LDB, INFO )

          IF ( INFO > 0 ) THEN
             RETURN
          END IF

          ! B(N+1:M,1:NRHS) = ZERO

          DO J = 1, NRHS
             DO I = N + 1, M
                B( I, J ) = ZERO
             END DO
          END DO

          ! B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

          CALL SORMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          SCLLEN = M

       END IF

    ELSE

       ! Compute LQ factorization of A

       CALL SGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

       ! workspace at least M, optimally M*NB.

       IF ( .NOT.TPSD ) THEN

          ! underdetermined system of equations A * X = B

          ! B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

          CALL STRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS, &
                       A, LDA, B, LDB, INFO )

          IF ( INFO > 0 ) THEN
             RETURN
          END IF

          ! B(M+1:N,1:NRHS) = 0

          DO J = 1, NRHS
             DO I = M + 1, N
                B( I, J ) = ZERO
             END DO
          END DO

          ! B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

          CALL SORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          SCLLEN = N

       ELSE

          ! overdetermined system min || A**T * X - B ||

          ! B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

          CALL SORMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          ! B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

          CALL STRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS, &
                       A, LDA, B, LDB, INFO )

          IF ( INFO > 0 ) THEN
             RETURN
          END IF

          SCLLEN = M

       END IF

    END IF

    ! Undo scaling

    IF ( IASCL == 1 ) THEN
       CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO )
    ELSE IF ( IASCL == 2 ) THEN
       CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO )
    END IF

    IF ( IBSCL == 1 ) THEN
       CALL SLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
    ELSE IF ( IBSCL == 2 ) THEN
       CALL SLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
    END IF

    WORK( 1 ) = REAL( WSIZE )

  END SUBROUTINE SGELS

!=====================================================================
!
!  SGEQR2
!
!  Purpose
!  =======
!
!  SGEQR2 computes a QR factorization of a real m by n matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) REAL array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!=====================================================================

  SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )

    IMPORT, ONLY: XERBLA, SLARF, SLARFG

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ONE = 1.0

    ! .. Local Scalars ..
    INTEGER :: I, K
    REAL    :: AII

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0
    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -4
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SGEQR2', -INFO )
       RETURN
    END IF

    K = MIN( M, N )

    DO I = 1, K

       ! Generate elementary reflector H(i) to annihilate A(i+1:m,i)

       CALL SLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, TAU( I ) )

       IF ( I < N ) THEN

          ! Apply H(i) to A(i:m,i+1:n) from the left
          AII = A( I, I )
          A( I, I ) = ONE
          CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                     A( I, I+1 ), LDA, WORK )
          A( I, I ) = AII

       END IF

    END DO

  END SUBROUTINE SGEQR2

!=====================================================================
!
!  SGEQRF
!
!  Purpose
!  =======
!
!  SGEQRF computes a QR factorization of a real M-by-N matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of min(m,n) elementary reflectors (see Further
!          Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) REAL array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v**T
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!=====================================================================

  SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

    IMPORT,  ONLY: XERBLA, ILAENV, SGEQR2, SLARFB, SLARFT

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Local Scalars ..
    LOGICAL :: LQUERY
    INTEGER :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0
    NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
    LWKOPT = N*NB
    WORK( 1 ) = REAL( LWKOPT )
    LQUERY = ( LWORK == -1 )

    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -4
    ELSE IF ( LWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
       INFO = -7
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SGEQRF', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    K = MIN( M, N )
    IF ( K == 0 ) THEN
       WORK( 1 ) = 1
       RETURN
    END IF

    NBMIN = 2
    NX = 0
    IWS = N

    IF ( NB > 1 .AND. NB < K ) THEN

       ! Determine when to cross over from blocked to unblocked code.

       NX = MAX( 0, ILAENV( 3, 'SGEQRF', ' ', M, N, -1, -1 ) )

       IF ( NX < K ) THEN

          ! Determine if workspace is large enough for blocked code.

          LDWORK = N
          IWS = LDWORK*NB

          IF ( LWORK < IWS ) THEN
             ! Not enough workspace to use optimal NB:  reduce NB and
             ! determine the minimum value of NB.
             NB = LWORK / LDWORK
             NBMIN = MAX( 2, ILAENV( 2, 'SGEQRF', ' ', M, N, -1, -1 ) )
          END IF

       END IF

    END IF

    IF ( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN

       ! Use blocked code initially

       DO I = 1, K - NX, NB
          IB = MIN( K-I+1, NB )

          ! Compute the QR factorization of the current block
          ! A(i:m,i:i+ib-1)

          CALL SGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )

          IF ( I+IB <= N ) THEN

             ! Form the triangular factor of the block reflector
             ! H = H(i) H(i+1) . . . H(i+ib-1)

             CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                          A( I, I ), LDA, TAU( I ), WORK, LDWORK )

             ! Apply H**T to A(i:m,i+ib:n) from the left

             CALL SLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', &
                          M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, &
                          A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK )
          END IF
       END DO

    ELSE

       I = 1

    END IF

    ! Use unblocked code to factor the last or only block.

    IF ( I <= K ) &
         CALL SGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

    WORK( 1 ) = REAL( IWS )

  END SUBROUTINE SGEQRF

!=====================================================================
!
!  SISNAN
!
!  Purpose
!  =======
!
!  SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!  future.
!
!  Arguments
!  =========
!
!  SIN     (input) REAL
!          Input to test for NaN.
!
!=====================================================================

  LOGICAL FUNCTION SISNAN( SSIN )

    USE IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! June 2010

    ! .. Scalar Arguments ..
    REAL, INTENT(IN) :: SSIN

    ! .. Executable Statements ..
    SISNAN = IEEE_IS_NAN( SSIN )

  END FUNCTION SISNAN

!=====================================================================
!
!  SLABAD
!
!  Purpose
!  =======
!
!  SLABAD takes as input the values computed by SLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by SLAMCH.  This subroutine is needed because
!  SLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) REAL
!          On entry, the underflow threshold as computed by SLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) REAL
!          On entry, the overflow threshold as computed by SLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!=====================================================================

  SUBROUTINE SLABAD( SMALL, LARGE )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    REAL :: LARGE, SMALL


    ! .. Executable Statements ..

    ! If it looks like we're on a Cray, take the square root of
    ! SMALL and LARGE to avoid overflow and underflow problems.

    IF ( LOG10( LARGE ) > 2000. ) THEN
       SMALL = SQRT( SMALL )
       LARGE = SQRT( LARGE )
    END IF

  END SUBROUTINE SLABAD

!=====================================================================
!
!  SLANGE
!
!  Purpose
!  =======
!
!  SLANGE  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real matrix A.
!
!  Description
!  ===========
!
!  SLANGE returns the value
!
!     SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in SLANGE as described
!          above.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          SLANGE is set to zero.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          SLANGE is set to zero.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
!  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
!=====================================================================

  REAL FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )

    IMPORT, ONLY: LSAME, XERBLA, SLASSQ

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN) :: NORM
    INTEGER,   INTENT(IN) :: LDA, M, N

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), WORK( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0

    ! .. Local Scalars ..
    INTEGER  :: I, J
    REAL :: SCALE, SUM, VALUE

    ! .. Executable Statements ..

    IF ( MIN( M, N ) == 0 ) THEN

       VALUE = ZERO

    ELSE IF ( LSAME( NORM, 'M' ) ) THEN

       ! Find max(abs(A(i,j))).
       VALUE = ZERO
       DO J = 1, N
          DO I = 1, M
             VALUE = MAX( VALUE, ABS( A( I, J ) ) )
          END DO
       END DO

    ELSE IF ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) THEN

       ! Find norm1(A).
       VALUE = ZERO
       DO J = 1, N
          SUM = ZERO
          DO I = 1, M
             SUM = SUM + ABS( A( I, J ) )
          END DO
          VALUE = MAX( VALUE, SUM )
       END DO

    ELSE IF ( LSAME( NORM, 'I' ) ) THEN

       ! Find normI(A).
       DO I = 1, M
          WORK( I ) = ZERO
       END DO
       DO J = 1, N
          DO I = 1, M
             WORK( I ) = WORK( I ) + ABS( A( I, J ) )
          END DO
       END DO
       VALUE = ZERO
       DO I = 1, M
          VALUE = MAX( VALUE, WORK( I ) )
       END DO

    ELSE IF ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN

       ! Find normF(A).
       SCALE = ZERO
       SUM = ONE
       DO J = 1, N
          CALL SLASSQ( M, A( 1, J ), 1, SCALE, SUM )
       END DO
       VALUE = SCALE*SQRT( SUM )

    END IF

    SLANGE = VALUE

  END FUNCTION SLANGE

!=====================================================================
!
!  SLAPY2
!
!  Purpose
!  =======
!
!  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) REAL
!  Y       (input) REAL
!          X and Y specify the values x and y.
!
!=====================================================================

  REAL FUNCTION SLAPY2( X, Y )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    REAL, INTENT(IN) :: X, Y

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL :: W, XABS, YABS, Z

    ! .. Executable Statements ..

    XABS = ABS( X )
    YABS = ABS( Y )
    W = MAX( XABS, YABS )
    Z = MIN( XABS, YABS )
    IF ( Z < EPS ) THEN
       SLAPY2 = W
    ELSE
       SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
    END IF

  END FUNCTION SLAPY2

!=====================================================================
!
!  SLARF
!
!  Purpose
!  =======
!
!  SLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v**T
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) REAL array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) REAL
!          The value tau in the representation of H.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!=====================================================================

  SUBROUTINE SLARF( SIDEC, M, N, V, INCV, TAU, C, LDC, WORK )

    IMPORT, ONLY: LSAME, XERBLA, ILASLR, ILASLC, SGEMV, SGER

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: SIDEC
    INTEGER,      INTENT(IN) :: INCV, LDC, M, N
    REAL,         INTENT(IN) :: TAU

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: C( LDC,* ), V( * ), WORK( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL   :: APPLYLEFT
    INTEGER   :: I, LASTV, LASTC
    CHARACTER :: SIDE

    ! .. Executable Statements ..

    SIDE = SIDEC(1:1)

    APPLYLEFT = LSAME( SIDE, 'L' )
    LASTV = 0
    LASTC = 0

    IF ( ABS(TAU) >= EPS ) THEN
       ! Set up variables for scanning V.  LASTV begins pointing to the end
       ! of V.
       IF ( APPLYLEFT ) THEN
          LASTV = M
       ELSE
          LASTV = N
       END IF
       IF ( INCV > 0 ) THEN
          I = 1 + (LASTV-1) * INCV
       ELSE
          I = 1
       END IF
       ! Look for the last non-zero row in V.
       DO WHILE( LASTV > 0 .AND. ABS( V( I ) ) < EPS )
          LASTV = LASTV - 1
          I = I - INCV
       END DO
       IF ( APPLYLEFT ) THEN
          ! Scan for the last non-zero column in C(1:lastv,:).
          LASTC = ILASLC(LASTV, N, C, LDC)
       ELSE
          ! Scan for the last non-zero row in C(:,1:lastv).
          LASTC = ILASLR(M, LASTV, C, LDC)
       END IF

    END IF

    ! Note that lastc == 0 renders the BLAS operations null; no special
    ! case is needed at this level.

    IF ( APPLYLEFT ) THEN

       ! Form  H * C
       IF ( LASTV > 0 ) THEN

          ! w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
          CALL SGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
                      ZERO, WORK, 1 )

          ! C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
          CALL SGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )

       END IF

    ELSE

       ! Form  C * H
       IF ( LASTV > 0 ) THEN

          ! w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
          CALL SGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, &
                      ZERO, WORK, 1 )

          ! C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
          CALL SGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
       END IF

    END IF

  END SUBROUTINE SLARF

!=====================================================================
!
!  SLARFB
!
!  Purpose
!  =======
!
!  SLARFB applies a real block reflector H or its transpose H**T to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H**T from the Left
!          = 'R': apply H or H**T from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H**T (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) REAL array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See Further Details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) REAL array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!=====================================================================

  SUBROUTINE SLARFB( SIDEC, TRANSC, DIRECTC, STOREVC, M, N, K, V, LDV, &
                     T, LDT, C, LDC, WORK, LDWORK )

    IMPORT, ONLY: LSAME, SGEMM, ILASLR, ILASLC, SCOPY, STRMM

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: DIRECTC, SIDEC, STOREVC, TRANSC
    INTEGER,      INTENT(IN) :: K, LDC, LDT, LDV, LDWORK, M, N

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: C( LDC,* ), T( LDT,* ), V( LDV,* ), WORK( LDWORK,* )

    ! .. Parameters ..
    REAL, PARAMETER :: ONE = 1.0

    ! .. Local Scalars ..
    CHARACTER :: TRANST
    INTEGER   :: I, J
    CHARACTER :: DIRECT, SIDE, STOREV, TRANS

    ! .. Executable Statements ..

    DIRECT = DIRECTC(1:1)
    SIDE   = SIDEC  (1:1)
    STOREV = STOREVC(1:1)
    TRANS  = TRANSC (1:1)

    ! Quick return if possible

    IF ( M <= 0 .OR. N <= 0 ) RETURN

    IF ( LSAME( TRANS, 'N' ) ) THEN
       TRANST = 'T'
    ELSE
       TRANST = 'N'
    END IF

    IF ( LSAME( STOREV, 'C' ) ) THEN

       IF ( LSAME( DIRECT, 'F' ) ) THEN

          ! Let  V =  ( V1 )    (first K rows)
          !           ( V2 )
          ! where  V1  is unit lower triangular.

          IF ( LSAME( SIDE, 'L' ) ) THEN

             ! Form  H * C  or  H**T * C  where  C = ( C1 )
             !                                       ( C2 )
             ! W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
             ! W := C1**T

             DO J = 1, K
                CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1

             CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, K, &
                         ONE, V, LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C2**T *V2
                CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V * W**T

             IF ( M > K ) THEN
                ! C2 := C2 - V2 * W**T
                CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                            -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                            C( K+1, 1 ), LDC )
             END IF

             ! W := W * V1**T

             CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                         ONE, V, LDV, WORK, LDWORK )

             ! C1 := C1 - W**T

             DO J = 1, K
                DO I = 1, N
                   C( J, I ) = C( J, I ) - WORK( I, J )
                END DO
             END DO

          ELSE IF ( LSAME( SIDE, 'R' ) ) THEN

             ! Form  C * H  or  C * H**T  where  C = ( C1  C2 )
             ! W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
             ! W := C1

             DO J = 1, K
                CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1

             CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, K, &
                         ONE, V, LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C2 * V2
                CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V**T

             IF ( N > K ) THEN
                ! C2 := C2 - W * V2**T
                CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                            C( 1, K+1 ), LDC )
             END IF

             ! W := W * V1**T

             CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                         ONE, V, LDV, WORK, LDWORK )

             ! C1 := C1 - W

             DO J = 1, K
                DO I = 1, M
                   C( I, J ) = C( I, J ) - WORK( I, J )
                END DO
             END DO
          END IF

       ELSE

          ! Let  V =  ( V1 )
          !           ( V2 )    (last K rows)
          ! where  V2  is unit upper triangular.

          IF ( LSAME( SIDE, 'L' ) ) THEN

             ! Form  H * C  or  H**T * C  where  C = ( C1 )
             !                                       ( C2 )
             ! W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
             ! W := C2**T

             DO J = 1, K
                CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2

             CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, K, &
                         ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C1**T*V1
                CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V * W**T

             IF ( M > K ) THEN
                ! C1 := C1 - V1 * W**T
                CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                            -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
             END IF

             ! W := W * V2**T

             CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C2 := C2 - W**T

             DO J = 1, K
                DO I = 1, N
                   C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
                END DO
             END DO

          ELSE IF ( LSAME( SIDE, 'R' ) ) THEN

             ! Form  C * H  or  C * H**T  where  C = ( C1  C2 )
             ! W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
             ! W := C2

             DO J = 1, K
                CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2

             CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, K, &
                         ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C1 * V1
                CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V**T

             IF ( N > K ) THEN
                ! C1 := C1 - W * V1**T
                CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
             END IF

             ! W := W * V2**T

             CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                         ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )

             ! C2 := C2 - W

             DO J = 1, K
                DO I = 1, N
                   C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
                END DO
             END DO
          END IF
       END IF

    ELSE IF ( LSAME( STOREV, 'R' ) ) THEN

       IF ( LSAME( DIRECT, 'F' ) ) THEN

          ! Let  V =  ( V1  V2 )    (V1: first K columns)
          ! where  V1  is unit upper triangular.

          IF ( LSAME( SIDE, 'L' ) ) THEN

             ! Form  H * C  or  H**T * C  where  C = ( C1 )
             !                                       ( C2 )
             ! W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
             ! W := C1**T

             DO J = 1, K
                CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1**T

             CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                         ONE, V, LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C2**T*V2**T
                CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                     C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V**T * W**T

             IF ( M > K ) THEN
                ! C2 := C2 - V2**T * W**T
                CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                            V( 1, K+1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC )
             END IF

             ! W := W * V1

             CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, K, &
                         ONE, V, LDV, WORK, LDWORK )

             ! C1 := C1 - W**T

             DO J = 1, K
                DO I = 1, N
                   C( J, I ) = C( J, I ) - WORK( I, J )
                END DO
             END DO

          ELSE IF ( LSAME( SIDE, 'R' ) ) THEN

             ! Form  C * H  or  C * H**T  where  C = ( C1  C2 )
             ! W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
             ! W := C1

             DO J = 1, K
                CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1**T

             CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                         ONE, V, LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C2 * V2**T
                CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                            ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V

             IF ( N > K ) THEN
                ! C2 := C2 - W * V2
                CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                            C( 1, K+1 ), LDC )
             END IF

             ! W := W * V1

             CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, K, &
                         ONE, V, LDV, WORK, LDWORK )

             ! C1 := C1 - W

             DO J = 1, K
                DO I = 1, M
                   C( I, J ) = C( I, J ) - WORK( I, J )
                END DO
             END DO

          END IF

       ELSE

          ! Let  V =  ( V1  V2 )    (V2: last K columns)
          ! where  V2  is unit lower triangular.

          IF ( LSAME( SIDE, 'L' ) ) THEN

             ! Form  H * C  or  H**T * C  where  C = ( C1 )
             !                                       ( C2 )
             ! W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
             ! W := C2**T

             DO J = 1, K
                CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2**T

             CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                         ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C1**T * V1**T
                CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                            C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V**T * W**T

             IF ( M > K ) THEN
                ! C1 := C1 - V1**T * W**T
                CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                            V, LDV, WORK, LDWORK, ONE, C, LDC )
             END IF

             ! W := W * V2

             CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, K, &
                         ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )

             ! C2 := C2 - W**T

             DO J = 1, K
                DO I = 1, N
                   C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
                END DO
             END DO

          ELSE IF ( LSAME( SIDE, 'R' ) ) THEN

             ! Form  C * H  or  C * H**T  where  C = ( C1  C2 )
             ! W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
             ! W := C2

             DO J = 1, K
                CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2**T

             CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                         ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C1 * V1**T
                CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
              END IF

             ! W := W * T  or  W * T**T

             CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V

             IF (N > K ) THEN
                ! C1 := C1 - W * V1
                CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
             END IF

             ! W := W * V2

             CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, K, &
                         ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )

             ! C1 := C1 - W

             DO J = 1, K
                DO I = 1, M
                   C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
                END DO
             END DO

          END IF

       END IF
    END IF

  END SUBROUTINE SLARFB

!=====================================================================
!
!  SLARFG
!
!  Purpose
!  =======
!
!  SLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H**T * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v**T ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) REAL
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) REAL array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) REAL
!          The value tau.
!
!=====================================================================

  SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )

    IMPORT, ONLY: SSCAL, SLAMCH, SLAPY2, SNRM2

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)    :: INCX, N
    REAL,    INTENT(INOUT) :: ALPHA, TAU

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: X( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER :: J, KNT
    REAL    :: BETA, RSAFMN, SAFMIN, XNORM

    ! .. Executable Statements ..

    IF ( N <= 1 ) THEN
       TAU = ZERO
       RETURN
    END IF

    XNORM = SNRM2( N-1, X, INCX )

    IF ( XNORM < EPS ) THEN

       TAU = ZERO

    ELSE

       ! general case

       BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
       SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
       KNT = 0

       IF ( ABS( BETA ) < SAFMIN ) THEN

          ! XNORM, BETA may be inaccurate; scale X and recompute them

          RSAFMN = ONE / SAFMIN

          DO WHILE ( ABS( BETA ) < SAFMIN )
             KNT = KNT + 1
             CALL SSCAL( N-1, RSAFMN, X, INCX )
             BETA = BETA*RSAFMN
             ALPHA = ALPHA*RSAFMN
          END DO

          ! New BETA is at most 1, at least SAFMIN

          XNORM = SNRM2( N-1, X, INCX )
          BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )

       END IF

       TAU = ( BETA-ALPHA ) / BETA
       CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )

       ! If ALPHA is subnormal, it may lose relative accuracy

       DO J = 1, KNT
          BETA = BETA*SAFMIN
       END DO
       ALPHA = BETA

    END IF

  END SUBROUTINE SLARFG

!=====================================================================
!
!  SLARFT
!
!  Purpose
!  =======
!
!  SLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V**T
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V**T * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) REAL array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) REAL array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )

!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!=====================================================================

  SUBROUTINE SLARFT( DIRECTC, STOREVC, N, K, V, LDV, TAU, T, LDT )

    IMPORT,  ONLY: LSAME, SGEMV, STRMV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: DIRECTC, STOREVC
    INTEGER,      INTENT(IN) :: K, LDT, LDV, N

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: T( LDT,* ), TAU( * ), V( LDV,* )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER   :: I, J, PREVLASTV, LASTV
    CHARACTER :: DIRECT, STOREV

    ! .. Executable Statements ..

    DIRECT = DIRECTC(1:1)
    STOREV = STOREVC(1:1)

    ! Quick return if possible

    IF ( N == 0 ) RETURN

    IF ( LSAME( DIRECT, 'F' ) ) THEN

       PREVLASTV = N
       DO I = 1, K

          PREVLASTV = MAX( I, PREVLASTV )

          IF ( ABS( TAU(I) ) < EPS ) THEN

             ! H(i)  =  I
             DO J = 1, I
                T( J, I ) = ZERO
             END DO

          ELSE

             ! general case

             IF ( LSAME( STOREV, 'C' ) ) THEN

                ! Skip any trailing zeros.
                DO LASTV = N, I+1, -1
                   IF ( ABS( V(LASTV,I) ) >= EPS ) EXIT
                END DO
                DO J = 1, I-1
                   T( J, I ) = -TAU( I ) * V( I , J )
                END DO
                J = MIN( LASTV, PREVLASTV )

                ! T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)

                CALL SGEMV( 'Transpose', J-I, I-1, -TAU( I ), V( I+1,1 ), LDV, &
                            V( I+1,I ), 1, ONE, T( 1,I ), 1 )


             ELSE

                ! Skip any trailing zeros.
                DO LASTV = N, I+1, -1
                   IF ( ABS( V(I,LASTV) ) >= EPS ) EXIT
                END DO
                DO J = 1, I-1
                   T( J, I ) = -TAU( I ) * V( J , I )
                END DO
                J = MIN( LASTV, PREVLASTV )

                ! T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T

                CALL SGEMV( 'No transpose', I-1, J-I, -TAU( I ), V( 1,I+1 ), &
                            LDV, V( I,I+1 ), LDV, ONE, T( 1,I ), 1 )

             END IF

             ! T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

             CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1,I ), 1 )

             T( I, I ) = TAU( I )
             IF ( I > 1 ) THEN
                PREVLASTV = MAX( PREVLASTV, LASTV )
             ELSE
                PREVLASTV = LASTV
             END IF

          END IF

       END DO

    ELSE

       PREVLASTV = 1
       DO I = K, 1, -1

          IF ( ABS( TAU(I) ) < EPS ) THEN

             ! H(i)  =  I
             DO J = I, K
                T( J, I ) = ZERO
             END DO

          ELSE

             ! general case

             IF ( I < K ) THEN
                IF ( LSAME( STOREV, 'C' ) ) THEN

                   ! Skip any leading zeros.
                   DO LASTV = 1, I-1
                      IF ( ABS( V(LASTV,I) ) >= EPS ) EXIT
                   END DO
                   DO J = I+1, K
                      T( J, I ) = -TAU( I ) * V( N-K+I , J )
                   END DO
                   J = MAX( LASTV, PREVLASTV )

                   ! T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)

                   CALL SGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ), V( J,I+1 ), &
                               LDV, V( J,I ), 1, ONE, T( I+1,I ), 1 )

                ELSE

                   ! Skip any leading zeros.
                   DO LASTV = 1, I-1
                      IF ( ABS( V(I,LASTV) ) >= EPS ) EXIT
                   END DO
                   DO J = I+1, K
                      T( J, I ) = -TAU( I ) * V( J, N-K+I )
                   END DO
                   J = MAX( LASTV, PREVLASTV )

                   ! T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T

                   CALL SGEMV( 'No transpose', K-I, N-K+I-J, -TAU( I ), V( I+1,J ), &
                               LDV, V( I,J ), LDV, ONE, T( I+1,I ), 1 )

                END IF

                ! T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

                CALL STRMV( 'Lower', 'No transpose', 'Non-unit', K-I, T( I+1,I+1 ), &
                            LDT, T( I+1,I ), 1 )

                IF ( I > 1 ) THEN
                   PREVLASTV = MIN( PREVLASTV, LASTV )
                ELSE
                   PREVLASTV = LASTV
                END IF

             END IF

             T( I, I ) = TAU( I )

          END IF
       END DO

    END IF

  END SUBROUTINE SLARFT

!=====================================================================
!
!  SLASCL
!
!  Purpose
!  =======
!
!  SLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU. See SGBTRF for storage details.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) REAL
!  CTO     (input) REAL
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!=====================================================================

  SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

    USE IEEE_ARITHMETIC, ONLY: IEEE_QUIET_EQ
    IMPORT,              ONLY: LSAME, XERBLA, SISNAN, SLAMCH

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.3.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2010

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: TYPE
    INTEGER,   INTENT(IN)  :: KL, KU, LDA, M, N
    INTEGER,   INTENT(OUT) :: INFO
    REAL,      INTENT(IN)  :: CFROM, CTO

    ! .. Array Arguments ..
    REAL,  INTENT(INOUT) :: A( LDA,* )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL :: DONE
    INTEGER :: I, ITYPE, J, K1, K2, K3, K4
    REAL    :: BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0

    IF ( LSAME( TYPE, 'G' ) ) THEN
       ITYPE = 0
    ELSE IF ( LSAME( TYPE, 'L' ) ) THEN
       ITYPE = 1
    ELSE IF ( LSAME( TYPE, 'U' ) ) THEN
       ITYPE = 2
    ELSE IF ( LSAME( TYPE, 'H' ) ) THEN
       ITYPE = 3
    ELSE IF ( LSAME( TYPE, 'B' ) ) THEN
       ITYPE = 4
    ELSE IF ( LSAME( TYPE, 'Q' ) ) THEN
       ITYPE = 5
    ELSE IF ( LSAME( TYPE, 'Z' ) ) THEN
       ITYPE = 6
    ELSE
       ITYPE = -1
    END IF

    IF ( ITYPE == -1 ) THEN
       INFO = -1
    ELSE IF ( ABS(CFROM) < EPS .OR. SISNAN(CFROM) ) THEN
       INFO = -4
    ELSE IF ( SISNAN(CTO) ) THEN
       INFO = -5
    ELSE IF ( M < 0 ) THEN
       INFO = -6
    ELSE IF ( N < 0 .OR. ( ITYPE == 4 .AND. N /= M ) .OR. &
                         ( ITYPE == 5 .AND. N /= M ) ) THEN
       INFO = -7
    ELSE IF ( ITYPE <= 3 .AND. LDA < MAX( 1, M ) ) THEN
       INFO = -9
    ELSE IF ( ITYPE >= 4 ) THEN
       IF ( KL < 0 .OR. KL > MAX( M-1, 0 ) ) THEN
          INFO = -2
       ELSE IF ( KU < 0 .OR. KU > MAX( N-1, 0 ) .OR. &
            ( ( ITYPE == 4 .OR. ITYPE == 5 ) .AND. KL /= KU ) ) THEN
          INFO = -3
       ELSE IF ( ( ITYPE == 4 .AND. LDA < KL+1 ) .OR. &
                 ( ITYPE == 5 .AND. LDA < KU+1 ) .OR. &
                 ( ITYPE == 6 .AND. LDA < 2*KL+KU+1 ) ) THEN
          INFO = -9
       END IF
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SLASCL', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( N == 0 .OR. M == 0 ) RETURN

    ! Get machine parameters

    SMLNUM = SLAMCH( 'S' )
    BIGNUM = ONE / SMLNUM

    CFROMC = CFROM
    CTOC = CTO

    DONE = .FALSE.
    DO WHILE (.NOT. DONE)

       CFROM1 = CFROMC*SMLNUM
       IF ( IEEE_QUIET_EQ(CFROM1,CFROMC) ) THEN
          ! CFROMC is an inf.  Multiply by a correctly signed zero for
          ! finite CTOC, or a NaN if CTOC is infinite.
          MUL = CTOC / CFROMC
          DONE = .TRUE.
          CTO1 = CTOC
       ELSE
          CTO1 = CTOC / BIGNUM
          IF ( IEEE_QUIET_EQ(CTO1,CTOC) ) THEN
             ! CTOC is either 0 or an inf.  In both cases, CTOC itself
             ! serves as the correct multiplication factor.
             MUL = CTOC
             DONE = .TRUE.
             CFROMC = ONE
          ELSE IF ( ABS( CFROM1 ) > ABS( CTOC ) .AND. ABS( CTOC ) >= EPS ) THEN
             MUL = SMLNUM
             DONE = .FALSE.
             CFROMC = CFROM1
          ELSE IF ( ABS( CTO1 ) > ABS( CFROMC ) ) THEN
             MUL = BIGNUM
             DONE = .FALSE.
             CTOC = CTO1
          ELSE
             MUL = CTOC / CFROMC
             DONE = .TRUE.
          END IF
       END IF

       IF ( ITYPE == 0 ) THEN

          ! Full matrix

          DO J = 1, N
             DO I = 1, M
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       ELSE IF ( ITYPE == 1 ) THEN

          ! Lower triangular matrix

          DO J = 1, N
             DO I = J, M
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       ELSE IF ( ITYPE == 2 ) THEN

          ! Upper triangular matrix

          DO J = 1, N
             DO I = 1, MIN( J, M )
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       ELSE IF ( ITYPE == 3 ) THEN

          ! Upper Hessenberg matrix

          DO J = 1, N
             DO I = 1, MIN( J+1, M )
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       ELSE IF ( ITYPE == 4 ) THEN

          ! Lower half of a symmetric band matrix

          K3 = KL + 1
          K4 = N + 1
          DO J = 1, N
             DO I = 1, MIN( K3, K4-J )
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       ELSE IF ( ITYPE == 5 ) THEN

          ! Upper half of a symmetric band matrix

          K1 = KU + 2
          K3 = KU + 1
          DO J = 1, N
             DO I = MAX( K1-J, 1 ), K3
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       ELSE IF ( ITYPE == 6 ) THEN

          ! Band matrix

          K1 = KL + KU + 2
          K2 = KL + 1
          K3 = 2*KL + KU + 1
          K4 = KL + KU + 1 + M
          DO J = 1, N
             DO I = MAX( K1-J, K2 ), MIN( K3, K4-J )
                A( I, J ) = A( I, J )*MUL
             END DO
          END DO

       END IF

    END DO

  END SUBROUTINE SLASCL

!=====================================================================
!
!  SLASET
!
!  Purpose
!  =======
!
!  SLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) REAL
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) REAL
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i /= j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!=====================================================================

  SUBROUTINE SLASET( UPLOC, M, N, ALPHA, BETA, A, LDA )

    IMPORT,  ONLY: LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: UPLOC
    INTEGER,      INTENT(IN) :: LDA, M, N
    REAL,         INTENT(IN) :: ALPHA, BETA

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* )

    ! .. Local Scalars ..
    INTEGER   :: I, J
    CHARACTER :: UPLO

    ! .. Executable Statements ..

    UPLO = UPLOC(1:1)

    IF ( LSAME( UPLO, 'U' ) ) THEN

       ! Set the strictly upper triangular or trapezoidal part of the
       ! array to ALPHA.

       DO J = 2, N
          DO I = 1, MIN( J-1, M )
             A( I, J ) = ALPHA
          END DO
       END DO

    ELSE IF ( LSAME( UPLO, 'L' ) ) THEN

       ! Set the strictly lower triangular or trapezoidal part of the
       ! array to ALPHA.

       DO J = 1, MIN( M, N )
          DO I = J + 1, M
             A( I, J ) = ALPHA
          END DO
       END DO

    ELSE

       ! Set the leading m-by-n submatrix to ALPHA.

       DO J = 1, N
          DO I = 1, M
             A( I, J ) = ALPHA
          END DO
       END DO
    END IF

    ! Set the first min(M,N) diagonal elements to BETA.

    DO I = 1, MIN( M, N )
       A( I, I ) = BETA
    END DO

  END SUBROUTINE SLASET

!=====================================================================
!
!  SLASSQ
!
!  Purpose
!  =======
!
!  SLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) REAL array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) REAL
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) REAL
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
!=====================================================================

  SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    INTEGER,  INTENT(IN)    :: INCX, N
    REAL, INTENT(INOUT) :: SCALE, SUMSQ

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: X( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER :: IX
    REAL    :: ABSXI

    ! .. Executable Statements ..

    IF ( N > 0 ) THEN
       DO IX = 1, 1 + ( N-1 )*INCX, INCX
          ABSXI = ABS( X( IX ) )
          IF ( ABSXI >=  EPS ) THEN
             IF ( SCALE < ABSXI ) THEN
                SUMSQ = 1. + SUMSQ*( SCALE / ABSXI )**2
                SCALE = ABSXI
             ELSE
                SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
             END IF
          END IF
       END DO
    END IF

  END SUBROUTINE SLASSQ

!=====================================================================
!
!  Purpose
!  =======
!
!  SORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q**T if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by SGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left
!          = 'R': apply Q or Q**T from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q**T (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) REAL array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          SGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGEQRF.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!=====================================================================

  SUBROUTINE SORM2R( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )

    IMPORT,  ONLY: LSAME, XERBLA, SLARF

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ONE = 1.0

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, NOTRAN
    INTEGER   :: I, I1, I2, I3, IC, JC, MI, NI, NQ
    REAL      :: AII
    CHARACTER :: SIDE, TRANS

    ! .. Executable Statements ..

    SIDE  = SIDEC (1:1)
    TRANS = TRANSC(1:1)

    ! Test the input arguments

    INFO   = 0
    LEFT   = LSAME( SIDE,  'L' )
    NOTRAN = LSAME( TRANS, 'N' )

    ! NQ is the order of Q

    IF ( LEFT ) THEN
       NQ = M
    ELSE
       NQ = N
    END IF

    IF ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
       INFO = -2
    ELSE IF ( M < 0 ) THEN
       INFO = -3
    ELSE IF ( N < 0 ) THEN
       INFO = -4
    ELSE IF ( K < 0 .OR. K > NQ ) THEN
       INFO = -5
    ELSE IF ( LDA < MAX( 1, NQ ) ) THEN
       INFO = -7
    ELSE IF ( LDC < MAX( 1, M ) ) THEN
       INFO = -10
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SORM2R', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) RETURN

    IF ( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
       I1 = 1
       I2 = K
       I3 = 1
    ELSE
       I1 = K
       I2 = 1
       I3 = -1
    END IF

    IF ( LEFT ) THEN
       NI = N
       JC = 1
    ELSE
       MI = M
       IC = 1
    END IF

    DO I = I1, I2, I3

       IF ( LEFT ) THEN
          ! H(i) is applied to C(i:m,1:n)
          MI = M - I + 1
          IC = I
       ELSE
          ! H(i) is applied to C(1:m,i:n)
          NI = N - I + 1
          JC = I
       END IF

       ! Apply H(i)
       AII = A( I, I )
       A( I, I ) = ONE
       CALL SLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), LDC, WORK )
       A( I, I ) = AII

    END DO

  END SUBROUTINE SORM2R

!=====================================================================
!
!  SORML2
!
!  Purpose
!  =======
!
!  SORML2 overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q**T if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by SGELQF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left
!          = 'R': apply Q or Q**T from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q**T (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) REAL array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          SGELQF in the first k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGELQF.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!=====================================================================

  SUBROUTINE SORML2( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )

    IMPORT,  ONLY: LSAME, XERBLA, SLARF

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL, PARAMETER :: ONE = 1.0

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, NOTRAN
    INTEGER   :: I, I1, I2, I3, IC, JC, MI, NI, NQ
    REAL      :: AII
    CHARACTER :: SIDE, TRANS

    ! .. Executable Statements ..

    SIDE  = SIDEC (1:1)
    TRANS = TRANSC(1:1)

    ! Test the input arguments

    INFO   = 0
    LEFT   = LSAME( SIDE,  'L' )
    NOTRAN = LSAME( TRANS, 'N' )

    ! NQ is the order of Q

    IF ( LEFT ) THEN
       NQ = M
    ELSE
       NQ = N
    END IF

    IF ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
       INFO = -2
    ELSE IF ( M < 0 ) THEN
       INFO = -3
    ELSE IF ( N < 0 ) THEN
       INFO = -4
    ELSE IF ( K < 0 .OR. K > NQ ) THEN
       INFO = -5
    ELSE IF ( LDA < MAX( 1, K ) ) THEN
       INFO = -7
    ELSE IF ( LDC < MAX( 1, M ) ) THEN
       INFO = -10
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SORML2', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) RETURN

    IF ( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
       I1 = 1
       I2 = K
       I3 = 1
    ELSE
       I1 = K
       I2 = 1
       I3 = -1
    END IF

    IF ( LEFT ) THEN
       NI = N
       JC = 1
    ELSE
       MI = M
       IC = 1
    END IF

    DO I = I1, I2, I3

       IF ( LEFT ) THEN
          ! H(i) is applied to C(i:m,1:n)
          MI = M - I + 1
          IC = I
       ELSE
          ! H(i) is applied to C(1:m,i:n)
          NI = N - I + 1
          JC = I
       END IF

       ! Apply H(i)

       AII = A( I, I )
       A( I, I ) = ONE
       CALL SLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ), C( IC, JC ), LDC, WORK )
       A( I, I ) = AII

    END DO

  END SUBROUTINE SORML2

!=====================================================================
!
!  SORMLQ
!
!  Purpose
!  =======
!
!  SORMLQ overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by SGELQF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) REAL array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          SGELQF in the first k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGELQF.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!=====================================================================

  SUBROUTINE SORMLQ( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, &
                     WORK, LWORK, INFO )

    IMPORT, ONLY: LSAME, XERBLA, ILAENV, SLARFB, SLARFT, SORML2

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, LWORK, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    INTEGER, PARAMETER :: NBMAX = 64
    INTEGER, PARAMETER :: LDT   = NBMAX+1

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, LQUERY, NOTRAN
    CHARACTER :: TRANST
    INTEGER   :: I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, LWKOPT, &
                 MI, NB, NBMIN, NI, NQ, NW
    CHARACTER :: SIDE, TRANS

    ! .. Local Arrays ..
    REAL :: T( LDT, NBMAX )

    ! .. Executable Statements ..

    SIDE  = SIDEC (1:1)
    TRANS = TRANSC(1:1)

    ! Test the input arguments

    INFO   = 0
    LEFT   = LSAME( SIDE,  'L' )
    NOTRAN = LSAME( TRANS, 'N' )
    LQUERY = ( LWORK == -1 )

    ! NQ is the order of Q and NW is the minimum dimension of WORK

    IF ( LEFT ) THEN
       NQ = M
       NW = N
    ELSE
       NQ = N
       NW = M
    END IF

    IF ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
       INFO = -2
    ELSE IF ( M < 0 ) THEN
       INFO = -3
    ELSE IF ( N < 0 ) THEN
       INFO = -4
    ELSE IF ( K < 0 .OR. K > NQ ) THEN
       INFO = -5
    ELSE IF ( LDA < MAX( 1, K ) ) THEN
       INFO = -7
    ELSE IF ( LDC < MAX( 1, M ) ) THEN
       INFO = -10
    ELSE IF ( LWORK < MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
       INFO = -12
    END IF

    IF ( INFO == 0 ) THEN
       ! Determine the block size.  NB may be at most NBMAX, where NBMAX
       ! is used to define the local array T.
       NB = MIN( NBMAX, ILAENV( 1, 'SORMLQ', SIDE // TRANS, M, N, K, -1 ) )
       LWKOPT = MAX( 1, NW )*NB
       WORK( 1 ) = REAL( LWKOPT )
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SORMLQ', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) THEN
       WORK( 1 ) = 1
       RETURN
    END IF

    NBMIN = 2
    LDWORK = NW

    IF ( NB > 1 .AND. NB < K ) THEN
       IWS = NW*NB
       IF ( LWORK < IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'SORMLQ', SIDE // TRANS, M, N, K, -1 ) )
       END IF
    ELSE
       IWS = NW
    END IF

    IF ( NB < NBMIN .OR. NB >= K ) THEN

       ! Use unblocked code
       CALL SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )

    ELSE

       ! Use blocked code

       IF ( ( LEFT .AND. NOTRAN ) .OR. &
            ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
          I1 = 1
          I2 = K
          I3 = NB
       ELSE
          I1 = ( ( K-1 ) / NB )*NB + 1
          I2 = 1
          I3 = -NB
       END IF

       IF ( LEFT ) THEN
          NI = N
          JC = 1
       ELSE
          MI = M
          IC = 1
       END IF

       IF ( NOTRAN ) THEN
          TRANST = 'T'
       ELSE
          TRANST = 'N'
       END IF

       DO I = I1, I2, I3
          IB = MIN( NB, K-I+1 )

          ! Form the triangular factor of the block reflector
          ! H = H(i) H(i+1) . . . H(i+ib-1)

          CALL SLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ), &
                       LDA, TAU( I ), T, LDT )

          IF ( LEFT ) THEN
             ! H or H**T is applied to C(i:m,1:n)
             MI = M - I + 1
             IC = I
          ELSE
             ! H or H**T is applied to C(1:m,i:n)
             NI = N - I + 1
             JC = I
          END IF

          ! Apply H or H**T
          CALL SLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB, &
                       A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK, LDWORK )
       END DO

    END IF

    WORK( 1 ) = REAL( LWKOPT )

  END SUBROUTINE SORMLQ

!=====================================================================
!
!  SORMQR
!
!  Purpose
!  =======
!
!  SORMQR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by SGEQRF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) REAL array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          SGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGEQRF.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!=====================================================================

  SUBROUTINE SORMQR( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, &
                     WORK, LWORK, INFO )

    IMPORT,  ONLY: LSAME, XERBLA, ILAENV, SLARFB, SLARFT, SORM2R

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, LWORK, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(OUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    INTEGER, PARAMETER :: NBMAX = 64
    INTEGER, PARAMETER :: LDT   = NBMAX+1

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, LQUERY, NOTRAN
    INTEGER   :: I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, LWKOPT, MI, &
                 NB, NBMIN, NI, NQ, NW
    CHARACTER :: SIDE, TRANS

    ! .. Local Arrays ..
    REAL :: T( LDT, NBMAX )

    ! .. Executable Statements ..

    SIDE  = SIDEC (1:1)
    TRANS = TRANSC(1:1)

    ! Test the input arguments

    INFO   = 0
    LEFT   = LSAME( SIDE,  'L' )
    NOTRAN = LSAME( TRANS, 'N' )
    LQUERY = ( LWORK == -1 )

    ! NQ is the order of Q and NW is the minimum dimension of WORK

    IF ( LEFT ) THEN
       NQ = M
       NW = N
    ELSE
       NQ = N
       NW = M
    END IF

    IF ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
       INFO = -2
    ELSE IF ( M < 0 ) THEN
       INFO = -3
    ELSE IF ( N < 0 ) THEN
       INFO = -4
    ELSE IF ( K < 0 .OR. K > NQ ) THEN
       INFO = -5
    ELSE IF ( LDA < MAX( 1, NQ ) ) THEN
       INFO = -7
    ELSE IF ( LDC < MAX( 1, M ) ) THEN
       INFO = -10
    ELSE IF ( LWORK < MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
       INFO = -12
    END IF

    IF ( INFO == 0 ) THEN
       ! Determine the block size.  NB may be at most NBMAX, where NBMAX
       ! is used to define the local array T.
       NB = MIN( NBMAX, ILAENV( 1, 'SORMQR', SIDE // TRANS, M, N, K, -1 ) )
       LWKOPT = MAX( 1, NW )*NB
       WORK( 1 ) = REAL( LWKOPT )
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'SORMQR', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) THEN
       WORK( 1 ) = 1
       RETURN
    END IF

    NBMIN = 2
    LDWORK = NW
    IF ( NB > 1 .AND. NB < K ) THEN
       IWS = NW*NB
       IF ( LWORK < IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'SORMQR', SIDE // TRANS, M, N, K, -1 ) )
       END IF
    ELSE
       IWS = NW
    END IF

    IF ( NB < NBMIN .OR. NB >= K ) THEN

       ! Use unblocked code
       CALL SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )

    ELSE

       ! Use blocked code

       IF ( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
          I1 = 1
          I2 = K
          I3 = NB
       ELSE
          I1 = ( ( K-1 ) / NB )*NB + 1
          I2 = 1
          I3 = -NB
       END IF

       IF ( LEFT ) THEN
          NI = N
          JC = 1
       ELSE
          MI = M
          IC = 1
       END IF

       DO I = I1, I2, I3
          IB = MIN( NB, K-I+1 )

          ! Form the triangular factor of the block reflector
          ! H = H(i) H(i+1) . . . H(i+ib-1)

          CALL SLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), LDA, &
                       TAU( I ), T, LDT )

          IF ( LEFT ) THEN
             ! H or H**T is applied to C(i:m,1:n)
             MI = M - I + 1
             IC = I
          ELSE
             ! H or H**T is applied to C(1:m,i:n)
             NI = N - I + 1
             JC = I
          END IF

          ! Apply H or H**T

          CALL SLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, IB, &
                       A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK, LDWORK )
       END DO

    END IF

    WORK( 1 ) = REAL( LWKOPT )

  END SUBROUTINE SORMQR

!=====================================================================
!
!  STRTRS
!
!  Purpose
!  =======
!
!  STRTRS solves a triangular system of the form
!
!     A * X = B  or  A**T * X = B,
!
!  where A is a triangular matrix of order N, and B is an N-by-NRHS
!  matrix.  A check is made to verify that A is nonsingular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
!          upper triangular part of the array A contains the upper
!          triangular matrix, and the strictly lower triangular part of
!          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
!          triangular part of the array A contains the lower triangular
!          matrix, and the strictly upper triangular part of A is not
!          referenced.  If DIAG = 'U', the diagonal elements of A are
!          also not referenced and are assumed to be 1.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, if INFO = 0, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, the i-th diagonal element of A is zero,
!               indicating that the matrix is singular and the solutions
!               X have not been computed.
!
!=====================================================================

  SUBROUTINE STRTRS( UPLOC, TRANSC, DIAGC, N, NRHS, A, LDA, B, LDB, INFO )

    IMPORT, ONLY: LSAME, XERBLA, STRSM

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: DIAGC, TRANSC, UPLOC
    INTEGER,      INTENT(IN)  :: LDA, LDB, N, NRHS
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A( LDA,* ), B( LDB,* )

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL   :: NOUNIT
    CHARACTER :: DIAG, TRANS, UPLO

    ! .. Executable Statements ..

    DIAG  = DIAGC (1:1)
    TRANS = TRANSC(1:1)
    UPLO  = UPLOC (1:1)
    ! Test the input parameters.

    INFO = 0
    NOUNIT = LSAME( DIAG, 'N' )

    IF ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. &
                   LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
       INFO = -2
    ELSE IF ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
       INFO = -3
    ELSE IF ( N < 0 ) THEN
       INFO = -4
    ELSE IF ( NRHS < 0 ) THEN
       INFO = -5
    ELSE IF ( LDA < MAX( 1, N ) ) THEN
       INFO = -7
    ELSE IF ( LDB < MAX( 1, N ) ) THEN
       INFO = -9
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'STRTRS', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( N == 0 ) RETURN

    ! Check for singularity.

    IF ( NOUNIT ) THEN
       DO INFO = 1, N
          IF ( ABS( A(INFO,INFO) ) < EPS ) RETURN
       END DO
    END IF
    INFO = 0

    ! Solve A * x = b  or  A**T * x = b.

    CALL STRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB )

  END SUBROUTINE STRTRS

!===================================================================================
!
!  SGEMM
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!  Purpose:
!  =============
!
! SGEMM  performs one of the matrix-matrix operations
!
!    C := alpha*op( A )*op( B ) + beta*C,
!
! where  op( X ) is one of
!
!    op( X ) = X   or   op( X ) = X**T,
!
! alpha and beta are scalars, and A, B and C are matrices, with op( A )
! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments:
!  ==========
!
! param[in] TRANSA

!          TRANSA is CHARACTER*1
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A**T.
!
!              TRANSA = 'C' or 'c',  op( A ) = A**T.
!
! param[in] TRANSB

!          TRANSB is CHARACTER*1
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B**T.
!
!              TRANSB = 'C' or 'c',  op( B ) = B**T.
!
! param[in] M

!          M is INTEGER
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!
! param[in] N

!          N is INTEGER
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!
! param[in] K

!          K is INTEGER
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!
! param[in] ALPHA

!          ALPHA is REAL
!           On entry, ALPHA specifies the scalar alpha.
!
! param[in] A

!          A is REAL array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!
! param[in] LDA

!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!
! param[in] B

!          B is REAL array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!
! param[in] LDB

!          LDB is INTEGER
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!
! param[in] BETA

!          BETA is REAL
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!
! param[in,out] C

!          C is REAL array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
! param[in] LDC

!          LDC is INTEGER
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!  Level 3 Blas routine.
!
! Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!=====================================================================

  SUBROUTINE SGEMM(TRANSAC,TRANSBC,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    IMPORT,  ONLY: LSAME, XERBLA
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level3 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011


    ! .. Scalar Arguments ..
    REAL,         INTENT(IN) :: ALPHA, BETA
    INTEGER,      INTENT(IN) :: K, LDA, LDB, LDC, M, N
    CHARACTER(*), INTENT(IN) :: TRANSAC, TRANSBC

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A(LDA,*),B(LDB,*),C(LDC,*)

    ! .. Local Scalars ..
    REAL      :: TEMP
    INTEGER   :: I,INFO,J,L,NCOLA,NROWA,NROWB
    LOGICAL   :: NOTA,NOTB
    CHARACTER :: TRANSA,TRANSB

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    TRANSA = TRANSAC(1:1)
    TRANSB = TRANSBC(1:1)

    ! Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    ! transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    ! and  columns of  A  and the  number of  rows  of  B  respectively.

    NOTA = LSAME(TRANSA,'N')
    NOTB = LSAME(TRANSB,'N')

    IF (NOTA) THEN
       NROWA = M
       NCOLA = K
    ELSE
       NROWA = K
       NCOLA = M
    END IF

    IF (NOTB) THEN
       NROWB = K
    ELSE
       NROWB = N
    END IF

    ! Test the input parameters.

    INFO = 0
    IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
                           (.NOT.LSAME(TRANSA,'T'))) THEN
       INFO = 1
    ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
                               (.NOT.LSAME(TRANSB,'T'))) THEN
       INFO = 2
    ELSE IF (M < 0) THEN
       INFO = 3
    ELSE IF (N < 0) THEN
       INFO = 4
    ELSE IF (K < 0) THEN
       INFO = 5
    ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 8
    ELSE IF (LDB < MAX(1,NROWB)) THEN
       INFO = 10
    ELSE IF (LDC < MAX(1,M)) THEN
       INFO = 13
    END IF

    IF (INFO /= 0) THEN
       CALL XERBLA('SGEMM ',INFO)
       RETURN
    END IF

    ! Quick return if possible.

    IF ( (M == 0) .OR. (N == 0) .OR. &
         ( ((ABS(ALPHA) < EPS) .OR. (K == 0)) .AND. &
           (ABS(BETA-ONE) < SPACING(ONE)) ) ) RETURN

    ! And if  alpha == zero.

    IF ( ABS(ALPHA) < EPS ) THEN

       IF ( ABS(BETA) < EPS ) THEN
          DO J = 1,N
             DO I = 1,M
                C(I,J) = ZERO
             END DO
          END DO
       ELSE
          DO J = 1,N
             DO I = 1,M
                C(I,J) = BETA*C(I,J)
             END DO
          END DO
       END IF

       RETURN
    END IF

    ! Start the operations.

    IF (NOTB) THEN

       IF (NOTA) THEN

          ! Form  C := alpha*A*B + beta*C.

          DO J = 1,N
             IF ( ABS(BETA) < EPS ) THEN
                DO I = 1,M
                   C(I,J) = ZERO
                END DO
             ELSE IF ( ABS(BETA-ONE) >= SPACING(ONE) ) THEN
                DO I = 1,M
                   C(I,J) = BETA*C(I,J)
                END DO
             END IF
             DO L = 1,K
                IF ( ABS(B(L,J)) >= EPS ) THEN
                   TEMP = ALPHA*B(L,J)
                   DO I = 1,M
                      C(I,J) = C(I,J) + TEMP*A(I,L)
                   END DO
                END IF
             END DO
          END DO

       ELSE

          ! Form  C := alpha*A**T*B + beta*C

          DO J = 1,N
             DO I = 1,M
                TEMP = ZERO
                DO L = 1,K
                   TEMP = TEMP + A(L,I)*B(L,J)
                END DO
                IF ( ABS(BETA) < EPS ) THEN
                   C(I,J) = ALPHA*TEMP
                ELSE
                   C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                END IF
             END DO
          END DO
       END IF

    ELSE

       IF (NOTA) THEN

          ! Form  C := alpha*A*B**T + beta*C

          DO J = 1,N
             IF ( ABS(BETA) < EPS ) THEN
                DO I = 1,M
                   C(I,J) = ZERO
                END DO
             ELSE IF ( ABS(BETA-ONE) >= SPACING(ONE) ) THEN
                DO I = 1,M
                   C(I,J) = BETA*C(I,J)
                END DO
             END IF
             DO L = 1,K
                IF ( ABS(B(J,L)) >= EPS ) THEN
                   TEMP = ALPHA*B(J,L)
                   DO I = 1,M
                      C(I,J) = C(I,J) + TEMP*A(I,L)
                   END DO
                END IF
             END DO
          END DO

       ELSE

          ! Form  C := alpha*A**T*B**T + beta*C

          DO J = 1,N
             DO I = 1,M
                TEMP = ZERO
                DO L = 1,K
                   TEMP = TEMP + A(L,I)*B(J,L)
                END DO
                IF ( ABS(BETA) < EPS ) THEN
                   C(I,J) = ALPHA*TEMP
                ELSE
                   C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                END IF
             END DO
          END DO

       END IF

    END IF

  END SUBROUTINE SGEMM

!=====================================================================
!
!  SGEMV
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*),Y(*)
!       ..
!
!  Purpose:
!  =============
!
! SGEMV  performs one of the matrix-vector operations
!
!    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!
! where alpha and beta are scalars, x and y are vectors and A is an
! m by n matrix.
!
!  Arguments:
!  ==========
!
! param[in] TRANS

!          TRANS is CHARACTER*1
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!
! param[in] M

!          M is INTEGER
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!
! param[in] N

!          N is INTEGER
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!
! param[in] ALPHA

!          ALPHA is REAL
!           On entry, ALPHA specifies the scalar alpha.
!
! param[in] A

!          A is REAL array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!
! param[in] LDA

!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!
! param[in] X

!          X is REAL array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!
! param[in] INCX

!          INCX is INTEGER
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!
! param[in] BETA

!          BETA is REAL
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!
! param[in,out] Y

!          Y is REAL array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
! param[in] INCY

!          INCY is INTEGER
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!  Level 2 Blas routine.
!  The vector and matrix arguments are not referenced when N = 0, or M = 0
!
! Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!=====================================================================

  SUBROUTINE SGEMV(TRANSC,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

    IMPORT,  ONLY: LSAME, XERBLA
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level2 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL,         INTENT(IN) :: ALPHA,BETA
    INTEGER,      INTENT(IN) :: INCX,INCY,LDA,M,N
    CHARACTER(*), INTENT(IN) :: TRANSC

    ! .. Array Arguments ..
    REAL,  INTENT(INOUT) :: A(LDA,*),X(*),Y(*)

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL      :: TEMP
    INTEGER   :: I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
    CHARACTER :: TRANS

    TRANS = TRANSC(1:1)

    ! Test the input parameters.

    INFO = 0
    IF ( .NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
         .NOT.LSAME(TRANS,'C') ) THEN
       INFO = 1
    ELSE IF (M < 0) THEN
       INFO = 2
    ELSE IF (N < 0) THEN
       INFO = 3
    ELSE IF (LDA < MAX(1,M)) THEN
       INFO = 6
    ELSE IF (INCX == 0) THEN
       INFO = 8
    ELSE IF (INCY == 0) THEN
       INFO = 11
    END IF

    IF (INFO /= 0) THEN
       CALL XERBLA('SGEMV ',INFO)
       RETURN
    END IF

    ! Quick return if possible.

    IF ( (M == 0) .OR. (N == 0) .OR. &
         ( (ABS(ALPHA) < EPS) .AND. (ABS(BETA-ONE) < SPACING(ONE)))) RETURN

    ! Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    ! up the start points in  X  and  Y.

    IF (LSAME(TRANS,'N')) THEN
       LENX = N
       LENY = M
    ELSE
       LENX = M
       LENY = N
    END IF

    IF (INCX > 0) THEN
       KX = 1
    ELSE
       KX = 1 - (LENX-1)*INCX
    END IF

    IF (INCY > 0) THEN
       KY = 1
    ELSE
       KY = 1 - (LENY-1)*INCY
    END IF

    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through A.

    ! First form  y := beta*y.

    IF ( ABS(BETA-ONE) >= SPACING(ONE) ) THEN
       IF (INCY == 1) THEN
          IF ( ABS(BETA) < EPS ) THEN
             DO I = 1,LENY
                Y(I) = ZERO
             END DO
          ELSE
             DO I = 1,LENY
                Y(I) = BETA*Y(I)
             END DO
          END IF
       ELSE
          IY = KY
          IF ( ABS(BETA) < EPS ) THEN
             DO I = 1,LENY
                Y(IY) = ZERO
                IY = IY + INCY
             END DO
          ELSE
             DO I = 1,LENY
                Y(IY) = BETA*Y(IY)
                IY = IY + INCY
             END DO
          END IF
       END IF
    END IF

    IF ( ABS(ALPHA) < EPS ) RETURN

    IF (LSAME(TRANS,'N')) THEN

       ! Form  y := alpha*A*x + y.

       JX = KX
       IF (INCY == 1) THEN
          DO J = 1,N
             IF ( ABS(X(JX)) >= EPS ) THEN
                TEMP = ALPHA*X(JX)
                DO I = 1,M
                   Y(I) = Y(I) + TEMP*A(I,J)
                END DO
             END IF
             JX = JX + INCX
          END DO
       ELSE
          DO J = 1,N
             IF ( ABS(X(JX)) >= EPS ) THEN
                TEMP = ALPHA*X(JX)
                IY = KY
                DO I = 1,M
                   Y(IY) = Y(IY) + TEMP*A(I,J)
                   IY = IY + INCY
                END DO
             END IF
             JX = JX + INCX
          END DO
       END IF
    ELSE

       ! Form  y := alpha*A**T*x + y.

       JY = KY
       IF (INCX == 1) THEN
          DO J = 1,N
             TEMP = ZERO
             DO I = 1,M
                TEMP = TEMP + A(I,J)*X(I)
             END DO
             Y(JY) = Y(JY) + ALPHA*TEMP
             JY = JY + INCY
          END DO
       ELSE
          DO J = 1,N
             TEMP = ZERO
             IX = KX
             DO I = 1,M
                TEMP = TEMP + A(I,J)*X(IX)
                IX = IX + INCX
             END DO
             Y(JY) = Y(JY) + ALPHA*TEMP
             JY = JY + INCY
          END DO
       END IF
    END IF

  END SUBROUTINE SGEMV

!=====================================================================
!
!  SGER
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*),Y(*)
!       ..
!
!  Purpose:
!  =============
!
! SGER   performs the rank 1 operation
!
!    A := alpha*x*y**T + A,
!
! where alpha is a scalar, x is an m element vector, y is an n element
! vector and A is an m by n matrix.
!
!  Arguments:
!  ==========
!
! param[in] M

!          M is INTEGER
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!
! param[in] N

!          N is INTEGER
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!
! param[in] ALPHA

!          ALPHA is REAL
!           On entry, ALPHA specifies the scalar alpha.
!
! param[in] X

!          X is REAL array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!
! param[in] INCX

!          INCX is INTEGER
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!
! param[in] Y

!          Y is REAL array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!
! param[in] INCY

!          INCY is INTEGER
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!
! param[in,out] A

!          A is REAL array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
! param[in] LDA

!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!  Level 2 Blas routine.
!
! Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!=====================================================================

  SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

    IMPORT,  ONLY: XERBLA
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level2 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL,    INTENT(IN) :: ALPHA
    INTEGER, INTENT(IN) :: INCX,INCY,LDA,M,N

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A(LDA,*),X(*),Y(*)

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL    :: TEMP
    INTEGER :: I,INFO,IX,J,JY,KX

    ! Test the input parameters.

    INFO = 0
    IF (M < 0) THEN
       INFO = 1
    ELSE IF (N < 0) THEN
       INFO = 2
    ELSE IF (INCX == 0) THEN
       INFO = 5
    ELSE IF (INCY == 0) THEN
       INFO = 7
    ELSE IF (LDA < MAX(1,M)) THEN
       INFO = 9
    END IF

    IF (INFO /= 0) THEN
       CALL XERBLA('SGER  ',INFO)
       RETURN
    END IF

    ! Quick return if possible.

    IF ((M == 0) .OR. (N == 0) .OR. (ABS(ALPHA) < EPS)) RETURN

    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through A.

    IF (INCY > 0) THEN
       JY = 1
    ELSE
       JY = 1 - (N-1)*INCY
    END IF

    IF (INCX == 1) THEN

       DO J = 1,N
          IF ( ABS(Y(JY)) >= EPS) THEN
             TEMP = ALPHA*Y(JY)
             DO I = 1,M
                A(I,J) = A(I,J) + X(I)*TEMP
             END DO
          END IF
          JY = JY + INCY
       END DO

    ELSE

       IF (INCX > 0) THEN
          KX = 1
       ELSE
          KX = 1 - (M-1)*INCX
       END IF
       DO J = 1,N
          IF ( ABS(Y(JY)) >= EPS ) THEN
             TEMP = ALPHA*Y(JY)
             IX = KX
             DO I = 1,M
                A(I,J) = A(I,J) + X(IX)*TEMP
                IX = IX + INCX
             END DO
          END IF
          JY = JY + INCY
       END DO

    END IF

  END SUBROUTINE SGER

!=====================================================================
!
!  SNRM2
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL X(*)
!       ..
!
!  Purpose:
!  =============
!
! SNRM2 returns the euclidean norm of a vector via the function
! name, so that
!
!    SNRM2 := sqrt( x'*x ).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
! This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to SLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!=====================================================================

  REAL FUNCTION SNRM2(N,X,INCX)

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: INCX, N

    ! .. Array Arguments ..
    REAL, INTENT(IN) :: X(*)

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL    :: ABSXI,NORM,SCALE,SSQ
    INTEGER :: IX

    IF (N < 1 .OR. INCX < 1) THEN

       NORM = ZERO

    ELSE IF (N == 1) THEN

       NORM = ABS(X(1))

    ELSE

       SCALE = ZERO
       SSQ = ONE

       ! The following loop is equivalent to this call to the LAPACK
       ! auxiliary routine:
       ! CALL SLASSQ( N, X, INCX, SCALE, SSQ )

       DO IX = 1,1 + (N-1)*INCX,INCX
          IF ( ABS(X(IX)) >= EPS) THEN
             ABSXI = ABS(X(IX))
             IF (SCALE < ABSXI) THEN
                SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                SCALE = ABSXI
             ELSE
                SSQ = SSQ + (ABSXI/SCALE)**2
             END IF
          END IF
       END DO
       NORM = SCALE*SQRT(SSQ)

    END IF

    SNRM2 = NORM

  END FUNCTION SNRM2

!=====================================================================
!
!  SSCAL
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSCAL(N,SA,SX,INCX)
!
!       .. Scalar Arguments ..
!       REAL SA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*)
!       ..
!
!  Purpose:
!  =============
!
!    scales a vector by a constant.
!    uses unrolled loops for increment equal to 1.
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx <= 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!=====================================================================

  SUBROUTINE SSCAL(N,SA,SX,INCX)

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL,    INTENT(IN) :: SA
    INTEGER, INTENT(IN) :: INCX, N

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: SX(*)

    ! .. Local Scalars ..
    INTEGER :: I,M,MP1,NINCX

    IF (N <= 0 .OR. INCX <= 0) RETURN

    IF (INCX == 1) THEN

       ! code for increment equal to 1

       ! clean-up loop
       M = MOD(N,5)
       IF (M /= 0) THEN
          DO I = 1,M
             SX(I) = SA*SX(I)
          END DO
          IF (N < 5) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,5
          SX(I) = SA*SX(I)
          SX(I+1) = SA*SX(I+1)
          SX(I+2) = SA*SX(I+2)
          SX(I+3) = SA*SX(I+3)
          SX(I+4) = SA*SX(I+4)
       END DO

    ELSE

       ! code for increment not equal to 1

       NINCX = N*INCX
       DO I = 1,NINCX,INCX
          SX(I) = SA*SX(I)
       END DO

    END IF

  END SUBROUTINE SSCAL

!=====================================================================
!
!  STRMM
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),B(LDB,*)
!       ..
!
!  Purpose:
!  =============
!
! STRMM  performs one of the matrix-matrix operations
!
!    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!    op( A ) = A   or   op( A ) = A**T.
!
!  Arguments:
!  ==========
!
! param[in] SIDE

!          SIDE is CHARACTER*1
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
! param[in] UPLO

!          UPLO is CHARACTER*1
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
! param[in] TRANSA

!          TRANSA is CHARACTER*1
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A**T.
!
!              TRANSA = 'C' or 'c'   op( A ) = A**T.
!
! param[in] DIAG

!          DIAG is CHARACTER*1
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
! param[in] M

!          M is INTEGER
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!
! param[in] N

!          N is INTEGER
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!
! param[in] ALPHA

!          ALPHA is REAL
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!
! param[in] A

!          A is REAL array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!
! param[in] LDA

!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!
! param[in,out] B

!          B is REAL array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
! param[in] LDB

!          LDB is INTEGER
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!  Level 3 Blas routine.
!
! Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!=====================================================================

  SUBROUTINE STRMM(SIDES,UPLOS,TRANSAS,DIAGS,M,N,ALPHA,A,LDA,B,LDB)

    IMPORT,  ONLY: LSAME, XERBLA
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level3 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL,         INTENT(IN) :: ALPHA
    INTEGER,      INTENT(IN) :: LDA,LDB,M,N
    CHARACTER(*), INTENT(IN) :: DIAGS,SIDES,TRANSAS,UPLOS

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A(LDA,*),B(LDB,*)

    ! .. Local Scalars ..
    REAL      :: TEMP
    INTEGER   :: I,INFO,J,K,NROWA
    LOGICAL   :: LSIDE,NOUNIT,UPPER
    CHARACTER :: DIAG,SIDE,TRANSA,UPLO

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    DIAG   = DIAGS  (1:1)
    SIDE   = SIDES  (1:1)
    TRANSA = TRANSAS(1:1)
    UPLO   = UPLOS  (1:1)

    ! Test the input parameters.

    LSIDE = LSAME(SIDE,'L')
    IF (LSIDE) THEN
       NROWA = M
    ELSE
       NROWA = N
    END IF
    NOUNIT = LSAME(DIAG,'N')
    UPPER = LSAME(UPLO,'U')

    INFO = 0
    IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
       INFO = 1
    ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
       INFO = 2
    ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
             (.NOT.LSAME(TRANSA,'T')) .AND. &
             (.NOT.LSAME(TRANSA,'C'))) THEN
       INFO = 3
    ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
       INFO = 4
    ELSE IF (M < 0) THEN
       INFO = 5
    ELSE IF (N < 0) THEN
       INFO = 6
    ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 9
    ELSE IF (LDB < MAX(1,M)) THEN
       INFO = 11
    END IF

    IF (INFO /= 0) THEN
       CALL XERBLA('STRMM ',INFO)
       RETURN
    END IF

    ! Quick return if possible.

    IF (M == 0 .OR. N == 0) RETURN

    ! And when  alpha == zero.

    IF ( ABS(ALPHA) < EPS ) THEN
       DO J = 1,N
          DO I = 1,M
             B(I,J) = ZERO
          END DO
       END DO
       RETURN
    END IF

    ! Start the operations.

    IF (LSIDE) THEN

       IF (LSAME(TRANSA,'N')) THEN

          ! Form  B := alpha*A*B.
          IF (UPPER) THEN
             DO J = 1,N
                DO K = 1,M
                   IF ( ABS( B(K,J) ) >= EPS ) THEN
                      TEMP = ALPHA*B(K,J)
                      DO I = 1,K - 1
                         B(I,J) = B(I,J) + TEMP*A(I,K)
                      END DO
                      IF (NOUNIT) TEMP = TEMP*A(K,K)
                      B(K,J) = TEMP
                   END IF
                END DO
             END DO
          ELSE
             DO J = 1,N
                DO K = M,1,-1
                   IF ( ABS( B(K,J) ) >= EPS ) THEN
                      TEMP = ALPHA*B(K,J)
                      B(K,J) = TEMP
                      IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                      DO I = K + 1,M
                         B(I,J) = B(I,J) + TEMP*A(I,K)
                      END DO
                   END IF
                END DO
             END DO
          END IF

       ELSE

          ! Form  B := alpha*A**T*B.
          IF (UPPER) THEN
             DO J = 1,N
                DO I = M,1,-1
                   TEMP = B(I,J)
                   IF (NOUNIT) TEMP = TEMP*A(I,I)
                   DO K = 1,I - 1
                      TEMP = TEMP + A(K,I)*B(K,J)
                   END DO
                   B(I,J) = ALPHA*TEMP
                END DO
             END DO
          ELSE
             DO J = 1,N
                DO I = 1,M
                   TEMP = B(I,J)
                   IF (NOUNIT) TEMP = TEMP*A(I,I)
                   DO K = I + 1,M
                      TEMP = TEMP + A(K,I)*B(K,J)
                   END DO
                   B(I,J) = ALPHA*TEMP
                END DO
             END DO
          END IF
       END IF

    ELSE

       IF (LSAME(TRANSA,'N')) THEN

          ! Form  B := alpha*B*A.
          IF (UPPER) THEN
             DO J = N,1,-1
                TEMP = ALPHA
                IF (NOUNIT) TEMP = TEMP*A(J,J)
                DO I = 1,M
                   B(I,J) = TEMP*B(I,J)
                END DO
                DO K = 1,J - 1
                   IF ( ABS( A(K,J) ) >= EPS ) THEN
                      TEMP = ALPHA*A(K,J)
                      DO I = 1,M
                         B(I,J) = B(I,J) + TEMP*B(I,K)
                      END DO
                   END IF
                END DO
             END DO
          ELSE
             DO J = 1,N
                TEMP = ALPHA
                IF (NOUNIT) TEMP = TEMP*A(J,J)
                DO I = 1,M
                   B(I,J) = TEMP*B(I,J)
                END DO
                DO K = J + 1,N
                   IF ( ABS( A(K,J) ) >= EPS ) THEN
                      TEMP = ALPHA*A(K,J)
                      DO I = 1,M
                         B(I,J) = B(I,J) + TEMP*B(I,K)
                      END DO
                   END IF
                END DO
             END DO
          END IF

       ELSE

          ! Form  B := alpha*B*A**T.
          IF (UPPER) THEN
             DO K = 1,N
                DO J = 1,K - 1
                   IF ( ABS( A(J,K) ) >= EPS ) THEN
                      TEMP = ALPHA*A(J,K)
                      DO I = 1,M
                         B(I,J) = B(I,J) + TEMP*B(I,K)
                      END DO
                   END IF
                END DO
                TEMP = ALPHA
                IF (NOUNIT) TEMP = TEMP*A(K,K)
                IF ( ABS(TEMP-ONE) >= SPACING(ONE) ) THEN
                   DO I = 1,M
                      B(I,K) = TEMP*B(I,K)
                   END DO
                END IF
             END DO
          ELSE
             DO K = N,1,-1
                DO J = K + 1,N
                   IF ( ABS( A(J,K) ) >= EPS ) THEN
                      TEMP = ALPHA*A(J,K)
                      DO I = 1,M
                         B(I,J) = B(I,J) + TEMP*B(I,K)
                      END DO
                   END IF
                END DO
                TEMP = ALPHA
                IF (NOUNIT) TEMP = TEMP*A(K,K)
                IF ( ABS(TEMP-ONE) >= SPACING(ONE) ) THEN
                   DO I = 1,M
                      B(I,K) = TEMP*B(I,K)
                   END DO
                END IF
             END DO
          END IF

       END IF

    END IF

  END SUBROUTINE STRMM

!=====================================================================
!
!  STRMV
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*)
!       ..
!
!  Purpose:
!  =============
!
! STRMV  performs one of the matrix-vector operations
!
!    x := A*x,   or   x := A**T*x,
!
! where x is an n element vector and  A is an n by n unit, or non-unit,
! upper or lower triangular matrix.
!
!  Arguments:
!  ==========
!
! param[in] UPLO

!          UPLO is CHARACTER*1
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
! param[in] TRANS

!          TRANS is CHARACTER*1
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A**T*x.
!
!              TRANS = 'C' or 'c'   x := A**T*x.
!
! param[in] DIAG

!          DIAG is CHARACTER*1
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
! param[in] N

!          N is INTEGER
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!
! param[in] A

!          A is REAL array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!
! param[in] LDA

!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!
! param[in,out] X

!          X is REAL array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
! param[in] INCX

!          INCX is INTEGER
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!  Level 2 Blas routine.
!  The vector and matrix arguments are not referenced when N = 0, or M = 0
!
! Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!=====================================================================

  SUBROUTINE STRMV(UPLOC,TRANSC,DIAGC,N,A,LDA,X,INCX)

    IMPORT,  ONLY: LSAME, XERBLA
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level2 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER,      INTENT(IN) :: INCX, LDA, N
    CHARACTER(*), INTENT(IN) :: DIAGC, TRANSC, UPLOC

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A(LDA,*),X(*)

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL      :: TEMP
    INTEGER   :: I,INFO,IX,J,JX,KX
    LOGICAL   :: NOUNIT
    CHARACTER :: DIAG,TRANS,UPLO

    DIAG  = DIAGC (1:1)
    TRANS = TRANSC(1:1)
    UPLO  = UPLOC (1:1)

    ! Test the input parameters.

    INFO = 0
    IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
       INFO = 1
    ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
             .NOT.LSAME(TRANS,'C')) THEN
       INFO = 2
    ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
       INFO = 3
    ELSE IF (N < 0) THEN
       INFO = 4
    ELSE IF (LDA < MAX(1,N)) THEN
       INFO = 6
    ELSE IF (INCX == 0) THEN
       INFO = 8
    END IF

    IF (INFO /= 0) THEN
       CALL XERBLA('STRMV ',INFO)
       RETURN
    END IF

    ! Quick return if possible.

    IF (N == 0) RETURN

    NOUNIT = LSAME(DIAG,'N')

    ! Set up the start point in X if the increment is not unity. This
    ! will be  ( N - 1 )*INCX  too small for descending loops.

    IF (INCX <= 0) THEN
       KX = 1 - (N-1)*INCX
    ELSE IF (INCX /= 1) THEN
       KX = 1
    END IF

    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through A.

    IF (LSAME(TRANS,'N')) THEN

       ! Form  x := A*x.

       IF (LSAME(UPLO,'U')) THEN
          IF (INCX == 1) THEN
             DO J = 1,N
                IF ( ABS( X(J) ) >= EPS ) THEN
                   TEMP = X(J)
                   DO I = 1,J - 1
                      X(I) = X(I) + TEMP*A(I,J)
                   END DO
                   IF (NOUNIT) X(J) = X(J)*A(J,J)
                END IF
             END DO
          ELSE
             JX = KX
             DO J = 1,N
                IF ( ABS( X(JX) ) >= EPS ) THEN
                   TEMP = X(JX)
                   IX = KX
                   DO I = 1,J - 1
                      X(IX) = X(IX) + TEMP*A(I,J)
                      IX = IX + INCX
                   END DO
                   IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                END IF
                JX = JX + INCX
             END DO
          END IF
       ELSE
          IF (INCX == 1) THEN
             DO J = N,1,-1
                IF ( ABS( X(J) ) >= EPS ) THEN
                   TEMP = X(J)
                   DO I = N,J + 1,-1
                      X(I) = X(I) + TEMP*A(I,J)
                   END DO
                   IF (NOUNIT) X(J) = X(J)*A(J,J)
                END IF
             END DO
          ELSE
             KX = KX + (N-1)*INCX
             JX = KX
             DO J = N,1,-1
                IF ( ABS( X(JX) ) >= EPS ) THEN
                   TEMP = X(JX)
                   IX = KX
                   DO I = N,J + 1,-1
                      X(IX) = X(IX) + TEMP*A(I,J)
                      IX = IX - INCX
                   END DO
                   IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                END IF
                JX = JX - INCX
             END DO
          END IF
       END IF

    ELSE

       ! Form  x := A**T*x.

       IF (LSAME(UPLO,'U')) THEN
          IF (INCX == 1) THEN
             DO J = N,1,-1
                TEMP = X(J)
                IF (NOUNIT) TEMP = TEMP*A(J,J)
                DO I = J - 1,1,-1
                   TEMP = TEMP + A(I,J)*X(I)
                END DO
                X(J) = TEMP
             END DO
          ELSE
             JX = KX + (N-1)*INCX
             DO J = N,1,-1
                TEMP = X(JX)
                IX = JX
                IF (NOUNIT) TEMP = TEMP*A(J,J)
                DO I = J - 1,1,-1
                   IX = IX - INCX
                   TEMP = TEMP + A(I,J)*X(IX)
                END DO
                X(JX) = TEMP
                JX = JX - INCX
             END DO
          END IF
       ELSE
          IF (INCX == 1) THEN
             DO J = 1,N
                TEMP = X(J)
                IF (NOUNIT) TEMP = TEMP*A(J,J)
                DO I = J + 1,N
                   TEMP = TEMP + A(I,J)*X(I)
                END DO
                X(J) = TEMP
             END DO
          ELSE
             JX = KX
             DO J = 1,N
                TEMP = X(JX)
                IX = JX
                IF (NOUNIT) TEMP = TEMP*A(J,J)
                DO I = J + 1,N
                   IX = IX + INCX
                   TEMP = TEMP + A(I,J)*X(IX)
                END DO
                X(JX) = TEMP
                JX = JX + INCX
             END DO
          END IF
       END IF

    END IF

  END SUBROUTINE STRMV

!=====================================================================
!
!  STRSM
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),B(LDB,*)
!       ..
!
!  Purpose:
!  =============
!
! STRSM  solves one of the matrix equations
!
!    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!    op( A ) = A   or   op( A ) = A**T.
!
! The matrix X is overwritten on B.
!
!  Arguments:
!  ==========
!
! param[in] SIDE

!          SIDE is CHARACTER*1
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
! param[in] UPLO

!          UPLO is CHARACTER*1
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
! param[in] TRANSA

!          TRANSA is CHARACTER*1
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A**T.
!
!              TRANSA = 'C' or 'c'   op( A ) = A**T.
!
! param[in] DIAG

!          DIAG is CHARACTER*1
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
! param[in] M

!          M is INTEGER
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!
! param[in] N

!          N is INTEGER
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!
! param[in] ALPHA

!          ALPHA is REAL
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!
! param[in] A

!          A is REAL array of DIMENSION ( LDA, k ),
!           where k is m when SIDE = 'L' or 'l'
!             and k is n when SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!
! param[in] LDA

!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!
! param[in,out] B

!          B is REAL array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
! param[in] LDB

!          LDB is INTEGER
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!  Further Details:
!=====================
!
!  Level 3 Blas routine.
!
! Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!=====================================================================

  SUBROUTINE STRSM(SIDEC,UPLOC,TRANSAC,DIAGC,M,N,ALPHA,A,LDA,B,LDB)

    IMPORT,  ONLY: LSAME, XERBLA
    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level3 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL,         INTENT(IN) :: ALPHA
    INTEGER,      INTENT(IN) :: LDA,LDB,M,N
    CHARACTER(*), INTENT(IN) :: DIAGC,SIDEC,TRANSAC,UPLOC

    ! .. Array Arguments ..
    REAL, INTENT(INOUT) :: A(LDA,*),B(LDB,*)

    ! .. Local Scalars ..
    REAL      :: TEMP
    INTEGER   :: I,INFO,J,K,NROWA
    LOGICAL   :: LSIDE,NOUNIT,UPPER
    CHARACTER :: DIAG,SIDE,TRANSA,UPLO

    ! .. Parameters ..
    REAL, PARAMETER :: ZERO = 0.0
    REAL, PARAMETER :: ONE  = 1.0
    REAL, PARAMETER :: EPS  = TINY( ZERO )

    DIAG   = DIAGC  (1:1)
    SIDE   = SIDEC  (1:1)
    TRANSA = TRANSAC(1:1)
    UPLO   = UPLOC  (1:1)

    ! Test the input parameters.

    LSIDE = LSAME(SIDE,'L')
    IF (LSIDE) THEN
       NROWA = M
    ELSE
       NROWA = N
    END IF

    NOUNIT = LSAME(DIAG,'N')
    UPPER = LSAME(UPLO,'U')

    INFO = 0
    IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
       INFO = 1
    ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
       INFO = 2
    ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
            (.NOT.LSAME(TRANSA,'T')) .AND. &
            (.NOT.LSAME(TRANSA,'C'))) THEN
       INFO = 3
    ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
       INFO = 4
    ELSE IF (M < 0) THEN
       INFO = 5
    ELSE IF (N < 0) THEN
       INFO = 6
    ELSE IF (LDA < MAX(1,NROWA)) THEN
       INFO = 9
    ELSE IF (LDB < MAX(1,M)) THEN
       INFO = 11
    END IF

    IF (INFO /= 0) THEN
       CALL XERBLA('STRSM ',INFO)
       RETURN
    END IF

    ! Quick return if possible.

    IF (M == 0 .OR. N == 0) RETURN

    ! And when  alpha == zero.

    IF (ABS(ALPHA) < EPS) THEN
       DO J = 1,N
          DO I = 1,M
             B(I,J) = ZERO
          END DO
       END DO
       RETURN
    END IF

    ! Start the operations.

    IF (LSIDE) THEN

       IF (LSAME(TRANSA,'N')) THEN

          ! Form  B := alpha*inv( A )*B.

          IF (UPPER) THEN

             DO J = 1,N
                IF (ABS(ALPHA-ONE) >= SPACING(ONE)) THEN
                   DO I = 1,M
                      B(I,J) = ALPHA*B(I,J)
                   END DO
                END IF
                DO K = M,1,-1
                   IF (ABS(B(K,J)) >= EPS) THEN
                      IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                      DO I = 1,K - 1
                         B(I,J) = B(I,J) - B(K,J)*A(I,K)
                      END DO
                   END IF
                END DO
             END DO

          ELSE

             DO J = 1,N
                IF (ABS(ALPHA-ONE) >= SPACING(ONE)) THEN
                   DO I = 1,M
                      B(I,J) = ALPHA*B(I,J)
                   END DO
                END IF
                DO K = 1,M
                   IF (ABS(B(K,J)) >= EPS) THEN
                      IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                      DO I = K + 1,M
                         B(I,J) = B(I,J) - B(K,J)*A(I,K)
                      END DO
                   END IF
                END DO
             END DO

          END IF

       ELSE

          ! Form  B := alpha*inv( A**T )*B.

          IF (UPPER) THEN
             DO J = 1,N
                DO I = 1,M
                   TEMP = ALPHA*B(I,J)
                   DO K = 1,I - 1
                      TEMP = TEMP - A(K,I)*B(K,J)
                   END DO
                   IF (NOUNIT) TEMP = TEMP/A(I,I)
                   B(I,J) = TEMP
                END DO
             END DO

          ELSE

             DO J = 1,N
                DO I = M,1,-1
                   TEMP = ALPHA*B(I,J)
                   DO K = I + 1,M
                      TEMP = TEMP - A(K,I)*B(K,J)
                   END DO
                   IF (NOUNIT) TEMP = TEMP/A(I,I)
                   B(I,J) = TEMP
                END DO
             END DO

          END IF

       END IF

    ELSE

       IF (LSAME(TRANSA,'N')) THEN

          ! Form  B := alpha*B*inv( A ).

          IF (UPPER) THEN

             DO J = 1,N
                IF (ABS(ALPHA-ONE) >= SPACING(ONE)) THEN
                   DO I = 1,M
                      B(I,J) = ALPHA*B(I,J)
                   END DO
                END IF
                DO K = 1,J - 1
                   IF (ABS(A(K,J)) >= EPS) THEN
                      DO I = 1,M
                         B(I,J) = B(I,J) - A(K,J)*B(I,K)
                      END DO
                   END IF
                END DO
                IF (NOUNIT) THEN
                   TEMP = ONE/A(J,J)
                   DO I = 1,M
                      B(I,J) = TEMP*B(I,J)
                   END DO
                END IF
             END DO

          ELSE

             DO J = N,1,-1
                IF (ABS(ALPHA-ONE) >= SPACING(ONE)) THEN
                   DO I = 1,M
                      B(I,J) = ALPHA*B(I,J)
                   END DO
                END IF
                DO K = J + 1,N
                   IF (ABS(A(K,J)) >= EPS) THEN
                      DO I = 1,M
                         B(I,J) = B(I,J) - A(K,J)*B(I,K)
                      END DO
                   END IF
                END DO
                IF (NOUNIT) THEN
                   TEMP = ONE/A(J,J)
                   DO I = 1,M
                      B(I,J) = TEMP*B(I,J)
                   END DO
                END IF
             END DO

          END IF

       ELSE

          ! Form  B := alpha*B*inv( A**T ).

          IF (UPPER) THEN

             DO K = N,1,-1
                IF (NOUNIT) THEN
                   TEMP = ONE/A(K,K)
                   DO I = 1,M
                      B(I,K) = TEMP*B(I,K)
                   END DO
                END IF
                DO J = 1,K - 1
                   IF (ABS(A(J,K)) >= EPS) THEN
                      TEMP = A(J,K)
                      DO I = 1,M
                         B(I,J) = B(I,J) - TEMP*B(I,K)
                      END DO
                   END IF
                END DO
                IF (ABS(ALPHA-ONE) >= SPACING(ONE)) THEN
                   DO I = 1,M
                      B(I,K) = ALPHA*B(I,K)
                   END DO
                END IF
             END DO

          ELSE

             DO K = 1,N
                IF (NOUNIT) THEN
                   TEMP = ONE/A(K,K)
                   DO I = 1,M
                      B(I,K) = TEMP*B(I,K)
                   END DO
                END IF
                DO J = K + 1,N
                   IF (ABS(A(J,K)) >= EPS) THEN
                      TEMP = A(J,K)
                      DO I = 1,M
                         B(I,J) = B(I,J) - TEMP*B(I,K)
                      END DO
                   END IF
                END DO
                IF (ABS(ALPHA-ONE) >= SPACING(ONE)) THEN
                   DO I = 1,M
                      B(I,K) = ALPHA*B(I,K)
                   END DO
                END IF
             END DO

          END IF

       END IF

    END IF

  END SUBROUTINE STRSM

End Module Lapack_Single
