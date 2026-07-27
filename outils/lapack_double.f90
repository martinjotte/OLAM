Module Lapack_Double

  use Lapack_Utils

  ! This module contains a subset of the Lapack/Blas double precision library routines
  ! used in OLAM

Contains

!==========================================================================
!
!  ILADLC scans a matrix for its last non-zero column.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download ILADLC + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlc.f
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADLC( M, N, A, LDA )
!
!        IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* )
!       ..
!
!  Purpose:
!  =============
!
! ILADLC scans A for its last non-zero column.
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
!          A is REAL(DP)         array, dimension (LDA,N)
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

  INTEGER FUNCTION ILADLC( M, N, A, LDA )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER,  INTENT(IN) :: M, N, LDA

    ! .. Array Arguments ..
    REAL(DP), INTENT(IN) :: A( LDA,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER :: I

    ! .. Executable Statements ..

    ! Quick test for the common case where one corner is non-zero.
    IF ( N == 0 ) THEN
       ILADLC = N
    ELSE IF ( ABS( A(1,N) ) >= EPS .OR. ABS( A(M,N) ) >= EPS ) THEN
       ILADLC = N
    ELSE
       ! Now scan each column from the end, returning with the first non-zero.
       DO ILADLC = N, 1, -1
          DO I = 1, M
             IF ( ABS( A(I,ILADLC) ) >= EPS ) RETURN
          END DO
       END DO
    END IF

  END FUNCTION ILADLC

!==========================================================================
!
!  ILADLR scans a matrix for its last non-zero row.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download ILADLR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladlr.f
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADLR( M, N, A, LDA )
!
!        IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* )
!       ..
!
!  Purpose:
!  =============
!
! ILADLR scans A for its last non-zero row.
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
!          A is REAL(DP)         array, dimension (LDA,N)
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

  INTEGER FUNCTION ILADLR( M, N, A, LDA )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: M, N, LDA

    ! .. Array Arguments ..
    REAL(DP), INTENT(IN) :: A( LDA,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER :: I, J

    ! .. Executable Statements ..

    ! Quick test for the common case where one corner is non-zero.
    IF ( M == 0 ) THEN
       ILADLR = M
    ELSE IF ( ABS( A(M,1) ) >= EPS .OR. ABS( A(M,N) ) >= EPS ) THEN
       ILADLR = M
    ELSE
       ! Scan up each column tracking the last zero row seen.
       ILADLR = 0
       DO J = 1, N
          I=M
          DO WHILE( (ABS( A(MAX(I,1),J) ) < EPS) .AND. (I >= 1) )
             I=I-1
          ENDDO
          ILADLR = MAX( ILADLR, I )
       END DO
    END IF

  END FUNCTION ILADLR

!==========================================================================
!
!  DLAMCH
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      REAL(DP)         FUNCTION DLAMCH( CMACH )
!
!       IMPLICIT NONE (EXTERNAL, TYPE)
!
!  Purpose:
!  =============
!
! DLAMCH determines double precision machine parameters.
!
!  Arguments:
!  ==========
!
! param[in] CMACH
!
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax
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

  REAL(DP) FUNCTION DLAMCH( CMACHS )

    USE IEEE_ARITHMETIC, ONLY: IEEE_QUIET_EQ
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: CMACHS

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    REAL(DP)  :: RND, EPS, SFMIN, SMALL, RMACH
    CHARACTER :: CMACH

    ! .. Executable Statements ..

    CMACH = CMACHS(1:1)

    ! Assume rounding, not chopping. Always.

    RND = ONE

    IF ( IEEE_QUIET_EQ( ONE, RND ) ) THEN
       EPS = EPSILON(ZERO) * 0.5_DP
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

    DLAMCH = RMACH

  END FUNCTION DLAMCH

!===================================================================================
!
!  DGESVD computes the singular value decomposition (SVD) for GE matrices
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGESVD + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvd.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU, JOBVT
!       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA,* ), S( * ), U( LDU,* ), &
!                          VT( LDVT,* ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGESVD computes the singular value decomposition (SVD) of a real
! M-by-N matrix A, optionally computing the left and/or right singular
! vectors. The SVD is written
!
!      A = U * SIGMA * transpose(V)
!
! where SIGMA is an M-by-N matrix which is zero except for its
! min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
! V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
! are the singular values of A; they are real and non-negative, and
! are returned in descending order.  The first min(m,n) columns of
! U and V are the left and right singular vectors of A.
!
! Note that the routine returns V**T, not V.
!
!  Arguments:
!  ==========
!
! param[in] JOBU
!
!          JOBU is CHARACTER*1
!          Specifies options for computing all or part of the matrix U:
!          = 'A':  all M columns of U are returned in array U:
!          = 'S':  the first min(m,n) columns of U (the left singular
!                  vectors) are returned in the array U;
!          = 'O':  the first min(m,n) columns of U (the left singular
!                  vectors) are overwritten on the array A;
!          = 'N':  no columns of U (no left singular vectors) are
!                  computed.
!
! param[in] JOBVT
!
!          JOBVT is CHARACTER*1
!          Specifies options for computing all or part of the matrix
!          V**T:
!          = 'A':  all N rows of V**T are returned in the array VT;
!          = 'S':  the first min(m,n) rows of V**T (the right singular
!                  vectors) are returned in the array VT;
!          = 'O':  the first min(m,n) rows of V**T (the right singular
!                  vectors) are overwritten on the array A;
!          = 'N':  no rows of V**T (no right singular vectors) are
!                  computed.
!
!          JOBVT and JOBU cannot both be 'O'.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the input matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the input matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!          if JOBU = 'O',  A is overwritten with the first min(m,n)
!                          columns of U (the left singular vectors,
!                          stored columnwise);
!          if JOBVT = 'O', A is overwritten with the first min(m,n)
!                          rows of V**T (the right singular vectors,
!                          stored rowwise);
!          if JOBU /= 'O' and JOBVT /= 'O', the contents of A
!                          are destroyed.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] S
!
!          S is DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A, sorted so that S(i) >= S(i+1).
!
! param[out] U
!
!          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
!          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
!          if JOBU = 'S', U contains the first min(m,n) columns of U
!          (the left singular vectors, stored columnwise);
!          if JOBU = 'N' or 'O', U is not referenced.
!
! param[in] LDU
!
!          LDU is INTEGER
!          The leading dimension of the array U.  LDU >= 1; if
!          JOBU = 'S' or 'A', LDU >= M.
!
! param[out] VT
!
!          VT is DOUBLE PRECISION array, dimension (LDVT,N)
!          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
!          V**T;
!          if JOBVT = 'S', VT contains the first min(m,n) rows of
!          V**T (the right singular vectors, stored rowwise);
!          if JOBVT = 'N' or 'O', VT is not referenced.
!
! param[in] LDVT
!
!          LDVT is INTEGER
!          The leading dimension of the array VT.  LDVT >= 1; if
!          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
!
! param[out] WORK
!
!          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
!          superdiagonal elements of an upper bidiagonal matrix B
!          whose diagonal is in S (not necessarily sorted). B
!          satisfies A = U * B * VT, so it has the same singular values
!          as A, and singular vectors related by U and VT.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The dimension of the array WORK.
!          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
!             - PATH 1  (M much larger than N, JOBU='N')
!             - PATH 1t (N much larger than M, JOBVT='N')
!          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if DBDSQR did not converge, INFO specifies how many
!                superdiagonals of an intermediate bidiagonal form B
!                did not converge to zero. See the description of WORK
!                above for details.
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! April 2012
!
!=====================================================================

  SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, &
                     VT, LDVT, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64

    IMPORT, ONLY: DBDSQR, DGEBRD, DGELQF, DGEMM,  DGEQRF, DLACPY, &
                  DLASCL, DLASET, DORGBR, DORGLQ, DORGQR, DORMBR, &
                  XERBLA, LSAME,  ILAENV, DLAMCH, DLANGE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK driver routine (version 3.4.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: JOBU, JOBVT
    INTEGER,   INTENT(IN)  :: LDA, LDU, LDVT, LWORK, M, N
    INTEGER,   INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), S( * ), U( LDU,* ), &
                               VT( LDVT,* ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    LOGICAL  :: LQUERY, WNTUA, WNTUAS, WNTUN, WNTUO, WNTUS, &
                WNTVA, WNTVAS, WNTVN, WNTVO, WNTVS
    INTEGER  :: BDSPAC, BLK, CHUNK, I, IE, IERR, IR, ISCL,     &
                ITAU, ITAUP, ITAUQ, IU, IWORK, LDWRKR, LDWRKU, &
                MAXWRK, MINMN, MINWRK, MNTHR, NCU, NCVT, NRU,  &
                NRVT, WRKBL
    INTEGER  :: LWORK_DGEQRF, LWORK_DORGQR_N, LWORK_DORGQR_M, &
                LWORK_DGEBRD, LWORK_DORGBR_P, LWORK_DORGBR_Q, &
                LWORK_DGELQF, LWORK_DORGLQ_N, LWORK_DORGLQ_M
    REAL(DP) :: ANRM, BIGNUM, EPS, SMLNUM

    ! .. Local Arrays ..
    REAL(DP) :: DUM( 1 )

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO   = 0
    MINMN  = MIN( M, N )
    WNTUA  = LSAME( JOBU, 'A' )
    WNTUS  = LSAME( JOBU, 'S' )
    WNTUAS = WNTUA .OR. WNTUS
    WNTUO  = LSAME( JOBU, 'O' )
    WNTUN  = LSAME( JOBU, 'N' )
    WNTVA  = LSAME( JOBVT, 'A' )
    WNTVS  = LSAME( JOBVT, 'S' )
    WNTVAS = WNTVA .OR. WNTVS
    WNTVO  = LSAME( JOBVT, 'O' )
    WNTVN  = LSAME( JOBVT, 'N' )
    LQUERY = ( LWORK == -1 )

    IF ( .NOT.( WNTUA .OR. WNTUS .OR. WNTUO .OR. WNTUN ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.( WNTVA .OR. WNTVS .OR. WNTVO .OR. WNTVN ) .OR. &
              ( WNTVO .AND. WNTUO ) ) THEN
       INFO = -2
    ELSE IF ( M < 0 ) THEN
       INFO = -3
    ELSE IF ( N < 0 ) THEN
       INFO = -4
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -6
    ELSE IF ( LDU < 1 .OR. ( WNTUAS .AND. LDU < M ) ) THEN
       INFO = -9
    ELSE IF ( LDVT < 1 .OR. ( WNTVA .AND. LDVT < N ) .OR. &
              ( WNTVS .AND. LDVT < MINMN ) ) THEN
       INFO = -11
    END IF

    ! Compute workspace
    ! (Note: Comments in the code beginning "Workspace:" describe the
    ! minimal amount of workspace needed at that point in the code,
    ! as well as the preferred amount for good performance.
    ! NB refers to the optimal block size for the immediately
    ! following subroutine, as returned by ILAENV.)

    IF ( INFO == 0 ) THEN
       MINWRK = 1
       MAXWRK = 1
       IF ( M >= N .AND. MINMN > 0 ) THEN

          ! Compute space needed for DBDSQR

          MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 )
          BDSPAC = 5*N
          ! Compute space needed for DGEQRF
          CALL DGEQRF( M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DGEQRF = INT( DUM(1) )
          ! Compute space needed for DORGQR
          CALL DORGQR( M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGQR_N = INT( DUM(1) )
          CALL DORGQR( M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGQR_M = INT( DUM(1) )
          ! Compute space needed for DGEBRD
          CALL DGEBRD( N, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
          LWORK_DGEBRD = INT( DUM(1) )
          ! Compute space needed for DORGBR P &
          CALL DORGBR( 'P', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGBR_P = INT( DUM(1) )
          ! Compute space needed for DORGBR Q
          CALL DORGBR( 'Q', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGBR_Q = INT( DUM(1) )

          IF ( M >= MNTHR ) THEN
             IF ( WNTUN ) THEN

                ! Path 1 (M much larger than N, JOBU='N')

                MAXWRK = N + LWORK_DGEQRF
                MAXWRK = MAX( MAXWRK, 3*N+LWORK_DGEBRD )
                IF ( WNTVO .OR. WNTVAS ) &
                     MAXWRK = MAX( MAXWRK, 3*N+LWORK_DORGBR_P )
                MAXWRK = MAX( MAXWRK, BDSPAC )
                MINWRK = MAX( 4*N, BDSPAC )

             ELSE IF ( WNTUO .AND. WNTVN ) THEN

                ! Path 2 (M much larger than N, JOBU='O', JOBVT='N')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_N )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = MAX( N*N+WRKBL, N*N+M*N+N )
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUO .AND. WNTVAS ) THEN

                ! Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
                ! 'A')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_N )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = MAX( N*N+WRKBL, N*N+M*N+N )
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUS .AND. WNTVN ) THEN

                ! Path 4 (M much larger than N, JOBU='S', JOBVT='N')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_N )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = N*N + WRKBL
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUS .AND. WNTVO ) THEN

                ! Path 5 (M much larger than N, JOBU='S', JOBVT='O')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_N )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = 2*N*N + WRKBL
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUS .AND. WNTVAS ) THEN

                ! Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
                ! 'A')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_N )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = N*N + WRKBL
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUA .AND. WNTVN ) THEN

                ! Path 7 (M much larger than N, JOBU='A', JOBVT='N')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_M )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = N*N + WRKBL
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUA .AND. WNTVO ) THEN

                ! Path 8 (M much larger than N, JOBU='A', JOBVT='O')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_M )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = 2*N*N + WRKBL
                MINWRK = MAX( 3*N+M, BDSPAC )

             ELSE IF ( WNTUA .AND. WNTVAS ) THEN

                ! Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
                ! 'A')

                WRKBL = N + LWORK_DGEQRF
                WRKBL = MAX( WRKBL, N+LWORK_DORGQR_M )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, 3*N+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = N*N + WRKBL
                MINWRK = MAX( 3*N+M, BDSPAC )

             END IF

          ELSE

             ! Path 10 (M at least N, but not much larger)

             CALL DGEBRD( M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
             LWORK_DGEBRD = INT( DUM(1) )
             MAXWRK = 3*N + LWORK_DGEBRD
             IF ( WNTUS .OR. WNTUO ) THEN
                CALL DORGBR( 'Q', M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR )
                LWORK_DORGBR_Q = INT( DUM(1) )
                MAXWRK = MAX( MAXWRK, 3*N+LWORK_DORGBR_Q )
             END IF
             IF ( WNTUA ) THEN
                CALL DORGBR( 'Q', M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
                LWORK_DORGBR_Q = INT( DUM(1) )
                MAXWRK = MAX( MAXWRK, 3*N+LWORK_DORGBR_Q )
             END IF
             IF ( .NOT.WNTVN ) THEN
                MAXWRK = MAX( MAXWRK, 3*N+LWORK_DORGBR_P )
             END IF
             MAXWRK = MAX( MAXWRK, BDSPAC )
             MINWRK = MAX( 3*N+M, BDSPAC )

          END IF

       ELSE IF ( MINMN > 0 ) THEN

          ! Compute space needed for DBDSQR

          MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 )
          BDSPAC = 5*M
          ! Compute space needed for DGELQF
          CALL DGELQF( M, N, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DGELQF = INT( DUM(1) )
          ! Compute space needed for DORGLQ
          CALL DORGLQ( N, N, M, DUM(1), N, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGLQ_N = INT( DUM(1) )
          CALL DORGLQ( M, N, M, A, LDA, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGLQ_M = INT( DUM(1) )
          ! Compute space needed for DGEBRD
          CALL DGEBRD( M, M, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
          LWORK_DGEBRD = INT( DUM(1) )
          ! Compute space needed for DORGBR P
          CALL DORGBR( 'P', M, M, M, A, N, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGBR_P = INT( DUM(1) )
          ! Compute space needed for DORGBR Q
          CALL DORGBR( 'Q', M, M, M, A, N, DUM(1), DUM(1), -1, IERR )
          LWORK_DORGBR_Q = INT( DUM(1) )

          IF ( N >= MNTHR ) THEN
             IF ( WNTVN ) THEN

                ! Path 1t(N much larger than M, JOBVT='N')

                MAXWRK = M + LWORK_DGELQF
                MAXWRK = MAX( MAXWRK, 3*M+LWORK_DGEBRD )
                IF ( WNTUO .OR. WNTUAS ) &
                     MAXWRK = MAX( MAXWRK, 3*M+LWORK_DORGBR_Q )
                MAXWRK = MAX( MAXWRK, BDSPAC )
                MINWRK = MAX( 4*M, BDSPAC )

             ELSE IF ( WNTVO .AND. WNTUN ) THEN

                ! Path 2t(N much larger than M, JOBU='N', JOBVT='O')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_M )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = MAX( M*M+WRKBL, M*M+M*N+M )
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVO .AND. WNTUAS ) THEN

                ! Path 3t(N much larger than M, JOBU='S' or 'A',
                ! JOBVT='O')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_M )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = MAX( M*M+WRKBL, M*M+M*N+M )
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVS .AND. WNTUN ) THEN

                ! Path 4t(N much larger than M, JOBU='N', JOBVT='S')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_M )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = M*M + WRKBL
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVS .AND. WNTUO ) THEN

                ! Path 5t(N much larger than M, JOBU='O', JOBVT='S')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_M )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = 2*M*M + WRKBL
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVS .AND. WNTUAS ) THEN

                ! Path 6t(N much larger than M, JOBU='S' or 'A',
                ! JOBVT='S')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_M )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = M*M + WRKBL
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVA .AND. WNTUN ) THEN

                ! Path 7t(N much larger than M, JOBU='N', JOBVT='A')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_N )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = M*M + WRKBL
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVA .AND. WNTUO ) THEN

                ! Path 8t(N much larger than M, JOBU='O', JOBVT='A')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_N )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = 2*M*M + WRKBL
                MINWRK = MAX( 3*M+N, BDSPAC )

             ELSE IF ( WNTVA .AND. WNTUAS ) THEN

                ! Path 9t(N much larger than M, JOBU='S' or 'A',
                ! JOBVT='A')

                WRKBL = M + LWORK_DGELQF
                WRKBL = MAX( WRKBL, M+LWORK_DORGLQ_N )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DGEBRD )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_P )
                WRKBL = MAX( WRKBL, 3*M+LWORK_DORGBR_Q )
                WRKBL = MAX( WRKBL, BDSPAC )
                MAXWRK = M*M + WRKBL
                MINWRK = MAX( 3*M+N, BDSPAC )

             END IF

          ELSE

             ! Path 10t(N greater than M, but not much larger)

             CALL DGEBRD( M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
             LWORK_DGEBRD = INT( DUM(1) )
             MAXWRK = 3*M + LWORK_DGEBRD
             IF ( WNTVS .OR. WNTVO ) THEN
                ! Compute space needed for DORGBR P
                CALL DORGBR( 'P', M, N, M, A, N, DUM(1), DUM(1), -1, IERR )
                LWORK_DORGBR_P = INT( DUM(1) )
                MAXWRK = MAX( MAXWRK, 3*M+LWORK_DORGBR_P )
             END IF
             IF ( WNTVA ) THEN
                CALL DORGBR( 'P', N, N, M, A, N, DUM(1), DUM(1), -1, IERR )
                LWORK_DORGBR_P = INT( DUM(1) )
                MAXWRK = MAX( MAXWRK, 3*M+LWORK_DORGBR_P )
             END IF
             IF ( .NOT.WNTUN ) THEN
                MAXWRK = MAX( MAXWRK, 3*M+LWORK_DORGBR_Q )
             END IF
             MAXWRK = MAX( MAXWRK, BDSPAC )
             MINWRK = MAX( 3*M+N, BDSPAC )

          END IF

       END IF

       MAXWRK = MAX( MAXWRK, MINWRK )
       WORK( 1 ) = REAL( MAXWRK, DP )

       IF ( LWORK < MINWRK .AND. .NOT.LQUERY ) THEN
          INFO = -13
       END IF

    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DGESVD', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 ) THEN
       RETURN
    END IF

    ! Get machine constants

    EPS = DLAMCH( 'P' )
    SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
    BIGNUM = ONE / SMLNUM

    ! Scale A if max element outside range [SMLNUM,BIGNUM]

    ANRM = DLANGE( 'M', M, N, A, LDA, DUM )
    ISCL = 0
    IF ( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN
       ISCL = 1
       CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR )
    ELSE IF ( ANRM > BIGNUM ) THEN
       ISCL = 1
       CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR )
    END IF

    IF ( M >= N ) THEN

       ! A has at least as many rows as columns. If A has sufficiently
       ! more rows than columns, first reduce using the QR
       ! decomposition (if sufficient workspace available)

       IF ( M >= MNTHR ) THEN

          IF ( WNTUN ) THEN

             ! Path 1 (M much larger than N, JOBU='N')
             ! No left singular vectors to be computed

             ITAU = 1
             IWORK = ITAU + N

             ! Compute A=Q*R
             ! (Workspace: need 2*N, prefer N+N*NB)

             CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )

             ! Zero out below R

             CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
             IE = 1
             ITAUQ = IE + N
             ITAUP = ITAUQ + N
             IWORK = ITAUP + N

             ! Bidiagonalize R in A
             ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

             CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                          WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
             NCVT = 0
             IF ( WNTVO .OR. WNTVAS ) THEN

                ! If right singular vectors desired, generate P'.
                ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )
                NCVT = N
             END IF
             IWORK = IE + N

             ! Perform bidiagonal QR iteration, computing right
             ! singular vectors of A in A if desired
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'U', N, NCVT, 0, 0, S, WORK( IE ), A, LDA, &
                          DUM, 1, DUM, 1, WORK( IWORK ), INFO )

             ! If right singular vectors desired in VT, copy them there

             IF ( WNTVAS ) CALL DLACPY( 'F', N, N, A, LDA, VT, LDVT )

          ELSE IF ( WNTUO .AND. WNTVN ) THEN

             ! Path 2 (M much larger than N, JOBU='O', JOBVT='N')
             ! N left singular vectors to be overwritten on A and
             ! no right singular vectors to be computed

             IF ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) THEN

                ! Sufficient workspace for a fast algorithm

                IR = 1

                IF ( LWORK >= MAX( WRKBL, LDA*N+N )+LDA*N ) THEN

                   ! WORK(IU) is LDA by N, WORK(IR) is LDA by N

                   LDWRKU = LDA
                   LDWRKR = LDA

                ELSE IF ( LWORK >= MAX( WRKBL, LDA*N+N )+N*N ) THEN

                   ! WORK(IU) is LDA by N, WORK(IR) is N by N

                   LDWRKU = LDA
                   LDWRKR = N

                ELSE

                   ! WORK(IU) is LDWRKU by N, WORK(IR) is N by N

                   LDWRKU = ( LWORK-N*N-N ) / N
                   LDWRKR = N

                END IF

                ITAU = IR + LDWRKR*N
                IWORK = ITAU + N

                ! Compute A=Q*R
                ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Copy R to WORK(IR) and zero out below it

                CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
                CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), LDWRKR )

                ! Generate Q in A
                ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                IE = ITAU
                ITAUQ = IE + N
                ITAUP = ITAUQ + N
                IWORK = ITAUP + N

                ! Bidiagonalize R in WORK(IR)
                ! (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ), &
                             WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate left vectors bidiagonalizing R
                ! (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, &
                             WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                IWORK = IE + N

                ! Perform bidiagonal QR iteration, computing left
                ! singular vectors of R in WORK(IR)
                ! (Workspace: need N*N+BDSPAC)

                CALL DBDSQR( 'U', N, 0, N, 0, S, WORK( IE ), DUM, 1, &
                             WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO )
                IU = IE + N

                ! Multiply Q in A by left singular vectors of R in
                ! WORK(IR), storing result in WORK(IU) and copying to A
                ! (Workspace: need N*N+2*N, prefer N*N+M*N+N)

                DO I = 1, M, LDWRKU
                   CHUNK = MIN( M-I+1, LDWRKU )
                   CALL DGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ), &
                               LDA, WORK( IR ), LDWRKR, ZERO, WORK( IU ), LDWRKU )
                   CALL DLACPY( 'F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA )
                END DO

             ELSE

                ! Insufficient workspace for a fast algorithm

                IE = 1
                ITAUQ = IE + N
                ITAUP = ITAUQ + N
                IWORK = ITAUP + N

                ! Bidiagonalize A
                ! (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)

                CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                             WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate left vectors bidiagonalizing A
                ! (Workspace: need 4*N, prefer 3*N+N*NB)

                CALL DORGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IWORK = IE + N

                ! Perform bidiagonal QR iteration, computing left
                ! singular vectors of A in A
                ! (Workspace: need BDSPAC)

                CALL DBDSQR( 'U', N, 0, M, 0, S, WORK( IE ), DUM, 1, A, LDA, &
                             DUM, 1, WORK( IWORK ), INFO )

             END IF

          ELSE IF ( WNTUO .AND. WNTVAS ) THEN

             ! Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
             ! N left singular vectors to be overwritten on A and
             ! N right singular vectors to be computed in VT

             IF ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) THEN

                ! Sufficient workspace for a fast algorithm

                IR = 1
                IF ( LWORK >= MAX( WRKBL, LDA*N+N )+LDA*N ) THEN

                   ! WORK(IU) is LDA by N and WORK(IR) is LDA by N

                   LDWRKU = LDA
                   LDWRKR = LDA

                ELSE IF ( LWORK >= MAX( WRKBL, LDA*N+N )+N*N ) THEN

                   ! WORK(IU) is LDA by N and WORK(IR) is N by N

                   LDWRKU = LDA
                   LDWRKR = N

                ELSE

                   ! WORK(IU) is LDWRKU by N and WORK(IR) is N by N

                   LDWRKU = ( LWORK-N*N-N ) / N
                   LDWRKR = N

                END IF

                ITAU = IR + LDWRKR*N
                IWORK = ITAU + N

                ! Compute A=Q*R
                ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Copy R to VT, zeroing out below it

                CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                IF ( N > 1 ) &
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT )

                ! Generate Q in A
                ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IE = ITAU
                ITAUQ = IE + N
                ITAUP = ITAUQ + N
                IWORK = ITAUP + N

                ! Bidiagonalize R in VT, copying result to WORK(IR)
                ! (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), &
                             WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                CALL DLACPY( 'L', N, N, VT, LDVT, WORK( IR ), LDWRKR )

                ! Generate left vectors bidiagonalizing R in WORK(IR)
                ! (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate right vectors bidiagonalizing R in VT
                ! (Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB)

                CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )
                IWORK = IE + N

                ! Perform bidiagonal QR iteration, computing left
                ! singular vectors of R in WORK(IR) and computing right
                ! singular vectors of R in VT
                ! (Workspace: need N*N+BDSPAC)

                CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), VT, LDVT, &
                             WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO )
                IU = IE + N

                ! Multiply Q in A by left singular vectors of R in
                ! WORK(IR), storing result in WORK(IU) and copying to A
                ! (Workspace: need N*N+2*N, prefer N*N+M*N+N)

                DO I = 1, M, LDWRKU
                   CHUNK = MIN( M-I+1, LDWRKU )
                   CALL DGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ), &
                               LDA, WORK( IR ), LDWRKR, ZERO, WORK( IU ), LDWRKU )
                   CALL DLACPY( 'F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA )
                END DO

             ELSE

                ! Insufficient workspace for a fast algorithm

                ITAU = 1
                IWORK = ITAU + N

                ! Compute A=Q*R
                ! (Workspace: need 2*N, prefer N+N*NB)

                CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Copy R to VT, zeroing out below it

                CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                IF ( N > 1 ) &
                     CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT )

                ! Generate Q in A
                ! (Workspace: need 2*N, prefer N+N*NB)

                CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )
                IE = ITAU
                ITAUQ = IE + N
                ITAUP = ITAUQ + N
                IWORK = ITAUP + N

                ! Bidiagonalize R in VT
                ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), &
                             WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Multiply Q in A by left vectors bidiagonalizing R
                ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                CALL DORMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), &
                             A, LDA, WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate right vectors bidiagonalizing R in VT
                ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )
                IWORK = IE + N

                ! Perform bidiagonal QR iteration, computing left
                ! singular vectors of A in A and computing right
                ! singular vectors of A in VT
                ! (Workspace: need BDSPAC)

                CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), VT, LDVT, &
                             A, LDA, DUM, 1, WORK( IWORK ), INFO )

             END IF

          ELSE IF ( WNTUS ) THEN

             IF ( WNTVN ) THEN

                ! Path 4 (M much larger than N, JOBU='S', JOBVT='N')
                ! N left singular vectors to be computed in U and
                ! no right singular vectors to be computed

                IF ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IR = 1
                   IF ( LWORK >= WRKBL+LDA*N ) THEN

                      ! WORK(IR) is LDA by N

                      LDWRKR = LDA
                   ELSE

                      ! WORK(IR) is N by N

                      LDWRKR = N
                   END IF
                   ITAU = IR + LDWRKR*N
                   IWORK = ITAU + N

                   ! Compute A=Q*R
                   ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Copy R to WORK(IR), zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), LDWRKR )

                   ! Generate Q in A
                   ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                   CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in WORK(IR)
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                   CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left vectors bidiagonalizing R in WORK(IR)
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                   CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of R in WORK(IR)
                   ! (Workspace: need N*N+BDSPAC)

                   CALL DBDSQR( 'U', N, 0, N, 0, S, WORK( IE ), DUM, 1, WORK( IR ), &
                                LDWRKR, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply Q in A by left singular vectors of R in
                   ! WORK(IR), storing result in U
                   ! (Workspace: need N*N)

                   CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, &
                               WORK( IR ), LDWRKR, ZERO, U, LDU )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DORGQR( M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Zero out below R in A

                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )

                   ! Bidiagonalize R in A
                   ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                   CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply Q in U by left vectors bidiagonalizing R
                   ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                   CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', N, 0, M, 0, S, WORK( IE ), DUM,   &
                                1, U, LDU, DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTVO ) THEN

                ! Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                ! N left singular vectors to be computed in U and
                ! N right singular vectors to be overwritten on A

                IF ( LWORK >= 2*N*N+MAX( 4*N, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+2*LDA*N ) THEN

                      ! WORK(IU) is LDA by N and WORK(IR) is LDA by N

                      LDWRKU = LDA
                      IR = IU + LDWRKU*N
                      LDWRKR = LDA

                   ELSE IF ( LWORK >= WRKBL+( LDA+N )*N ) THEN

                      ! WORK(IU) is LDA by N and WORK(IR) is N by N

                      LDWRKU = LDA
                      IR = IU + LDWRKU*N
                      LDWRKR = N

                   ELSE

                      ! WORK(IU) is N by N and WORK(IR) is N by N

                      LDWRKU = N
                      IR = IU + LDWRKU*N
                      LDWRKR = N

                   END IF

                   ITAU = IR + LDWRKR*N
                   IWORK = ITAU + N

                   ! Compute A=Q*R
                   ! (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), &
                               WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Copy R to WORK(IU), zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU )

                   ! Generate Q in A
                   ! (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)

                   CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in WORK(IU), copying result to WORK(IR)
                   ! (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+2*N*NB)

                   CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR )

                   ! Generate left bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)

                   CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, &
                                WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in WORK(IR)
                   ! (Workspace: need 2*N*N+4*N-1,
                   ! prefer 2*N*N+3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of R in WORK(IU) and computing
                   ! right singular vectors of R in WORK(IR)
                   ! (Workspace: need 2*N*N+BDSPAC)

                   CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), WORK( IR ), LDWRKR, &
                                WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply Q in A by left singular vectors of R in
                   ! WORK(IU), storing result in U
                   ! (Workspace: need N*N)

                   CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, WORK( IU ), &
                               LDWRKU, ZERO, U, LDU )

                   ! Copy right singular vectors of R to A
                   ! (Workspace: need N*N)

                   CALL DLACPY( 'F', N, N, WORK( IR ), LDWRKR, A, LDA )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DORGQR( M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Zero out below R in A

                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )

                   ! Bidiagonalize R in A
                   ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                   CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply Q in U by left vectors bidiagonalizing R
                   ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                   CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, &
                                LDU, WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right vectors bidiagonalizing R in A
                   ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U and computing right
                   ! singular vectors of A in A
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), A, &
                                LDA, U, LDU, DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTVAS ) THEN

                ! Path 6 (M much larger than N, JOBU='S', JOBVT='S'
                ! or 'A')
                ! N left singular vectors to be computed in U and
                ! N right singular vectors to be computed in VT

                IF ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+LDA*N ) THEN

                      ! WORK(IU) is LDA by N

                      LDWRKU = LDA

                   ELSE

                      ! WORK(IU) is N by N

                      LDWRKU = N

                   END IF

                   ITAU = IU + LDWRKU*N
                   IWORK = ITAU + N

                   ! Compute A=Q*R
                   ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy R to WORK(IU), zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU )

                   ! Generate Q in A
                   ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                   CALL DORGQR( M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in WORK(IU), copying result to VT
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                   CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, VT, LDVT )

                   ! Generate left bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                   CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in VT
                   ! (Workspace: need N*N+4*N-1,
                   ! prefer N*N+3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of R in WORK(IU) and computing
                   ! right singular vectors of R in VT
                   ! (Workspace: need N*N+BDSPAC)

                   CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), VT, LDVT, &
                                WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply Q in A by left singular vectors of R in
                   ! WORK(IU), storing result in U
                   ! (Workspace: need N*N)

                   CALL DGEMM( 'N', 'N', M, N, N, ONE, A, LDA, WORK( IU ), &
                               LDWRKU, ZERO, U, LDU )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DORGQR( M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy R to VT, zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                   IF ( N > 1 ) &
                        CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in VT
                   ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                   CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply Q in U by left bidiagonalizing vectors
                   ! in VT
                   ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                   CALL DORMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), &
                                U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in VT
                   ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U and computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), VT, &
                                LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO )

                END IF

             END IF

          ELSE IF ( WNTUA ) THEN

             IF ( WNTVN ) THEN

                ! Path 7 (M much larger than N, JOBU='A', JOBVT='N')
                ! M left singular vectors to be computed in U and
                ! no right singular vectors to be computed

                IF ( LWORK >= N*N+MAX( N+M, 4*N, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IR = 1

                   IF ( LWORK >= WRKBL+LDA*N ) THEN

                      ! WORK(IR) is LDA by N

                      LDWRKR = LDA

                   ELSE

                      ! WORK(IR) is N by N

                      LDWRKR = N

                   END IF

                   ITAU = IR + LDWRKR*N
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Copy R to WORK(IR), zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), LDWRKR )

                   ! Generate Q in U
                   ! (Workspace: need N*N+N+M, prefer N*N+N+M*NB)

                   CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in WORK(IR)
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                   CALL DGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in WORK(IR)
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                   CALL DORGBR( 'Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of R in WORK(IR)
                   ! (Workspace: need N*N+BDSPAC)

                   CALL DBDSQR( 'U', N, 0, N, 0, S, WORK( IE ), DUM, 1, WORK( IR ), &
                                LDWRKR, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply Q in U by left singular vectors of R in
                   ! WORK(IR), storing result in A
                   ! (Workspace: need N*N)

                   CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, WORK( IR ), &
                               LDWRKR, ZERO, A, LDA )

                   ! Copy left singular vectors of A from A to U

                   CALL DLACPY( 'F', M, N, A, LDA, U, LDU )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need N+M, prefer N+M*NB)

                   CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Zero out below R in A

                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )

                   ! Bidiagonalize R in A
                   ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                   CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply Q in U by left bidiagonalizing vectors
                   ! in A
                   ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                   CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', N, 0, M, 0, S, WORK( IE ), DUM, &
                                1, U, LDU, DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTVO ) THEN

                ! Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                ! M left singular vectors to be computed in U and
                ! N right singular vectors to be overwritten on A

                IF ( LWORK >= 2*N*N+MAX( N+M, 4*N, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+2*LDA*N ) THEN

                      ! WORK(IU) is LDA by N and WORK(IR) is LDA by N

                      LDWRKU = LDA
                      IR = IU + LDWRKU*N
                      LDWRKR = LDA

                   ELSE IF ( LWORK >= WRKBL+( LDA+N )*N ) THEN

                      ! WORK(IU) is LDA by N and WORK(IR) is N by N

                      LDWRKU = LDA
                      IR = IU + LDWRKU*N
                      LDWRKR = N

                   ELSE

                      ! WORK(IU) is N by N and WORK(IR) is N by N

                      LDWRKU = N
                      IR = IU + LDWRKU*N
                      LDWRKR = N

                   END IF

                   ITAU = IR + LDWRKR*N
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)

                   CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy R to WORK(IU), zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in WORK(IU), copying result to
                   ! WORK(IR)
                   ! (Workspace: need 2*N*N+4*N,
                   ! prefer 2*N*N+3*N+2*N*NB)

                   CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR )

                   ! Generate left bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)

                   CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in WORK(IR)
                   ! (Workspace: need 2*N*N+4*N-1,
                   ! prefer 2*N*N+3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of R in WORK(IU) and computing
                   ! right singular vectors of R in WORK(IR)
                   ! (Workspace: need 2*N*N+BDSPAC)

                   CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), WORK( IR ), LDWRKR, &
                                WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply Q in U by left singular vectors of R in
                   ! WORK(IU), storing result in A
                   ! (Workspace: need N*N)

                   CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, WORK( IU ), &
                               LDWRKU, ZERO, A, LDA )

                   ! Copy left singular vectors of A from A to U

                   CALL DLACPY( 'F', M, N, A, LDA, U, LDU )

                   ! Copy right singular vectors of R from WORK(IR) to A

                   CALL DLACPY( 'F', N, N, WORK( IR ), LDWRKR, A, LDA )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need N+M, prefer N+M*NB)

                   CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Zero out below R in A

                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )

                   ! Bidiagonalize R in A
                   ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                   CALL DGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply Q in U by left bidiagonalizing vectors
                   ! in A
                   ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                   CALL DORMBR( 'Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, &
                                LDU, WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in A
                   ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U and computing right
                   ! singular vectors of A in A
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), A, LDA, &
                                U, LDU, DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTVAS ) THEN

                ! Path 9 (M much larger than N, JOBU='A', JOBVT='S'
                ! or 'A')
                ! M left singular vectors to be computed in U and
                ! N right singular vectors to be computed in VT

                IF ( LWORK >= N*N+MAX( N+M, 4*N, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+LDA*N ) THEN

                      ! WORK(IU) is LDA by N

                      LDWRKU = LDA

                   ELSE

                      ! WORK(IU) is N by N

                      LDWRKU = N

                   END IF

                   ITAU = IU + LDWRKU*N
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need N*N+N+M, prefer N*N+N+M*NB)

                   CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy R to WORK(IU), zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in WORK(IU), copying result to VT
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                   CALL DGEBRD( N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', N, N, WORK( IU ), LDWRKU, VT, LDVT )

                   ! Generate left bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                   CALL DORGBR( 'Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in VT
                   ! (Workspace: need N*N+4*N-1,
                   ! prefer N*N+3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of R in WORK(IU) and computing
                   ! right singular vectors of R in VT
                   ! (Workspace: need N*N+BDSPAC)

                   CALL DBDSQR( 'U', N, N, N, 0, S, WORK( IE ), VT, LDVT, WORK( IU ), &
                                LDWRKU, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply Q in U by left singular vectors of R in
                   ! WORK(IU), storing result in A
                   ! (Workspace: need N*N)

                   CALL DGEMM( 'N', 'N', M, N, N, ONE, U, LDU, WORK( IU ), &
                               LDWRKU, ZERO, A, LDA )

                   ! Copy left singular vectors of A from A to U

                   CALL DLACPY( 'F', M, N, A, LDA, U, LDU )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + N

                   ! Compute A=Q*R, copying result to U
                   ! (Workspace: need 2*N, prefer N+N*NB)

                   CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, N, A, LDA, U, LDU )

                   ! Generate Q in U
                   ! (Workspace: need N+M, prefer N+M*NB)

                   CALL DORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy R from A to VT, zeroing out below it

                   CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
                   IF ( N > 1 ) &
                        CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT )
                   IE = ITAU
                   ITAUQ = IE + N
                   ITAUP = ITAUQ + N
                   IWORK = ITAUP + N

                   ! Bidiagonalize R in VT
                   ! (Workspace: need 4*N, prefer 3*N+2*N*NB)

                   CALL DGEBRD( N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply Q in U by left bidiagonalizing vectors
                   ! in VT
                   ! (Workspace: need 3*N+M, prefer 3*N+M*NB)

                   CALL DORMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), U, LDU, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in VT
                   ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                   CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + N

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U and computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', N, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, &
                                DUM, 1, WORK( IWORK ), INFO )

                END IF

             END IF

          END IF

       ELSE

          ! M < MNTHR

          ! Path 10 (M at least N, but not much larger)
          ! Reduce to bidiagonal form without QR decomposition

          IE = 1
          ITAUQ = IE + N
          ITAUP = ITAUQ + N
          IWORK = ITAUP + N

          ! Bidiagonalize A
          ! (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)

          CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                       WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

          IF ( WNTUAS ) THEN

             ! If left singular vectors desired in U, copy result to U
             ! and generate left bidiagonalizing vectors in U
             ! (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB)

             CALL DLACPY( 'L', M, N, A, LDA, U, LDU )
             IF ( WNTUS ) NCU = N
             IF ( WNTUA ) NCU = M
             CALL DORGBR( 'Q', M, NCU, N, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )
          END IF

          IF ( WNTVAS ) THEN

             ! If right singular vectors desired in VT, copy result to
             ! VT and generate right bidiagonalizing vectors in VT
             ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

             CALL DLACPY( 'U', N, N, A, LDA, VT, LDVT )
             CALL DORGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )
          END IF

          IF ( WNTUO ) THEN

             ! If left singular vectors desired in A, generate left
             ! bidiagonalizing vectors in A
             ! (Workspace: need 4*N, prefer 3*N+N*NB)

             CALL DORGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )
          END IF

          IF ( WNTVO ) THEN

             ! If right singular vectors desired in A, generate right
             ! bidiagonalizing vectors in A
             ! (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

             CALL DORGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )
          END IF

          IWORK = IE + N
          IF ( WNTUAS .OR. WNTUO ) NRU = M
          IF ( WNTUN ) NRU = 0
          IF ( WNTVAS .OR. WNTVO ) NCVT = N
          IF ( WNTVN ) NCVT = 0

          IF ( ( .NOT.WNTUO ) .AND. ( .NOT.WNTVO ) ) THEN

             ! Perform bidiagonal QR iteration, if desired, computing
             ! left singular vectors in U and computing right singular
             ! vectors in VT
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'U', N, NCVT, NRU, 0, S, WORK( IE ), VT, &
                          LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO )

          ELSE IF ( ( .NOT.WNTUO ) .AND. WNTVO ) THEN

             ! Perform bidiagonal QR iteration, if desired, computing
             ! left singular vectors in U and computing right singular
             ! vectors in A
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'U', N, NCVT, NRU, 0, S, WORK( IE ), A, LDA, &
                          U, LDU, DUM, 1, WORK( IWORK ), INFO )
          ELSE

             ! Perform bidiagonal QR iteration, if desired, computing
             ! left singular vectors in A and computing right singular
             ! vectors in VT
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'U', N, NCVT, NRU, 0, S, WORK( IE ), VT, &
                          LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO )
          END IF

       END IF

    ELSE

       ! A has more columns than rows. If A has sufficiently more
       ! columns than rows, first reduce using the LQ decomposition (if
       ! sufficient workspace available)

       IF ( N >= MNTHR ) THEN

          IF ( WNTVN ) THEN

             ! Path 1t(N much larger than M, JOBVT='N')
             ! No right singular vectors to be computed

             ITAU = 1
             IWORK = ITAU + M

             ! Compute A=L*Q
             ! (Workspace: need 2*M, prefer M+M*NB)

             CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR )

             ! Zero out above L

             CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )
             IE = 1
             ITAUQ = IE + M
             ITAUP = ITAUQ + M
             IWORK = ITAUP + M

             ! Bidiagonalize L in A
             ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

             CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                          WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

             IF ( WNTUO .OR. WNTUAS ) THEN

                ! If left singular vectors desired, generate Q
                ! (Workspace: need 4*M, prefer 3*M+M*NB)

                CALL DORGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
             END IF

             IWORK = IE + M
             NRU = 0
             IF ( WNTUO .OR. WNTUAS ) NRU = M

             ! Perform bidiagonal QR iteration, computing left singular
             ! vectors of A in A if desired
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'U', M, 0, NRU, 0, S, WORK( IE ), DUM, 1, A, LDA, DUM, &
                          1, WORK( IWORK ), INFO )

             ! If left singular vectors desired in U, copy them there

             IF ( WNTUAS ) CALL DLACPY( 'F', M, M, A, LDA, U, LDU )

          ELSE IF ( WNTVO .AND. WNTUN ) THEN

             ! Path 2t(N much larger than M, JOBU='N', JOBVT='O')
             ! M right singular vectors to be overwritten on A and
             ! no left singular vectors to be computed

             IF ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) THEN

                ! Sufficient workspace for a fast algorithm

                IR = 1

                IF ( LWORK >= MAX( WRKBL, LDA*N+M )+LDA*M ) THEN

                   ! WORK(IU) is LDA by N and WORK(IR) is LDA by M

                   LDWRKU = LDA
                   CHUNK = N
                   LDWRKR = LDA

                ELSE IF ( LWORK >= MAX( WRKBL, LDA*N+M )+M*M ) THEN

                   ! WORK(IU) is LDA by N and WORK(IR) is M by M

                   LDWRKU = LDA
                   CHUNK = N
                   LDWRKR = M

                ELSE

                   ! WORK(IU) is M by CHUNK and WORK(IR) is M by M

                   LDWRKU = M
                   CHUNK = ( LWORK-M*M-M ) / M
                   LDWRKR = M

                END IF

                ITAU = IR + LDWRKR*M
                IWORK = ITAU + M

                ! Compute A=L*Q
                ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                            LWORK-IWORK+1, IERR )

                ! Copy L to WORK(IR) and zero out above it

                CALL DLACPY( 'L', M, M, A, LDA, WORK( IR ), LDWRKR )
                CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IR+LDWRKR ), LDWRKR )

                ! Generate Q in A
                ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IE = ITAU
                ITAUQ = IE + M
                ITAUP = ITAUQ + M
                IWORK = ITAUP + M

                ! Bidiagonalize L in WORK(IR)
                ! (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                CALL DGEBRD( M, M, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), &
                             WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate right vectors bidiagonalizing L
                ! (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)

                CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )
                IWORK = IE + M

                ! Perform bidiagonal QR iteration, computing right
                ! singular vectors of L in WORK(IR)
                ! (Workspace: need M*M+BDSPAC)

                CALL DBDSQR( 'U', M, M, 0, 0, S, WORK( IE ), WORK( IR ), LDWRKR, DUM, &
                             1, DUM, 1, WORK( IWORK ), INFO )
                IU = IE + M

                ! Multiply right singular vectors of L in WORK(IR) by Q
                ! in A, storing result in WORK(IU) and copying to A
                ! (Workspace: need M*M+2*M, prefer M*M+M*N+M)

                DO I = 1, N, CHUNK
                   BLK = MIN( N-I+1, CHUNK )
                   CALL DGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IR ), LDWRKR, A( 1, I ), &
                               LDA, ZERO, WORK( IU ), LDWRKU )
                   CALL DLACPY( 'F', M, BLK, WORK( IU ), LDWRKU, A( 1, I ), LDA )
                END DO

             ELSE

                ! Insufficient workspace for a fast algorithm

                IE = 1
                ITAUQ = IE + M
                ITAUP = ITAUQ + M
                IWORK = ITAUP + M

                ! Bidiagonalize A
                ! (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

                CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate right vectors bidiagonalizing A
                ! (Workspace: need 4*M, prefer 3*M+M*NB)

                CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IWORK = IE + M

                ! Perform bidiagonal QR iteration, computing right
                ! singular vectors of A in A
                ! (Workspace: need BDSPAC)

                CALL DBDSQR( 'L', M, N, 0, 0, S, WORK( IE ), A, LDA, DUM, 1, DUM, &
                             1, WORK( IWORK ), INFO )

             END IF

          ELSE IF ( WNTVO .AND. WNTUAS ) THEN

             ! Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
             ! M right singular vectors to be overwritten on A and
             ! M left singular vectors to be computed in U

             IF ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) THEN

                ! Sufficient workspace for a fast algorithm

                IR = 1

                IF ( LWORK >= MAX( WRKBL, LDA*N+M )+LDA*M ) THEN

                   ! WORK(IU) is LDA by N and WORK(IR) is LDA by M

                   LDWRKU = LDA
                   CHUNK = N
                   LDWRKR = LDA

                ELSE IF ( LWORK >= MAX( WRKBL, LDA*N+M )+M*M ) THEN

                   ! WORK(IU) is LDA by N and WORK(IR) is M by M

                   LDWRKU = LDA
                   CHUNK = N
                   LDWRKR = M

                ELSE

                   ! WORK(IU) is M by CHUNK and WORK(IR) is M by M

                   LDWRKU = M
                   CHUNK = ( LWORK-M*M-M ) / M
                   LDWRKR = M

                END IF

                ITAU = IR + LDWRKR*M
                IWORK = ITAU + M

                ! Compute A=L*Q
                ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )

                ! Copy L to U, zeroing about above it

                CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU )

                ! Generate Q in A
                ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IE = ITAU
                ITAUQ = IE + M
                ITAUP = ITAUQ + M
                IWORK = ITAUP + M

                ! Bidiagonalize L in U, copying result to WORK(IR)
                ! (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )
                CALL DLACPY( 'U', M, M, U, LDU, WORK( IR ), LDWRKR )

                ! Generate right vectors bidiagonalizing L in WORK(IR)
                ! (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)

                CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate left vectors bidiagonalizing L in U
                ! (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)

                CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IWORK = IE + M

                ! Perform bidiagonal QR iteration, computing left
                ! singular vectors of L in U, and computing right
                ! singular vectors of L in WORK(IR)
                ! (Workspace: need M*M+BDSPAC)

                CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), WORK( IR ), LDWRKR, U, &
                             LDU, DUM, 1, WORK( IWORK ), INFO )
                IU = IE + M

                ! Multiply right singular vectors of L in WORK(IR) by Q
                ! in A, storing result in WORK(IU) and copying to A
                ! (Workspace: need M*M+2*M, prefer M*M+M*N+M))

                DO I = 1, N, CHUNK
                   BLK = MIN( N-I+1, CHUNK )
                   CALL DGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IR ), LDWRKR, A( 1, I ), &
                               LDA, ZERO, WORK( IU ), LDWRKU )
                   CALL DLACPY( 'F', M, BLK, WORK( IU ), LDWRKU, A( 1, I ), LDA )
                END DO

             ELSE

                ! Insufficient workspace for a fast algorithm

                ITAU = 1
                IWORK = ITAU + M

                ! Compute A=L*Q
                ! (Workspace: need 2*M, prefer M+M*NB)

                CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )

                ! Copy L to U, zeroing out above it

                CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU )

                ! Generate Q in A
                ! (Workspace: need 2*M, prefer M+M*NB)

                CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IE = ITAU
                ITAUQ = IE + M
                ITAUP = ITAUQ + M
                IWORK = ITAUP + M

                ! Bidiagonalize L in U
                ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Multiply right vectors bidiagonalizing L by Q in A
                ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                CALL DORMBR( 'P', 'L', 'T', M, N, M, U, LDU, WORK( ITAUP ), A, LDA, &
                             WORK( IWORK ), LWORK-IWORK+1, IERR )

                ! Generate left vectors bidiagonalizing L in U
                ! (Workspace: need 4*M, prefer 3*M+M*NB)

                CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                             LWORK-IWORK+1, IERR )
                IWORK = IE + M

                ! Perform bidiagonal QR iteration, computing left
                ! singular vectors of A in U and computing right
                ! singular vectors of A in A
                ! (Workspace: need BDSPAC)

                CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), A, LDA, U, LDU, DUM, &
                             1, WORK( IWORK ), INFO )

             END IF

          ELSE IF ( WNTVS ) THEN

             IF ( WNTUN ) THEN

                ! Path 4t(N much larger than M, JOBU='N', JOBVT='S')
                ! M right singular vectors to be computed in VT and
                ! no left singular vectors to be computed

                IF ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IR = 1

                   IF ( LWORK >= WRKBL+LDA*M ) THEN

                      ! WORK(IR) is LDA by M

                      LDWRKR = LDA

                   ELSE

                      ! WORK(IR) is M by M

                      LDWRKR = M

                   END IF

                   ITAU = IR + LDWRKR*M
                   IWORK = ITAU + M

                   ! Compute A=L*Q
                   ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to WORK(IR), zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, WORK( IR ), LDWRKR )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IR+LDWRKR ), LDWRKR )

                   ! Generate Q in A
                   ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                   CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in WORK(IR)
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                   CALL DGEBRD( M, M, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right vectors bidiagonalizing L in
                   ! WORK(IR)
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)

                   CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing right
                   ! singular vectors of L in WORK(IR)
                   ! (Workspace: need M*M+BDSPAC)

                   CALL DBDSQR( 'U', M, M, 0, 0, S, WORK( IE ), WORK( IR ), LDWRKR, &
                                DUM, 1, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply right singular vectors of L in WORK(IR) by
                   ! Q in A, storing result in VT
                   ! (Workspace: need M*M)

                   CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IR ), LDWRKR, A, LDA, &
                               ZERO, VT, LDVT )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + M

                   ! Compute A=L*Q
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy result to VT

                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DORGLQ( M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Zero out above L in A

                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )

                   ! Bidiagonalize L in A
                   ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                   CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply right vectors bidiagonalizing L by Q in VT
                   ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                   CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', M, N, 0, 0, S, WORK( IE ), VT, LDVT, DUM, 1, &
                                DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTUO ) THEN

                ! Path 5t(N much larger than M, JOBU='O', JOBVT='S')
                ! M right singular vectors to be computed in VT and
                ! M left singular vectors to be overwritten on A

                IF ( LWORK >= 2*M*M+MAX( 4*M, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+2*LDA*M ) THEN

                      ! WORK(IU) is LDA by M and WORK(IR) is LDA by M

                      LDWRKU = LDA
                      IR = IU + LDWRKU*M
                      LDWRKR = LDA

                   ELSE IF ( LWORK >= WRKBL+( LDA+M )*M ) THEN

                      ! WORK(IU) is LDA by M and WORK(IR) is M by M

                      LDWRKU = LDA
                      IR = IU + LDWRKU*M
                      LDWRKR = M

                   ELSE

                      ! WORK(IU) is M by M and WORK(IR) is M by M

                      LDWRKU = M
                      IR = IU + LDWRKU*M
                      LDWRKR = M

                   END IF

                   ITAU = IR + LDWRKR*M
                   IWORK = ITAU + M

                   ! Compute A=L*Q
                   ! (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to WORK(IU), zeroing out below it

                   CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU )

                   ! Generate Q in A
                   ! (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)

                   CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in WORK(IU), copying result to
                   ! WORK(IR)
                   ! (Workspace: need 2*M*M+4*M,
                   ! prefer 2*M*M+3*M+2*M*NB)

                   CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR )

                   ! Generate right bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need 2*M*M+4*M-1,
                   ! prefer 2*M*M+3*M+(M-1)*NB)

                   CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in WORK(IR)
                   ! (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of L in WORK(IR) and computing
                   ! right singular vectors of L in WORK(IU)
                   ! (Workspace: need 2*M*M+BDSPAC)

                   CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, &
                                WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply right singular vectors of L in WORK(IU) by
                   ! Q in A, storing result in VT
                   ! (Workspace: need M*M)

                   CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, A, LDA, &
                               ZERO, VT, LDVT )

                   ! Copy left singular vectors of L to A
                   ! (Workspace: need M*M)

                   CALL DLACPY( 'F', M, M, WORK( IR ), LDWRKR, A, LDA )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DORGLQ( M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Zero out above L in A

                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )

                   ! Bidiagonalize L in A
                   ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                   CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply right vectors bidiagonalizing L by Q in VT
                   ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                   CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors of L in A
                   ! (Workspace: need 4*M, prefer 3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, compute left
                   ! singular vectors of A in A and compute right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, LDVT, A, LDA, DUM, &
                                1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTUAS ) THEN

                ! Path 6t(N much larger than M, JOBU='S' or 'A',
                ! JOBVT='S')
                ! M right singular vectors to be computed in VT and
                ! M left singular vectors to be computed in U

                IF ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+LDA*M ) THEN

                      ! WORK(IU) is LDA by N

                      LDWRKU = LDA

                   ELSE

                      ! WORK(IU) is LDA by M

                      LDWRKU = M

                   END IF

                   ITAU = IU + LDWRKU*M
                   IWORK = ITAU + M

                   ! Compute A=L*Q
                   ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to WORK(IU), zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU )

                   ! Generate Q in A
                   ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                   CALL DORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in WORK(IU), copying result to U
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                   CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, U, LDU )

                   ! Generate right bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need M*M+4*M-1,
                   ! prefer M*M+3*M+(M-1)*NB)

                   CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in U
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of L in U and computing right
                   ! singular vectors of L in WORK(IU)
                   ! (Workspace: need M*M+BDSPAC)

                   CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, U, &
                                LDU, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply right singular vectors of L in WORK(IU) by
                   ! Q in A, storing result in VT
                   ! (Workspace: need M*M)

                   CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, A, LDA, &
                               ZERO, VT, LDVT )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DORGLQ( M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to U, zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU )

                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in U
                   ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                   CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply right bidiagonalizing vectors in U by Q
                   ! in VT
                   ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                   CALL DORMBR( 'P', 'L', 'T', M, N, M, U, LDU, WORK( ITAUP ), VT, LDVT, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in U
                   ! (Workspace: need 4*M, prefer 3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U and computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, &
                                1, WORK( IWORK ), INFO )

                END IF

             END IF

          ELSE IF ( WNTVA ) THEN

             IF ( WNTUN ) THEN

                ! Path 7t(N much larger than M, JOBU='N', JOBVT='A')
                ! N right singular vectors to be computed in VT and
                ! no left singular vectors to be computed

                IF ( LWORK >= M*M+MAX( N+M, 4*M, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IR = 1

                   IF ( LWORK >= WRKBL+LDA*M ) THEN

                      ! WORK(IR) is LDA by M

                      LDWRKR = LDA

                   ELSE

                      ! WORK(IR) is M by M

                      LDWRKR = M

                   END IF

                   ITAU = IR + LDWRKR*M
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Copy L to WORK(IR), zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, WORK( IR ), LDWRKR )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IR+LDWRKR ), LDWRKR )

                   ! Generate Q in VT
                   ! (Workspace: need M*M+M+N, prefer M*M+M+N*NB)

                   CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in WORK(IR)
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                   CALL DGEBRD( M, M, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate right bidiagonalizing vectors in WORK(IR)
                   ! (Workspace: need M*M+4*M-1,
                   ! prefer M*M+3*M+(M-1)*NB)

                   CALL DORGBR( 'P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing right
                   ! singular vectors of L in WORK(IR)
                   ! (Workspace: need M*M+BDSPAC)

                   CALL DBDSQR( 'U', M, M, 0, 0, S, WORK( IE ), WORK( IR ), LDWRKR, &
                                DUM, 1, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply right singular vectors of L in WORK(IR) by
                   ! Q in VT, storing result in A
                   ! (Workspace: need M*M)

                   CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IR ), LDWRKR, VT, LDVT, &
                               ZERO, A, LDA )

                   ! Copy right singular vectors of A from A to VT

                   CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                               LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need M+N, prefer M+N*NB)

                   CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Zero out above L in A

                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )

                   ! Bidiagonalize L in A
                   ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                   CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), &
                        WORK( ITAUQ ), WORK( ITAUP ), &
                        WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply right bidiagonalizing vectors in A by Q
                   ! in VT
                   ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                   CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', M, N, 0, 0, S, WORK( IE ), VT, LDVT, DUM, 1, &
                                DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTUO ) THEN

                ! Path 8t(N much larger than M, JOBU='O', JOBVT='A')
                ! N right singular vectors to be computed in VT and
                ! M left singular vectors to be overwritten on A

                IF ( LWORK >= 2*M*M+MAX( N+M, 4*M, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+2*LDA*M ) THEN

                      ! WORK(IU) is LDA by M and WORK(IR) is LDA by M

                      LDWRKU = LDA
                      IR = IU + LDWRKU*M
                      LDWRKR = LDA

                   ELSE IF ( LWORK >= WRKBL+( LDA+M )*M ) THEN

                      ! WORK(IU) is LDA by M and WORK(IR) is M by M

                      LDWRKU = LDA
                      IR = IU + LDWRKU*M
                      LDWRKR = M

                   ELSE

                      ! WORK(IU) is M by M and WORK(IR) is M by M

                      LDWRKU = M
                      IR = IU + LDWRKU*M
                      LDWRKR = M

                   END IF

                   ITAU = IR + LDWRKR*M
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)

                   CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to WORK(IU), zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU )

                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in WORK(IU), copying result to
                   ! WORK(IR)
                   ! (Workspace: need 2*M*M+4*M,
                   ! prefer 2*M*M+3*M+2*M*NB)

                   CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR )

                   ! Generate right bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need 2*M*M+4*M-1,
                   ! prefer 2*M*M+3*M+(M-1)*NB)

                   CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in WORK(IR)
                   ! (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUQ ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of L in WORK(IR) and computing
                   ! right singular vectors of L in WORK(IU)
                   ! (Workspace: need 2*M*M+BDSPAC)

                   CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, &
                                WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply right singular vectors of L in WORK(IU) by
                   ! Q in VT, storing result in A
                   ! (Workspace: need M*M)

                   CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, VT, &
                               LDVT, ZERO, A, LDA )

                   ! Copy right singular vectors of A from A to VT

                   CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )

                   ! Copy left singular vectors of A from WORK(IR) to A

                   CALL DLACPY( 'F', M, M, WORK( IR ), LDWRKR, A, LDA )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need M+N, prefer M+N*NB)

                   CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Zero out above L in A

                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )

                   ! Bidiagonalize L in A
                   ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                   CALL DGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply right bidiagonalizing vectors in A by Q
                   ! in VT
                   ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                   CALL DORMBR( 'P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, &
                                LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in A
                   ! (Workspace: need 4*M, prefer 3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in A and computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, LDVT, A, LDA, &
                                DUM, 1, WORK( IWORK ), INFO )

                END IF

             ELSE IF ( WNTUAS ) THEN

                ! Path 9t(N much larger than M, JOBU='S' or 'A',
                ! JOBVT='A')
                ! N right singular vectors to be computed in VT and
                ! M left singular vectors to be computed in U

                IF ( LWORK >= M*M+MAX( N+M, 4*M, BDSPAC ) ) THEN

                   ! Sufficient workspace for a fast algorithm

                   IU = 1

                   IF ( LWORK >= WRKBL+LDA*M ) THEN

                      ! WORK(IU) is LDA by M

                      LDWRKU = LDA

                   ELSE

                      ! WORK(IU) is M by M

                      LDWRKU = M

                   END IF

                   ITAU = IU + LDWRKU*M
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need M*M+M+N, prefer M*M+M+N*NB)

                   CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to WORK(IU), zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, WORK( IU ), LDWRKU )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU )

                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in WORK(IU), copying result to U
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                   CALL DGEBRD( M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'L', M, M, WORK( IU ), LDWRKU, U, LDU )

                   ! Generate right bidiagonalizing vectors in WORK(IU)
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)

                   CALL DORGBR( 'P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), &
                                WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in U
                   ! (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of L in U and computing right
                   ! singular vectors of L in WORK(IU)
                   ! (Workspace: need M*M+BDSPAC)

                   CALL DBDSQR( 'U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, &
                                U, LDU, DUM, 1, WORK( IWORK ), INFO )

                   ! Multiply right singular vectors of L in WORK(IU) by
                   ! Q in VT, storing result in A
                   ! (Workspace: need M*M)

                   CALL DGEMM( 'N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, VT, &
                               LDVT, ZERO, A, LDA )

                   ! Copy right singular vectors of A from A to VT

                   CALL DLACPY( 'F', M, N, A, LDA, VT, LDVT )

                ELSE

                   ! Insufficient workspace for a fast algorithm

                   ITAU = 1
                   IWORK = ITAU + M

                   ! Compute A=L*Q, copying result to VT
                   ! (Workspace: need 2*M, prefer M+M*NB)

                   CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )

                   ! Generate Q in VT
                   ! (Workspace: need M+N, prefer M+N*NB)

                   CALL DORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )

                   ! Copy L to U, zeroing out above it

                   CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
                   CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU )

                   IE = ITAU
                   ITAUQ = IE + M
                   ITAUP = ITAUQ + M
                   IWORK = ITAUP + M

                   ! Bidiagonalize L in U
                   ! (Workspace: need 4*M, prefer 3*M+2*M*NB)

                   CALL DGEBRD( M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), &
                                WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Multiply right bidiagonalizing vectors in U by Q
                   ! in VT
                   ! (Workspace: need 3*M+N, prefer 3*M+N*NB)

                   CALL DORMBR( 'P', 'L', 'T', M, N, M, U, LDU, WORK( ITAUP ), VT, &
                                LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR )

                   ! Generate left bidiagonalizing vectors in U
                   ! (Workspace: need 4*M, prefer 3*M+M*NB)

                   CALL DORGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                                LWORK-IWORK+1, IERR )
                   IWORK = IE + M

                   ! Perform bidiagonal QR iteration, computing left
                   ! singular vectors of A in U and computing right
                   ! singular vectors of A in VT
                   ! (Workspace: need BDSPAC)

                   CALL DBDSQR( 'U', M, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, &
                                DUM, 1, WORK( IWORK ), INFO )

                END IF

             END IF

          END IF

       ELSE

          ! N < MNTHR

          ! Path 10t(N greater than M, but not much larger)
          ! Reduce to bidiagonal form without LQ decomposition

          IE = 1
          ITAUQ = IE + M
          ITAUP = ITAUQ + M
          IWORK = ITAUP + M

          ! Bidiagonalize A
          ! (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

          CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), &
                       WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR )

          IF ( WNTUAS ) THEN

             ! If left singular vectors desired in U, copy result to U
             ! and generate left bidiagonalizing vectors in U
             ! (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)

             CALL DLACPY( 'L', M, M, A, LDA, U, LDU )
             CALL DORGBR( 'Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )

          END IF

          IF ( WNTVAS ) THEN

             ! If right singular vectors desired in VT, copy result to
             ! VT and generate right bidiagonalizing vectors in VT
             ! (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB)

             CALL DLACPY( 'U', M, N, A, LDA, VT, LDVT )
             IF ( WNTVA ) NRVT = N
             IF ( WNTVS ) NRVT = M
             CALL DORGBR( 'P', NRVT, N, M, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )
          END IF

          IF ( WNTUO ) THEN

             ! If left singular vectors desired in A, generate left
             ! bidiagonalizing vectors in A
             ! (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)

             CALL DORGBR( 'Q', M, M, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), &
                  LWORK-IWORK+1, IERR )

          END IF

          IF ( WNTVO ) THEN

             ! If right singular vectors desired in A, generate right
             ! bidiagonalizing vectors in A
             ! (Workspace: need 4*M, prefer 3*M+M*NB)

             CALL DORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), &
                          LWORK-IWORK+1, IERR )
          END IF

          IWORK = IE + M
          IF ( WNTUAS .OR. WNTUO ) NRU = M
          IF ( WNTUN ) NRU = 0
          IF ( WNTVAS .OR. WNTVO ) NCVT = N
          IF ( WNTVN ) NCVT = 0

          IF ( ( .NOT.WNTUO ) .AND. ( .NOT.WNTVO ) ) THEN

             ! Perform bidiagonal QR iteration, if desired, computing
             ! left singular vectors in U and computing right singular
             ! vectors in VT
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'L', M, NCVT, NRU, 0, S, WORK( IE ), VT, LDVT, U, LDU, &
                          DUM, 1, WORK( IWORK ), INFO )

          ELSE IF ( ( .NOT.WNTUO ) .AND. WNTVO ) THEN

             ! Perform bidiagonal QR iteration, if desired, computing
             ! left singular vectors in U and computing right singular
             ! vectors in A
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'L', M, NCVT, NRU, 0, S, WORK( IE ), A, LDA, U, LDU, &
                          DUM, 1, WORK( IWORK ), INFO )

          ELSE

             ! Perform bidiagonal QR iteration, if desired, computing
             ! left singular vectors in A and computing right singular
             ! vectors in VT
             ! (Workspace: need BDSPAC)

             CALL DBDSQR( 'L', M, NCVT, NRU, 0, S, WORK( IE ), VT, LDVT, A, LDA, &
                          DUM, 1, WORK( IWORK ), INFO )

          END IF

       END IF

    END IF

    ! If DBDSQR failed to converge, copy unconverged superdiagonals
    ! to WORK( 2:MINMN )

    IF ( INFO /= 0 ) THEN
       IF ( IE > 2 ) THEN
          DO I = 1, MINMN - 1
             WORK( I+1 ) = WORK( I+IE-1 )
          END DO
       END IF
       IF ( IE < 2 ) THEN
          DO I = MINMN - 1, 1, -1
             WORK( I+1 ) = WORK( I+IE-1 )
          END DO
       END IF
    END IF

    ! Undo scaling if necessary

    IF ( ISCL == 1 ) THEN
       IF ( ANRM > BIGNUM ) &
            CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR )
       IF ( INFO /= 0 .AND. ANRM > BIGNUM ) &
            CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN-1, 1, WORK( 2 ), MINMN, IERR )
       IF ( ANRM < SMLNUM ) &
            CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR )
       IF ( INFO /= 0 .AND. ANRM < SMLNUM ) &
            CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN-1, 1, WORK( 2 ), MINMN, IERR )
    END IF

    ! Return optimal workspace in WORK(1)

    WORK( 1 ) = REAL( MAXWRK, DP )

  END SUBROUTINE DGESVD

!==========================================================================
!
!  DBDSQR
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DBDSQR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dbdsqr.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
!                          LDU, C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
!       ..
!       .. Array Arguments ..
!       REAL(DP)           C( LDC,* ), D( * ), E( * ), U( LDU,* ),
!                          VT( LDVT,* ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DBDSQR computes the singular values and, optionally, the right and/or
! left singular vectors from the singular value decomposition (SVD) of
! a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
! zero-shift QR algorithm.  The SVD of B has the form
!
!    B = Q * S * P**T
!
! where S is the diagonal matrix of singular values, Q is an orthogonal
! matrix of left singular vectors, and P is an orthogonal matrix of
! right singular vectors.  If left singular vectors are requested, this
! subroutine actually returns U*Q instead of Q, and, if right singular
! vectors are requested, this subroutine returns P**T*VT instead of
! P**T, for given real input matrices U and VT.  When U and VT are the
! orthogonal matrices that reduce a general matrix A to bidiagonal
! form:  A = U*B*VT, as computed by DGEBRD, then
!
!    A = (U*Q) * S * (P**T*VT)
!
! is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
! for a given real input matrix C.
!
! See "Computing  Small Singular Values of Bidiagonal Matrices With
! Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
! LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
! no. 5, pp. 873-912, Sept 1990) and
! "Accurate singular values and differential qd algorithms," by
! B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
! Department, University of California at Berkeley, July 1992
! for a detailed description of the algorithm.
!
!  Arguments:
!  ==========
!
! param[in] UPLO
!
!          UPLO is CHARACTER*1
!          = 'U':  B is upper bidiagonal;
!          = 'L':  B is lower bidiagonal.
!
! param[in] N
!
!          N is INTEGER
!          The order of the matrix B.  N >= 0.
!
! param[in] NCVT
!
!          NCVT is INTEGER
!          The number of columns of the matrix VT. NCVT >= 0.
!
! param[in] NRU
!
!          NRU is INTEGER
!          The number of rows of the matrix U. NRU >= 0.
!
! param[in] NCC
!
!          NCC is INTEGER
!          The number of columns of the matrix C. NCC >= 0.
!
! param[in,out] D
!
!          D is REAL(DP)         array, dimension (N)
!          On entry, the n diagonal elements of the bidiagonal matrix B.
!          On exit, if INFO=0, the singular values of B in decreasing
!          order.
!
! param[in,out] E
!
!          E is REAL(DP)         array, dimension (N-1)
!          On entry, the N-1 offdiagonal elements of the bidiagonal
!          matrix B.
!          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
!          will contain the diagonal and superdiagonal elements of a
!          bidiagonal matrix orthogonally equivalent to the one given
!          as input.
!
! param[in,out] VT
!
!          VT is REAL(DP)         array, dimension (LDVT, NCVT)
!          On entry, an N-by-NCVT matrix VT.
!          On exit, VT is overwritten by P**T * VT.
!          Not referenced if NCVT = 0.
!
! param[in] LDVT
!
!          LDVT is INTEGER
!          The leading dimension of the array VT.
!          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
!
! param[in,out] U
!
!          U is REAL(DP)         array, dimension (LDU, N)
!          On entry, an NRU-by-N matrix U.
!          On exit, U is overwritten by U * Q.
!          Not referenced if NRU = 0.
!
! param[in] LDU
!
!          LDU is INTEGER
!          The leading dimension of the array U.  LDU >= max(1,NRU).
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC, NCC)
!          On entry, an N-by-NCC matrix C.
!          On exit, C is overwritten by Q**T * C.
!          Not referenced if NCC = 0.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C.
!          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (4*N)
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  If INFO = -i, the i-th argument had an illegal value
!          > 0:
!             if NCVT = NRU = NCC = 0,
!                = 1, a split was marked by a positive value in E
!                = 2, current block of Z not diagonalized after 30*N
!                     iterations (in inner while loop)
!                = 3, termination criterion of outer while loop not met
!                     (program created more than N unreduced blocks)
!             else NCVT = NRU = NCC = 0,
!                   the algorithm did not converge; D and E contain the
!                   elements of a bidiagonal matrix which is orthogonally
!                   similar to the input matrix B;  if INFO = i, i
!                   elements of E have not converged to zero.
!
!  Internal Parameters:
!  =========================
!
!  TOLMUL  REAL(DP)        , default = max(10,min(100,EPS**(-1/8)))
!          TOLMUL controls the convergence criterion of the QR loop.
!          If it is positive, TOLMUL*EPS is the desired relative
!             precision in the computed singular values.
!          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
!             desired absolute accuracy in the computed singular
!             values (corresponds to relative accuracy
!             abs(TOLMUL*EPS) in the largest singular value.
!          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
!             between 10 (for fast convergence) and .1/EPS
!             (for there to be some accuracy in the results).
!          Default is to lose at either one eighth or 2 of the
!             available decimal digits in each computed singular value
!             (whichever is smaller).
!
!  MAXITR  INTEGER, default = 6
!          MAXITR controls the maximum number of passes of the
!          algorithm through its inner loop. The algorithms stops
!          (and so fails to converge) if the number of passes
!          through the inner loop exceeds MAXITR*N**2.
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

  SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, &
                     LDU, C, LDC, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DLARTG, DLAS2, DLASQ1, DLASR, DLASV2, DROT, &
                               DSCAL,  DSWAP, XERBLA, LSAME, DLAMCH

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: UPLO
    INTEGER,   INTENT(IN)  :: LDC, LDU, LDVT, N, NCC, NCVT, NRU
    INTEGER,   INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: C( LDC,* ), D( * ), E( * ), U( LDU,* ), &
                               VT( LDVT,* ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO   = 0.0_DP
    REAL(DP), PARAMETER :: ONE    = 1.0_DP
    REAL(DP), PARAMETER :: NEGONE = -1.0_DP
    REAL(DP), PARAMETER :: HNDRTH = 0.01_DP
    REAL(DP), PARAMETER :: TEN    = 10.0_DP
    REAL(DP), PARAMETER :: HNDRD  = 100.0_DP
    REAL(DP), PARAMETER :: MEIGTH = -0.125_DP
    REAL(DP), PARAMETER :: ZEPS   = TINY( ZERO )
    INTEGER,  PARAMETER :: MAXITR = 6

    ! .. Local Scalars ..
    LOGICAL  :: LOWER, ROTATE
    INTEGER  :: I, IDIR, ISUB, ITER, J, LL, LLL, M, MAXIT, NM1, &
                NM12, NM13, OLDLL, OLDM
    REAL(DP) :: ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU, OLDCS, &
                OLDSN, R, SHIFT, SIGMN, SIGMX, SINL, SINR, SLL, SMAX, &
                SMIN, SMINL, SMINOA, SN, THRESH, TOL, TOLMUL, UNFL, DN

    ! .. Executable Statements ..

    ! Test the input parameters.

    INFO = 0
    LOWER = LSAME( UPLO, 'L' )

    IF ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( NCVT < 0 ) THEN
       INFO = -3
    ELSE IF ( NRU < 0 ) THEN
       INFO = -4
    ELSE IF ( NCC < 0 ) THEN
       INFO = -5
    ELSE IF ( ( NCVT == 0 .AND. LDVT < 1 ) .OR. &
              ( NCVT > 0 .AND. LDVT < MAX( 1, N ) ) ) THEN
       INFO = -9
    ELSE IF ( LDU < MAX( 1, NRU ) ) THEN
       INFO = -11
    ELSE IF ( ( NCC == 0 .AND. LDC < 1 ) .OR. &
              ( NCC > 0 .AND. LDC < MAX( 1, N ) ) ) THEN
       INFO = -13
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DBDSQR', -INFO )
       RETURN
    END IF

    IF ( N == 0 ) RETURN
    IF ( N > 1 ) THEN

       ! ROTATE is true if any singular vectors desired, false otherwise

       ROTATE = ( NCVT > 0 ) .OR. ( NRU > 0 ) .OR. ( NCC > 0 )

       ! If no singular vectors desired, use qd algorithm

       IF ( .NOT.ROTATE ) THEN

          CALL DLASQ1( N, D, E, WORK, INFO )

          ! If INFO equals 2, dqds didn't finish, try to finish
          IF ( INFO /= 2 ) RETURN
          INFO = 0

       END IF

       NM1 = N - 1
       NM12 = NM1 + NM1
       NM13 = NM12 + NM1
       IDIR = 0

       ! Get machine constants

       EPS = DLAMCH( 'Epsilon' )
       UNFL = DLAMCH( 'Safe minimum' )

       ! If matrix lower bidiagonal, rotate to be upper bidiagonal
       ! by applying Givens rotations on the left

       IF ( LOWER ) THEN

          DO I = 1, N - 1
             CALL DLARTG( D( I ), E( I ), CS, SN, R )
             D( I ) = R
             E( I ) = SN*D( I+1 )
             D( I+1 ) = CS*D( I+1 )
             WORK( I ) = CS
             WORK( NM1+I ) = SN
          END DO

          ! Update singular vectors if desired

          IF ( NRU > 0 ) &
               CALL DLASR( 'R', 'V', 'F', NRU, N, WORK( 1 ), WORK( N ), U, LDU )
          IF ( NCC > 0 ) &
               CALL DLASR( 'L', 'V', 'F', N, NCC, WORK( 1 ), WORK( N ), C, LDC )

       END IF

       ! Compute singular values to relative accuracy TOL
       ! (By setting TOL to be negative, algorithm will compute
       ! singular values to absolute accuracy ABS(TOL)*norm(input matrix))

       TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
       TOL = TOLMUL*EPS

       ! Compute approximate maximum, minimum singular values

       SMAX = ZERO
       DO I = 1, N
          SMAX = MAX( SMAX, ABS( D( I ) ) )
       END DO
       DO I = 1, N - 1
          SMAX = MAX( SMAX, ABS( E( I ) ) )
       END DO
       SMINL = ZERO

       DN = REAL( N,DP )
       IF ( TOL >= ZERO ) THEN

          ! Relative accuracy desired
          SMINOA = ABS( D( 1 ) )
          IF ( SMINOA >= ZEPS ) THEN
             MU = SMINOA
             DO I = 2, N
                MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
                SMINOA = MIN( SMINOA, MU )
                IF ( ABS(SMINOA) < ZEPS ) EXIT
             END DO
          END IF
          SMINOA = SMINOA / SQRT( DN )
          THRESH = MAX( TOL*SMINOA, REAL(MAXITR,DP)*DN*DN*UNFL )

       ELSE

          ! Absolute accuracy desired
          THRESH = MAX( ABS( TOL )*SMAX, REAL(MAXITR,DP)*DN*DN*UNFL )

       END IF

       ! Prepare for main iteration loop for the singular values
       ! (MAXIT is the maximum number of passes through the inner
       ! loop permitted before nonconvergence signalled.)

       MAXIT = MAXITR*N*N
       ITER = 0
       OLDLL = -1
       OLDM = -1

       ! M points to last element of unconverged part of matrix
       M = N

       ! Begin main iteration loop

       DO WHILE (.TRUE.)

          ! Check for convergence

          IF ( M <= 1 ) EXIT

          ! Check if Maximum number of iterations exceeded

          IF ( ITER > MAXIT ) THEN

             ! Failure to converge
             INFO = 0
             DO I = 1, N - 1
                IF ( ABS( E(I) ) >= ZEPS ) INFO = INFO + 1
             END DO

             RETURN
          ENDIF

          ! Find diagonal block of matrix to work on

          IF ( TOL < ZERO .AND. ABS( D( M ) ) <= THRESH ) D( M ) = ZERO

          SMAX = ABS( D( M ) )
          SMIN = SMAX

          DO LLL = 1, M - 1
             LL = M - LLL
             ABSS = ABS( D( LL ) )
             ABSE = ABS( E( LL ) )
             IF ( TOL < ZERO .AND. ABSS <= THRESH ) D( LL ) = ZERO
             IF ( ABSE <= THRESH ) EXIT
             SMIN = MIN( SMIN, ABSS )
             SMAX = MAX( SMAX, ABSS, ABSE )
          END DO

          IF ( LLL == M ) THEN
             LL = 0
          ELSE
             ! Matrix splits since E(LL) = 0
             E( LL ) = ZERO
             IF ( LL == M-1 ) THEN
                ! Convergence of bottom singular value, return to top of loop
                M = M - 1
                CYCLE
             END IF
          ENDIF

          LL = LL + 1

          ! E(LL) through E(M-1) are nonzero, E(LL-1) is zero

          IF ( LL == M-1 ) THEN

             ! 2 by 2 block, handle separately

             CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR, &
                         COSR, SINL, COSL )
             D( M-1 ) = SIGMX
             E( M-1 ) = ZERO
             D( M ) = SIGMN

             ! Compute singular vectors, if desired

             IF ( NCVT > 0 ) &
                  CALL DROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT, COSR, SINR )
             IF ( NRU > 0 ) &
                  CALL DROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
             IF ( NCC > 0 ) &
                  CALL DROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL, SINL )
             M = M - 2

             CYCLE
          END IF

          ! If working on new submatrix, choose shift direction
          ! (from larger end diagonal element towards smaller)

          IF ( LL > OLDM .OR. M < OLDLL ) THEN
             IF ( ABS( D( LL ) ) >= ABS( D( M ) ) ) THEN

                ! Chase bulge from top (big end) to bottom (small end)
                IDIR = 1

             ELSE

                ! Chase bulge from bottom (big end) to top (small end)
                IDIR = 2

             END IF
          END IF

          ! Apply convergence tests

          IF ( IDIR == 1 ) THEN

             ! Run convergence test in forward direction
             ! First apply standard test to bottom of matrix

             IF ( ( ABS( E( M-1 ) ) <= ABS( TOL )*ABS( D( M ) ) ) .OR. &
                  ( TOL < ZERO .AND. ABS( E( M-1 ) ) <= THRESH ) ) THEN
                E( M-1 ) = ZERO
                CYCLE
             END IF

             IF ( TOL >= ZERO ) THEN

                ! If relative accuracy desired,
                ! apply convergence criterion forward

                MU = ABS( D( LL ) )
                SMINL = MU
                DO LLL = LL, M - 1
                   IF ( ABS( E( LLL ) ) <= TOL*MU ) THEN
                      E( LLL ) = ZERO
                      CYCLE
                   END IF
                   MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
                   SMINL = MIN( SMINL, MU )
                END DO

             END IF

          ELSE

             ! Run convergence test in backward direction
             ! First apply standard test to top of matrix

             IF ( ( ABS( E( LL ) ) <= ABS( TOL )*ABS( D( LL ) ) ) .OR. &
                  ( TOL < ZERO .AND. ABS( E( LL ) ) <= THRESH ) ) THEN
                E( LL ) = ZERO
                CYCLE
             END IF

             IF ( TOL >= ZERO ) THEN

                ! If relative accuracy desired,
                ! apply convergence criterion backward

                MU = ABS( D( M ) )
                SMINL = MU
                DO LLL = M - 1, LL, -1
                   IF ( ABS( E( LLL ) ) <= TOL*MU ) THEN
                      E( LLL ) = ZERO
                      CYCLE
                   END IF
                   MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
                   SMINL = MIN( SMINL, MU )
                END DO
             END IF
          END IF
          OLDLL = LL
          OLDM = M

          ! Compute shift.  First, test if shifting would ruin relative
          ! accuracy, and if so set the shift to zero.

          IF ( TOL >= ZERO .AND. REAL(N,DP)*TOL*( SMINL / SMAX ) <= &
               MAX( EPS, HNDRTH*TOL ) ) THEN

             ! Use a zero shift to avoid loss of relative accuracy
             SHIFT = ZERO

          ELSE

             ! Compute the shift from 2-by-2 block at end of matrix

             IF ( IDIR == 1 ) THEN
                SLL = ABS( D( LL ) )
                CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
             ELSE
                SLL = ABS( D( M ) )
                CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
             END IF

             ! Test if shift negligible, and if so set to zero

             IF ( SLL > ZERO ) THEN
                IF ( ( SHIFT / SLL )**2 < EPS ) SHIFT = ZERO
             END IF

          END IF

          ! Increment iteration count

          ITER = ITER + M - LL

          ! If SHIFT = 0, do simplified QR iteration

          IF ( ABS(SHIFT) < ZEPS ) THEN

             IF ( IDIR == 1 ) THEN

                ! Chase bulge from top to bottom
                ! Save cosines and sines for later singular vector updates

                CS = ONE
                OLDCS = ONE
                DO I = LL, M - 1
                   CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
                   IF ( I > LL ) E( I-1 ) = OLDSN*R
                   CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
                   WORK( I-LL+1 ) = CS
                   WORK( I-LL+1+NM1 ) = SN
                   WORK( I-LL+1+NM12 ) = OLDCS
                   WORK( I-LL+1+NM13 ) = OLDSN
                END DO
                H = D( M )*CS
                D( M ) = H*OLDCS
                E( M-1 ) = H*OLDSN

                ! Update singular vectors

                IF ( NCVT > 0 ) &
                     CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                                 WORK( N ), VT( LL, 1 ), LDVT )
                IF ( NRU > 0 ) &
                     CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                                 WORK( NM13+1 ), U( 1, LL ), LDU )
                IF ( NCC > 0 ) &
                     CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                                 WORK( NM13+1 ), C( LL, 1 ), LDC )

                ! Test convergence

                IF ( ABS( E( M-1 ) ) <= THRESH ) E( M-1 ) = ZERO

             ELSE

                ! Chase bulge from bottom to top
                ! Save cosines and sines for later singular vector updates

                CS = ONE
                OLDCS = ONE
                DO I = M, LL + 1, -1
                   CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
                   IF ( I < M ) E( I ) = OLDSN*R
                   CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
                   WORK( I-LL ) = CS
                   WORK( I-LL+NM1 ) = -SN
                   WORK( I-LL+NM12 ) = OLDCS
                   WORK( I-LL+NM13 ) = -OLDSN
                END DO
                H = D( LL )*CS
                D( LL ) = H*OLDCS
                E( LL ) = H*OLDSN

                ! Update singular vectors

                IF ( NCVT > 0 ) &
                     CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                                 WORK( NM13+1 ), VT( LL, 1 ), LDVT )
                IF ( NRU > 0 ) &
                     CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                                 WORK( N ), U( 1, LL ), LDU )
                IF ( NCC > 0 ) &
                     CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                                 WORK( N ), C( LL, 1 ), LDC )

                ! Test convergence

                IF ( ABS( E( LL ) ) <= THRESH ) E( LL ) = ZERO

             END IF

          ELSE

             ! Use nonzero shift

             IF ( IDIR == 1 ) THEN

                ! Chase bulge from top to bottom
                ! Save cosines and sines for later singular vector updates

                F = ( ABS( D( LL ) )-SHIFT )*( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
                G = E( LL )
                DO I = LL, M - 1
                   CALL DLARTG( F, G, COSR, SINR, R )
                   IF ( I > LL ) E( I-1 ) = R
                   F = COSR*D( I ) + SINR*E( I )
                   E( I ) = COSR*E( I ) - SINR*D( I )
                   G = SINR*D( I+1 )
                   D( I+1 ) = COSR*D( I+1 )
                   CALL DLARTG( F, G, COSL, SINL, R )
                   D( I ) = R
                   F = COSL*E( I ) + SINL*D( I+1 )
                   D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
                   IF ( I < M-1 ) THEN
                      G = SINL*E( I+1 )
                      E( I+1 ) = COSL*E( I+1 )
                   END IF
                   WORK( I-LL+1 ) = COSR
                   WORK( I-LL+1+NM1 ) = SINR
                   WORK( I-LL+1+NM12 ) = COSL
                   WORK( I-LL+1+NM13 ) = SINL
                END DO
                E( M-1 ) = F

                ! Update singular vectors

                IF ( NCVT > 0 ) &
                     CALL DLASR( 'L', 'V', 'F', M-LL+1, NCVT, WORK( 1 ), &
                                 WORK( N ), VT( LL, 1 ), LDVT )
                IF ( NRU > 0 ) &
                     CALL DLASR( 'R', 'V', 'F', NRU, M-LL+1, WORK( NM12+1 ), &
                                 WORK( NM13+1 ), U( 1, LL ), LDU )
                IF ( NCC > 0 ) &
                     CALL DLASR( 'L', 'V', 'F', M-LL+1, NCC, WORK( NM12+1 ), &
                                 WORK( NM13+1 ), C( LL, 1 ), LDC )

                ! Test convergence

                IF ( ABS( E( M-1 ) ) <= THRESH ) E( M-1 ) = ZERO

             ELSE

                ! Chase bulge from bottom to top
                ! Save cosines and sines for later singular vector updates

                F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT / D( M ) )
                G = E( M-1 )
                DO I = M, LL + 1, -1
                   CALL DLARTG( F, G, COSR, SINR, R )
                   IF ( I < M ) E( I ) = R
                   F = COSR*D( I ) + SINR*E( I-1 )
                   E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
                   G = SINR*D( I-1 )
                   D( I-1 ) = COSR*D( I-1 )
                   CALL DLARTG( F, G, COSL, SINL, R )
                   D( I ) = R
                   F = COSL*E( I-1 ) + SINL*D( I-1 )
                   D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
                   IF ( I > LL+1 ) THEN
                      G = SINL*E( I-2 )
                      E( I-2 ) = COSL*E( I-2 )
                   END IF
                   WORK( I-LL ) = COSR
                   WORK( I-LL+NM1 ) = -SINR
                   WORK( I-LL+NM12 ) = COSL
                   WORK( I-LL+NM13 ) = -SINL
                END DO
                E( LL ) = F

                ! Test convergence

                IF ( ABS( E( LL ) ) <= THRESH ) E( LL ) = ZERO

                ! Update singular vectors if desired

                IF ( NCVT > 0 ) &
                     CALL DLASR( 'L', 'V', 'B', M-LL+1, NCVT, WORK( NM12+1 ), &
                                 WORK( NM13+1 ), VT( LL, 1 ), LDVT )
                IF ( NRU > 0 ) &
                     CALL DLASR( 'R', 'V', 'B', NRU, M-LL+1, WORK( 1 ), &
                                 WORK( N ), U( 1, LL ), LDU )
                IF ( NCC > 0 ) &
                     CALL DLASR( 'L', 'V', 'B', M-LL+1, NCC, WORK( 1 ), &
                                 WORK( N ), C( LL, 1 ), LDC )
             END IF

          END IF

          ! QR iteration finished, go back and check convergence

       END DO

    END IF

    ! All singular values converged, so make them positive

    DO I = 1, N
       IF ( D( I ) < ZERO ) THEN
          D( I ) = -D( I )

          ! Change sign of singular vectors, if desired
          IF ( NCVT > 0 ) CALL DSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
       END IF
    END DO

    ! Sort the singular values into decreasing order (insertion sort on
    ! singular values, but only one transposition per singular vector)

    DO I = 1, N - 1

       ! Scan for smallest D(I)
       ISUB = 1
       SMIN = D( 1 )
       DO J = 2, N + 1 - I
          IF ( D( J ) <= SMIN ) THEN
             ISUB = J
             SMIN = D( J )
          END IF
       END DO

       IF ( ISUB /= N+1-I ) THEN
          ! Swap singular values and vectors
          D( ISUB ) = D( N+1-I )
          D( N+1-I ) = SMIN
          IF ( NCVT > 0 ) &
               CALL DSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ), LDVT )
          IF ( NRU > 0 ) &
               CALL DSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
          IF ( NCC > 0 ) &
               CALL DSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
       END IF

    END DO

  END SUBROUTINE DBDSQR

!==========================================================================
!
!  DGEBRD
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGEBRD + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebrd.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), D( * ), E( * ), TAUP( * ), &
!                          TAUQ( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGEBRD reduces a general real M-by-N matrix A to upper or lower
! bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!
! If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows in the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns in the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the M-by-N general matrix to be reduced.
!          On exit,
!          if m >= n, the diagonal and the first superdiagonal are
!            overwritten with the upper bidiagonal matrix B; the
!            elements below the diagonal, with the array TAUQ, represent
!            the orthogonal matrix Q as a product of elementary
!            reflectors, and the elements above the first superdiagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors;
!          if m < n, the diagonal and the first subdiagonal are
!            overwritten with the lower bidiagonal matrix B; the
!            elements below the first subdiagonal, with the array TAUQ,
!            represent the orthogonal matrix Q as a product of
!            elementary reflectors, and the elements above the diagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors.
!          See Further Details.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] D
!
!          D is REAL(DP)         array, dimension (min(M,N))
!          The diagonal elements of the bidiagonal matrix B:
!          D(i) = A(i,i).
!
! param[out] E
!
!          E is REAL(DP)         array, dimension (min(M,N)-1)
!          The off-diagonal elements of the bidiagonal matrix B:
!          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
!
! param[out] TAUQ
!
!          TAUQ is REAL(DP)         array dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix Q. See Further Details.
!
! param[out] TAUP
!
!          TAUP is REAL(DP)         array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix P. See Further Details.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The length of the array WORK.  LWORK >= max(1,M,N).
!          For optimum performance LWORK >= (M+N)*NB, where NB
!          is the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!  =====================
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!  If m >= n,
!
!     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
!  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n,
!
!     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The contents of A on exit are illustrated by the following examples:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!    (  v1  v2  v3  v4  v5 )
!
!  where d and e denote diagonal and off-diagonal elements of B, vi
!  denotes an element of the vector defining H(i), and ui an element of
!  the vector defining G(i).
!
!=====================================================================

  SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DGEBD2, DGEMM, DLABRD, XERBLA, ILAENV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), D( * ), E( * ), TAUP( * ), &
                               TAUQ( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ONE = 1.0_DP

    ! .. Local Scalars ..
    LOGICAL  :: LQUERY
    INTEGER  :: I, IINFO, J, LDWRKX, LDWRKY, LWKOPT, MINMN, NB, NBMIN, NX, LWS

    ! .. Executable Statements ..

    ! Test the input parameters

    INFO = 0
    NB = MAX( 1, ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 ) )
    LWKOPT = ( M+N )*NB
    WORK( 1 ) = REAL( LWKOPT, DP )
    LQUERY = ( LWORK == -1 )

    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -4
    ELSE IF ( LWORK < MAX( 1, M, N ) .AND. .NOT.LQUERY ) THEN
       INFO = -10
    END IF

    IF ( INFO < 0 ) THEN
       CALL XERBLA( 'DGEBRD', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    MINMN = MIN( M, N )
    IF ( MINMN == 0 ) THEN
       WORK( 1 ) = 1
       RETURN
    END IF

    LWS = MAX( M, N )
    LDWRKX = M
    LDWRKY = N

    IF ( NB > 1 .AND. NB < MINMN ) THEN

       ! Set the crossover point NX.

       NX = MAX( NB, ILAENV( 3, 'DGEBRD', ' ', M, N, -1, -1 ) )

       ! Determine when to switch from blocked to unblocked code.

       IF ( NX < MINMN ) THEN
          LWS = LWKOPT
          IF ( LWORK < LWS ) THEN

             ! Not enough work space for the optimal NB, consider using
             ! a smaller block size.
             NBMIN = ILAENV( 2, 'DGEBRD', ' ', M, N, -1, -1 )
             IF ( LWORK >= ( M+N )*NBMIN ) THEN
                NB = LWORK / ( M+N )
             ELSE
                NB = 1
                NX = MINMN
             END IF

          END IF
       END IF

    ELSE

       NX = MINMN

    END IF

    DO I = 1, MINMN - NX, NB

       ! Reduce rows and columns i:i+nb-1 to bidiagonal form and return
       ! the matrices X and Y which are needed to update the unreduced
       ! part of the matrix

       CALL DLABRD( M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), &
                    TAUP( I ), WORK, LDWRKX, WORK( LDWRKX*NB+1 ), LDWRKY )

       ! Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
       ! of the form  A := A - V*Y**T - X*U**T

       CALL DGEMM( 'No transpose', 'Transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, A( I+NB, I ), &
                   LDA, WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE, A( I+NB, I+NB ), LDA )

       CALL DGEMM( 'No transpose', 'No transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, WORK( NB+1 ), &
                   LDWRKX, A( I, I+NB ), LDA, ONE, A( I+NB, I+NB ), LDA )

       ! Copy diagonal and off-diagonal elements of B back into A

       IF ( M >= N ) THEN
          DO J = I, I + NB - 1
             A( J, J ) = D( J )
             A( J, J+1 ) = E( J )
          END DO
       ELSE
          DO J = I, I + NB - 1
             A( J, J ) = D( J )
             A( J+1, J ) = E( J )
          END DO
       END IF
    END DO

    ! Use unblocked code to reduce the remainder of the matrix

    CALL DGEBD2( M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), &
                 TAUP( I ), WORK, IINFO )

    WORK( 1 ) = REAL( LWS, DP )

  END SUBROUTINE DGEBRD

!==========================================================================
!
!  DGELQF
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGELQF + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelqf.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGELQF computes an LQ factorization of a real M-by-N matrix A:
! A = L * Q.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and below the diagonal of the array
!          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
!          lower triangular if m <= n); the elements above the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] TAU
!
!          TAU is REAL(DP)         array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!  =====================
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

  SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DGELQ2, DLARFB, DLARFT, XERBLA, ILAENV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Local Scalars ..
    LOGICAL :: LQUERY
    INTEGER :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0
    NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
    LWKOPT = M*NB
    WORK( 1 ) = REAL( LWKOPT, DP )
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
       CALL XERBLA( 'DGELQF', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    K = MIN( M, N )
    IF ( K == 0 ) THEN
       WORK( 1 ) = 1._DP
       RETURN
    END IF

    NBMIN = 2
    NX = 0
    IWS = M

    IF ( NB > 1 .AND. NB < K ) THEN

       ! Determine when to cross over from blocked to unblocked code.

       NX = MAX( 0, ILAENV( 3, 'DGELQF', ' ', M, N, -1, -1 ) )

       IF ( NX < K ) THEN

          ! Determine if workspace is large enough for blocked code.

          LDWORK = M
          IWS = LDWORK*NB

          IF ( LWORK < IWS ) THEN
             ! Not enough workspace to use optimal NB:  reduce NB and
             ! determine the minimum value of NB.
             NB = LWORK / LDWORK
             NBMIN = MAX( 2, ILAENV( 2, 'DGELQF', ' ', M, N, -1, -1 ) )
          END IF

       END IF

    END IF

    IF ( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN

       ! Use blocked code initially

       DO I = 1, K - NX, NB
          IB = MIN( K-I+1, NB )

          ! Compute the LQ factorization of the current block
          ! A(i:i+ib-1,i:n)

          CALL DGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

          IF ( I+IB <= M ) THEN

             ! Form the triangular factor of the block reflector
             ! H = H(i) H(i+1) . . . H(i+ib-1)

             CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), LDA, &
                          TAU( I ), WORK, LDWORK )

             ! Apply H to A(i+ib:m,i:n) from the right

             CALL DLARFB( 'Right', 'No transpose', 'Forward', 'Rowwise', M-I-IB+1, &
                          N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), &
                          LDA, WORK( IB+1 ), LDWORK )
          END IF
       END DO

    ELSE

       I = 1

    END IF

    ! Use unblocked code to factor the last or only block.

    IF ( I <= K ) &
         CALL DGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

    WORK( 1 ) = REAL( IWS, DP )

  END SUBROUTINE DGELQF

!==========================================================================
!
!  DGEQRF
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGEQRF + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqrf.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGEQRF computes a QR factorization of a real M-by-N matrix A:
! A = Q * R.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of min(m,n) elementary reflectors (see Further
!          Details).
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] TAU
!
!          TAU is REAL(DP)         array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!  =====================
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

  SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DGEQR2, DLARFB, DLARFT, XERBLA, ILAENV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Local Scalars ..
    LOGICAL :: LQUERY
    INTEGER :: I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX

    ! .. Executable Statements ..

    ! Test the input arguments
    INFO = 0
    NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
    LWKOPT = N*NB
    WORK( 1 ) = REAL( LWKOPT, DP )
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
       CALL XERBLA( 'DGEQRF', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    K = MIN( M, N )
    IF ( K == 0 ) THEN
       WORK( 1 ) = 1._DP
       RETURN
    END IF

    NBMIN = 2
    NX = 0
    IWS = N

    IF ( NB > 1 .AND. NB < K ) THEN

       ! Determine when to cross over from blocked to unblocked code.

       NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )

       IF ( NX < K ) THEN

          ! Determine if workspace is large enough for blocked code.

          LDWORK = N
          IWS = LDWORK*NB

          IF ( LWORK < IWS ) THEN
             ! Not enough workspace to use optimal NB:  reduce NB and
             ! determine the minimum value of NB.
             NB = LWORK / LDWORK
             NBMIN = MAX( 2, ILAENV( 2, 'DGEQRF', ' ', M, N, -1, -1 ) )
          END IF

       END IF

    END IF

    IF ( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN

       ! Use blocked code initially

       DO I = 1, K - NX, NB
          IB = MIN( K-I+1, NB )

          ! Compute the QR factorization of the current block
          ! A(i:m,i:i+ib-1)

          CALL DGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )

          IF ( I+IB <= N ) THEN

             ! Form the triangular factor of the block reflector
             ! H = H(i) H(i+1) . . . H(i+ib-1)

             CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, A( I, I ), &
                          LDA, TAU( I ), WORK, LDWORK )

             ! Apply H**T to A(i:m,i+ib:n) from the left

             CALL DLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', &
                          M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, &
                          A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK )
          END IF
       END DO

    ELSE

       I = 1

    END IF

    ! Use unblocked code to factor the last or only block.

    IF ( I <= K ) &
         CALL DGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, IINFO )

    WORK( 1 ) = REAL( IWS, DP )

  END SUBROUTINE DGEQRF

!==========================================================================
!
!  DORGQR
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORGQR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORGQR generates an M-by-N real matrix Q with orthonormal columns,
! which is defined as the first N columns of a product of K elementary
! reflectors of order M
!
!       Q  =  H(1) H(2) . . . H(k)
!
! as returned by DGEQRF.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
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

  SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DLARFB, DLARFT, DORG2R, XERBLA, ILAENV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: K, LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP

    ! .. Local Scalars ..
    LOGICAL :: LQUERY
    INTEGER :: I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO = 0
    NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
    LWKOPT = MAX( 1, N )*NB
    WORK( 1 ) = REAL( LWKOPT, DP )
    LQUERY = ( LWORK == -1 )
    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 .OR. N > M ) THEN
       INFO = -2
    ELSE IF ( K < 0 .OR. K > N ) THEN
       INFO = -3
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -5
    ELSE IF ( LWORK < MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
       INFO = -8
    END IF
    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORGQR', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( N <= 0 ) THEN
       WORK( 1 ) = 1._DP
       RETURN
    END IF

    NBMIN = 2
    NX = 0
    IWS = N

    IF ( NB > 1 .AND. NB < K ) THEN

       ! Determine when to cross over from blocked to unblocked code.

       NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
       IF ( NX < K ) THEN

          ! Determine if workspace is large enough for blocked code.

          LDWORK = N
          IWS = LDWORK*NB
          IF ( LWORK < IWS ) THEN

             ! Not enough workspace to use optimal NB:  reduce NB and
             ! determine the minimum value of NB.

             NB = LWORK / LDWORK
             NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
          END IF
       END IF

    END IF

    IF ( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN

       ! Use blocked code after the last block.
       ! The first kk columns are handled by the block method.

       KI = ( ( K-NX-1 ) / NB )*NB
       KK = MIN( K, KI+NB )

       ! Set A(1:kk,kk+1:n) to zero.

       DO J = KK + 1, N
          DO I = 1, KK
             A( I, J ) = ZERO
          END DO
       END DO

    ELSE

       KK = 0

    END IF

    ! Use unblocked code for the last or only block.

    IF ( KK < N ) THEN
       CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO )
    ENDIF

    IF ( KK > 0 ) THEN

       ! Use blocked code

       DO I = KI + 1, 1, -NB
          IB = MIN( NB, K-I+1 )

          IF ( I+IB <= N ) THEN

             ! Form the triangular factor of the block reflector
             ! H = H(i) H(i+1) . . . H(i+ib-1)

             CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                          A( I, I ), LDA, TAU( I ), WORK, LDWORK )

             ! Apply H to A(i:m,i+ib:n) from the left

             CALL DLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise', &
                          M-I+1, N-I-IB+1, IB, A( I, I ), LDA, WORK, LDWORK, &
                          A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK )
          END IF

          ! Apply H to rows i:m of current block

          CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )

          ! Set rows 1:i-1 of current block to zero

          DO J = I, I + IB - 1
             DO L = 1, I - 1
                A( L, J ) = ZERO
             END DO
          END DO

       END DO
    END IF

    WORK( 1 ) = REAL( IWS, DP )

  END SUBROUTINE DORGQR

!=====================================================================
!
!  DLACPY copies all or part of one two-dimensional array to another.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLACPY + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), B( LDB,* )
!       ..
!
!  Purpose:
!  =============
!
! DLACPY copies all or part of a two-dimensional matrix A to another
! matrix B.
!
!  Arguments:
!  ==========
!
! param[in] UPLO
!
!          UPLO is CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] B
!
!          B is REAL(DP)         array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
! param[in] LDB
!
!          LDB is INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
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

  SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN) :: UPLO
    INTEGER,   INTENT(IN) :: LDA, LDB, M, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), B( LDB,* )

    ! .. Local Scalars ..
    INTEGER :: I, J

    ! .. Executable Statements ..

    IF ( LSAME( UPLO, 'U' ) ) THEN
       DO J = 1, N
          DO I = 1, MIN( J, M )
             B( I, J ) = A( I, J )
          END DO
       END DO
    ELSE IF ( LSAME( UPLO, 'L' ) ) THEN
       DO J = 1, N
          DO I = J, M
             B( I, J ) = A( I, J )
          END DO
       END DO
    ELSE
       DO J = 1, N
          DO I = 1, M
             B( I, J ) = A( I, J )
          END DO
       END DO
    END IF

  END SUBROUTINE DLACPY

!==========================================================================
!
!  DLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm,
! or the largest absolute value of any element of a general rectangular matrix.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLANGE + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlange.f
!
!  Definition:
!  ===========
!
!       REAL(DP)         FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
!
!        IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DLANGE  returns the value of the one norm,  or the Frobenius norm, or
! the  infinity norm,  or the  element of  largest absolute value  of a
! real matrix A.
!
! \return DLANGE
!    DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!             (
!             ( norm1(A),         NORM = '1', 'O' or 'o'
!             (
!             ( normI(A),         NORM = 'I' or 'i'
!             (
!             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
! where  norm1  denotes the  one norm of a matrix (maximum column sum),
! normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
! normF  denotes the  Frobenius norm of a matrix (square root of sum of
! squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!
!  Arguments:
!  ==========
!
! param[in] NORM
!
!          NORM is CHARACTER*1
!          Specifies the value to be returned in DLANGE as described
!          above.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          DLANGE is set to zero.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          DLANGE is set to zero.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          The m by n matrix A.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK)),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
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

  REAL(DP) FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, DLASSQ, DISNAN

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN) :: NORM
    INTEGER,   INTENT(IN) :: LDA, M, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    INTEGER  :: I, J
    REAL(DP) :: SCALE, SUM, VALUE

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
          CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
       END DO
       VALUE = SCALE*SQRT( SUM )

    END IF

    DLANGE = VALUE

  END FUNCTION DLANGE

!==========================================================================
!
!  DLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASCL + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlascl.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       INTEGER            INFO, KL, KU, LDA, M, N
!       REAL(DP)           CFROM, CTO
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* )
!       ..
!
!  Purpose:
!  =============
!
! DLASCL multiplies the M by N real matrix A by the real scalar
! CTO/CFROM.  This is done without over/underflow as long as the final
! result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
! A may be full, upper triangular, lower triangular, upper Hessenberg,
! or banded.
!
!  Arguments:
!  ==========
!
! param[in] TYPE
!
!          TYPE is CHARACTER*1
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
!                  bandwidth KU. See DGBTRF for storage details.
!
! param[in] KL
!
!          KL is INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
! param[in] KU
!
!          KU is INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
! param[in] CFROM
!
!          CFROM is REAL(DP)
!
! param[in] CTO
!
!          CTO is REAL(DP)
!
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] INFO
!
!          INFO is INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
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

  SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    USE IEEE_ARITHMETIC, ONLY: IEEE_QUIET_EQ
    IMPORT,              ONLY: LSAME, XERBLA, DISNAN, DLAMCH

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: TYPE
    INTEGER,   INTENT(IN)  :: KL, KU, LDA, M, N
    INTEGER,   INTENT(OUT) :: INFO
    REAL(DP),  INTENT(IN)  :: CFROM, CTO

    ! .. Array Arguments ..
    REAL(DP),  INTENT(INOUT) :: A( LDA,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL  :: DONE
    INTEGER  :: I, ITYPE, J, K1, K2, K3, K4
    REAL(DP) :: BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM

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
    ELSE IF ( ABS(CFROM) < EPS .OR. DISNAN(CFROM) ) THEN
       INFO = -4
    ELSE IF ( DISNAN(CTO) ) THEN
       INFO = -5
    ELSE IF ( M < 0 ) THEN
       INFO = -6
    ELSE IF ( (N < 0) .OR. (ITYPE == 4 .AND. N /= M ) .OR. &
                           (ITYPE == 5 .AND. N /= M ) ) THEN
       INFO = -7
    ELSE IF ( ITYPE <= 3 .AND. LDA < MAX( 1, M ) ) THEN
       INFO = -9
    ELSE IF ( ITYPE >= 4 ) THEN
       IF ( KL < 0 .OR. KL > MAX( M-1, 0 ) ) THEN
          INFO = -2
       ELSE IF ( ( KU < 0 ) .OR. ( KU > MAX( N-1, 0 ) ) .OR. &
                 ( (ITYPE == 4 .OR. ITYPE == 5) .AND. KL /= KU ) ) THEN
          INFO = -3
       ELSE IF ( ( ITYPE == 4 .AND. LDA < KL+1 ) .OR. &
                 ( ITYPE == 5 .AND. LDA < KU+1 ) .OR. &
                 ( ITYPE == 6 .AND. LDA < 2*KL+KU+1 ) ) THEN
          INFO = -9
       END IF
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DLASCL', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( N == 0 .OR. M == 0 ) RETURN

    ! Get machine parameters

    SMLNUM = DLAMCH( 'S' )
    BIGNUM = ONE / SMLNUM

    CFROMC = CFROM
    CTOC = CTO

    DONE = .FALSE.
    DO WHILE (.NOT. DONE)

       CFROM1 = CFROMC*SMLNUM
       IF ( IEEE_QUIET_EQ( CFROM1, CFROMC ) ) THEN
          ! CFROMC is an inf.  Multiply by a correctly signed zero for
          ! finite CTOC, or a NaN if CTOC is infinite.
          MUL = CTOC / CFROMC
          DONE = .TRUE.
          CTO1 = CTOC
       ELSE
          CTO1 = CTOC / BIGNUM
          IF ( IEEE_QUIET_EQ( CTO1, CTOC ) ) THEN
             ! CTOC is either 0 or an inf.  In both cases, CTOC itself
             ! serves as the correct multiplication factor.
             MUL = CTOC
             DONE = .TRUE.
             CFROMC = ONE
          ELSE IF ( ABS( CFROM1 ) > ABS( CTOC ) .AND. ABS(CTOC) >= EPS ) THEN
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

  END SUBROUTINE DLASCL

!==========================================================================
!
!  DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASET + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       REAL(DP)           ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* )
!       ..
!
!  Purpose:
!  =============
!
! DLASET initializes an m-by-n matrix A to BETA on the diagonal and
! ALPHA on the offdiagonals.
!
!  Arguments:
!  ==========
!
! param[in] UPLO
!
!          UPLO is CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in] ALPHA
!
!          ALPHA is REAL(DP)
!          The constant to which the offdiagonal elements are to be set.
!
! param[in] BETA
!
!          BETA is REAL(DP)
!          The constant to which the diagonal elements are to be set.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i /= j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
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

  SUBROUTINE DLASET( UPLOC, M, N, ALPHA, BETA, A, LDA )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: UPLOC
    INTEGER,      INTENT(IN) :: LDA, M, N
    REAL(DP),     INTENT(IN) :: ALPHA, BETA

    ! .. Array Arguments ..
    REAL(DP),  INTENT(INOUT) :: A( LDA,* )

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

  END SUBROUTINE DLASET

!==========================================================================
!
!  DORGBR
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORGBR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgbr.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          VECT
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORGBR generates one of the real orthogonal matrices Q or P**T
! determined by DGEBRD when reducing a real matrix A to bidiagonal
! form: A = Q * B * P**T.  Q and P**T are defined as products of
! elementary reflectors H(i) or G(i) respectively.
!
! If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
! is of order M:
! if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
! columns of Q, where m >= n >= k;
! if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
! M-by-M matrix.
!
! If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
! is of order N:
! if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
! rows of P**T, where n >= m >= k;
! if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
! an N-by-N matrix.
!
!  Arguments:
!  ==========
!
! param[in] VECT
!
!          VECT is CHARACTER*1
!          Specifies whether the matrix Q or the matrix P**T is
!          required, as defined in the transformation applied by DGEBRD:
!          = 'Q':  generate Q;
!          = 'P':  generate P**T.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix Q or P**T to be returned.
!          M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix Q or P**T to be returned.
!          N >= 0.
!          If VECT = 'Q', M >= N >= min(M,K);
!          if VECT = 'P', N >= M >= min(N,K).
!
! param[in] K
!
!          K is INTEGER
!          If VECT = 'Q', the number of columns in the original M-by-K
!          matrix reduced by DGEBRD.
!          If VECT = 'P', the number of rows in the original K-by-N
!          matrix reduced by DGEBRD.
!          K >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DGEBRD.
!          On exit, the M-by-N matrix Q or P**T.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension
!                                (min(M,K)) if VECT = 'Q'
!                                (min(N,K)) if VECT = 'P'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i) or G(i), which determines Q or P**T, as
!          returned by DGEBRD in its array argument TAUQ or TAUP.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
!          For optimum performance LWORK >= min(M,N)*NB, where NB
!          is the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! April 2012
!
!=====================================================================

  SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, ILAENV, DORGLQ, DORGQR

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: VECT
    INTEGER,   INTENT(IN)  :: K, LDA, LWORK, M, N
    INTEGER,   INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    LOGICAL :: LQUERY, WANTQ
    INTEGER :: I, IINFO, J, LWKOPT, MN

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO = 0
    WANTQ = LSAME( VECT, 'Q' )
    MN = MIN( M, N )
    LQUERY = ( LWORK == -1 )

    IF ( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
       INFO = -1
    ELSE IF ( M < 0 ) THEN
       INFO = -2
    ELSE IF ( N < 0 .OR. ( WANTQ .AND. ( N > M .OR. N < MIN( M, K ) ) ) .OR. &
                    ( .NOT.WANTQ .AND. ( M > N .OR. M < MIN( N, K ) ) ) ) THEN
       INFO = -3
    ELSE IF ( K < 0 ) THEN
       INFO = -4
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -6
    ELSE IF ( LWORK < MAX( 1, MN ) .AND. .NOT.LQUERY ) THEN
       INFO = -9
    END IF

    IF ( INFO == 0 ) THEN
       WORK( 1 ) = 1
       IF ( WANTQ ) THEN
          IF ( M >= K ) THEN
             CALL DORGQR( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
          ELSE
             IF ( M > 1 ) THEN
                CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, -1, IINFO )
             END IF
          END IF
       ELSE
          IF ( K < N ) THEN
             CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, -1, IINFO )
          ELSE
             IF ( N > 1 ) THEN
                CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, -1, IINFO )
             END IF
          END IF
       END IF
       LWKOPT = INT( WORK( 1 ) )
       LWKOPT = MAX (LWKOPT, MN)
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORGBR', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       WORK( 1 ) = REAL( LWKOPT, DP )
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 ) THEN
       WORK( 1 ) = 1
       RETURN
    END IF

    IF ( WANTQ ) THEN

       ! Form Q, determined by a call to DGEBRD to reduce an m-by-k
       ! matrix

       IF ( M >= K ) THEN

          ! If m >= k, assume m >= n >= k

          CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )

       ELSE

          ! If m < k, assume m = n

          ! Shift the vectors which define the elementary reflectors one
          ! column to the right, and set the first row and column of Q
          ! to those of the unit matrix

          DO J = M, 2, -1
             A( 1, J ) = ZERO
             DO I = J + 1, M
                A( I, J ) = A( I, J-1 )
             END DO
          END DO

          A( 1, 1 ) = ONE
          DO I = 2, M
             A( I, 1 ) = ZERO
          END DO

          IF ( M > 1 ) THEN
             ! Form Q(2:m,2:m)
             CALL DORGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
          END IF

       END IF

    ELSE

       ! Form P**T, determined by a call to DGEBRD to reduce a k-by-n
       ! matrix

       IF ( K < N ) THEN

          ! If k < n, assume k <= m <= n
          CALL DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )

       ELSE

          ! If k >= n, assume m = n

          ! Shift the vectors which define the elementary reflectors one
          ! row downward, and set the first row and column of P**T to
          ! those of the unit matrix

          A( 1, 1 ) = ONE
          DO I = 2, N
             A( I, 1 ) = ZERO
          END DO

          DO J = 2, N
             DO I = J - 1, 2, -1
                A( I, J ) = A( I-1, J )
             END DO
             A( 1, J ) = ZERO
          END DO

          IF ( N > 1 ) THEN
             ! Form P**T(2:n,2:n)
             CALL DORGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
          END IF

       END IF

    END IF

    WORK( 1 ) = REAL( LWKOPT, DP )

  END SUBROUTINE DORGBR

!==========================================================================
!
!  DORGLQ
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORGLQ + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorglq.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
! which is defined as the first M rows of a product of K elementary
! reflectors of order N
!
!       Q  =  H(k) . . . H(2) H(1)
!
! as returned by DGELQF.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix Q. N >= M.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. M >= K >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the i-th row must contain the vector which defines
!          the elementary reflector H(i), for i = 1,2,...,k, as returned
!          by DGELQF in the first k rows of its array argument A.
!          On exit, the M-by-N matrix Q.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
!          The dimension of the array WORK. LWORK >= max(1,M).
!          For optimum performance LWORK >= M*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
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

  SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: XERBLA, DLARFB, DLARFT, DORGL2, ILAENV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: K, LDA, LWORK, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP

    ! .. Local Scalars ..
    LOGICAL :: LQUERY
    INTEGER :: I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO = 0
    NB = ILAENV( 1, 'DORGLQ', ' ', M, N, K, -1 )
    LWKOPT = MAX( 1, M )*NB
    WORK( 1 ) = REAL( LWKOPT, DP )
    LQUERY = ( LWORK == -1 )

    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < M ) THEN
       INFO = -2
    ELSE IF ( K < 0 .OR. K > M ) THEN
       INFO = -3
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -5
    ELSE IF ( LWORK < MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
       INFO = -8
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORGLQ', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( M <= 0 ) THEN
       WORK( 1 ) = 1._DP
       RETURN
    END IF

    NBMIN = 2
    NX = 0
    IWS = M

    IF ( NB > 1 .AND. NB < K ) THEN

       ! Determine when to cross over from blocked to unblocked code.

       NX = MAX( 0, ILAENV( 3, 'DORGLQ', ' ', M, N, K, -1 ) )
       IF ( NX < K ) THEN

          ! Determine if workspace is large enough for blocked code.

          LDWORK = M
          IWS = LDWORK*NB
          IF ( LWORK < IWS ) THEN

             ! Not enough workspace to use optimal NB:  reduce NB and
             ! determine the minimum value of NB.

             NB = LWORK / LDWORK
             NBMIN = MAX( 2, ILAENV( 2, 'DORGLQ', ' ', M, N, K, -1 ) )
          END IF
       END IF
    END IF

    IF ( NB >= NBMIN .AND. NB < K .AND. NX < K ) THEN

       ! Use blocked code after the last block.
       ! The first kk rows are handled by the block method.

       KI = ( ( K-NX-1 ) / NB )*NB
       KK = MIN( K, KI+NB )

       ! Set A(kk+1:m,1:kk) to zero.

       DO J = 1, KK
          DO I = KK + 1, M
             A( I, J ) = ZERO
          END DO
       END DO

    ELSE

       KK = 0

    END IF

    ! Use unblocked code for the last or only block.

    IF ( KK < M ) THEN
       CALL DORGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO )
    ENDIF

    IF ( KK > 0 ) THEN

       ! Use blocked code

       DO I = KI + 1, 1, -NB
          IB = MIN( NB, K-I+1 )

          IF ( I+IB <= M ) THEN

             ! Form the triangular factor of the block reflector
             ! H = H(i) H(i+1) . . . H(i+ib-1)

             CALL DLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), &
                          LDA, TAU( I ), WORK, LDWORK )

             ! Apply H**T to A(i+ib:m,i:n) from the right

             CALL DLARFB( 'Right', 'Transpose', 'Forward', 'Rowwise', M-I-IB+1, &
                          N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), &
                          LDA, WORK( IB+1 ), LDWORK )
          END IF

          ! Apply H**T to columns i:n of current block

          CALL DORGL2( IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )

          ! Set columns 1:i-1 of current block to zero

          DO J = 1, I - 1
             DO L = I, I + IB - 1
                A( L, J ) = ZERO
             END DO
          END DO

       END DO

    END IF

    WORK( 1 ) = REAL( IWS, DP )

  END SUBROUTINE DORGLQ

!==========================================================================
!
!  DORMBR
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORMBR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormbr.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, VECT
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
! with
!                 SIDE = 'L'     SIDE = 'R'
! TRANS = 'N':      Q * C          C * Q
! TRANS = 'T':      Q**T * C       C * Q**T
!
! If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
! with
!                 SIDE = 'L'     SIDE = 'R'
! TRANS = 'N':      P * C          C * P
! TRANS = 'T':      P**T * C       C * P**T
!
! Here Q and P**T are the orthogonal matrices determined by DGEBRD when
! reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
! P**T are defined as products of elementary reflectors H(i) and G(i)
! respectively.
!
! Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
! order of the orthogonal matrix Q or P**T that is applied.
!
! If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
! if nq >= k, Q = H(1) H(2) . . . H(k);
! if nq < k, Q = H(1) H(2) . . . H(nq-1).
!
! If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
! if k < nq, P = G(1) G(2) . . . G(k);
! if k >= nq, P = G(1) G(2) . . . G(nq-1).
!
!  Arguments:
!  ==========
!
! param[in] VECT
!
!          VECT is CHARACTER*1
!          = 'Q': apply Q or Q**T;
!          = 'P': apply P or P**T.
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': apply Q, Q**T, P or P**T from the Left;
!          = 'R': apply Q, Q**T, P or P**T from the Right.
!
! param[in] TRANS
!
!          TRANS is CHARACTER*1
!          = 'N':  No transpose, apply Q  or P;
!          = 'T':  Transpose, apply Q**T or P**T.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C. N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          If VECT = 'Q', the number of columns in the original
!          matrix reduced by DGEBRD.
!          If VECT = 'P', the number of rows in the original
!          matrix reduced by DGEBRD.
!          K >= 0.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension
!                                (LDA,min(nq,K)) if VECT = 'Q'
!                                (LDA,nq)        if VECT = 'P'
!          The vectors which define the elementary reflectors H(i) and
!          G(i), whose products determine the matrices Q and P, as
!          returned by DGEBRD.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.
!          If VECT = 'Q', LDA >= max(1,nq);
!          if VECT = 'P', LDA >= max(1,min(nq,K)).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (min(nq,K))
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i) or G(i) which determines Q or P, as returned
!          by DGEBRD in the array argument TAUQ or TAUP.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
!          or P*C or P**T*C or C*P or C*P**T.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
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
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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

  SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C, &
                     LDC, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, ILAENV, DORMLQ, DORMQR

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: SIDE, TRANS, VECT
    INTEGER,   INTENT(IN)  :: K, LDA, LDC, LWORK, M, N
    INTEGER,   INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Local Scalars ..
    LOGICAL   :: APPLYQ, LEFT, LQUERY, NOTRAN
    CHARACTER :: TRANST
    INTEGER   :: I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO = 0
    APPLYQ = LSAME( VECT, 'Q' )
    LEFT = LSAME( SIDE, 'L' )
    NOTRAN = LSAME( TRANS, 'N' )
    LQUERY = ( LWORK == -1 )

    ! NQ is the order of Q or P and NW is the minimum dimension of WORK

    IF ( LEFT ) THEN
       NQ = M
       NW = N
    ELSE
       NQ = N
       NW = M
    END IF
    IF ( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
       INFO = -1
    ELSE IF ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
       INFO = -2
    ELSE IF ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
       INFO = -3
    ELSE IF ( M < 0 ) THEN
       INFO = -4
    ELSE IF ( N < 0 ) THEN
       INFO = -5
    ELSE IF ( K < 0 ) THEN
       INFO = -6
    ELSE IF ( ( APPLYQ .AND. LDA < MAX( 1, NQ ) ) .OR. &
         ( .NOT.APPLYQ .AND. LDA < MAX( 1, MIN( NQ, K ) ) ) ) THEN
       INFO = -8
    ELSE IF ( LDC < MAX( 1, M ) ) THEN
       INFO = -11
    ELSE IF ( LWORK < MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
       INFO = -13
    END IF

    IF ( INFO == 0 ) THEN
       IF ( APPLYQ ) THEN
          IF ( LEFT ) THEN
             NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M-1, N, M-1, -1 )
          ELSE
             NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N-1, N-1, -1 )
          END IF
       ELSE
          IF ( LEFT ) THEN
             NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M-1, N, M-1, -1 )
          ELSE
             NB = ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N-1, N-1, -1 )
          END IF
       END IF
       LWKOPT = MAX( 1, NW )*NB
       WORK( 1 ) = REAL( LWKOPT, DP )
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORMBR', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    WORK( 1 ) = 1._DP
    IF ( M == 0 .OR. N == 0 ) RETURN

    IF ( APPLYQ ) THEN

       ! Apply Q

       IF ( NQ >= K ) THEN

          ! Q was determined by a call to DGEBRD with nq >= k
          CALL DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO )

       ELSE IF ( NQ > 1 ) THEN

          ! Q was determined by a call to DGEBRD with nq < k

          IF ( LEFT ) THEN
             MI = M - 1
             NI = N
             I1 = 2
             I2 = 1
          ELSE
             MI = M
             NI = N - 1
             I1 = 1
             I2 = 2
          END IF
          CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), &
                       LDC, WORK, LWORK, IINFO )
       END IF

    ELSE

       ! Apply P

       IF ( NOTRAN ) THEN
          TRANST = 'T'
       ELSE
          TRANST = 'N'
       END IF

       IF ( NQ > K ) THEN

          ! P was determined by a call to DGEBRD with nq > k
          CALL DORMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, IINFO )

       ELSE IF ( NQ > 1 ) THEN

          ! P was determined by a call to DGEBRD with nq <= k

          IF ( LEFT ) THEN
             MI = M - 1
             NI = N
             I1 = 2
             I2 = 1
          ELSE
             MI = M
             NI = N - 1
             I1 = 1
             I2 = 2
          END IF
          CALL DORMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C( I1, I2 ), &
                       LDC, WORK, LWORK, IINFO )
       END IF

    END IF

   WORK( 1 ) = REAL( LWKOPT, DP )

  END SUBROUTINE DORMBR

!==========================================================================
!
!  DISNAN tests input for NaN.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DISNAN + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
!
!        IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       REAL(DP)           DIN
!       ..
!
!  Purpose:
!  =============
!
! DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
! otherwise.  To be replaced by the Fortran 2003 intrinsic in the
! future.
!
!  Arguments:
!  ==========
!
! param[in] DIN
!
!          DIN is REAL(DP)
!          Input to test for NaN.
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

  LOGICAL FUNCTION DISNAN( DIN )

    USE IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN) :: DIN

    ! .. Executable Statements ..

    DISNAN = IEEE_IS_NAN( DIN )

  END FUNCTION DISNAN

!=====================================================================
!
!  DORGL2
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORGL2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgl2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORGL2 generates an m by n real matrix Q with orthonormal rows,
! which is defined as the first m rows of a product of k elementary
! reflectors of order n
!
!       Q  =  H(k) . . . H(2) H(1)
!
! as returned by DGELQF.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix Q. N >= M.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. M >= K >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the i-th row must contain the vector which defines
!          the elementary reflector H(i), for i = 1,2,...,k, as returned
!          by DGELQF in the first k rows of its array argument A.
!          On exit, the m-by-n matrix Q.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (M)
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
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

  SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DSCAL, XERBLA, DLARF

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: K, LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    INTEGER :: I, J, L

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO = 0
    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < M ) THEN
       INFO = -2
    ELSE IF ( K < 0 .OR. K > M ) THEN
       INFO = -3
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -5
    END IF
    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORGL2', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( M <= 0 ) RETURN

    IF ( K < M ) THEN

       ! Initialise rows k+1:m to rows of the unit matrix

       DO J = 1, N
          DO L = K + 1, M
             A( L, J ) = ZERO
          END DO
          IF ( J > K .AND. J <= M ) A( J, J ) = ONE
       END DO
    END IF

    DO I = K, 1, -1

       ! Apply H(i) to A(i:m,i:n) from the right

       IF ( I < N ) THEN
          IF ( I < M ) THEN
             A( I, I ) = ONE
             CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), &
                         A( I+1, I ), LDA, WORK )
          END IF
          CALL DSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
       END IF

       A( I, I ) = ONE - TAU( I )

       ! Set A(i,1:i-1) to zero

       DO L = 1, I - 1
          A( I, L ) = ZERO
       END DO

    END DO

  END SUBROUTINE DORGL2

!==========================================================================
!
!  DLARFT forms the triangular factor T of a block reflector H = I - vtvH
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLARFT + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarft.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       INTEGER            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           T( LDT,* ), TAU( * ), V( LDV,* )
!       ..
!
!  Purpose:
!  =============
!
! DLARFT forms the triangular factor T of a real block reflector H
! of order n, which is defined as a product of k elementary reflectors.
!
! If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
! If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
! If STOREV = 'C', the vector which defines the elementary reflector
! H(i) is stored in the i-th column of the array V, and
!
!    H  =  I - V * T * V**T
!
! If STOREV = 'R', the vector which defines the elementary reflector
! H(i) is stored in the i-th row of the array V, and
!
!    H  =  I - V**T * T * V
!
!  Arguments:
!  ==========
!
! param[in] DIRECT
!
!          DIRECT is CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
! param[in] STOREV
!
!          STOREV is CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
! param[in] N
!
!          N is INTEGER
!          The order of the block reflector H. N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
! param[in] V
!
!          V is REAL(DP)         array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
! param[in] LDV
!
!          LDV is INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
! param[out] T
!
!          T is REAL(DP)         array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
! param[in] LDT
!
!          LDT is INTEGER
!          The leading dimension of the array T. LDT >= K.
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
!  Further Details:
!=====================
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored.
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

  SUBROUTINE DLARFT( DIRECTC, STOREVC, N, K, V, LDV, TAU, T, LDT )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, DGEMV, DTRMV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: DIRECTC, STOREVC
    INTEGER,      INTENT(IN) :: K, LDT, LDV, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: T( LDT,* ), TAU( * ), V( LDV,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

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

                CALL DGEMV( 'Transpose', J-I, I-1, -TAU( I ), V( I+1,1 ), LDV, &
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

                CALL DGEMV( 'No transpose', I-1, J-I, -TAU( I ), V( 1,I+1 ), &
                            LDV, V( I,I+1 ), LDV, ONE, T( 1,I ), 1 )

             END IF

             ! T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

             CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1,I ), 1 )

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

                   CALL DGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ), V( J,I+1 ), &
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

                   CALL DGEMV( 'No transpose', K-I, N-K+I-J, -TAU( I ), V( I+1,J ), &
                               LDV, V( I,J ), LDV, ONE, T( I+1,I ), 1 )

                END IF

                ! T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

                CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, T( I+1,I+1 ), &
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

  END SUBROUTINE DLARFT

!==========================================================================
!
!  DLARFB applies a block reflector or its transpose to a general rectangular matrix.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLARFB + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           C( LDC,* ), T( LDT,* ), V( LDV,* ), &
!                          WORK( LDWORK,* )
!       ..
!
!  Purpose:
!  =============
!
! DLARFB applies a real block reflector H or its transpose H**T to a
! real m by n matrix C, from either the left or the right.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': apply H or H**T from the Left
!          = 'R': apply H or H**T from the Right
!
! param[in] TRANS
!
!          TRANS is CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H**T (Transpose)
!
! param[in] DIRECT
!
!          DIRECT is CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
! param[in] STOREV
!
!          STOREV is CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C.
!
! param[in] K
!
!          K is INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
! param[in] V
!
!          V is REAL(DP)         array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See Further Details.
!
! param[in] LDV
!
!          LDV is INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
! param[in] T
!
!          T is REAL(DP)         array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
! param[in] LDT
!
!          LDT is INTEGER
!          The leading dimension of the array T. LDT >= K.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (LDWORK,K)
!
! param[in] LDWORK
!
!          LDWORK is INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  Authors:
!  ========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! June 2013
!
!  Further Details:
!=====================
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

  SUBROUTINE DLARFB( SIDEC, TRANSC, DIRECTC, STOREVC, M, N, K, V, LDV, &
                     T, LDT, C, LDC, WORK, LDWORK )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, DGEMM, DCOPY, DTRMM

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.5.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! June 2013

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: DIRECTC, SIDEC, STOREVC, TRANSC
    INTEGER,      INTENT(IN) :: K, LDC, LDT, LDV, LDWORK, M, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: C( LDC,* ), T( LDT,* ), V( LDV,* ), WORK( LDWORK,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ONE = 1.0_DP

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
                CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1

             CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                         K, ONE, V, LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C2**T * V2
                CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V * W**T

             IF ( M > K ) THEN
                ! C2 := C2 - V2 * W**T
                CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                            -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, &
                            ONE, C( K+1, 1 ), LDC )
             END IF

             ! W := W * V1**T

             CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
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
                CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1

             CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, K, &
                        ONE, V, LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C2 * V2
                CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V**T

             IF ( N > K ) THEN
                ! C2 := C2 - W * V2**T
                CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                            C( 1, K+1 ), LDC )
             END IF

             ! W := W * V1**T

             CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
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
                CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2

             CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, K, &
                         ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C1**T * V1
                CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V * W**T

             IF ( M > K ) THEN
                ! C1 := C1 - V1 * W**T
                CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                            -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
             END IF

             ! W := W * V2**T

             CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                         ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )

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
                CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2

             CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, K, &
                         ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C1 * V1
                CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V**T

             IF ( N > K ) THEN
                ! C1 := C1 - W * V1**T
                CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
             END IF

             ! W := W * V2**T

             CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                         ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )

             ! C2 := C2 - W

             DO J = 1, K
                DO I = 1, M
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
                CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1**T

             CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                         ONE, V, LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C2**T * V2**T
                CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                     C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V**T * W**T

             IF ( M > K ) THEN
                ! C2 := C2 - V2**T * W**T
                CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                            V( 1, K+1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC )
             END IF

             ! W := W * V1

             CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, K, &
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
                CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V1**T

             CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                         ONE, V, LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C2 * V2**T
                CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                            ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V

             IF ( N > K ) THEN
                ! C2 := C2 - W * V2
                CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                            C( 1, K+1 ), LDC )
             END IF

             ! W := W * V1

             CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, K, &
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
                CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2**T

             CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                         ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )

             IF ( M > K ) THEN
                ! W := W + C1**T * V1**T
                CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                            C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T**T  or  W * T

             CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - V**T * W**T

             IF ( M > K ) THEN
                ! C1 := C1 - V1**T * W**T
                CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                            V, LDV, WORK, LDWORK, ONE, C, LDC )
             END IF

             ! W := W * V2

             CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, K, &
                         ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )

             ! C2 := C2 - W**T

             DO J = 1, K
                DO I = 1, N
                   C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
                END DO
             END DO

          ELSE IF ( LSAME( SIDE, 'R' ) ) THEN

             ! Form  C * H  or  C * H'  where  C = ( C1  C2 )
             ! W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
             ! W := C2

             DO J = 1, K
                CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
             END DO

             ! W := W * V2**T

             CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                         ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )

             IF ( N > K ) THEN
                ! W := W + C1 * V1**T
                CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
             END IF

             ! W := W * T  or  W * T**T

             CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                         ONE, T, LDT, WORK, LDWORK )

             ! C := C - W * V

             IF ( N > K ) THEN
                ! C1 := C1 - W * V1
                CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
             END IF

             ! W := W * V2

             CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, K, &
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

  END SUBROUTINE DLARFB

!==========================================================================
!
!  DGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorithm.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGELQ2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelq2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGELQ2 computes an LQ factorization of a real m by n matrix A:
! A = L * Q.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and below the diagonal of the array
!          contain the m by min(m,n) lower trapezoidal matrix L (L is
!          lower triangular if m <= n); the elements above the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] TAU
!
!          TAU is REAL(DP)         array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (M)
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
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
!  Further Details:
!=====================
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

  SUBROUTINE DGELQ2( M, N, A, LDA, TAU, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: XERBLA, DLARF, DLARFG

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ONE = 1.0_DP

    ! .. Local Scalars ..
    INTEGER  :: I, K
    REAL(DP) :: AII

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
       CALL XERBLA( 'DGELQ2', -INFO )
       RETURN
    END IF

    K = MIN( M, N )

    DO I = 1, K

       ! Generate elementary reflector H(i) to annihilate A(i,i+1:n)

       CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, TAU( I ) )

       IF ( I < M ) THEN

          ! Apply H(i) to A(i+1:m,i:n) from the right
          AII = A( I, I )
          A( I, I ) = ONE
          CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), &
                      A( I+1, I ), LDA, WORK )
          A( I, I ) = AII

       END IF
    END DO

  END SUBROUTINE DGELQ2

!==========================================================================
!
!  DLASSQ updates a sum of squares represented in scaled form.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASSQ + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL(DP)           SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       REAL(DP)           X( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASSQ  returns the values  scl  and  smsq  such that
!
!    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
! where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
! assumed to be non-negative and  scl  returns the value
!
!    scl = max( scale, abs( x( i ) ) ).
!
! scale and sumsq must be supplied in SCALE and SUMSQ and
! scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
! The routine makes only one pass through the vector x.
!
!  Arguments:
!  ==========
!
! param[in] N
!
!          N is INTEGER
!          The number of elements to be used from the vector X.
!
! param[in] X
!
!          X is REAL(DP)         array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
! param[in] INCX
!
!          INCX is INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
! param[in,out] SCALE
!
!          SCALE is REAL(DP)
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
! param[in,out] SUMSQ
!
!          SUMSQ is REAL(DP)
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
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

  SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER,  INTENT(IN)    :: INCX, N
    REAL(DP), INTENT(INOUT) :: SCALE, SUMSQ

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: X( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: IX
    REAL(DP) :: ABSXI

    ! .. Executable Statements ..

    IF ( N > 0 ) THEN
       DO IX = 1, 1 + ( N-1 )*INCX, INCX
          ABSXI = ABS( X( IX ) )
          IF ( ABSXI >= EPS ) THEN
             IF ( SCALE < ABSXI ) THEN
                SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                SCALE = ABSXI
             ELSE
                SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
             END IF
          END IF
       END DO
    END IF

  END SUBROUTINE DLASSQ

!==========================================================================
!
!  DORMQR
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORMQR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormqr.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORMQR overwrites the general real M-by-N matrix C with
!
!                 SIDE = 'L'     SIDE = 'R'
! TRANS = 'N':      Q * C          C * Q
! TRANS = 'T':      Q**T * C       C * Q**T
!
! where Q is a real orthogonal matrix defined as the product of k
! elementary reflectors
!
!       Q = H(1) H(2) . . . H(k)
!
! as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
! if SIDE = 'R'.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
! param[in] TRANS
!
!          TRANS is CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C. N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
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
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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

  SUBROUTINE DORMQR( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, &
                     WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, ILAENV, DLARFB, DLARFT, DORM2R

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, LWORK, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(OUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    INTEGER, PARAMETER :: NBMAX = 64
    INTEGER, PARAMETER :: LDT   = NBMAX+1

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, LQUERY, NOTRAN
    INTEGER   :: I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, LWKOPT, MI, &
                 NB, NBMIN, NI, NQ, NW
    CHARACTER :: SIDE, TRANS

    ! .. Local Arrays ..
    REAL(DP)  :: T( LDT, NBMAX )

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
       NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K, -1 ) )
       LWKOPT = MAX( 1, NW )*NB
       WORK( 1 ) = REAL( LWKOPT, DP )
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORMQR', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) THEN
       WORK( 1 ) = 1._DP
       RETURN
    END IF

    NBMIN = 2
    LDWORK = NW
    IF ( NB > 1 .AND. NB < K ) THEN
       IWS = NW*NB
       IF ( LWORK < IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K, -1 ) )
       END IF
    ELSE
       IWS = NW
    END IF

    IF ( NB < NBMIN .OR. NB >= K ) THEN

       ! Use unblocked code
       CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )

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

          CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), LDA, &
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

          CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, IB, &
                       A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK, LDWORK )
       END DO

    END IF

    WORK( 1 ) = REAL( LWKOPT, DP )

  END SUBROUTINE DORMQR

!==========================================================================
!
!  DORMLQ
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORMLQ + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormlq.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORMLQ overwrites the general real M-by-N matrix C with
!
!                 SIDE = 'L'     SIDE = 'R'
! TRANS = 'N':      Q * C          C * Q
! TRANS = 'T':      Q**T * C       C * Q**T
!
! where Q is a real orthogonal matrix defined as the product of k
! elementary reflectors
!
!       Q = H(k) . . . H(2) H(1)
!
! as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
! if SIDE = 'R'.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
! param[in] TRANS
!
!          TRANS is CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C. N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGELQF in the first k rows of its array argument A.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! param[in] LWORK
!
!          LWORK is INTEGER
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
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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

  SUBROUTINE DORMLQ( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, &
                     WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, ILAENV, DLARFB, DLARFT, DORML2

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.0)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, LWORK, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP) , INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

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
    REAL(DP)  :: T( LDT, NBMAX )

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
       NB = MIN( NBMAX, ILAENV( 1, 'DORMLQ', SIDE // TRANS, M, N, K, -1 ) )
       LWKOPT = MAX( 1, NW )*NB
       WORK( 1 ) = REAL( LWKOPT, DP )
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORMLQ', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) THEN
       WORK( 1 ) = 1._DP
       RETURN
    END IF

    NBMIN = 2
    LDWORK = NW

    IF ( NB > 1 .AND. NB < K ) THEN
       IWS = NW*NB
       IF ( LWORK < IWS ) THEN
          NB = LWORK / LDWORK
          NBMIN = MAX( 2, ILAENV( 2, 'DORMLQ', SIDE // TRANS, M, N, K, -1 ) )
       END IF
    ELSE
       IWS = NW
    END IF

    IF ( NB < NBMIN .OR. NB >= K ) THEN

       ! Use unblocked code
       CALL DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )

    ELSE

       ! Use blocked code
       IF ( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
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

          CALL DLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ), &
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
          CALL DLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB, &
                       A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK, LDWORK )
       END DO

    END IF

    WORK( 1 ) = REAL( LWKOPT, DP )

  END SUBROUTINE DORMLQ

!==========================================================================
!
!  DGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorithm.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGEQR2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGEQR2 computes a QR factorization of a real m by n matrix A:
! A = Q * R.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] TAU
!
!          TAU is REAL(DP)         array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (N)
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
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
!  Further Details:
!=====================
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
  SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, DLARF, DLARFG

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ONE = 1.0_DP

    ! .. Local Scalars ..
    INTEGER  :: I, K
    REAL(DP) :: AII

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
       CALL XERBLA( 'DGEQR2', -INFO )
       RETURN
    END IF

    K = MIN( M, N )

    DO I = 1, K

       ! Generate elementary reflector H(i) to annihilate A(i+1:m,i)

       CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, TAU( I ) )

       IF ( I < N ) THEN

          ! Apply H(i) to A(i:m,i+1:n) from the left
          AII = A( I, I )
          A( I, I ) = ONE
          CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                      A( I, I+1 ), LDA, WORK )
          A( I, I ) = AII

       END IF

    END DO

  END SUBROUTINE DGEQR2

!==========================================================================
!
!  DORG2R generates all or part of the orthogonal matrix Q from a QR
! factorization determined by sgeqrf (unblocked algorithm).
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORG2R + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2r.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORG2R generates an m by n real matrix Q with orthonormal columns,
! which is defined as the first n columns of a product of k elementary
! reflectors of order m
!
!       Q  =  H(1) H(2) . . . H(k)
!
! as returned by DGEQRF.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (N)
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
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

  SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, DSCAL, DLARF

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: K, LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    INTEGER :: I, J, L

    ! .. Executable Statements ..

    ! Test the input arguments

    INFO = 0
    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 .OR. N > M ) THEN
       INFO = -2
    ELSE IF ( K < 0 .OR. K > N ) THEN
       INFO = -3
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -5
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DORG2R', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( N <= 0 ) RETURN

    ! Initialise columns k+1:n to columns of the unit matrix

    DO J = K + 1, N
       DO L = 1, M
          A( L, J ) = ZERO
       END DO
       A( J, J ) = ONE
    END DO

    DO I = K, 1, -1

       ! Apply H(i) to A(i:m,i:n) from the left
       IF ( I < N ) THEN
          A( I, I ) = ONE
          CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                      A( I, I+1 ), LDA, WORK )
       END IF
       IF ( I < M ) CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
       A( I, I ) = ONE - TAU( I )

       ! Set A(1:i-1,i) to zero
       DO L = 1, I - 1
          A( L, I ) = ZERO
       END DO

    END DO

  END SUBROUTINE DORG2R

!==========================================================================
!
!  DLABRD reduces the first nb rows and columns of a general matrix to a bidiagonal form.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLABRD + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabrd.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,
!                          LDY )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDX, LDY, M, N, NB
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), D( * ), E( * ), TAUP( * ), &
!                          TAUQ( * ), X( LDX,* ), Y( LDY,* )
!       ..
!
!  Purpose:
!  =============
!
! DLABRD reduces the first NB rows and columns of a real general
! m by n matrix A to upper or lower bidiagonal form by an orthogonal
! transformation Q**T * A * P, and returns the matrices X and Y which
! are needed to apply the transformation to the unreduced part of A.
!
! If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
! bidiagonal form.
!
! This is an auxiliary routine called by DGEBRD
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows in the matrix A.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns in the matrix A.
!
! param[in] NB
!
!          NB is INTEGER
!          The number of leading rows and columns of A to be reduced.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the m by n general matrix to be reduced.
!          On exit, the first NB rows and columns of the matrix are
!          overwritten; the rest of the array is unchanged.
!          If m >= n, elements on and below the diagonal in the first NB
!            columns, with the array TAUQ, represent the orthogonal
!            matrix Q as a product of elementary reflectors; and
!            elements above the diagonal in the first NB rows, with the
!            array TAUP, represent the orthogonal matrix P as a product
!            of elementary reflectors.
!          If m < n, elements below the diagonal in the first NB
!            columns, with the array TAUQ, represent the orthogonal
!            matrix Q as a product of elementary reflectors, and
!            elements on and above the diagonal in the first NB rows,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors.
!          See Further Details.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] D
!
!          D is REAL(DP)         array, dimension (NB)
!          The diagonal elements of the first NB rows and columns of
!          the reduced matrix.  D(i) = A(i,i).
!
! param[out] E
!
!          E is REAL(DP)         array, dimension (NB)
!          The off-diagonal elements of the first NB rows and columns of
!          the reduced matrix.
!
! param[out] TAUQ
!
!          TAUQ is REAL(DP)         array dimension (NB)
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix Q. See Further Details.
!
! param[out] TAUP
!
!          TAUP is REAL(DP)         array, dimension (NB)
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix P. See Further Details.
!
! param[out] X
!
!          X is REAL(DP)         array, dimension (LDX,NB)
!          The m-by-nb matrix X required to update the unreduced part
!          of A.
!
! param[in] LDX
!
!          LDX is INTEGER
!          The leading dimension of the array X. LDX >= max(1,M).
!
! param[out] Y
!
!          Y is REAL(DP)         array, dimension (LDY,NB)
!          The n-by-nb matrix Y required to update the unreduced part
!          of A.
!
! param[in] LDY
!
!          LDY is INTEGER
!          The leading dimension of the array Y. LDY >= max(1,N).
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
!  Further Details:
!=====================
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors.
!
!  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
!  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
!  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
!  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
!  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The elements of the vectors v and u together form the m-by-nb matrix
!  V and the nb-by-n matrix U**T which are needed, with X and Y, to apply
!  the transformation to the unreduced part of the matrix, using a block
!  update of the form:  A := A - V*Y**T - X*U**T.
!
!  The contents of A on exit are illustrated by the following examples
!  with nb = 2:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
!    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
!    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
!    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
!    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
!    (  v1  v2  a   a   a  )
!
!  where a denotes an element of the original matrix which is unchanged,
!  vi denotes an element of the vector defining H(i), and ui an element
!  of the vector defining G(i).
!
!=====================================================================

  SUBROUTINE DLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y, LDY )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DSCAL, DGEMV, DLARFG

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: LDA, LDX, LDY, M, N, NB

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), D( * ), E( * ), TAUP( * ), &
                               TAUQ( * ), X( LDX,* ), Y( LDY,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    INTEGER :: I

    ! .. Executable Statements ..

    ! Quick return if possible

    IF ( M <= 0 .OR. N <= 0 ) RETURN

    IF ( M >= N ) THEN

       ! Reduce to upper bidiagonal form

       DO I = 1, NB

          ! Update A(i:m,i)

          CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, A( I, 1 ), &
                      LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )
          CALL DGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ), &
                      LDX, A( 1, I ), 1, ONE, A( I, I ), 1 )

          ! Generate reflection Q(i) to annihilate A(i+1:m,i)

          CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                       TAUQ( I ) )
          D( I ) = A( I, I )

          IF ( I < N ) THEN

             A( I, I ) = ONE

             ! Compute Y(i+1:n,i)

             CALL DGEMV( 'Transpose', M-I+1, N-I, ONE, A( I, I+1 ), &
                         LDA, A( I, I ), 1, ZERO, Y( I+1, I ), 1 )
             CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, A( I, 1 ), LDA, &
                         A( I, I ), 1, ZERO, Y( 1, I ), 1 )
             CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), &
                         LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
             CALL DGEMV( 'Transpose', M-I+1, I-1, ONE, X( I, 1 ), LDX, &
                         A( I, I ), 1, ZERO, Y( 1, I ), 1 )
             CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ), &
                         LDA, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
             CALL DSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )

             ! Update A(i,i+1:n)

             CALL DGEMV( 'No transpose', N-I, I, -ONE, Y( I+1, 1 ), &
                         LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )
             CALL DGEMV( 'Transpose', I-1, N-I, -ONE, A( 1, I+1 ), &
                         LDA, X( I, 1 ), LDX, ONE, A( I, I+1 ), LDA )

             ! Generate reflection P(i) to annihilate A(i,i+2:n)

             CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ), &
                          LDA, TAUP( I ) )
             E( I ) = A( I, I+1 )
             A( I, I+1 ) = ONE

             ! Compute X(i+1:m,i)

             CALL DGEMV( 'No transpose', M-I, N-I, ONE, A( I+1, I+1 ), &
                         LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )
             CALL DGEMV( 'Transpose', N-I, I, ONE, Y( I+1, 1 ), LDY, &
                         A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
             CALL DGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ), &
                         LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
             CALL DGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), &
                         LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
             CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), &
                         LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
             CALL DSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )

          END IF

       END DO

    ELSE

       ! Reduce to lower bidiagonal form

       DO I = 1, NB

          ! Update A(i,i:n)

          CALL DGEMV( 'No transpose', N-I+1, I-1, -ONE, Y( I, 1 ), &
                      LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )
          CALL DGEMV( 'Transpose', I-1, N-I+1, -ONE, A( 1, I ), LDA, &
                      X( I, 1 ), LDX, ONE, A( I, I ), LDA )

          ! Generate reflection P(i) to annihilate A(i,i+1:n)

          CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, &
                       TAUP( I ) )

          D( I ) = A( I, I )

          IF ( I < M ) THEN

             A( I, I ) = ONE

             ! Compute X(i+1:m,i)

             CALL DGEMV( 'No transpose', M-I, N-I+1, ONE, A( I+1, I ), &
                         LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )
             CALL DGEMV( 'Transpose', N-I+1, I-1, ONE, Y( I, 1 ), LDY, &
                         A( I, I ), LDA, ZERO, X( 1, I ), 1 )
             CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ), &
                         LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
             CALL DGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ), &
                         LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )
             CALL DGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ), &
                         LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
             CALL DSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )

             ! Update A(i+1:m,i)

             CALL DGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ), &
                         LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )
             CALL DGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ), &
                         LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 )

             ! Generate reflection Q(i) to annihilate A(i+2:m,i)

             CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1, &
                          TAUQ( I ) )
             E( I ) = A( I+1, I )
             A( I+1, I ) = ONE

             ! Compute Y(i+1:n,i)

             CALL DGEMV( 'Transpose', M-I, N-I, ONE, A( I+1, I+1 ), &
                         LDA, A( I+1, I ), 1, ZERO, Y( I+1, I ), 1 )
             CALL DGEMV( 'Transpose', M-I, I-1, ONE, A( I+1, 1 ), LDA, &
                         A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )
             CALL DGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ), &
                         LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
             CALL DGEMV( 'Transpose', M-I, I, ONE, X( I+1, 1 ), LDX, &
                         A( I+1, I ), 1, ZERO, Y( 1, I ), 1 )
             CALL DGEMV( 'Transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, &
                         Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
             CALL DSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )

          END IF

       END DO

    END IF

  END SUBROUTINE DLABRD

!==========================================================================
!
!  DGEBD2 reduces a general matrix to bidiagonal form using an unblocked algorithm.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DGEBD2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgebd2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), D( * ), E( * ), TAUP( * ), &
!                          TAUQ( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DGEBD2 reduces a real general m by n matrix A to upper or lower
! bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!
! If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
!
!  Arguments:
!  ==========
!
! param[in] M
!
!          M is INTEGER
!          The number of rows in the matrix A.  M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns in the matrix A.  N >= 0.
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          On entry, the m by n general matrix to be reduced.
!          On exit,
!          if m >= n, the diagonal and the first superdiagonal are
!            overwritten with the upper bidiagonal matrix B; the
!            elements below the diagonal, with the array TAUQ, represent
!            the orthogonal matrix Q as a product of elementary
!            reflectors, and the elements above the first superdiagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors;
!          if m < n, the diagonal and the first subdiagonal are
!            overwritten with the lower bidiagonal matrix B; the
!            elements below the first subdiagonal, with the array TAUQ,
!            represent the orthogonal matrix Q as a product of
!            elementary reflectors, and the elements above the diagonal,
!            with the array TAUP, represent the orthogonal matrix P as
!            a product of elementary reflectors.
!          See Further Details.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! param[out] D
!
!          D is REAL(DP)         array, dimension (min(M,N))
!          The diagonal elements of the bidiagonal matrix B:
!          D(i) = A(i,i).
!
! param[out] E
!
!          E is REAL(DP)         array, dimension (min(M,N)-1)
!          The off-diagonal elements of the bidiagonal matrix B:
!          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
!          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
!
! param[out] TAUQ
!
!          TAUQ is REAL(DP)         array dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix Q. See Further Details.
!
! param[out] TAUP
!
!          TAUP is REAL(DP)         array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors which
!          represent the orthogonal matrix P. See Further Details.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (max(M,N))
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
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
!  Further Details:
!=====================
!
!  The matrices Q and P are represented as products of elementary
!  reflectors:
!
!  If m >= n,
!
!     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
!  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  If m < n,
!
!     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
!
!  Each H(i) and G(i) has the form:
!
!     H(i) = I - tauq * v * v**T  and G(i) = I - taup * u * u**T
!
!  where tauq and taup are real scalars, and v and u are real vectors;
!  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
!  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
!  tauq is stored in TAUQ(i) and taup in TAUP(i).
!
!  The contents of A on exit are illustrated by the following examples:
!
!  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
!
!    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
!    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
!    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
!    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
!    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
!    (  v1  v2  v3  v4  v5 )
!
!  where d and e denote diagonal and off-diagonal elements of B, vi
!  denotes an element of the vector defining H(i), and ui an element of
!  the vector defining G(i).
!
!=====================================================================

  SUBROUTINE DGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: XERBLA, DLARF, DLARFG

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: LDA, M, N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), D( * ), E( * ), TAUP( * ), &
                               TAUQ( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP

    ! .. Local Scalars ..
    INTEGER            I

    ! .. Executable Statements ..

    ! Test the input parameters

    INFO = 0
    IF ( M < 0 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = -4
    END IF

    IF ( INFO < 0 ) THEN
       CALL XERBLA( 'DGEBD2', -INFO )
       RETURN
    END IF

    IF ( M >= N ) THEN

       ! Reduce to upper bidiagonal form

       DO I = 1, N

          ! Generate elementary reflector H(i) to annihilate A(i+1:m,i)

          CALL DLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                       TAUQ( I ) )
          D( I ) = A( I, I )
          A( I, I ) = ONE

          ! Apply H(i) to A(i:m,i+1:n) from the left

          IF ( I < N ) &
               CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAUQ( I ), &
                           A( I, I+1 ), LDA, WORK )
          A( I, I ) = D( I )

          IF ( I < N ) THEN

             ! Generate elementary reflector G(i) to annihilate
             ! A(i,i+2:n)

             CALL DLARFG( N-I, A( I, I+1 ), A( I, MIN( I+2, N ) ), LDA, TAUP( I ) )
             E( I ) = A( I, I+1 )
             A( I, I+1 ) = ONE

             ! Apply G(i) to A(i+1:m,i+1:n) from the right

             CALL DLARF( 'Right', M-I, N-I, A( I, I+1 ), LDA, &
                         TAUP( I ), A( I+1, I+1 ), LDA, WORK )
             A( I, I+1 ) = E( I )

          ELSE

             TAUP( I ) = ZERO

          END IF

       END DO

    ELSE

       ! Reduce to lower bidiagonal form

       DO I = 1, M

          ! Generate elementary reflector G(i) to annihilate A(i,i+1:n)

          CALL DLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, TAUP( I ) )
          D( I ) = A( I, I )
          A( I, I ) = ONE

          ! Apply G(i) to A(i+1:m,i:n) from the right

          IF ( I < M ) &
               CALL DLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, &
                           TAUP( I ), A( I+1, I ), LDA, WORK )
          A( I, I ) = D( I )

          IF ( I < M ) THEN

             ! Generate elementary reflector H(i) to annihilate
             ! A(i+2:m,i)

             CALL DLARFG( M-I, A( I+1, I ), A( MIN( I+2, M ), I ), 1, TAUQ( I ) )
             E( I ) = A( I+1, I )
             A( I+1, I ) = ONE

             ! Apply H(i) to A(i+1:m,i+1:n) from the left

             CALL DLARF( 'Left', M-I, N-I, A( I+1, I ), 1, TAUQ( I ), &
                         A( I+1, I+1 ), LDA, WORK )
             A( I+1, I ) = E( I )

          ELSE

             TAUQ( I ) = ZERO

          END IF

       END DO

    END IF

  END SUBROUTINE DGEBD2

!==========================================================================
!
!  DLARTG generates a plane rotation with real cosine and real sine.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLARTG + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartg.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARTG( F, G, CS, SN, R )
!
!       .. Scalar Arguments ..
!       REAL(DP)           CS, F, G, R, SN
!       ..
!
!  Purpose:
!  =============
!
! DLARTG generate a plane rotation so that
!
!    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!    [ -SN  CS  ]     [ G ]     [ 0 ]
!
! This is a slower, more accurate version of the BLAS1 routine DROTG,
! with the following other differences:
!    F and G are unchanged on return.
!    If G=0, then CS=1 and SN=0.
!    If F=0 and (G /= 0), then CS=0 and SN=1 without doing any
!       floating point operations (saves work in DBDSQR when
!       there are zeros on the diagonal).
!
! If F exceeds G in magnitude, CS will be positive.
!
!  Arguments:
!  ==========
!
! param[in] F
!
!          F is REAL(DP)
!          The first component of vector to be rotated.
!
! param[in] G
!
!          G is REAL(DP)
!          The second component of vector to be rotated.
!
! param[out] CS
!
!          CS is REAL(DP)
!          The cosine of the rotation.
!
! param[out] SN
!
!          SN is REAL(DP)
!          The sine of the rotation.
!
! param[out] R
!
!          R is REAL(DP)
!          The nonzero component of the rotated vector.
!
!  This version has a few statements commented out for thread safety
!  (machine parameters are computed on each entry). 10 feb 03, SJH.
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

  SUBROUTINE DLARTG( F, G, CS, SN, R )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN)  :: F, G
    REAL(DP), INTENT(OUT) :: CS, SN, R

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO   = 0.0_DP
    REAL(DP), PARAMETER :: ONE    = 1.0_DP
    REAL(DP), PARAMETER :: TWO    = 2.0_DP
    REAL(DP), PARAMETER :: EPS    = TINY( ZERO )
    REAL(DP), PARAMETER :: SAFMIN = REAL( RADIX(ZERO), DP)**MAX( MINEXPONENT(ZERO)-1, &
                                                               1-MAXEXPONENT(ZERO) )
    REAL(DP), PARAMETER :: SAFMAX = ONE / SAFMIN
    REAL(DP), PARAMETER :: RTMIN  = SQRT( SAFMIN )
    REAL(DP), PARAMETER :: RTMAX  = SQRT( SAFMAX/TWO )

    ! .. Local Scalars ..
    REAL(DP) :: D, F1, FS, G1, GS, U

    ! .. Executable Statements ..

    F1 = ABS( F )
    G1 = ABS( G )
    IF ( ABS(G) < EPS ) THEN
       CS = ONE
       SN = ZERO
       R = F
    ELSE IF ( ABS(F) < EPS ) THEN
       CS = ZERO
       SN = SIGN( ONE, G )
       R = G
    ELSE IF ( F1 > RTMIN .AND. F1 < RTMAX .AND. &
              G1 > RTMIN .AND. G1 < RTMAX ) THEN
       D  = SQRT( F*F + G*G )
       CS = F1 / D
       R  = SIGN( D, F )
       SN = G / R
    ELSE
       U = MIN( SAFMAX, MAX( SAFMIN, F1, G1 ) )
       FS = F / U
       GS = G / U
       D  = SQRT( FS*FS + GS*GS )
       CS = ABS( FS ) / D
       R  = SIGN( D, F )
       SN = GS / R
       R  = R*U
    END IF

  END SUBROUTINE DLARTG

!==========================================================================
!
!  DLASR applies a sequence of plane rotations to a general rectangular matrix.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASR + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasr.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, PIVOT, SIDE
!       INTEGER            LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), C( * ), S( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASR applies a sequence of plane rotations to a real matrix A,
! from either the left or the right.
!
! When SIDE = 'L', the transformation takes the form
!
!    A := P*A
!
! and when SIDE = 'R', the transformation takes the form
!
!    A := A*P**T
!
! where P is an orthogonal matrix consisting of a sequence of z plane
! rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
! and P**T is the transpose of P.
!
! When DIRECT = 'F' (Forward sequence), then
!
!    P = P(z-1) * ... * P(2) * P(1)
!
! and when DIRECT = 'B' (Backward sequence), then
!
!    P = P(1) * P(2) * ... * P(z-1)
!
! where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
!
!    R(k) = (  c(k)  s(k) )
!         = ( -s(k)  c(k) ).
!
! When PIVOT = 'V' (Variable pivot), the rotation is performed
! for the plane (k,k+1), i.e., P(k) has the form
!
!    P(k) = (  1                                            )
!           (       ...                                     )
!           (              1                                )
!           (                   c(k)  s(k)                  )
!           (                  -s(k)  c(k)                  )
!           (                                1              )
!           (                                     ...       )
!           (                                            1  )
!
! where R(k) appears as a rank-2 modification to the identity matrix in
! rows and columns k and k+1.
!
! When PIVOT = 'T' (Top pivot), the rotation is performed for the
! plane (1,k+1), so P(k) has the form
!
!    P(k) = (  c(k)                    s(k)                 )
!           (         1                                     )
!           (              ...                              )
!           (                     1                         )
!           ( -s(k)                    c(k)                 )
!           (                                 1             )
!           (                                      ...      )
!           (                                             1 )
!
! where R(k) appears in rows and columns 1 and k+1.
!
! Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
! performed for the plane (k,z), giving P(k) the form
!
!    P(k) = ( 1                                             )
!           (      ...                                      )
!           (             1                                 )
!           (                  c(k)                    s(k) )
!           (                         1                     )
!           (                              ...              )
!           (                                     1         )
!           (                 -s(k)                    c(k) )
!
! where R(k) appears in rows and columns k and z.  The rotations are
! performed without ever forming P(k) explicitly.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P**T
!
! param[in] PIVOT
!
!          PIVOT is CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
! param[in] DIRECT
!
!          DIRECT is CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
! param[in] C
!
!          C is REAL(DP)         array, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          The cosines c(k) of the plane rotations.
!
! param[in] S
!
!          S is REAL(DP)         array, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          The sines s(k) of the plane rotations.  The 2-by-2 plane
!          rotation part of the matrix P(k), R(k), has the form
!          R(k) = (  c(k)  s(k) )
!                 ( -s(k)  c(k) ).
!
! param[in,out] A
!
!          A is REAL(DP)         array, dimension (LDA,N)
!          The M-by-N matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P**T if SIDE = 'L'.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
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

  SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN) :: DIRECT, PIVOT, SIDE
    INTEGER,   INTENT(IN) :: LDA, M, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), C( * ), S( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS   = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: I, INFO, J
    REAL(DP) :: CTEMP, STEMP, TEMP

    ! .. Executable Statements ..

    ! Test the input parameters

    INFO = 0
    IF ( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
       INFO = 1
    ELSE IF ( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, 'T' ) &
                                         .OR. LSAME( PIVOT, 'B' ) ) ) THEN
       INFO = 2
    ELSE IF ( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) THEN
       INFO = 3
    ELSE IF ( M < 0 ) THEN
       INFO = 4
    ELSE IF ( N < 0 ) THEN
       INFO = 5
    ELSE IF ( LDA < MAX( 1, M ) ) THEN
       INFO = 9
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DLASR ', INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( ( M == 0 ) .OR. ( N == 0 ) ) RETURN

    IF ( LSAME( SIDE, 'L' ) ) THEN

       ! Form  P * A

       IF ( LSAME( PIVOT, 'V' ) ) THEN

          IF ( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, M - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, N
                      TEMP = A( J+1, I )
                      A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                      A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                   END DO
                END IF
             END DO

          ELSE IF ( LSAME( DIRECT, 'B' ) ) THEN

             DO J = M - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, N
                      TEMP = A( J+1, I )
                      A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                      A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
                   END DO
                END IF
             END DO
          END IF

       ELSE IF ( LSAME( PIVOT, 'T' ) ) THEN

          IF ( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 2, M
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                      A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                   END DO
                END IF
             END DO
          ELSE IF ( LSAME( DIRECT, 'B' ) ) THEN
             DO J = M, 2, -1
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                      A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
                   END DO
                END IF
             END DO
          END IF

       ELSE IF ( LSAME( PIVOT, 'B' ) ) THEN

          IF ( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, M - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                      A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                   END DO
                END IF
             END DO
          ELSE IF ( LSAME( DIRECT, 'B' ) ) THEN
             DO J = M - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, N
                      TEMP = A( J, I )
                      A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                      A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
                   END DO
                END IF
             END DO
          END IF

       END IF

    ELSE IF ( LSAME( SIDE, 'R' ) ) THEN

       ! Form A * P**T

       IF ( LSAME( PIVOT, 'V' ) ) THEN

          IF ( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, N - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J+1 )
                      A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                      A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                   END DO
                END IF
             END DO
          ELSE IF ( LSAME( DIRECT, 'B' ) ) THEN
             DO J = N - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J+1 )
                      A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                      A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
                   END DO
                END IF
             END DO
          END IF

       ELSE IF ( LSAME( PIVOT, 'T' ) ) THEN

          IF ( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 2, N
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                      A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                   END DO
                END IF
             END DO
          ELSE IF ( LSAME( DIRECT, 'B' ) ) THEN
             DO J = N, 2, -1
                CTEMP = C( J-1 )
                STEMP = S( J-1 )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                      A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
                   END DO
                END IF
             END DO
          END IF

       ELSE IF ( LSAME( PIVOT, 'B' ) ) THEN

          IF ( LSAME( DIRECT, 'F' ) ) THEN
             DO J = 1, N - 1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                      A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                   END DO
                END IF
             END DO
          ELSE IF ( LSAME( DIRECT, 'B' ) ) THEN
             DO J = N - 1, 1, -1
                CTEMP = C( J )
                STEMP = S( J )
                IF ( (ABS(CTEMP-ONE) >= SPACING(ONE)) .OR. &
                     (ABS(STEMP) >= EPS) ) THEN
                   DO I = 1, M
                      TEMP = A( I, J )
                      A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                      A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
                   END DO
                END IF
             END DO
          END IF

       END IF

    END IF

  END SUBROUTINE DLASR

!==========================================================================
!
!  DLASV2 computes the singular value decomposition of a 2-by-2 triangular matrix.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASV2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasv2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )
!
!       .. Scalar Arguments ..
!       REAL(DP)           CSL, CSR, F, G, H, SNL, SNR, SSMAX, SSMIN
!       ..
!
!  Purpose:
!  =============
!
! DLASV2 computes the singular value decomposition of a 2-by-2
! triangular matrix
!    [  F   G  ]
!    [  0   H  ].
! On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
! smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
! right singular vectors for abs(SSMAX), giving the decomposition
!
!    [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
!    [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
!
!  Arguments:
!  ==========
!
! param[in] F
!
!          F is REAL(DP)
!          The (1,1) element of the 2-by-2 matrix.
!
! param[in] G
!
!          G is REAL(DP)
!          The (1,2) element of the 2-by-2 matrix.
!
! param[in] H
!
!          H is REAL(DP)
!          The (2,2) element of the 2-by-2 matrix.
!
! param[out] SSMIN
!
!          SSMIN is REAL(DP)
!          abs(SSMIN) is the smaller singular value.
!
! param[out] SSMAX
!
!          SSMAX is REAL(DP)
!          abs(SSMAX) is the larger singular value.
!
! param[out] SNL
!
!          SNL is REAL(DP)
!
! param[out] CSL
!
!          CSL is REAL(DP)
!          The vector (CSL, SNL) is a unit left singular vector for the
!          singular value abs(SSMAX).
!
! param[out] SNR
!
!          SNR is REAL(DP)
!
! param[out] CSR
!
!          CSR is REAL(DP)
!          The vector (CSR, SNR) is a unit right singular vector for the
!          singular value abs(SSMAX).
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
!  Further Details:
!=====================
!
!  Any input parameter may be aliased with any output parameter.
!
!  Barring over/underflow and assuming a guard digit in subtraction, all
!  output quantities are correct to within a few units in the last
!  place (ulps).
!
!  In IEEE arithmetic, the code works correctly if one matrix element is
!  infinite.
!
!  Overflow will not occur unless the largest singular value itself
!  overflows or is within a few ulps of overflow. (On machines with
!  partial overflow, like the Cray, overflow may occur if the largest
!  singular value is within a factor of 2 of overflow.)
!
!  Underflow is harmless if underflow is gradual. Otherwise, results
!  may correspond to a matrix modified by perturbations of size near
!  the underflow threshold.
!
!=====================================================================

  SUBROUTINE DLASV2( F, G, H, SSMIN, SSMAX, SNR, CSR, SNL, CSL )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DLAMCH

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN)  :: F, G, H
    REAL(DP), INTENT(OUT) :: CSL, CSR, SNL, SNR, SSMAX, SSMIN

    ! .. Parameters ..

    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: HALF = 0.5_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: TWO  = 2.0_DP
    REAL(DP), PARAMETER :: FOUR = 4.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL  :: GASMAL, SWAP
    INTEGER  :: PMAX
    REAL(DP) :: A, CLT, CRT, D, FA, FT, GA, GT, HA, HT, L, M, MM, R, S, &
                SLT, SRT, T, TEMP, TSIGN, TT

    ! .. Executable Statements ..

    FT = F
    FA = ABS( FT )
    HT = H
    HA = ABS( H )

    ! PMAX points to the maximum absolute element of matrix
    ! PMAX = 1 if F largest in absolute values
    ! PMAX = 2 if G largest in absolute values
    ! PMAX = 3 if H largest in absolute values

    PMAX = 1
    SWAP = ( HA > FA )
    IF ( SWAP ) THEN
       PMAX = 3
       TEMP = FT
       FT = HT
       HT = TEMP
       TEMP = FA
       FA = HA
       HA = TEMP

       ! Now FA >= HA

    END IF
    GT = G
    GA = ABS( GT )

    IF ( GA < EPS ) THEN

       ! Diagonal matrix
       SSMIN = HA
       SSMAX = FA
       CLT = ONE
       CRT = ONE
       SLT = ZERO
       SRT = ZERO

    ELSE

       GASMAL = .TRUE.

       IF ( GA > FA ) THEN
          PMAX = 2
          IF ( FA < GA * DLAMCH( 'EPS' ) ) THEN

             ! Case of very large GA

             GASMAL = .FALSE.
             SSMAX = GA
             IF ( HA > ONE ) THEN
                SSMIN = FA / ( GA / HA )
             ELSE
                SSMIN = ( FA / GA )*HA
             END IF
             CLT = ONE
             SLT = HT / GT
             SRT = ONE
             CRT = FT / GT
          END IF
       END IF

       IF ( GASMAL ) THEN

          ! Normal case

          D = FA - HA

          IF ( FA <= EPS .or. D >= HUGE( ONE ) ) THEN
             ! Copes with infinite F or H
             L = ONE
          ELSE
             L = D / FA
          END IF

          ! Note that 0 <= L <= 1

          M = GT / FT

          ! Note that abs(M) <= 1/macheps

          T = TWO - L

          ! Note that T >= 1

          MM = M*M
          TT = T*T
          S = SQRT( TT+MM )

          ! Note that 1 <= S <= 1 + 1/macheps

          IF ( L < EPS ) THEN
             R = ABS( M )
          ELSE
             R = SQRT( L*L+MM )
          END IF

          ! Note that 0 <= R <= 1 + 1/macheps

          A = HALF*( S+R )

          ! Note that 1 <= A <= 1 + abs(M)

          SSMIN = HA / A
          SSMAX = FA*A

          IF ( MM < EPS ) THEN
             ! Note that M is very tiny
             IF ( L < EPS ) THEN
                T = SIGN( TWO, FT )*SIGN( ONE, GT )
             ELSE
                T = GT / SIGN( D, FT ) + M / T
             END IF
          ELSE
             T = ( M / ( S+T )+M / ( R+L ) )*( ONE+A )
          END IF

          L = SQRT( T*T+FOUR )
          CRT = TWO / L
          SRT = T / L
          CLT = ( CRT+SRT*M ) / A
          SLT = ( HT / FT )*SRT / A

       END IF

    END IF

    IF ( SWAP ) THEN
       CSL = SRT
       SNL = CRT
       CSR = SLT
       SNR = CLT
    ELSE
       CSL = CLT
       SNL = SLT
       CSR = CRT
       SNR = SRT
    END IF

    ! Correct signs of SSMAX and SSMIN

    IF ( PMAX == 1 ) &
         TSIGN = SIGN( ONE, CSR )*SIGN( ONE, CSL )*SIGN( ONE, F )
    IF ( PMAX == 2 ) &
         TSIGN = SIGN( ONE, SNR )*SIGN( ONE, CSL )*SIGN( ONE, G )
    IF ( PMAX == 3 ) &
         TSIGN = SIGN( ONE, SNR )*SIGN( ONE, SNL )*SIGN( ONE, H )
    SSMAX = SIGN( SSMAX, TSIGN )
    SSMIN = SIGN( SSMIN, TSIGN*SIGN( ONE, F )*SIGN( ONE, H ) )

  END SUBROUTINE DLASV2

!==========================================================================
!
!  DROT
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       REAL(DP)         C,S
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL(DP)         DX(*),DY(*)
!       ..
!
!  Purpose:
!  =============
!
!    DROT applies a plane rotation.
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
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!=====================================================================

  SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN) :: C,S
    INTEGER,  INTENT(IN) :: INCX,INCY,N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: DX(*), DY(*)

    ! .. Local Scalars ..
    REAL(DP) :: DTEMP
    INTEGER  :: I,IX,IY

    IF (N <= 0) RETURN

    IF (INCX == 1 .AND. INCY == 1) THEN

       ! code for both increments equal to 1
       DO I = 1,N
          DTEMP = C*DX(I) + S*DY(I)
          DY(I) = C*DY(I) - S*DX(I)
          DX(I) = DTEMP
       END DO

    ELSE

       ! code for unequal increments or equal increments not equal
       ! to 1
       IX = 1
       IY = 1
       IF (INCX < 0) IX = (-N+1)*INCX + 1
       IF (INCY < 0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          DTEMP = C*DX(IX) + S*DY(IY)
          DY(IY) = C*DY(IY) - S*DX(IX)
          DX(IX) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
       END DO
    END IF

  END SUBROUTINE DROT

!==========================================================================
!
!  DLAS2 computes singular values of a 2-by-2 triangular matrix.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLAS2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlas2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )
!
!       .. Scalar Arguments ..
!       REAL(DP)           F, G, H, SSMAX, SSMIN
!       ..
!
!  Purpose:
!  =============
!
! DLAS2  computes the singular values of the 2-by-2 matrix
!    [  F   G  ]
!    [  0   H  ].
! On return, SSMIN is the smaller singular value and SSMAX is the
! larger singular value.
!
!  Arguments:
!  ==========
!
! param[in] F
!
!          F is REAL(DP)
!          The (1,1) element of the 2-by-2 matrix.
!
! param[in] G
!
!          G is REAL(DP)
!          The (1,2) element of the 2-by-2 matrix.
!
! param[in] H
!
!          H is REAL(DP)
!          The (2,2) element of the 2-by-2 matrix.
!
! param[out] SSMIN
!
!          SSMIN is REAL(DP)
!          The smaller singular value.
!
! param[out] SSMAX
!
!          SSMAX is REAL(DP)
!          The larger singular value.
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
!  Further Details:
!=====================
!
!  Barring over/underflow, all output quantities are correct to within
!  a few units in the last place (ulps), even in the absence of a guard
!  digit in addition/subtraction.
!
!  In IEEE arithmetic, the code works correctly if one matrix element is
!  infinite.
!
!  Overflow will not occur unless the largest singular value itself
!  overflows, or is within a few ulps of overflow. (On machines with
!  partial overflow, like the Cray, overflow may occur if the largest
!  singular value is within a factor of 2 of overflow.)
!
!  Underflow is harmless if underflow is gradual. Otherwise, results
!  may correspond to a matrix modified by perturbations of size near
!  the underflow threshold.
!
!=====================================================================

  SUBROUTINE DLAS2( F, G, H, SSMIN, SSMAX )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN)  :: F, G, H
    REAL(DP), INTENT(OUT) :: SSMAX, SSMIN

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: TWO  = 2.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL(DP) :: AS, AT, AU, C, FA, FHMN, FHMX, GA, HA

    ! .. Executable Statements ..

    FA = ABS( F )
    GA = ABS( G )
    HA = ABS( H )
    FHMN = MIN( FA, HA )
    FHMX = MAX( FA, HA )

    IF ( ABS(FHMN) < EPS ) THEN

       SSMIN = ZERO
       IF ( ABS(FHMX) < EPS ) THEN
          SSMAX = GA
       ELSE
          SSMAX = MAX( FHMX, GA ) * &
                  SQRT( ONE + ( MIN( FHMX, GA ) / MAX( FHMX, GA ) )**2 )
       END IF

    ELSE

       IF ( GA < FHMX ) THEN
          AS    = ONE + FHMN / FHMX
          AT    = ( FHMX-FHMN ) / FHMX
          AU    = ( GA / FHMX )**2
          C     = TWO / ( SQRT( AS*AS+AU )+SQRT( AT*AT+AU ) )
          SSMIN = FHMN*C
          SSMAX = FHMX / C
       ELSE
          AU = FHMX / GA
          IF ( ABS(AU) < EPS ) THEN
             ! Avoid possible harmful underflow if exponent range
             ! asymmetric (true SSMIN may not underflow even if
             ! AU underflows)
             SSMIN = ( FHMN*FHMX ) / GA
             SSMAX = GA
          ELSE
             AS    = ONE + FHMN / FHMX
             AT    = ( FHMX-FHMN ) / FHMX
             C     = ONE / ( SQRT( ONE+( AS*AU )**2 ) &
                           + SQRT( ONE+( AT*AU )**2 ) )
             SSMIN = ( FHMN*C )*AU
             SSMIN = SSMIN + SSMIN
             SSMAX = GA / ( C+C )
          END IF
       END IF

    END IF

  END SUBROUTINE DLAS2

!==========================================================================
!
!  DSWAP
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL(DP)         DX(*),DY(*)
!       ..
!
!  Purpose:
!  =============
!
!    interchanges two vectors.
!    uses unrolled loops for increments equal one.
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
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!=====================================================================

  SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: INCX, INCY, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: DX(*), DY(*)

    ! .. Local Scalars ..
    REAL(DP) :: DTEMP
    INTEGER  :: I,IX,IY,M,MP1

    IF (N <= 0) RETURN

    IF (INCX == 1 .AND. INCY == 1) THEN

       ! code for both increments equal to 1

       ! clean-up loop
       M = MOD(N,3)
       IF (M /= 0) THEN
          DO I = 1,M
             DTEMP = DX(I)
             DX(I) = DY(I)
             DY(I) = DTEMP
          END DO
          IF (N < 3) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,3
          DTEMP   = DX(I)
          DX(I)   = DY(I)
          DY(I)   = DTEMP
          DTEMP   = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP   = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
       END DO

    ELSE

       ! code for unequal increments or equal increments not equal to 1

       IX = 1
       IY = 1
       IF (INCX < 0) IX = (-N+1)*INCX + 1
       IF (INCY < 0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
       END DO

    END IF

  END SUBROUTINE DSWAP

!==========================================================================
!
!  DLARF applies an elementary reflector to a general rectangular matrix.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLARF + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarf.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       INTEGER            INCV, LDC, M, N
!       REAL(DP)           TAU
!       ..
!       .. Array Arguments ..
!       REAL(DP)           C( LDC,* ), V( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DLARF applies a real elementary reflector H to a real m by n matrix
! C, from either the left or the right. H is represented in the form
!
!       H = I - tau * v * v**T
!
! where tau is a real scalar and v is a real vector.
!
! If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C.
!
! param[in] V
!
!          V is REAL(DP)         array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
! param[in] INCV
!
!          INCV is INTEGER
!          The increment between elements of v. INCV <> 0.
!
! param[in] TAU
!
!          TAU is REAL(DP)
!          The value tau in the representation of H.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
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

  SUBROUTINE DLARF( SIDEC, M, N, V, INCV, TAU, C, LDC, WORK )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, DGEMV, DGER, ILADLR, ILADLC

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: SIDEC
    INTEGER,      INTENT(IN) :: INCV, LDC, M, N
    REAL(DP),     INTENT(IN) :: TAU

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: C( LDC,* ), V( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

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
       DO WHILE( LASTV > 0 .AND. ABS( V(I) ) < EPS )
          LASTV = LASTV - 1
          I = I - INCV
       END DO
       IF ( APPLYLEFT ) THEN
          ! Scan for the last non-zero column in C(1:lastv,:).
          LASTC = ILADLC(LASTV, N, C, LDC)
       ELSE
          ! Scan for the last non-zero row in C(:,1:lastv).
          LASTC = ILADLR(M, LASTV, C, LDC)
       END IF

    END IF

    ! Note that lastc == 0 renders the BLAS operations null; no special
    ! case is needed at this level.

    IF ( APPLYLEFT ) THEN

       ! Form  H * C
       IF ( LASTV > 0 ) THEN

          ! w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
          CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
                      ZERO, WORK, 1 )

          ! C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
          CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )

       END IF

    ELSE

       ! Form  C * H
       IF ( LASTV > 0 ) THEN

          ! w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
          CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, &
                      ZERO, WORK, 1 )

          ! C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
          CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
       END IF

    END IF

  END SUBROUTINE DLARF

!==========================================================================
!
!  DCOPY
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL(DP)         DX(*),DY(*)
!       ..
!
!  Purpose:
!  =============
!
!    DCOPY copies a vector, x, to a vector, y.
!    uses unrolled loops for increments equal to one.
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
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!=====================================================================

  SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: INCX, INCY, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: DX(*), DY(*)

    ! .. Local Scalars ..
    INTEGER :: I,IX,IY,M,MP1

    IF (N <= 0) RETURN

    IF (INCX == 1 .AND. INCY == 1) THEN

       ! code for both increments equal to 1

       ! clean-up loop
       M = MOD(N,7)
       IF (M /= 0) THEN
          DO I = 1,M
             DY(I) = DX(I)
          END DO
          IF (N < 7) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,7
          DY(I)   = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
       END DO

    ELSE

       ! code for unequal increments or equal increments not equal to 1

       IX = 1
       IY = 1
       IF (INCX < 0) IX = (-N+1)*INCX + 1
       IF (INCY < 0) IY = (-N+1)*INCY + 1
       DO I = 1,N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
       END DO

    END IF

  END SUBROUTINE DCOPY

!==========================================================================
!
!  DTRMM
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       REAL(DP)         ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       REAL(DP)         A(LDA,*),B(LDB,*)
!       ..
!
!  Purpose:
!  =============
!
! DTRMM  performs one of the matrix-matrix operations
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
!
!          SIDE is CHARACTER*1
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
! param[in] UPLO
!
!          UPLO is CHARACTER*1
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
! param[in] TRANSA
!
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
!
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
!
!          M is INTEGER
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!
! param[in] N
!
!          N is INTEGER
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!
! param[in] ALPHA
!
!          ALPHA is REAL(DP)        .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!
! param[in] A
!
!           A is REAL(DP)         array of DIMENSION ( LDA, k ), where k is m
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
!
!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!
! param[in,out] B
!
!          B is REAL(DP)         array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
! param[in] LDB
!
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

  SUBROUTINE DTRMM(SIDES,UPLOS,TRANSAS,DIAGS,M,N,ALPHA,A,LDA,B,LDB)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level3 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL(DP),     INTENT(IN) :: ALPHA
    INTEGER,      INTENT(IN) :: LDA,LDB,M,N
    CHARACTER(*), INTENT(IN) :: DIAGS,SIDES,TRANSAS,UPLOS

    ! .. Array Arguments ..
    REAL(DP),  INTENT(INOUT) :: A(LDA,*), B(LDB,*)

    ! .. Local Scalars ..
    REAL(DP)  :: TEMP
    INTEGER   :: I,INFO,J,K,NROWA
    LOGICAL   :: LSIDE,NOUNIT,UPPER
    CHARACTER :: DIAG,SIDE,TRANSA,UPLO

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

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
       CALL XERBLA('DTRMM ',INFO)
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

  END SUBROUTINE DTRMM

!==========================================================================
!
!  DGEMV
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       REAL(DP)         ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       REAL(DP)         A(LDA,*),X(*),Y(*)
!       ..
!
!  Purpose:
!  =============
!
! DGEMV  performs one of the matrix-vector operations
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
!
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
!
!          M is INTEGER
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!
! param[in] N
!
!          N is INTEGER
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!
! param[in] ALPHA
!
!          ALPHA is REAL(DP)        .
!           On entry, ALPHA specifies the scalar alpha.
!
! param[in] A
!
!          A is REAL(DP)         array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!
! param[in] LDA
!
!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!
! param[in] X
!
!          X is REAL(DP)         array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!
! param[in] INCX
!
!          INCX is INTEGER
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!
! param[in] BETA
!
!          BETA is REAL(DP)        .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!
! param[in,out] Y
!
!          Y is REAL(DP)         array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
! param[in] INCY
!
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

  SUBROUTINE DGEMV(TRANSC,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level2 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL(DP),     INTENT(IN) :: ALPHA,BETA
    INTEGER,      INTENT(IN) :: INCX,INCY,LDA,M,N
    CHARACTER(*), INTENT(IN) :: TRANSC

    ! .. Array Arguments ..
    REAL(DP),  INTENT(INOUT) :: A(LDA,*),X(*),Y(*)

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL(DP)  :: TEMP
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
       CALL XERBLA('DGEMV ',INFO)
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

  END SUBROUTINE DGEMV

!==========================================================================
!
!  DTRMV
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       REAL(DP)         A(LDA,*),X(*)
!       ..
!
!  Purpose:
!  =============
!
! DTRMV  performs one of the matrix-vector operations
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
!
!          UPLO is CHARACTER*1
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
! param[in] TRANS
!
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
!
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
!
!          N is INTEGER
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!
! param[in] A
!
!          A is REAL(DP)         array of DIMENSION ( LDA, n ).
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
!
!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!
! param[in,out] X
!
!          X is REAL(DP)         array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
! param[in] INCX
!
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

  SUBROUTINE DTRMV(UPLOC,TRANSC,DIAGC,N,A,LDA,X,INCX)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level2 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER,      INTENT(IN) :: INCX, LDA, N
    CHARACTER(*), INTENT(IN) :: DIAGC, TRANSC, UPLOC

    ! .. Array Arguments ..
    REAL(DP),  INTENT(INOUT) :: A(LDA,*), X(*)

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL(DP)  :: TEMP
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
       CALL XERBLA('DTRMV ',INFO)
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

  END SUBROUTINE DTRMV

!==========================================================================
!
!  DLARFG generates an elementary reflector (Householder matrix).
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLARFG + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL(DP)           ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       REAL(DP)           X( * )
!       ..
!
!  Purpose:
!  =============
!
! DLARFG generates a real elementary reflector H of order n, such
! that
!
!       H * ( alpha ) = ( beta ),   H**T * H = I.
!           (   x   )   (   0  )
!
! where alpha and beta are scalars, and x is an (n-1)-element real
! vector. H is represented in the form
!
!       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!                     ( v )
!
! where tau is a real scalar and v is a real (n-1)-element
! vector.
!
! If the elements of x are all zero, then tau = 0 and H is taken to be
! the unit matrix.
!
! Otherwise  1 <= tau <= 2.
!
!  Arguments:
!  ==========
!
! param[in] N
!
!          N is INTEGER
!          The order of the elementary reflector.
!
! param[in,out] ALPHA
!
!          ALPHA is REAL(DP)
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
! param[in,out] X
!
!          X is REAL(DP)         array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
! param[in] INCX
!
!          INCX is INTEGER
!          The increment between elements of X. INCX > 0.
!
! param[out] TAU
!
!          TAU is REAL(DP)
!          The value tau.
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

  SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DSCAL, DLAMCH, DLAPY2, DNRM2

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER,  INTENT(IN)    :: INCX, N
    REAL(DP), INTENT(INOUT) :: ALPHA, TAU

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: X( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: J, KNT
    REAL(DP) :: BETA, RSAFMN, SAFMIN, XNORM

    ! .. Executable Statements ..

    IF ( N <= 1 ) THEN
       TAU = ZERO
       RETURN
    END IF

    XNORM = DNRM2( N-1, X, INCX )

    IF ( XNORM < EPS ) THEN

       TAU = ZERO

    ELSE

       ! general case

       BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
       SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
       KNT = 0

       IF ( ABS( BETA ) < SAFMIN ) THEN

          ! XNORM, BETA may be inaccurate; scale X and recompute them

          RSAFMN = ONE / SAFMIN

          DO WHILE ( ABS( BETA ) < SAFMIN )
             KNT = KNT + 1
             CALL DSCAL( N-1, RSAFMN, X, INCX )
             BETA = BETA*RSAFMN
             ALPHA = ALPHA*RSAFMN
          END DO

          ! New BETA is at most 1, at least SAFMIN

          XNORM = DNRM2( N-1, X, INCX )
          BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )

       END IF

       TAU = ( BETA-ALPHA ) / BETA
       CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )

       ! If ALPHA is subnormal, it may lose relative accuracy

       DO J = 1, KNT
          BETA = BETA*SAFMIN
       END DO
       ALPHA = BETA

    END IF

  END SUBROUTINE DLARFG

!==========================================================================
!
!  DORM2R multiplies a general matrix by the orthogonal matrix from a
! QR factorization determined by sgeqrf (unblocked algorithm).
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORM2R + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORM2R overwrites the general real m by n matrix C with
!
!       Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!
!       C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!       C * Q**T if SIDE = 'R' and TRANS = 'T',
!
! where Q is a real orthogonal matrix defined as the product of k
! elementary reflectors
!
!       Q = H(1) H(2) . . . H(k)
!
! as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
! if SIDE = 'R'.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': apply Q or Q**T from the Left
!          = 'R': apply Q or Q**T from the Right
!
! param[in] TRANS
!
!          TRANS is CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q**T (Transpose)
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C. N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
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

  SUBROUTINE DORM2R( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, DLARF

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP),   INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ONE = 1.0_DP

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, NOTRAN
    INTEGER   :: I, I1, I2, I3, IC, JC, MI, NI, NQ
    REAL(DP)  :: AII
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
       CALL XERBLA( 'DORM2R', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( M == 0 .OR. N == 0 .OR. K == 0 ) RETURN

    IF ( (LEFT .AND. .NOT.NOTRAN) .OR. (.NOT.LEFT .AND. NOTRAN) ) THEN
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
       CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), LDC, WORK )
       A( I, I ) = AII

    END DO

  END SUBROUTINE DORM2R

!==========================================================================
!
!  DORML2 multiplies a general matrix by the orthogonal matrix from a
! LQ factorization determined by sgelqf (unblocked algorithm).
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DORML2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorml2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DORML2 overwrites the general real m by n matrix C with
!
!       Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!
!       C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!       C * Q**T if SIDE = 'R' and TRANS = 'T',
!
! where Q is a real orthogonal matrix defined as the product of k
! elementary reflectors
!
!       Q = H(k) . . . H(2) H(1)
!
! as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n
! if SIDE = 'R'.
!
!  Arguments:
!  ==========
!
! param[in] SIDE
!
!          SIDE is CHARACTER*1
!          = 'L': apply Q or Q**T from the Left
!          = 'R': apply Q or Q**T from the Right
!
! param[in] TRANS
!
!          TRANS is CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q**T (Transpose)
!
! param[in] M
!
!          M is INTEGER
!          The number of rows of the matrix C. M >= 0.
!
! param[in] N
!
!          N is INTEGER
!          The number of columns of the matrix C. N >= 0.
!
! param[in] K
!
!          K is INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
! param[in] A
!
!          A is REAL(DP)         array, dimension
!                               (LDA,M) if SIDE = 'L',
!                               (LDA,N) if SIDE = 'R'
!          The i-th row must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGELQF in the first k rows of its array argument A.
!          A is modified by the routine but restored on exit.
!
! param[in] LDA
!
!          LDA is INTEGER
!          The leading dimension of the array A. LDA >= max(1,K).
!
! param[in] TAU
!
!          TAU is REAL(DP)         array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGELQF.
!
! param[in,out] C
!
!          C is REAL(DP)         array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
! param[in] LDC
!
!          LDC is INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
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

  SUBROUTINE DORML2( SIDEC, TRANSC, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, DLARF

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN)  :: SIDEC, TRANSC
    INTEGER,      INTENT(IN)  :: K, LDA, LDC, M, N
    INTEGER,      INTENT(OUT) :: INFO

    ! .. array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), C( LDC,* ), TAU( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ONE = 1.0_DP

    ! .. Local Scalars ..
    LOGICAL   :: LEFT, NOTRAN
    INTEGER   :: I, I1, I2, I3, IC, JC, MI, NI, NQ
    REAL(DP)  :: AII
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
       CALL XERBLA( 'DORML2', -INFO )
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
       CALL DLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ), C( IC, JC ), LDC, WORK )
       A( I, I ) = AII

    END DO

  END SUBROUTINE DORML2

!==========================================================================
!
!  DLASRT sorts numbers in increasing or decreasing order.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASRT + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasrt.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASRT( ID, N, D, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          ID
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           D( * )
!       ..
!
!  Purpose:
!  =============
!
! Sort the numbers in D in increasing order (if ID = 'I') or
! in decreasing order (if ID = 'D' ).
!
! Use Quick Sort, reverting to Insertion sort on arrays of
! size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments:
!  ==========
!
! param[in] ID
!
!          ID is CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
! param[in] N
!
!          N is INTEGER
!          The length of the array D.
!
! param[in,out] D
!
!          D is REAL(DP)         array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
! param[out] INFO
!
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
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

  SUBROUTINE DLASRT( ID, N, D, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: ID
    INTEGER,   INTENT(IN)  :: N
    INTEGER,   INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: D( * )

    ! .. Parameters ..
    INTEGER, PARAMETER :: SELECT = 20

    ! .. Local Scalars ..
    INTEGER  :: DIR, ENDD, I, J, START, STKPNT
    REAL(DP) :: D1, D2, D3, DMNMX, TMP

    ! .. Local Arrays ..
    INTEGER  :: STACK( 2, 32 )

    ! .. Executable Statements ..

    ! Test the input paramters.

    INFO = 0
    DIR = -1
    IF ( LSAME( ID, 'D' ) ) THEN
       DIR = 0
    ELSE IF ( LSAME( ID, 'I' ) ) THEN
       DIR = 1
    END IF

    IF ( DIR == -1 ) THEN
       INFO = -1
    ELSE IF ( N < 0 ) THEN
       INFO = -2
    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DLASRT', -INFO )
       RETURN
    END IF

    ! Quick return if possible

    IF ( N <= 1 ) RETURN

    STKPNT = 1
    STACK( 1, 1 ) = 1
    STACK( 2, 1 ) = N

    DO WHILE (STKPNT > 0)

       START = STACK( 1, STKPNT )
       ENDD = STACK( 2, STKPNT )
       STKPNT = STKPNT - 1

       IF ( ENDD-START <= SELECT .AND. ENDD-START > 0 ) THEN

          ! Do Insertion sort on D( START:ENDD )

          IF ( DIR == 0 ) THEN

             ! Sort into decreasing order

             DO I = START + 1, ENDD
                DO J = I, START + 1, -1
                   IF ( D( J ) > D( J-1 ) ) THEN
                      DMNMX = D( J )
                      D( J ) = D( J-1 )
                      D( J-1 ) = DMNMX
                   ELSE
                      EXIT
                   END IF
                END DO
             END DO

          ELSE

             ! Sort into increasing order

             DO I = START + 1, ENDD
                DO J = I, START + 1, -1
                   IF ( D( J ) < D( J-1 ) ) THEN
                      DMNMX = D( J )
                      D( J ) = D( J-1 )
                      D( J-1 ) = DMNMX
                   ELSE
                      EXIT
                   END IF
                END DO
             END DO

          END IF

       ELSE IF ( ENDD-START > SELECT ) THEN

          ! Partition D( START:ENDD ) and stack parts, largest one first

          ! Choose partition entry as median of 3

          D1 = D( START )
          D2 = D( ENDD )
          I = ( START+ENDD ) / 2
          D3 = D( I )
          IF ( D1 < D2 ) THEN
             IF ( D3 < D1 ) THEN
                DMNMX = D1
             ELSE IF ( D3 < D2 ) THEN
                DMNMX = D3
             ELSE
                DMNMX = D2
             END IF
          ELSE
             IF ( D3 < D2 ) THEN
                DMNMX = D2
             ELSE IF ( D3 < D1 ) THEN
                DMNMX = D3
             ELSE
                DMNMX = D1
             END IF
          END IF

          IF ( DIR == 0 ) THEN

             ! Sort into decreasing order

             I = START
             J = ENDD
             DO WHILE (.TRUE.)
                DO WHILE ( D( J ) < DMNMX )
                   J = J - 1
                END DO
                DO WHILE ( D( I ) > DMNMX )
                   I = I + 1
                END DO
                IF ( I < J ) THEN
                   TMP = D( I )
                   D( I ) = D( J )
                   D( J ) = TMP
                ELSE
                   EXIT
                END IF
             END DO

             IF ( J-START > ENDD-J-1 ) THEN
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = START
                STACK( 2, STKPNT ) = J
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = J + 1
                STACK( 2, STKPNT ) = ENDD
             ELSE
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = J + 1
                STACK( 2, STKPNT ) = ENDD
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = START
                STACK( 2, STKPNT ) = J
             END IF

          ELSE

             ! Sort into increasing order

             I = START
             J = ENDD
             DO WHILE (.TRUE.)
                DO WHILE ( D( J ) > DMNMX )
                   J = J - 1
                END DO
                DO WHILE ( D( I ) < DMNMX )
                   I = I + 1
                END DO
                IF ( I < J ) THEN
                   TMP = D( I )
                   D( I ) = D( J )
                   D( J ) = TMP
                ELSE
                   EXIT
                END IF
             END DO

             IF ( J-START > ENDD-J-1 ) THEN
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = START
                STACK( 2, STKPNT ) = J
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = J + 1
                STACK( 2, STKPNT ) = ENDD
             ELSE
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = J + 1
                STACK( 2, STKPNT ) = ENDD
                STKPNT = STKPNT + 1
                STACK( 1, STKPNT ) = START
                STACK( 2, STKPNT ) = J
             END IF
          END IF
       END IF

    END DO

  END SUBROUTINE DLASRT

!==========================================================================
!
!  DGER
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       REAL(DP)         ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       REAL(DP)         A(LDA,*),X(*),Y(*)
!       ..
!
!  Purpose:
!  =============
!
! DGER   performs the rank 1 operation
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
!
!          M is INTEGER
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!
! param[in] N
!
!          N is INTEGER
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!
! param[in] ALPHA
!
!          ALPHA is REAL(DP)        .
!           On entry, ALPHA specifies the scalar alpha.
!
! param[in] X
!
!          X is REAL(DP)         array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!
! param[in] INCX
!
!          INCX is INTEGER
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!
! param[in] Y
!
!          Y is REAL(DP)         array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!
! param[in] INCY
!
!          INCY is INTEGER
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!
! param[in,out] A
!
!          A is REAL(DP)         array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
! param[in] LDA
!
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

  SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level2 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN) :: ALPHA
    INTEGER,  INTENT(IN) :: INCX,INCY,LDA,M,N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A(LDA,*),X(*),Y(*)

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL(DP) :: TEMP
    INTEGER  :: I,INFO,IX,J,JY,KX

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
       CALL XERBLA('DGER  ',INFO)
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

  END SUBROUTINE DGER

!==========================================================================
!
!  DNRM2
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL(DP)         FUNCTION DNRM2(N,X,INCX)
!
!        IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL(DP)         X(*)
!       ..
!
!  Purpose:
!  =============
!
! DNRM2 returns the euclidean norm of a vector via the function
! name, so that
!
!    DNRM2 := sqrt( x'*x )
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
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.
!
!=====================================================================

  REAL(DP) FUNCTION DNRM2(N,X,INCX)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) ::  INCX,N

    ! .. Array Arguments ..
    REAL(DP), INTENT(IN) :: X(*)

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL(DP) :: ABSXI,NORM,SCALE,SSQ
    INTEGER  :: IX

    IF (N < 1 .OR. INCX < 1) THEN

       NORM = ZERO

    ELSE IF (N == 1) THEN

       NORM = ABS(X(1))

    ELSE

       SCALE = ZERO
       SSQ = ONE

       ! The following loop is equivalent to this call to the LAPACK
       ! auxiliary routine:
       ! CALL DLASSQ( N, X, INCX, SCALE, SSQ )

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

    DNRM2 = NORM

  END FUNCTION DNRM2

!==========================================================================
!
!  DLAPY2 returns sqrt(x2+y2).
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLAPY2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f
!
!  Definition:
!  ===========
!
!       REAL(DP)         FUNCTION DLAPY2( X, Y )
!
!        IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       REAL(DP)           X, Y
!       ..
!
!  Purpose:
!  =============
!
! DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
! overflow.
!
!  Arguments:
!  ==========
!
! param[in] X
!
!          X is REAL(DP)
!
! param[in] Y
!
!          Y is REAL(DP)
!          X and Y specify the values x and y.
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

  REAL(DP) FUNCTION DLAPY2( X, Y )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN) :: X, Y

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    REAL(DP) :: W, XABS, YABS, Z

    ! .. Executable Statements ..

    XABS = ABS( X )
    YABS = ABS( Y )
    W = MAX( XABS, YABS )
    Z = MIN( XABS, YABS )
    IF ( Z < EPS ) THEN
       DLAPY2 = W
    ELSE
       DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
    END IF

  END FUNCTION DLAPY2

!==========================================================================
!
!  DLASQ1 computes the singular values of a real square bidiagonal matrix. Used by sbdsqr.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASQ1 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq1.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ1( N, D, E, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           D( * ), E( * ), WORK( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASQ1 computes the singular values of a real N-by-N bidiagonal
! matrix with diagonal D and off-diagonal E. The singular values
! are computed to high relative accuracy, in the absence of
! denormalization, underflow and overflow. The algorithm was first
! presented in
!
! "Accurate singular values and differential qd algorithms" by K. V.
! Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
! 1994,
!
! and the present implementation is described in "An implementation of
! the dqds Algorithm (Positive Case)", LAPACK Working Note.
!
!  Arguments:
!  ==========
!
! param[in] N
!
!          N is INTEGER
!        The number of rows and columns in the matrix. N >= 0.
!
! param[in,out] D
!
!          D is REAL(DP)         array, dimension (N)
!        On entry, D contains the diagonal elements of the
!        bidiagonal matrix whose SVD is desired. On normal exit,
!        D contains the singular values in decreasing order.
!
! param[in,out] E
!
!          E is REAL(DP)         array, dimension (N)
!        On entry, elements E(1:N-1) contain the off-diagonal elements
!        of the bidiagonal matrix whose SVD is desired.
!        On exit, E is overwritten.
!
! param[out] WORK
!
!          WORK is REAL(DP)         array, dimension (4*N)
!
! param[out] INFO
!
!          INFO is INTEGER
!        = 0: successful exit
!        < 0: if INFO = -i, the i-th argument had an illegal value
!        > 0: the algorithm failed
!             = 1, a split was marked by a positive value in E
!             = 2, current block of Z not diagonalized after 100*N
!                  iterations (in inner while loop)  On exit D and E
!                  represent a matrix with the same singular values
!                  which the calling subroutine could use to finish the
!                  computation, or even feed back into DLASQ1
!             = 3, termination criterion of outer while loop not met
!                  (program created more than N unreduced blocks)
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

  SUBROUTINE DLASQ1( N, D, E, WORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DLAMCH, DCOPY, DLAS2, DLASCL, DLASQ2, DLASRT, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: D( * ), E( * ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: I, IINFO
    REAL(DP) :: SMALL, SCALE, SAFMIN, SIGMN, SIGMX

    ! .. Executable Statements ..

    INFO = 0
    IF ( N < 0 ) THEN
       INFO = -2
       CALL XERBLA( 'DLASQ1', -INFO )
       RETURN
    ELSE IF ( N == 0 ) THEN
       RETURN
    ELSE IF ( N == 1 ) THEN
       D( 1 ) = ABS( D( 1 ) )
       RETURN
    ELSE IF ( N == 2 ) THEN
       CALL DLAS2( D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX )
       D( 1 ) = SIGMX
       D( 2 ) = SIGMN
       RETURN
    END IF

    ! Estimate the largest singular value.

    SIGMX = ZERO
    DO I = 1, N - 1
       D( I ) = ABS( D( I ) )
       SIGMX = MAX( SIGMX, ABS( E( I ) ) )
    END DO
    D( N ) = ABS( D( N ) )

    ! Early return if SIGMX is zero (matrix is already diagonal).

    IF ( ABS(SIGMX) < EPS ) THEN
       CALL DLASRT( 'D', N, D, IINFO )
       RETURN
    END IF

    DO I = 1, N
       SIGMX = MAX( SIGMX, D( I ) )
    END DO

    ! Copy D and E into WORK (in the Z format) and scale (squaring the
    ! input data makes scaling by a power of the radix pointless).

    SMALL = DLAMCH( 'Precision' )
    SAFMIN = DLAMCH( 'Safe minimum' )
    SCALE = SQRT( SMALL / SAFMIN )
    CALL DCOPY( N, D, 1, WORK( 1 ), 2 )
    CALL DCOPY( N-1, E, 1, WORK( 2 ), 2 )
    CALL DLASCL( 'G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, IINFO )

    ! Compute the q's and e's.

    DO I = 1, 2*N - 1
       WORK( I ) = WORK( I )**2
    END DO
    WORK( 2*N ) = ZERO

    CALL DLASQ2( N, WORK, INFO )

    IF ( INFO == 0 ) THEN

       DO I = 1, N
          D( I ) = SQRT( WORK( I ) )
       END DO
       CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )

    ELSE IF ( INFO == 2 ) THEN

       ! Maximum number of iterations exceeded.  Move data from WORK
       ! into D and E so the calling subroutine can try to finish

       DO I = 1, N
          D( I ) = SQRT( WORK( 2*I-1 ) )
          E( I ) = SQRT( WORK( 2*I ) )
       END DO
       CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO )
       CALL DLASCL( 'G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO )

    END IF

  END SUBROUTINE DLASQ1

!==========================================================================
!
!  DLASQ2 computes all the eigenvalues of the symmetric positive definite
! tridiagonal matrix associated with the qd Array Z to high relative accuracy.
! Used by sbdsqr and sstegr.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASQ2 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq2.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ2( N, Z, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL(DP)           Z( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASQ2 computes all the eigenvalues of the symmetric positive
! definite tridiagonal matrix associated with the qd array Z to high
! relative accuracy are computed to high relative accuracy, in the
! absence of denormalization, underflow and overflow.
!
! To see the relation of Z to the tridiagonal matrix, let L be a
! unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
! let U be an upper bidiagonal matrix with 1's above and diagonal
! Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
! symmetric tridiagonal to which it is similar.
!
! Note : DLASQ2 defines a logical variable, IEEE, which is true
! on machines which follow ieee-754 floating-point standard in their
! handling of infinities and NaNs, and false otherwise. This variable
! is passed to DLASQ3.
!
!  Arguments:
!  ==========
!
! param[in] N
!
!          N is INTEGER
!        The number of rows and columns in the matrix. N >= 0.
!
! param[in,out] Z
!
!          Z is REAL(DP)         array, dimension ( 4*N )
!        On entry Z holds the qd array. On exit, entries 1 to N hold
!        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
!        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If
!        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )
!        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of
!        shifts that failed.
!
! param[out] INFO
!
!          INFO is INTEGER
!        = 0: successful exit
!        < 0: if the i-th argument is a scalar and had an illegal
!             value, then INFO = -i, if the i-th argument is an
!             array and the j-entry had an illegal value, then
!             INFO = -(i*100+j)
!        > 0: the algorithm failed
!              = 1, a split was marked by a positive value in E
!              = 2, current block of Z not diagonalized after 100*N
!                   iterations (in inner while loop).  On exit Z holds
!                   a qd array with the same eigenvalues as the given Z.
!              = 3, termination criterion of outer while loop not met
!                   (program created more than N unreduced blocks)
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
!  Further Details:
!=====================
!
!  Local Variables: I0:N0 defines a current unreduced segment of Z.
!  The shifts are accumulated in SIGMA. Iteration count is in ITER.
!  Ping-pong is controlled by PP (alternates between 0 and 1).
!
!=====================================================================

  SUBROUTINE DLASQ2( N, Z, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: XERBLA, DLASQ3, DLASRT, DLAMCH, ILAENV

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN)  :: N
    INTEGER, INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: Z( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: CBIAS  = 1.5_DP
    REAL(DP), PARAMETER :: ZERO   = 0.0_DP
    REAL(DP), PARAMETER :: HALF   = 0.5_DP
    REAL(DP), PARAMETER :: ONE    = 1.0_DP
    REAL(DP), PARAMETER :: TWO    = 2.0_DP
    REAL(DP), PARAMETER :: FOUR   = 4.0_DP
    REAL(DP), PARAMETER :: HUNDRD = 100.0_DP
    REAL(DP), PARAMETER :: EPS    = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL  :: IEEE
    INTEGER  :: I0, I1, I4, IINFO, IPN4, ITER, IWHILA, IWHILB, &
                K, KMIN, N0, N1, NBIG, NDIV, NFAIL, PP, SPLT, TTYPE
    REAL(DP) :: D, DEE, DEEMIN, DESIG, DMIN, DMIN1, DMIN2, DN,   &
                DN1, DN2, E, EMAX, EMIN, SMALL, G, OLDEMN, QMAX, &
                QMIN, S, SAFMIN, SIGMA, T, TAU, TEMP, TOL, TOL2, &
                TRACE, ZMAX, TEMPE, TEMPQ

    ! .. Executable Statements ..

    ! Test the input arguments.
    ! (in case DLASQ2 is not called by DLASQ1)

    INFO = 0
    SMALL = DLAMCH( 'Precision' )
    SAFMIN = DLAMCH( 'Safe minimum' )
    TOL = SMALL*HUNDRD
    TOL2 = TOL**2

    IF ( N < 0 ) THEN

       INFO = -1
       CALL XERBLA( 'DLASQ2', 1 )
       RETURN

    ELSE IF ( N == 0 ) THEN

       RETURN

    ELSE IF ( N == 1 ) THEN

       ! 1-by-1 case.
       IF ( Z( 1 ) < ZERO ) THEN
          INFO = -201
          CALL XERBLA( 'DLASQ2', 2 )
       END IF
       RETURN

    ELSE IF ( N == 2 ) THEN

       ! 2-by-2 case.
       IF ( Z( 2 ) < ZERO .OR. Z( 3 ) < ZERO ) THEN
          INFO = -2
          CALL XERBLA( 'DLASQ2', 2 )
          RETURN
       ELSE IF ( Z( 3 ) > Z( 1 ) ) THEN
          D = Z( 3 )
          Z( 3 ) = Z( 1 )
          Z( 1 ) = D
       END IF
       Z( 5 ) = Z( 1 ) + Z( 2 ) + Z( 3 )
       IF ( Z( 2 ) > Z( 3 )*TOL2 ) THEN
          T = HALF*( ( Z( 1 )-Z( 3 ) )+Z( 2 ) )
          S = Z( 3 )*( Z( 2 ) / T )
          IF ( S <= T ) THEN
             S = Z( 3 )*( Z( 2 ) / ( T*( ONE+SQRT( ONE+S / T ) ) ) )
          ELSE
             S = Z( 3 )*( Z( 2 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
          END IF
          T = Z( 1 ) + ( S+Z( 2 ) )
          Z( 3 ) = Z( 3 )*( Z( 1 ) / T )
          Z( 1 ) = T
       END IF
       Z( 2 ) = Z( 3 )
       Z( 6 ) = Z( 2 ) + Z( 1 )
       RETURN

    END IF

    ! Check for negative data and compute sums of q's and e's.

    Z( 2*N ) = ZERO
    EMIN = Z( 2 )
    QMAX = ZERO
    ZMAX = ZERO
    D = ZERO
    E = ZERO

    DO K = 1, 2*( N-1 ), 2
       IF ( Z( K ) < ZERO ) THEN
          INFO = -( 200+K )
          CALL XERBLA( 'DLASQ2', 2 )
          RETURN
       ELSE IF ( Z( K+1 ) < ZERO ) THEN
          INFO = -( 200+K+1 )
          CALL XERBLA( 'DLASQ2', 2 )
          RETURN
       END IF
       D = D + Z( K )
       E = E + Z( K+1 )
       QMAX = MAX( QMAX, Z( K ) )
       EMIN = MIN( EMIN, Z( K+1 ) )
       ZMAX = MAX( QMAX, ZMAX, Z( K+1 ) )
    END DO

    IF ( Z( 2*N-1 ) < ZERO ) THEN
       INFO = -( 200+2*N-1 )
       CALL XERBLA( 'DLASQ2', 2 )
       RETURN
    END IF

    D = D + Z( 2*N-1 )
    QMAX = MAX( QMAX, Z( 2*N-1 ) )
    ZMAX = MAX( QMAX, ZMAX )

    ! Check for diagonality.

    IF ( ABS(E) < EPS ) THEN
       DO K = 2, N
          Z( K ) = Z( 2*K-1 )
       END DO
       CALL DLASRT( 'D', N, Z, IINFO )
       Z( 2*N-1 ) = D
       RETURN
    END IF

    TRACE = D + E

    ! Check for zero data.

    IF ( ABS(TRACE) < EPS ) THEN
       Z( 2*N-1 ) = ZERO
       RETURN
    END IF

    ! Check whether the machine is IEEE conformable.

    IEEE = ILAENV( 10, 'DLASQ2', 'N', 1, 2, 3, 4 ) == 1 .AND. &
           ILAENV( 11, 'DLASQ2', 'N', 1, 2, 3, 4 ) == 1

    ! Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).

    DO K = 2*N, 2, -2
       Z( 2*K ) = ZERO
       Z( 2*K-1 ) = Z( K )
       Z( 2*K-2 ) = ZERO
       Z( 2*K-3 ) = Z( K-1 )
    END DO

    I0 = 1
    N0 = N

    ! Reverse the qd-array, if warranted.

    IF ( CBIAS*Z( 4*I0-3 ) < Z( 4*N0-3 ) ) THEN
       IPN4 = 4*( I0+N0 )
       DO I4 = 4*I0, 2*( I0+N0-1 ), 4
          TEMP = Z( I4-3 )
          Z( I4-3 ) = Z( IPN4-I4-3 )
          Z( IPN4-I4-3 ) = TEMP
          TEMP = Z( I4-1 )
          Z( I4-1 ) = Z( IPN4-I4-5 )
          Z( IPN4-I4-5 ) = TEMP
       END DO
    END IF

    ! Initial split checking via dqd and Li's test.

    PP = 0

    DO K = 1, 2

       D = Z( 4*N0+PP-3 )
       DO I4 = 4*( N0-1 ) + PP, 4*I0 + PP, -4
          IF ( Z( I4-1 ) <= TOL2*D ) THEN
             Z( I4-1 ) = -ZERO
             D = Z( I4-3 )
          ELSE
             D = Z( I4-3 )*( D / ( D+Z( I4-1 ) ) )
          END IF
       END DO

       ! dqd maps Z to ZZ plus Li's test.

       EMIN = Z( 4*I0+PP+1 )
       D = Z( 4*I0+PP-3 )
       DO I4 = 4*I0 + PP, 4*( N0-1 ) + PP, 4
          Z( I4-2*PP-2 ) = D + Z( I4-1 )
          IF ( Z( I4-1 ) <= TOL2*D ) THEN
             Z( I4-1 ) = -ZERO
             Z( I4-2*PP-2 ) = D
             Z( I4-2*PP ) = ZERO
             D = Z( I4+1 )
          ELSE IF ( SAFMIN*Z( I4+1 ) < Z( I4-2*PP-2 ) .AND. &
                    SAFMIN*Z( I4-2*PP-2 ) < Z( I4+1 ) ) THEN
             TEMP = Z( I4+1 ) / Z( I4-2*PP-2 )
             Z( I4-2*PP ) = Z( I4-1 )*TEMP
             D = D*TEMP
          ELSE
             Z( I4-2*PP ) = Z( I4+1 )*( Z( I4-1 ) / Z( I4-2*PP-2 ) )
             D = Z( I4+1 )*( D / Z( I4-2*PP-2 ) )
          END IF
          EMIN = MIN( EMIN, Z( I4-2*PP ) )
       END DO
       Z( 4*N0-PP-2 ) = D

       ! Now find qmax.

       QMAX = Z( 4*I0-PP-2 )
       DO I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
          QMAX = MAX( QMAX, Z( I4 ) )
       END DO

       ! Prepare for the next iteration on K.

       PP = 1 - PP
    END DO

    ! Initialise variables to pass to DLASQ3.

    TTYPE = 0
    DMIN1 = ZERO
    DMIN2 = ZERO
    DN    = ZERO
    DN1   = ZERO
    DN2   = ZERO
    G     = ZERO
    TAU   = ZERO

    ITER = 2
    NFAIL = 0
    NDIV = 2*( N0-I0 )

    ALOOP: DO IWHILA = 1, N + 1
       IF ( N0 < 1 ) EXIT ALOOP

       ! While array unfinished do

       ! E(N0) holds the value of SIGMA when submatrix in I0:N0
       ! splits from the rest of the array, but is negated.

       DESIG = ZERO
       IF ( N0 == N ) THEN
          SIGMA = ZERO
       ELSE
          SIGMA = -Z( 4*N0-1 )
       END IF
       IF ( SIGMA < ZERO ) THEN
          INFO = 1
          RETURN
       END IF

       ! Find last unreduced submatrix's top index I0, find QMAX and
       ! EMIN. Find Gershgorin-type bound if Q's much greater than E's.

       EMAX = ZERO
       IF ( N0 > I0 ) THEN
          EMIN = ABS( Z( 4*N0-5 ) )
       ELSE
          EMIN = ZERO
       END IF
       QMIN = Z( 4*N0-3 )
       QMAX = QMIN
       DO I4 = 4*N0, 8, -4
          IF ( Z( I4-5 ) <= ZERO ) EXIT
          IF ( QMIN >= FOUR*EMAX ) THEN
             QMIN = MIN( QMIN, Z( I4-3 ) )
             EMAX = MAX( EMAX, Z( I4-5 ) )
          END IF
          QMAX = MAX( QMAX, Z( I4-7 )+Z( I4-5 ) )
          EMIN = MIN( EMIN, Z( I4-5 ) )
       END DO

       I0 = I4 / 4
       PP = 0

       IF ( N0-I0 > 1 ) THEN
          DEE = Z( 4*I0-3 )
          DEEMIN = DEE
          KMIN = I0
          DO I4 = 4*I0+1, 4*N0-3, 4
             DEE = Z( I4 )*( DEE /( DEE+Z( I4-2 ) ) )
             IF ( DEE <= DEEMIN ) THEN
                DEEMIN = DEE
                KMIN = ( I4+3 )/4
             END IF
          END DO
          IF ( (KMIN-I0)*2 < N0-KMIN .AND. DEEMIN <= HALF*Z(4*N0-3) ) THEN
             IPN4 = 4*( I0+N0 )
             PP = 2
             DO I4 = 4*I0, 2*( I0+N0-1 ), 4
                TEMP = Z( I4-3 )
                Z( I4-3 ) = Z( IPN4-I4-3 )
                Z( IPN4-I4-3 ) = TEMP
                TEMP = Z( I4-2 )
                Z( I4-2 ) = Z( IPN4-I4-2 )
                Z( IPN4-I4-2 ) = TEMP
                TEMP = Z( I4-1 )
                Z( I4-1 ) = Z( IPN4-I4-5 )
                Z( IPN4-I4-5 ) = TEMP
                TEMP = Z( I4 )
                Z( I4 ) = Z( IPN4-I4-4 )
                Z( IPN4-I4-4 ) = TEMP
             END DO
          END IF
       END IF

       ! Put -(initial shift) into DMIN.

       DMIN = -MAX( ZERO, QMIN-TWO*SQRT( QMIN )*SQRT( EMAX ) )

       ! Now I0:N0 is unreduced.
       ! PP = 0 for ping, PP = 1 for pong.
       ! PP = 2 indicates that flipping was applied to the Z array and
       ! and that the tests for deflation upon entry in DLASQ3
       ! should not be performed.

       NBIG = 100*( N0-I0+1 )
       BLOOP: DO IWHILB = 1, NBIG
          IF ( I0 > N0 ) CYCLE ALOOP

          ! While submatrix unfinished take a good dqds step.

          CALL DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
                       ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
                       DN2, G, TAU )

          PP = 1 - PP

          ! When EMIN is very small check for splits.

          IF ( PP == 0 .AND. N0-I0 >= 3 ) THEN
             IF ( Z( 4*N0 )   <= TOL2*QMAX .OR. &
                  Z( 4*N0-1 ) <= TOL2*SIGMA ) THEN
                SPLT = I0 - 1
                QMAX = Z( 4*I0-3 )
                EMIN = Z( 4*I0-1 )
                OLDEMN = Z( 4*I0 )
                DO I4 = 4*I0, 4*( N0-3 ), 4
                   IF ( Z( I4 )   <= TOL2*Z( I4-3 ) .OR. &
                        Z( I4-1 ) <= TOL2*SIGMA ) THEN
                      Z( I4-1 ) = -SIGMA
                      SPLT = I4 / 4
                      QMAX = ZERO
                      EMIN = Z( I4+3 )
                      OLDEMN = Z( I4+4 )
                   ELSE
                      QMAX = MAX( QMAX, Z( I4+1 ) )
                      EMIN = MIN( EMIN, Z( I4-1 ) )
                      OLDEMN = MIN( OLDEMN, Z( I4 ) )
                   END IF
                END DO
                Z( 4*N0-1 ) = EMIN
                Z( 4*N0 ) = OLDEMN
                I0 = SPLT + 1
             END IF
          END IF

       END DO BLOOP

       INFO = 2

       ! Maximum number of iterations exceeded, restore the shift
       ! SIGMA and place the new d's and e's in a qd array.
       ! This might need to be done for several blocks

       I1 = I0
       N1 = N0

       DO WHILE (.TRUE.)

          TEMPQ = Z( 4*I0-3 )
          Z( 4*I0-3 ) = Z( 4*I0-3 ) + SIGMA
          DO K = I0+1, N0
             TEMPE = Z( 4*K-5 )
             Z( 4*K-5 ) = Z( 4*K-5 ) * (TEMPQ / Z( 4*K-7 ))
             TEMPQ = Z( 4*K-3 )
             Z( 4*K-3 ) = Z( 4*K-3 ) + SIGMA + TEMPE - Z( 4*K-5 )
          END DO

          ! Prepare to do this on the previous block if there is one

          IF ( I1 > 1 ) THEN
             N1 = I1-1
             DO WHILE( ( I1 >= 2 ) .AND. ( Z(4*I1-5) >= ZERO ) )
                I1 = I1 - 1
             END DO
             SIGMA = -Z(4*N1-1)
          ELSE
             EXIT
          END IF

       END DO

       DO K = 1, N
          Z( 2*K-1 ) = Z( 4*K-3 )

          ! Only the block 1..N0 is unfinished.  The rest of the e's
          ! must be essentially zero, although sometimes other data
          ! has been stored in them.

          IF ( K < N0 ) THEN
             Z( 2*K ) = Z( 4*K-1 )
          ELSE
             Z( 2*K ) = 0
          END IF
       END DO

       RETURN

    END DO ALOOP

    IF (IWHILA > N + 1) THEN
       INFO = 3
       RETURN
    ENDIF

    ! Move q's to the front.

    DO K = 2, N
       Z( K ) = Z( 4*K-3 )
    END DO

    ! Sort and compute sum of eigenvalues.

    CALL DLASRT( 'D', N, Z, IINFO )

    E = ZERO
    DO K = N, 1, -1
       E = E + Z( K )
    END DO

    ! Store trace, sum(eigenvalues) and information on performance.

    Z( 2*N+1 ) = TRACE
    Z( 2*N+2 ) = E
    Z( 2*N+3 ) = REAL( ITER,DP )
    Z( 2*N+4 ) = REAL( NDIV,DP ) / REAL( N**2,DP )
    Z( 2*N+5 ) = HUNDRD * REAL( NFAIL,DP ) / REAL( ITER,DP )

  END SUBROUTINE DLASQ2

!==========================================================================
!
!  DLASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASQ3 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq3.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
!                          ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,
!                          DN2, G, TAU )
!
!       .. Scalar Arguments ..
!       LOGICAL            IEEE
!       INTEGER            I0, ITER, N0, NDIV, NFAIL, PP
!       REAL(DP)           DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, &
!                          QMAX, SIGMA, TAU
!       ..
!       .. Array Arguments ..
!       REAL(DP)           Z( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
! In case of failure it changes shifts, and tries again until output
! is positive.
!
!  Arguments:
!  ==========
!
! param[in] I0
!
!          I0 is INTEGER
!         First index.
!
! param[in,out] N0
!
!          N0 is INTEGER
!         Last index.
!
! param[in] Z
!
!          Z is REAL(DP)         array, dimension ( 4*N )
!         Z holds the qd array.
!
! param[in,out] PP
!
!          PP is INTEGER
!         PP=0 for ping, PP=1 for pong.
!         PP=2 indicates that flipping was applied to the Z array
!         and that the initial tests for deflation should not be
!         performed.
!
! param[out] DMIN
!
!          DMIN is REAL(DP)
!         Minimum value of d.
!
! param[out] SIGMA
!
!          SIGMA is REAL(DP)
!         Sum of shifts used in current segment.
!
! param[in,out] DESIG
!
!          DESIG is REAL(DP)
!         Lower order part of SIGMA
!
! param[in] QMAX
!
!          QMAX is REAL(DP)
!         Maximum value of q.
!
! param[out] NFAIL
!
!          NFAIL is INTEGER
!         Number of times shift was too big.
!
! param[out] ITER
!
!          ITER is INTEGER
!         Number of iterations.
!
! param[out] NDIV
!
!          NDIV is INTEGER
!         Number of divisions.
!
! param[in] IEEE
!
!          IEEE is LOGICAL
!         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).
!
! param[in,out] TTYPE
!
!          TTYPE is INTEGER
!         Shift type.
!
! param[in,out] DMIN1
!
!          DMIN1 is REAL(DP)
!
! param[in,out] DMIN2
!
!          DMIN2 is REAL(DP)
!
! param[in,out] DN
!
!          DN is REAL(DP)
!
! param[in,out] DN1
!
!          DN1 is REAL(DP)
!
! param[in,out] DN2
!
!          DN2 is REAL(DP)
!
! param[in,out] G
!
!          G is REAL(DP)
!
! param[in,out] TAU
!
!          TAU is REAL(DP)
!
!         These are passed as arguments in order to save their values
!         between calls to DLASQ3.
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

  SUBROUTINE DLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL, &
                     ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1, &
                     DN2, G, TAU )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DLASQ4, DLASQ5, DLASQ6, DISNAN, DLAMCH

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    LOGICAL,  INTENT(IN)    :: IEEE
    INTEGER,  INTENT(INOUT) :: I0, ITER, N0, NDIV, NFAIL, PP
    REAL(DP), INTENT(INOUT) :: DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, &
                               QMAX, SIGMA, TAU

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: Z( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: CBIAS  = 1.5_DP
    REAL(DP), PARAMETER :: ZERO   = 0.0_DP
    REAL(DP), PARAMETER :: QURTR  = 0.25_DP
    REAL(DP), PARAMETER :: HALF   = 0.5_DP
    REAL(DP), PARAMETER :: ONE    = 1.0_DP
    REAL(DP), PARAMETER :: TWO    = 2.0_DP
    REAL(DP), PARAMETER :: HUNDRD = 100.0_DP
    REAL(DP), PARAMETER :: EPS    = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: IPN4, J4, N0IN, NN, TTYPE
    REAL(DP) :: SMALL, S, T, TEMP, TOL, TOL2

    ! .. Executable Statements ..

    N0IN = N0
    SMALL = DLAMCH( 'Precision' )
    TOL = SMALL*HUNDRD
    TOL2 = TOL**2

    ! Check for deflation.

10  CONTINUE

    IF ( N0 < I0 ) RETURN
    IF ( N0 == I0 ) GO TO 20
    NN = 4*N0 + PP
    IF ( N0 == I0+1 ) GO TO 40

    ! Check whether E(N0-1) is negligible, 1 eigenvalue.

    IF ( Z( NN-5 ) > TOL2*( SIGMA+Z( NN-3 ) ) .AND. &
         Z( NN-2*PP-4 ) > TOL2*Z( NN-7 ) ) GO TO 30

20  CONTINUE

    Z( 4*N0-3 ) = Z( 4*N0+PP-3 ) + SIGMA
    N0 = N0 - 1
    GO TO 10

    ! Check  whether E(N0-2) is negligible, 2 eigenvalues.

30  CONTINUE

    IF ( Z( NN-9 ) > TOL2*SIGMA .AND. &
         Z( NN-2*PP-8 ) > TOL2*Z( NN-11 ) ) GO TO 50

40  CONTINUE

    IF ( Z( NN-3 ) > Z( NN-7 ) ) THEN
       S = Z( NN-3 )
       Z( NN-3 ) = Z( NN-7 )
       Z( NN-7 ) = S
    END IF
    T = HALF*( ( Z( NN-7 )-Z( NN-3 ) )+Z( NN-5 ) )
    IF ( Z( NN-5 ) > Z( NN-3 )*TOL2 .AND. ABS(T) > EPS ) THEN
       S = Z( NN-3 )*( Z( NN-5 ) / T )
       IF ( S <= T ) THEN
          S = Z( NN-3 )*( Z( NN-5 ) / &
               ( T*( ONE+SQRT( ONE+S / T ) ) ) )
       ELSE
          S = Z( NN-3 )*( Z( NN-5 ) / ( T+SQRT( T )*SQRT( T+S ) ) )
       END IF
       T = Z( NN-7 ) + ( S+Z( NN-5 ) )
       Z( NN-3 ) = Z( NN-3 )*( Z( NN-7 ) / T )
       Z( NN-7 ) = T
    END IF
    Z( 4*N0-7 ) = Z( NN-7 ) + SIGMA
    Z( 4*N0-3 ) = Z( NN-3 ) + SIGMA
    N0 = N0 - 2
    GO TO 10

50  CONTINUE

    IF ( PP == 2 ) PP = 0

    ! Reverse the qd-array, if warranted.

    IF ( DMIN <= ZERO .OR. N0 < N0IN ) THEN
       IF ( CBIAS*Z( 4*I0+PP-3 ) < Z( 4*N0+PP-3 ) ) THEN
          IPN4 = 4*( I0+N0 )
          DO J4 = 4*I0, 2*( I0+N0-1 ), 4
             TEMP = Z( J4-3 )
             Z( J4-3 ) = Z( IPN4-J4-3 )
             Z( IPN4-J4-3 ) = TEMP
             TEMP = Z( J4-2 )
             Z( J4-2 ) = Z( IPN4-J4-2 )
             Z( IPN4-J4-2 ) = TEMP
             TEMP = Z( J4-1 )
             Z( J4-1 ) = Z( IPN4-J4-5 )
             Z( IPN4-J4-5 ) = TEMP
             TEMP = Z( J4 )
             Z( J4 ) = Z( IPN4-J4-4 )
             Z( IPN4-J4-4 ) = TEMP
          END DO
          IF ( N0-I0 <= 4 ) THEN
             Z( 4*N0+PP-1 ) = Z( 4*I0+PP-1 )
             Z( 4*N0-PP ) = Z( 4*I0-PP )
          END IF
          DMIN2 = MIN( DMIN2, Z( 4*N0+PP-1 ) )
          Z( 4*N0+PP-1 ) = MIN( Z( 4*N0+PP-1 ), Z( 4*I0+PP-1 ), &
               Z( 4*I0+PP+3 ) )
          Z( 4*N0-PP ) = MIN( Z( 4*N0-PP ), Z( 4*I0-PP ), &
               Z( 4*I0-PP+4 ) )
          QMAX = MAX( QMAX, Z( 4*I0+PP-3 ), Z( 4*I0+PP+1 ) )
          DMIN = -ZERO
       END IF
    END IF

    ! Choose a shift.

    CALL DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1, &
         DN2, TAU, TTYPE, G )

    ! Call dqds until DMIN > 0.

70  CONTINUE

    CALL DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, &
         DN1, DN2, IEEE, SMALL )

    NDIV = NDIV + ( N0-I0+2 )
    ITER = ITER + 1

    ! Check status.

    IF ( DMIN >= ZERO .AND. DMIN1 >= ZERO ) THEN

       ! Success.

       GO TO 90

    ELSE IF ( DMIN < ZERO .AND. DMIN1 > ZERO .AND.  &
         Z( 4*( N0-1 )-PP ) < TOL*( SIGMA+DN1 ) .AND. &
         ABS( DN ) < TOL*SIGMA ) THEN

       ! Convergence hidden by negative DN.

       Z( 4*( N0-1 )-PP+2 ) = ZERO
       DMIN = ZERO
       GO TO 90

    ELSE IF ( DMIN < ZERO ) THEN

       ! TAU too big. Select new TAU and try again.

       NFAIL = NFAIL + 1
       IF ( TTYPE < -22 ) THEN

          ! Failed twice. Play it safe.

          TAU = ZERO
       ELSE IF ( DMIN1 > ZERO ) THEN

          ! Late failure. Gives excellent shift.

          TAU = ( TAU+DMIN )*( ONE-TWO*SMALL )
          TTYPE = TTYPE - 11
       ELSE

          ! Early failure. Divide by 4.

          TAU = QURTR*TAU
          TTYPE = TTYPE - 12
       END IF
       GO TO 70

    ELSE IF ( DISNAN( DMIN ) ) THEN

       ! NaN.

       IF ( ABS(TAU) < EPS ) THEN
          GO TO 80
       ELSE
          TAU = ZERO
          GO TO 70
       END IF

    ELSE

       ! Possible underflow. Play it safe.

       GO TO 80

    END IF

    ! Risk of underflow.

80  CONTINUE

    CALL DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2 )
    NDIV = NDIV + ( N0-I0+2 )
    ITER = ITER + 1
    TAU = ZERO

90  CONTINUE

    IF ( TAU < SIGMA ) THEN
       DESIG = DESIG + TAU
       T = SIGMA + DESIG
       DESIG = DESIG - ( T-SIGMA )
    ELSE
       T = SIGMA + TAU
       DESIG = SIGMA - ( T-TAU ) + DESIG
    END IF
    SIGMA = T

  END SUBROUTINE DLASQ3

!==========================================================================
!
!  DLASQ4 computes an approximation to the smallest eigenvalue using values
! of d from the previous transform. Used by sbdsqr.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASQ4 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq4.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,
!                          DN1, DN2, TAU, TTYPE, G )
!
!       .. Scalar Arguments ..
!       INTEGER            I0, N0, N0IN, PP, TTYPE
!       REAL(DP)           DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU
!       ..
!       .. Array Arguments ..
!       REAL(DP)           Z( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASQ4 computes an approximation TAU to the smallest eigenvalue
! using values of d from the previous transform.
!
!  Arguments:
!  ==========
!
! param[in] I0
!
!          I0 is INTEGER
!        First index.
!
! param[in] N0
!
!          N0 is INTEGER
!        Last index.
!
! param[in] Z
!
!          Z is REAL(DP)         array, dimension ( 4*N )
!        Z holds the qd array.
!
! param[in] PP
!
!          PP is INTEGER
!        PP=0 for ping, PP=1 for pong.
!
! param[in] N0IN
!
!          N0IN is INTEGER
!        The value of N0 at start of EIGTEST.
!
! param[in] DMIN
!
!          DMIN is REAL(DP)
!        Minimum value of d.
!
! param[in] DMIN1
!
!          DMIN1 is REAL(DP)
!        Minimum value of d, excluding D( N0 ).
!
! param[in] DMIN2
!
!          DMIN2 is REAL(DP)
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
! param[in] DN
!
!          DN is REAL(DP)
!        d(N)
!
! param[in] DN1
!
!          DN1 is REAL(DP)
!        d(N-1)
!
! param[in] DN2
!
!          DN2 is REAL(DP)
!        d(N-2)
!
! param[out] TAU
!
!          TAU is REAL(DP)
!        This is the shift.
!
! param[out] TTYPE
!
!          TTYPE is INTEGER
!        Shift type.
!
! param[in,out] G
!
!          G is REAL
!        G is passed as an argument in order to save its value between
!        calls to DLASQ4.
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
!  Further Details:
!=====================
!
!  CNST1 = 9/16
!
!=====================================================================

  SUBROUTINE DLASQ4( I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, &
                     DN1, DN2, TAU, TTYPE, G )

    use IEEE_ARITHMETIC, ONLY: IEEE_QUIET_EQ
    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER,  INTENT(INOUT) :: I0, N0, N0IN, PP, TTYPE
    REAL(DP), INTENT(INOUT) :: DMIN, DMIN1, DMIN2, DN, DN1, DN2, G, TAU

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: Z( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: CNST1  = 0.563_DP
    REAL(DP), PARAMETER :: CNST2  = 1.010_DP
    REAL(DP), PARAMETER :: CNST3  = 1.050_DP
    REAL(DP), PARAMETER :: ZERO   = 0.0_DP
    REAL(DP), PARAMETER :: QURTR  = 0.25_DP
    REAL(DP), PARAMETER :: HALF   = 0.5_DP
    REAL(DP), PARAMETER :: ONE    = 1.0_DP
    REAL(DP), PARAMETER :: TWO    = 2.0_DP
    REAL(DP), PARAMETER :: THREE  = 3.0_DP
    REAL(DP), PARAMETER :: HUNDRD = 100.0_DP
    REAL(DP), PARAMETER :: THIRD  = ONE / THREE
    REAL(DP), PARAMETER :: EPS    = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: I4, NN, NP
    REAL(DP) :: A2, B1, B2, GAM, GAP1, GAP2, S

    ! .. Executable Statements ..

    ! A negative DMIN forces the shift to take that absolute value
    ! TTYPE records the type of shift.

    IF ( DMIN <= ZERO ) THEN
       TAU = -DMIN
       TTYPE = -1
       RETURN
    END IF

    NN = 4*N0 + PP

    IF ( N0IN == N0 ) THEN

       ! No eigenvalues deflated.

       IF ( IEEE_QUIET_EQ(DMIN,DN) .OR. IEEE_QUIET_EQ(DMIN,DN1) ) THEN

          B1 = SQRT( Z( NN-3 ) )*SQRT( Z( NN-5 ) )
          B2 = SQRT( Z( NN-7 ) )*SQRT( Z( NN-9 ) )
          A2 = Z( NN-7 ) + Z( NN-5 )

          ! Cases 2 and 3.

          IF ( IEEE_QUIET_EQ(DMIN,DN) .AND. IEEE_QUIET_EQ(DMIN1,DN1) ) THEN

             GAP2 = DMIN2 - A2 - DMIN2*QURTR
             IF ( GAP2 > ZERO .AND. GAP2 > B2 ) THEN
                GAP1 = A2 - DN - ( B2 / GAP2 )*B2
             ELSE
                GAP1 = A2 - DN - ( B1+B2 )
             END IF
             IF ( GAP1 > ZERO .AND. GAP1 > B1 ) THEN
                S = MAX( DN-( B1 / GAP1 )*B1, HALF*DMIN )
                TTYPE = -2
             ELSE
                S = ZERO
                IF ( DN > B1 ) S = DN - B1
                IF ( A2 > ( B1+B2 ) ) S = MIN( S, A2-( B1+B2 ) )
                S = MAX( S, THIRD*DMIN )
                TTYPE = -3
             END IF

          ELSE

             ! Case 4.

             TTYPE = -4
             S = QURTR*DMIN
             IF ( IEEE_QUIET_EQ(DMIN,DN) ) THEN
                GAM = DN
                A2 = ZERO
                IF ( Z( NN-5 ) > Z( NN-7 ) ) RETURN
                B2 = Z( NN-5 ) / Z( NN-7 )
                NP = NN - 9
             ELSE
                NP = NN - 2*PP
                B2 = Z( NP-2 )
                GAM = DN1
                IF ( Z( NP-4 ) > Z( NP-2 ) ) RETURN
                A2 = Z( NP-4 ) / Z( NP-2 )
                IF ( Z( NN-9 ) > Z( NN-11 ) ) RETURN
                B2 = Z( NN-9 ) / Z( NN-11 )
                NP = NN - 13
             END IF

             ! Approximate contribution to norm squared from I < NN-1.

             A2 = A2 + B2
             DO I4 = NP, 4*I0 - 1 + PP, -4
                IF ( ABS(B2) < EPS ) EXIT
                B1 = B2
                IF ( Z( I4 ) > Z( I4-2 ) ) RETURN
                B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                A2 = A2 + B2
                IF ( HUNDRD*MAX( B2, B1 ) < A2 .OR. CNST1 < A2 ) EXIT
             END DO
             A2 = CNST3*A2

             ! Rayleigh quotient residual bound.

             IF ( A2 < CNST1 ) S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )

          END IF

       ELSE IF ( IEEE_QUIET_EQ(DMIN,DN2) ) THEN

          ! Case 5.

          TTYPE = -5
          S = QURTR*DMIN

          ! Compute contribution to norm squared from I > NN-2.

          NP = NN - 2*PP
          B1 = Z( NP-2 )
          B2 = Z( NP-6 )
          GAM = DN2
          IF ( Z( NP-8 ) > B2 .OR. Z( NP-4 ) > B1 ) RETURN
          A2 = ( Z( NP-8 ) / B2 )*( ONE+Z( NP-4 ) / B1 )

          ! Approximate contribution to norm squared from I < NN-2.

          IF ( N0-I0 > 2 ) THEN
             B2 = Z( NN-13 ) / Z( NN-15 )
             A2 = A2 + B2
             DO I4 = NN - 17, 4*I0 - 1 + PP, -4
                IF ( ABS(B2) < EPS ) EXIT
                B1 = B2
                IF ( Z( I4 ) > Z( I4-2 ) ) RETURN
                B2 = B2*( Z( I4 ) / Z( I4-2 ) )
                A2 = A2 + B2
                IF ( HUNDRD*MAX( B2, B1 ) < A2 .OR. CNST1 < A2 ) EXIT
             END DO
             A2 = CNST3*A2
          END IF

          IF ( A2 < CNST1 ) S = GAM*( ONE-SQRT( A2 ) ) / ( ONE+A2 )

       ELSE

          ! Case 6, no information to guide us.

          IF ( TTYPE == -6 ) THEN
             G = G + THIRD*( ONE-G )
          ELSE IF ( TTYPE == -18 ) THEN
             G = QURTR*THIRD
          ELSE
             G = QURTR
          END IF
          S = G*DMIN
          TTYPE = -6

       END IF

    ELSE IF ( N0IN == ( N0+1 ) ) THEN

       ! One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.

       IF ( IEEE_QUIET_EQ(DMIN1,DN1) .AND. IEEE_QUIET_EQ(DMIN2,DN2) ) THEN

          ! Cases 7 and 8.

          TTYPE = -7
          S = THIRD*DMIN1
          IF ( Z( NN-5 ) > Z( NN-7 ) ) RETURN
          B1 = Z( NN-5 ) / Z( NN-7 )
          B2 = B1
          IF ( ABS(B2) >= EPS ) THEN
             DO I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
                A2 = B1
                IF ( Z( I4 ) > Z( I4-2 ) ) RETURN
                B1 = B1*( Z( I4 ) / Z( I4-2 ) )
                B2 = B2 + B1
                IF ( HUNDRD*MAX( B1, A2 ) < B2 ) EXIT
             END DO
          END IF
          B2 = SQRT( CNST3*B2 )
          A2 = DMIN1 / ( ONE+B2**2 )
          GAP2 = HALF*DMIN2 - A2
          IF ( GAP2 > ZERO .AND. GAP2 > B2*A2 ) THEN
             S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
          ELSE
             S = MAX( S, A2*( ONE-CNST2*B2 ) )
             TTYPE = -8
          END IF

       ELSE

          ! Case 9.

          S = QURTR*DMIN1
          IF ( IEEE_QUIET_EQ(DMIN1,DN1) ) S = HALF*DMIN1
          TTYPE = -9

       END IF

    ELSE IF ( N0IN == ( N0+2 ) ) THEN

       ! Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.

       ! Cases 10 and 11.

       IF ( IEEE_QUIET_EQ(DMIN2,DN2) .AND. TWO*Z( NN-5 ) < Z( NN-7 ) ) THEN
          TTYPE = -10
          S = THIRD*DMIN2
          IF ( Z( NN-5 ) > Z( NN-7 ) ) RETURN
          B1 = Z( NN-5 ) / Z( NN-7 )
          B2 = B1
          IF ( ABS(B2) >= EPS ) THEN
             DO I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
                IF ( Z( I4 ) > Z( I4-2 ) ) RETURN
                B1 = B1*( Z( I4 ) / Z( I4-2 ) )
                B2 = B2 + B1
                IF ( HUNDRD*B1 < B2 ) EXIT
             END DO
          ENDIF
          B2 = SQRT( CNST3*B2 )
          A2 = DMIN2 / ( ONE+B2**2 )
          GAP2 = Z( NN-7 ) + Z( NN-9 ) - SQRT( Z( NN-11 ) )*SQRT( Z( NN-9 ) ) - A2
          IF ( GAP2 > ZERO .AND. GAP2 > B2*A2 ) THEN
             S = MAX( S, A2*( ONE-CNST2*A2*( B2 / GAP2 )*B2 ) )
          ELSE
             S = MAX( S, A2*( ONE-CNST2*B2 ) )
          END IF
       ELSE
          S = QURTR*DMIN2
          TTYPE = -11
       END IF

    ELSE IF ( N0IN > ( N0+2 ) ) THEN

       ! Case 12, more than two eigenvalues deflated. No information.

       S = ZERO
       TTYPE = -12

    END IF

    TAU = S

  END SUBROUTINE DLASQ4

!==========================================================================
!
!  DLASQ5 computes one dqds transform in ping-pong form. Used by sbdsqr and sstegr.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASQ5 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq5.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN,
!                          DNM1, DNM2, IEEE, EPS )
!
!       .. Scalar Arguments ..
!       LOGICAL            IEEE
!       INTEGER            I0, N0, PP
!       REAL(DP)           DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS
!       ..
!       .. Array Arguments ..
!       REAL(DP)           Z( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASQ5 computes one dqds transform in ping-pong form, one
! version for IEEE machines another for non IEEE machines.
!
!  Arguments:
!  ==========
!
! param[in] I0
!
!          I0 is INTEGER
!        First index.
!
! param[in] N0
!
!          N0 is INTEGER
!        Last index.
!
! param[in] Z
!
!          Z is REAL(DP)         array, dimension ( 4*N )
!        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!        an extra argument.
!
! param[in] PP
!
!          PP is INTEGER
!        PP=0 for ping, PP=1 for pong.
!
! param[in] TAU
!
!          TAU is REAL(DP)
!        This is the shift.
!
! param[in] SIGMA
!
!          SIGMA is REAL(DP)
!        This is the accumulated shift up to this step.
!
! param[out] DMIN
!
!          DMIN is REAL(DP)
!        Minimum value of d.
!
! param[out] DMIN1
!
!          DMIN1 is REAL(DP)
!        Minimum value of d, excluding D( N0 ).
!
! param[out] DMIN2
!
!          DMIN2 is REAL(DP)
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
! param[out] DN
!
!          DN is REAL(DP)
!        d(N0), the last value of d.
!
! param[out] DNM1
!
!          DNM1 is REAL(DP)
!        d(N0-1).
!
! param[out] DNM2
!
!          DNM2 is REAL(DP)
!        d(N0-2).
!
! param[in] IEEE
!
!          IEEE is LOGICAL
!        Flag for IEEE or non IEEE arithmetic.
!
! param[in] EPS
!
!          EPS is REAL(DP)
!        This is the value of epsilon used.
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

  SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, &
                     DN, DNM1, DNM2, IEEE, SMALL )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    LOGICAL,  INTENT(IN)    :: IEEE
    INTEGER , INTENT(INOUT) :: I0, N0, PP
    REAL(DP), INTENT(INOUT) :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, &
                               SIGMA, SMALL

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: Z( * )

    ! .. Parameter ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: HALF = 0.5_DP
    REAL(DP), PARAMETER :: EPS = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: J4, J4P2
    REAL(DP) :: D, EMIN, TEMP, DTHRESH

    ! .. Executable Statements ..

    IF ( ( N0-I0-1 ) <= 0 ) RETURN

    DTHRESH = SMALL*(SIGMA+TAU)
    IF ( TAU < DTHRESH*HALF ) TAU = ZERO

    IF ( ABS(TAU) > EPS ) THEN

       J4 = 4*I0 + PP - 3
       EMIN = Z( J4+4 )
       D = Z( J4 ) - TAU
       DMIN = D
       DMIN1 = -Z( J4 )

       IF ( IEEE ) THEN

          ! Code for IEEE arithmetic.

          IF ( PP == 0 ) THEN
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-2 ) = D + Z( J4-1 )
                TEMP = Z( J4+1 ) / Z( J4-2 )
                D = D*TEMP - TAU
                DMIN = MIN( DMIN, D )
                Z( J4 ) = Z( J4-1 )*TEMP
                EMIN = MIN( Z( J4 ), EMIN )
             END DO
          ELSE
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-3 ) = D + Z( J4 )
                TEMP = Z( J4+2 ) / Z( J4-3 )
                D = D*TEMP - TAU
                DMIN = MIN( DMIN, D )
                Z( J4-1 ) = Z( J4 )*TEMP
                EMIN = MIN( Z( J4-1 ), EMIN )
             END DO
          END IF

          ! Unroll last two steps.

          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*( N0-2 ) - PP
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM2 + Z( J4P2 )
          Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
          DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
          DMIN = MIN( DMIN, DNM1 )

          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM1 + Z( J4P2 )
          Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
          DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
          DMIN = MIN( DMIN, DN )

       ELSE

          ! Code for non IEEE arithmetic.

          IF ( PP == 0 ) THEN
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-2 ) = D + Z( J4-1 )
                IF ( D < ZERO ) THEN
                   RETURN
                ELSE
                   Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                   D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
                END IF
                DMIN = MIN( DMIN, D )
                EMIN = MIN( EMIN, Z( J4 ) )
             END DO
          ELSE
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-3 ) = D + Z( J4 )
                IF ( D < ZERO ) THEN
                   RETURN
                ELSE
                   Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                   D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
                END IF
                DMIN = MIN( DMIN, D )
                EMIN = MIN( EMIN, Z( J4-1 ) )
             END DO
          END IF

          ! Unroll last two steps.

          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*( N0-2 ) - PP
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM2 + Z( J4P2 )
          IF ( DNM2 < ZERO ) THEN
             RETURN
          ELSE
             Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
             DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
          END IF
          DMIN = MIN( DMIN, DNM1 )

          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM1 + Z( J4P2 )
          IF ( DNM1 < ZERO ) THEN
             RETURN
          ELSE
             Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
             DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
          END IF
          DMIN = MIN( DMIN, DN )

       END IF

    ELSE

       ! This is the version that sets d's to zero if they are small enough
       J4 = 4*I0 + PP - 3
       EMIN = Z( J4+4 )
       D = Z( J4 ) - TAU
       DMIN = D
       DMIN1 = -Z( J4 )

       IF ( IEEE ) THEN

          ! Code for IEEE arithmetic.

          IF ( PP == 0 ) THEN
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-2 ) = D + Z( J4-1 )
                TEMP = Z( J4+1 ) / Z( J4-2 )
                D = D*TEMP - TAU
                IF ( D < DTHRESH ) D = ZERO
                DMIN = MIN( DMIN, D )
                Z( J4 ) = Z( J4-1 )*TEMP
                EMIN = MIN( Z( J4 ), EMIN )
             END DO
          ELSE
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-3 ) = D + Z( J4 )
                TEMP = Z( J4+2 ) / Z( J4-3 )
                D = D*TEMP - TAU
                IF ( D < DTHRESH ) D = ZERO
                DMIN = MIN( DMIN, D )
                Z( J4-1 ) = Z( J4 )*TEMP
                EMIN = MIN( Z( J4-1 ), EMIN )
             END DO
          END IF

          ! Unroll last two steps.

          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*( N0-2 ) - PP
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM2 + Z( J4P2 )
          Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
          DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
          DMIN = MIN( DMIN, DNM1 )

          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM1 + Z( J4P2 )
          Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
          DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
          DMIN = MIN( DMIN, DN )

       ELSE

          ! Code for non IEEE arithmetic.

          IF ( PP == 0 ) THEN
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-2 ) = D + Z( J4-1 )
                IF ( D < ZERO ) THEN
                   RETURN
                ELSE
                   Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                   D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
                END IF
                IF ( D < DTHRESH) D = ZERO
                DMIN = MIN( DMIN, D )
                EMIN = MIN( EMIN, Z( J4 ) )
             END DO
          ELSE
             DO J4 = 4*I0, 4*( N0-3 ), 4
                Z( J4-3 ) = D + Z( J4 )
                IF ( D < ZERO ) THEN
                   RETURN
                ELSE
                   Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                   D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
                END IF
                IF ( D < DTHRESH) D = ZERO
                DMIN = MIN( DMIN, D )
                EMIN = MIN( EMIN, Z( J4-1 ) )
             END DO
          END IF

          ! Unroll last two steps.

          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*( N0-2 ) - PP
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM2 + Z( J4P2 )
          IF ( DNM2 < ZERO ) THEN
             RETURN
          ELSE
             Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
             DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
          END IF
          DMIN = MIN( DMIN, DNM1 )

          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z( J4-2 ) = DNM1 + Z( J4P2 )
          IF ( DNM1 < ZERO ) THEN
             RETURN
          ELSE
             Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
             DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
          END IF
          DMIN = MIN( DMIN, DN )

       END IF

    END IF

    Z( J4+2 ) = DN
    Z( 4*N0-PP ) = EMIN

  END SUBROUTINE DLASQ5

!==========================================================================
!
!  DLASQ6 computes one dqd transform in ping-pong form. Used by sbdsqr and sstegr.
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download DLASQ6 + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasq6.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,
!                          DNM1, DNM2 )
!
!       .. Scalar Arguments ..
!       INTEGER            I0, N0, PP
!       REAL(DP)           DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!       ..
!       .. Array Arguments ..
!       REAL(DP)           Z( * )
!       ..
!
!  Purpose:
!  =============
!
! DLASQ6 computes one dqd (shift equal to zero) transform in
! ping-pong form, with protection against underflow and overflow.
!
!  Arguments:
!  ==========
!
! param[in] I0
!
!          I0 is INTEGER
!        First index.
!
! param[in] N0
!
!          N0 is INTEGER
!        Last index.
!
! param[in] Z
!
!          Z is REAL(DP)         array, dimension ( 4*N )
!        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
!        an extra argument.
!
! param[in] PP
!
!          PP is INTEGER
!        PP=0 for ping, PP=1 for pong.
!
! param[out] DMIN
!
!          DMIN is REAL(DP)
!        Minimum value of d.
!
! param[out] DMIN1
!
!          DMIN1 is REAL(DP)
!        Minimum value of d, excluding D( N0 ).
!
! param[out] DMIN2
!
!          DMIN2 is REAL(DP)
!        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
!
! param[out] DN
!
!          DN is REAL(DP)
!        d(N0), the last value of d.
!
! param[out] DNM1
!
!          DNM1 is REAL(DP)
!        d(N0-1).
!
! param[out] DNM2
!
!          DNM2 is REAL(DP)
!        d(N0-2).
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

  SUBROUTINE DLASQ6( I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, &
                     DNM1, DNM2 )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: DLAMCH

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK computational routine (version 3.4.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! September 2012

    ! .. Scalar Arguments ..
    INTEGER,  INTENT(INOUT) :: I0, N0, PP
    REAL(DP), INTENT(INOUT) :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: Z( * )

    ! .. Parameter ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: EPS = TINY( ZERO )

    ! .. Local Scalars ..
    INTEGER  :: J4, J4P2
    REAL(DP) :: D, EMIN, SAFMIN, TEMP

    ! .. Executable Statements ..

    IF ( ( N0-I0-1 ) <= 0 ) RETURN

    SAFMIN = DLAMCH( 'Safe minimum' )
    J4 = 4*I0 + PP - 3
    EMIN = Z( J4+4 )
    D = Z( J4 )
    DMIN = D

    IF ( PP == 0 ) THEN
       DO J4 = 4*I0, 4*( N0-3 ), 4
          Z( J4-2 ) = D + Z( J4-1 )
          IF ( ABS( Z(J4-2) ) < EPS ) THEN
             Z( J4 ) = ZERO
             D = Z( J4+1 )
             DMIN = D
             EMIN = ZERO
          ELSE IF ( SAFMIN*Z( J4+1 ) < Z( J4-2 ) .AND. &
                    SAFMIN*Z( J4-2 ) < Z( J4+1 ) ) THEN
             TEMP = Z( J4+1 ) / Z( J4-2 )
             Z( J4 ) = Z( J4-1 )*TEMP
             D = D*TEMP
          ELSE
             Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
             D = Z( J4+1 )*( D / Z( J4-2 ) )
          END IF
          DMIN = MIN( DMIN, D )
          EMIN = MIN( EMIN, Z( J4 ) )
       END DO
    ELSE
       DO J4 = 4*I0, 4*( N0-3 ), 4
          Z( J4-3 ) = D + Z( J4 )
          IF ( ABS( Z(J4-3) ) < EPS ) THEN
             Z( J4-1 ) = ZERO
             D = Z( J4+2 )
             DMIN = D
             EMIN = ZERO
          ELSE IF ( SAFMIN*Z( J4+2 ) < Z( J4-3 ) .AND. &
                    SAFMIN*Z( J4-3 ) < Z( J4+2 ) ) THEN
             TEMP = Z( J4+2 ) / Z( J4-3 )
             Z( J4-1 ) = Z( J4 )*TEMP
             D = D*TEMP
          ELSE
             Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
             D = Z( J4+2 )*( D / Z( J4-3 ) )
          END IF
          DMIN = MIN( DMIN, D )
          EMIN = MIN( EMIN, Z( J4-1 ) )
       END DO
    END IF

    ! Unroll last two steps.

    DNM2 = D
    DMIN2 = DMIN
    J4 = 4*( N0-2 ) - PP
    J4P2 = J4 + 2*PP - 1
    Z( J4-2 ) = DNM2 + Z( J4P2 )
    IF ( ABS( Z(J4-2) ) < EPS ) THEN
       Z( J4 ) = ZERO
       DNM1 = Z( J4P2+2 )
       DMIN = DNM1
       EMIN = ZERO
    ELSE IF ( SAFMIN*Z( J4P2+2 ) < Z( J4-2 ) .AND. &
              SAFMIN*Z( J4-2 ) < Z( J4P2+2 ) ) THEN
       TEMP = Z( J4P2+2 ) / Z( J4-2 )
       Z( J4 ) = Z( J4P2 )*TEMP
       DNM1 = DNM2*TEMP
    ELSE
       Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
       DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) )
    END IF
    DMIN = MIN( DMIN, DNM1 )

    DMIN1 = DMIN
    J4 = J4 + 4
    J4P2 = J4 + 2*PP - 1
    Z( J4-2 ) = DNM1 + Z( J4P2 )
    IF ( ABS( Z(J4-2) ) < EPS ) THEN
       Z( J4 ) = ZERO
       DN = Z( J4P2+2 )
       DMIN = DN
       EMIN = ZERO
    ELSE IF ( SAFMIN*Z( J4P2+2 ) < Z( J4-2 ) .AND. &
              SAFMIN*Z( J4-2 ) < Z( J4P2+2 ) ) THEN
       TEMP = Z( J4P2+2 ) / Z( J4-2 )
       Z( J4 ) = Z( J4P2 )*TEMP
       DN = DNM1*TEMP
    ELSE
       Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
       DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) )
    END IF
    DMIN = MIN( DMIN, DN )

    Z( J4+2 ) = DN
    Z( 4*N0-PP ) = EMIN

  END SUBROUTINE DLASQ6

!=====================================================================
!
!  DGELS
!
!  Purpose
!  =======
!
!  DGELS solves overdetermined or underdetermined real linear systems
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
!  A       (input/output) REAL(DP)         array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!            if M >= N, A is overwritten by details of its QR
!                       factorization as returned by DGEQRF;
!            if M <  N, A is overwritten by details of its LQ
!                       factorization as returned by DGELQF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) REAL(DP)         array, dimension (LDB,NRHS)
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
!  WORK    (workspace/output) REAL(DP)         array, dimension (MAX(1,LWORK))
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

  SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: ILAENV, DLABAD, DLAMCH, DLANGE, DGELQF, DGEQRF, &
                               DLASCL, DLASET, DORMLQ, DORMQR, DTRTRS, XERBLA, LSAME

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK driver routine (version 3.3.1)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! April 2011

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: TRANS
    INTEGER,   INTENT(IN)  :: LDA, LDB, LWORK, M, N, NRHS
    INTEGER,   INTENT(OUT) :: INFO

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), B( LDB,* ), WORK( * )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

    ! .. Local Scalars ..
    LOGICAL  :: LQUERY, TPSD
    INTEGER  :: BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
    REAL(DP) :: ANRM, BIGNUM, BNRM, SMLNUM

    ! .. Local Arrays ..
    REAL(DP) :: RWORK( 1 )

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
    ELSE IF ( LWORK < MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY ) THEN
       INFO = -10
    END IF

    ! Figure out optimal block size

    IF ( INFO == 0 .OR. INFO == -10 ) THEN

       TPSD = .TRUE.
       IF ( LSAME( TRANS, 'N' ) ) TPSD = .FALSE.

       IF ( M >= N ) THEN
          NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
          IF ( TPSD ) THEN
             NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LN', M, NRHS, N, -1 ) )
          ELSE
             NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N, -1 ) )
          END IF
       ELSE
          NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
          IF ( TPSD ) THEN
             NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M, -1 ) )
          ELSE
             NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LN', N, NRHS, M, -1 ) )
          END IF
       END IF

       WSIZE = MAX( 1, MN+MAX( MN, NRHS )*NB )
       WORK( 1 ) = REAL( WSIZE, DP )

    END IF

    IF ( INFO /= 0 ) THEN
       CALL XERBLA( 'DGELS ', -INFO )
       RETURN
    ELSE IF ( LQUERY ) THEN
       RETURN
    END IF

    ! Quick return if possible

    IF ( MIN( M, N, NRHS ) == 0 ) THEN
       CALL DLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
       RETURN
    END IF

    ! Get machine parameters

    SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
    BIGNUM = ONE / SMLNUM
    CALL DLABAD( SMLNUM, BIGNUM )

    ! Scale A, B if max element outside range [SMLNUM,BIGNUM]

    ANRM = DLANGE( 'M', M, N, A, LDA, RWORK )
    IASCL = 0

    IF ( ANRM > ZERO .AND. ANRM < SMLNUM ) THEN

       ! Scale matrix norm up to SMLNUM
       CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
       IASCL = 1

    ELSE IF ( ANRM > BIGNUM ) THEN

       ! Scale matrix norm down to BIGNUM
       CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
       IASCL = 2

    ELSE IF ( ANRM < EPS ) THEN

       ! Matrix all zero. Return zero solution.
       CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
       WORK( 1 ) = REAL( WSIZE, DP )
       RETURN

    END IF

    BROW = M
    IF ( TPSD ) BROW = N
    BNRM = DLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
    IBSCL = 0

    IF ( BNRM > ZERO .AND. BNRM < SMLNUM ) THEN

       ! Scale matrix norm up to SMLNUM
       CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO )
       IBSCL = 1

    ELSE IF ( BNRM > BIGNUM ) THEN

       ! Scale matrix norm down to BIGNUM
       CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO )
       IBSCL = 2

    END IF

    IF ( M >= N ) THEN

       ! compute QR factorization of A

       CALL DGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

       ! workspace at least N, optimally N*NB

       IF ( .NOT.TPSD ) THEN

          ! Least-Squares Problem min || A * X - B ||

          ! B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

          CALL DORMQR( 'Left', 'Transpose', M, NRHS, N, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          ! B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

          CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS, &
                       A, LDA, B, LDB, INFO )

          IF ( INFO > 0 ) THEN
             RETURN
          END IF

          SCLLEN = N

       ELSE

          ! Overdetermined system of equations A**T * X = B

          ! B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

          CALL DTRTRS( 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
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

          CALL DORMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          SCLLEN = M

       END IF

    ELSE

       ! Compute LQ factorization of A

       CALL DGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

       ! workspace at least M, optimally M*NB.

       IF ( .NOT.TPSD ) THEN

          ! underdetermined system of equations A * X = B

          ! B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

          CALL DTRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS, &
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

          CALL DORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          SCLLEN = N

       ELSE

          ! overdetermined system min || A**T * X - B ||

          ! B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

          CALL DORMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA, &
                       WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

          ! workspace at least NRHS, optimally NRHS*NB

          ! B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

          CALL DTRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS, &
                       A, LDA, B, LDB, INFO )

          IF ( INFO > 0 ) THEN
             RETURN
          END IF

          SCLLEN = M

       END IF

    END IF

    ! Undo scaling

    IF ( IASCL == 1 ) THEN
       CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO )
    ELSE IF ( IASCL == 2 ) THEN
       CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO )
    END IF

    IF ( IBSCL == 1 ) THEN
       CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
    ELSE IF ( IBSCL == 2 ) THEN
       CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
    END IF

    WORK( 1 ) = REAL( WSIZE, DP )

  END SUBROUTINE DGELS

!=====================================================================
!
!  DLABAD
!
!  Purpose
!  =======
!
!  DLABAD takes as input the values computed by DLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by DLAMCH.  This subroutine is needed because
!  DLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) REAL(DP)
!          On entry, the underflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) REAL(DP)
!          On entry, the overflow threshold as computed by DLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!=====================================================================

  SUBROUTINE DLABAD( SMALL, LARGE )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! LAPACK auxiliary routine (version 3.2)
    ! LAPACK is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2006

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(INOUT) :: LARGE, SMALL

    ! .. Executable Statements ..

    ! If it looks like we're on a Cray, take the square root of
    ! SMALL and LARGE to avoid overflow and underflow problems.

    IF ( LOG10( LARGE ) > 2000.D0 ) THEN
       SMALL = SQRT( SMALL )
       LARGE = SQRT( LARGE )
    END IF

  END SUBROUTINE DLABAD

!=====================================================================
!
!  DTRTRS
!
!  Purpose
!  =======
!
!  DTRTRS solves a triangular system of the form
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
!  A       (input) REAL(DP)         array, dimension (LDA,N)
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
!  B       (input/output) REAL(DP)         array, dimension (LDB,NRHS)
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

  SUBROUTINE DTRTRS( UPLOC, TRANSC, DIAGC, N, NRHS, A, LDA, B, LDB, INFO )

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA, DTRSM

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
    REAL(DP), INTENT(INOUT) :: A( LDA,* ), B( LDB,* )

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

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
       CALL XERBLA( 'DTRTRS', -INFO )
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

    CALL DTRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB )

  END SUBROUTINE DTRTRS

!=====================================================================
!
!  DTRSM
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A**T.
!
!  The matrix X is overwritten on B.
!
!  Arguments
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A**T.
!
!              TRANSA = 'C' or 'c'   op( A ) = A**T.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(DP)        .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL(DP)         array of DIMENSION ( LDA, k ), where k is m
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
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL(DP)         array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
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

  SUBROUTINE DTRSM(SIDEC,UPLOC,TRANSAC,DIAGC,M,N,ALPHA,A,LDA,B,LDB)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! .. Scalar Arguments ..
    REAL(DP),     INTENT(IN) :: ALPHA
    INTEGER,      INTENT(IN) :: LDA,LDB,M,N
    CHARACTER(*), INTENT(IN) :: DIAGC,SIDEC,TRANSAC,UPLOC

    ! .. Array Arguments ..
    REAL(DP),  INTENT(INOUT) :: A(LDA,*),B(LDB,*)

    ! .. Local Scalars ..
    REAL(DP)  :: TEMP
    INTEGER   :: I,INFO,J,K,NROWA
    LOGICAL   :: LSIDE,NOUNIT,UPPER
    CHARACTER :: DIAG,SIDE,TRANSA,UPLO

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

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
       CALL XERBLA('DTRSM ',INFO)
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

  END SUBROUTINE DTRSM

!=====================================================================
!
!  DSCAL
!
!  Purpose
!  =======
!
!     DSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!
!  Further Details
!  ===============
!
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx <= 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
!=====================================================================

  SUBROUTINE DSCAL(N,DA,DX,INCX)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level1 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL(DP), INTENT(IN) :: DA
    INTEGER,  INTENT(IN) :: INCX, N

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: DX(*)

    ! .. Local Scalars ..
    INTEGER :: I,M,MP1,NINCX

    IF (N <= 0 .OR. INCX <= 0) RETURN

    IF (INCX == 1) THEN

       ! code for increment equal to 1

       ! clean-up loop
       M = MOD(N,5)
       IF (M /= 0) THEN
          DO I = 1,M
             DX(I) = DA*DX(I)
          END DO
          IF (N < 5) RETURN
       END IF
       MP1 = M + 1
       DO I = MP1,N,5
          DX(I)   = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
       END DO

    ELSE

       ! code for increment not equal to 1

       NINCX = N*INCX
       DO I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
       END DO

    END IF

  END SUBROUTINE DSCAL

!=====================================================================
!
!  DGEMM
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X**T,
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Arguments
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A**T.
!
!              TRANSA = 'C' or 'c',  op( A ) = A**T.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B**T.
!
!              TRANSB = 'C' or 'c',  op( B ) = B**T.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(DP)        .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(DP)         array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL(DP)         array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - REAL(DP)        .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - REAL(DP)         array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!!*           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!  Further Details
!  ===============
!
!  Level 3 Blas routine.
!
!     Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!=====================================================================

  SUBROUTINE DGEMM(TRANSAC,TRANSBC,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    USE ISO_FORTRAN_ENV, ONLY: DP => REAL64
    IMPORT,              ONLY: LSAME, XERBLA

    IMPLICIT NONE (EXTERNAL, TYPE)

    ! Reference BLAS level3 routine (version 3.4.0)
    ! Reference BLAS is a software package provided by Univ. of Tennessee,
    ! Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd
    ! November 2011

    ! .. Scalar Arguments ..
    REAL(DP),     INTENT(IN) :: ALPHA, BETA
    INTEGER,      INTENT(IN) :: K, LDA, LDB, LDC, M, N
    CHARACTER(*), INTENT(IN) :: TRANSAC, TRANSBC

    ! .. Array Arguments ..
    REAL(DP), INTENT(INOUT) :: A(LDA,*),B(LDB,*),C(LDC,*)

    ! .. Local Scalars ..
    REAL(DP)  :: TEMP
    INTEGER   :: I,INFO,J,L,NCOLA,NROWA,NROWB
    LOGICAL   :: NOTA,NOTB
    CHARACTER :: TRANSA,TRANSB

    ! .. Parameters ..
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
    REAL(DP), PARAMETER :: ONE  = 1.0_DP
    REAL(DP), PARAMETER :: EPS  = TINY( ZERO )

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
       CALL XERBLA('DGEMM ',INFO)
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

  END SUBROUTINE DGEMM

End Module Lapack_Double
