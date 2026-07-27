Module Lapack_Utils

  ! This module contains Lapack utililty routines used by both the single and
  ! double precision libraries

Contains

!=====================================================================
!
!  LSAME
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      LOGICAL FUNCTION LSAME( CA, CB )
!
!      IMPLICIT NONE (EXTERNAL, TYPE)
!
!     .. Scalar Arguments ..
!      CHARACTER          CA, CB
!     ..
!
!  Purpose:
!  =============
!
! LSAME returns .TRUE. if CA is the same letter as CB regardless of
! case.
!
!  Arguments:
!  ==========
!
! param[in] CA
!
! param[in] CB
!
!          CA and CB specify the single characters to be compared.
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

  LOGICAL FUNCTION LSAME( CA, CB )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    !  -- LAPACK auxiliary routine (version 3.4.0) --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     November 2011

    ! .. Scalar Arguments ..
    CHARACTER, INTENT(IN) :: CA, CB

    ! .. Local Scalars ..
    INTEGER :: INTA, INTB, ZCODE

    ! .. Executable Statements ..

    ! Test if the characters are equal

    LSAME = CA == CB
    IF ( LSAME ) RETURN

    ! Now test for equivalence if both characters are alphabetic.

    ZCODE = ICHAR( 'Z' )

    ! Use 'Z' rather than 'A' so that ASCII can be detected on Prime
    ! machines, on which ICHAR returns a value with bit 8 set.
    ! ICHAR('A') on Prime machines returns 193 which is the same as
    ! ICHAR('A') on an EBCDIC machine.

    INTA = ICHAR( CA )
    INTB = ICHAR( CB )

    IF ( ZCODE == 90 .OR. ZCODE == 122 ) THEN

       ! ASCII is assumed - ZCODE is the ASCII code of either lower or
       ! upper case 'Z'.

       IF ( INTA >= 97 .AND. INTA <= 122 ) INTA = INTA - 32
       IF ( INTB >= 97 .AND. INTB <= 122 ) INTB = INTB - 32

    ELSE IF ( ZCODE == 233 .OR. ZCODE == 169 ) THEN

       ! EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
       ! upper case 'Z'.

       IF ( INTA >= 129 .AND. INTA <= 137 .OR. &
            INTA >= 145 .AND. INTA <= 153 .OR. &
            INTA >= 162 .AND. INTA <= 169 ) INTA = INTA + 64
       IF ( INTB >= 129 .AND. INTB <= 137 .OR. &
            INTB >= 145 .AND. INTB <= 153 .OR. &
            INTB >= 162 .AND. INTB <= 169 ) INTB = INTB + 64

    ELSE IF ( ZCODE == 218 .OR. ZCODE == 250 ) THEN

       ! ASCII is assumed, on Prime machines - ZCODE is the ASCII code
       ! plus 128 of either lower or upper case 'Z'.

       IF ( INTA >= 225 .AND. INTA <= 250 ) INTA = INTA - 32
       IF ( INTB >= 225 .AND. INTB <= 250 ) INTB = INTB - 32
    END IF
    LSAME = INTA == INTB

  END FUNCTION LSAME

!=====================================================================
!
!  XERBLA
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download XERBLA + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f
!
!  Definition:
!  ===========
!
!       SUBROUTINE XERBLA( SRNAME, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!
! Purpose:
!  =============
!
! XERBLA  is an error handler for the LAPACK routines.
! It is called by an LAPACK routine if an input parameter has an
! invalid value.  A message is printed and execution stops.
!
! Installers may consider modifying the STOP statement in order to
! call system-specific exception-handling facilities.
!
!  Arguments:
!  ==========
!
! param[in] SRNAME
!
!          SRNAME is CHARACTER*(*)
!          The name of the routine which called XERBLA.
!
! param[in] INFO
!
!          INFO is INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
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
! \ingroup auxOTHERauxiliary

!=====================================================================

  SUBROUTINE XERBLA( SRNAME, INFO )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    !  -- LAPACK auxiliary routine (version 3.4.0) --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     November 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: SRNAME
    INTEGER,      INTENT(IN) :: INFO

    CHARACTER(80) :: FORMAT = &
         "(' ** On entry to ', A, ' parameter number ', I2, ' had an illegal value')"

    ! .. Executable Statements ..

    WRITE(*,FORMAT) TRIM(SRNAME), INFO

    STOP

  END SUBROUTINE XERBLA

!=====================================================================
!
!  IEEECK
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download IEEECK + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f
!
!  Definition:
!===========
!
!       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!       .. Scalar Arguments ..
!       INTEGER            ISPEC
!       REAL               ONE, ZERO
!       ..
!
! Purpose:
!=============
!
! IEEECK is called from the ILAENV to verify that Infinity and
! possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments:
!==========
!
! param[in] ISPEC
!
!          ISPEC is INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
! param[in] ZERO
!
!          ZERO is REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
! param[in] ONE
!
!          ONE is REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!  Authors:
!========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
!=====================================================================

  INTEGER FUNCTION IEEECK( ISPEC, ZERO, ONE )

    USE IEEE_ARITHMETIC, ONLY: IEEE_QUIET_NE, IEEE_IS_NAN
    IMPORT,              NONE

    IMPLICIT NONE (EXTERNAL, TYPE)

    !  -- LAPACK auxiliary routine (version 3.4.0) --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     November 2011

    ! .. Scalar Arguments ..
    INTEGER, INTENT(IN) :: ISPEC
    REAL,    INTENT(IN) :: ONE, ZERO

    ! .. Local Scalars ..
    REAL :: NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, NEGZRO, NEWZRO, POSINF

    ! .. Executable Statements ..
    IEEECK = 1

    POSINF = ONE / ZERO
    IF ( POSINF <= ONE ) THEN
       IEEECK = 0
       RETURN
    END IF

    NEGINF = -ONE / ZERO
    IF ( NEGINF >= ZERO ) THEN
       IEEECK = 0
       RETURN
    END IF

    NEGZRO = ONE / ( NEGINF+ONE )
    IF ( IEEE_QUIET_NE( NEGZRO,ZERO ) ) THEN
       IEEECK = 0
       RETURN
    END IF

    NEGINF = ONE / NEGZRO
    IF ( NEGINF >= ZERO ) THEN
       IEEECK = 0
       RETURN
    END IF

    NEWZRO = NEGZRO + ZERO
    IF ( IEEE_QUIET_NE( NEWZRO,ZERO ) ) THEN
       IEEECK = 0
       RETURN
    END IF

    POSINF = ONE / NEWZRO
    IF ( POSINF <= ONE ) THEN
       IEEECK = 0
       RETURN
    END IF

    NEGINF = NEGINF*POSINF
    IF ( NEGINF >= ZERO ) THEN
       IEEECK = 0
       RETURN
    END IF

    POSINF = POSINF*POSINF
    IF ( POSINF <= ONE ) THEN
       IEEECK = 0
       RETURN
    END IF

    ! Return if we were only asked to check infinity arithmetic

    IF ( ISPEC == 0 ) RETURN

    NAN1 = POSINF + NEGINF

    NAN2 = POSINF / NEGINF

    NAN3 = POSINF / POSINF

    NAN4 = POSINF*ZERO

    NAN5 = NEGINF*NEGZRO

    NAN6 = NAN5*ZERO

    IF ( .NOT. IEEE_IS_NAN(NAN1)) THEN
       IEEECK = 0
       RETURN
    END IF

    IF ( .NOT. IEEE_IS_NAN(NAN2) ) THEN
       IEEECK = 0
       RETURN
    END IF

    IF ( .NOT. IEEE_IS_NAN(NAN3) ) THEN
       IEEECK = 0
       RETURN
    END IF

    IF ( .NOT. IEEE_IS_NAN(NAN4) ) THEN
       IEEECK = 0
       RETURN
    END IF

    IF ( .NOT. IEEE_IS_NAN(NAN5) ) THEN
       IEEECK = 0
       RETURN
    END IF

    IF ( .NOT. IEEE_IS_NAN(NAN6) ) THEN
       IEEECK = 0
       RETURN
    END IF

  END FUNCTION IEEECK

!=====================================================================
!
!  ILAENV
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download ILAENV + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f
!
!  Definition:
!===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!
! Purpose:
!=============
!
! ILAENV is called from the LAPACK routines to choose problem-dependent
! parameters for the local environment.  See ISPEC for a description of
! the parameters.
!
! ILAENV returns an INTEGER
! if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
! if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
! This version provides a set of parameters which should give good,
! but not optimal, performance on many of the currently available
! computers.  Users are encouraged to modify this subroutine to set
! the tuning parameters for their particular machine using the option
! and problem size information in the arguments.
!
! This routine will not function correctly if it is converted to all
! lower case.  Converting it to all upper case is allowed.
!
!  Arguments:
!==========
!
! param[in] ISPEC
!
!          ISPEC is INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines (DEPRECATED)
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR method
!               for nonsymmetric eigenvalue problems (DEPRECATED)
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!          12 <= ISPEC <= 16:
!               xHSEQR or one of its subroutines,
!               see IPARMQ for detailed explanation
!
! param[in] NAME
!
!          NAME is CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
! param[in] OPTS
!
!          OPTS is CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
! param[in] N1
!
!          N1 is INTEGER
!
! param[in] N2
!
!          N2 is INTEGER
!
! param[in] N3
!
!          N3 is INTEGER
!
! param[in] N4
!
!          N4 is INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
!  Authors:
!========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
! Further Details:
!=====================
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF ( NB <= 1 ) NB = MAX( 1, N )
!
!=====================================================================

  INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )

    IMPORT,  ONLY: IEEECK, IPARMQ
    IMPLICIT NONE (EXTERNAL, TYPE)

    !  -- LAPACK auxiliary routine (version 3.4.0) --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     November 2011

    ! .. Scalar Arguments ..
    CHARACTER(*), INTENT(IN) :: NAME, OPTS
    INTEGER,      INTENT(IN) :: ISPEC, N1, N2, N3, N4

    ! .. Local Scalars ..
    INTEGER   :: I, IC, IZ, NB, NBMIN, NX
    LOGICAL   :: CNAME, SNAME
    CHARACTER :: C1*1, C2*2, C4*2, C3*3, SUBNAM*6

    ! .. Executable Statements ..

    SELECT CASE (ISPEC)

    CASE (1:3)

       ! Convert NAME to upper case if the first character is lower case.

       ILAENV = 1
       SUBNAM = NAME
       IC = ICHAR( SUBNAM( 1: 1 ) )
       IZ = ICHAR( 'Z' )

       IF ( IZ == 90 .OR. IZ == 122 ) THEN

          ! ASCII character set

          IF ( IC >= 97 .AND. IC <= 122 ) THEN
             SUBNAM( 1: 1 ) = CHAR( IC-32 )
             DO I = 2, 6
                IC = ICHAR( SUBNAM( I: I ) )
                IF ( IC >= 97 .AND. IC <= 122 ) SUBNAM( I:I ) = CHAR( IC-32 )
             END DO
          END IF

       ELSE IF ( IZ == 233 .OR. IZ == 169 ) THEN

          ! EBCDIC character set

          IF ( ( IC >= 129 .AND. IC <= 137 ) .OR. &
               ( IC >= 145 .AND. IC <= 153 ) .OR. &
               ( IC >= 162 .AND. IC <= 169 ) ) THEN
             SUBNAM( 1: 1 ) = CHAR( IC+64 )
             DO I = 2, 6
                IC = ICHAR( SUBNAM( I: I ) )
                IF ( ( IC >= 129 .AND. IC <= 137 ) .OR. &
                     ( IC >= 145 .AND. IC <= 153 ) .OR. &
                     ( IC >= 162 .AND. IC <= 169 ) ) SUBNAM( I:I ) = CHAR( IC+64 )
             END DO
          END IF

       ELSE IF ( IZ == 218 .OR. IZ == 250 ) THEN

          ! Prime machines:  ASCII+128

          IF ( IC >= 225 .AND. IC <= 250 ) THEN
             SUBNAM( 1: 1 ) = CHAR( IC-32 )
             DO I = 2, 6
                IC = ICHAR( SUBNAM( I: I ) )
                IF ( IC >= 225 .AND. IC <= 250 ) &
                     SUBNAM( I: I ) = CHAR( IC-32 )
             END DO
          END IF
       END IF

       C1 = SUBNAM( 1: 1 )
       SNAME = C1 == 'S' .OR. C1 == 'D'
       CNAME = C1 == 'C' .OR. C1 == 'Z'
       IF ( .NOT.( CNAME .OR. SNAME ) ) RETURN
       C2 = SUBNAM( 2: 3 )
       C3 = SUBNAM( 4: 6 )
       C4 = C3( 2: 3 )

       IF (ISPEC == 1) THEN

          ! ISPEC = 1:  block size

          ! In these examples, separate code is provided for setting NB for
          ! real and complex.  We assume that NB will take the same value in
          ! single or double precision.

          NB = 1

          IF ( C2 == 'GE' ) THEN
             IF ( C3 == 'TRF' ) THEN
                IF ( SNAME ) THEN
                   NB = 64
                ELSE
                   NB = 64
                END IF
             ELSE IF ( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. &
                       C3 == 'QLF' ) THEN
                IF ( SNAME ) THEN
                   NB = 32
                ELSE
                   NB = 32
                END IF
             ELSE IF ( C3 == 'HRD' ) THEN
                IF ( SNAME ) THEN
                   NB = 32
                ELSE
                   NB = 32
                END IF
             ELSE IF ( C3 == 'BRD' ) THEN
                IF ( SNAME ) THEN
                   NB = 32
                ELSE
                   NB = 32
                END IF
             ELSE IF ( C3 == 'TRI' ) THEN
                IF ( SNAME ) THEN
                   NB = 64
                ELSE
                   NB = 64
                END IF
             END IF
          ELSE IF ( C2 == 'PO' ) THEN
             IF ( C3 == 'TRF' ) THEN
                IF ( SNAME ) THEN
                   NB = 64
                ELSE
                   NB = 64
                END IF
             END IF
          ELSE IF ( C2 == 'SY' ) THEN
             IF ( C3 == 'TRF' ) THEN
                IF ( SNAME ) THEN
                   NB = 64
                ELSE
                   NB = 64
                END IF
             ELSE IF ( SNAME .AND. C3 == 'TRD' ) THEN
                NB = 32
             ELSE IF ( SNAME .AND. C3 == 'GST' ) THEN
                NB = 64
             END IF
          ELSE IF ( CNAME .AND. C2 == 'HE' ) THEN
             IF ( C3 == 'TRF' ) THEN
                NB = 64
             ELSE IF ( C3 == 'TRD' ) THEN
                NB = 32
             ELSE IF ( C3 == 'GST' ) THEN
                NB = 64
             END IF
          ELSE IF ( SNAME .AND. C2 == 'OR' ) THEN
             IF ( C3( 1: 1 ) == 'G' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NB = 32
                END IF
             ELSE IF ( C3( 1: 1 ) == 'M' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NB = 32
                END IF
             END IF
          ELSE IF ( CNAME .AND. C2 == 'UN' ) THEN
             IF ( C3( 1: 1 ) == 'G' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NB = 32
                END IF
             ELSE IF ( C3( 1: 1 ) == 'M' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NB = 32
                END IF
             END IF
          ELSE IF ( C2 == 'GB' ) THEN
             IF ( C3 == 'TRF' ) THEN
                IF ( SNAME ) THEN
                   IF ( N4 <= 64 ) THEN
                      NB = 1
                   ELSE
                      NB = 32
                   END IF
                ELSE
                   IF ( N4 <= 64 ) THEN
                      NB = 1
                   ELSE
                      NB = 32
                   END IF
                END IF
             END IF
          ELSE IF ( C2 == 'PB' ) THEN
             IF ( C3 == 'TRF' ) THEN
                IF ( SNAME ) THEN
                   IF ( N2 <= 64 ) THEN
                      NB = 1
                   ELSE
                      NB = 32
                   END IF
                ELSE
                   IF ( N2 <= 64 ) THEN
                      NB = 1
                   ELSE
                      NB = 32
                   END IF
                END IF
             END IF
          ELSE IF ( C2 == 'TR' ) THEN
             IF ( C3 == 'TRI' ) THEN
                IF ( SNAME ) THEN
                   NB = 64
                ELSE
                   NB = 64
                END IF
             END IF
          ELSE IF ( C2 == 'LA' ) THEN
             IF ( C3 == 'UUM' ) THEN
                IF ( SNAME ) THEN
                   NB = 64
                ELSE
                   NB = 64
                END IF
             END IF
          ELSE IF ( SNAME .AND. C2 == 'ST' ) THEN
             IF ( C3 == 'EBZ' ) THEN
                NB = 1
             END IF
          END IF
          ILAENV = NB

       ELSE IF (ISPEC == 2) THEN

          ! ISPEC = 2:  minimum block size

          NBMIN = 2
          IF ( C2 == 'GE' ) THEN
             IF ( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. C3 == &
                  'QLF' ) THEN
                IF ( SNAME ) THEN
                   NBMIN = 2
                ELSE
                   NBMIN = 2
                END IF
             ELSE IF ( C3 == 'HRD' ) THEN
                IF ( SNAME ) THEN
                   NBMIN = 2
                ELSE
                   NBMIN = 2
                END IF
             ELSE IF ( C3 == 'BRD' ) THEN
                IF ( SNAME ) THEN
                   NBMIN = 2
                ELSE
                   NBMIN = 2
                END IF
             ELSE IF ( C3 == 'TRI' ) THEN
                IF ( SNAME ) THEN
                   NBMIN = 2
                ELSE
                   NBMIN = 2
                END IF
             END IF
          ELSE IF ( C2 == 'SY' ) THEN
             IF ( C3 == 'TRF' ) THEN
                IF ( SNAME ) THEN
                   NBMIN = 8
                ELSE
                   NBMIN = 8
                END IF
             ELSE IF ( SNAME .AND. C3 == 'TRD' ) THEN
                NBMIN = 2
             END IF
          ELSE IF ( CNAME .AND. C2 == 'HE' ) THEN
             IF ( C3 == 'TRD' ) THEN
                NBMIN = 2
             END IF
          ELSE IF ( SNAME .AND. C2 == 'OR' ) THEN
             IF ( C3( 1: 1 ) == 'G' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NBMIN = 2
                END IF
             ELSE IF ( C3( 1: 1 ) == 'M' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NBMIN = 2
                END IF
             END IF
          ELSE IF ( CNAME .AND. C2 == 'UN' ) THEN
             IF ( C3( 1: 1 ) == 'G' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NBMIN = 2
                END IF
             ELSE IF ( C3( 1: 1 ) == 'M' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NBMIN = 2
                END IF
             END IF
          END IF
          ILAENV = NBMIN

       ELSE

          ! ISPEC = 3:  crossover point

          NX = 0
          IF ( C2 == 'GE' ) THEN
             IF ( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. C3 == 'QLF' ) THEN
                IF ( SNAME ) THEN
                   NX = 128
                ELSE
                   NX = 128
                END IF
             ELSE IF ( C3 == 'HRD' ) THEN
                IF ( SNAME ) THEN
                   NX = 128
                ELSE
                   NX = 128
                END IF
             ELSE IF ( C3 == 'BRD' ) THEN
                IF ( SNAME ) THEN
                   NX = 128
                ELSE
                   NX = 128
                END IF
             END IF
          ELSE IF ( C2 == 'SY' ) THEN
             IF ( SNAME .AND. C3 == 'TRD' ) THEN
                NX = 32
             END IF
          ELSE IF ( CNAME .AND. C2 == 'HE' ) THEN
             IF ( C3 == 'TRD' ) THEN
                NX = 32
             END IF
          ELSE IF ( SNAME .AND. C2 == 'OR' ) THEN
             IF ( C3( 1: 1 ) == 'G' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NX = 128
                END IF
             END IF
          ELSE IF ( CNAME .AND. C2 == 'UN' ) THEN
             IF ( C3( 1: 1 ) == 'G' ) THEN
                IF ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. C4 == 'QL' .OR. &
                     C4 == 'HR' .OR. C4 == 'TR' .OR. C4 == 'BR' ) THEN
                   NX = 128
                END IF
             END IF
          END IF
          ILAENV = NX

       ENDIF

    CASE (4)

       ! ISPEC = 4:  number of shifts (used by xHSEQR)

       ILAENV = 6

    CASE (5)

       ! ISPEC = 5:  minimum column dimension (not used)

       ILAENV = 2

    CASE (6)

       ! ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)

       ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )

    CASE (7)

       ! ISPEC = 7:  number of processors (not used)

       ILAENV = 1

    CASE (8)

       ! ISPEC = 8:  crossover point for multishift (used by xHSEQR)

       ILAENV = 50

    CASE (9)

       ! ISPEC = 9:  maximum size of the subproblems at the bottom of the
       !             computation tree in the divide-and-conquer algorithm
       !             (used by xGELSD and xGESDD)

       ILAENV = 25

    CASE (10)

       ! ISPEC = 10: ieee NaN arithmetic can be trusted not to trap

       ILAENV = 1
       IF ( ILAENV == 1 ) THEN
          ILAENV = IEEECK( 1, 0.0, 1.0 )
       END IF

    CASE (11)

       ! ISPEC = 11: infinity arithmetic can be trusted not to trap

       ILAENV = 1
       IF ( ILAENV == 1 ) THEN
          ILAENV = IEEECK( 0, 0.0, 1.0 )
       END IF

    CASE (12:16)

       ! 12 <= ISPEC <= 16: xHSEQR or one of its subroutines.

       ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )

    CASE DEFAULT

       ! Invalid value for ISPEC

       ILAENV = -1

    END SELECT

  END FUNCTION ILAENV

!=====================================================================
!
!  IPARMQ
!
!=========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
! Download IPARMQ + dependencies
! http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f
!
!  Definition:
!===========
!
!       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!       IMPLICIT NONE (EXTERNAL, TYPE)
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, ISPEC, LWORK, N
!       CHARACTER          NAME*( * ), OPTS*( * )
!
! Purpose:
!=============
!
!      This program sets problem and machine dependent parameters
!      useful for xHSEQR and its subroutines. It is called whenever
!      ILAENV is called with 12 <= ISPEC <= 16
!
!  Arguments:
!==========
!
! param[in] ISPEC
!
!          ISPEC is integer scalar
!              ISPEC specifies which tunable parameter IPARMQ should
!              return.
!
!              ISPEC=12: (INMIN)  Matrices of order nmin or less
!                        are sent directly to xLAHQR, the implicit
!                        double shift QR algorithm.  NMIN must be
!                        at least 11.
!
!              ISPEC=13: (INWIN)  Size of the deflation window.
!                        This is best set greater than or equal to
!                        the number of simultaneous shifts NS.
!                        Larger matrices benefit from larger deflation
!                        windows.
!
!              ISPEC=14: (INIBL) Determines when to stop nibbling and
!                        invest in an (expensive) multi-shift QR sweep.
!                        If the aggressive early deflation subroutine
!                        finds LD converged eigenvalues from an order
!                        NW deflation window and LD > (NW*NIBBLE)/100,
!                        then the next QR sweep is skipped and early
!                        deflation is applied immediately to the
!                        remaining active diagonal block.  Setting
!                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!                        multi-shift QR sweep whenever early deflation
!                        finds a converged eigenvalue.  Setting
!                        IPARMQ(ISPEC=14) greater than or equal to 100
!                        prevents TTQRE from skipping a multi-shift
!                        QR sweep.
!
!              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!                        a multi-shift QR iteration.
!
!              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!                        following meanings.
!                        0:  During the multi-shift QR sweep,
!                            xLAQR5 does not accumulate reflections and
!                            does not use matrix-matrix multiply to
!                            update the far-from-diagonal matrix
!                            entries.
!                        1:  During the multi-shift QR sweep,
!                            xLAQR5 and/or xLAQRaccumulates reflections and uses
!                            matrix-matrix multiply to update the
!                            far-from-diagonal matrix entries.
!                        2:  During the multi-shift QR sweep.
!                            xLAQR5 accumulates reflections and takes
!                            advantage of 2-by-2 block structure during
!                            matrix-matrix multiplies.
!                        (If xTRMM is slower than xGEMM, then
!                        IPARMQ(ISPEC=16)=1 may be more efficient than
!                        IPARMQ(ISPEC=16)=2 despite the greater level of
!                        arithmetic work implied by the latter choice.)
!
! param[in] NAME
!
!          NAME is character string
!               Name of the calling subroutine
!
! param[in] OPTS
!
!          OPTS is character string
!               This is a concatenation of the string arguments to
!               TTQRE.
!
! param[in] N
!
!          N is integer scalar
!               N is the order of the Hessenberg matrix H.
!
! param[in] ILO
!
!          ILO is INTEGER
!
! param[in] IHI
!
!          IHI is INTEGER
!               It is assumed that H is already upper triangular
!               in rows and columns 1:ILO-1 and IHI+1:N.
!
! param[in] LWORK
!
!          LWORK is integer scalar
!               The amount of workspace available.
!
!  Authors:
!========
!
! Univ. of Tennessee
! Univ. of California Berkeley
! Univ. of Colorado Denver
! NAG Ltd.
!
! November 2011
!
! Further Details:
!=====================
!
!       Little is known about how best to choose these parameters.
!       It is possible to use different values of the parameters
!       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!
!       It is probably best to choose different parameters for
!       different matrices and different parameters at different
!       times during the iteration, but this has not been
!       implemented --- yet.
!
!       The best choices of most of the parameters depend
!       in an ill-understood way on the relative execution
!       rate of xLAQR3 and xLAQR5 and on the nature of each
!       particular eigenvalue problem.  Experiment may be the
!       only practical way to determine which choices are most
!       effective.
!
!       Following is a list of default values supplied by IPARMQ.
!       These defaults may be adjusted in order to attain better
!       performance in any particular computational environment.
!
!       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!                        Default: 75. (Must be at least 11.)
!
!       IPARMQ(ISPEC=13) Recommended deflation window size.
!                        This depends on ILO, IHI and NS, the
!                        number of simultaneous shifts returned
!                        by IPARMQ(ISPEC=15).  The default for
!                        (IHI-ILO+1) <= 500 is NS.  The default
!                        for (IHI-ILO+1) > 500 is 3*NS/2.
!
!       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!
!       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!                        a multi-shift QR iteration.
!
!                        If IHI-ILO+1 is ...
!
!                        greater than      ...but less    ... the
!                        or equal to ...      than        default is
!
!                                0               30       NS =   2+
!                               30               60       NS =   4+
!                               60              150       NS =  10
!                              150              590       NS =  **
!                              590             3000       NS =  64
!                             3000             6000       NS = 128
!                             6000             infinity   NS = 256
!
!                    (+)  By default matrices of this order are
!                         passed to the implicit double shift routine
!                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!                         values of NS are used only in case of a rare
!                         xLAHQR failure.
!
!                    (**) The asterisks (**) indicate an ad-hoc
!                         function increasing from 10 to 64.
!
!       IPARMQ(ISPEC=16) Select structured matrix multiply.
!                        (See ISPEC=16 above for details.)
!                        Default: 3.
!
!=====================================================================

  INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )

    IMPORT,  NONE
    IMPLICIT NONE (EXTERNAL, TYPE)

    !  -- LAPACK auxiliary routine (version 3.4.0) --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     November 2011

    ! .. Scalar Arguments ..
    INTEGER,      INTENT(IN) :: IHI, ILO, ISPEC, LWORK, N
    CHARACTER(*), INTENT(IN) :: NAME, OPTS

    ! .. Parameters ..
    INTEGER, PARAMETER :: INMIN  = 12, INWIN  = 13,  INIBL  = 14,  &
                          ISHFTS = 15, IACC22 = 16
    INTEGER, PARAMETER :: NMIN   = 75, K22MIN = 14,  KACMIN = 14,  &
                          NIBBLE = 14, KNWSWP = 500
    REAL,    PARAMETER :: TWO = 2.0

    ! .. Local Scalars ..
    INTEGER            :: NH, NS

    ! .. Executable Statements ..
    IF ( ( ISPEC == ISHFTS ) .OR. ( ISPEC == INWIN ) .OR. &
         ( ISPEC == IACC22 ) ) THEN

       !==== Set the number simultaneous shifts ====

       NH = IHI - ILO + 1
       NS = 2
       IF ( NH >= 30 ) &
            NS = 4
       IF ( NH >= 60 ) &
            NS = 10
       IF ( NH >= 150 ) &
            NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
       IF ( NH >= 590 ) &
            NS = 64
       IF ( NH >= 3000 ) &
            NS = 128
       IF ( NH >= 6000 ) &
            NS = 256
       NS = MAX( 2, NS-MOD( NS, 2 ) )
    END IF

    IF ( ISPEC == INMIN ) THEN

       !===== Matrices of order smaller than NMIN get sent
       !        .     to xLAHQR, the classic double shift algorithm.
       !        .     This must be at least 11. ====

       IPARMQ = NMIN

    ELSE IF ( ISPEC == INIBL ) THEN

       !==== INIBL: skip a multi-shift qr iteration and
       !        .    whenever aggressive early deflation finds
       !        .    at least (NIBBLE*(window size)/100) deflations. ====

       IPARMQ = NIBBLE

    ELSE IF ( ISPEC == ISHFTS ) THEN

       !==== NSHFTS: The number of simultaneous shifts =====

       IPARMQ = NS

    ELSE IF ( ISPEC == INWIN ) THEN

       !==== NW: deflation window size.  ====

       IF ( NH <= KNWSWP ) THEN
          IPARMQ = NS
       ELSE
          IPARMQ = 3*NS / 2
       END IF

    ELSE IF ( ISPEC == IACC22 ) THEN

       !==== IACC22: Whether to accumulate reflections
       !        .     before updating the far-from-diagonal elements
       !        .     and whether to use 2-by-2 block structure while
       !        .     doing it.  A small amount of work could be saved
       !        .     by making this choice dependent also upon the
       !        .     NH=IHI-ILO+1.

       IPARMQ = 0
       IF ( NS >= KACMIN ) IPARMQ = 1
       IF ( NS >= K22MIN ) IPARMQ = 2

    ELSE

       !===== invalid value of ispec =====
       IPARMQ = -1

    END IF

    ! Unused arguments (removes compiler warnings)
    IF (.FALSE.) THEN
       WRITE(*,*) N, LWORK, NAME, OPTS, " ARE UNUSED"
    END IF

  END FUNCTION IPARMQ

End Module Lapack_Utils
