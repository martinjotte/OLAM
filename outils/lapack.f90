!> \brief \b IEEECK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download IEEECK + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
! 
!       .. Scalar Arguments ..
!       INTEGER            ISPEC
!       REAL               ONE, ZERO
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> IEEECK is called from the ILAENV to verify that Infinity and
!> possibly NaN arithmetic is safe (i.e. will not trap).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies whether to test just for inifinity arithmetic
!>          or whether to test for infinity and NaN arithmetic.
!>          = 0: Verify infinity arithmetic only.
!>          = 1: Verify infinity and NaN arithmetic.
!> \endverbatim
!>
!> \param[in] ZERO
!> \verbatim
!>          ZERO is REAL
!>          Must contain the value 0.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!> \endverbatim
!>
!> \param[in] ONE
!> \verbatim
!>          ONE is REAL
!>          Must contain the value 1.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!>
!>  RETURN VALUE:  INTEGER
!>          = 0:  Arithmetic failed to produce the correct answers
!>          = 1:  Arithmetic produced the correct answers
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                         NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 ) &
         RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*ZERO
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END
!> \brief \b ILAENV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILAENV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
! 
!       .. Scalar Arguments ..
!       CHARACTER*( * )    NAME, OPTS
!       INTEGER            ISPEC, N1, N2, N3, N4
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAENV is called from the LAPACK routines to choose problem-dependent
!> parameters for the local environment.  See ISPEC for a description of
!> the parameters.
!>
!> ILAENV returns an INTEGER
!> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
!> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!>
!> This version provides a set of parameters which should give good,
!> but not optimal, performance on many of the currently available
!> computers.  Users are encouraged to modify this subroutine to set
!> the tuning parameters for their particular machine using the option
!> and problem size information in the arguments.
!>
!> This routine will not function correctly if it is converted to all
!> lower case.  Converting it to all upper case is allowed.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies the parameter to be returned as the value of
!>          ILAENV.
!>          = 1: the optimal blocksize; if this value is 1, an unblocked
!>               algorithm will give the best performance.
!>          = 2: the minimum block size for which the block routine
!>               should be used; if the usable block size is less than
!>               this value, an unblocked routine should be used.
!>          = 3: the crossover point (in a block routine, for N less
!>               than this value, an unblocked routine should be used)
!>          = 4: the number of shifts, used in the nonsymmetric
!>               eigenvalue routines (DEPRECATED)
!>          = 5: the minimum column dimension for blocking to be used;
!>               rectangular blocks must have dimension at least k by m,
!>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!>          = 6: the crossover point for the SVD (when reducing an m by n
!>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!>               this value, a QR factorization is used first to reduce
!>               the matrix to a triangular form.)
!>          = 7: the number of processors
!>          = 8: the crossover point for the multishift QR method
!>               for nonsymmetric eigenvalue problems (DEPRECATED)
!>          = 9: maximum size of the subproblems at the bottom of the
!>               computation tree in the divide-and-conquer algorithm
!>               (used by xGELSD and xGESDD)
!>          =10: ieee NaN arithmetic can be trusted not to trap
!>          =11: infinity arithmetic can be trusted not to trap
!>          12 <= ISPEC <= 16:
!>               xHSEQR or one of its subroutines,
!>               see IPARMQ for detailed explanation
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is CHARACTER*(*)
!>          The name of the calling subroutine, in either upper case or
!>          lower case.
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is CHARACTER*(*)
!>          The character options to the subroutine NAME, concatenated
!>          into a single character string.  For example, UPLO = 'U',
!>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!>          be specified as OPTS = 'UTN'.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!> \endverbatim
!>
!> \param[in] N3
!> \verbatim
!>          N3 is INTEGER
!> \endverbatim
!>
!> \param[in] N4
!> \verbatim
!>          N4 is INTEGER
!>          Problem dimensions for the subroutine NAME; these may not all
!>          be required.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following conventions have been used when calling ILAENV from the
!>  LAPACK routines:
!>  1)  OPTS is a concatenation of all of the character options to
!>      subroutine NAME, in the same order that they appear in the
!>      argument list for NAME, even if they are not used in determining
!>      the value of the parameter specified by ISPEC.
!>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!>      that they appear in the argument list for NAME.  N1 is used
!>      first, N2 second, and so on, and unused problem dimensions are
!>      passed a value of -1.
!>  3)  The parameter value returned by ILAENV is checked for validity in
!>      the calling subroutine.  For example, ILAENV is used to retrieve
!>      the optimal blocksize for STRTRI as follows:
!>
!>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!>      IF( NB.LE.1 ) NB = MAX( 1, N )
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
!     ..
!     .. Executable Statements ..
!
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
              130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
   10 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 ) &
                  SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
             ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                   ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                   ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: &
                   I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 ) &
                  SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
!
      GO TO ( 50, 60, 70 )ISPEC
!
   50 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
                  C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
   60 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
             'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
   70 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ. &
             'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. &
                'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) &
                 THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
   80 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
   90 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  100 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  110 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  120 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
  130 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
      ILAENV = 25
      RETURN
!
  140 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
!
  150 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
!
  160 CONTINUE
!
!     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
!
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
!
!     End of ILAENV
!
      END
!> \brief \b ILASLC scans a matrix for its last non-zero column.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILASLC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
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
!       REAL               A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILASLC scans A for its last non-zero column.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILASLC( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL             ZERO
      PARAMETER ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER I
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILASLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILASLC = N
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO ILASLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILASLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
!> \brief \b ILASLR scans a matrix for its last non-zero row.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILASLR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaslr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaslr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaslr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILASLR( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       INTEGER            M, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILASLR scans A for its last non-zero row.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date September 2012
!
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION ILASLR( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     September 2012
!
!     .. Scalar Arguments ..
      INTEGER            M, N, LDA
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL             ZERO
      PARAMETER ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILASLR = M
      ELSEIF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILASLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILASLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILASLR = MAX( ILASLR, I )
         END DO
      END IF
      RETURN
      END
!> \brief \b IPARMQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download IPARMQ + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
! 
!       .. Scalar Arguments ..
!       INTEGER            IHI, ILO, ISPEC, LWORK, N
!       CHARACTER          NAME*( * ), OPTS*( * )
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      This program sets problem and machine dependent parameters
!>      useful for xHSEQR and its subroutines. It is called whenever 
!>      ILAENV is called with 12 <= ISPEC <= 16
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is integer scalar
!>              ISPEC specifies which tunable parameter IPARMQ should
!>              return.
!>
!>              ISPEC=12: (INMIN)  Matrices of order nmin or less
!>                        are sent directly to xLAHQR, the implicit
!>                        double shift QR algorithm.  NMIN must be
!>                        at least 11.
!>
!>              ISPEC=13: (INWIN)  Size of the deflation window.
!>                        This is best set greater than or equal to
!>                        the number of simultaneous shifts NS.
!>                        Larger matrices benefit from larger deflation
!>                        windows.
!>
!>              ISPEC=14: (INIBL) Determines when to stop nibbling and
!>                        invest in an (expensive) multi-shift QR sweep.
!>                        If the aggressive early deflation subroutine
!>                        finds LD converged eigenvalues from an order
!>                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
!>                        then the next QR sweep is skipped and early
!>                        deflation is applied immediately to the
!>                        remaining active diagonal block.  Setting
!>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
!>                        multi-shift QR sweep whenever early deflation
!>                        finds a converged eigenvalue.  Setting
!>                        IPARMQ(ISPEC=14) greater than or equal to 100
!>                        prevents TTQRE from skipping a multi-shift
!>                        QR sweep.
!>
!>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
!>                        a multi-shift QR iteration.
!>
!>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
!>                        following meanings.
!>                        0:  During the multi-shift QR sweep,
!>                            xLAQR5 does not accumulate reflections and
!>                            does not use matrix-matrix multiply to
!>                            update the far-from-diagonal matrix
!>                            entries.
!>                        1:  During the multi-shift QR sweep,
!>                            xLAQR5 and/or xLAQRaccumulates reflections and uses
!>                            matrix-matrix multiply to update the
!>                            far-from-diagonal matrix entries.
!>                        2:  During the multi-shift QR sweep.
!>                            xLAQR5 accumulates reflections and takes
!>                            advantage of 2-by-2 block structure during
!>                            matrix-matrix multiplies.
!>                        (If xTRMM is slower than xGEMM, then
!>                        IPARMQ(ISPEC=16)=1 may be more efficient than
!>                        IPARMQ(ISPEC=16)=2 despite the greater level of
!>                        arithmetic work implied by the latter choice.)
!> \endverbatim
!>
!> \param[in] NAME
!> \verbatim
!>          NAME is character string
!>               Name of the calling subroutine
!> \endverbatim
!>
!> \param[in] OPTS
!> \verbatim
!>          OPTS is character string
!>               This is a concatenation of the string arguments to
!>               TTQRE.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer scalar
!>               N is the order of the Hessenberg matrix H.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>               It is assumed that H is already upper triangular
!>               in rows and columns 1:ILO-1 and IHI+1:N.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer scalar
!>               The amount of workspace available.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>       Little is known about how best to choose these parameters.
!>       It is possible to use different values of the parameters
!>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
!>
!>       It is probably best to choose different parameters for
!>       different matrices and different parameters at different
!>       times during the iteration, but this has not been
!>       implemented --- yet.
!>
!>
!>       The best choices of most of the parameters depend
!>       in an ill-understood way on the relative execution
!>       rate of xLAQR3 and xLAQR5 and on the nature of each
!>       particular eigenvalue problem.  Experiment may be the
!>       only practical way to determine which choices are most
!>       effective.
!>
!>       Following is a list of default values supplied by IPARMQ.
!>       These defaults may be adjusted in order to attain better
!>       performance in any particular computational environment.
!>
!>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
!>                        Default: 75. (Must be at least 11.)
!>
!>       IPARMQ(ISPEC=13) Recommended deflation window size.
!>                        This depends on ILO, IHI and NS, the
!>                        number of simultaneous shifts returned
!>                        by IPARMQ(ISPEC=15).  The default for
!>                        (IHI-ILO+1).LE.500 is NS.  The default
!>                        for (IHI-ILO+1).GT.500 is 3*NS/2.
!>
!>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
!>
!>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
!>                        a multi-shift QR iteration.
!>
!>                        If IHI-ILO+1 is ...
!>
!>                        greater than      ...but less    ... the
!>                        or equal to ...      than        default is
!>
!>                                0               30       NS =   2+
!>                               30               60       NS =   4+
!>                               60              150       NS =  10
!>                              150              590       NS =  **
!>                              590             3000       NS =  64
!>                             3000             6000       NS = 128
!>                             6000             infinity   NS = 256
!>
!>                    (+)  By default matrices of this order are
!>                         passed to the implicit double shift routine
!>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
!>                         values of NS are used only in case of a rare
!>                         xLAHQR failure.
!>
!>                    (**) The asterisks (**) indicate an ad-hoc
!>                         function increasing from 10 to 64.
!>
!>       IPARMQ(ISPEC=16) Select structured matrix multiply.
!>                        (See ISPEC=16 above for details.)
!>                        Default: 3.
!> \endverbatim
!>
!  =====================================================================
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
!
!  ================================================================
!     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14, &
                         ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, &
                         NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
!     ..
!     .. Local Scalars ..
      INTEGER            NH, NS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
!     ..
!     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. &
          ( ISPEC.EQ.IACC22 ) ) THEN
!
!        ==== Set the number simultaneous shifts ====
!
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 ) &
            NS = 4
         IF( NH.GE.60 ) &
            NS = 10
         IF( NH.GE.150 ) &
            NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 ) &
            NS = 64
         IF( NH.GE.3000 ) &
            NS = 128
         IF( NH.GE.6000 ) &
            NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
!
      IF( ISPEC.EQ.INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to xLAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      ELSE IF( ISPEC.EQ.INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = NS
!
      ELSE IF( ISPEC.EQ.INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
!
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
         IPARMQ = 0
         IF( NS.GE.KACMIN ) &
            IPARMQ = 1
         IF( NS.GE.K22MIN ) &
            IPARMQ = 2
!
      ELSE
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      END IF
!
!     ==== End of IPARMQ ====
!
      END
!> \brief \b LSAME
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!      LOGICAL FUNCTION LSAME( CA, CB )
!
!     .. Scalar Arguments ..
!      CHARACTER          CA, CB
!     ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAME returns .TRUE. if CA is the same letter as CB regardless of
!> case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CA
!> \verbatim
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CA and CB specify the single characters to be compared.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME ) &
         RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
             INTA.GE.145 .AND. INTA.LE.153 .OR. &
             INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
             INTB.GE.145 .AND. INTB.LE.153 .OR. &
             INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END
!> \brief \b SLAMCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!      REAL             FUNCTION SLAMCH( CMACH )
!
!     .. Scalar Arguments ..
!      CHARACTER          CMACH
!     ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAMCH determines single precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          Specifies the value to be returned by SLAMCH:
!>          = 'E' or 'e',   SLAMCH := eps
!>          = 'S' or 's ,   SLAMCH := sfmin
!>          = 'B' or 'b',   SLAMCH := base
!>          = 'P' or 'p',   SLAMCH := eps*base
!>          = 'N' or 'n',   SLAMCH := t
!>          = 'R' or 'r',   SLAMCH := rnd
!>          = 'M' or 'm',   SLAMCH := emin
!>          = 'U' or 'u',   SLAMCH := rmin
!>          = 'L' or 'l',   SLAMCH := emax
!>          = 'O' or 'o',   SLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      REAL             FUNCTION SLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      REAL               RND, EPS, SFMIN, SMALL, RMACH
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, &
                         MINEXPONENT, RADIX, TINY
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
      RND = ONE
!
      IF( ONE.EQ.RND ) THEN
         EPS = EPSILON(ZERO) * 0.5
      ELSE
         EPS = EPSILON(ZERO)
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         SFMIN = TINY(ZERO)
         SMALL = ONE / HUGE(ZERO)
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = EPS * RADIX(ZERO)
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = DIGITS(ZERO)
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = MINEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = tiny(zero)
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = MAXEXPONENT(ZERO)
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = HUGE(ZERO)
      ELSE
         RMACH = ZERO
      END IF
!
      SLAMCH = RMACH
      RETURN
!
!     End of SLAMCH
!
      END
!***********************************************************************
!> \brief \b SLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date November 2011
!> \ingroup auxOTHERauxiliary
!>
!> \param[in] A
!> \verbatim
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          The values A and B.
!> \endverbatim
!>
!
      REAL             FUNCTION SLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      REAL               A, B
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B
!
      RETURN
!
!     End of SLAMC3
!
      END
!
!***********************************************************************
!> \brief \b XERBLA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download XERBLA + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/xerbla.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/xerbla.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/xerbla.f"> 
!> [TXT]</a>
!> \endhtmlonly 
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
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> XERBLA  is an error handler for the LAPACK routines.
!> It is called by an LAPACK routine if an input parameter has an
!> invalid value.  A message is printed and execution stops.
!>
!> Installers may consider modifying the STOP statement in order to
!> call system-specific exception-handling facilities.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the routine which called XERBLA.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The position of the invalid parameter in the parameter list
!>          of the calling routine.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
!     ..
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
!     ..
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
      END
      SUBROUTINE SGELQ2( M, N, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      REAL               AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELQ2', -INFO )
         RETURN
      END IF
!
      K = MIN( M, N )
!
      DO 10 I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
!
         CALL SLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, &
                      TAU( I ) )
         IF( I.LT.M ) THEN
!
!           Apply H(i) to A(i+1:m,i:n) from the right
!
            AII = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), &
                        A( I+1, I ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
!
!     End of SGELQ2
!
      END
      SUBROUTINE SGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, &
                         NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGELQ2, SLARFB, SLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'SGELQF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SGELQF', ' ', M, N, -1, &
                       -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!
!           Compute the LQ factorization of the current block
!           A(i:i+ib-1,i:n)
!
            CALL SGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
            IF( I+IB.LE.M ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL SLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ), &
                            LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i+ib:m,i:n) from the right
!
               CALL SLARFB( 'Right', 'No transpose', 'Forward', &
                            'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), &
                            LDA, WORK, LDWORK, A( I+IB, I ), LDA, &
                            WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K ) &
         CALL SGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                      IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SGELQF
!
      END
      SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, &
                        INFO )
!
!  -- LAPACK driver routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
!     ..
!     .. Local Arrays ..
      REAL               RWORK( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANGE
      EXTERNAL           LSAME, ILAENV, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGELQF, SGEQRF, SLABAD, SLASCL, SLASET, SORMLQ, &
                         SORMQR, STRTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'T' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN + MAX( MN, NRHS ) ) .AND. &
         .NOT.LQUERY ) THEN
         INFO = -10
      END IF
!
!     Figure out optimal block size
!
      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN
!
         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) ) &
            TPSD = .FALSE.
!
         IF( M.GE.N ) THEN
            NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LN', M, NRHS, N, &
                    -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LT', M, NRHS, N, &
                    -1 ) )
            END IF
         ELSE
            NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LT', N, NRHS, M, &
                    -1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LN', N, NRHS, M, &
                    -1 ) )
            END IF
         END IF
!
         WSIZE = MAX( 1, MN + MAX( MN, NRHS )*NB )
         WORK( 1 ) = REAL( WSIZE )
!
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL SLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RETURN
      END IF
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
!
!     Scale A, B if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = SLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL SLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         GO TO 50
      END IF
!
      BROW = M
      IF( TPSD ) &
         BROW = N
      BNRM = SLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, &
                      INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, &
                      INFO )
         IBSCL = 2
      END IF
!
      IF( M.GE.N ) THEN
!
!        compute QR factorization of A
!
         CALL SGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, &
                      INFO )
!
!        workspace at least N, optimally N*NB
!
         IF( .NOT.TPSD ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
!
            CALL SORMQR( 'Left', 'Transpose', M, NRHS, N, A, LDA, &
                         WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, &
                         INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
            CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS, &
                         A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
            SCLLEN = N
!
         ELSE
!
!           Overdetermined system of equations A**T * X = B
!
!           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
!
            CALL STRTRS( 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                         A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
!           B(N+1:M,1:NRHS) = ZERO
!
            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
!
            CALL SORMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA, &
                         WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, &
                         INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
            SCLLEN = M
!
         END IF
!
      ELSE
!
!        Compute LQ factorization of A
!
         CALL SGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, &
                      INFO )
!
!        workspace at least M, optimally M*NB.
!
         IF( .NOT.TPSD ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
            CALL STRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS, &
                         A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
!           B(M+1:N,1:NRHS) = 0
!
            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
!
!           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
!
            CALL SORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA, &
                         WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, &
                         INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
            SCLLEN = N
!
         ELSE
!
!           overdetermined system min || A**T * X - B ||
!
!           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
!
            CALL SORMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA, &
                         WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, &
                         INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
!
            CALL STRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS, &
                         A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
            SCLLEN = M
!
         END IF
!
      END IF
!
!     Undo scaling
!
      IF( IASCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, &
                      INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, &
                      INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, &
                      INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, &
                      INFO )
      END IF
!
   50 CONTINUE
      WORK( 1 ) = REAL( WSIZE )
!
      RETURN
!
!     End of SGELS
!
      END
      SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      REAL               AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQR2', -INFO )
         RETURN
      END IF
!
      K = MIN( M, N )
!
      DO 10 I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         CALL SLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                      TAU( I ) )
         IF( I.LT.N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            AII = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
!
!     End of SGEQR2
!
      END
      SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, &
                         NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEQR2, SLARFB, SLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'SGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SGEQRF', ' ', M, N, -1, &
                       -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!
!           Compute the QR factorization of the current block
!           A(i:m,i:i+ib-1)
!
            CALL SGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK, &
                         IINFO )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H**T to A(i:m,i+ib:n) from the left
!
               CALL SLARFB( 'Left', 'Transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K ) &
         CALL SGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK, &
                      IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SGEQRF
!
      END
      LOGICAL FUNCTION SISNAN( SIN )
!
!  -- LAPACK auxiliary routine (version 3.2.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2010
!
!     .. Scalar Arguments ..
      REAL               SIN
!     ..
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
!  =====================================================================
!
!  .. External Functions ..
      LOGICAL SLAISNAN
      EXTERNAL SLAISNAN
!  ..
!  .. Executable Statements ..
      SISNAN = SLAISNAN(SIN,SIN)
      RETURN
      END
      SUBROUTINE SLABAD( SMALL, LARGE )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      REAL               LARGE, SMALL
!     ..
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
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE ).GT.2000. ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of SLABAD
!
      END
      LOGICAL FUNCTION SLAISNAN( SIN1, SIN2 )
!
!  -- LAPACK auxiliary routine (version 3.2.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2010
!
!     .. Scalar Arguments ..
      REAL               SIN1, SIN2
!     ..
!
!  Purpose
!  =======
!
!  This routine is not for general use.  It exists solely to avoid
!  over-optimization in SISNAN.
!
!  SLAISNAN checks for NaNs by comparing its two arguments for
!  inequality.  NaN is the only floating-point value where NaN != NaN
!  returns .TRUE.  To check for NaNs, pass the same variable as both
!  arguments.
!
!  A compiler must assume that the two arguments are
!  not the same variable, and the test will not be optimized away.
!  Interprocedural or whole-program optimization may delete this
!  test.  The ISNAN functions will be replaced by the correct
!  Fortran 03 intrinsic once the intrinsic is widely available.
!
!  Arguments
!  =========
!
!  SIN1     (input) REAL
!
!  SIN2     (input) REAL
!          Two numbers to compare for inequality.
!
!  =====================================================================
!
!  .. Executable Statements ..
      SLAISNAN = (SIN1.NE.SIN2)
      RETURN
      END
      REAL             FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
!     ..
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
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL               SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL SLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      SLANGE = VALUE
      RETURN
!
!     End of SLANGE
!
      END
      REAL             FUNCTION SLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      REAL               X, Y
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      REAL               W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         SLAPY2 = W
      ELSE
         SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of SLAPY2
!
      END
      SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      REAL               TAU
!     ..
!     .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILASLR, ILASLC
      EXTERNAL           LSAME, ILASLR, ILASLC
!     ..
!     .. Executable Statements ..
!
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILASLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILASLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
!
!        Form  H * C
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
!
            CALL SGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
                 ZERO, WORK, 1 )
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
!
            CALL SGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
            CALL SGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, &
                 V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
!
            CALL SGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of SLARF
!
      END
      SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                         T, LDT, C, LDC, WORK, LDWORK )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                         WORK( LDWORK, * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J, LASTV, LASTC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILASLR, ILASLC
      EXTERNAL           LSAME, ILASLR, ILASLC
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMM, STRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) &
         RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILASLR( M, K, V, LDV ) )
               LASTC = ILASLC( LASTV, N, C, LDC )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C1**T
!
               DO 10 J = 1, K
                  CALL SCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2**T *V2
!
                  CALL SGEMM( 'Transpose', 'No transpose', &
                       LASTC, K, LASTV-K, &
                       ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**T
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - V2 * W**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', &
                       LASTV-K, LASTC, K, &
                       -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                       C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**T
!
               DO 30 J = 1, K
                  DO 20 I = 1, LASTC
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILASLR( N, K, V, LDV ) )
               LASTC = ILASLR( M, LASTV, C, LDC )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL SCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL SGEMM( 'No transpose', 'No transpose', &
                       LASTC, K, LASTV-K, &
                       ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**T
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - W * V2**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', &
                       LASTC, LASTV-K, K, &
                       -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                       C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, LASTC
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILASLR( M, K, V, LDV ) )
               LASTC = ILASLC( LASTV, N, C, LDC )
!
!              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
!
!              W := C2**T
!
               DO 70 J = 1, K
                  CALL SCOPY( LASTC, C( LASTV-K+J, 1 ), LDC, &
                       WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1**T*V1
!
                  CALL SGEMM( 'Transpose', 'No transpose', &
                       LASTC, K, LASTV-K, ONE, C, LDC, V, LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**T
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - V1 * W**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', &
                       LASTV-K, LASTC, K, -ONE, V, LDV, WORK, LDWORK, &
                       ONE, C, LDC )
               END IF
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
!
!              C2 := C2 - W**T
!
               DO 90 J = 1, K
                  DO 80 I = 1, LASTC
                     C( LASTV-K+J, I ) = C( LASTV-K+J, I ) - WORK(I, J)
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILASLR( N, K, V, LDV ) )
               LASTC = ILASLR( M, LASTV, C, LDC )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL SCOPY( LASTC, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL SGEMM( 'No transpose', 'No transpose', &
                       LASTC, K, LASTV-K, ONE, C, LDC, V, LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**T
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - W * V1**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', &
                       LASTC, LASTV-K, K, -ONE, WORK, LDWORK, V, LDV, &
                       ONE, C, LDC )
               END IF
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, LASTC
                     C( I, LASTV-K+J ) = C( I, LASTV-K+J ) - WORK(I, J)
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILASLC( K, M, V, LDV ) )
               LASTC = ILASLC( LASTV, N, C, LDC )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C1**T
!
               DO 130 J = 1, K
                  CALL SCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2**T*V2**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', &
                       LASTC, K, LASTV-K, &
                       ONE, C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**T * W**T
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - V2**T * W**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', &
                       LASTV-K, LASTC, K, &
                       -ONE, V( 1, K+1 ), LDV, WORK, LDWORK, &
                       ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**T
!
               DO 150 J = 1, K
                  DO 140 I = 1, LASTC
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILASLC( K, N, V, LDV ) )
               LASTC = ILASLR( M, LASTV, C, LDC )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL SCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2 * V2**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', &
                       LASTC, K, LASTV-K, &
                       ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL SGEMM( 'No transpose', 'No transpose', &
                       LASTC, LASTV-K, K, &
                       -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, &
                       ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, LASTC
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**T * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILASLC( K, M, V, LDV ) )
               LASTC = ILASLC( LASTV, N, C, LDC )
!
!              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
!
!              W := C2**T
!
               DO 190 J = 1, K
                  CALL SCOPY( LASTC, C( LASTV-K+J, 1 ), LDC, &
                       WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1**T * V1**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', &
                       LASTC, K, LASTV-K, ONE, C, LDC, V, LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**T * W**T
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - V1**T * W**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', &
                       LASTV-K, LASTC, K, -ONE, V, LDV, WORK, LDWORK, &
                       ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
!
!              C2 := C2 - W**T
!
               DO 210 J = 1, K
                  DO 200 I = 1, LASTC
                     C( LASTV-K+J, I ) = C( LASTV-K+J, I ) - WORK(I, J)
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILASLC( K, N, V, LDV ) )
               LASTC = ILASLR( M, LASTV, C, LDC )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL SCOPY( LASTC, C( 1, LASTV-K+J ), 1, &
                       WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', &
                    LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1 * V1**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', &
                       LASTC, K, LASTV-K, ONE, C, LDC, V, LDV, &
                       ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', &
                    LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL SGEMM( 'No transpose', 'No transpose', &
                       LASTC, LASTV-K, K, -ONE, WORK, LDWORK, V, LDV, &
                       ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, LASTC
                     C( I, LASTV-K+J ) = C( I, LASTV-K+J ) &
                          - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of SLARFB
!
      END
      SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               ALPHA, TAU
!     ..
!     .. Array Arguments ..
      REAL               X( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      REAL               SLAMCH, SLAPY2, SNRM2
      EXTERNAL           SLAMCH, SLAPY2, SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           SSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = SNRM2( N-1, X, INCX )
!
      IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL SSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = SNRM2( N-1, X, INCX )
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
!
      RETURN
!
!     End of SLARFG
!
      END
      SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
      IMPLICIT NONE
!
!  -- LAPACK auxiliary routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      REAL               T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
      REAL               VII
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, STRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO 20 I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!
!                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
!
                  CALL SGEMV( 'Transpose', J-I+1, I-1, -TAU( I ), &
                              V( I, 1 ), LDV, V( I, I ), 1, ZERO, &
                              T( 1, I ), 1 )
               ELSE
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
!
                  CALL SGEMV( 'No transpose', I-1, J-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                              T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                           LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
   20    CONTINUE
      ELSE
         PREVLASTV = 1
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
!
                     CALL SGEMV( 'Transpose', N-K+I-J+1, K-I, -TAU( I ), &
                                 V( J, I+1 ), LDV, V( J, I ), 1, ZERO, &
                                 T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
!
                     CALL SGEMV( 'No transpose', K-I, N-K+I-J+1, &
                          -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, &
                          ZERO, T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL STRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                              T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of SLARFT
!
      END
      SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 3.3.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2010
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL               CFROM, CTO
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH, SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. SISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( SISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
               ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                  ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
                   THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                  ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                  ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of SLASCL
!
      END
      SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      REAL               ALPHA, BETA
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * )
!     ..
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
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of SLASET
!
      END
      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      REAL               X( * )
!     ..
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
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      REAL               ABSXI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of SLASSQ
!
      END
      SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      REAL               AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i)
!
         AII = A( I, I )
         A( I, I ) = ONE
         CALL SLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), &
                     LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of SORM2R
!
      END
      SUBROUTINE SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      REAL               AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORML2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
         RETURN
!
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i)
!
         AII = A( I, I )
         A( I, I ) = ONE
         CALL SLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ), &
                     C( IC, JC ), LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of SORML2
!
      END
      SUBROUTINE SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), &
                         WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
                         LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      REAL               T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORML2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'SORMLQ', SIDE // TRANS, M, N, K, &
                   -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF 
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORMLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'SORMLQ', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL SLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ), &
                         LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H**T is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H**T is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H**T
!
            CALL SLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB, &
                         A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK, &
                         LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of SORMLQ
!
      END
      SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), &
                         WORK( * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, &
                         LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      REAL               T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORM2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'SORMQR', SIDE // TRANS, M, N, K, &
              -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'SORMQR', SIDE // TRANS, M, N, K, &
                    -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                      IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
             ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL SLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                         LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H**T is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H**T is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H**T
!
            CALL SLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                         IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, &
                         WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of SORMQR
!
      END
      SUBROUTINE STRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, &
                         INFO )
!
!  -- LAPACK routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
!     ..
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. &
               LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STRTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
         RETURN
!
!     Check for singularity.
!
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO ) &
               RETURN
   10    CONTINUE
      END IF
      INFO = 0
!
!     Solve A * x = b  or  A**T * x = b.
!
      CALL STRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, &
                  LDB )
!
      RETURN
!
!     End of STRTRS
!
      END
