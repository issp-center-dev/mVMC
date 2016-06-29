      SUBROUTINE DSKTD2( UPLO, MODE, N, A, LDA, E, TAU, INFO )
*
*  -- Written on 10/22/2010
*     Michael Wimmer, Universiteit Leiden
*
*  Based on  the LAPACK routine ZHETD2 (www.netlib.org)
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  DSKTD2 reduces a real skew-symmetric matrix A to skew-symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q^T * A * Q = T. Alternatively, the routine can also compute
*  a partial tridiagonal form useful for computing the Pfaffian.
*
*  This routine uses unblocked code.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          skew-symmetric matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  MODE    (input) CHARACTER*1
*          = 'N':  A is fully tridiagonalized
*          = 'P':  A is partially tridiagonalized for Pfaffian computation
*                  (details see below)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the skew-symmetric matrix A.
*            If UPLO = 'U', the leading N-by-N upper triangular part
*              of A contains the upper triangular part of the matrix A,
*              and the strictly lower triangular part of A is not referenced.
*            If UPLO = 'L', the leading N-by-N lower triangular part
*              of A contains the lower triangular part of the matrix A,
*              and the strictly upper triangular part of A is not referenced.
*          On exit, if MODE = 'N':
*            If UPLO = 'U', the diagonal and first superdiagonal
*              of A are overwritten by the corresponding elements of the
*              tridiagonal matrix T, and the elements above the first
*              superdiagonal, with the array TAU, represent the unitary
*              matrix Q as a product of elementary reflectors;
*            If UPLO = 'L', the diagonal and first subdiagonal of A are over-
*              written by the corresponding elements of the tridiagonal
*              matrix T, and the elements below the first subdiagonal, with
*              the array TAU, represent the unitary matrix Q as a product
*              of elementary reflectors.
*            See Further Details, also for information about MODE = 'P'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*          If MODE = 'P', only the entries at i odd are well-defined
*          (see Further Details)
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The normal use for DSKTD2 is to compute the tridiagonal form of
*  a skew-symmetric matrix under an orthogonal similarity transformation,
*  and chosen by setting MODE = 'N' ("normal" mode). The other
*  use of DSKTD2 is the computation the Pfaffian of a skew-symmetric matrix,
*  which only requires a partial tridiagonalization, this mode is chosen
*  by setting MODE = 'P' ("Pfaffian" mode).
*
*  Normal mode (MODE = 'N'):
*  ========================
*
*  The routine computes a tridiagonal matrix T and an orthogonal Q such
*  that A = Q * T * Q^T .
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  0   e   v2  v3  v4 )              (  0                  )
*    (      0   e   v3  v4 )              (  e   0              )
*    (          0   e   v4 )              (  v1  e   0          )
*    (              0   e  )              (  v1  v2  e   0      )
*    (                  0  )              (  v1  v2  v3  e   0  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  The LAPACK routine DORGTR can be used to form the transformation
*  matrix explicitely, and DORMTR can be used to multiply another
*  matrix without forming the transformation.
*
*  Pfaffian mode (MODE = 'P'):
*  ==========================
*
*  For computing the Pfaffian, it is enough to bring A into a partial
*  tridiagonal form. In particular, assuming n even, it is enough to
*  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
*  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
*  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
*  only the off-diagonal entries in the odd columns (if UPLO = 'L')
*  or in the even columns (if UPLU = 'U') are properly computed by DSKTD2.
*
*  A is brought into this special form pT using an orthogonal matrix Q:
*  A = Q * pT * Q^T
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) H(n-3) . . . H(3) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v^T
*
*  where tau is a real scalar, and v is a real vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(3) . . . H(n-3) H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v^T
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 6:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  0   e   x   v3  x   v5 )        (  0                      )
*    (      0   x   v3  x   v5 )        (  e   0                  )
*    (          0   e   x   v5 )        (  v1  x   0              )
*    (              0   x   v5 )        (  v1  x   e   0          )
*    (                  0   e  )        (  v1  x   v3  x   0      )
*    (                      0  )        (  v1  x   v3  x   e   0  )
*
*  where d and e denote diagonal and off-diagonal elements of T, vi
*  denotes an element of the vector defining H(i), and x denotes an
*  element not computed by DSKTD2.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0,
     $                   ZERO = 0.0D+0)
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER, NORMAL
      INTEGER            I, STEP
      DOUBLE PRECISION   ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DSKMV, DSKR2, DLARFG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NORMAL = LSAME( MODE, 'N' )

      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NORMAL .AND. .NOT.LSAME( MODE, 'P' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( .NOT.NORMAL .AND. MOD(N,2).NE.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKTD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( .NOT. NORMAL ) THEN
         STEP = 2
      ELSE
         STEP = 1
      END IF

      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         A( N, N ) = ZERO
         DO 10 I = N - 1, 1, -STEP
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            ALPHA = A( I, I+1 )
            CALL DLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i-step+1,1:i-step+1)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
*
               CALL DSKMV( UPLO, I, TAUI, A, LDA,
     $                     A( 1, I+1 ),1, ZERO, TAU, 1 )

*
*              Apply the transformation as a rank-2 update:
*                 A := A + v * w^T - w * v^T
*
               CALL DSKR2( UPLO, I-STEP+1, ONE, A( 1, I+1 ), 1, TAU, 1,
     $                     A, LDA )
*
            ELSE
               A( I, I ) = ZERO
            END IF
            A( I, I+1 ) = E( I )
            TAU( I ) = TAUI
   10    CONTINUE
      ELSE
*
*        Reduce the lower triangle of A
*
         A( 1, 1 ) = ZERO
         DO 20 I = 1, N - 1, STEP
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            ALPHA = A( I+1, I )
            CALL DLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+step:n,i+step:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau^* * A * v^*  storing y in TAU(i:n-1)
*
               CALL DSKMV( UPLO, N-I, TAUI,
     $                     A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A + v * x^T - x * v^T
*
               CALL DSKR2( UPLO, N-I-STEP+1, ONE, A( I+STEP, I ), 1,
     $                     TAU( I+STEP-1 ), 1,
     $                     A( I+STEP, I+STEP ), LDA )
*
            ELSE
               A( I+1, I+1 ) = ZERO
            END IF
            A( I+1, I ) = E( I )
            TAU( I ) = TAUI
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DSKTD2
*
      END
