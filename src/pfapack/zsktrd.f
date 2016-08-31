      SUBROUTINE ZSKTRD( UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK,
     $                   INFO )
*
*  -- Written on 10/22/2010
*     Michael Wimmer
*
*  -- derived from LAPACK routine ZHETRD (version 3.2, www.netlib.org) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      DOUBLE COMPLEX     A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZSKTRD reduces a complex skew-symmetric matrix A to real skew-symmetric
*  tridiagonal form T by a unitary congruence transformation:
*  Q^H * A * Q^* = T. Alternatively, the routine can also compute
*  a partial tridiagonal form useful for computing the Pfaffian.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  MODE    (input) CHARACTER*1
*          = 'N':  A is fully tridiagonalized
*          = 'P':  A is partially tridiagonalized for Pfaffian computation
*                  (details see below)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
*
*  A       (input/output) DOUBLE COMPLEX array, dimension (LDA,N)
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
*  TAU     (output) DOUBLE COMPLEX array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE COMPLEX array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The normal use for ZSKTRD is to compute the tridiagonal form of
*  a skew-symmetric matrix under an orthogonal similarity transformation,
*  and chosen by setting MODE = 'N' ("normal" mode). The other
*  use of ZSKTRD is the computation the Pfaffian of a skew-symmetric matrix,
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
*  where tau is a complex scalar, and v is a complex vector with
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
*  where tau is a complex scalar, and v is a complex vector with
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
*  The LAPACK routine ZUNGTR can be used to form the transformation
*  matrix explicitely, and ZUNMTR can be used to multiply another
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
*  or in the even columns (if UPLU = 'U') are properly computed by ZSKTRD.
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
*  where tau is a complex scalar, and v is a complex vector with
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
*  where tau is a complex scalar, and v is a complex vector with
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
*  element not computed by ZSKTRD.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE COMPLEX     CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, NORMAL
      INTEGER            I, IINFO, IWS, J, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX, STEP, NPANEL, NXPANEL
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSKR2K, ZSKTD2, ZLASKTRD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NORMAL = LSAME( MODE, 'N' )

      LQUERY = ( LWORK.EQ.-1 )
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
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.
*        NOTE: I keep 'ZHETRD" here, as of course ILAENV has no information
*              about ZSKTRD ...

         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSKTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         NX = MAX( NB, ILAENV( 3, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'ZHETRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN .OR. NB.LE.1 )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*

      IF( .NOT.NORMAL ) THEN
         STEP = 2
      ELSE
         STEP = 1
      END IF

      NPANEL = NB * STEP
      NXPANEL = NX * STEP

      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         DO 20 I = N, NXPANEL + NPANEL, -NPANEL
*
*           Reduce columns i-npanel+1:i to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL ZLASKTRD( UPLO, MODE, I, NB, A, LDA, E, TAU, WORK,
     $                     LDWORK )
*
*           Update the unreduced submatrix A(1:i-npanel,1:i-npanel), using an
*           update of the form:  A := A + V*W^T - W*V^T
*
            CALL ZSKR2K( UPLO, 'No transpose', I-NPANEL, NB, CONE,
     $                   A( 1, I-NPANEL+STEP ), LDA*STEP,
     $                   WORK, LDWORK, CONE, A, LDA )
*
*           Copy superdiagonal elements back into A
*
            DO 10 J = I-NPANEL+1+STEP-1, I, STEP
               A( J-1, J ) = E( J-1 )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL ZSKTD2( UPLO, MODE, I, A, LDA, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NXPANEL, NPANEL
*
*           Reduce columns i:i+npanel-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL ZLASKTRD( UPLO, MODE, N-I+1, NB, A( I, I ), LDA,
     $                     E( I ), TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+npanel:n,i+npanel:n), using
*           an update of the form:  A := A + V*W^T - W*V^T
*
            CALL ZSKR2K( UPLO, 'No transpose', N-I-NPANEL+1, NB, CONE,
     $                   A( I+NPANEL, I ), LDA*STEP,
     $                   WORK( NPANEL+1 ), LDWORK, CONE,
     $                   A( I+NPANEL, I+NPANEL ), LDA )
*
*           Copy subdiagonal elements back into A
*
            DO 30 J = I, I + NPANEL - 1, STEP
               A( J+1, J ) = E( J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL ZSKTD2( UPLO, MODE, N-I+1, A( I, I ), LDA, E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZSKTRD
*
      END
