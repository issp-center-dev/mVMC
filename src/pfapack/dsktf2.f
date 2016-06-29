      SUBROUTINE DSKTF2( UPLO, MODE, N, A, LDA, IPIV, INFO )
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*
*
*  Purpose
*  =======
*
*  DSKTF2 computes the factorization of a skew-symmetric matrix A
*  using the Parlett-Reid algorithm:
*
*     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
*
*  where U (or L) unit upper (lower) triangular matrix (^T denotes
*  the transpose), T is a skew-symmetric tridiagonal matrix and P
*  is a permutation matrix. In addition to being unit triangular,
*  U(1:n-1,n)=0 and L(2:n,1)=0.
*  Instead of a full tridiagonalization, DFSKTF2 can also compute a
*  partial tridiagonal form for computing the Pfaffian.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored:
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
*          On entry, the symmetric matrix A.
*             If UPLO = 'U', the leading n-by-n upper triangular part
*                of A contains the upper triangular part of the matrix A,
*                and the strictly lower triangular part of A is not referenced.
*             If UPLO = 'L', the leading n-by-n lower triangular part
*                of A contains the lower triangular part of the matrix A,
*                and the strictly upper triangular part of A is not referenced.
*
*          On exit, the tridiagonal matrix T and the multipliers used
*          to obtain the factor U or L (see below for further details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Information about the permutation matrix P: row and column
*          i are interchanged with IPIV(i). If UPLO = 'U', those
*          interchanges are done in the order i = N ... 1, if UPLO = 'L'
*          in the order i = 1 ... N.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the off-diagonal entry in the k-th row
*                            (UPLO = 'U') or k-th column (UPLO = 'L')
*                            is exactly zero.
*
*  Further Details
*  ===============
*
*
*  The normal use for DSKTD2 is to compute the U T U^T or L T L^T
*  decomposition of a skew-symmetric matrix with pivoting. This mode
*  is chosen by setting MODE = 'N' ("normal" mode). The other
*  use of DSKTD2 is the computation the Pfaffian of a skew-symmetric matrix,
*  which only requires a partial computation of T, this mode is chosen
*  by setting MODE = 'P' ("Pfaffian" mode).
*
*  Normal mode (MODE = 'N'):
*  ========================
*
*  If UPLO = 'U', the U*T*U^T decomposition of A is computed. U is a
*  upper triangular unit matrix with the additional constraint
*  U(1:n-1,n) = 0, and T a tridiagonal matrix. The upper diagonal
*  of T is stored on exit in A(i,i+1) for i = 1 .. n-1. The column
*  U(1:i-1, i) is stored in A(1:i-1,i+1).
*
*  If UPLO = 'L', the L*T*L^T decomposition of A is computed. L is a
*  lower triangular unit matrix with the additional constraint
*  L(2:n,1) = 0, and T a tridiagonal matrix. The lower diagonal
*  of T is stored on exit in A(i+1,i) for i = 1 .. n-1. The column
*  L(i+1:n, i) is stored in A(i+1:n,i-1).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  0   e   u2  u3  u4 )              (  0                  )
*    (      0   e   u3  u4 )              (  e   0              )
*    (          0   e   u4 )              (  l2  e   0          )
*    (              0   e  )              (  l2  l3  e   0      )
*    (                  0  )              (  l2  l3  l4  e   0  )
*
*  where e denotes the off-diagonal elements of T, and ui (li)
*  denotes an element of the i-th column of U (L).
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
*  or in the even columns (if UPLU = 'U') are properly computed by DSKTF2.
*
*  If UPLO = 'U', the U*pT*U^T decomposition of A is computed. U is a
*  upper triangular unit matrix with the additional constraint
*  U(1:i-1,i) = 0 for even i, and pT a partially tridiagonal matrix.
*  The entries in the odd rows of the upper diagonal of pT are stored
*  on exit in A(i,i+1) for i odd. The column U(1:i-1, i) for odd i
*  is stored in A(1:i-1,i+1).
*
*  If UPLO = 'L', the L*pT*L^T decomposition of A is computed. L is a
*  lower triangular unit matrix with the additional constraint
*  L(i+1:n,i) = 0 for odd i, and pT a partially tridiagonal matrix.
*  The entries in odd columns in the lower diagonal of pT are stored
*  on exit in A(i+1,i) for i odd. The column L(i+1:n, i) for i odd
*  is stored in A(i+1:n,i-1).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 6:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  0   e   x   u3  x   u5 )              (  0                    )
*    (      0   x   u3  x   u5 )              (  e   0                )
*    (          0   e   x   u5 )              (  l2  x   0            )
*    (              0   x   u5 )              (  l2  x   e   0        )
*    (                  0   e  )              (  l2  x   l4  x   0    )
*    (                      0  )              (  l2  x   l4  x   e  0 )
*
*  where e denotes the off-diagonal elements of T, ui (li)
*  denotes an element of the i-th column of U (L), and x denotes an
*  element not computed by DSKTF2.
*
*  =====================================================================
*

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     .. Local Scalars ..
      LOGICAL            UPPER, NORMAL
      INTEGER            K, KK, KP, STEP
      DOUBLE PRECISION   COLMAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSKR2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
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
      ELSE IF( .NOT.NORMAL .AND. MOD(N,2).EQ.1 ) THEN
*     If STEP == 2, we need an even-dimensional matrix
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKTF2', -INFO )
         RETURN
      END IF

*     Quick return if possible
      IF( N .EQ. 0 ) RETURN

      IF( NORMAL ) THEN
         STEP = 1
      ELSE
         STEP = 2
      END IF

      IF( UPPER ) THEN
*     Factorize A as U * T * U^T using the upper triangle of A
         IPIV( N ) = N

         DO 10 K=N, 2, -1
*     Either all columns or only the even ones (MODE = 'P')
            IF( MOD(K, STEP) .EQ. 0) THEN
*     Find the pivot
               KP = IDAMAX(K-1, A( 1, K ), 1)
               COLMAX = ABS( A( KP, K ) )

               IF( COLMAX.EQ.ZERO ) THEN
*     The column is completely zero - do nothing
                  IF( INFO.EQ.0 ) THEN
                     INFO = K
                  END IF
                  KP = K-1
               END IF

*     swap rows and columns K+1 and IMAX in the
*     full sub-matrix A(1:N,1:N)
               KK = K-1

               IF( KP .NE. KK ) THEN
                  CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1)
                  CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1,
     $                 A( KP, KP+1 ), LDA )

                  CALL DSWAP( N-K+1, A( KK, K), LDA, A( KP, K), LDA)

                  CALL DSCAL(KK-KP, -ONE, A(KP, KK), 1)
                  CALL DSCAL(KK-KP-1, -ONE, A(KP, KP+1), LDA)
               END IF

*     Update the leading submatrix A(1:K-2, 1:K-2) in a rank 2 update
*     (The column/row K-1 is not affected by the update)
               IF( COLMAX .NE. ZERO ) THEN
                  CALL DSKR2( UPLO, K-2, ONE/A( K-1,K ), A( 1, K ), 1,
     $                 A( 1, K-1 ), 1, A( 1, 1 ), LDA )

*     Store L(k+1) in A(k)
                  CALL DSCAL(K-2, ONE/A( K-1, K ), A(1, K), 1)
               END IF
*     Store Pivot
               IPIV( K-1 ) = KP
            ELSE
               IPIV( K-1 ) = K-1
            END IF
10      CONTINUE

      ELSE
*     Factorize A as L * T * L^T using the lower triangle of A

         IPIV( 1 ) = 1

         DO 20 K=1, N-1
*     Either all columns or only the odd ones (MODE = 'P')
            IF( MOD(K, STEP).EQ.1 .OR. STEP.EQ.1 ) THEN
*     Find the pivot
               KP = K + IDAMAX(N-K, A( K+1, K ), 1)
               COLMAX = ABS( A( KP, K ) )

               IF( COLMAX.EQ.ZERO ) THEN
*     The column is completely zero - do nothing
                  IF( INFO.EQ.0 ) THEN
                     INFO = K
                  END IF
                  KP = K+1
               END IF

*     swap rows and columns K+1 and IMAX in the
*     full matrix A(1:N,1:N)
               KK = K+1

               IF( KP .NE. KK ) THEN
                  IF( KP.LT.N ) THEN
                     CALL DSWAP( N-KP, A( KP+1, KK ), 1,
     $                    A( KP+1, KP ),1 )
                  END IF

                  CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1,
     $                 A( KP, KK+1 ), LDA )

                  CALL DSWAP( K, A( KK, 1), LDA, A( KP, 1), LDA)

                  CALL DSCAL(KP-KK, -ONE, A(KK+1, KK), 1)
                  CALL DSCAL(KP-KK-1, -ONE, A(KP, KK+1), LDA)
               END IF

*     Update the trailing submatrix A(K+2:N, K+2:N) in a rank 2 update
*     (The column/row K+1 is not affected by the update)
               IF( COLMAX .NE. ZERO .AND. K+2 .LE. N) THEN
                  CALL DSKR2( UPLO, N-K-1, ONE/A( K+1,K ),
     $                 A( K+2, K ), 1, A( K+2, K+1 ), 1,
     $                 A( K+2, K+2 ), LDA )

*     Store L(k+1) in A(k)
                  CALL DSCAL(N-K-1, ONE/A( K+1, K ), A(K+2, K), 1)
               END IF

*     Store Pivot
               IPIV( K+1 ) = KP
            ELSE
               IPIV( K+1 ) = K+1
            END IF
 20      CONTINUE

      END IF

      END
