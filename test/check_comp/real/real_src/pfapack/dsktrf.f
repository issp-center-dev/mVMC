      SUBROUTINE DSKTRF( UPLO, MODE, N, A, LDA, IPIV, WORK, LWORK, INFO)
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSKTRF computes the factorization of a skew-symmetric matrix A
*  using the Parlett-Reid algorithm:
*
*     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
*
*  where U (or L) unit upper (lower) triangular matrix (^T denotes
*  the transpose), T is a skew-symmetric tridiagonal matrix and P
*  is a permutation matrix. In addition to being unit triangular,
*  U(1:n-1,n)=0 and L(2:n,1)=0.
*  Instead of a full tridiagonalization, DSKTRF can also compute a
*  partial tridiagonal form for computing the Pfaffian.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
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
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of WORK.  LWORK >=1.  For best performance
*          LWORK >= N*NB, where NB is the block size returned by ILAENV
*          (at the moment, uses the same block size as DSYTRF from Lapack).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
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
*  The normal use for DSKTRF is to compute the U T U^T or L T L^T
*  decomposition of a skew-symmetric matrix with pivoting. This mode
*  is chosen by setting MODE = 'N' ("normal" mode). The other
*  use of DSKTRF is the computation the Pfaffian of a skew-symmetric matrix,
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
*  or in the even columns (if UPLU = 'U') are properly computed by DSKTRF.
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
*  element not computed by DSKTRF.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, NORMAL
      INTEGER            IINFO, J, K, K2, PIV, LWKOPT, NB, NBMIN, NPANEL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASKTRF, DSKTF2, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
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
      ELSE IF( .NOT.NORMAL .AND. MOD(N,2).EQ.1 ) THEN
*     If MODE = 'P', we need an even-dimensional matrix
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size
*        Note: Obviously, since Lapack does not know about DSKTRF, I query the
*        similar routine 'DSYTRF'
         NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKTRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      NBMIN=NB
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IF( LWORK .LT. N*NB ) THEN
            NB = MAX( LWORK / N, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         NB = N
      END IF

      IF( NB.LT.NBMIN ) NB = N
*     Note: NB = N means do not use blocked code

*     Quick return if possible
      IF( N .EQ. 0 ) RETURN

      IF( LSAME( MODE, 'N' ) ) THEN
         NPANEL = NB
      ELSE
         NPANEL = MIN(NB*2, N)
      END IF

      IF( UPPER ) THEN
*
*     Factorize A as L*T*L^T using the lower triangle of A

         IPIV( N ) = N
*

*     Loop throgh the system in steps of NPANEL
         DO 10 K = N, MAX(NPANEL,1), -NPANEL
*
            IF( K.GE.NPANEL*2 ) THEN
*
*     Factorize columns k-nb*step+1:k of A and use blocked code to
*     update columns 1:k-nb*step
*
               CALL DLASKTRF( UPLO, MODE, K, NB, A, LDA,
     $                        IPIV, WORK, N, IINFO )

               K2 = K-NPANEL
            ELSE
*
*     Use unblocked code to factorize columns 1:k of A
*
*     IPIV( K ) is overwritten by DSKTF2, need to restore it later
               PIV = IPIV( K )

               CALL DSKTF2( UPLO, MODE, K, A, LDA, IPIV, IINFO )

               IPIV( K ) = PIV

               K2 = 1
            END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $           INFO = IINFO

*     Perform the missing row interchanges in the trailing N-K columns
            IF( K .LT. N ) THEN
               DO 20 J=K-1, K2, -1
                  CALL DSWAP( N-K, A( J, K+1 ), LDA,
     $                 A( IPIV( J ), K+1 ), LDA )
 20            CONTINUE
            END IF

 10   CONTINUE
*
      ELSE
*
*        Factorize A as L*T*L^T using the lower triangle of A

         IPIV( 1 ) = 1
*
*        Loop throgh the system in steps of NPANEL
         DO 30 K = 1, MIN(N-NPANEL+1, N-1), NPANEL
*
            IF( K.LE.N-NPANEL*2+1 ) THEN
*
*     Factorize columns k:k+nb-1 of A and use blocked code to
*     update columns k+nb-1:n
*
               CALL DLASKTRF( UPLO, MODE, N-K+1, NB, A( K, K ), LDA,
     $                        IPIV( K ), WORK, N, IINFO )

               K2 = K + NPANEL
            ELSE
*
*     Use unblocked code to factorize columns k:n of A
*
*     IPIV( K ) is overwritten by DSKTF2, need to restore it later
               PIV = IPIV( K )

               CALL DSKTF2( UPLO, MODE, N-K+1, A( K, K ), LDA,
     $                      IPIV( K ), IINFO )

               IPIV( K ) = PIV

               K2 = N
            END IF
*
*     Set INFO on the first occurrence of a zero pivot
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $           INFO = IINFO + K - 1
*
*     Adjust IPIV
*
            DO 40 J = K+1, K2
               IPIV( J ) = IPIV( J ) + K - 1
 40         CONTINUE

*     Perform the missing row interchanges in the leading K-1 columns
            IF( K .GT. 1 ) THEN
               DO 50 J=K+1, K2
                  CALL DSWAP( K-1, A( J, 1 ), LDA,
     $                 A( IPIV( J ), 1 ), LDA )
 50            CONTINUE
            END IF


 30   CONTINUE
*
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DSKTRF
*
      END
