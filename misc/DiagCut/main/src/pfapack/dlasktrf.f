      SUBROUTINE DLASKTRF( UPLO, MODE, N, NB, A, LDA, IPIV, W, LDW,
     $                     INFO)
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            INFO, LDA, LDW, N, NB, STEP
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
      DOUBLE PRECISION   W( LDW, * )
*
*  Purpose
*  =======
*
*  DLASKTRF computes NB steps of the factorization of a skew-symmetric
*  matrix A using the Parlett-Reid algorithm:
*
*     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
*
*  If UPLO = 'U', DLASKTRF reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', DLASKTRF reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  Alternatively, the routine can also be used to compute a partial
*  tridiagonal form
*
*  This is an auxiliary routine called by DSKTRF.
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
*  NB      (input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the skew-symmetric matrix A.
*            If UPLO = 'U', the leading n-by-n upper triangular part
*               of A contains the upper triangular part of the matrix A,
*               and the strictly lower triangular part of A is not referenced.
*            If UPLO = 'L', the leading n-by-n lower triangular part
*               of A contains the lower triangular part of the matrix A,
*               and the strictly upper triangular part of A is not referenced.
*          On exit:
*          if UPLO = 'U', the last NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements above the diagonal
*            represent the upper unit triangular matrix U. If MODE = 'P',
*            only the even columns of the last NB columns contain
*            meaningful values.
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements below the diagonal
*            represent the lower unit triangular matrix L. If MODE = 'P',
*            only the  odd columns of the first NB columns contain
*            meaningful values.
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
*  W       (workspace) DOUBLE PRECISION array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A. (The update is performed in this routine)
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     .. Local Scalars ..
      INTEGER            K, KK, KP, NPANEL, WK
      DOUBLE PRECISION   COLMAX, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX, DISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, DSKR2K, DCOPY, DGEMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     No safety checks, it's an internal function

      INFO = 0

      IF( LSAME( MODE, 'P' ) ) THEN
         STEP = 2
      ELSE
         STEP = 1
      END IF

*     double the amount of panels before the block update, if STEP == 2
      NPANEL = NB * STEP

      IF( LSAME( UPLO, 'U' ) ) THEN

*     Factorize A as U * T * U^T using the upper triangle of A

*     Consider the last NB columns of the matrix:
*     First compute the transformations and apply them to the last NB
*     columns (only), then apply the accumulated transformations to the rest
*     of the matrix

         WK = 0
         DO 10 K=N, MAX(N-NPANEL+1, 2), -1
*
*     Update A(1:K,K) with all the accumulated transformations
*     (if K<N: K=N,N-1 is updated when computing the Gauss vector,
*     K=N-1 is not affected by the K=N transform)
*
            KK = K-1

            IF( K .LT. N) THEN

               IF( WK .GT. 0 ) THEN
                  A( K, K ) = ZERO
                  CALL DGEMV( 'N', K, WK, +ONE, A( 1, N-(WK-1)*STEP ),
     $                 LDA*STEP, W( K, NB-WK+1 ), LDW, ONE,
     $                 A( 1, K ), 1 )
                  CALL DGEMV( 'N', K, WK, -ONE, W( 1, NB-WK+1 ),
     $                 LDW, A( K, N-(WK-1)*STEP ), LDA*STEP, ONE,
     $                 A( 1, K ), 1 )
                  A( K, K ) = ZERO
               END IF

*     Store the (updated) column K in W(:,NB-WK+1)
               IF( MOD(K, STEP).EQ.1 .OR. STEP.EQ.1 ) THEN
                  WK = WK + 1
                  CALL DCOPY(K, A(1, K), 1, W(1, NB-WK+1 ), 1)
               END IF
            END IF

            IF( MOD(K, STEP) .EQ. 0) THEN
*     For STEP == 1, process every column, but if
*     STEP == 2, do only things for the even columns

*     Find the pivot
               KP = IDAMAX(K-1, A( 1, K ), 1)
               COLMAX = ABS( A( KP, K ) )

               IF( COLMAX.EQ.ZERO ) THEN
*     The column is completely zero - do nothing
                  IF( INFO.EQ.0 ) THEN
                     INFO = K
                  END IF
                  KP = KK
               END IF

*     swap rows and columns K+1 and KP in the
*     full matrix A(1:N,1:N)
*     Also, swap the first K-1 columns of W

               IF( KP .NE. KK ) THEN
                  CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ),1)
                  CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1,
     $                 A( KP, KP+1 ), LDA )

                  CALL DSWAP( N-K+1, A( KK, K), LDA, A( KP, K), LDA)

                  CALL DSCAL(KK-KP, -ONE, A(KP, KK), 1)
                  CALL DSCAL(KK-KP-1, -ONE, A(KP, KP+1), LDA)

*     Swap in the last columns of W
                  IF( WK .GT. 0 ) THEN
                     CALL DSWAP( WK, W( KK, NB-WK+1 ), LDW,
     $                    W( KP, NB-WK+1 ), LDW)
                  END IF
               END IF

*     (The column/row K+1 is not affected by the update)
               IF( COLMAX .NE. ZERO ) THEN
*     Store L(k+1) in A(k)
                  CALL DSCAL(K-2, ONE/A( K-1, K ), A(1, K), 1)
               END IF

*     Delay the update of the trailing submatrix (done at the beginning
*     of each loop for every column of the first NB columns, and then
*     finally in a rank-3 update for the last N-NB block)

*     Store Pivot
               IPIV( K-1 ) = KP

            ELSE
*     STEP == 2 and an even column, do nothing
               IPIV( K-1 ) = K-1
            END IF
 10      CONTINUE

*     Now update the leading A(1:N-NB, 1:N-NB) submatrix in a level 3 update,
*     if necessary
         IF( N-NPANEL+1 .GT. 2 ) THEN

*     For this we have to set the N-NB,N-NB+1 entry to zero,
*     but restore it later
            T = A( N-NPANEL, N-NPANEL+1 )
            A( N-NPANEL, N-NPANEL+1 ) = ZERO

            IF( WK .LT. NB) THEN
*     Store the column N-NB in W, and update it in W
               CALL DCOPY(N-NPANEL, A(1, N-NPANEL), 1,
     $              W(1, 1), 1)

               W( N-NPANEL, 1 ) = ZERO
               CALL DGEMV( 'N', N-NPANEL, WK, +ONE,
     $              A( 1, N-(WK-1)*STEP ), LDA*STEP,
     $              W( N-NPANEL, NB-WK+1 ), LDW, ONE, W( 1, 1 ), 1 )
               CALL DGEMV( 'N', N-NPANEL, WK, -ONE, W( 1, NB-WK+1 ),
     $              LDW, A( N-NPANEL, N-(WK-1)*STEP ), LDA*STEP, ONE,
     $              W( 1, 1 ), 1 )
               W( N-NPANEL, 1 ) = ZERO

               WK = WK + 1
            END IF

*     Now do the rank-3 update
            CALL DSKR2K( UPLO, "N", N-NPANEL, NB, ONE,
     $           A(1, N-(WK-1)*STEP), LDA*STEP,
     $           W(1,1), LDW, ONE, A(1, 1), LDA)

            A( N-NPANEL, N-NPANEL+1 ) = T
         END IF

      ELSE

*     Factorize A as L * T * L^T using the lower triangle of A

*     Consider the first NB columns of the matrix:
*     First compute the transformations and apply them to the first NB
*     columns (only), then apply the accumulated transformations to the rest
*     of the matrix

         WK = 0
         DO 30 K=1, MIN(NPANEL, N-1)

*
*     Update A(K+1:n,K+1) with all the accumulated transformations
*     (if K>2: K=1,2 is updated when computing the Gauss vector,
*     K=2 is not affected by the K=1 transform)
*
            KK = K+1

            IF( K .GT. 1) THEN

               IF( WK .GT. 0) THEN
                  A( K, K ) = ZERO
                  CALL DGEMV( 'N', N-K+1, WK, +ONE, A( K, 1 ),
     $                 LDA*STEP, W( K, 1 ), LDW, ONE, A( K, K ), 1 )
                  CALL DGEMV( 'N', N-K+1, WK, -ONE, W( K, 1 ),
     $                 LDW, A( K, 1 ), LDA*STEP, ONE, A( K, K ), 1 )
                  A( K, K ) = ZERO
               END IF

*     Store the (updated) column K in W(:,WK)
               IF( MOD(K, STEP) . EQ. 0 ) THEN
                  WK = WK + 1
                  CALL DCOPY(N-K+1, A(K, K), 1, W(K, WK), 1)
               END IF
            END IF

            IF( MOD(K, STEP) .EQ. 1 .OR. STEP .EQ. 1) THEN
*     For STEP == 1, process every column, but if
*     STEP == 2, do only things for the odd columns

*     Find the pivot
               KP = K + IDAMAX(N-K, A( K+1, K ), 1)
               COLMAX = ABS( A( KP, K ) )

               IF( COLMAX.EQ.ZERO ) THEN
*     The column is completely zero - do nothing
                  IF( INFO.EQ.0 ) THEN
                     INFO = K
                  END IF
                  KP = KK
               END IF

*     swap rows and columns K+1 and KP in the
*     full matrix A(1:N,1:N)
*     Also, swap the first K-1 columns of W

               IF( KP .NE. KK ) THEN
                  IF( KP.LT.N ) THEN
                     CALL DSWAP( N-KP, A( KP+1, KK ), 1,
     $                           A( KP+1, KP ),1 )
                  END IF

                  CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1,
     $                        A( KP, KK+1 ), LDA )

                  CALL DSWAP( K, A( KK, 1), LDA, A( KP, 1), LDA)

*     The matrix is anti-symmetric, hence swaps beyond the diagonal need a -1

                  CALL DSCAL( KP-KK, -ONE, A(KK+1, KK), 1 )
                  CALL DSCAL( KP-KK-1, -ONE, A(KP, KK+1), LDA )

*     Swap in the first columns of W
                  CALL DSWAP( WK, W( KK, 1 ), LDW, W( KP, 1 ), LDW)
               END IF

*     (The column/row K+1 is not affected by the update)
               IF( COLMAX .NE. ZERO .AND. K .LE. N-2) THEN
*     Store L(k+1) in A(k)
                  CALL DSCAL( N-K-1, ONE/A( K+1, K ), A(K+2, K), 1 )
               END IF

*     Delay the update of the trailing submatrix (done at the beginning
*     of each loop for every column of the first NB columns, and then
*     finally in a rank-3 update for the last N-NB block)

*     Store Pivot
               IPIV( K+1 ) = KP

            ELSE
*     STEP == 2 and an even column, do nothing
               IPIV(K+1) = K+1
            END IF
 30      CONTINUE

*     Now update the trailing A(NB+1:N, NB+1:N) submatrix in a level 3 update,
*     if necessary

         IF( NPANEL .LT. N-1) THEN
*     For this we have to set the NB+1,NB entry to zero, but restore it later
            T = A( NPANEL+1, NPANEL )
            A( NPANEL+1, NPANEL ) = ZERO

            IF( WK .LT. NB) THEN
*     Store the column NB+1 in W, and update it in W
               CALL DCOPY(N-NPANEL, A(NPANEL+1, NPANEL+1), 1,
     $              W(NPANEL+1, NB), 1)

               W( NPANEL+1, NB ) = ZERO
               CALL DGEMV( 'N', N-NPANEL, NB-1, +ONE, A( NPANEL+1, 1 ),
     $              LDA*STEP, W( NPANEL+1, 1 ), LDW,
     $              ONE, W( NPANEL+1, NB ), 1 )
               CALL DGEMV( 'N', N-NPANEL, NB-1, -ONE, W( NPANEL+1, 1 ),
     $              LDW, A( NPANEL+1, 1 ), LDA*STEP, ONE,
     $              W( NPANEL+1, NB ), 1 )
               W( NPANEL+1, NB ) = ZERO
            END IF

*     Now do the rank-3 update
            CALL DSKR2K( UPLO, "N", N-NPANEL, NB, ONE,
     $           A(NPANEL+1,1), LDA*STEP,
     $           W(NPANEL+1,1), LDW, ONE, A(NPANEL+1, NPANEL+1), LDA)

            A(NPANEL+1, NPANEL)=T
         END IF

      END IF

      END
