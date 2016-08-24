      SUBROUTINE DLASKTRD( UPLO, MODE, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- Written on 10/22/2010
*     Michael Wimmer, Universiteit Leiden
*     Derived from the LAPACK routine ZLATRD (www.netlib.org)
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MODE
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASKTRD reduces NB rows and columns of a real skew-symmetric matrix A to
*  skew-symmetric tridiagonal form by an orthogonal similarity
*  transformation Q^T * A * Q, and returns the matrices V and W which are
*  needed to apply the transformation to the unreduced part of A.
*
*  If UPLO = 'U', DLASKTRD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', DLASKTRD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  Alternatively, the routine can also be used to compute a partial
*  tridiagonal form
*
*  This is an auxiliary routine called by DSKTRD.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          skew-symmetric matrix A is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  MODE    (input) CHARACTER*1
*          = 'N':  A is fully tridiagonalized
*          = 'P':  A is partially tridiagonalized for Pfaffian computation
*
*  N       (input) INTEGER
*          The order of the matrix A. N must be even if MODE = 'P'.
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
*            with the array TAU, represent the unitary matrix Q as a
*            product of elementary reflectors. If MODE = 'P' only the
*            even columns of the last NB columns contain meaningful values.
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements below the diagonal
*            with the array TAU, represent the  unitary matrix Q as a
*            product of elementary reflectors. If MODE = 'P' only the
*            odd columns of the first NB columns contain meaningful values.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*          elements of the last NB columns of the reduced matrix;
*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
*          the first NB columns of the reduced matrix.
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors, stored in
*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*          See Further Details.
*
*  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0,
     $                   ONE = 1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER            I, NW, NW2, STEP, NPANEL
      DOUBLE PRECISION   ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DSKMV, DLARFG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*

      IF( LSAME( MODE, 'P' ) ) THEN
         STEP = 2
      ELSE
         STEP = 1
      END IF
      NPANEL = NB * STEP

      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NPANEL columns of upper triangle
*
         NW=0
*         IW = NB + 1

         DO 10 I = N, MAX(N - NPANEL + 1, 2), -1

            NW2 = NW - MOD(I,STEP)
            IF( NW2 .GT. 0 ) THEN
*
*              Update A(1:i,i)
*
               A( I, I ) = ZERO
               CALL DGEMV( 'No transpose', I, NW2, +ONE,
     $                     A( 1, N-(NW2-1)*STEP ), LDA*STEP,
     $                     W( I, NB-NW2+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'No transpose', I, NW2, -ONE,
     $                     W( 1, NB-NW2+1 ), LDW,
     $                     A( I, N-(NW2-1)*STEP ), LDA*STEP,
     $                     ONE, A( 1, I ), 1 )
               A( I, I ) = ZERO
            END IF

*     In the Pfaffian mode, only zero all even columns
            IF( STEP.EQ.2 .AND. MOD( I, STEP ).EQ.1 ) THEN
               GOTO 10
            END IF

            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               ALPHA = A( I-1, I )
               CALL DLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = ALPHA
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL DSKMV( 'Upper', I-1, TAU(I-1), A, LDA,
     $                     A( 1, I ), 1, ZERO, W( 1, NB-NW ), 1 )

               IF( NW .GT. 0 ) THEN
                  CALL DGEMV( 'Transpose', I-1, NW, ONE,
     $                        W( 1, NB-NW+1 ), LDW, A( 1, I ), 1, ZERO,
     $                        W( I+1, NB-NW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, NW, TAU(I-1),
     $                        A( 1, N-(NW-1)*STEP ), LDA*STEP,
     $                        W( I+1, NB-NW ), 1,
     $                        ONE, W( 1, NB-NW ), 1)
                  CALL DGEMV( 'Transpose', I-1, NW, ONE,
     $                        A( 1, N-(NW-1)*STEP ), LDA*STEP,
     $                        A( 1, I ), 1, ZERO, W( I+1, NB-NW ), 1 )
                  CALL DGEMV( 'No transpose', I-1, NW, -TAU(I-1),
     $                        W( 1, NB-NW+1 ), LDW,
     $                        W( I+1, NB-NW ), 1,
     $                        ONE, W( 1, NB-NW ), 1 )
               END IF

*     One more complete entry in W
               NW = NW + 1

            END IF

*     Note: setting A(I-1,I) back to alpha happens in the calling routine

   10    CONTINUE
      ELSE
*
*        Reduce first NPANEL columns of lower triangle
*
         NW = 0

         DO 20 I = 1, MIN(NPANEL, N-1)
*
*           Update A(i:n,i)
*
            NW2 = NW - MOD(I+1,STEP)
            IF( NW2 .GT. 0 ) THEN
               A( I, I ) = ZERO
               CALL DGEMV( 'No transpose', N-I+1, NW2, +ONE, A( I, 1 ),
     $                     LDA*STEP, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
               CALL DGEMV( 'No transpose', N-I+1, NW2, -ONE, W( I, 1 ),
     $                     LDW, A( I, 1 ), LDA*STEP, ONE, A( I, I ), 1 )
               A( I, I ) = ZERO
            END IF

*     In the Pfaffian mode, only zero all odd columns
            IF( STEP.EQ.2 .AND. MOD( I, STEP ).EQ.0 ) THEN
               GOTO 20
            END IF

            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               ALPHA = A( I+1, I )
               CALL DLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*              This is given by tau A^(i)*v^*=tau(A*v^* + VW^T v^* - WV^T v^*)

               CALL DSKMV( 'Lower', N-I, TAU( I ),
     $                     A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, NW+1 ), 1 )

               IF( NW .GT. 0 ) THEN
                  CALL DGEMV( 'Transpose', N-I, NW, ONE,
     $                 W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO,
     $                 W( 1, NW+1 ), 1 )
                  CALL DGEMV( 'No transpose', N-I, NW, TAU( I ),
     $                 A( I+1, 1 ), LDA*STEP, W( 1, NW+1 ), 1,
     $                 ONE, W( I+1, NW+1 ), 1 )
                  CALL DGEMV( 'Transpose', N-I, NW, ONE,
     $                 A( I+1, 1 ), LDA*STEP, A( I+1, I ), 1, ZERO,
     $                 W( 1, NW+1 ), 1 )
                  CALL DGEMV( 'No transpose', N-I, NW, -TAU( I ),
     $                 W( I+1, 1 ), LDW, W( 1, NW+1 ), 1,
     $                 ONE, W( I+1, NW+1 ), 1 )
               END IF

*     One more complete entry in W
               NW = NW + 1
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASKTRD
*
      END
