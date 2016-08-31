      SUBROUTINE DSKPFA( UPLO, MTHD, N, A, LDA, PFAFF,
     $                   IWORK, WORK, LWORK, INFO)
*
* -- Written on 10/22/2010
*    Michael Wimmer, Universiteit Leiden
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, MTHD
      INTEGER            INFO, LDA, LWORK, N
      DOUBLE PRECISION   PFAFF
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSKPFA computes the Pfaffian of a real skew-symmetric matrix.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  MTHD    (input) CHARACTER*1
*          = 'P': Compute Pfaffian using Parlett-Reid algorithm (recommended)
*          = 'H': Compute Pfaffian using Householder reflections
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the skew-symmetric matrix A.
*             If UPLO = 'U', the upper triangular part of A contains
*                the upper triangular part of the matrix A, and the
*                strictly lower triangular part of A is not referenced.
*             If UPLO = 'L', the lower triangular part of A contains
*                the lower triangular part of the matrix A, and the
*                strictly upper triangular part of A is not referenced.
*          If the matrix size is odd, A is not referenced. If the matrix
*          size is even, A is overwritten by values generated during
*          the computation.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  PFAFF   (output) DOUBLE PRECISION
*          The value of the Pfaffian.
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*          Not referenced if MTHD = 'H'.
*
*  WORK    (workspace) DOUBLE PRECISION array,
*             dimension (MAX(1, LWORK)), if MTHD = 'P';
*             dimension (MAX(2*N-1,LWORK)), if MTHD = 'H'.
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If MTHD = 'P', LWORK >= 1,
*          If MTHD = 'H', LWORK >= 2*N-1.
*
*          For optimum performance LWORK >= N*NB for MTHD = 'P' or
*          LWORK >= N*NB+2*N-2 for MTHD = 'H', where NB is the
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
*  The Pfaffian is computed by bringing the skew-symmetric matrix A into
*  a partial tridiagonal form pT, either by computing a partial L pT L^T
*  decomposition (MTHD = 'P'), or by a by a unitary congruence transformation
*  Q^H * A * Q^* = pT (MTHD = 'H').
*  These transformations are computed by the routines DSKTRF or DSKTRD,
*  respectively (for further details see there).
*
*
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0 )
      PARAMETER          ( ZERO = 0.0D+0 )

      INTEGER            I

*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER, LTL

*     .. External Subroutines ..
      EXTERNAL           XERBLA, DSKTRD, DSKTRF
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LTL = LSAME( MTHD, 'P' )
      LQUERY = ( LWORK.EQ.-1 )

      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LTL .AND. .NOT.LSAME( MTHD, 'H' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      ELSE IF( MOD(N,2).NE.1 .AND. .NOT.LTL .AND.
     $         LWORK.LT.2*N-1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF

      IF( INFO.EQ.0 .AND. LQUERY) THEN
         IF( MOD(N,2).EQ.1 ) THEN
            WORK(1) = 1
         ELSE IF( LTL ) THEN
*     Defer workspace query to DSKTRF
            CALL DSKTRF( UPLO, "P", N, A, LDA, IWORK, WORK, LWORK,
     $                   INFO )
         ELSE
*     Defer workspace query to DSKTRD
            CALL DSKTRD( UPLO, "P", N, A, LDA, WORK, WORK, WORK,
     $                   LWORK, INFO)
            WORK(1) = WORK(1) + 2*N - 2
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSKPFA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

*     Quick return if possible
      IF( N.EQ.0 ) THEN
         PFAFF = ONE
         RETURN
      ELSE IF( MOD(N,2).EQ.1 ) THEN
         PFAFF = ZERO
         RETURN
      END IF

      IF( LTL ) THEN
*     Compute tridiagonal form
         CALL DSKTRF( UPLO, "P", N, A, LDA, IWORK, WORK, LWORK, INFO )

         IF( UPPER ) THEN
            PFAFF = ONE

            DO 10 I = 1, N-1, 2
               PFAFF = PFAFF * A( I, I+1 )

*     Accumulate the determinant of the permutations
               IF( IWORK( I ) .NE. I ) PFAFF = -PFAFF
 10         CONTINUE

         ELSE
            PFAFF = ONE

            DO 20 I = 1, N-1, 2
               PFAFF = PFAFF * (-A( I+1, I ))

*     Accumulate the determinant of the permutations
               IF( IWORK( I+1 ) .NE. I+1 ) PFAFF = -PFAFF
 20         CONTINUE

         END IF
      ELSE

*     Reduce to tridiagonal form
         CALL DSKTRD(UPLO, "P", N, A, LDA, WORK(1), WORK(N),
     $               WORK( 2*N-1 ), LWORK-2*N+2, INFO)

         PFAFF = ONE

         IF( UPPER ) THEN
*     Multiply every other entry on the superdiagonal
            DO 30 I = 1, N-1, 2
               PFAFF = PFAFF * WORK( I )

*     Accumulate the determinant of the Householder reflection
*     (which in the real case can only be +1 or -1)
               IF (WORK( N-1+I ) .GT. ZERO) PFAFF = -PFAFF
 30         CONTINUE

         ELSE

*     Multiply every other entry on the superdiagonal
            DO 40 I = 1, N-1, 2
               PFAFF = PFAFF * (-WORK( I ))

*     Accumulate the determinant of the Householder reflection
*     (which in the real case can only be +1 or -1)
               IF (WORK( N-1+I ) .GT. ZERO) PFAFF = -PFAFF
 40         CONTINUE

         END IF

*     Shift optimal workspace size to first position in the WORK array
         WORK( 1 ) = WORK( 2*N-1 ) + 2*N-2
      END IF

      RETURN

*     end of DSKPFA

      END
