!     ------------------------------------------------------------------
!     [DSYTRFI]   Calcula el determinante y la inversa de una matriz
!                 simétrica a partir de la factorización UDU'.
!     ------------------------------------------------------------------
      SUBROUTINE DSYTRFI(N, A, AI, DET)
      IMPLICIT NONE
      INTEGER N
      REAL(8)  A(N,N), AI(N,N), DET
!     Variables locales
      INTEGER IPIV(N), LWORK, INFO, i
      REAL(8) tmp(1)
      REAL(8), ALLOCATABLE ::  WORK(:)
!     ------------------------------------------------------------------
      AI = A
!     Factorizar matriz
!     DSYTRF computes the factorization of a real symmetric matrix A using
!           the Bunch-Kaufman diagonal pivoting method.  The form of the
!           factorization is A = U*D*U**T  or  A = L*D*L**T
!     The leading N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.
      LWORK = -1    ! Determine the block size
      CALL DSYTRF( 'U', N, AI, N, IPIV, tmp, LWORK, INFO )
      LWORK = NINT(tmp(1))
      ALLOCATE(WORK(LWORK), STAT=i)
      IF (i.NE.0) CALL Error(i, 'DSYTRFI: ALLOCATE')
      CALL DSYTRF( 'U', N, AI, N, IPIV, WORK, LWORK, INFO )
      IF (INFO.NE.0) CALL Error(INFO, 'DSYTRFI: DSYTRF')
!     Calcular el determinante
!     If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!           interchanged and D(k,k) is a 1-by-1 diagonal block.
!     If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!           columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!           is a 2-by-2 diagonal block.
      DET = 1.0d0
      DO i = 1, N
         IF (IPIV(i).GT.0) THEN
            DET = DET * AI(i,i)
         ELSE IF((i.GT.1).AND.(IPIV(i).LT.0).AND.(IPIV(i).EQ.IPIV(i-1))) THEN
            DET = DET * ( AI(i,i) * AI(i-1,i-1) - AI(i-1,i) * AI(i,i-1) )
         END IF
      END DO
!     Invertir la matriz
!     DSYTRI computes the inverse of a real symmetric indefinite matrix
!           A using the factorization computed by DSYTRF.
      CALL DSYTRI( 'U', N, AI, N, IPIV, WORK, INFO )
      IF (INFO.NE.0) CALL Error(INFO, 'DSYTRFI: DSYTRI')
!     Liberar memoria
      DEALLOCATE (WORK,STAT=i)
      IF (i.NE.0) CALL Error(i, 'DSYTRFI: DEALLOCATE')
      RETURN
      END SUBROUTINE DSYTRFI



