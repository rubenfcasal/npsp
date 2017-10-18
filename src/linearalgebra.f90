!   --------------------------------------------------------------------
    subroutine DPOSV_R(N, NRHS, A, B, INFO)
!   --------------------------------------------------------------------
!   Interface to the LAPACK routine DPOSV:
!   DPOSV computes the solution to a real system of linear equations
!     A * X = B,
!   where A is an N-by-N symmetric positive definite matrix and X and B
!   are N-by-NRHS matrices.
!   The Cholesky decomposition is used to factor A as
!     A = L * L**T,  if UPLO = 'L',
!   where L is a lower triangular matrix. The factored form of A is then
!   used to solve the system of equations A * X = B.
!   Parameters
!   [in]	N INTEGER
!         The number of linear equations, i.e., the order of the
!         matrix A.  N >= 0.
!   [in]	NRHS INTEGER
!         The number of right hand sides, i.e., the number of columns
!         of the matrix B.  NRHS >= 0.
!   [in,out]	A DOUBLE PRECISION array, dimension (LDA,N)
!         On entry, the symmetric matrix A. The
!         leading N-by-N lower triangular part of A contains the lower
!         triangular part of the matrix A, and the strictly upper
!         triangular part of A is not referenced.
!
!         On exit, if INFO = 0, the factor L from the Cholesky
!         factorization A = L*L**T.
!   [in,out]	B DOUBLE PRECISION array, dimension (LDB,NRHS)
!         On entry, the N-by-NRHS right hand side matrix B.
!         On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!   [out]	INFO INTEGER
!         = 0:  successful exit
!         < 0:  if INFO = -i, the i-th argument had an illegal value
!         > 0:  if INFO = i, the leading minor of order i of A is not
!               positive definite, so the factorization could not be
!               completed, and the solution has not been computed.
!   --------------------------------------------------------------------
    implicit none
    integer N, NRHS
!   input/output:
    integer INFO
    real(8) A(N, N), B(N, NRHS)
    call DPOSV('L', N, NRHS, A, N, B, N, INFO)
    return
    end subroutine DPOSV_R


!   --------------------------------------------------------------------
    subroutine DNRM2_R(nd, x, nx, z, nz, dist)
!   Interface to the BLAS routine DNRM2
!   --------------------------------------------------------------------
    implicit none
    integer nd, nx, nz, i, j
    real(8)  x(nd, nx), z(nd, nz)
    real(8), external :: DNRM2   ! DNRM2(N,X,INCX)
!   output:
    real(8)  dist(nz, nx)
!       Recorrer datos
        do j = 1, nx
            do i = 1, nz
                dist(i, j) = DNRM2(nd, x(1:nd, j) - z(1:nd, i), 1)
            end do
        end do  ! Recorrer datos
    return
    end subroutine DNRM2_R


