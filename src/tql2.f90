!     ------------------------------------------------------------------
!     [tql2]  Modificación de la subrutina del paquete EISPACK 
!             (descargada de http://www.netlib.org/eispack/)
!             para el cálculo de autovalores (y autovectores) de una
!             matriz tridiagonal simétrica.
!
!     NOTA: Se eliminó la ordenación de autovalores y autovectores.
!     ------------------------------------------------------------------
      SUBROUTINE TQL2(Nm,N,D,E,Z,Ierr)
      IMPLICIT NONE
      INTEGER i, j, k, l, m, N, ii, l1, l2, Nm, mml, Ierr
      DOUBLE PRECISION D(N), E(N), Z(Nm,N)
      DOUBLE PRECISION c, c2, c3, dl1, el1, f, g, h, p, r, s, &
                     & s2, tst1, tst2, PYTHAG
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306 (1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240 (1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      Ierr = 0
      IF ( N==1 ) GOTO 400
!
      DO i = 2, N
         E(i-1) = E(i)
      ENDDO
!
      f = 0.0D0
      tst1 = 0.0D0
      E(N) = 0.0D0
!
      DO l = 1, N
         j = 0
         h = DABS(D(l)) + DABS(E(l))
         IF ( tst1<h ) tst1 = h
!     .......... look for small sub-diagonal element ..........
         DO m = l, N
            tst2 = tst1 + DABS(E(m))
            IF ( tst2==tst1 ) GOTO 50
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
         ENDDO
!
 50      IF ( m==l ) GOTO 200
 100     IF ( j==30 ) GOTO 300
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = D(l)
         p = (D(l1)-g)/(2.0D0*E(l))
         r = PYTHAG(p,1.0D0)
         D(l) = E(l)/(p+DSIGN(r,p))
         D(l1) = E(l)*(p+DSIGN(r,p))
         dl1 = D(l1)
         h = g - D(l)
         IF ( l2>N ) GOTO 150
!
         DO i = l2, N
            D(i) = D(i) - h
         ENDDO
!
 150     f = f + h
!     .......... ql transformation ..........
         p = D(m)
         c = 1.0D0
         c2 = c
         c3 = c
         el1 = E(l1)
         s = 0.0D0
         s2 = 0.0D0
         mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
         DO ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c*E(i)
            h = c*p
            r = PYTHAG(p,E(i))
            E(i+1) = s*r
            s = E(i)/r
            c = p/r
            p = c*D(i) - s*g
            D(i+1) = h + s*(c*g+s*D(i))
!     .......... form vector ..........
            DO k = 1, N
               h = Z(k,i+1)
               Z(k,i+1) = s*Z(k,i) + c*h
               Z(k,i) = c*Z(k,i) - s*h
            ENDDO
!
         ENDDO
!
         p = -s*s2*c3*el1*E(l)/dl1
         E(l) = s*p
         D(l) = c*p
         tst2 = tst1 + DABS(E(l))
         IF ( tst2>tst1 ) GOTO 100
 200     D(l) = D(l) + f
      ENDDO
!     .......... order eigenvalues and eigenvectors ..........
!      do 300 ii = 2, n
!         i = ii - 1
!         k = i
!         p = d(i)
!
!         do 260 j = ii, n
!            if (d(j) .ge. p) go to 260
!            k = j
!            p = d(j)
!  260    continue
!
!         if (k .eq. i) go to 300
!         d(k) = d(i)
!         d(i) = p
!
!         do 280 j = 1, n
!            p = z(j,i)
!            z(j,i) = z(j,k)
!            z(j,k) = p
!  280    continue
!
!  300 continue
!
      GOTO 400
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 300  Ierr = l
 400  CONTINUE
      END

!     ------------------------------------------------------------------
!     [pythag]    Función utilizada por la subrutina tql2
!                 (paquete EISPACK).
!     ------------------------------------------------------------------
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!     ------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      IMPLICIT NONE
      DOUBLE PRECISION A, B
!
      DOUBLE PRECISION p, r, s, t, u
      p = DMAX1(DABS(A),DABS(B))
      IF ( p==0.0D0 ) GOTO 200
      r = (DMIN1(DABS(A),DABS(B))/p)**2
 100  CONTINUE
      t = 4.0D0 + r
      IF ( t==4.0D0 ) GOTO 200
      s = r/t
      u = 1.0D0 + 2.0D0*s
      p = u*p
      r = (s/u)**2*r
      GOTO 100
 200  PYTHAG = p
      CONTINUE
      END
 
