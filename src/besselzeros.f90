!-----------------------------------------------------------------------
!   [besselzeros.f90]   Utilidades para discretización de la distribución
!                       espectral de un variograma multidimensional
!
!   Interfaces con R:
!       disc_sbv      ptos de discretización para un modelo de variograma de 
!                     Shapiro-Botha (R "disc_sbv")
!
!   Autor: (c) Ruben Fernandez-Casal                Creacion: Abr 2002
!   Revisiones: Mar 2013
!-----------------------------------------------------------------------

!     ------------------------------------------------------------------
!     [disc_sbv]  Obtiene los ptos de discretización de la función de 
!                 distribución espectral para un modelo de Shapiro-Botha 
!                 extendido. Basado en el artículo:
!                 Gorsich y Genton (2001) "On the discretization of 
!                 nonparametric covariogram estimators"
!
!     PARÁMETROS:
!         nx = nº de nodos                                            (I)
!         x(nx) = nodos                                               (O)
!         dim = dimensión correspondiente                             (I)
!         ((dim-2.0)/2.0 orden de la función de Bessel)
!         rango = máximo salto                                        (I)
!     ------------------------------------------------------------------
      SUBROUTINE disc_sbv(nx, x, dim, rango)
      IMPLICIT NONE
!     Parámetros
      INTEGER nx, dim, i
      REAL(8) x(nx), rango, a
!     En el caso infinito se toman equidistantes
      IF (dim.LE.0) THEN
          DO i = 1, nx
              x(i) = i*0.3
          END DO
      ELSE
          a =(dim-2.0)/2.0
          CALL besselzeros(nx, a, x)
      END IF
!     Se reescala por máximo salto
      DO i = 1,nx
          x(i) = x(i)/rango
      END DO
      RETURN
      END SUBROUTINE disc_sbv


!     ------------------------------------------------------------------
!     [Besselzeros]   Calcula los ceros de una función de Bessel de orden
!                     real utilizando el algoritmo propuesto por J.S. Ball
!                     (2000) "Automatic computation of zeros of Bessel
!                     functions and other special functions", 1458-1464,
!                     J. Sci. Comput.
!
!     PARÁMETROS:
!         nt = nº de ceros                                            (I)
!         a = orden de la función de Bessel                           (I)
!         c(nt) = ceros (ordenados)                                   (O)
!     ------------------------------------------------------------------
      SUBROUTINE besselzeros(nt, a, c)
      IMPLICIT NONE
      INTEGER nt
      REAL(8) a, c(nt)
!     Variables locales
      INTEGER i, j, nmax, ierr       
      REAL(8) fl, a1, aa
      REAL(8), ALLOCATABLE :: e(:), d(:), z(:,:)
!     Asignar memoria a variables locales
      nmax = 2*nt
      ALLOCATE (e(nmax+1), d(nmax), z(nmax, nmax),STAT=i)
      IF (i.NE.0) CALL Error(i,'besselzeros: ALLOCATE')
!     Valores iniciales
      z=0.0D0
      do 10 i=1,nmax
          fl=dfloat(i)+.5d0*a+.5d0
          d(i)=.125d0/fl/(fl-1.d0)
          e(i+1)=.125d0/fl/dsqrt(4.d0*fl**2-1.d0)
          z(i,i)=1.0D0
   10 continue
!     input e(i)=alpha, d(i)= beta; tql2 returns eigenvalues
!     of symmetric tridiagonal matrix in d(i); unsorted
      call tql2(nmax,nmax,d,e,z,ierr)
!     form zeros by inverse of g(x)
      do 20 i=1,nmax
          d(i)=1.d0/dsqrt(d(i))
   20 continue
!     sort zeros
      do 30 i=1,nmax
          a1=d(i)
          do 25 j=i+1,nmax
              if(d(j).lt.a1)then
              aa=d(j)
              d(j)=a1
              a1=aa
              endif
   25     continue
          d(i)=a1
   30 continue
      do 40 i=1,nt
          c(i) = d(i)
   40 continue
!     Liberar memoria de variables locales
      DEALLOCATE (e,d,z,STAT=i)
      IF (i.NE.0) CALL Error(i,'besselzeros: DEALLOCATE')
      RETURN
      END SUBROUTINE Besselzeros

