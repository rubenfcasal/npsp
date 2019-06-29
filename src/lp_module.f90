!-----------------------------------------------------------------------
!   Due to the CRAN policy requirement of using only Fortran 90/95, 
!   the implementation of additional grid types is postponed until Fortran 
!   compilers used at CRAN (specially in the case of Mac OS X) support the 
!   required Fortran 2003 features (mainly type-bound procedures).
!-----------------------------------------------------------------------
!   [lp_module.f90]   Modulo para la estimacion lineal/polinomica local
!                       multidimensional
!   Interfaces con R (en "locpol.bin.R"):
!       lp_raw        Estimador y rejilla binning ("locpol.default")
!       lp_bin        Estimador a partir de una rejilla binning ("locpol.bin.data")
!       lp_data_grid  Estimador a partir de datos en una rejilla ("locpol.data.grid", "np.den.bin.den")  
!       predict_locpol    Estimacion y matriz hat en observaciones ("predict.locpol.bin")
!
!   Autor: (c) Ruben Fernandez-Casal    
!   Fecha revision: Abr 2007, Oct 2012, Jun 2013, Oct 2013, Mar 2014, Dic 2015
!-----------------------------------------------------------------------

!   --------------------------------------------------------------------
    subroutine lp_raw( nd, nbin, ntbin, x, ny, y,                               &
   &                    bin_min, bin_max, bin_med, bin_y, bin_w,                &
   &                    h, lpe, degree, ideriv, deriv, ihat, hatlp,             &
   &                    ncv, rm, rss, nrl0)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "locpol.default"
!       Devuelve estimador y rejilla binning
!
!   IMPORTANTE:
!       - Se supone que lpe contiene NA's en la entrada
!       - Se supone que deriv contiene NA's en la entrada
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ntbin, ny
    real(8)  x(nd,ny), y(ny)
    real(8)  bin_min(nd),  bin_max(nd), bin_med, bin_y(ntbin), bin_w(ntbin)
    type(grid_bin) :: bin
    integer degree, ideriv, ihat, ncv, nrl0
    real(8)  h(ND,ND), lpe(ntbin), rm, rss, deriv(ntbin, nd), hatlp(ntbin, ntbin)
    integer NDelCV(ND)
    real(8), external :: KTWMD
!   --------------------------------------------------------------------
!       call bin%set_bin(nd, nbin, x, ny, y) ! Establece la rejilla binning (lineal)
        call set_grid_bin(bin, nd, nbin, x, ny, y)
!       Estimacion y obtencion matriz Hat
        NDelCV = ncv
!       SUBROUTINE lp(bin, h, FNucMD, GetE, lpe, degree, GetDERIV, deriv, ldderiv,  &
!   &                GetHAT, hatlp, ldhatlp, NDelCV, RMNP, RSSNP, nrl0)  
        call lp(  bin, h, KTWMD, .true., lpe, degree,                           &
                    ideriv == 1, deriv, ntbin, ihat == 1, hatlp, ntbin,         &
   &                NDelCV, rm, rss, nrl0)
        bin_min = bin%min
        bin_max = bin%max
        bin_med = bin%med
        bin_y(1:bin%ngrid) = bin%y
        bin_w(1:bin%ngrid) = bin%w
!       call bin%end_bin
        call end_grid_bin(bin)

    return
    end subroutine lp_raw


!   --------------------------------------------------------------------
    subroutine lp_bin( nd, nbin, ntbin, bin_min, bin_max, bin_med, bin_y, bin_w,  &
   &                    h, lpe, degree, ideriv, deriv, ihat, hatlp,               &
   &                    ncv, rm, rss, nrl0)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "locpol.bin.data"
!       Devuelve estimador a partir de una rejilla binning
!       Establece type(grid_bin)::bin a partir de parametros
!       (calculos manteniendo datos)
!
!   IMPORTANTE:
!       - Se supone que lpe contiene NA's en la entrada
!       - Se supone que deriv contiene NA's en la entrada
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ntbin
    real(8)  bin_min(nd),  bin_max(nd), bin_med, bin_y(ntbin), bin_w(ntbin)
    type(grid_bin) :: bin
    integer degree, ideriv, ihat, ncv, nrl0
    real(8)  h(nd,nd), lpe(ntbin), rm, rss, deriv(ntbin, nd), hatlp(ntbin,ntbin)
    integer NDelCV(nd)
    real(8), external :: KTWMD
!   --------------------------------------------------------------------
!       Establecer la rejilla
!       call bin%set(nd, nbin, bin_min, bin_max)
        call set_grid(bin, nd, nbin, bin_min, bin_max)
!       Asignar memoria rejilla binning
        allocate(bin%y(bin%ngrid), bin%w(bin%ngrid))
        bin%med = bin_med
        bin%y = bin_y(1:bin%ngrid)
        bin%w = bin_w(1:bin%ngrid)
        bin%ny = INT(SUM(bin%w))
!       Estimacion y obtencion matriz Hat
        NDelCV = ncv
        call lp(bin, h, KTWMD, .true., lpe, degree,                             &
                    ideriv == 1, deriv, ntbin, ihat == 1, hatlp, ntbin,         &
   &                NDelCV, rm, rss, nrl0)
!       call bin%end_bin
        call end_grid_bin(bin)
    return
    end subroutine lp_bin


!   --------------------------------------------------------------------
    subroutine lp_data_grid( nd, nbin, ntbin, bin_min, bin_max, bin_med, bin_y, &
   &                    h, lpe, degree, ideriv, deriv, ihat, hatlp,             &
   &                    ncv, rm, rss, nrl0)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "locpol.data.grid" 
!       Interfaz para la rutina de R "np.den.bin.den"    
!       Devuelve estimador a partir de datos en una rejilla regular (bin_w = 1)
!       Establece type(grid_bin)::bin a partir de parametros
!       (calculos manteniendo datos)
!
!   IMPORTANTE:
!       - Se supone que lpe contiene NA's en la entrada (locpol.data.grid)
!         o 0's para estimacion densidad (np.den.bin.den)
!       - bin_med se utiliza para el calculo de las medidas de error 
!         cuando NRL <  NINDRL. Debe ser 0 para estimacion densidad 
!       - Se supone que deriv contiene NA's en la entrada
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ntbin
    real(8)  bin_min(nd),  bin_max(nd), bin_med, bin_y(ntbin)
    type(grid_bin) :: bin
    integer degree, ideriv, ihat, ncv, nrl0
    real(8)  h(nd,nd), lpe(ntbin), rm, rss, deriv(ntbin, nd), hatlp(ntbin,ntbin)
    integer NDelCV(nd)
    real(8), external :: KTWMD
!   --------------------------------------------------------------------
!       Establecer la rejilla
!       call bin%set_bin(nd, nbin, x, ny, y)
!       call bin%set(nd, nbin, bin_min, bin_max)
        call set_grid(bin, nd, nbin, bin_min, bin_max)
!       Asignar memoria rejilla binning
        allocate(bin%y(bin%ngrid), bin%w(bin%ngrid))
        bin%med = bin_med
        bin%y = bin_y(1:bin%ngrid)
        bin%w = 1.0d0
        bin%ny = bin%ngrid
!       Estimacion y obtencion matriz Hat
        NDelCV = ncv
!       SUBROUTINE lp(bin, h, FNucMD, GetE, lpe, degree, GetDERIV, deriv, ldderiv,  &
!   &                GetHAT, hatlp, ldhatlp, NDelCV, RMNP, RSSNP, nrl0)
        call lp(bin, h, KTWMD, .true., lpe, degree,                             &
                    ideriv == 1, deriv, ntbin, ihat == 1, hatlp, ntbin,         &
   &                NDelCV, rm, rss, nrl0)
!       call bin%end_bin
        call end_grid_bin(bin)
    return
    end subroutine lp_data_grid




!   --------------------------------------------------------------------
!   [KEpanMD] Nucleo de Epanechnikov multidimensional producto
!   --------------------------------------------------------------------
      real(8) FUNCTION KEpanMD(u, ndim)
      IMPLICIT NONE
      integer ndim, i
      real(8) u(ndim), tmp
      real(8), external :: KEpan
!   --------------------------------------------------------------------
      tmp = 1.0d0
      DO i = 1, ndim
          tmp = tmp*KEpan(u(i))
      END DO
      KEpanMD = tmp
      RETURN
      END FUNCTION KEpanMD


!   --------------------------------------------------------------------
!   [KEpan] Nucleo de Epanechnikov
!   --------------------------------------------------------------------
      real(8) FUNCTION KEpan(u)
      IMPLICIT NONE
      real(8) u, tmp
!   --------------------------------------------------------------------
      tmp = 1.0d0 - u**2.0d0
      IF (tmp <= 0.0d0) THEN
          KEpan = 0.0d0
      ELSE
          KEpan = 3.0d0*tmp/4.0d0
      ENDIF
      RETURN
      END FUNCTION KEpan


!   --------------------------------------------------------------------
!   [KTWMD] Nucleo Triweight multidimensional producto
!   --------------------------------------------------------------------
      real(8) FUNCTION KTWMD(u, ndim)
      IMPLICIT NONE
      integer ndim, i
      real(8) u(ndim), tmp
      real(8), external :: KTW
!   --------------------------------------------------------------------
      tmp = 1.0d0
      DO i = 1, ndim
          tmp = tmp*KTW(u(i))
      END DO
      KTWMD = tmp
      RETURN
      END FUNCTION KTWMD


!   --------------------------------------------------------------------
!   [KTW] Nucleo Triweight
!   --------------------------------------------------------------------
!     \int u^2K(u)du = \frac {1}{9}
!     \int K(u)^2 du = \frac {350}{429}
      real(8) FUNCTION KTW(u)
      IMPLICIT NONE
      real(8) u, tmp
!   --------------------------------------------------------------------
      tmp = 1.0d0 - u**2.0d0
      IF (tmp <= 0.0d0) THEN
          KTW = 0.0d0
      ELSE
          KTW = (35.0d0/32.0d0)*tmp**3
      ENDIF
      RETURN
      END FUNCTION KTW




!   --------------------------------------------------------------------
    SUBROUTINE lp(bin, h, FNucMD, GetE, lpe, degree, GetDERIV, deriv, ldderiv,  &
   &                GetHAT, hatlp, ldhatlp, NDelCV, RMNP, RSSNP, nrl0)
!   --------------------------------------------------------------------
!       Obtiene la estimacion no parametrica lineal local multidimensional
!       a partir de la rejilla binning (lineal)
!
!   Parametros:
!       type(grid_bin) :: bin
!       GetE =.true. para obtener estimaciones
!           (en caso contrario no se accede a lpe)
!       lpe(NBx) = estimaciones polinomico locales en nodos binning
!       degree = grado del polinomio local
!               (0 =Nadaraya-Watson, 1 = lineal local, 2 = cuadratico local)
!       GetHAT  =.true. para construir matriz hat (de datos binning)
!           (en caso contrario no se accede a hatlp)
!       hatlp(NBx,NBx) = matriz hat de datos binning lpe = hatlp %*% bin%y
!       GetDERIV  =.true. para construir matriz de derivadas parciales
!           (en caso contrario no se accede a deriv)
!       deriv(NBx,NDimBM) = matriz de derivadas parciales
!       NDelCV(NDimBM)= Num de ptos a eliminar del vecindario en cada dimension
!               (para validacion cruzada)
!           NDelCV(j)=0 para no eliminar ningun pto
!           NDelCV(j)=1 para eliminar el pto en la posicion de estimacion
!               (validacion cruzada tradicional)
!           NDelCV(j)=k+1 para eliminar el rango [-k,k] en dimension j
!               (p.e. CV con dependencia)
!       nrl0 = numero de nodos binning (con peso no nulo) sin datos suficientes
!
!   PENDIENTE:
!       - Eliminar calculo de deth en DSYTRFI?
!       - Incluir Epsilon como parametro ?
!             tmp = bin%w(i) * Bk(j)  ! bin%w = Peso/frecuencia nodo binning
!             IF (tmp < Epsilon**2.0d0) CYCLE
!       - Incluir Tolerance para evitar extrapolaciones alejadas de los datos?
!             Calcular WSum en el bucle
!             IF (deth * WSum < Tolerance) CYCLE  
!       - Devolver rejilla binning (grid_bin) + rejilla datos
!
!   IMPORTANTE:
!       - Se supone que lpe contiene NA's en la entrada
!         (o 0's para estimacion densidad)
!       - Se supone que deriv contiene NA's en la entrada
!       - Los nodos con peso negativo bin%w(i) < 0.0d0 son ignorados
!
!   Autor: (c) Ruben Fernandez-Casal
!   Fecha revision: Mar 2008, Sep 2012, Jun 2013, Mar 2014, Dic 2015
!   --------------------------------------------------------------------
    use grid_module
    use linreg_module
    implicit none
    type(grid_bin) :: bin
    integer :: nd, ny
    integer nbin(bin%ndim), ngrid
!   --------------------------------------------------------------------
    integer degree, ldhatlp, ldderiv, NDelCV(bin%ndim), nrl0
    real(8)  h(bin%ndim, bin%ndim), lpe(bin%ngrid), RMNP, RSSNP,                     &
            hatlp(ldhatlp, bin%ngrid), deriv(ldderiv, bin%ndim) 
    real(8), external :: FNucMD
    logical GetE, GetHAT, GetDERIV
!   Variables locales
    real(8), PARAMETER :: Epsilon = 1.0D-6
    integer indb(bin%ndim), ii0(bin%ndim), ii(bin%ndim), iinc(bin%ndim), nindb,     &
   &        indy(bin%ngrid), i, i0, i1, j
    real(8)  IH(bin%ndim, bin%ndim), deth, d(bin%ndim), t(bin%ndim), Bk(bin%ngrid),  &
   &        VHAT(bin%ngrid), WSNP, WSum, tmp
    logical LOUT, LDELCV
    integer, ALLOCATABLE :: iindb(:,:)
!   --------------------------------------------------------------------
        nd = bin%ndim
        ngrid = bin%ngrid     ! ntotbm
        nbin = bin%n
        ny = bin%ny
!       Valores iniciales
!       lpe = 0.0D0
        IF (GetHAT) hatlp = 0.0D0
        nrl0 = 0
        RMNP = 0.0D0
        RSSNP = 0.0D0
        WSNP = 0.0D0
!       Calcular el determinante y la inversa de la matriz simetrica
!       Only the upper triangular part of h is referenced.
        CALL DSYTRFI(nd, h, IH, deth)
        IF (DABS(deth) < Epsilon) &
   &                    CALL error(1, 'lp: the bandwith matrix is singular.')
!       Recorrer rejilla y evaluar nucleo
        indb = 1
        bin%ii = 1
!       CALL bin%set_ind(bin%ii)
        bin%igrid = ind(bin, bin%ii)
        DO i = 1, ngrid
!           bin%igrid=bin%ind(ii) ! ii0(1:nd)=iBM(1:nd,i0)
!           Indice multidimensional correspondiente al indice unidimensional i0
            DO j = 1, nd
                d(j) = (bin%ii(j)-1.0) * bin%lag(j)
            END DO           
!           CALL DMURRV(nd, nd, IH, nd, nd, d, 1, nd, t)
            t = matmul(IH, d)
!           tmp = FNucMD(t, nd)/deth  ! problems with large values
            tmp = FNucMD(t, nd)
            Bk(i) = tmp
!           Calcular rango datos estimacion
            IF (tmp > Epsilon) THEN
                DO j = 1, nd
                    IF (bin%ii(j) > indb(j)) indb(j) = bin%ii(j)
                END DO
            END IF
!           CALL bin%incii()    ! Incrementa el indice (uni y multi dimensional)
            CALL incii(bin)
        END DO  !   DO i = 1, ngrid
        indb = indb - 1
        nindb = 1
        DO j = 1, nd
          nindb = nindb*(2*indb(j)+1)
        END DO
        NINDRL = 1 + nd * degree
!       Verificar problemas vecindarios 
!           PENDENTE: NON DATOS EN DIMENSION I
!           PENDENTE: indb(J)-NDelCV(J) < 0 --> WARNING
!           PENDENTE: opcion salir sin hacer nada / warning / stop ?
        IF ( MAXVAL(indb - NDelCV) < 0 )   & 
!               Los vecindarios estan vacios
   &            CALL error(2, 'lp: there is no data in the neighborhoods.')
        IF (nindb < NINDRL)    &
!               No hay suficientes datos para el ajuste propuesto
   &            CALL error(3, 'lp: there is not enough data in neighborhoods.')
!       Preparar rejilla MD con posiciones relativas: iindb
!           iindb indice multidimensional rejilla -indb:indb        
        ALLOCATE (iindb(nd,nindb))
!           It would be advisable to verify that the memory allocation was correct...            
!           , STAT = i) IF (i /= 0) CALL error(-1, 'lp: ALLOCATE (iindb(nd,nindb)).')
!           warning: ... may be used uninitialized in this function [-Wmaybe-uninitialized]
        ii = -indb
        DO i = 1, nindb
          iindb(1:nd,i) = ii(1:nd)
          DO j = 1, nd
              ii(j) = ii(j)+1
              IF (ii(j) <= indb(j)) EXIT
              ii(j) = -indb(j)
          END DO
        END DO
!       Asignar memoria para ajuste lineal
        CALL ModRegLinInit(nindb, NINDRL)
!       Recorrer rejilla binning
        bin%ii = 1
!       CALL bin%set_ind(bin%ii)
        bin%igrid = ind(bin, bin%ii)
        DO i0 = 1, ngrid
!           Ignorar datos enmascarados
            IF (bin%w(i0) < 0.0d0) THEN
!               CALL bin%incii()    ! Incrementa el indice (uni y multi dimensional)
                CALL incii(bin)
                CYCLE
            END IF
!           Indice multidimensional correspondiente al indice unidimensional i0
            ii0 = bin%ii
!           Seleccionar puntos estimacion
            NRL = 0
            DO i1 = 1, nindb
                LOUT = .FALSE.
                LDELCV = .TRUE.
                DO j = 1, nd
                    ii(j) = ii0(j) + iindb(j,i1)
                    IF ((ii(j) < 1).OR.(ii(j) > nbin(j))) THEN
                        LOUT = .TRUE.
                        EXIT
                    END IF
                    iinc(j) = IABS(iindb(j,i1))+1   ! OJO: iinc=1 en salto 0
                    IF (iinc(j) > NDelCV(j)) LDELCV = .FALSE.
                END DO
                IF (LOUT.OR.LDELCV) CYCLE
!               i = bin%ind(ii)                 ! Indice unidimensional equivalente
                i = ind(bin, ii)
!               j = bin%ind(iinc)               ! iiBM(iinc)
                j = ind(bin, iinc)
                tmp = bin%w(i) * Bk(j)          ! bin%w = Peso/frecuencia nodo binning
!               IF (tmp < Epsilon) CYCLE
                IF (tmp < Epsilon**2.0d0) CYCLE
                NRL = NRL+1
                IF (degree > 0) THEN
                    tmp = DSQRT(tmp) 
                    DO j = 1, nd
                        XRL(NRL, j + 1) = iindb(j, i1) * bin%lag(j) * tmp
                    END DO
                    IF (degree > 1) THEN
                        DO j = 1, nd
                            XRL(NRL, j + nd + 1) = XRL(NRL, j + 1)**2.0d0
                        END DO
                    END IF
                END IF    
                XRL(NRL,1) = tmp
                YRL(NRL) = bin%y(i) * tmp
                IF (GetHAT) indy(NRL) = i
            END DO  ! DO i1 = 1, nindb
!           Estimacion
            IF (NRL <  NINDRL) THEN
                IF (bin%w(i0) > 0.0d0) THEN
                    BRL(1) = bin%med    ! Media global para medidas de error                    
                    nrl0 = nrl0 + 1
                    IF (GetHAT) hatlp(i0, 1:ngrid) = 1/ngrid  ! Aprox. (considerar bin%w?)
                END IF    
            ELSE
                IF (degree == 0) THEN
!                   Media ponderada (pesos bin%w(i) * Bk(j))
                    WSum = SUM(XRL(1:NRL,1))
                    BRL(1) = SUM(YRL(1:NRL))/WSum
                    IF (GetE) lpe(i0) = BRL(1)
                    IF (GetHAT) VHAT(1:NRL) = XRL(1:NRL,1)/WSum 
                ELSE 
!                   Calcular estimacion lineal
                    CALL ModRegLinRL       ! OJO: ModRegLinRL FALLA SI NRL < NINDRL
!                   IF (GetE.AND.(RANKRL = NINDRL)) lpe(i0) = BRL(1)
                    IF (GetE) lpe(i0) = BRL(1)
!                   IF (GetHAT.AND.(RANKRL = NINDRL)) CALL GetVHatLP(VHAT)
                    IF (GetHAT) CALL GetVHatLP(VHAT)
!                   PENDIENTE:Controlar la posibilidad de error o warning
!                   LError = (INFORL /= 0) ! Actualmente genera error en ModRegLinRL
                    IF (GetDERIV) THEN
                        DO i = 2, RANKRL
                            IF ( JPVTRL( i ) <= nd + 1 )                          &
   &                            deriv( i0, JPVTRL( i ) - 1 ) = BRL( i )
                        END DO                    
                    END IF
                END IF
                IF (GetHAT) THEN
                    DO i = 1, NRL
                        hatlp(i0,indy(i)) = VHAT(i)
                    END DO
                END IF
            END IF
!           Medidas error
            IF (bin%w(i0) > 0.0d0) THEN
                tmp = bin%y(i0) - BRL(1)
!               RMNP = RMNP + bin%w(i0) * tmp/DBLE(bin%ny)
                WSNP = WSNP + bin%w(i0)
                RMNP = RMNP + bin%w(i0) * tmp
                RSSNP = RSSNP + bin%w(i0) * tmp * tmp
            END IF
!           CALL bin%incii()    ! Incrementa el indice (uni y multi dimensional)
            CALL incii(bin)
        END DO  ! DO i0 = 1, ngrid
        IF (WSNP > 0.0d0) RMNP = RMNP/WSNP
!       Verificar missing values
!        IF (nrl0 > 0) call rwarn('Not enough data in some neighborhoods.')
!   &        call warning('missing values were generated', &
!   &        'lp: some neighborhoods are empty, nrl0 =', nrl0)
!           PENDENTE:  Controlar warning
!       Liberar memoria
        CALL ModRegLinExit
        DEALLOCATE (iindb, STAT = i)
        IF (i /= 0) CALL error(-2, 'lp: DEALLOCATE(iindb).')
    RETURN
    END SUBROUTINE lp



!   ----------------------------------------------------------------
    subroutine predict_locpol_bin(g, lpe, gethat, hatlp, x, lpy, hatlpy)
!   --------------------------------------------------------------------
!       Devuelve las estimaciones y matriz hat correspondiente a las observaciones
!       NOTA: x deberia ser el empleado para construir la rejilla binning type(grid_bin) :: g
!
!   Parametros:
!       type(grid_bin) :: g
!       lpe(NBx) = estimaciones polinomico locales en nodos binning
!       gethat  =.true. para construir matriz hat (de datos originales)
!           (en caso contrario no se accede a hatlpy)
!       hatlp(NBx,NBx) = matriz hat de datos binning lpe = hatlp %*% bin%y
!       x(NBx, ny) = covariables/coordenadas de datos originales 
!       lpy(ny) = estimaciones polinomico locales en datos originales
!       hatlpy(ny,ny) = matriz hat de datos originales lpy = hatlpy %*% y
!
!   PENDIENTE:
!         - Warning si extrapolacion
!         - ALTERNATIVA PARA SPARSE MATRIX PBIN
!   ----------------------------------------------------------------
    use grid_module
    implicit none
    type(grid_bin) :: g
    real(8)  x(g%ndim,g%ny), lpe(g%ngrid), hatlp(g%ngrid, g%ngrid), lpy(g%ny), hatlpy(g%ny, g%ny)
    logical gethat
!
    integer nd, ny, ngrid
    integer ii(g%ndim), iib(g%ndim), ib, i, j, k
    real(8)  w(2,g%ndim), tmp
    integer niinc, iinc(g%ndim, 2**g%ndim)
    integer, allocatable :: ibin(:, :)
    real(8), allocatable  :: ibinw(:, :), phat(:, :)
    real(8), parameter :: epsilon = 1.0D-7
!       ------------------------------------------------------------
        nd = g%ndim
        ny = g%ny
        ngrid = g%ngrid
        niinc = 2**nd
        if (gethat) then
!           Asignar memoria matriz temporal
            allocate(ibin(ny, 2**nd), ibinw(ny, 2**nd), phat(ny, ngrid))
!           It would be advisable to verify that the memory allocation was correct...            
!           , stat = i) if (i /= 0) call error(i, 'predict_locpol_bin: ALLOCATE')
!           warning: ... may be used uninitialized in this function [-Wmaybe-uninitialized]
        else
!           This would not be really necessary, only to avoid        
!           warning: ... may be used uninitialized in this function [-Wmaybe-uninitialized]
            allocate(ibin(1, 1), ibinw(1, 1), phat(1, 1))
        end if    
!       Indice rejilla actualizacion
        ii = 0
        do i = 1, niinc
            do j = 1, nd-1
                if (ii(j) > 1) then
                    ii(j) = 0
                    ii(j+1) = ii(j+1)+1
                else
                    exit
                end if
            end do
            iinc(1:nd, i) = ii(1:nd)
            ii(1) = ii(1) + 1
        end do
!       Recorrer datos
        lpy = 0.0d0
        if (gethat) phat = 0.0d0
        do i = 1, ny
            do j = 1, nd
                iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!               Extrapolacion:
                if (iib(j) < 1) iib(j) = 1
                if (iib(j) >= g%n(j)) iib(j) = g%n(j) - 1
!               calculo de los pesos
                w(2, j) = (x(j, i) - g%min(j) - (iib(j)-1)*g%lag(j)) / g%lag(j)
                w(1, j) = 1.0d0 - w(2, j)
            end do
!           Actualizar valores
            do k = 1, niinc
                tmp = 1.0d0
                do j = 1, nd
                    ii(j) = iib(j) + iinc(j, k)
                    tmp = tmp * w(iinc(j, k) + 1, j)
                end do
                ! Cuidado puede ocurrir ibinw(i, k) = tmp = 0  y g%w(ib) = 0
!               ib = g%ind(ii)
                ib = ind(g, ii)
                lpy(i) = lpy(i) + tmp * lpe(ib)
!--------------------------------------
!               Opcional si matriz hat:               
                if (.not.gethat) cycle
!               Obtencion de phat = pbin2y * hatlp
                do j = 1, ngrid
                    phat(i, j) = phat(i, j) + tmp * hatlp(ib, j)
                end do
                ibin(i, k) = ib
                ibinw(i, k) = tmp
!               pbin2y(i, ibin(i, k)) = ibinw(i, k)
            end do
        end do
!--------------------------------------
!       Opcional si matriz hat:               
        if (.not.gethat) return        
!       Obtencion de hatlpy = phat * py2bin
!       Matriz de proyeccion binning: py2bin(ibin(j, k), j) = ibinw(j, k) / g%w(ibin(j, k))
!       NOTA: g%w(ib) se podria calcular en el codigo anterior...
        hatlpy = 0.0d0
        do j = 1, ny
            do k = 1, niinc
                ib = ibin(j, k)
                if (g%w(ib) > epsilon) then 
                    tmp = ibinw(j, k) / g%w(ib)   
                    do i = 1, ny
                        hatlpy(i, j) = hatlpy(i, j) + phat(i, ib) * tmp
                    end do
                end if        
            end do
        end do
!       liberar memoria
       deallocate (ibin, ibinw, phat, stat = i)
       if (i /= 0) call error(i, 'predict_locpol_bin: DEALLOCATE')
    return
    end subroutine predict_locpol_bin



!   --------------------------------------------------------------------
    subroutine predict_locpol( nd, nbin, bin_min, bin_max, bin_med, bin_y, bin_w,    &
   &                    ngrid, lpe, ihat, hatlp, x, ny, lpy, hatlpy)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "predict.locpol.bin"
!       Devuelve las estimaciones y matriz hat correspondientes a las observaciones
!       Establece type(grid_bin):: bin a partir de parametros
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ngrid, ny, ihat
    real(8)  bin_min(nd),  bin_max(nd), bin_med, bin_y(ngrid), bin_w(ngrid)
    type(grid_bin) :: bin
    real(8)  x(nd,ny), lpe(ngrid), hatlp(ngrid, ngrid), lpy(ny), hatlpy(ny, ny)
!       ----------------------------------------------------------------
!       Establecer la rejilla
!       call bin%set_bin(nd, nbin, x, ny, y)
!       call bin%set(nd, nbin, bin_min, bin_max)
        call set_grid(bin, nd, nbin, bin_min, bin_max)
!       Asignar memoria rejilla binning
        allocate(bin%y(bin%ngrid), bin%w(bin%ngrid))
        bin%med = bin_med
        bin%y = bin_y(1:bin%ngrid)
        bin%w = bin_w(1:bin%ngrid)
        bin%ny = ny
!       Estimacion y obtencion matriz Hat
        call predict_locpol_bin(bin, lpe, ihat == 1, hatlp, x, lpy, hatlpy)
!       call bin%end_bin
        call end_grid_bin(bin)
    return
    end subroutine predict_locpol