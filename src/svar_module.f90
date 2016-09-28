!-----------------------------------------------------------------------
!   [svar_module.f90]   Modulo para la estimacion lineal local de un
!                       semivariograma multidimensional
!
!   Interfaces con R:
!       svar_iso_np     Estimador np del svar y rejilla binning (R "svarisonp")
!       svar_iso_bin    Rejilla binning y estimación clásica/robusta (R "svariso")
!
!   Autor: (c) Ruben Fernandez-Casal                Creacion: Ago 2012
!   Revision: Jun 2013, Nov 2013
!-----------------------------------------------------------------------


!   --------------------------------------------------------------------
    subroutine set_bin_svar_iso(g, nd, x, ny, y, nlags, minlag, maxlag, itipo)
!   --------------------------------------------------------------------
!       Establece la rejilla binning (lineal) para la estimación np de un
!       semivariograma isotrópico.
!           g       = rejilla binning (type(grid_bin))
!           nlag    = número de saltos
!           minlag  = mínimo salto (si <0 se toma el valor por defecto maxlag/nlags)
!           maxlag  = máximo salto
!           itipo   = Tipo de estimador a calcular
!                   0 = promedio de las diferencias al cuadrado (equivalente al estimador clásico)
!                   1 = promedio de las diferencias absolutas
!                   2 = reescalado del promedio de las diferencias absolutas  (equivalente al estimador robusto)
!       NOTAS:
!           - Se ignoran los saltos fuera del rango [g%min-g%lag, g%max+g%lag]
!           - g%ny = Nº total de pares = (ny*(ny-1))/2 si se consideraran todos los posibles
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    type(grid_bin) :: g
    integer nd, ny, nlags, itipo
    real(8)  x(nd,ny), y(ny), minlag, maxlag, lag
    integer ib, i, j
    real(8)  wl, wr, tmp
    real(8), external :: DNRM2   ! DNRM2(N,X,INCX)
!   --------------------------------------------------------------------
        if (minlag < 0) minlag = maxlag/nlags
!       Establecer rejilla binning
        call set_grid1d(g, nlags, minlag, maxlag)
        lag = g%lag(1)
!       Asignar memoria rejilla binning
        allocate(g%y(g%ngrid), g%w(g%ngrid))
        g%y = 0.0d0
        g%w = 0.0d0
        g%ny = 0
!       Recorrer datos
        do i = 1, ny-1
            do j = i+1, ny
                tmp = DNRM2(nd, x(1:nd, i) - x(1:nd, j), 1)
                ib = 1 + INT( (tmp - minlag)/lag )
                if ((ib < 0).or.(ib > nlags)) cycle
!               Calculo de los pesos
!               wr = (tmp-minlag-(ib-1)*lag)/lag
!               wl = 1.0d0-wr
                wl = (minlag + ib*lag - tmp)/lag
                wr = 1.0d0 - wl
!               Actualizar valores
                if (itipo > 0) then
                    tmp = DSQRT(DABS( y(i) - y(j) ))   ! Robusto
                else
                    tmp = 0.5d0 * (y(i) - y(j))**2.0d0   ! Clásico
                end if
                if (ib > 0) then
                    g%y(ib) = g%y(ib) + wl*tmp
                    g%w(ib) = g%w(ib) + wl
                end if
                if (ib < nlags) then
                    g%y(ib+1) = g%y(ib+1) + wr*tmp
                    g%w(ib+1) = g%w(ib+1) + wr
                end if
            end do
        end do  ! Recorrer datos
!       Promediar y calcular valor medio
        g%med = 0.0d0
        tmp = SUM(g%w)
        g%ny = INT(tmp)
        do i = 1, g%ngrid
            if (g%w(i) > 0.d0) then
                g%med = g%med + g%y(i)/tmp
                g%y(i) = g%y(i)/g%w(i)
                if (itipo == 2) g%y(i) = 0.5d0 *(g%y(i)**4d0)/(0.457d0+0.494d0/g%w(i))
            end if
        end do
    return
    end subroutine set_bin_svar_iso


!   --------------------------------------------------------------------
    subroutine svar_iso_np( nd, x, ny, y, nlags, minlag, maxlag,                &
   &                        bin_lag, bin_med, bin_y, bin_w,                     &
   &                        h, lpe, degree, ideriv, deriv, ihat, hatlp,         &
   &                        ndelcv, rm, rss, nrl0)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "svarisonp"
!       Devuelve estimador np y rejilla binning
!
!   PENDENTE:
!       - Implementar itipo (actualmente = 0)
!       - Opción/nueva rutina para establecer rejilla binning a partir de parámetros
!         (calculos manteniendo datos)
!   OJO: rejilla unidimensional (e.g. ndelcv(1))
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, ny, nlags, ndelcv, degree, ideriv, ihat, nrl0
    real(8)  x(nd,ny), y(ny), minlag, maxlag
    real(8)  bin_lag, bin_med, bin_y(nlags), bin_w(nlags)
    real(8)  h, lpe(nlags), rm, rss, deriv(nlags), hatlp(nlags,nlags)
    type(grid_bin) :: bin
    real(8), EXTERNAL :: KTWMD
!   --------------------------------------------------------------------
!       subroutine set_bin_svar_iso(g, nd, x, ny, y, nlags, minlag, maxlag, itipo)
        call set_bin_svar_iso(bin, nd, x, ny, y, nlags, minlag, maxlag, 0)    ! Establece la rejilla binning (lineal)
!       Estimación y obtención matriz Hat
!       ndelcv = NCV
        call lp(  bin, h, KTWMD, .true., lpe, degree,                           &
                    ideriv == 1, deriv, nlags, ihat == 1, hatlp, nlags,         &
   &                ndelcv, rm, rss, nrl0)
        where ( lpe < 0 )   lpe = 0.0d0     ! Cuidado con hatlp
!       Copiar resultados
        bin_lag = bin%lag(1)
        bin_med = bin%med
        bin_y = bin%y
        bin_w = bin%w
!       call bin%end_bin
        call end_grid_bin(bin)
    return
    end subroutine svar_iso_np


!   --------------------------------------------------------------------
    subroutine svar_iso_bin(nd, x, ny, y, nlags, minlag, maxlag, itipo,         &
   &                        bin_lag, bin_med, bin_y, bin_w)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "svariso"
!       Devuelve la rejilla binning (lineal) para la estimación np de un
!       semivariograma isotrópico
!       Se puede emplear para estimación clásica/robusta.
!           itipo   = Tipo de estimador a calcular
!                   0 = promedio de las diferencias al cuadrado (equivalente al estimador clásico)
!                   1 = promedio de las diferencias absolutas
!                   2 = reescalado del promedio de las diferencias absolutas  (equivalente al estimador robusto)
!
!   Autor: (c) Ruben Fernandez-Casal    Ultima revision: Ago 2012
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, ny, nlags, itipo
    real(8)  x(nd,ny), y(ny), minlag, maxlag
    type(grid_bin) :: bin
    real(8)  bin_lag, bin_med, bin_y(nlags), bin_w(nlags)
        call set_bin_svar_iso(bin, nd, x, ny, y, nlags, minlag, maxlag, itipo)    ! Establece la rejilla binning (lineal)
        bin_lag = bin%lag(1)
        bin_med = bin%med
        bin_y = bin%y
        bin_w = bin%w
!       call bin%end_bin
        call end_grid_bin(bin)
    return
    end subroutine svar_iso_bin

