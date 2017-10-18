!-----------------------------------------------------------------------
    module grid_module
!-----------------------------------------------------------------------
!   Due to the CRAN policy requirement of using only Fortran 90/95, 
!   the implementation of additional grid types is postponed until Fortran 
!   compilers used at CRAN (specially in the case of Mac OS X) support the 
!   required Fortran 2003 features (mainly type-bound procedures).
!-----------------------------------------------------------------------
!   Modulo clases rejilla (grid) y rejilla binning (grid_bin)
!       type(grid) :: g
!           g%ndim          = Nº de dimensiones
!           g%n(1:g%ndim)   = Nº de nodos en cada dimensión
!           g%ngrid         = Nº total de nodos (= PRODUCT(n) rejilla estándar)
!           g%min(1:g%ndim) = Mínimo de la rejilla binning en cada dimensión
!           g%max(1:g%ndim) = Máximo      ""
!           g%lag(1:g%ndim) = Espaciado   ""
!           g%i             = indice unidimensional
!           g%ii(1:g%ndim)  = indice multidimensional
!           procedure :: set => set_grid, end => end_grid, ind, set_ind, incii
!       type(grid_den) :: bin
!           bin%w(1:bin%ngrid)  = Peso/frecuencia nodo binning (unidim.)
!           bin%ny              = Nº de datos (suma de los pesos binning)
!           procedure :: set_bin_den, end_bin_den
!       type(grid_bin) :: bin
!           bin%y(1:bin%ngrid)  = Valor nodo binning (indice unidimensional)
!           bin%med             = Media ponderada de bin%y
!                                 (media de los datos & binning; para estimación aditiva)
!           procedure :: set_bin => set_grid_bin, end_bin => end_grid_bin
!       set_grid1d
!
!   Interfaces con R:
!       binning_r           Rejilla binning ("binning" en "bin.data.R")
!       interp_data_grid    Interpolacion lineal de una rejilla ("interp.grid.par" en "interp.R")
!
!   PENDENTE:
!       - Crear clase grid_den intermedia entre grid y grid_bin
!         (renombrar %ny %nw? %nx?)  
!       - Implementar ndim como parámetro LEN de type grid
!       - Crear clase grid_data (intermedia entre grid y grid_bin?)
!       - procedure :: end type-bound / final
!
!   Autor: (c) Ruben Fernandez-Casal       
!   Fecha revision: Oct 2012, Nov 2013
!-----------------------------------------------------------------------
        implicit none

!       ----------------------------------------------------------------
        type grid_bin
!       ----------------------------------------------------------------
!       Rejilla binning multidimensional
!       ----------------------------------------------------------------
!       type grid
            integer ndim, ngrid, igrid
            integer, allocatable :: n(:), ii(:)
!           n(g%ndim) = NBM(1..NDimBM)= Nº nodos de la rejilla binning
            real(8), allocatable :: min(:), max(:), lag(:)
!       type, extends(grid) :: grid_den
            integer ny    !NBMc
            real(8), allocatable :: w(:)
!       type, extends(grid_den) :: grid_bin
            real(8) :: med
            real(8), allocatable :: y(:)
        end type

!   --------------------------------------------------------------------
    contains
!   --------------------------------------------------------------------

!       ----------------------------------------------------------------
        subroutine set_grid(g, ndim, n, min, max)
!       ----------------------------------------------------------------
!       Establece la rejilla
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer ndim, n(ndim)
        real(8) min(ndim), max(ndim)
!       ----------------------------------------------------------------
            g%ndim = ndim
            allocate(g%n(ndim), g%ii(ndim), g%min(ndim), g%max(ndim), g%lag(ndim))
            g%n = n
            g%ngrid = PRODUCT(n)
            g%min = min
            g%max = max
            g%lag = (max - min)/(n - 1.0d0)
        return
        end subroutine set_grid


!       ----------------------------------------------------------------
        subroutine set_grid1d(g, n, min, max)
!       Establece una rejilla unidimensional
!       Avoid rank mismatch in argument(rank-1 and scalar)
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer n
        real(8) min, max
!       ----------------------------------------------------------------
            g%ndim = 1
            allocate(g%n(1), g%ii(1), g%min(1), g%max(1), g%lag(1))
            g%n(1) = n
            g%ngrid = n
            g%min(1) = min
            g%max(1) = max
            g%lag(1) = (max - min)/(n - 1.0d0)
        return
        end subroutine set_grid1d


!       ----------------------------------------------------------------
        subroutine end_grid(g)
!       Libera memoria
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
!       ----------------------------------------------------------------
            deallocate(g%n, g%ii, g%min, g%max, g%lag)
        return
        end subroutine end_grid


!       ----------------------------------------------------------------
        integer function ind(g, ii)
!       Devuelve el indice unidimensional equivalente a uno
!       multidimensional. Por ejemplo, en el caso tridimensional:
!       iibm  = ix + (iy-1)*nbmx + (iz-1)*nbmx*nbmy
!             = ix + ((iy-1) + (iz-1)*nbmy)*nbmx
!       (se evalua empleando una regla tipo Horner)
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer ii(g%ndim)
        integer i, k
!       ----------------------------------------------------------------
            i = 0
            do k = g%ndim, 2, -1 
                i = g%n(k-1)*(i + ii(k) - 1)
            end do
            ind = i + ii(1)
        return
        end function ind


!       ----------------------------------------------------------------
        subroutine set_ind(g, ii)
!       Establece el indice unidimensional y multidimensional, a partir
!       de un indice multidimensional
!       PENDIENTE: VERIFICAR RANGO INDICES
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer ii(g%ndim)
            g%ii = ii
            g%igrid = ind(g, ii)
        return
        end subroutine set_ind


!       ----------------------------------------------------------------
        subroutine incii(g)
!       Incrementa el indice (uni y multi dimensional)
!       NO VERIFICA RANGOS
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer j
!           ------------------------------------------------------------
            g%igrid = g%igrid + 1
            do j = 1, g%ndim
                g%ii(j) = g%ii(j) + 1
                if (g%ii(j) <= g%n(j)) return
                g%ii(j) = 1
            end do
            return
        end subroutine incii


!       ----------------------------------------------------------------
        subroutine set_bin_den(g, nd, nbin, x, ny)
!       Establece la rejilla binning (lineal) para densidad
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer nd, nbin(nd), ny
        real(8)  x(nd,ny)
!
        integer ii(nd), iib(nd), ib, i, j, k
        real(8)  minx(nd), maxx(nd), w(2,nd), tmp
        integer niinc, iinc(nd, 2**nd)
!       real(8), external :: dmach
!           ------------------------------------------------------------
!           Calcular dimensiones rejilla binning
            minx = x(1:nd, 1)
            maxx = minx
            do j = 1, nd
                do i = 2, ny
                    if (minx(j) > x(j, i)) then
                        minx(j) = x(j, i)
                    else if (maxx(j) < x(j, i)) then
                        maxx(j) = x(j, i)
                    end if
                end do
!               Expandir un poco                
                tmp = MAX(maxx(j)-minx(j), DABS(minx(j)))* (1.0d2*epsilon(1.0d0)) 
                minx(j) = minx(j) - tmp
                maxx(j) = maxx(j) + tmp
            end do
!           Establecer rejilla
            call set_grid(g, nd, nbin, minx, maxx)
!           Asignar memoria rejilla binning
            allocate(g%w(g%ngrid))
            g%ny = ny
            g%w = 0.0d0
!           Indice rejilla actualización
            niinc = 2**nd
            ii = 0
            do i = 1, niinc
                do j = 1, nd - 1
                    if (ii(j) > 1) then
                        ii(j) = 0
                        ii(j+1) = ii(j+1) + 1
                    else
                        exit
                    end if
                end do
                iinc(1:nd, i) = ii(1:nd)
                ii(1) = ii(1) + 1
            end do
!           Recorrer datos
            do i = 1, ny
                do j = 1, nd
                    iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!                   calculo de los pesos
                    w(2, j) = (x(j, i) - g%min(j) - (iib(j)-1)*g%lag(j)) / g%lag(j)
                    w(1, j) = 1.0d0 - w(2, j)
                end do
!               Actualizar valores
                do k = 1, niinc
                    tmp = 1.0d0
                    do j = 1, nd
                        ii(j) = iib(j) + iinc(j, k)
                        tmp = tmp * w(iinc(j, k) + 1, j)
                    end do
                    ! if (tmp < epsilon) cycle
!                   ib = g%ind(ii)
                    ib = ind(g, ii)
                    g%w(ib) = g%w(ib) + tmp
                end do
            end do
        return
        end subroutine set_bin_den


!       ----------------------------------------------------------------
        subroutine end_bin_den(g)
!       Libera memoria
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
!           call g%end
            call end_grid(g)
            deallocate(g%w)
        return
        end subroutine end_bin_den

        
!       ----------------------------------------------------------------
        subroutine set_grid_bin(g, nd, nbin, x, ny, y)
!       Establece la rejilla binning (lineal)
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
        integer nd, nbin(nd), ny
        real(8)  x(nd,ny), y(ny)
!
        integer ii(nd), iib(nd), ib, i, j, k
        real(8)  minx(nd), maxx(nd), w(2,nd), tmp
        integer niinc, iinc(nd, 2**nd)
!       real(8), external :: dmach
!           ------------------------------------------------------------
!           Calcular dimensiones rejilla binning
            minx = x(1:nd, 1)
            maxx = minx
            do j = 1, nd
                do i = 2, ny
                    if (minx(j) > x(j, i)) then
                        minx(j) = x(j, i)
                    else if (maxx(j) < x(j, i)) then
                        maxx(j) = x(j, i)
                    end if
                end do
!               Expandir un poco                
                tmp = MAX(maxx(j)-minx(j), DABS(minx(j)))* (1.0d2*epsilon(1.0d0)) 
                minx(j) = minx(j) - tmp
                maxx(j) = maxx(j) + tmp
            end do
!           Establecer rejilla
            call set_grid(g, nd, nbin, minx, maxx)
!           Asignar memoria rejilla binning
            allocate(g%y(g%ngrid), g%w(g%ngrid))
            g%ny = ny
            g%y = 0.0d0
            g%w = 0.0d0
!           Indice rejilla actualización
            niinc = 2**nd
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
!           Recorrer datos
            do i = 1, ny
                do j = 1, nd
                    iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!                   calculo de los pesos
                    w(2, j) = (x(j, i) - g%min(j) - (iib(j)-1)*g%lag(j)) / g%lag(j)
                    w(1, j) = 1.0d0 - w(2, j)
                end do
!               Actualizar valores
                do k = 1, niinc
                    tmp = 1.0d0
                    do j = 1, nd
                        ii(j) = iib(j) + iinc(j, k)
                        tmp = tmp * w(iinc(j, k) + 1, j)
                    end do
                    ! if (tmp < epsilon) cycle
!                   ib = g%ind(ii)
                    ib = ind(g, ii)
                    g%y(ib) = g%y(ib) + tmp * y(i)
                    g%w(ib) = g%w(ib) + tmp
                end do
            end do
!           promediar y calcular valor medio
            g%med = 0.0d0
            do i = 1, g%ngrid
                if (g%w(i) > 0.d0) then
                    g%med = g%med + g%y(i)/dble(g%ny)
                    g%y(i) = g%y(i) / g%w(i)
                end if
            end do
        return
        end subroutine set_grid_bin


!       ----------------------------------------------------------------
        subroutine end_grid_bin(g)
!       Libera memoria
!       ----------------------------------------------------------------
        implicit none
        type(grid_bin) :: g
!           call g%end
            call end_grid(g)
            deallocate(g%y, g%w)
        return
        end subroutine end_grid_bin
        
!   --------------------------------------------------------------------
    end module grid_module
!   --------------------------------------------------------------------



!   --------------------------------------------------------------------
    subroutine binning_r(nd, nbin, x, ny, y, bin_min, bin_max, bin_med, bin_y, bin_w)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "binning"
!       Devuelve la rejilla binning type(grid_bin)%set_bin(nd, nbin, x, ny, y)
!
!   Autor: (c) Ruben Fernandez-Casal    Ultima revision: Jun 2012
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ny
    real(8)  x(nd,ny), y(ny)
    real(8)  bin_min(nd), bin_max(nd), bin_med, bin_y(*), bin_w(*)
    type(grid_bin) :: bin
!       call bin%set_bin(nd, nbin, x, ny, y) ! Establece la rejilla binning (lineal)
        call set_grid_bin(bin, nd, nbin, x, ny, y)
        bin_min = bin%min
        bin_max = bin%max
        bin_med = bin%med
        bin_y(1:bin%ngrid) = bin%y
        bin_w(1:bin%ngrid) = bin%w
!       call bin%end_bin
        call end_grid_bin(bin)
    return
    end subroutine binning_r


!   --------------------------------------------------------------------
    subroutine bin_den(nd, nbin, x, ny, bin_min, bin_max, bin_w)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "bin.den"
!       Devuelve la rejilla binning type(grid_den)%set_bin_den(nd, nbin, x, ny, y)
!
!   Autor: (c) Ruben Fernandez-Casal    Ultima revision: Nov 2013
!   --------------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ny
    real(8)  x(nd,ny)
    real(8)  bin_min(nd), bin_max(nd), bin_w(*)
    type(grid_bin) :: bin
!       call bin%set_bin_den(nd, nbin, x, ny) ! Establece la rejilla binning (lineal)
        call set_bin_den(bin, nd, nbin, x, ny)
        bin_min = bin%min
        bin_max = bin%max
        bin_w(1:bin%ngrid) = bin%w
!       call bin%end_bin_den
        call end_bin_den(bin)
    return
    end subroutine bin_den

!   --------------------------------------------------------------------
    subroutine interp_data_grid( nd, nbin, bin_min, bin_max,     &
   &                    ngrid, gy, x, ny, y)
!   --------------------------------------------------------------------
!       Interfaz para la rutina de R "interp.grid.par"
!       Interpolacion lineal de una rejilla
!       Establece type(grid) :: g a partir de parámetros
!   --------------------------------------------------------------------
!       PENDENTE: FALLA CON DATOS MISSING
!   ----------------------------------------------------------------
    use grid_module
    implicit none
    integer nd, nbin(nd), ngrid, ny
    real(8)  bin_min(nd), bin_max(nd)
    type(grid_bin) :: g
    real(8)  x(nd,ny), gy(ngrid), y(ny)
!       ----------------------------------------------------------------
!       Establecer la rejilla
!       call g%set(nd, nbin, bin_min, bin_max)
        call set_grid(g, nd, nbin, bin_min, bin_max)
!       Interpolación
        call interp_grid(g, gy, x, ny, y)
!       call g%end
        call end_grid(g)
    return
    end subroutine interp_data_grid



!   ----------------------------------------------------------------
    subroutine interp_grid(g, gy, x, ny, y)
!       Interpola linealmente valores en una rejilla
!   ----------------------------------------------------------------
!       PENDENTE: FALLA CON DATOS MISSING
!   ----------------------------------------------------------------
    use grid_module
    implicit none
    type(grid_bin) :: g
    integer ny
    real(8)  x(g%ndim, ny), gy(g%ngrid), y(ny)
!
    integer nd, ii(g%ndim), iib(g%ndim), ib, i, j, k
    real(8)  w(2,g%ndim), tmp
    integer niinc, iinc(g%ndim, 2**g%ndim)
!       ------------------------------------------------------------
        nd = g%ndim
        niinc = 2**nd
!       Indice rejilla actualización
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
        y = 0.0d0
        do i = 1, ny
            do j = 1, nd
                iib(j) = 1 + int((x(j, i) - g%min(j)) / g%lag(j))
!               Extrapolación:
                if (iib(j) < 1) iib(j) = 1
                if (iib(j) >= g%n(j)) iib(j) = g%n(j) - 1
!               Calculo de los pesos
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
!               ib = g%ind(ii)
                ib = ind(g, ii)
                y(i) = y(i) + tmp * gy(ib)
            end do
        end do
    return
    end subroutine interp_grid




!   --------------------------------------------------------------------
!   NOTAS SOBRE VARIABLES PARA CONVERSIÓN DE CÓDIGO SKUET
!       type(grid_bin) :: bin
!           bin%min             !   MinBM(1..NDimBM)= Mínimo de la rejilla binning
!                                                       en cada dimensión
!           bin%max             !   MaxBM(1..NDimBM)= Máximo ""
!           bin%lag             !   LBM(1..NDimBM)= Espaciado ""
!           bin%y(1:bin%ngrid)  !   BMy(1..NTotBM)= Valor nodo binning
!                                                       (indice unidimensional)
!           bin%w(1:bin%ngrid)  !   BMc(1..NTotBM)= Peso/frecuencia nodo binning
!                                                       (indice unidimensional)
!           bin%med             !   MedBMy= Media ponderada de BMy
!                                                       (media de los datos binning)
!       iBM(1:ND,i)= indice multidimensional correspondiente al indice unidimensional i
!   --------------------------------------------------------------------
!           print *, g%ndim
!           print *, g%igrid
!           print *, g%ii
!           print *, g%n
!           print *, g%ngrid



