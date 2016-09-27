!-----------------------------------------------------------------------
!     [linreg_module.f90]   Modulo para la regresion lineal en la estimación
!                           lineal local multidimensional
!
!     (c) Ruben Fernandez-Casal       Fecha revision: Abr 2007, Abr 2012, Jun 2013
!-----------------------------------------------------------------------
      MODULE linreg_module
      IMPLICIT NONE
!     Variables principales
      INTEGER LDXRL, NINDRL, NRL
      REAL(8), ALLOCATABLE :: XRL(:,:), YRL(:), BRL(:)
!     Variables secundarias
      REAL(8)  RCONDRL
      PARAMETER (RCONDRL = 1.0D-6)
      INTEGER RANKRL, LWKRL, INFORL
      REAL(8), ALLOCATABLE :: HATRL(:,:), RRL(:,:), WKRL(:)
      INTEGER, ALLOCATABLE :: JPVTRL(:)
!
      LOGICAL ModRegLinIni   !Modulo inicializado
      DATA ModRegLinIni /.FALSE./
!     ------------------------------------------------------------------
      CONTAINS

!         --------------------------------------------------------------
!         [ModRegLinInit]   Asignar memoria variables datos.
!         --------------------------------------------------------------
          SUBROUTINE ModRegLinInit(NMaxOBS, NMaxIND)
          INTEGER NMaxOBS, NMaxIND, i
!         NOTA: SE DEBE VERIFICAR NMaxOBS > NMaxIND
!         --------------------------------------------------------------
          IF (ModRegLinIni) CALL ModRegLinExit
!         Valores iniciales
          LDXRL = NMaxOBS
          NINDRL = NMaxIND
!         Asignar memoria
          ALLOCATE (XRL(LDXRL, NINDRL), YRL(LDXRL), BRL(LDXRL),           &
     &            RRL(LDXRL, NINDRL), JPVTRL(NINDRL), STAT = i)
          IF (i /= 0) CALL Error(i)
!         Obtener espacio lapack
          LWKRL = -1
          CALL DGELSYR( LDXRL, NINDRL, 1, RRL, LDXRL, BRL, LDXRL, JPVTRL,  & 
     &                  RCONDRL, RANKRL, BRL(1), LWKRL, INFORL )
          LWKRL = NINT(BRL(1))
          ALLOCATE(WKRL(LWKRL), STAT = i)
          IF (i /= 0) CALL Error(i)
          ModRegLinIni = .TRUE.
          RETURN
          END SUBROUTINE ModRegLinInit



!         --------------------------------------------------------------
!         [ModRegLinExit]   Liberar variables datos.
!         --------------------------------------------------------------
          SUBROUTINE ModRegLinExit
          INTEGER i
!         --------------------------------------------------------------
!         Liberar memoria
          DEALLOCATE (XRL, YRL, BRL, RRL, WKRL, JPVTRL, STAT = i)
          IF (i /= 0) CALL Error(i)
          IF (ALLOCATED(HATRL)) THEN
              DEALLOCATE(HATRL, STAT = i)
              IF (i /= 0) CALL Error(i)
          END IF
          ModRegLinIni = .FALSE.
          RETURN
          END SUBROUTINE ModRegLinExit


!         --------------------------------------------------------------
!         [ModRegLinRL]   Regresión lineal
!         --------------------------------------------------------------
          SUBROUTINE ModRegLinRL
          REAL(8), EXTERNAL :: DDOT
!         --------------------------------------------------------------
!         Generar error si NRL <  NINDRL
          IF( NRL <  NINDRL ) CALL Error(1, 'ModRegLinRL: NRL <  NINDRL') 
!         Establecer valores
          RRL(1:NRL, 1:NINDRL) = XRL(1:NRL, 1:NINDRL)
          BRL(1:NRL) = YRL(1:NRL)
          JPVTRL = 0
          JPVTRL(1) = 1 ! the 1ST column of X to the front of XP
!         Realizar regresión
          CALL DGELSYR( NRL, NINDRL, 1, RRL, LDXRL, BRL, LDXRL, JPVTRL,  & 
     &                  RCONDRL, RANKRL, WKRL, LWKRL, INFORL )
!         B = 0.0d0
!         DO i = 1, RANKRL
!           B( JPVTRL( i ) ) = BRL( i )
!         END DO
          RETURN
          END SUBROUTINE ModRegLinRL


!         --------------------------------------------------------------
!         [GetVHatLP] Obtiene el vector de la matrix hat de la estimación
!                     polinómico local correspondiente a la última regresión lineal
!         --------------------------------------------------------------
          SUBROUTINE GetVHatLP(VHAT)
          REAL(8)  VHAT(NRL)
!         Variables locales
          INTEGER i
!         --------------------------------------------------------------
          VHAT = 0.0d0
          VHAT(1) = 1.0d0
!         Obtener Z * VHAT
          IF( RANKRL.LT.NINDRL ) THEN
              CALL DORMRZ( 'Left', 'No transpose', NINDRL, 1, RANKRL, NINDRL-RANKRL,  &
     &                RRL, LDXRL, WKRL( NINDRL+1 ), VHAT, NRL, WKRL( 2*NINDRL+1 ),    &
     &                LWKRL-2*NINDRL, INFORL )
          END IF
!         Obtener (R^-1)^t * Z * VHAT
          CALL DTRSV('Upper', 'Transpose', 'Non-unit', RANKRL, RRL, LDXRL, VHAT, 1)
!         Obtener Q * (R^-1)^t * Z * VHAT
          CALL DORMQR( 'Left', 'No transpose', NRL, 1, NINDRL, RRL, LDXRL,            &
     &        WKRL( 1 ), VHAT, NRL, WKRL( 2*NINDRL+1 ), LWKRL-2*NINDRL, INFORL )
!         Multiplicar por los pesos almacenados en XRL(i, 1)
          DO i = 1, NRL
            VHAT(i) = VHAT(i)*XRL(i,1)
          END DO
          RETURN
          END SUBROUTINE GetVHatLP

!     ------------------------------------------------------------------
      END MODULE linreg_module
!     ------------------------------------------------------------------



!     ------------------------------------------------------------------
!     [RegLin]   Regresión lineal
!     ------------------------------------------------------------------
      SUBROUTINE RegLin(N, NIND, X, LDX, Y, B, INFO)
      USE linreg_module
      IMPLICIT NONE
      INTEGER LDX, NIND, N, INFO, i
      REAL(8)  X(LDX,NIND), Y(N), B(NIND)
!     ------------------------------------------------------------------
!     Establecer valores
      CALL ModRegLinInit(N, NIND)
      NRL = N
      XRL(1:NRL, 1:NINDRL) = X(1:NRL, 1:NINDRL)
      YRL(1:NRL) = Y(1:NRL)
!     Realizar regresión
      CALL ModRegLinRL
!     Resultados
      INFO = INFORL
      IF( INFO.GT.0 ) THEN
          CALL Error(INFO,'RegLin')
      END IF
      B = 0.0d0
      DO i = 1, RANKRL
       B( JPVTRL( i ) ) = BRL( i )
      END DO
      CALL ModRegLinExit
      RETURN
      END SUBROUTINE RegLin



