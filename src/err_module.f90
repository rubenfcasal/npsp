!-----------------------------------------------------------------------
    module err_module
!-----------------------------------------------------------------------
!     Modulo para manejar los codigos de error.
!
!     AVISO: Este modulo es una versión preliminar y será modificado.
!     (para que todas las rutinas controlen los errores)
!
!   Autor: (c) Ruben Fernandez-Casal    
!   Fecha revision: Mar 2002, Jul 2012
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PRIVATE :: errcode, erraction
      DATA erraction / 0 /

!   --------------------------------------------------------------------
    contains
!   --------------------------------------------------------------------

!         --------------------------------------------------------------
!         [SetErrAct] Establecer acción al producirse un error:
!                     iact = -2 -> No parar ni imprimir
!                     iact = -1 -> No parar
!                     iact =  0 -> Acción por defecto
!                     iact =  1 -> Defecto + Imprimir
!                     iact =  2 -> Defecto + Imprimir y parar
!         --------------------------------------------------------------
          SUBROUTINE SetErrAct(iact)
          INTEGER iact
          erraction = iact
          END SUBROUTINE SetErrAct

!         --------------------------------------------------------------
!         [GetErrAct] Devuelve acción en caso de error:
!         --------------------------------------------------------------
          INTEGER FUNCTION GetErrAct()
          GetErrAct = erraction
          END FUNCTION GetErrAct

!         --------------------------------------------------------------
!         [SetErrCode]    Establece el código de error
!         --------------------------------------------------------------
          SUBROUTINE SetErrCode(icode)
          INTEGER icode
          errcode = icode
          END SUBROUTINE SetErrCode

!         --------------------------------------------------------------
!         [GetErrCode]    Devuelve el código de error
!         --------------------------------------------------------------
          INTEGER FUNCTION GetErrCode()
          GetErrCode = errcode
          END FUNCTION GetErrCode

!   --------------------------------------------------------------------
    end module Err_module
!   --------------------------------------------------------------------

!     ------------------------------------------------------------------
!     [Error]  Se llama a esta función si se produce un error
!     ------------------------------------------------------------------
      subroutine error(i, label)
      integer i
      character(*) label
!     ------------------------------------------------------------------
!     Fortran optional arguments produces LTO warnings from gcc9  
!     character(*), optional :: label
!     Print error
!     if (present(label)) WRITE(*,*) label, i
!     STOP 'CRITICAL ERROR'
!     if (present(label)) call intpr(label, -1, i, 1)
!     ------------------------------------------------------------------
      call intpr(label, -1, i, 1)
      call rexit('npsp shared library: CRITICAL ERROR')
      return
      end subroutine error


!     ------------------------------------------------------------------
!     [warning]  Se llama a esta función si se produce un warning
!     ------------------------------------------------------------------
      subroutine warning(label1, label2, i)
      character(*) label1
      character(*), optional :: label2
      integer, optional :: i
!     ------------------------------------------------------------------
!     Escribir error
!     if (present(label2)) WRITE(*,*) label2, i
!     WRITE(*,*) 'WARNING'
!     return
!     ------------------------------------------------------------------
      if (present(label2)) call intpr(label2, -1, i, 1)
      call rwarn(label1)
      return
      end subroutine warning
