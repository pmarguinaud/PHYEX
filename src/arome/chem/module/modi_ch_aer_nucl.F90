!     ######spl
MODULE MODI_CH_AER_NUCL
!!    ################################ 
!!
INTERFACE
  !!
  SUBROUTINE CH_AER_NUCL(ZRH,ZT,ZCONC,ZJ,ZAL,KVECNPT)
  IMPLICIT NONE
  !!
REAL, DIMENSION(:), INTENT(INOUT) :: ZJ,ZAL
REAL, DIMENSION(:), INTENT(IN)    :: ZRH,ZT
REAL, DIMENSION(:), INTENT(INOUT) :: ZCONC
INTEGER, INTENT(IN)               :: KVECNPT
  !!
  !!
  END SUBROUTINE CH_AER_NUCL
  !!
END INTERFACE
!!
END MODULE MODI_CH_AER_NUCL
