!     ######spl
     MODULE MODI_CH_AER_INIT
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_INIT(PCHEM, PAERO, PRHODREF)
IMPLICIT NONE
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)   :: PCHEM, PAERO
REAL,       DIMENSION(:,:,:),    INTENT(IN)      :: PRHODREF
END SUBROUTINE CH_AER_INIT
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_INIT
