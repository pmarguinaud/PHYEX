!     ######spl
     MODULE MODI_SALT_FILTER
!!   ##############################
!!
INTERFACE
!
SUBROUTINE SALT_FILTER(PSV, PRHODREF)

IMPLICIT NONE

REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF

END SUBROUTINE SALT_FILTER
!!
END INTERFACE
!!
END MODULE MODI_SALT_FILTER
