!     ######spl
     MODULE MODI_AER_VELGRAV
!!   ########################
!!
INTERFACE
!!
SUBROUTINE AER_VELGRAV(PRG, PABST, KMODE, PMU, &
                       PVGG, PDPG,PTEMP, PCOR, PDENSITY_AER)
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(IN) :: PRG
REAL, DIMENSION(:),   INTENT(IN) :: PTEMP
REAL, DIMENSION(:),   INTENT(IN) :: PABST   
INTEGER,              INTENT(IN) :: KMODE
REAL,                 INTENT(IN) :: PDENSITY_AER
REAL, DIMENSION(:),   INTENT(OUT) :: PMU
REAL, DIMENSION(:,:), INTENT(OUT) :: PVGG, PDPG
REAL, DIMENSION(:,:), INTENT(OUT) :: PCOR
END SUBROUTINE AER_VELGRAV
!!
END INTERFACE
!!
END MODULE MODI_AER_VELGRAV
