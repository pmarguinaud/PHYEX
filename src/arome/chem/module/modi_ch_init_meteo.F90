!     ######spl
      MODULE MODI_CH_INIT_METEO
!!    #########################
!!
!
INTERFACE
SUBROUTINE CH_INIT_METEO(TPM)
USE MODD_CH_M9,      ONLY: METEOTRANSTYPE
TYPE(METEOTRANSTYPE), INTENT(OUT) :: TPM  ! the meteo variables
END SUBROUTINE CH_INIT_METEO
END INTERFACE
END MODULE MODI_CH_INIT_METEO
