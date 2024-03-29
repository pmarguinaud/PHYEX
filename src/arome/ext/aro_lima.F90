!     ######spl
      SUBROUTINE  ARO_LIMA(KKA,KKU,KKL,KLON,KLEV,  KRR, KSV, KTCOUNT, KSPLITR, KSPLITG, &
                                  PTSTEP, PDZZ, PRHODJ, PRHODREF, PEXNREF,&
                                  PPABSM, PW_NU, PTHT, PRT, PSVT, &
                                  PTHS, PRS, PSVS, PEVAP,  &
                                  PINPRR,PINPRS,                 &
                                  PINPRG,PINPRH,PFPR,     &
                                  YDDDH, YDLDDH, YDMDDH    )

      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##########################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation in LIMA
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    Vie et al., 2015 GMD
!!
!!    AUTHOR
!!    ------
!!    B. Vie 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/09/13
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CONF
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_RAIN_ICE_DESCR
!
USE MODD_PARAM_ICE
USE MODD_PARAM_LIMA
USE MODD_NSV
!
USE MODD_BUDGET
USE MODI_BUDGET_DDH
!
USE MODI_LIMA
!
USE MODI_LIMA_WARM
USE MODI_LIMA_COLD
USE MODI_LIMA_MIXED
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
INTEGER,                  INTENT(IN)   :: KKA  !near ground array index
INTEGER,                  INTENT(IN)   :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN)   :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSV      ! Number of LIMA variables
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
INTEGER,                  INTENT(IN)   :: KSPLITG  ! Number of small time step
                                       ! integrations for  graupel/hail sedimendation
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PDZZ     ! Height (z)
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PPABSM  ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PW_NU   ! w for CCN activation
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT):: PRT   ! Moist variables at time t
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(INOUT):: PSVT   ! LIMA variables at time t
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(INOUT) :: PSVS   ! LIMA variable sources
REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT) :: PEVAP ! Rain evap profile
!
!

REAL, DIMENSION(KLON,1), INTENT(INOUT)     :: PINPRR! Rain instant precip
REAL, DIMENSION(KLON,1), INTENT(INOUT)     :: PINPRS! Snow instant precip
REAL, DIMENSION(KLON,1), INTENT(INOUT)     :: PINPRG! Graupel instant precip
REAL, DIMENSION(KLON,1), INTENT(INOUT)     :: PINPRH! Hail instant precip
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PFPR ! upper-air precip
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
!
!*       0.2   Declarations of local variables :

CHARACTER(LEN=4)            :: HCLOUD   ! kind of cloud
                                                   ! paramerization

INTEGER                   :: KMI      ! Model index

!
INTEGER :: JRR, JL           ! Loop index for the moist and scalar variables
!
!
!
REAL, DIMENSION(KLON,1,KLEV):: ZT,ZLV,ZLS,ZCPH
REAL, DIMENSION(KLON,1,KLEV):: ZCOR,ZDUM3DC,ZDUM3DR,ZDUM3DS,ZDUM3DG,ZDUM3DH
REAL, DIMENSION(KLON,1,KLEV):: &
   & ZRAINFR, ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC
REAL, DIMENSION(KLON,1):: ZINPRC    ! surf cloud sedimentation
                                    ! for the correction of negative rv
REAL, DIMENSION(KLON,1):: ZINPRI    ! surf cloud ice sedimentation
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR

LOGICAL :: LL_RRR_BUDGET
!
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_LIMA',0,ZHOOK_HANDLE)

HCLOUD='LIMA'
KMI=1
ZINPRC=0.
ZDUM3DC=0.
ZDUM3DR=0.
ZDUM3DS=0.
ZDUM3DG=0.
ZDUM3DH=0.
PINPRH=0.


!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!
!  complete the vertical boundaries
!
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
ZT(:,:,:)= PTHT(:,:,:)*PEXNREF(:,:,:)
ZLV(:,:,:)=XLVTT +(XCPV-XCL) *(ZT(:,:,:)-XTT)
ZLS(:,:,:)=XLSTT +(XCPV-XCI) *(ZT(:,:,:)-XTT)
ZCPH(:,:,:)=XCPD +XCPV*2.*PTSTEP*PRS(:,:,:,1)
!

!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for 1-moment precipitating species (Rood 87)
!
DO JRR = 3,KRR
   SELECT CASE (JRR)
   CASE(5,6,7) ! snow, graupel and hail
        WHERE (PRS(:,:,:,JRR) < 1.E-15 )
           PRS(:,:,:,JRR) = 0.
        END WHERE
   END SELECT
END DO

!
!*       3.2    Correct negative values
!
! Correction where rc<0
     IF (LWARM_LIMA) THEN
        WHERE (PRS(:,:,:,2) < 1.E-15 .OR. PSVS(:,:,:,NSV_LIMA_NC) < 1.E-15)
           PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
           PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
                ZCPH(:,:,:) / PEXNREF(:,:,:)
           PRS(:,:,:,2)  = 0.0
           PSVS(:,:,:,NSV_LIMA_NC) = 0.0
        END WHERE
     END IF
! Correction where rr<0
     IF (LWARM_LIMA .AND. LRAIN_LIMA) THEN
        WHERE (PRS(:,:,:,3) < 1.E-15 .OR. PSVS(:,:,:,NSV_LIMA_NR) < 1.E-15)
           PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,3)
           PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,3) * ZLV(:,:,:) /  &
                ZCPH(:,:,:) / PEXNREF(:,:,:)
           PRS(:,:,:,3)  = 0.0
           PSVS(:,:,:,NSV_LIMA_NR) = 0.0
        END WHERE
     END IF
! Correction of IFN concentrations where ri<0 or Ni<0
!     IF (LCOLD_LIMA) THEN
!        DO JMOD = 1, NMOD_IFN 
!           WHERE (PRS(:,:,:,4) < 0. .OR. PSVS(:,:,:,NSV_LIMA_NI) < 0.) ! ri or Ni < 0.
!              PSVS(:,:,:,NSV_LIMA_IFN_FREE+JMOD-1) =               &
!                   PSVS(:,:,:,NSV_LIMA_IFN_FREE+JMOD-1) + &
!                   PSVS(:,:,:,NSV_LIMA_IFN_NUCL+JMOD-1)     ! N_IF =N_IF+N_IN
!              PSVS(:,:,:,NSV_LIMA_IFN_NUCL+JMOD-1) = 0.0             ! N_IN =0.
!           END WHERE
!        ENDDO
!     END IF
! Correction where ri<0
     IF (LCOLD_LIMA) THEN
        WHERE (PRS(:,:,:,4) < 1.E-15 .OR. PSVS(:,:,:,NSV_LIMA_NI) < 1.E-15)
           PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
           PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLS(:,:,:) /  &
                ZCPH(:,:,:) / PEXNREF(:,:,:)
           PRS(:,:,:,4)  = 0.0
           PSVS(:,:,:,NSV_LIMA_NI) = 0.0
        END WHERE
     END IF
!
     PSVS(:,:,:,:) = MAX( 0.0,PSVS(:,:,:,:) )
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------

LL_RRR_BUDGET = (LBUDGET_RV).OR.(LBUDGET_RC).OR.(LBUDGET_RR).OR.(LBUDGET_RI) &
             & .OR.(LBUDGET_RS).OR.(LBUDGET_RG).OR.(LBUDGET_RH)       

IF (LL_RRR_BUDGET) THEN      
  DO JRR=1,KRR
     CALL BUDGET_DDH (PRS(:,:,:,JRR) * PRHODJ(:,:,:), JRR+5,'NEGA_BU_RRR',YDDDH,YDLDDH, YDMDDH)
  END DO 
END IF
IF (LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:,:)  * PRHODJ(:,:,:),4,'NEGA_BU_RTH',YDDDH, YDLDDH, YDMDDH)
IF (LBUDGET_SV) THEN
   CALL BUDGET_DDH (PSVS(:,:,:,NSV_LIMA_NC)*PRHODJ(:,:,:),12+NSV_LIMA_NC,'NEGA_BU_RSV',YDDDH, YDLDDH, YDMDDH)
   CALL BUDGET_DDH (PSVS(:,:,:,NSV_LIMA_NR)*PRHODJ(:,:,:),12+NSV_LIMA_NR,'NEGA_BU_RSV',YDDDH, YDLDDH, YDMDDH)
   CALL BUDGET_DDH (PSVS(:,:,:,NSV_LIMA_NI)*PRHODJ(:,:,:),12+NSV_LIMA_NI,'NEGA_BU_RSV',YDDDH, YDLDDH, YDMDDH)
   IF (NMOD_CCN.GE.1) THEN
      DO JL=1, NMOD_CCN
         CALL BUDGET_DDH ( PSVS(:,:,:,NSV_LIMA_CCN_FREE+JL-1)* &
              PRHODJ(:,:,:),12+NSV_LIMA_CCN_FREE+JL-1,'NEGA_BU_RSV',YDDDH,YDLDDH, YDMDDH) 
      END DO
   END IF
   IF (NMOD_IFN.GE.1) THEN
      DO JL=1, NMOD_IFN
         CALL BUDGET_DDH ( PSVS(:,:,:,NSV_LIMA_IFN_FREE+JL-1)* &
              PRHODJ(:,:,:),12+NSV_LIMA_IFN_FREE+JL-1,'NEGA_BU_RSV',YDDDH,YDLDDH, YDMDDH) 
      END DO
   END IF
END IF


!
!
!-------------------------------------------------------------------------------
!

!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!*                   Compute the explicit microphysical sources
!
!
!
IF (LPTSPLIT) THEN
   CALL LIMA (PTSTEP=2*PTSTEP, HFMFILE='DUMMY', OCLOSE_OUT=.FALSE.,                    &
              PRHODREF=PRHODREF, PEXNREF=PEXNREF, PZZ=PDZZ,                         &
              PRHODJ=PRHODJ, PPABST=PPABSM,                                 &
              NCCN=NMOD_CCN, NIFN=NMOD_IFN, NIMM=NMOD_IMM,                   &
              PTHM=PTHT, PTHT=PTHT, PRT=PRT, PSVT=PSVT, PW_NU=PW_NU,                  &
              PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                &
              PINPRC=ZINPRC, PINPRR=PINPRR, PINPRI=ZINPRI, PINPRS=PINPRS, PINPRG=PINPRG, PINPRH=PINPRH, &
              PEVAP3D=PEVAP, KSPLITR=KSPLITR, KSPLITG=KSPLITG, YDDDH=YDDDH, YDLDDH=YDLDDH, YDMDDH=YDMDDH                                   )
ELSE
   IF (LWARM_LIMA) CALL LIMA_WARM(OACTIT=LACTIT_LIMA, OSEDC=LSEDC_LIMA, ORAIN=LRAIN_LIMA, KSPLITR=KSPLITR, PTSTEP=2*PTSTEP, KMI=KMI, &
                             HFMFILE='DUMMY', HLUOUT='DUMMY', OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ,     &
                             PRHODREF=PRHODREF, PEXNREF=PEXNREF, PW_NU=PW_NU, PPABSM=PPABSM, PPABST=PPABSM,         &
                             PTHM=PTHT, PRCM=PRT(:,:,:,2),                  &
                             PTHT=PTHT, PRT=PRT, PSVT=PSVT,                                   &
                             PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                   &
                             PINPRC=ZINPRC,PINPRR=PINPRR, PINPRR3D=ZDUM3DR, PEVAP3D=PEVAP,YDDDH=YDDDH, YDLDDH=YDLDDH, YDMDDH=YDMDDH         )
   !
   IF (LCOLD_LIMA) CALL LIMA_COLD(OSEDI=LSEDI_LIMA, OHHONI=LHHONI_LIMA, KSPLITG=KSPLITG, PTSTEP=2*PTSTEP, KMI=KMI,               &
                             HFMFILE='DUMMY', HLUOUT='DUMMY', OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ,     &
                             PRHODREF=PRHODREF, PEXNREF=PEXNREF, PPABST=PPABSM, PW_NU=PW_NU,                 &
                             PTHM=PTHT, PPABSM=PPABSM,                                      &
                             PTHT=PTHT, PRT=PRT, PSVT=PSVT,                                   &
                             PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                   &
                             PINPRS=PINPRS, PINPRG=PINPRG, PINPRH=PINPRH, YDDDH=YDDDH, YDLDDH=YDLDDH, YDMDDH=YDMDDH)
   !
   IF (LWARM_LIMA .AND. LCOLD_LIMA) CALL LIMA_MIXED(OSEDI=LSEDI_LIMA, OHHONI=LHHONI_LIMA, KSPLITG=KSPLITG, PTSTEP=2*PTSTEP, KMI=KMI,     &
                                          HFMFILE='DUMMY', HLUOUT='DUMMY', OCLOSE_OUT=.FALSE., KRR=KRR, PZZ=PDZZ, PRHODJ=PRHODJ, &
                                          PRHODREF=PRHODREF, PEXNREF=PEXNREF, PPABST=PPABSM, PW_NU=PW_NU,             &
                                          PTHM=PTHT, PPABSM=PPABSM,                                  &
                                          PTHT=PTHT, PRT=PRT, PSVT=PSVT,                               &
                                          PTHS=PTHS, PRS=PRS, PSVS=PSVS,YDDDH=YDDDH, YDLDDH=YDLDDH, YDMDDH=YDMDDH                                )
   
ENDIF
!add ZINPRC in PINPRR
PINPRR=PINPRR+ZINPRC
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_LIMA',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_LIMA
