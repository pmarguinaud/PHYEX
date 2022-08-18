!     ##################
      MODULE DDH_MIX
!     ##################


!!    PURPOSE
!!    -------
!       module for type structure in  DDH

!!**  IMPLICIT ARGUMENTS
!!    ------------------

!!    AUTHOR
!!    ------
!!    O.Riviere   *Meteo France*
!!
!!    Modifications :
!!    --------------
!!    F.Voitus     18/09/17 : New DDH flexible structure for OpenMP
  
!-------------------------------------------------------------------------------


USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK,DR_HOOK

IMPLICIT NONE

SAVE


!*       1.   TYPE DEFINITION


TYPE DDHFLEX
  CHARACTER(LEN=11)::CNAME !name of field
  CHARACTER(LEN=1)::CFLUX !'F' if flux 'V' if variable 'T' if tendency
  CHARACTER(LEN=3)::CMOD  ! 'ARP','ARO': name of model
  LOGICAL:: LKDDH !TRUE if to be stored into DDH

! rfield has to be a pointer because allocatable not allowed in structure type
  REAL(KIND=JPRB),DIMENSION(:,:),POINTER:: RFIELD=>NULL()  ! value of retrieved field
  INTEGER(KIND=JPIM):: NFIELDIND! position of flux in ddh array
END TYPE DDHFLEX


! structure for vertical profile
TYPE TYP_FIELD3D
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL() ! values (domains x vertical)
  CHARACTER(LEN=16)         :: CNAME             ! name
  CHARACTER(LEN=1)          :: CTYPE             ! indicate wether it's a flux, variable or tendency
END TYPE TYP_FIELD3D

TYPE TYP_FIELD2D
  REAL(KIND=JPRB), POINTER :: RVAL(:)=>NULL() ! values (domains)
  CHARACTER(LEN=16)        :: CNAME           ! name
  CHARACTER(LEN=1)         :: CTYPE           ! indicate wether it's a flux, variable or tendency
END TYPE TYP_FIELD2D

TYPE TYP_BUVAR
  REAL(KIND=JPRB), POINTER  :: RVAL(:,:)=>NULL() ! values (domains)
  CHARACTER(LEN=2)          :: CNAME             ! name
END TYPE TYP_BUVAR


! structure containing fields, weights for horizontal averaging, etc.
TYPE TYP_DDH
  INTEGER(KIND=JPIM)          :: NLEV                  ! vertical dimension
  INTEGER(KIND=JPIM)          :: NPROMA                ! horizontal dimension                                                    
  INTEGER(KIND=JPIM)          :: KST                   ! first point
  INTEGER(KIND=JPIM)          :: KEND                  ! last point

  ! 3D field info
  INTEGER(KIND=JPIM)          :: NFIELDS3D             ! number of fields
  INTEGER(KIND=JPIM)          :: NFIELDS3D_OFFSET      ! offset value for the fields when 
                                                       ! putting in the global arrays
  INTEGER(KIND=JPIM)          :: NFIELDS3D_AUTO        ! number of fields that should be 
                                                       ! automatically allocated 
                                                       ! (maximum of cpg and cpglag)
  TYPE(TYP_FIELD3D), POINTER  :: YFIELD3D(:)=>NULL()   ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVAL3D(:,:,:)=>NULL() ! auxiliary array for fast 
                                                       ! (automatic) allocation
                                                       ! (vertical x fields x domains)
  ! 2D field info
  INTEGER(KIND=JPIM)          :: NFIELDS2D             ! number of fields
  INTEGER(KIND=JPIM)          :: NFIELDS2D_OFFSET      ! offset value for the fields when 
                                                       ! putting in the global arrays
  INTEGER(KIND=JPIM)          :: NFIELDS2D_AUTO        ! number of fields that should be 
                                                       ! automatically allocated 
                                                       ! (maximum of cpg and cpglag)
  TYPE(TYP_FIELD2D), POINTER  :: YFIELD2D(:)=>NULL()   ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVAL2D(:,:)=>NULL()   ! auxiliary array for fast 
                                                       ! (automatic) allocation
                                                       ! (vertical x fields x domains)
  ! horizontal info
  REAL(KIND=JPRB), POINTER    :: WEIGHT(:)=>NULL()     ! weights inside one NPROMA block
  INTEGER(KIND=JPIM), POINTER :: NDDHI(:)=>NULL()      ! cfr. KDDHI in cpg_dia


  TYPE(TYP_BUVAR), POINTER    :: YVARMULT(:)=>NULL()    ! array of fields
  REAL(KIND=JPRB), POINTER    :: RVARSM(:,:,:,:)=>NULL()! auxiliary array for fast 
                                                        ! (automatic) allocation

END TYPE TYP_DDH



!*       2.   DECLARATIONS

REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE,TARGET::RDDH_FIELD
! ddh_field is the target of rfield

TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: RDDH_DESCR ! array of several structure type
INTEGER(KIND=JPIM):: NFLEVGDDH=0   ! number of processed field
INTEGER(KIND=JPIM):: NPROMADDH=0   ! number of processed field
INTEGER(KIND=JPIM):: NFIELDDDH=0   ! number of processed field
INTEGER(KIND=JPIM):: NTOTFIELD=0 !nb of fields in structure
INTEGER(KIND=JPIM):: NTOTVAR=0 ! nb of var. in structure
INTEGER(KIND=JPIM):: NTOTFSL=0 ! nb of SL fields. in structure
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDH ! mask for operations in ddh operators
CHARACTER(LEN=13)::CNAMEINI ! name of first field
LOGICAL:: LALLOC=.FALSE.  !if TRUE allocation of arrrays completed


REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE,TARGET::RDDHSURF_FIELD
! surface ddh_field is the target of rfield
TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: RDDHSURF_DESCR ! array of several structure type
INTEGER(KIND=JPIM):: NSURFDDH=0   ! number of processed surface field
INTEGER(KIND=JPIM):: NTOTSURF=0 !nb of surface fields in structure
INTEGER(KIND=JPIM):: NTOTSVAR=0 ! nb of surface var. in structure
INTEGER(KIND=JPIM):: NTOTSVFS=0 ! nb of surface free style var. in structure
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDHVAS ! mask for operations in ddh operators
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDHVFS ! mask for operations in ddh operators
LOGICAL,DIMENSION(:),ALLOCATABLE:: LMASKDDHSURF ! mask for operations in ddh operators
CHARACTER(LEN=13)::CNAMESURFINI ! name of first surface field
LOGICAL:: LALLOCSURF=.FALSE.  !if TRUE allocation of arrrays completed


CONTAINS

!-----------------------------------------------
!*       3.   SUBROUTINES
SUBROUTINE ADD_FIELD_3D(YDLDDH,PMAT,CDNAME,CDFLUX,CDMOD,LDINST,LDDH,LDSLFD)
USE YOMLDDH , ONLY : TLDDH
TYPE(TLDDH)                   ,INTENT(IN)   :: YDLDDH
REAL(KIND=JPRB),DIMENSION(:,:),INTENT(IN)   :: PMAT
CHARACTER(LEN=*)              ,INTENT(IN)   :: CDNAME
CHARACTER(LEN=1)              ,INTENT(IN)   :: CDFLUX
CHARACTER(LEN=3)              ,INTENT(IN)   :: CDMOD
LOGICAL                       ,INTENT(IN)   :: LDINST,LDDH
LOGICAL              ,OPTIONAL,INTENT(IN)   :: LDSLFD

REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE::ZDDH_FIELD_BKP
INTEGER(KIND=JPIM)::IDIM1,IDIM2,IDIM3,ITEST
TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: YLDDH_DESCR_BKP
LOGICAL::LLSLFD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE ADD_FIELD_3D

SUBROUTINE ADD_FIELD_2D(YDLDDH,PMAT,CDNAME,CDFLUX,CDMOD,LDINST,LDDH)
USE YOMLDDH , ONLY : TLDDH
TYPE(TLDDH)                 ,INTENT(IN)   :: YDLDDH
REAL(KIND=JPRB),DIMENSION(:),INTENT(IN)   :: PMAT
CHARACTER(LEN=*)            ,INTENT(IN)   :: CDNAME
CHARACTER(LEN=1)            ,INTENT(IN)   :: CDFLUX
CHARACTER(LEN=3)            ,INTENT(IN)   :: CDMOD
LOGICAL                     ,INTENT(IN)   :: LDINST,LDDH
REAL(KIND=JPRB),DIMENSION(:,:,:),ALLOCATABLE::ZSDDH_FIELD_BKP
INTEGER(KIND=JPIM)::IDIM1,IDIM2,IDIM3,ITEST
TYPE(DDHFLEX),ALLOCATABLE,DIMENSION(:):: YSLDDH_DESCR_BKP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE ADD_FIELD_2D

SUBROUTINE RESET_DDHFLEX

REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE RESET_DDHFLEX

SUBROUTINE MK_MASK
INTEGER(KIND=JPIM)::JFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE MK_MASK

SUBROUTINE MK_MASKSURF
INTEGER(KIND=JPIM)::JFIELD

REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE MK_MASKSURF

SUBROUTINE SUDDHFLEX(KPROMA,KLEV)
INTEGER(KIND=JPIM)::KPROMA,KLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE SUDDHFLEX

SUBROUTINE NEW_ADD_FIELD_3D(YDMDDH,PVAL,CDNAME,YDDDH,CDTYPE)

USE YOMMDDH , ONLY : TMDDH
  IMPLICIT NONE

 TYPE(TMDDH)      , INTENT(IN)    :: YDMDDH
  REAL(KIND=JPRB) , INTENT(IN)    :: PVAL(:,:)          ! numerical values
  CHARACTER(LEN=*), INTENT(IN)    :: CDNAME             ! name of the field
  TYPE(TYP_DDH)   , INTENT(INOUT) :: YDDDH              ! DDH superstructure
  CHARACTER(LEN=1), OPTIONAL, INTENT(IN)   :: CDTYPE   ! type of the field (flux, tendency or variable)
  TYPE(TYP_FIELD3D), ALLOCATABLE   :: YLFIELDS_BKP(:)  ! array of fields, necessary for reallocation
  TYPE(TYP_FIELD3D), POINTER       :: YLFIELD          ! pointer to the added field

END SUBROUTINE NEW_ADD_FIELD_3D


SUBROUTINE NEW_ADD_FIELD_2D(YDMDDH,PVAL,CDNAME,YDDDH,CDTYPE)

USE YOMMDDH , ONLY : TMDDH
  IMPLICIT NONE

 TYPE(TMDDH)      , INTENT(IN)    :: YDMDDH
  REAL(KIND=JPRB) , INTENT(IN)    :: PVAL(:)            ! numerical values
  CHARACTER(LEN=*), INTENT(IN)    :: CDNAME             ! name of the field
  TYPE(TYP_DDH)   , INTENT(INOUT) :: YDDDH              ! DDH superstructure
  CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: CDTYPE    ! type of the field (flux, tendency or variable)
  TYPE(TYP_FIELD2D), ALLOCATABLE   :: YLFIELDS_BKP(:)  ! array of fields, necessary for reallocation
  TYPE(TYP_FIELD2D), POINTER       :: YLFIELD          ! pointer to the added field

  INTEGER(KIND=JPIM) :: JROF, JFLD

END SUBROUTINE NEW_ADD_FIELD_2D

SUBROUTINE SETDDH(YDMDDH,YDDDH,KST,KEND,KPROMA,KLEV,KITER,PXW,KDDHI,PAUX3D,PAUX2D,PAUXSM)

USE YOMMDDH , ONLY : TMDDH
  IMPLICIT NONE

  TYPE(TMDDH)       , INTENT(IN)       :: YDMDDH
  TYPE(TYP_DDH)     , INTENT(INOUT)    :: YDDDH       ! DDH superstructure
  INTEGER(KIND=JPIM), INTENT(IN)       :: KST         ! index of first point in this NPROMA block
  INTEGER(KIND=JPIM), INTENT(IN)       :: KEND        ! index of first point in this NPROMA block
  INTEGER(KIND=JPIM), INTENT(IN)       :: KPROMA      ! horizontal dimension
  INTEGER(KIND=JPIM), INTENT(IN)       :: KLEV        ! vertical dimension
  INTEGER(KIND=JPIM), INTENT(IN)       :: KITER       ! Nb of Corrector iterations

  REAL(KIND=JPRB), TARGET         :: PXW(KPROMA) ! averaging weights
  INTEGER(KIND=JPIM), TARGET      :: KDDHI(:)    ! DDH domain index
  REAL(KIND=JPRB), TARGET         :: PAUX3D(:,:,:)   ! automatically allocated array for 3d fields
  REAL(KIND=JPRB), TARGET         :: PAUX2D(:,:)     ! automatically allocated array for 2d fields
  REAL(KIND=JPRB), TARGET         :: PAUXSM(:,:,:,:)  ! automatically allocated array for variables

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE SETDDH

SUBROUTINE CLEANDDH(YDTDDH,YDDDH,KSTG,KITER)

USE YOMMDDH , ONLY : TMDDH
USE YOMTDDH , ONLY : TTDDH
  IMPLICIT NONE

  TYPE(TTDDH)   , INTENT(INOUT) :: YDTDDH
  TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH      ! DDH superstructure
  INTEGER(KIND=JPIM), INTENT(IN) :: KSTG     ! global index of first point
  INTEGER(KIND=JPIM), INTENT(IN) :: KITER    ! Nbe of corrector iterations

  INTEGER(KIND=JPIM)           :: JFLD

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE CLEANDDH

SUBROUTINE STOREDDH(YDMDDH,YDTDDH,YDDDH)

USE YOMMDDH , ONLY : TMDDH
USE YOMTDDH , ONLY : TTDDH
  IMPLICIT NONE

 TYPE(TMDDH)   , INTENT(INOUT) :: YDMDDH
 TYPE(TTDDH)   , INTENT(IN)    :: YDTDDH
  TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH      ! DDH superstructure
  INTEGER(KIND=JPIM)           :: JFLD, IPLST, IPLEND

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

END SUBROUTINE STOREDDH

END MODULE DDH_MIX

