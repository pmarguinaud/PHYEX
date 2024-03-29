!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ##########################
      MODULE MODD_WET_DEP_PARAM
!     ##########################
!
!!****  *MODD_WET_DEP_PARAM* - declaration of some microphysical factors
!!                              extensively used in the warm and cold schemes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare some precomputed
!     microphysical paramters directly used in routine RAIN_ICE.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_RAIN_ICE_PARAM)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95
!!       J.-P. Pinty   29/11/02 add ICE4
!!       M. Mokhtari 01-02-2011 Adaptation of MODD_RAIN_ICE_PARAM under MADD_WET_DEP_PARAM for Aladin
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL,DIMENSION(2),SAVE :: XFSEDC                 ! Constants for sedimentation fluxes of C
REAL,SAVE :: XFSEDR,XEXSEDR,                   & ! Constants for sedimentation
             XFSEDI,XEXCSEDI,XEXRSEDI,         & ! fluxes of R, I, S and G
             XFSEDS,XEXSEDS,                   &
             XFSEDG,XEXSEDG
!
REAL,SAVE :: XNU10,XALPHA1,XBETA1,             & ! Constants for heterogeneous
             XNU20,XALPHA2,XBETA2,             & ! ice nucleation : HEN
             XMNU0                               ! mass of nucleated ice crystal
!
REAL,SAVE :: XALPHA3,XBETA3,                   & ! Constants for homogeneous
             XHON                                ! ice nucleation : HON
!
REAL,SAVE :: XSCFAC,                           & ! Constants for raindrop
             X0EVAR,X1EVAR,XEX0EVAR,XEX1EVAR,  & ! evaporation: EVA and for
             X0DEPI,X2DEPI,                    & ! deposition : DEP on I,
             X0DEPS,X1DEPS,XEX0DEPS,XEX1DEPS,  & !                  on S and
             X0DEPG,X1DEPG,XEX0DEPG,XEX1DEPG     !                  on G
!
REAL,SAVE :: XTIMAUTI,XTEXAUTI,                & ! Constants for pristine ice
             XCRIAUTI                            ! autoconversion : AUT
!
REAL,SAVE :: XCOLIS,XCOLEXIS,                  & ! Constants for snow
             XFIAGGS,                          & ! aggregation : AGG
             XEXIAGGS
!
REAL,SAVE :: XTIMAUTC,                         & ! Constants for cloud droplet
             XCRIAUTC                            ! autoconversion : AUT
!
REAL,SAVE :: XFCACCR,                          & ! Constants for cloud droplet
             XEXCACCR                            ! accretion on raindrops : ACC
!
REAL,SAVE :: XDCSLIM,XCOLCS,                   & ! Constants for the riming of
             XEXCRIMSS,XCRIMSS,                & ! the aggregates : RIM
             XEXCRIMSG,XCRIMSG,                & !
             XEXSRIMCG,XSRIMCG,                & !
             XGAMINC_BOUND_MIN,                & ! Min val. of Lbda_s for RIM
             XGAMINC_BOUND_MAX,                & ! Max val. of Lbda_s for RIM
             XRIMINTP1,XRIMINTP2                 ! Csts for lin. interpol. of
                                                 ! the tab. incomplete Gamma law
INTEGER,SAVE :: NGAMINC                          ! Number of tab. Lbda_s
REAL, DIMENSION(:), SAVE, ALLOCATABLE          &
                       :: XGAMINC_RIM1,        & ! Tab. incomplete Gamma funct.
                          XGAMINC_RIM2           ! for XDS+2 and for XBS
!
REAL,SAVE :: XFRACCSS,                         & ! Constants for the accretion
             XLBRACCS1,XLBRACCS2,XLBRACCS3,    & ! raindrops onto the aggregates
             XFSACCRG,                         & ! ACC (processes RACCSS and
             XLBSACCR1,XLBSACCR2,XLBSACCR3,    & !                SACCRG)
             XACCLBDAS_MIN,                    & ! Min val. of Lbda_s for ACC
             XACCLBDAS_MAX,                    & ! Max val. of Lbda_s for ACC
             XACCLBDAR_MIN,                    & ! Min val. of Lbda_r for ACC
             XACCLBDAR_MAX,                    & ! Max val. of Lbda_r for ACC
             XACCINTP1S,XACCINTP2S,            & ! Csts for bilin. interpol. of
             XACCINTP1R,XACCINTP2R               !   Lbda_s and Lbda_r in the
                                                 ! XKER_RACCSS and XKER_SACCRG
                                                 !            tables
INTEGER,SAVE :: NACCLBDAS,                     & ! Number of Lbda_s values and
                NACCLBDAR                        !   of Lbda_r values in the
                                                 ! XKER_RACCSS and XKER_SACCRG
                                                 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
                         :: XKER_RACCSS,       & ! Normalized kernel for RACCSS
                            XKER_RACCS,        & ! Normalized kernel for RACCS
                            XKER_SACCRG          ! Normalized kernel for SACCRG
REAL,SAVE :: XFSCVMG                             ! Melting-conversion factor of
                                                 ! the aggregates
!
REAL,SAVE :: XCOLIR,                           & ! Constants for rain contact
             XEXRCFRI,XRCFRI,                  & ! freezing : CFR
             XEXICFRR,XICFRR                     !
!
REAL,SAVE :: XFCDRYG,                          & ! Constants for the dry growth
             XCOLIG,XCOLEXIG,XFIDRYG,          & ! of the graupeln : DRY
             XCOLSG,XCOLEXSG,XFSDRYG,          & !   processes RCDRYG
             XLBSDRYG1,XLBSDRYG2,XLBSDRYG3,    & !             RIDRYG
             XFRDRYG,                          & !             RSDRYG
             XLBRDRYG1,XLBRDRYG2,XLBRDRYG3,    & !             RRDRYG
             XDRYLBDAR_MIN,                    & ! Min val. of Lbda_r for DRY
             XDRYLBDAR_MAX,                    & ! Max val. of Lbda_r for DRY
             XDRYLBDAS_MIN,                    & ! Min val. of Lbda_s for DRY
             XDRYLBDAS_MAX,                    & ! Max val. of Lbda_s for DRY
             XDRYLBDAG_MIN,                    & ! Min val. of Lbda_g for DRY
             XDRYLBDAG_MAX,                    & ! Max val. of Lbda_g for DRY
             XDRYINTP1R,XDRYINTP2R,            & ! Csts for bilin. interpol. of
             XDRYINTP1S,XDRYINTP2S,            & ! Lbda_r, Lbda_s and Lbda_g in
             XDRYINTP1G,XDRYINTP2G               ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
INTEGER,SAVE :: NDRYLBDAR,                     & ! Number of Lbda_r,
                NDRYLBDAS,                     & !        of Lbda_s and
                NDRYLBDAG                        !        of Lbda_g values in
                                                 ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
                         :: XKER_SDRYG,        & ! Normalized kernel for SDRYG
                            XKER_RDRYG           ! Normalized kernel for RDRYG
!
! addition of Hail category
!
REAL,SAVE :: XFSEDH,XEXSEDH                      ! Constants for sedimentation
!
!
REAL,SAVE :: X0DEPH,X1DEPH,XEX0DEPH,XEX1DEPH     ! Constants for deposition
!
REAL,SAVE :: XFWETH,XFSWETH,                   & ! Constants for the wet growth
             XLBSWETH1,XLBSWETH2,XLBSWETH3,    & ! of the hailstones : WET
             XFGWETH,                          & !   processes RSWETH
             XLBGWETH1,XLBGWETH2,XLBGWETH3,    & !             RGWETH
             XWETLBDAS_MIN,                    & ! Min val. of Lbda_s for WET
             XWETLBDAS_MAX,                    & ! Max val. of Lbda_s for WET
             XWETLBDAG_MIN,                    & ! Min val. of Lbda_g for WET
             XWETLBDAG_MAX,                    & ! Max val. of Lbda_g for WET
             XWETLBDAH_MIN,                    & ! Min val. of Lbda_h for WET
             XWETLBDAH_MAX,                    & ! Max val. of Lbda_h for WET
             XWETINTP1S,XWETINTP2S,            & ! Csts for bilin. interpol. of
             XWETINTP1G,XWETINTP2G,            & ! Lbda_r, Lbda_s and Lbda_g in
             XWETINTP1H,XWETINTP2H               ! the XKER_SWETH and XKER_GWETH
                                                 !            tables
INTEGER,SAVE :: NWETLBDAS,                     & ! Number of Lbda_s,
                NWETLBDAG,                     & !        of Lbda_g and
                NWETLBDAH                        !        of Lbda_h values in
                                                 ! the XKER_SWETH and XKER_GWETH
                                                 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
                         :: XKER_SWETH,        & ! Normalized kernel for SWETH
                            XKER_GWETH           ! Normalized kernel for GWETH
!
END MODULE MODD_WET_DEP_PARAM
