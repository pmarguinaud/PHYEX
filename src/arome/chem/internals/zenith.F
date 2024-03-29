!     ######spl
      SUBROUTINE zenith(lat,long,idate,ut,azim,zen)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      REAL lat,long
      REAL ut
      INTEGER idate

! output:
      REAL azim, zen

! local:
      REAL lbut,lzut
      REAL rlt
      REAL d, tz, rdecl, eqr, eqh, zpt
      REAL csz, zr, caz, raz
      REAL sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz

      INTEGER iiyear, imth, iday, ijd
      INTEGER imn(12)

      INTEGER i

! program constants:

      REAL pi, dr
      PARAMETER(pi=3.1415926535898)
      PARAMETER (dr=pi/180.D0)
!_______________________________________________________________________

      DATA imn/31,28,31,30,31,30,31,31,30,31,30,31/
!_______________________________________________________________________

! convert to radians

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ZENITH',0,ZHOOK_HANDLE)
      rlt = lat*dr

! parse date

      iiyear = idate/10000
      imth = (idate - iiyear*10000)/100
      iday = idate - iiyear*10000 - imth*100

! identify and correct leap years

      IF (MOD(iiyear,4) .EQ. 0) THEN
         imn(2) = 29
      ELSE
         imn(2) = 28
      ENDIF

! compute current (Julian) day of year IJD = 1 to 365

      ijd = 0
      DO 30, i = 1, imth - 1
         ijd = ijd + imn(i)
   30 CONTINUE
      ijd = ijd + iday

! calculate decimal Julian day from start of year:

      d = FLOAT(ijd-1) + ut/24.

! Equation 3.8 for "day-angle"

      tz = 2.*pi*d/365.

! Calculate sine and cosine from addition theoremes for
! better performance;  the computation of sin2tz,
! sin3tz, cos2tz and cos3tz is about 5-6 times faster
! than the evaluation of the intrinsic functions
!
! It is SIN(x+y) = SIN(x)*COS(y)+COS(x)*SIN(y)
! and   COS(x+y) = COS(x)*COS(y)-SIN(x)*SIN(y)
!
! sintz  = SIN(tz)      costz  = COS(tz)
! sin2tz = SIN(2.*tz)   cos2tz = SIN(2.*tz)
! sin3tz = SIN(3.*tz)   cos3tz = COS(3.*tz)
!
      sintz = SIN(tz)
      costz = COS(tz)
      sin2tz = 2.*sintz*costz
      cos2tz = costz*costz-sintz*sintz
      sin3tz = sintz*cos2tz + costz*sin2tz
      cos3tz = costz*cos2tz - sintz*sin2tz

! Equation 3.7 for declination in radians

      rdecl = 0.006918 - 0.399912*costz  + 0.070257*sintz               &
     &                 - 0.006758*cos2tz + 0.000907*sin2tz              &
     &                 - 0.002697*cos3tz + 0.001480*sin3tz

! Equation 3.11 for Equation of time  in radians

      eqr   = 0.000075 + 0.001868*costz  - 0.032077*sintz               &
     &                 - 0.014615*cos2tz - 0.040849*sin2tz

! convert equation of time to hours:

      eqh = eqr*24./(2.*pi)

! calculate local hour angle (hours):

      lbut = 12. - eqh - long*24./360

! convert to angle from UT

      lzut = 15.*(ut - lbut)
      zpt = lzut*dr

! Equation 2.4 for cosine of zenith angle

      csz = SIN(rlt)*SIN(rdecl) + COS(rlt)*COS(rdecl)*COS(zpt)
      zr = ACOS(csz)
      zen = zr/dr

!   calc local solar azimuth

      caz = (SIN(rdecl) - SIN(rlt)*COS(zr))/(COS(rlt)*SIN(zr))
      raz = ACOS(caz)
      azim = raz/dr
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('ZENITH',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE zenith
