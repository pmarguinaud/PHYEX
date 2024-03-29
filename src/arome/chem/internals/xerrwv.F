!     ######spl
      SUBROUTINE XERRWV (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################################
      REAL R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
!
!KS:  changed to adapt to Fortran90
!     CHARACTER*1 MSG(NMES)
      CHARACTER*(*) MSG
!-----------------------------------------------------------------------
! Subroutines XERRWV, XSETF, XSETUN, and the two function routines
! MFLGSV and LUNSAV, as given here, constitute a simplified version of
! the SLATEC error handling package.
! Written by A. C. Hindmarsh and P. N. Brown at LLNL.
! Version of 13 April, 1989.
! This version is in single precision.
!
! All arguments are input arguments.
!
! MSG    = The message (character array).
! NMES   = The length of MSG (number of characters).
! NERR   = The error number (not used).
! LEVEL  = The error level..
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
! NI     = Number of integers (0, 1, or 2) to be printed with message.
! I1,I2  = Integers to be printed, depending on NI.
! NR     = Number of reals (0, 1, or 2) to be printed with message.
! R1,R2  = Reals to be printed, depending on NR.
!
! Note..  this routine is machine-dependent and specialized for use
! in limited context, in the following ways..
! 1. The argument MSG is assumed to be of type CHARACTER, and
!    the message is printed with a format of (1X,80A1).
! 2. The message is assumed to take only one line.
!    Multi-line messages are generated by repeated calls.
! 3. If LEVEL = 2, control passes to the statement   STOP
!    to abort the run.  This statement may be machine-dependent.
! 4. R1 and R2 are assumed to be in single precision and are printed
!    in E21.13 format.
!
! For a different default logical unit number, change the data
! statement in function routine LUNSAV.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWV.. None
! Function routines called by XERRWV.. MFLGSV, LUNSAV
!-----------------------------------------------------------------------
!
      INTEGER I, LUNIT, LUNSAV, MESFLG, MFLGSV
!
! Get message print flag and logical unit number. ----------------------
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('XERRWV',0,ZHOOK_HANDLE)
      MESFLG = MFLGSV (0,.FALSE.)
      LUNIT = LUNSAV (0,.FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
! Write the message. ---------------------------------------------------
!KS:  changed to adapt to Fortran90
!     WRITE (LUNIT,10) (MSG(I),I=1,NMES)
!10   FORMAT(1X,80A1)
      WRITE (LUNIT,'(A)') MSG
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',E21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',E21.13,3X,'R2 =',E21.13)
! Abort the run if LEVEL = 2. ------------------------------------------
 100  CONTINUE
      IF (LEVEL .NE. 2)  THEN
        IF (LHOOK) CALL DR_HOOK('XERRWV',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      STOP
      IF (LHOOK) CALL DR_HOOK('XERRWV',1,ZHOOK_HANDLE)
      ENDSUBROUTINE XERRWV
