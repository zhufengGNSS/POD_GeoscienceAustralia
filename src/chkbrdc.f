      SUBROUTINE chkbrdc(EPH,AVE,STD)
CC
CC NAME       :  chkbrdc_gps
CC
CC PURPOSE    :  CHECK ONE EPHEMERIS MESSAGE FOR GARBAGE. SET STATUS
CC               ACCORDING TO THE RESULT OF THE CHECK
CC
CC PARAMETERS :
CC         IN :  EPH    : EPHEMERIDES INFORMATION              R*8(20)
CC                        EPH(I):
CC                          I: EPHEMERIDE ELEMENT
CC                          EPH(1) : GPS-WEEK
CC                          EPH(2) : T0E
CC                          EPH(3) : A
CC                          EPH(4) : E
CC                          EPH(5) : I
CC                          EPH(6) : R.A. OF ASCENDING NODE
CC                          EPH(7) : PERIGEE
CC                          EPH(8) : MEAN ANOMALY (T0E)
CC                          EPH(9) : DN (CORRECTION TO MEAN MOTION)
CC                          EPH(10): RATE OF NODE
CC                          EPH(11): CUS
CC                          EPH(12): CUC
CC                          EPH(13): CRS
CC                          EPH(14): CRC
CC                          EPH(15): CIS
CC                          EPH(16): CIC
CC                          EPH(17): AODE
CC                          EPH(18): IDOT
CC                          EPH(19): NOT USED
CC                          EPH(20): NOT USED
CC                          EPH(21): GPS WEEK OF THE NEXT EPHEMERIDE
CC                              :        :
CC        OUT :  STATUS : STATUS OF THE MESSAGE (RESULT OF    CH*8
CC                          THE CHECKS PERFORMED)
CC
      IMPLICIT NONE
C
C DECLARATIONS INSTEAD OF IMPLICIT
C --------------------------------
      INTEGER*4 IWEEK, N
C
      REAL*8    A    , DN   , E    , ODOT , PER  , T0E  , XI   , XM0  ,
     1          XNODE, PI
C
      REAL*8      EPH(20),AVE(8),STD(8)
      CHARACTER*8 STATUS
C
C
C INITIALIZE STATUS
C -----------------
      STATUS=' '
      PI = 4*ATAN(1.D0)
      N = 4
C      
C SET EPHEMERIS ELEMENTS
C ----------------------
      IWEEK=IDNINT(EPH(1))
      T0E  =EPH(2)
      A    =EPH(3)
      E    =EPH(4)
      XI   =EPH(5)*180/PI
      XNODE=EPH(6)*180/PI
      PER  =EPH(7)*180/PI
      XM0  =EPH(8)*180/PI
      DN   =EPH(9)
      ODOT =EPH(10)*180/PI
C
C CHECK GPS WEEK
C --------------
      IF (IWEEK == 0) RETURN
      IF(IWEEK.LE.0.OR.IWEEK.GT.2000) THEN
        STATUS='BAD WEEK'
        PRINT*,'STATUS =', STATUS, 'GPSWEEK =', IWEEK
        GOTO 100
      ENDIF
C
C CHECK T0E
C ---------
!      IF(T0E.LT.0.D0 .OR. T0E.GT.7*86400.D0 .OR.
!     1   (DABS(DMOD(T0E,100.D0)-3.D0).GT.1.D-5 .AND.
!     2   (DMOD(T0E+1.D-5,100.D0)).GT.2.D-5)) THEN
!        STATUS='BAD T0E'
!        PRINT*,'STATUS =', STATUS
!        GOTO 100
!      ENDIF
C
C CHECK A
C -------
      IF (A == 0.d0) RETURN
!      IF(A.LT.26.0D6.OR.A.GT.27.0D6) THEN
      IF(ABS(A-AVE(1))>N*STD(1))THEN
        STATUS='BAD A'
        PRINT*,'STATUS =',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(A-AVE(1))
        GOTO 100
      ENDIF
C
C CHECK E
C -------
!      IF(E.LT.0.D0.OR.E.GT.0.1D0) THEN
      IF(ABS(E-AVE(2))>N*STD(2))THEN
        STATUS='BAD E'
        PRINT*,'STATUS =',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(E-AVE(2))
        GOTO 100
      ENDIF
C
C CHECK I
C -------
      !IF((XI.GT.65.D0.OR.XI.LT.60.D0).AND.(XI.GT.58.OR.XI.LT.50)) THEN
      IF(ABS(XI-AVE(3)*180/PI)>N*STD(3)*180/PI)THEN
        STATUS='BAD I'
        PRINT*,'STATUS =',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(XI-AVE(3)*180/PI)
        GOTO 100
      ENDIF
C
C CHECK NODE
C ----------
      !IF(DABS(XNODE).LT.1.D-6) THEN
      IF(ABS(XNODE-AVE(4)*180/PI)>N*STD(4)*180/PI)THEN
        STATUS='BAD NODE'
        PRINT*,'STATUS =',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(XNODE-AVE(4)*180/PI)
        GOTO 100
      ENDIF
C
C CHECK PERIGEE
C -------------
      !IF(DABS(PER).LT.1.D-6) THEN
      IF(ABS(PER-AVE(5)*180/PI)>N*STD(5)*180/PI)THEN
        STATUS='BAD PERI'
        PRINT*,'STATUS=',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(PER-AVE(5)*180/PI)
        GOTO 100
      ENDIF
C
C CHECK MEAN ANOMALY
C ------------------
      !IF(DABS(XM0).LT.1.D-6) THEN
      IF(ABS(XM0-AVE(6)*180/PI)>N*STD(6)*180/PI)THEN
        STATUS='BAD M0'
        PRINT*,'STATUS=',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(XM0-AVE(6)*180/PI)
        GOTO 100
      ENDIF
C
C CHECK CORRECTION TO MEAN MOTION
C -------------------------------
      !IF(DN.LT.0.01D-8.OR.DN.GT.0.70D-8) THEN
      IF(ABS(DN-AVE(7))>N*STD(7))THEN
        STATUS='BAD DN'
        PRINT*,'STATUS =',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)', ABS(DN-AVE(7))
        GOTO 100
      ENDIF
C
C CHECK RATE OF NODE
C ------------------
      !IF(ODOT.LT.-0.60D-6.OR.ODOT.GT.-0.25D-6) THEN
      IF(ABS(ODOT-AVE(8)*180/PI)>N*STD(8)*180/PI)THEN
        STATUS='BAD ODOT'
        PRINT*,'STATUS=',STATUS
        PRINT*,'ABS(OBSERVED-MEAN)',ABS(ODOT-AVE(8)*180/PI)
        GOTO 100
      ENDIF
C
C END
C ---
100   CONTINUE
      RETURN
      END SUBROUTINE

