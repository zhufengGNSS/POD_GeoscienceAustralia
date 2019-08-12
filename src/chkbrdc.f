      SUBROUTINE chkbrdc(EPH1,EPH2)
CC
CC NAME       :  chkbrdc
CC
CC PURPOSE    :  CHECK TWO BROADCAST MESSAGES FOR CONSISTENCY
CC
CC PARAMETERS :
CC         IN :  EPH1   : EPHEMERIDES INFORMATION (FIRST SET)  R*8(20)
CC                        EPH1(I):
CC                          I: EPHEMERIDE ELEMENT
CC                          EPH1(1) : GPS-WEEK
CC                          EPH1(2) : T0E
CC                          EPH1(3) : A
CC                          EPH1(4) : E
CC                          EPH1(5) : I
CC                          EPH1(6) : R.A. OF ASCENDING NODE
CC                          EPH1(7) : PERIGEE
CC                          EPH1(8) : MEAN ANOMALY (T0E)
CC                          EPH1(9) : DN (CORRECTION TO MEAN MOTION)
CC                          EPH1(10): RATE OF NODE
CC                          EPH1(11): CUS
CC                          EPH1(12): CUC
CC                          EPH1(13): CRS
CC                          EPH1(14): CRC
CC                          EPH1(15): CIS
CC                          EPH1(16): CIC
CC                          EPH1(17): AODE
CC                          EPH1(18): IDOT
CC                          EPH1(19): NOT USED
CC                          EPH1(20): NOT USED
CC                          EPH1(21): GPS WEEK OF THE NEXT EPHEMERIDE
CC                              :        :
CC               EPH2   : EPHEMERIDES INFORMATION (SECOND SET) R*8(20)
CC                          (SEE ARRAY EPH1)
CC
CC
C*
      IMPLICIT NONE
C
C DECLARATIONS INSTEAD OF IMPLICIT
C --------------------------------
      INTEGER*4 IWEEK1, IWEEK2
C
      REAL*8    A1    , A2    , DA    , DE    , DI    , DM0   , DN1   ,
     1          DNODE , DPER  , DT    , E1    , E2    , HOURS , PER1  ,
     2          PER2  , PER21 , T0E1  , T0E2  , XI1   , XI2   , XIDOT1,
     3          XM01  , XM02  , XM021 , XN    , XNDOT1, XNOD21, XNODE1,
     4          XNODE2, XNODEW
C
CCC       IMPLICIT REAL*8(A-H,O-Z)
CCC       IMPLICIT INTEGER*4 (I-N)
      REAL*8      EPH1(20),EPH2(20)
      CHARACTER*8 STATUS
      REAL*8 GM, OMEGA, PI
C
      PI = 4*ATAN(1.D0)
      OMEGA  = 7292115.1467D-11 ! RAD/SEC
      GM     = 398.6004415D12   !  M**3/SEC**2

C
C
C FIRST SET OF ELEMENTS
C ---------------------
      IWEEK1=IDNINT(EPH1(1))
      T0E1  =EPH1(2)
      A1    =EPH1(3)
      E1    =EPH1(4)
      XI1   =EPH1(5)*180.D0/PI
      XNODE1=EPH1(6)*180.D0/PI
      PER1  =EPH1(7)*180.D0/PI
      XM01  =EPH1(8)*180.D0/PI
      DN1   =EPH1(9)
      XNDOT1=EPH1(10)
      XIDOT1=EPH1(18)
C
C SECOND SET OF ELEMENTS
C ----------------------
      IWEEK2=IDNINT(EPH2(1))
      T0E2  =EPH2(2)
      A2    =EPH2(3)
      E2    =EPH2(4)
      XI2   =EPH2(5)*180.D0/PI
      XNODE2=EPH2(6)*180.D0/PI
C
C CORRECT NODE FOR SYSTEM CHANGE AT END OF GPS WEEK
      IF(IWEEK1.NE.IWEEK2) THEN
        XNODEW=DMOD((IWEEK2-IWEEK1)*86400.D0*7.D0*OMEGA*180.D0
     1              /PI,360.D0)
        XNODE2=XNODE2+XNODEW
      ENDIF
      PER2  =EPH2(7)*180.D0/PI
      XM02  =EPH2(8)*180.D0/PI
      IF(XM02.LE.-180.D0) XM02=XM02+360.D0
      IF(XM02.GT.+180.D0) XM02=XM02-360.D0
C
C COMPUTE CHANGE IN FIRST ELEMENTS IN THE TIME INTERVAL (T0E1,T0E2)
C -----------------------------------------------------------------
      DT=(IWEEK2-IWEEK1)*86400.D0*7+T0E2-T0E1
      HOURS=DT/3600.D0
      XI1=XI1+XIDOT1*DT*180.D0/PI
      XNODE1=XNODE1+XNDOT1*DT*180.D0/PI
      XN=DSQRT(GM/A1**3)+DN1
      XM01=XM01+XN*DT*180.D0/PI
      IF(XM01.LE.-180.D0) XM01=XM01+360.D0
      IF(XM01.GT.+180.D0) XM01=XM01-360.D0
C
C CHECK DIFFERENCE IN A
C ---------------------
      DA=100.D0+10.D0*HOURS
      IF(DABS(A2-A1).GT.DA) THEN
        PRINT*,'BAD DA'
        GOTO 100
      ENDIF
C
C CHECK DIFFERENCE IN E
C ---------------------
      DE=5.D-6+0.5D-6*HOURS
      IF(DABS(E2-E1).GT.DE) THEN
        PRINT*,'BAD DE'
        GOTO 100
      ENDIF
C
C CHECK DIFFERENCE IN I
C ---------------------
      DI=0.0005D0+0.0003D0*HOURS
      IF(DABS(XI2-XI1).GT.DI) THEN
        PRINT*,'BAD DI'
        GOTO 100
      ENDIF
C
C CHECK DIFFERENCE IN NODE
C ------------------------
      XNOD21=DABS(XNODE2-XNODE1)
      IF (XNOD21.GT.100.D0) XNOD21=DABS(XNOD21-360.D0)
      DNODE=0.0002D0+0.0005D0*HOURS
      IF(XNOD21.GT.DNODE) THEN
        PRINT*,'BAD DNOD'
        GOTO 100
      ENDIF
C
C CHECK DIFFERENCE IN PERIGEE
C ---------------------------
      PER21=DABS(PER2-PER1)
      IF (PER21.GT.100.D0) PER21=DABS(PER21-360.D0)
      DPER=2.0D0+0.01D0*HOURS
      IF (EPH1(4).LT.0.001) DPER=10.0D0+0.01D0*HOURS
      IF(PER21.GT.DPER) THEN
        PRINT*,'BAD DPER'
        GOTO 100
      ENDIF
C
C CHECK DIFFERENCE IN MEAN ANOMALY
C --------------------------------
      XM021=DABS(XM02-XM01)
      IF(XM021.GT.100.D0) XM021=DABS(XM021-360.D0)
      DM0=2.00D0+0.08D0*HOURS
      IF (EPH1(4).LT.0.001) DM0=10.0D0+0.01D0*HOURS
      IF(XM021.GT.DM0) THEN
        PRINT*,'BAD DM0'
        GOTO 100
      ENDIF
C
C END
C ---
100   CONTINUE
      RETURN
      END SUBROUTINE

