       SUBROUTINE glnorbint(EPH,TEPO,EPOEPH,X)
CC
CC NAME       :  glnorbint
CC
CC PURPOSE    :  COMPUTE GLONASS SATELLITE POSITION AND SATELLITE
CC               CLOCK CORRECTION FOR SPECIFIED EPOCH USING
CC               NUMERICAL INTEGRATION
CC
CC PARAMETERS :
CC         IN :  EPH    : NAVIGATION MESSAGE RECORD            R*8(16)
CC               TEPO   : EPOCH FOR NEW SATELLITE POSITION     R*8
CC               EPOEPH : EPOCH FOR GIVEN SATELLITE POSITION   R*8
CC        OUT :  X      : SATELITE POSITION FOR TEPO           R*8(3)
CC                        X(1)=X
CC                        X(2)=Y
CC                        X(3)=Z
CC               DTSATC : SATELLITE CLOCK CORRECTION FOR TEPO  R*8
CC                        (SEC)
CC
CC
C*
      IMPLICIT NONE
C
C DECLARATIONS INSTEAD OF IMPLICIT
C --------------------------------
      INTEGER*4 I     , IRCODE, ISTEP , LFNERR, NRSTEP
C
      REAL*8    DTSATC, EPOEPH, TEPO  , X0    , X1    , X2    , XI    ,
     1          XIH   , XII   , XIII  , XIIII , XK1   , XK2   , XK3   ,
     2          XK4   , XKK   , XVI   , XVIH  , XVII  , XVIII , XVIIII,
     3          Y0    , Y1    , Y2    , YI    , YIH   , YII   , YIII  ,
     4          YIIII , YK1   , YK2   , YK3   , YK4   , YKK   , YVI   ,
     5          YVIH  , YVII  , YVIII , YVIIII, Z0    , Z1    , Z2    ,
     6          ZI    , ZIH   , ZII   , ZIII  , ZIIII , ZK1   , ZK2   ,
     7          ZK3   , ZK4   , ZKK   , ZVI   , ZVIH  , ZVII  , ZVIII ,
     8          ZVIIII, ACC
C
C
C PARAMETERS
C MAXINT : MAXIMUM NUMBER OF INTEGRATION STEPS
C
C DECLARATIONS
C ------------
      REAL*8       EPH(*),X(3)
      REAL*8       H,HSTEP(500),INTERV
C
C COMMON BLOCKS
C -------------
C
      DATA INTERV/10.D0/
C INTERV : INTERVAL FOR ONE INETRAGTION STEP (SEC)
C
C X COORDINATE IN METER
C ---------------------
      X0=EPH(5)*1000
      X1=EPH(6)*1000
      X2=EPH(7)*1000
C
C Y COORDINATE IN METER
C ---------------------
      Y0=EPH(9) *1000
      Y1=EPH(10)*1000
      Y2=EPH(11)*1000
C
C Z COORDINATE IN METER
C ---------------------
      Z0=EPH(13)*1000
      Z1=EPH(14)*1000
      Z2=EPH(15)*1000
C
C SATELLITE CLOCK
C ---------------
!      IF (EPH(2).NE.0D0.OR.EPH(3).NE.0D0) THEN
!        DTSATC=EPH(2)+EPH(3)*(TEPO-EPOEPH)*86400.D0
!      ENDIF
C
C NUMERICAL INTEGRATION
C ---------------------
C
C FORMULA ACCORDING GLONASS INTERFACE CONTROL DOCUMENT
C
C RUNGE-KUTTA-METHOD FOR INTEGRATION, NOTATION ACCORDING BRONSTEIN
C
C INITIALIZATION
C --------------
      XI=X0
      YI=Y0
      ZI=Z0
      XVI=X1
      YVI=Y1
      ZVI=Z1
C
C DEFINE INTEGRATION INTERVAL
C ---------------------------
C     NRSTEP = NUMBER OF ITERATION STEPS
C     HSTEP   = INTERVAL FOR EACH ITERATION STEP (SECONDS)
      NRSTEP=IABS(IDINT((TEPO-EPOEPH)/INTERV))+1
C      PRINT*,'TEPO-EPOEPH =', TEPO-EPOEPH
C      PRINT*,' NUMBER OF ITERATION STEPS =', NRSTEP

      DO 100 I=1,NRSTEP-1
        HSTEP(I)=INTERV
        IF(TEPO-EPOEPH.LT.0.D0)HSTEP(I)=HSTEP(I)*(-1.D0)
100   CONTINUE
      HSTEP(NRSTEP)=DMOD((TEPO-EPOEPH),INTERV)
C      PRINT*,'HSTEP(NRSTEP) =',HSTEP(NRSTEP), HSTEP(I) 
C
C LOOP FOR ALL INTEGRATION STEPS
C ------------------------------
      DO 3000 ISTEP=1,NRSTEP
        H=HSTEP(ISTEP)
C
C INTEGRATION FOR VELOCITY
C ------------------------
C FORMULA TO COMPUTE NEW VELOCITY:
C V(I+1)=V(I)+DV/DT*DT
C
        CALL glnacc(1,XI,YI,ZI,XVI,YVI,ACC)
        XK1=ACC+X2
        CALL glnacc(2,XI,YI,ZI,XVI,YVI,ACC)
        YK1=ACC+Y2
        CALL glnacc(3,XI,YI,ZI,XVI,YVI,ACC)
        ZK1=ACC+Z2
        XVII=XVI+XK1*H/2
        YVII=YVI+YK1*H/2
        ZVII=ZVI+ZK1*H/2
        XII=XI+XVII*H/2
        YII=YI+YVII*H/2
        ZII=ZI+ZVII*H/2
C
        CALL glnacc(1,XII,YII,ZII,XVII,YVII,ACC)
        XK2=ACC+X2
        CALL glnacc(2,XII,YII,ZII,XVII,YVII,ACC)
        YK2=ACC+Y2
        CALL glnacc(3,XII,YII,ZII,XVII,YVII,ACC)
        ZK2=ACC+Z2
        XVIII=XVI+XK2*H/2
        YVIII=YVI+YK2*H/2
        ZVIII=ZVI+ZK2*H/2
        XIII=XI+XVIII*H/2
        YIII=YI+YVIII*H/2
        ZIII=ZI+ZVIII*H/2
C
        CALL glnacc(1,XIII,YIII,ZIII,XVIII,YVIII,ACC)
        XK3=ACC+X2
        CALL glnacc(2,XIII,YIII,ZIII,XVIII,YVIII,ACC)
        YK3=ACC+Y2
        CALL glnacc(3,XIII,YIII,ZIII,XVIII,YVIII,ACC)
        ZK3=ACC+Z2
        XVIIII=XVI+XK3*H
        YVIIII=YVI+YK3*H
        ZVIIII=ZVI+ZK3*H
        XIIII=XI+XVIIII*H
        YIIII=YI+YVIIII*H
        ZIIII=ZI+ZVIIII*H
C
        CALL glnacc(1,XIIII,YIIII,ZIIII,XVIIII,YVIIII,ACC)
        XK4=ACC+X2
        CALL glnacc(2,XIIII,YIIII,ZIIII,XVIIII,YVIIII,ACC)
        YK4=ACC+Y2
        CALL glnacc(3,XIIII,YIIII,ZIIII,XVIIII,YVIIII,ACC)
        ZK4=ACC+Z2
        XKK=(XK1+2*XK2+2*XK3+XK4)/6
        YKK=(YK1+2*YK2+2*YK3+YK4)/6
        ZKK=(ZK1+2*ZK2+2*ZK3+ZK4)/6
        XVIH=XVI+XKK*H
        YVIH=YVI+YKK*H
        ZVIH=ZVI+ZKK*H
C
C INTEGRATION FOR POSITION
C ------------------------
C FORMULA TO COMPUTE NEW POSITION:
C P(I+1)=P(I)+DP/DT*DT
C
        XK1=XVI
        YK1=YVI
        ZK1=ZVI
C
        XK2=XVII
        YK2=YVII
        ZK2=ZVII
C
        XK3=XVIII
        YK3=YVIII
        ZK3=ZVIII
C
        XK4=XVIIII
        YK4=YVIIII
        ZK4=ZVIIII
C
        XKK=(XK1+2*XK2+2*XK3+XK4)/6
        YKK=(YK1+2*YK2+2*YK3+YK4)/6
        ZKK=(ZK1+2*ZK2+2*ZK3+ZK4)/6
        XIH=XI+XKK*H
        YIH=YI+YKK*H
        ZIH=ZI+ZKK*H
C
C
C INITIALIZATION OF NEW INTEGRATION STEP
C --------------------------------------
        XVI=XVIH
        YVI=YVIH
        ZVI=ZVIH
        XI=XIH
        YI=YIH
        ZI=ZIH
C
C NEXT INTEGRATION STEP
3000  CONTINUE
C
C WRITE RESULTS
C -------------
      X(1)=XI
      X(2)=YI
      X(3)=ZI

      RETURN
      END SUBROUTINE

