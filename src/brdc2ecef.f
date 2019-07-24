C*
      SUBROUTINE brdc2ecef(T,EPH,X)
CC
CC NAME       :  brdc2ecef
CC
CC PURPOSE    :  COMPUTE EARTH FIXED COORDINATES X(I),I=1,2,3 IN
CC               WGS-84 SYSTEM FOR A SATELLITE WITH BROADCAST EPHEMERIS
CC               PARAMETERS GIVEN IN ARRAY EPH(I),I=1,2,..,20 
CC
CC PARAMETERS :
CC         IN :  T     :  ACTUAL TIME IN GPS-TIME                  
CC               EPH   :  ARRAY WITH 20 ELEMENTS CONTAINING THE    
CC                        BROADCAST PARAMETERS
CC        OUT :  X(I),I=1,2,3:  EARTH FIXED SATELLITE COORDINATES 
CC
CC
CC
C*
      IMPLICIT NONE

      INTEGER*4 I
C
      REAL*8    A    , CP   , DI   , DR   , DT   , DU   , E    , EX   ,
     1          PHI  , R    , SP   , T    , U    , V    , XI   , XM   ,
     2          XN   , XNODE, XP   , YP
C
      REAL*8 EPH(20),X(3)
      REAL*8 GM, OMEGA

      OMEGA  = 7292115.1467D-11 ! RAD/SEC
      GM     = 398.6004415D12   !  M**3/SEC**2
    
C
C COMPUTE SATELLITE POSITION
C --------------------------
      A=EPH(3)
      XN=DSQRT(GM/A**3)
      DT=T-EPH(2)
      IF(DT.GT. 302400.D0)DT=DT-604800.D0
      IF(DT.LT.-302400.D0)DT=DT+604800.D0
      XN=XN+EPH(9)
      XM=EPH(8)+XN*DT
      EX=XM
      E=EPH(4)
      DO I=1,10
        EX=XM+E*DSIN(EX)
      ENDDO
      V=2.D0*DATAN(DSQRT((1.D0+E)/(1.D0-E))*DTAN(EX/2.D0))
      PHI=V+EPH(7)
      CP=DCOS(2.D0*PHI)
      SP=DSIN(2.D0*PHI)
      DU=EPH(11)*SP+EPH(12)*CP
      DR=EPH(14)*CP+EPH(13)*SP
      DI=EPH(16)*CP+EPH(15)*SP
      R=A*(1-E*DCOS(EX))+DR
      U=PHI+DU
      XI=EPH(5)+DI+EPH(18)*DT
      XP=R*DCOS(U)
      YP=R*DSIN(U)
      XNODE=EPH(6)+(EPH(10)-OMEGA)*DT-OMEGA*EPH(2)
      X(1)=XP*DCOS(XNODE)-YP*DCOS(XI)*DSIN(XNODE)
      X(2)=XP*DSIN(XNODE)+YP*DCOS(XI)*DCOS(XNODE)
      X(3)=YP*DSIN(XI)
      RETURN
      END SUBROUTINE
