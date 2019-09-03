      SUBROUTINE glnacc(ITYP,X,Y,Z,VX,VY,ACC)
CC
CC NAME       : glnacc
CC
CC
CC PURPOSE    :  COMPUTES GLONASS ACCELERATION FOR SATELLITE POSITIONS 
CC               PROPAGATION IN ECEF SYSTEM (REFER TO BERNESE GLDVDT.f)
CC
CC PARAMETERS :
CC         IN :  ITYP   : =1 FOR X-COMPONENT                  I*4
CC                      : =2 FOR Y-COMPONENT                  I*4
CC                      : =3 FOR Z-COMPONENT                  I*4
CC               X      : X-COORDINATE OF SATELLITE POSITION  R*8
CC               Y      : Y-COORDINATE OF SATELLITE POSITION  R*8
CC               Z      : Z-COORDINATE OF SATELLITE POSITION  R*8
CC               VX     : VELOCITY IN X DIRECTION             R*8
CC               VY     : VELOCITY IN Y DIRECTION             R*8
CC        OUT :  ACC    : ACCELERATION                        R*8
CC
CC AUTHOR     :  H.HABRICH
CC
CC CREATED    :  13-FEB-1997

C*
      IMPLICIT NONE
      REAL*8    X,Y,Z,VX,VY,GM,AE,OMEGA,C20,ACC,T1,T2,T3,T4,R,D
      INTEGER*4 ITYP
C
      GM=398.60044D12
      AE=6378136.D0
      OMEGA=7292115.D-11
      C20=-1082.63D-6
C
      D=X*X+Y*Y+Z*Z
      R=DSQRT(D)
C
      IF(ITYP.EQ.1) THEN
        T1=-GM/(R*R*R)*X
        T2=3.D0/2.D0*C20*(GM*AE*AE)/(R*R*R*R*R)*X*(1-(5*Z*Z)/(R*R))
        T3=OMEGA*OMEGA*X
        T4=2*OMEGA*VY
        ACC=T1+T2+T3+T4
      ELSEIF(ITYP.EQ.2)THEN
        T1=-GM/(R*R*R)*Y
        T2=3.D0/2.D0*C20*(GM*AE*AE)/(R*R*R*R*R)*Y*(1-(5*Z*Z)/(R*R))
        T3=OMEGA*OMEGA*Y
        T4=-2*OMEGA*VX
        ACC=T1+T2+T3+T4
      ELSEIF(ITYP.EQ.3)THEN
        T1=-GM/(R*R*R)*Z
        T2=3.D0/2.D0*C20*(GM*AE*AE)/(R*R*R*R*R)*Z*(3-(5*Z*Z)/(R*R))
        T3=0.
        T4=0.
        ACC=T1+T2+T3+T4
      ENDIF
      RETURN


      END 
