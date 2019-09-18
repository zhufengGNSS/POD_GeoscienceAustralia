      SUBROUTINE chkbrdc(EPH,AVE,STD)
CC
CC NAME       :  chkbrdc
CC
CC PURPOSE    :  Check the broadcast reference epoch by using 3 to 4 sigma critierion
CC               
CC
CC PARAMETERS :
CC         IN :  EPH: Broadcast elements 
CC                           
CC        OUT :  Print out the status of dynamic elements 
CC
CC AUTHOR     :  Tzupang Tseng 
CC
CC
CC CREATED    :  02-08-2019 
CC
CC CHANGES    :  02-08-2019  T. Tseng: Introduce the STD values for quality control

      USE mdl_precision
      IMPLICIT NONE
C
C DECLARATIONS INSTEAD OF IMPLICIT
C --------------------------------
      INTEGER*4  ::N
      REAL(KIND = prec_q)::PI
      REAL(KIND = prec_q)::axis,cmm,ecc,ron,per,inc,ma,node
      REAL(KIND = prec_q)::EPH(20),AVE(8),STD(8)

      PI = 4*ATAN(1.D0)
      N = 4
C      
      axis  =EPH(3)
      ecc   =EPH(4)
      inc   =EPH(5)*180/PI
      node  =EPH(6)*180/PI
      per   =EPH(7)*180/PI
      ma    =EPH(8)*180/PI
      cmm   =EPH(9)
      ron   =EPH(10)*180/PI
C
C CHECK SEMI-MAJOR AXIS
C ---------------------
      IF (axis == 0.d0) RETURN
      IF(ABS(axis-AVE(1))>N*STD(1))
     1 PRINT*,'BAD SEMI-MAJOR AXIS, ABS(OBSERVED-MEAN)',
     2 ABS(axis-AVE(1))
     
C
C CHECK ECCENTRICITY
C ------------------
      IF(ABS(ecc-AVE(2))>N*STD(2))
     1 PRINT*,'BAD ECCENTRICITY, ABS(OBSERVED-MEAN)',
     2 ABS(ecc-AVE(2))
      
C
C CHECK INCLINATION
C ------------------
      IF(ABS(inc-AVE(3)*180/PI)>N*STD(3)*180/PI)
     1 PRINT*,'BAD INCLINATION, ABS(OBSERVED-MEAN)',
     2 ABS(inc-AVE(3)*180/PI)
     
C
C CHECK NODE
C ----------
      IF(ABS(node-AVE(4)*180/PI)>N*STD(4)*180/PI)
     1 PRINT*,'BAD NODE, ABS(OBSERVED-MEAN)',
     2 ABS(node-AVE(4)*180/PI)
     
C
C CHECK PERIGEE
C -------------
      IF(ABS(per-AVE(5)*180/PI)>N*STD(5)*180/PI)
     1 PRINT*,'BAD PERIGEE, ABS(OBSERVED-MEAN)',
     2 ABS(per-AVE(5)*180/PI)
      
C
C CHECK MEAN ANOMALY
C ------------------
      IF(ABS(ma-AVE(6)*180/PI)>N*STD(6)*180/PI)
     1 PRINT*,'BAD MEAN ANOMALY, ABS(OBSERVED-MEAN)',
     2 ABS(ma-AVE(6)*180/PI)
     
C
C CHECK CORRECTION TO MEAN MOTION
C -------------------------------
      IF(ABS(cmm-AVE(7))>N*STD(7))
     1 PRINT*,'BAD CORRECTION TO MEAN MOTION, ABS(OBSERVED-MEAN)', 
     2 ABS(cmm-AVE(7))
      
C
C CHECK RATE OF NODE
C ------------------
      IF(ABS(ron-AVE(8)*180/PI)>N*STD(8)*180/PI) 
     1 PRINT*,'BAD RATE OF NODE, ABS(OBSERVED-MEAN)',
     2 ABS(ron-AVE(8)*180/PI)
     
      END SUBROUTINE

