C*
      SUBROUTINE brdc2ecef(GPST,EPH,ECEFPOS)
CC
CC NAME       :  brdc2ecef
CC
CC PURPOSE    :  Convert the broadcast elements to ECEF coordinates  
CC
CC PARAMETERS :
CC         IN :  GPST     :  GPS time                   
CC               EPH      :  broadcast parameters 
CC               EPH(2)   :  toe
CC               EPH(3)   :  semi-major axis
CC               EPH(4)   :  eccentricity
CC               EPH(5)   :  inclination
CC               EPH(6)   :  right ascension of ascending node
CC               EPH(7)   :  perigee
CC               EPH(8)   :  mean anomaly
CC               EPH(9)   :  correction to mean motion
CC               EPH(10)  :  rate of node
CC               EPH(11)  :  cus
CC               EPH(12)  :  cuc
CC               EPH(13)  :  crs
CC               EPH(14)  :  crc
CC               EPH(15)  :  cis
CC               EPH(16)  :  cic
CC               EPH(18)  :  idot

CC        OUT :  ECEFPOS  :  ECEF coordinates  
CC
CC Author: Tzupang Tseng
CC
CC Create: 18-09-2019
CC
CC Corpyright: Geoscience Australia

C*
      USE mdl_precision
      IMPLICIT NONE

      INTEGER :: I
      REAL (KIND = prec_d) :: GPST
      REAL (KIND = prec_d) :: EPH(20),ECEFPOS(3)
      REAL (KIND = prec_d) :: axis, toe, DT
      REAL (KIND = prec_d) :: u0, u
      REAL (KIND = prec_d) :: Et, e
      REAL (KIND = prec_d) :: v 
      REAL (KIND = prec_d) :: cus,cuc,crs,crc,cis,cic
      REAL (KIND = prec_d) :: w0, inc0, idot 
      REAL (KIND = prec_d) :: w, it, r
      REAL (KIND = prec_d) :: ascnode0, ascnode,dot_ascnode
      REAL (KIND = prec_d) :: GM, OMEGA

      OMEGA  = 7292115.1467D-11 ! RAD/SEC
      GM     = 398.6004415D12   ! M**3/SEC**2
    
C
C Convert the broadcast to ECEF coordinates 
C -----------------------------------------

!  Time elapse since toe
! --------------------------
      toe = EPH(2)
      DT=GPST-toe
      IF(DT.GT. 302400.D0)DT=DT-604800.D0
      IF(DT.LT.-302400.D0)DT=DT+604800.D0

! Mean anomaly at GPST
! --------------------
      axis = EPH(3)
      u0   = EPH(8)
      u    =u0+(sqrt(GM/axis**3)+EPH(9))*DT

! Iteration solution of eccentric anomaly
! ---------------------------------------
      Et=u
      e=EPH(4)
      DO I=1,10
        Et=u+e*sin(Et)
      ENDDO

! True anomaly
! -------------
      v = atan(sqrt(1.d0-e**2)*sin(Et)/(cos(Et)-e))

! Correct for orbital perturbations
! ----------------------------------
      cus =  EPH(11)
      cuc =  EPH(12)
      crs =  EPH(13)
      crc =  EPH(14)
      cis =  EPH(15)
      cic =  EPH(16)

      w0   = EPH(7) + v
      inc0 = EPH(5)
      idot = EPH(18)
! Argument of perigee
      w = w0               + cuc*cos(2.D0*w0)+cus*sin(2.D0*w0)
! Radial distance
      r=axis*(1-e*cos(Et)) + crc*cos(2.D0*w0)+crs*sin(2.D0*w0)
! Inclination
      it = inc0 + idot*DT  + cic*cos(2.D0*w0)+cis*sin(2.D0*w0)

! Compute the right ascension of ascending node
! ---------------------------------------------
      ascnode0     = EPH(6)
      dot_ascnode  = EPH(10)
      ascnode=ascnode0+(dot_ascnode-OMEGA)*DT-OMEGA*toe

! Convert the Keplerian parameters to ECEF coordinates
! ------------------------------------------------------
      ECEFPOS(1)=r*cos(w)*cos(ascnode)-r*sin(w)*cos(it)*sin(ascnode)
      ECEFPOS(2)=r*cos(w)*sin(ascnode)+r*sin(w)*cos(it)*cos(ascnode)
      ECEFPOS(3)=r*sin(w)*sin(it)
 
      END SUBROUTINE
