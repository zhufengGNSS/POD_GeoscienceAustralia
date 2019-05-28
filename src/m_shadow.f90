MODULE m_shadow


! ----------------------------------------------------------------------
! MODULE: m_shadow
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the m_shadow subroutine
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng, Geosceince Australia, Australia 
!               
!               
! Created:      16-04-2019
! ----------------------------------------------------------------------

      IMPLICIT NONE
      !SAVE


Contains


SUBROUTINE shadow ( r_sat, r_sun, r_moon, lambda )


! ----------------------------------------------------------------------
! SUBROUTINE: shadow.f90
! ----------------------------------------------------------------------
! Purpose:
! This subrutine is used to compute the coefficient for scaling the solar radiation pressure 
! from penumbra to umbra. The concept of this model is to use the ratio of the celestial radius 
! over the distance (the celestial to the satellite) as a threshold for the eclipsing judgement. 
! Such a ratio forms an apparent angle (in radian) viewed from the satellite. 
! By comparing with different apparent angles as viewed from the satellite, the satellite can be 
! judged in the sun light, penumbra or umbra.
!
! This concept can be easily understood through three illustrations: a sideview,
! a look-down view from the satellite and a sideview of the vector yy.    
! ----------------------------------------------------------------------
! Input arguments:
! - r_sat        : satellite position vector (m)
! - r_sun        : Sun position vector wrt the earth
! - r_moon       : Moon position vector wrt the earth
! 
! Output arguments:
! - lambda       : Coefficient for scaling the SRP-induced acceleration acting
!                  on the eclipsed GNSS satellites
!                  
!                  lambda = 1     : In SUN LIGHT AREA
!                         = 0     : In UMBRA AREA (full eclipse)
!                  0 < lambda < 1 : In PENUMBRA AREA 
!- ECLTYP        : ECLIPSE TYPE = 'E' : eclipsed by the earth shadow
!                                 'M' : eclipsed by the moon  shadow 
!
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	16-04-2019
!
! Changes:    
! 
! Copyright:  GEOSCIENCE AUSTRALIA
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_get_lambda
      IMPLICIT NONE

! ----------------------------------------------------------------------

      REAL (KIND = prec_q), DIMENSION(3),INTENT(IN) :: r_sat,r_moon
      REAL (KIND = prec_q), DIMENSION(3),INTENT(IN) :: r_sun
      REAL (KIND = prec_q), INTENT(OUT) :: lambda
! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      CHARACTER (KIND = 1) :: ECLTYP
      REAL (KIND = prec_q) :: AU,Pi
      REAL (KIND = prec_q) :: sunrad, moonrad, ertrad
      REAL (KIND = prec_q) :: Dsun, Dmoonsun
      REAL (KIND = prec_q) :: Ds,sclfa
      REAL (KIND = prec_q) :: AB, CD
      REAL (KIND = prec_q) :: rp, rs, xx
      REAL (KIND = prec_q), DIMENSION(3) :: ud 
      REAL (KIND = prec_q), DIMENSION(3) :: yy
      REAL (KIND = prec_q), DIMENSION(3) :: Rmoonsun, Rsatmoon
      INTEGER              :: i,j,k
      INTEGER (KIND = 4) :: idbprt,idebug
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Numerical Constants
      Pi = 4*atan(1.0d0)
      sunrad = 696.d6         !Units = m
      moonrad = 1738.d3       !Units = m
      ertrad = 6378.13630d3   !Units = m
      lambda = 1.0d0
! ---------------------------------------------------------------------
! The distance between the earth and sun
      Dsun=sqrt(r_sun(1)**2+r_sun(2)**2+r_sun(3)**2)

! The unit vector of the satellite wrt the sun (SUN->SAT)
      Ds=sqrt((r_sun(1)-r_sat(1))**2+(r_sun(2)-r_sat(2))**2+(r_sun(3)-r_sat(3))**2)
      ud(1)=(r_sat(1)-r_sun(1))/Ds
      ud(2)=(r_sat(2)-r_sun(2))/Ds
      ud(3)=(r_sat(3)-r_sun(3))/Ds

! The satellite is in the sunlight area.
IF(Ds.le.Dsun) return

! Project the satellite position vector onto the unit vector SUN->SAT
      
      CALL productdot(r_sat,ud,AB)

! A vector resulted from the cross product of the position vector EARTH->SAT and the unit vector SUN->SAT 

      CALL productcross(r_sat,ud,yy)

! rs  : apparent angle of sun as viewed from satellite (radians)
! rp  : apparent angle of eclipsing body as viewed from satellite (radians)
! xx  : apparent separation of the center of the Sun and eclipsing body (radians)
!     rs+rp <= xx : no eclipse
!    
!     rp-rs >= xx : full eclipse

     xx=sqrt(yy(1)**2+yy(2)**2+yy(3)**2)/AB
     rs=sunrad/Ds 
     rp=ertrad/AB 

! Check the earth shadow  

      CALL get_lambda(rs, rp, xx, lambda, idbprt, idebug) 

    IF (lambda .lt. 1.0d0) THEN
    ECLTYP = 'E'
!    PRINT*,' ECLTYP =  ', ECLTYP
    RETURN

    ELSE
! Check the lunar eclipse
          DO i=1,3
          Rmoonsun(i) = r_moon(i) - r_sun(i)
          Rsatmoon(i) = r_sat(i) - r_moon(i)
          END DO 

          Dmoonsun = dsqrt(Rmoonsun(1)**2+Rmoonsun(2)**2+Rmoonsun(3)**2)

          IF(Ds.le.Dmoonsun) return

          CALL productdot(Rsatmoon,ud,CD)
          CALL productcross (Rsatmoon,ud,yy)
          xx=sqrt(yy(1)**2+yy(2)**2+yy(3)**2)/CD

          rs=sunrad/Ds 
          rp=moonrad/CD 
          CALL get_lambda(rs, rp, xx, lambda, idbprt, idebug)

          IF( lambda.lt.1.d0 ) ECLTYP = 'M'
!          PRINT*,' ECLTYP =  ', ECLTYP

    END IF


END SUBROUTINE

END MODULE 
