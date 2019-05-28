
SUBROUTINE pdsrp (GM,prnnum,BLKNUM,srpid,r,v,r_sun,beta,del_u,Asrp )


! ----------------------------------------------------------------------
! SUBROUTINE: pdsrp.f90
! ----------------------------------------------------------------------
! Purpose:
! This subrutine is used to compute the partial derivative of force
! w.r.t. SRP parameters. 
! a box-wing model (srpid=2) and a ECOM (D2B1) model (srpid=3) 
! ----------------------------------------------------------------------
! Input arguments:
! - prnnum       : satellite PRN number 
! - srpid        : =2: box-wing model;
!                  =3: ECOM (D2B1)model
! - r            : satellite position vector (m)
! - v            : satellite velocity vector
! - r_sun        : Sun position vector
! 
! Output arguments:
! - Asrp         : Partial derivative output 
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	08-Jan-2018
! 
! Copyright:  GEOSCIENCE AUSTRALIA
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      INTEGER                           :: srpid
      INTEGER                           :: prnnum,BLKNUM
      REAL (KIND = prec_q), DIMENSION(3) :: r,v,r_sun

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: Ps,AU,Pi
      REAL (KIND = prec_q) :: Cr,ab,Lab
      REAL (KIND = prec_q) :: ANG,GM,E,Edot
      REAL (KIND = prec_q) :: Ds,sclfa

      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: surforce(4,3)
      REAL (KIND = prec_q) :: Asrp(3,9)
      REAL (KIND = prec_q), DIMENSION(3) :: er,ed,ey,eb,ex
      REAL (KIND = prec_q), DIMENSION(3) :: fsrp
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(9) :: srpcoef
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q), DIMENSION(4) :: AREA1,REFL1,DIFU1,ABSP1
      INTEGER              :: zta
      INTEGER              :: i,j,k
! ----------------------------------------------------------------------
! Satellite physical informaiton
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: X_SIDE,Z_SIDE
      REAL (KIND = prec_q) :: MASS,AREA
      REAL (KIND = prec_q) :: A_SOLAR
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: u_sat,i_sat,omega_sat

! ----------------------------------------------------------------------
! Sun-related variables
! ----------------------------------------------------------------------
       REAL (KIND = prec_q) :: u_sun,beta,del_u
       REAL (KIND = prec_q), DIMENSION(3) :: r_sun1,r_sun2

! ----------------------------------------------------------------------
! Numerical Constants
      Ps = 4.56D-6 ! (Nm^-2)
      AU = 1.4959787066d11 ! (m)
      Pi = 4*atan(1.0d0)
! ---------------------------------------------------------------------



! The unit vector er SAT->EARTH
      er(1)=-r(1)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(2)=-r(2)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(3)=-r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2)

! The unit vector ed SAT->SUN (where is just opposite to the solar radiation vector)
      Ds=sqrt((r_sun(1)-r(1))**2+(r_sun(2)-r(2))**2+(r_sun(3)-r(3))**2)
      ed(1)=((r_sun(1)-r(1))/Ds)
      ed(2)=((r_sun(2)-r(2))/Ds)
      ed(3)=((r_sun(3)-r(3))/Ds)

! The unit vector ey = er x ed, parallel to the rotation axis of solar panel
      CALL cross_product (er,ed,ey)
      ey(1)=ey(1)/sqrt(ey(1)**2+ey(2)**2+ey(3)**2)
      ey(2)=ey(2)/sqrt(ey(1)**2+ey(2)**2+ey(3)**2)
      ey(3)=ey(3)/sqrt(ey(1)**2+ey(2)**2+ey(3)**2)

! The unit vector eb = ed x ey
      CALL cross_product (ed,ey,eb)
      eb(1)=eb(1)/sqrt(eb(1)**2+eb(2)**2+eb(3)**2)
      eb(2)=eb(2)/sqrt(eb(1)**2+eb(2)**2+eb(3)**2)
      eb(3)=eb(3)/sqrt(eb(1)**2+eb(2)**2+eb(3)**2)
!print *, eb(1),eb(2), eb(3)


! computation of the satellite argument of latitude and orbit inclination
      CALL kepler_z2k (r,v,GM,kepler)
      u_sat = kepler(9)*Pi/180.d0
      i_sat = kepler(3)*Pi/180.d0
      omega_sat = kepler(4)*Pi/180.d0

! compute the sun position in the satellite orbit plane by rotating big Omega_sat and i_sat,
! allowing us for the computation of u_sun and sun elevation angles (beta)
!   sun                 sun
!  r    = R(big_w,i) x r
!   kep                 ICF
!
      do i=1,3
         do j=1,3
          R33(i,j)=1.0d0
         end do
      end do
      R33(1,3)=0.0d0
      R33(2,3)=0.0d0
      R33(3,1)=0.0d0
      R33(3,2)=0.0d0


      do i=1,3
         do j=1,3
          R11(i,j)=1.0d0
         end do
      end do
      R11(1,2)=0.0d0
      R11(1,3)=0.0d0
      R11(2,1)=0.0d0
      R11(3,1)=0.0d0


! rotate big omega to make the X-axis of the celestial frame consistent with the
! direction of right ascension of ascending node
      CALL R3(omega_sat,R33)

! rotate inclination to make the XY plane of the celestial frame consistent with
! the orbit plane
      CALL R1(i_sat,R11)

! convert the sun position from the celestial frame to the orbit plane
      CALL matrix_Rr(R33,r_sun,r_sun1)
      CALL matrix_Rr(R11,r_sun1,r_sun2)

!compute the sun argument of lattitude and beta angle
      u_sun = atan2(r_sun2(2),r_sun2(1)) ! in rad
!      beta  = atan2(r_sun2(3),sqrt(r_sun2(1)**2+r_sun2(2)**2)) ! in rad

      del_u = u_sat - u_sun ! in rad
      if (del_u*180/Pi .gt.360.0d0) then
      del_u=del_u-2*Pi
      else if(del_u*180/Pi .lt.0.0d0) then
      del_u=del_u+2*Pi
      end if


!write (*,*) beta*180.0d0/Pi, del_u*180.0d0/Pi


     if (srpid .eq. 2) then

! Box-wing model
! ================================

! GPS Block IIR types
      if (prnnum.ge.11 .and. prnnum.le.23 .or. prnnum.eq.2 .or. &
       prnnum.eq.5 .or. prnnum.eq.7 .or. prnnum.eq.28 .or. &
       prnnum.eq.29 .or. prnnum.eq.31) then
      MASS   = 1100.0D0
! Block IIF types
      else
      MASS   = 1555.3D0
      end if

! The unit vector of x side (ex = ey x er), which is always illuminated by the
! sun, parallel to the satellite flight/velocity direction

!     CALL cross_product (ey,er,ex)

! The unit vector ex 
      ex(1)=v(1)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ex(2)=v(2)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ex(3)=v(3)/sqrt(v(1)**2+v(2)**2+v(3)**2)


! angles between ed and each surface(k)
! k=1: +X
!   2: +Y
!   3: +Z
!   4: solar panels

     cosang(1)=ed(1)*ex(1)+ed(2)*ex(2)+ed(3)*ex(3)
     cosang(2)=ed(1)*ey(1)+ed(2)*ey(2)+ed(3)*ey(3)
     cosang(3)=ed(1)*er(1)+ed(2)*er(2)+ed(3)*er(3)
     cosang(4)=ed(1)*ed(1)+ed(2)*ed(2)+ed(3)*ed(3)

!write(*,*) acos(cosang(1))*180/pi, del_u*180.0d0/Pi

      CALL surfprop(BLKNUM,AREA1,REFL1,DIFU1,ABSP1)

!      write(*,*)'BLK=', BLKNUM, AREA1,REFL1,DIFU1,ABSP1

      CALL productdot(r_sun,r,ab)
      Lab = sqrt(r(1)**2+r(2)**2+r(3)**2)*sqrt(r_sun(1)**2+r_sun(2)**2+r_sun(3)**2)
      E   = Pi - acos(ab/Lab)
      Edot = -(0.01*(Pi/180)*sin(beta)*cos(del_u)+0.09*cos(beta)*sin(del_u))/sin(E)
 
!
! partial derivative outputs
!----------------------------
! Asrp(i,k), k=1-9
! k = 1 : (1+uv+2/3*(v-uv)) for solar panels
!     2 : (1-uv)            for the absorption and diffusion of +X
!     3 : uv                for the reflection               of +X  
!     4 : (1-uv)            for the absorption and diffusion of +Z
!     5 : uv                for the reflection               of +Z
!     6 : (1-uv)            for the absorption and diffusion of -Z
!     7 : uv                for the reflection               of -Z
!     8 :                   for the Y bias
!     9 :                   for the solar panel rotation lag
!----------------------------------------------------------------------

      do i=1,3
      Asrp(i,1)=-AREA1(4)*(Ps/MASS)*ed(i)
! Judge the positive surface or negative surface facing to the Sun
      IF(cosang(1).GE.0D0)THEN
      Asrp(i,2)=-AREA1(1)*(Ps/MASS)*abs(cosang(1))*(ed(i)+(2/3)*ex(i))
      Asrp(i,3)=-AREA1(1)*(Ps/MASS)*2*cosang(1)**2*ex(i)
       ELSEIF(cosang(1).LT.0D0)THEN
      Asrp(i,2)=-AREA1(1)*(Ps/MASS)*abs(cosang(1))*(ed(i)+(2/3)*(-ex(i)))
      Asrp(i,3)=-AREA1(1)*(Ps/MASS)*2*cosang(1)**2*(-ex(i))
      END IF

       IF(cosang(3).GE.0D0)THEN
      Asrp(i,4)=-AREA1(3)*(Ps/MASS)*abs(cosang(3))*(ed(i)+(2/3)*er(i))
      Asrp(i,5)=-AREA1(3)*(Ps/MASS)*2*cosang(3)**2*er(i)
      Asrp(i,6)= 0.0d0
      Asrp(i,7)= 0.0d0
       ELSEIF(cosang(3).LT.0D0)THEN
      Asrp(i,4)= 0.0d0
      Asrp(i,5)= 0.0d0
      Asrp(i,6)=-AREA1(3)*(Ps/MASS)*abs(cosang(3))*(ed(i)+(2/3)*(-er(i)))
      Asrp(i,7)=-AREA1(3)*(Ps/MASS)*2*cosang(3)**2*(-er(i))
       END IF

      Asrp(i,8)=-1.0d0*ey(i)
      if (Edot.lt.0.0d0) then
      Asrp(i,9)=AREA1(4)*(Ps/MASS)*2*((DIFU1(4)/3)+REFL1(4))*eb(i)
      else
      Asrp(i,9)=-AREA1(4)*(Ps/MASS)*2*((DIFU1(4)/3)+REFL1(4))*eb(i)
      end if

     if (abs(beta) .lt. 14.0d0 .and.  abs(del_u*180.0d0/Pi-180) .lt. 14.0d0)then

      Asrp(i,1)=0.0d0
      Asrp(i,2)=0.0d0
      Asrp(i,3)=0.0d0
      Asrp(i,4)=0.0d0
      Asrp(i,5)=0.0d0
      Asrp(i,6)=0.0d0
      Asrp(i,7)=0.0d0
      Asrp(i,8)=0.0d0
      Asrp(i,9)=0.0d0

     end if

      end do


! The following is used to output the partial derivative results of the box-wing
! model. 
! -------------------------------------------------------------------------------
!        i=1

!        Asrp(i,1)=-AREA1(4)*(Ps/MASS)
!        IF(cosang(1).GE.0D0)THEN
!        Asrp(i,2)=-AREA1(1)*(Ps/MASS)*(cosang(1))
!        Asrp(i,3)=-AREA1(1)*(Ps/MASS)*2*cosang(1)**2

!        ELSEIF(cosang(1).LT.0D0)THEN
!        Asrp(i,2)= AREA1(1)*(Ps/MASS)*abs(cosang(1))
!        Asrp(i,3)= AREA1(1)*(Ps/MASS)*2*cosang(1)**2

!        ENDIF

!         IF(cosang(3).GE.0D0)THEN
!         Asrp(i,4)=-AREA1(3)*(Ps/MASS)*(cosang(3))
!         Asrp(i,5)=-AREA1(3)*(Ps/MASS)*2*cosang(3)**2
!         Asrp(i,6)=0.0d0
!         Asrp(i,7)=0.0d0
!         ELSEIF(cosang(3).LT.0D0)THEN
!         Asrp(i,4)=0.0d0
!         Asrp(i,5)=0.0d0
!         Asrp(i,6)= AREA1(3)*(Ps/MASS)*abs(cosang(3))
!         Asrp(i,7)= AREA1(3)*(Ps/MASS)*2*cosang(3)**2
!         END IF


!      Asrp(i,8)=-1.0d0
!      if (Edot.lt.0.0d0) then
!      Asrp(i,9)=AREA1(4)*(Ps/MASS)*2*((DIFU1(4)/3)+REFL1(4))
!      else
!      Asrp(i,9)=-AREA1(4)*(Ps/MASS)*2*((DIFU1(4)/3)+REFL1(4))
!      end if

!     if (abs(beta) .lt. 14.0d0 .and.  abs(del_u*180.0d0/Pi-180) .lt. 14.0d0)then

!      Asrp(1,1)=0.0d0
!      Asrp(1,2)=0.0d0 
!      Asrp(1,3)=0.0d0
!      Asrp(1,4)=0.0d0
!      Asrp(1,5)=0.0d0
!      Asrp(1,6)=0.0d0
!      Asrp(1,7)=0.0d0
!      Asrp(1,8)=0.0d0
!      Asrp(1,9)=0.0d0

!     end if
!
!---------------------------------------------------------------------------------------------

!      write(*,*)'Asrp', Asrp 
! -------------------------
! end of the box-wing model
! =======================================================================


     else if (srpid .eq. 3) then
! ECOM (D2B1) model
! ***************************************
! a=D(u)*ed+Y(u)*ey+B(u)*eb
!
! D(u)=D0+D2c*cos2*(del_u)+D2s*sin2*(del_u)
!        +D4c*cos4*(del_u)+D4s*sin4*(del_u)
!
! Y(u)=Y0
!
! B(u)=B0+B1c*cos1*(del_u)+B1s*sin1*(del_u)
! ***************************************
! da/dD0  = 1*ed
! da/dY0  = 1*ey
! da/dB0  = 1*eb
! da/dD2c = cos2u*ed
! da/dD4c = cos4u*ed
! da/dD2s = sin2u*ed
! da/dD4s = sin4u*ed
! da/dB1c = cos1u*eb 
! da/dB1s = sin1u*eb
!=========================================

! A scaling factor is applied to ECOM model as suggested by Bernese
!******************************************************************

      sclfa=(AU/Ds)**2

      do i=1,3
      
      Asrp(i,1)=-sclfa*1.0d0*ed(i)
      Asrp(i,2)=-sclfa*1.0d0*ey(i)
      Asrp(i,3)=-sclfa*1.0d0*eb(i)
      Asrp(i,4)=-sclfa*DCOS(2*del_u)*ed(k)
      Asrp(i,5)=-sclfa*DCOS(4*del_u)*ed(k)
      Asrp(i,6)=-sclfa*DSIN(2*del_u)*ed(k)
      Asrp(i,7)=-sclfa*DSIN(4*del_u)*ed(k)
      Asrp(i,8)=-sclfa*DCOS(1*del_u)*eb(k)
      Asrp(i,9)=-sclfa*DSIN(1*del_u)*eb(k)

      
      end do

     if (abs(beta) .lt. 14.0d0 .and.  abs(del_u*180.0d0/Pi-180) .lt. 14.0d0)then

      Asrp(i,1)=0.0d0
      Asrp(i,2)=0.0d0
      Asrp(i,3)=0.0d0
      Asrp(i,4)=0.0d0
      Asrp(i,5)=0.0d0
      Asrp(i,6)=0.0d0
      Asrp(i,7)=0.0d0
      Asrp(i,8)=0.0d0
      Asrp(i,9)=0.0d0

     end if



     end if

! end of ECOM (D2B1) model
! ==================================================================

END
