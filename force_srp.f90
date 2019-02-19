
SUBROUTINE force_srp (GM,prnnum,satsvn,eclpf,srpid,r,v,r_sun,fx,fy,fz )


! ----------------------------------------------------------------------
! SUBROUTINE: force_srp.f90
! ----------------------------------------------------------------------
! Purpose:
! Acceleration due to the solar radiation pressure 
! Computation of SRP acceleration using various SRP models, 
! such as a simply cannonball model (srpid=1),a box-wing model (srpid=2) 
! and a ECOM model (srpid=3) 
! ----------------------------------------------------------------------
! Input arguments:
! - GM           : the earth gravitational constant  
! - prnnum       : satellite PRN number
! - satsvn       : satellite SVN
! - eclpf        : =0: in the sun illumination; 
!                  =1: in the satellite eclipse
! - srpid        : =1: a simply cannonball model; 
!                  =2: box-wing model;
!                  =3: ECOM model
! - r            : satellite position vector (m)
! - v            : satellite velocity vector
! - r_sun        : Sun position vector
! 
! Output arguments:
! - fx,fy,fz:		Acceleration's cartesian components (m)
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	01-10-2017
! 
! Changes:      10-10-2018 Tzupang Tseng: modify the Cannonball model to be 
!                                            similar to the box-wing model in
!                                            preparation for BDS with the ON
!                                            attitude mode.
!               11-12-2018 Tzupang Tseng: make the force matrix dynamic
!               23-01-2019 Tzupang Tseng: add the ECOM2 model
!               31-01-2019 Tzupang Tseng: change the definition of ey by
!                                         dividing the length of ey
!
! Copyright:  GEOSCIENCE AUSTRALIA, AUSTRALIA
! ----------------------------------------------------------------------

      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_satinfo
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      INTEGER                           :: eclpf,srpid,ECOM
      INTEGER                           :: prnnum,BLKNUM
      REAL (KIND = prec_q), DIMENSION(3) :: r,v,r_sun
      REAL (KIND = prec_q)               :: fx,fy,fz
      INTEGER                            :: satsvn

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: Ps,AU,Pi
      REAL (KIND = prec_q) :: Cr
      REAL (KIND = prec_q) :: ANG,GM
      REAL (KIND = prec_q) :: Ds,sclfa

      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: surforce(4,3)
      REAL (KIND = prec_q), DIMENSION(3) :: er,ed,ey,eb,ex,en,ev,ez
      REAL (KIND = prec_q), DIMENSION(3) :: yy
      REAL (KIND = prec_q), DIMENSION(3) :: fsrp
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(:), ALLOCATABLE :: srpcoef 
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q), DIMENSION(4) :: AREA1,REFL1,DIFU1,ABSP1
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
      INTEGER              :: zta
      INTEGER              :: i,j,k,m
      INTEGER              :: N_param, PD_Param_ID
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
      Cr = 1.4 ! SRP coefficient ranges from 1.3 to 1.5
      AU = 1.496d11 ! (m)
      Pi = 4*atan(1.0d0)
! ---------------------------------------------------------------------

! initialize the SRP force
! -------------------------
     DO i=1,3
     fsrp(i)=0.0d0
     END DO

! GPS constellation
! -----------------
      if(prnnum.le.100)then
! IIF
         if(satsvn.ge.62.and.satsvn.le.73) then
         MASS   = 1633.0d0
         Z_SIDE = 5.05D0
         X_SIDE = 4.55D0
         A_SOLAR= 22.25D0
! IIR
         else
         MASS   = 1080.0d0
         Z_SIDE = 4.25D0
         X_SIDE = 4.11D0
         A_SOLAR= 13.92D0

         end if

! GLONASS constellation
! ---------------------
      else if (prnnum .gt. 100 .and. prnnum .le. 200) then
         Z_SIDE = 1.6620D0
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-K
         if(satsvn.eq.801.or.satsvn.eq.802.or.satsvn.eq.855)then
         MASS = 995.0d0
! GLONASS-M
         else
         MASS   = 1415.0d0
         end if

! GALILEO constellation
! ---------------------
      else if (prnnum .gt. 200 .and. prnnum .le. 300) then
         MASS   = 700.0d0
         Z_SIDE = 3.002D0
         X_SIDE = 1.323D0
         A_SOLAR= 11.0D0

! BDS constellation
! -----------------
      else if (prnnum .gt. 300 .and. prnnum .le. 400) then
         Z_SIDE = 3.96D0
         X_SIDE = 4.5D0
         A_SOLAR= 22.44D0

! BDS MEO
         if(satsvn.ge.12.and.satsvn.le.15)then
         MASS   = 800.0d0
! BDS IGSO
         elseif(satsvn.ge.7.and.satsvn.le.10.or.satsvn.eq.5.or.satsvn.eq.17)then
         MASS = 1400.0d0
         end if

! QZSS constellation
! ------------------
      else if (prnnum .gt. 400 .and. prnnum .le. 500) then
!         if(satsvn.eq.1)then
         MASS = 2000.0d0
         Z_SIDE = 6.00D0
         X_SIDE = 12.2D0
         A_SOLAR= 40.0D0
!         end if
      end if


! The unit vector ez SAT->EARTH
      er(1)=r(1)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(2)=r(2)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(3)=r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      ez(1)=-er(1)
      ez(2)=-er(2)
      ez(3)=-er(3)

! The unit vector ed SAT->SUN (where is just opposite to the solar radiation vector)
      Ds=sqrt((r_sun(1)-r(1))**2+(r_sun(2)-r(2))**2+(r_sun(3)-r(3))**2)
      ed(1)=((r_sun(1)-r(1))/Ds)
      ed(2)=((r_sun(2)-r(2))/Ds)
      ed(3)=((r_sun(3)-r(3))/Ds)
! The unit vector ey = ez x ed/|ez x ed|, parallel to the rotation axis of solar panel
      CALL cross_product (ez,ed,yy)
      ey(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
! The unit vector eb = ed x ey
      CALL cross_product (ed,ey,eb)

! The unit vector of x side, which is always illuminated by the
! sun.

      CALL cross_product (ey,ez,ex)

! the orbit normal vector
!------------------------

      ev(1)=v(1)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ev(2)=v(2)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ev(3)=v(3)/sqrt(v(1)**2+v(2)**2+v(3)**2)

     Call cross_product(er,ev,en)

! computation of the factor zta related to the Sun illumination
      if (eclpf.eq.0) then
         zta=1 ! in the sun light area
      else
         zta=0 ! in the eclipse
         fx=0.0d0
         fy=0.0d0
         fz=0.0d0
      return
      end if

! computation of the satellite argument of latitude and orbit inclination
      CALL kepler_z2k (r,v,GM,kepler)
      u_sat = kepler(9)*Pi/180.d0
      i_sat = kepler(3)*Pi/180.d0
      omega_sat = kepler(4)*Pi/180.d0

! compute the sun position in the satellite orbit plane by rotating omega_sat and i_sat,
! allowing us for the computation of u_sun and sun elevation angles (beta)
      do i=1,3
         do j=1,3
          R33(i,j)=1.0d0
         end do
      end do
      R33(1,3)=0.0d0
      R33(2,3)=0.0d0
      R33(3,1)=0.0d0
      R33(3,2)=0.0d0
      R33(3,3)=1.0d0

      do i=1,3
         do j=1,3
          R11(i,j)=1.0d0
         end do
      end do
      R11(1,1)=1.0d0
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

! convert the sun position in the celestial frame to the orbit plane
      CALL matrix_Rr(R33,r_sun,r_sun1)
      CALL matrix_Rr(R11,r_sun1,r_sun2)

      beta  = atan2(r_sun2(3),sqrt(r_sun2(1)**2+r_sun2(2)**2)) ! in rad
!write (*,*) beta*180.0d0/Pi

      u_sun = atan2(r_sun2(2),r_sun2(1)) ! in rad

      del_u = u_sat - u_sun

      if (del_u*180/Pi .gt.360.0d0) then
      del_u=del_u-2*Pi
      else if(del_u*180/Pi .lt.0.0d0) then
      del_u=del_u+2*Pi
      end if
!print*,'del_u=',del_u*180/Pi

!========================================
! angles between ed and each surface(k)
! k=1: +X
!   2: +Y
!   3: +Z
!   4: solar panels
!=======================================
     cosang(1)=ed(1)*ex(1)+ed(2)*ex(2)+ed(3)*ex(3)
     cosang(2)=ed(1)*ey(1)+ed(2)*ey(2)+ed(3)*ey(3)
     cosang(3)=ed(1)*ez(1)+ed(2)*ez(2)+ed(3)*ez(3)
     cosang(4)=ed(1)*ed(1)+ed(2)*ed(2)+ed(3)*ed(3)

!print*, 'COSANG1=', cosang(1)

      IF (srpid == 1) THEN
! A simply cannonball model
! *********************************
!  a=-zta*Cr*(A/m)P(1AU/r)^2*Dr(i)
! ********************************

! The main surface area face toward to the Sun using the SAT->SUN and SAT->EARTH
! vectors
     ANG=acos(ed(1)*ez(1)+ed(2)*ez(2)+ed(3)*ez(3))*180.0d0/Pi

     if (abs(ANG) .lt. 30.0d0) then
     AREA=Z_SIDE+A_SOLAR
     else
     AREA=X_SIDE+A_SOLAR
     end if
! Cartesian counterparts (fx,fy,fz) of acceleration fr
      fx = -zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(1)
      fy = -zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(2)
      fz = -zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(3) 

!      fx = -zta*Cr/MASS*Ps*(AU/Ds)**2*ed(1)*(X_SIDE*cosang(1)+Z_SIDE*cosang(3)+A_SOLAR*cosang(4))
!      fy = -zta*Cr/MASS*Ps*(AU/Ds)**2*ed(2)*(X_SIDE*cosang(1)+Z_SIDE*cosang(3)+A_SOLAR*cosang(4))
!      fz = -zta*Cr/MASS*Ps*(AU/Ds)**2*ed(3)*(X_SIDE*cosang(1)+Z_SIDE*cosang(3)+A_SOLAR*cosang(4))

!      if (abs(ANG) .le. 14 ) then
!         fx=0.0d0
!         fy=0.0d0
!         fz=0.0d0
!      end if

! end of the cannonball model
!--------------------------------------------------------------------------


! Box-wing model
! **************
     ELSE IF (srpid == 2) THEN

!      CALL surfprop(BLKNUM,AREA1,REFL1,DIFU1,ABSP1)

! Forces caused by different interactions between optical properties and the
! satellite surface

!write(*,*)Ps/MASS
     do k=1,4
     if (k .eq. 1) then
! Judge the positive surface or negative surface facing to the Sun
        IF(cosang(k).GE.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*ex(i))+2*cosang(k)**2*REFL1(k)*ex(i)
      end do
        ELSEIF(cosang(k).LT.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*(-ex(i)))+2*cosang(k)**2*REFL1(k)*(-ex(i))
      end do
        END IF

     elseif (k .eq. 2) then
        IF(cosang(k).GE.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*ey(i))+2*cosang(k)**2*REFL1(k)*ey(i)
      end do
        ELSEIF(cosang(k).LT.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*(-ey(i)))+2*cosang(k)**2*REFL1(k)*(-ey(i))
      end do
        END IF

     elseif (k .eq. 3) then
        IF(cosang(k).GE.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*ez(i))+2*cosang(k)**2*REFL1(k)*ez(i)
      end do
        ELSEIF(cosang(k).LT.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*(-ez(i)))+2*cosang(k)**2*REFL1(k)*(-ez(i))
      end do
        END IF

     elseif (k .eq. 4) then
      do i=1,3
! Here the normal of the soalr panels is assumed to be parallel to ed
      surforce(k,i)=AREA1(k)*ed(i)+REFL1(k)*ed(i)+(2/3)*DIFU1(k)*ed(i)
      end do


     end if
     end do

! Compute the total forces
!
     do i=1,3
     fsrp(i)=0.0d0
     end do

     do i=1,3
        do k=1,4
     fsrp(i)=fsrp(i)+(Ps/MASS)*surforce(k,i)
        end do
     end do

! forces in inertial frame
!-------------------------
     fx=-fsrp(1)
     fy=-fsrp(2)
     fz=-fsrp(3)

! end of the box-wing model
! =======================================================================


     ELSE IF (srpid == 3) THEN
! A scaling factor 
!*************************************************************************
     sclfa=(AU/Ds)**2

ALLOCATE (srpcoef(NPARAM_glb), STAT = AllocateStatus)

! ECOM1 model
! ***********************************************************************
      IF (ECOM_param_glb == 1 ) then
      PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ed(i)
        END DO
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i)
        END DO
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i)
        END DO
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*ed(i)
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*ed(i)
        END DO
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*ey(i)
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*ey(i)
        END DO
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*eb(i)
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*eb(i)
        END DO

End If

! ECOM2 model
! **********************************************************************

     ELSE IF (ECOM_param_glb == 2 ) then
      PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ed(i)
        END DO
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i)
        END DO
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i)
        END DO
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(2*del_u)*ed(i)
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(2*del_u)*ed(i)
        END DO
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(4*del_u)*ed(i)
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(4*del_u)*ed(i)
        END DO
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*eb(i)
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*eb(i)
        END DO

End If
      END IF 


     fx=-fsrp(1)
     fy=-fsrp(2)
     fz=-fsrp(3)

! Implement BW + ECOM model (To be worked on)
! ***********************************************************************
! IF (SRP_MOD == 4) THEN
!      fxo =-Cr/MASS*Ps*sclfa*ed(1)*(X_SIDE*cosang(1)+Z_SIDE*cosang(3)+A_SOLAR*cosang(4))
!      fyo =-Cr/MASS*Ps*sclfa*ed(2)*(X_SIDE*cosang(1)+Z_SIDE*cosang(3)+A_SOLAR*cosang(4))
!      fzo =-Cr/MASS*Ps*sclfa*ed(3)*(X_SIDE*cosang(1)+Z_SIDE*cosang(3)+A_SOLAR*cosang(4))
!     fx=-(fsrp(1) + fxo)
!     fy=-(fsrp(2) + fyo)
!     fz=-(fsrp(3) + fzo)

! END IF
     END IF 

! end of ECOM-based model
! ---------------------------------------------------------

END
