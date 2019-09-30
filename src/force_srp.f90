
SUBROUTINE force_srp (lambda, eBX_ecl, GM, prnnum,  srpid, r, v, r_sun, fx, fy, fz)


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
! - srpid        : =1: a simply cannonball model; 
!                  =2: box-wing model;
!                  =3: ECOM model
! - r            : satellite position vector (m)
! - v            : satellite velocity vector
! - r_sun        : Sun position vector
! - lambda       : shadow coefficient
! - eBX_ecl      : dynamic ex of satellite body frame
! 
! Output arguments:
! - fx,fy,fz:		Acceleration's cartesian components (m)
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	01-10-2017
! 
! Changes:      10-10-2018 Tzupang Tseng: modify the Cannonball model to be 
!                                         similar to the box-wing model in
!                                         preparation for BDS with the ON
!                                         attitude mode.
!               11-12-2018 Tzupang Tseng: make the force matrix dynamic
!               23-01-2019 Tzupang Tseng: add the ECOM2 model
!               31-01-2019 Tzupang Tseng: change the definition of ey by
!                                         dividing the length of ey
!               20-02-2019 Tzupang Tseng: create a functionality for switching
!                                         on and off some particular coefficients in ECOM models
!               22-03-2019 Tzupang Tseng: set a simple condition for the eclipsed satellites
!                                         where only the D0 accelerations are setup to zero 
!               02-05-2019 Tzupang Tseng: use the coefficient from the shadow.f90 for scaling the SRP effect
!               30-08-2019 Tzupang Tseng: implement the orbit-normal attitude for the BDS SRP estimation (still need to be tested)
!               03-09-2019 Tzupang Tseng: use a simple box-wing model as a  priori SRP value where the ECOM is 
!                                         used to adjust the box-wing model (BOX-WING + ECOM)
!
! Copyright:  GEOSCIENCE AUSTRALIA, AUSTRALIA
! ----------------------------------------------------------------------

      USE mdl_precision
      USE mdl_num
      USE mdl_param
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q), INTENT (IN) :: lambda
      REAL (KIND = prec_d) , Dimension(3), INTENT(IN) :: eBX_ecl
      INTEGER                           :: srpid,ECOM
      INTEGER                           :: prnnum
      REAL (KIND = prec_q), DIMENSION(3) :: r,v,r_sun
      REAL (KIND = prec_q)               :: fx,fy,fz

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: Ps,AU,Pi
      REAL (KIND = prec_q) :: Cr
      REAL (KIND = prec_q) :: ANG,GM
      REAL (KIND = prec_q) :: Ds,sclfa
      REAL (KIND = prec_q) :: fxo,fyo,fzo

      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: surforce(4,3)
      REAL (KIND = prec_q), DIMENSION(3) :: er,ed,ey,eb,ex,en,ev,ez,ECOM_ey
      REAL (KIND = prec_q), DIMENSION(3) :: yy
      REAL (KIND = prec_q), DIMENSION(3) :: fsrp
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(:), ALLOCATABLE :: srpcoef 
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q), DIMENSION(4) :: AREA1,REFL1,DIFU1,ABSP1
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
      INTEGER              :: zta
      INTEGER              :: ex_i
      INTEGER              :: i,j,k,m
      INTEGER              :: N_param, PD_Param_ID
      INTEGER              :: att_ON
      INTEGER              :: flag_BW
! ----------------------------------------------------------------------
! Satellite physical informaiton
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: X_SIDE,Z_SIDE
      REAL (KIND = prec_q) :: AREA
      REAL (KIND = prec_q) :: A_SOLAR
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: u_sat,i_sat,omega_sat
      REAL (KIND = 8)      :: II, KN, U
      REAL (KIND = prec_q) :: F0,alpha
! ----------------------------------------------------------------------
! Sun-related variables
! ----------------------------------------------------------------------
       REAL (KIND = prec_q) :: u_sun
       REAL (KIND = prec_q) :: beta,del_u
       REAL (KIND = prec_q), DIMENSION(3) :: r_sun1,r_sun2
! ----------------------------------------------------------------------
       INTEGER*4 BLKNUM,SVN,REFF,ERM,ANT,GRD,MONTH
       REAL*8  ACCEL(3),SUN(3)
       REAL*8  YSAT(6)
! ---------------------------------------------------------------------- 
! Numerical Constants
      Ps = 4.5567D-6 ! (Nm^-2)
      Cr = 1.5 ! SRP coefficient ranges from 1.3 to 1.5
      AU = 1.496d11 ! (m)
      Pi = 4*atan(1.0d0)
    ex_i = 0 ! change the definition of the unit vector ex 
             ! ex_i = 0 (default)
             !      = 1 (using dynamic ex vector from attitude routine)
  att_ON = 0 ! att_ON = 1 : use the orbit-normal attitude for BDS satellite
             !              when the beta < 4 deg
             !        = 0 : use the yaw-steering attitude for BDS satellite 
             !              for all beta angles 
 flag_BW = 0 !flag_BW = 1 : use the simple box-wing model as a priori SRP value
             !        = 2 : use the box-wing model from repro3 routines as a
             !              priori SRP value
             !        = 0 : use the constant f0 as a priori SRP value
             !        = any numbers : directly estimate the SRP parameters
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
         if(SVNID.ge.62.and.SVNID.le.73) then
         Z_SIDE = 5.05D0
         X_SIDE = 4.55D0
         A_SOLAR= 22.25D0
         F0 = 16.7d-5
! IIR
         else
         Z_SIDE = 4.25D0
         X_SIDE = 4.11D0
         A_SOLAR= 13.92D0
         F0 = 11.15d-5
         end if

! GLONASS constellation
! ---------------------
      else if (prnnum .gt. 100 .and. prnnum .le. 200) then
         Z_SIDE = 1.6620D0
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0
! GLONASS-K
         if(SVNID.eq.801.or.SVNID.eq.802.or.SVNID.eq.855)then
         F0 = 10.0d-5
! GLONASS-M
         else
         F0 = 20.9d-5
         end if

! GALILEO constellation
! ---------------------
      else if (prnnum .gt. 200 .and. prnnum .le. 300) then
         Z_SIDE = 3.002D0
         X_SIDE = 1.323D0
         A_SOLAR= 11.0D0
         F0 = 8.35d-5
! BDS constellation
! -----------------
      else if (prnnum .gt. 300 .and. prnnum .le. 400) then
         Z_SIDE = 3.96D0
         X_SIDE = 4.5D0
         A_SOLAR= 22.44D0
! BDS MEO
         if(SVNID.ge.12.and.SVNID.le.15)then
         F0 = 8.35d-5
! BDS IGSO
         elseif(SVNID.ge.7.and.SVNID.le.10.or.SVNID.eq.5.or.SVNID.eq.17)then
         F0 = 50.1d-5
         end if
! QZSS constellation
! ------------------
      else if (prnnum .gt. 400 .and. prnnum .le. 500) then
         if(SVNID.eq.1)then
         Z_SIDE = 6.00D0
         X_SIDE = 12.2D0
         A_SOLAR= 40.0D0
         F0 = 50.1d-5 ! Assumed to be the same with BDS/IGSO
         end if
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

      CALL productcross (ez,ed,yy)
      ey(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

! The unit vector of x side, which is always illuminated by the sun.

      CALL productcross (ey,ez,ex)

! Using the BODY-X univector from Kouba routine to redefine the ey for the satellite eclipsed
IF(ex_i .gt. 0) THEN

      ex(1)=eBX_ecl(1)/sqrt(eBX_ecl(1)**2+eBX_ecl(2)**2+eBX_ecl(3)**2)
      ex(2)=eBX_ecl(2)/sqrt(eBX_ecl(1)**2+eBX_ecl(2)**2+eBX_ecl(3)**2)
      ex(3)=eBX_ecl(3)/sqrt(eBX_ecl(1)**2+eBX_ecl(2)**2+eBX_ecl(3)**2)

      CALL productcross (ez,ex,ey)
END IF
! The unit vector eb = ed x ey

      CALL productcross (ed,ey,eb)

! the orbit normal vector
!------------------------

      ev(1)=v(1)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ev(2)=v(2)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ev(3)=v(3)/sqrt(v(1)**2+v(2)**2+v(3)**2)

     CALL productcross (er,ev,en)

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

!print*,'del_u=, lambda=',del_u*180/Pi, lambda

! Implement the orbit-normal attitude for BDS satellites when the beat < 4 deg
! ----------------------------------------------------------------------------
     if (att_ON == 1) then
     if(prnnum .gt. 300 .and. prnnum .le. 400) then
        if (abs(beta*180.0d0/Pi) < 4.d0) then
!        PRINT*,'The orbit-normal attitude is applied.'
!        PRINT*,'ed_YS =',ed, sqrt(ed(1)**2+ed(2)**2+ed(3)**2)
!        print*,'ey_YS =',ey, sqrt(ey(1)**2+ey(2)**2+ey(3)**2)
!        print*,'eb_YS =',eb, sqrt(eb(1)**2+eb(2)**2+eb(3)**2)
        CALL productcross (ez,ev,yy)
        ey(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        ey(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        ey(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

        CALL productcross (ed,ey,yy)
        eb(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        eb(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        eb(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        
        CALL productcross (ey,eb,yy)
        ed(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        ed(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
        ed(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
!        PRINT*,'ed_ON =',ed, sqrt(ed(1)**2+ed(2)**2+ed(3)**2)
!        print*,'ey_ON =',ey, sqrt(ey(1)**2+ey(2)**2+ey(3)**2)
!        print*,'eb_ON =',eb, sqrt(eb(1)**2+eb(2)**2+eb(3)**2)
!
        end if
     end if
     end if

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

! A scaling factor is applied to ECOM model
!******************************************************************
      sclfa=(AU/Ds)**2

! SIMPLE BOX-WING 
      if (flag_BW == 1 .or. srpid == 1) then
         fxo=sclfa*Ps/MASS*(X_SIDE*cosang(1)*ex(1)+Z_SIDE*cosang(3)*ez(1)+1*A_SOLAR*cosang(4)*ed(1))
         fyo=sclfa*Ps/MASS*(X_SIDE*cosang(1)*ex(2)+Z_SIDE*cosang(3)*ez(2)+1*A_SOLAR*cosang(4)*ed(2))
         fzo=sclfa*Ps/MASS*(X_SIDE*cosang(1)*ex(3)+Z_SIDE*cosang(3)*ez(3)+1*A_SOLAR*cosang(4)*ed(3))
         alpha = sqrt(fxo**2+fyo**2+fzo**2)

! BOX-WING model from the repro3 routine
! --------------------------------------
      else if (flag_BW == 2 .or. srpid == 2) then
         REFF = 0     
         YSAT(1:3) = r
         YSAT(4:6) = v
         CALL SRPFBOXW(REFF,YSAT,R_SUN,BLKID,SVNID,ACCEL)
         alpha = sclfa*sqrt(ACCEL(1)**2+ACCEL(2)**2+ACCEL(3)**2)
      
      else if (flag_BW == 0) then
         alpha = F0/MASS
         alpha = F0/MASS
      else
         alpha = 1.d0
      end if

     IF (srpid == 3) THEN 

ALLOCATE (srpcoef(NPARAM_glb), STAT = AllocateStatus)

! ECOM1 model
! ***********************************************************************
      IF (ECOM_param_glb == 1 ) then
      PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ed(i)*alpha
        END DO
!print*,'ECOM1-caused accelerations'
!print*,'D0'
Else 
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i)*alpha
        END DO
!print*,'Y0'
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i)*alpha
        END DO
!print*,'B0'
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*ed(i)*alpha
        END DO
!print*,'DC'
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*ed(i)*alpha
        END DO
!print*,'DS'
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*ey(i)*alpha
        END DO
!print*,'YC'
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*ey(i)*alpha
        END DO
!print*,'YS'
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*eb(i)*alpha
        END DO
!print*,'BC'
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*eb(i)*alpha
        END DO
!print*,'BS'
Else
        PD_Param_ID = PD_Param_ID
End If

! ECOM2 model
! **********************************************************************

     ELSE IF (ECOM_param_glb == 2 ) then
      PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ed(i)*alpha
        END DO
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i)*alpha
        END DO
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i)*alpha
        END DO
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(2*del_u)*ed(i)*alpha
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(2*del_u)*ed(i)*alpha
        END DO
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(4*del_u)*ed(i)*alpha
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(4*del_u)*ed(i)*alpha
        END DO
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*eb(i)*alpha
        END DO
! S term
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        DO i=1,3
        fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*eb(i)*alpha
        END DO
Else
        PD_Param_ID = PD_Param_ID
End If
      END IF 


     fx=-fsrp(1)
     fy=-fsrp(2)
     fz=-fsrp(3)



! use the shadow coefficient for scaling the SRP effect
!-------------------------------------------------------
IF (lambda .lt. 1)THEN
!print*,'beta=, lambda=, del_u=', beta*180/Pi, lambda, del_u*180/Pi
   DO i=1,3
   fsrp(i)=0.0d0
   END DO
   srpcoef(1) = lambda*srpcoef(1)
   IF (ECOM_param_glb == 1 ) then
      DO i=1,3
      fsrp(i) = fsrp(i) +   srpcoef(1)*sclfa*ed(i)*alpha              &
                        +   srpcoef(2)*sclfa*ey(i)*alpha              &
                        +   srpcoef(3)*sclfa*eb(i)*alpha              &
                        +   srpcoef(4)*sclfa*DCOS(del_u)*ed(i)*alpha  &
                        +   srpcoef(5)*sclfa*DSIN(del_u)*ed(i)*alpha  &
                        +   srpcoef(6)*sclfa*DCOS(del_u)*ey(i)*alpha  &
                        +   srpcoef(7)*sclfa*DSIN(del_u)*ey(i)*alpha  &
                        +   srpcoef(8)*sclfa*DCOS(del_u)*eb(i)*alpha  &
                        +   srpcoef(9)*sclfa*DSIN(del_u)*eb(i)*alpha  
      END DO
   ELSE 
      DO i=1,3
      fsrp(i) = fsrp(i) +   srpcoef(1)*sclfa*ed(i)*alpha              &
                        +   srpcoef(2)*sclfa*ey(i)*alpha              &
                        +   srpcoef(3)*sclfa*eb(i)*alpha              &
                        +   srpcoef(4)*sclfa*DCOS(2*del_u)*ed(i)*alpha  &
                        +   srpcoef(5)*sclfa*DSIN(2*del_u)*ed(i)*alpha  &
                        +   srpcoef(6)*sclfa*DCOS(4*del_u)*ed(i)*alpha  &
                        +   srpcoef(7)*sclfa*DSIN(4*del_u)*ed(i)*alpha  &
                        +   srpcoef(8)*sclfa*DCOS(del_u)*eb(i)*alpha  &
                        +   srpcoef(9)*sclfa*DSIN(del_u)*eb(i)*alpha
      END DO
   END IF
      fx=-fsrp(1)
      fy=-fsrp(2)
      fz=-fsrp(3)
END IF
!----------------------------------------------------------------------------------------------


     END IF 

! end of ECOM-based model
! ---------------------------------------------------------

END
