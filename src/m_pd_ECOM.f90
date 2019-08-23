MODULE m_pd_ECOM


! ----------------------------------------------------------------------
! MODULE: m_pd_ECOM
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the m_pd_ECOM subroutine
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng, Geosceince Australia, Australia 
!               
!               
! Created:      05-12-2018
! ----------------------------------------------------------------------

      IMPLICIT NONE
      !SAVE


Contains


SUBROUTINE pd_ECOM (lambda, eBX_ecl, GM, prnnum, eclipsf, r, v, r_sun, Asrp)


! ----------------------------------------------------------------------
! SUBROUTINE: pd_ECOM.f90
! ----------------------------------------------------------------------
! Purpose:
! This subrutine is used to compute the partial derivative of force
! w.r.t. ECOM SRP parameters.  
! ----------------------------------------------------------------------
! Input arguments:
! - prnnum       : satellite PRN number 
! - r            : satellite position vector (m)
! - v            : satellite velocity vector
! - r_sun        : Sun position vector
! - lambda       : shadow coefficient
! - eBX_ecl      : dynamic ex of satellite body frame
! 
! Output arguments:
! - Asrp         : Partial derivative output 
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	08-Jan-2018
!
! Changes:    11-12-2018  Tzupang Tseng : make the PD matrix dynamic
!             31-01-2019  Tzupang Tseng : change the definition of ey by
!                                         dividing the length of ey
!             20-02-2019  Tzupang Tseng : create a function for switching on and
!                                         off some particular coefficients in ECOM models
!             22-03-2019  Tzupang Tseng : set a simple condition for the eclipsed satellites
!                                         where only D0 partials are setup to zero 
!             02-05-2019  Tzupang Tseng : use the coefficient from the shadow.f90 for scaling
!                                         the SRP effect
! 
! Copyright:  GEOSCIENCE AUSTRALIA
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int4),INTENT(IN)         :: prnnum
      REAL (KIND = prec_d) , Dimension(3), INTENT(IN) :: eBX_ecl
      REAL (KIND = prec_q), DIMENSION(3),INTENT(IN) :: r,v
      REAL (KIND = prec_q), DIMENSION(3),INTENT(IN) :: r_sun
      INTEGER (KIND = 4), INTENT(IN) :: eclipsf
      REAL (KIND = prec_q),INTENT(IN) :: GM
      REAL (KIND = prec_q),INTENT(IN) :: lambda
! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: AU,Pi
      REAL (KIND = prec_q) :: Cr,ab,Lab
      REAL (KIND = prec_q) :: ANG,E,Edot
      REAL (KIND = prec_q) :: Ds,sclfa

      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
     ! REAL (KIND = prec_q) :: Asrp(3,9)
      REAL (KIND = prec_q), DIMENSION(3) :: er,ed,ey,eb,ex,ez,ev,en
      REAL (KIND = prec_q), DIMENSION(3) :: yy
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      INTEGER              :: i,j,k,ECOM
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: u_sat,i_sat,omega_sat
      INTEGER              :: ex_i
! ----------------------------------------------------------------------
! Sun-related variables
! ----------------------------------------------------------------------
       REAL (KIND = prec_q) :: u_sun,beta,del_u, dang
       REAL (KIND = prec_q), DIMENSION(3) :: r_sun1,r_sun2

      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: N_param, PD_Param_ID
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE,INTENT(OUT) :: Asrp
! ----------------------------------------------------------------------
      REAL (KIND = 8)      :: II, KN, U

! ----------------------------------------------------------------------
! Numerical Constants
      AU = 1.4959787066d11 ! (m)
      Pi = 4*atan(1.0d0)
    ex_i = 0 ! change the definition of the unit vector ex
             ! ex_i = 0 (default)
             !      = 1 (using dynamic ex vector from attitude routine)
! ---------------------------------------------------------------------
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
!      CALL XYZELE(GM, r, v, II, U, KN)
!      u_sat = U
!      i_sat = II
!      omega_sat = KN

! compute the sun position in the satellite orbit plane by rotating big Omega_sat and i_sat,
! allowing us for the computation of u_sun and sun elevation angles (beta)
!   sun                 sun
!  r    = R(big_w,i) x r
!   kep                 ICF
!
      DO i=1,3
         DO j=1,3
          R33(i,j)=1.0d0
         END DO
      END DO 

      R33(1,3)=0.0d0
      R33(2,3)=0.0d0
      R33(3,1)=0.0d0
      R33(3,2)=0.0d0

      DO i=1,3
         DO j=1,3
          R11(i,j)=1.0d0
         END DO 
      END DO 

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
      beta  = atan2(r_sun2(3),sqrt(r_sun2(1)**2+r_sun2(2)**2)) ! in rad

      del_u = u_sat - u_sun ! in rad
      IF (del_u*180/Pi .GT.360.0d0) THEN
      del_u=del_u-2*Pi
      ELSE IF(del_u*180/Pi .LT.0.0d0) THEN
      del_u=del_u+2*Pi
      END IF 

! ----------------------------------------------------------------------
! Partial derivatives w.r.t. unknown parameters

!ALLOCATE (Asrp(3,N_param), STAT = AllocateStatus)
IF (ECOM_param_glb /= 0) THEN
! Bias partial derivatives matrix allocation
PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
ELSE
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
ELSE
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
ELSE
        PD_Param_ID = PD_Param_ID
End IF

! CPR partial derivatives matrix allocation

If (ECOM_CPR_glb(1) == 1) THEN
        PD_Param_ID = PD_Param_ID + 2
ELSE
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(2) == 1) THEN
        PD_Param_ID = PD_Param_ID + 2
ELSE
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(3) == 1) THEN
        PD_Param_ID = PD_Param_ID + 2
ELSE
        PD_Param_ID = PD_Param_ID
End If

IF (NPARAM_glb /= PD_Param_ID) THEN
PRINT*, 'THE NUMBER OF FORCE PARAMETERS IS NOT CONSISTENT'
PRINT*,           'NPARAM_glb  =', NPARAM_glb
PRINT*,           'PD_Param_ID =', PD_Param_ID
END IF

END IF


ALLOCATE (Asrp(3,PD_Param_ID), STAT = AllocateStatus)


! A scaling factor is applied to ECOM model
!******************************************************************
sclfa=(AU/Ds)**2

! ECOM1 model
! ***************************************

     IF (ECOM_param_glb == 1) THEN
!print*,'ECOM1 partials'
PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*1.0d0*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*1.0d0*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*1.0d0*ed(3)
!print*,'D0=',PD_Param_ID
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*1.0d0*ey(1)
        Asrp (2,PD_Param_ID) = -sclfa*1.0d0*ey(2)
        Asrp (3,PD_Param_ID) = -sclfa*1.0d0*ey(3)
!print*,'Y0=',PD_Param_ID
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*1.0d0*eb(1)
        Asrp (2,PD_Param_ID) = -sclfa*1.0d0*eb(2)
        Asrp (3,PD_Param_ID) = -sclfa*1.0d0*eb(3)
!print*,'B0=',PD_Param_ID
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DCOS(del_u)*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*DCOS(del_u)*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*DCOS(del_u)*ed(3)
!print*,'DC=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DSIN(del_u)*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*DSIN(del_u)*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*DSIN(del_u)*ed(3)
!print*,'DS=',PD_Param_ID
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DCOS(del_u)*ey(1)
        Asrp (2,PD_Param_ID) = -sclfa*DCOS(del_u)*ey(2)
        Asrp (3,PD_Param_ID) = -sclfa*DCOS(del_u)*ey(3)
!print*,'YC=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DSIN(del_u)*ey(1)
        Asrp (2,PD_Param_ID) = -sclfa*DSIN(del_u)*ey(2)
        Asrp (3,PD_Param_ID) = -sclfa*DSIN(del_u)*ey(3)
!print*,'YS=',PD_Param_ID
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        Asrp(1,PD_Param_ID) = -sclfa*DCOS(del_u)*eb(1)
        Asrp(2,PD_Param_ID) = -sclfa*DCOS(del_u)*eb(2)
        Asrp(3,PD_Param_ID) = -sclfa*DCOS(del_u)*eb(3)
!print*,'BC=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        Asrp(1,PD_Param_ID) = -sclfa*DSIN(del_u)*eb(1)
        Asrp(2,PD_Param_ID) = -sclfa*DSIN(del_u)*eb(2)
        Asrp(3,PD_Param_ID) = -sclfa*DSIN(del_u)*eb(3)
!print*,'BS=',PD_Param_ID
Else
        PD_Param_ID = PD_Param_ID
End If


! ECOM2 model
! ***************************************
     ELSE IF (ECOM_param_glb == 2) THEN

PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then

        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*1.0d0*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*1.0d0*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*1.0d0*ed(3)
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(2) == 1) Then

        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*1.0d0*ey(1)
        Asrp (2,PD_Param_ID) = -sclfa*1.0d0*ey(2)
        Asrp (3,PD_Param_ID) = -sclfa*1.0d0*ey(3)
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_Bias_glb(3) == 1) Then

        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*1.0d0*eb(1)
        Asrp (2,PD_Param_ID) = -sclfa*1.0d0*eb(2)
        Asrp (3,PD_Param_ID) = -sclfa*1.0d0*eb(3)
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DCOS(2*del_u)*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*DCOS(2*del_u)*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*DCOS(2*del_u)*ed(3)
! S term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DSIN(2*del_u)*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*DSIN(2*del_u)*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*DSIN(2*del_u)*ed(3)
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DCOS(4*del_u)*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*DCOS(4*del_u)*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*DCOS(4*del_u)*ed(3)
! S term
        PD_Param_ID = PD_Param_ID + 1
        Asrp (1,PD_Param_ID) = -sclfa*DSIN(4*del_u)*ed(1)
        Asrp (2,PD_Param_ID) = -sclfa*DSIN(4*del_u)*ed(2)
        Asrp (3,PD_Param_ID) = -sclfa*DSIN(4*del_u)*ed(3)
Else
        PD_Param_ID = PD_Param_ID
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        Asrp(1,PD_Param_ID) = -sclfa*DCOS(del_u)*eb(1)
        Asrp(2,PD_Param_ID) = -sclfa*DCOS(del_u)*eb(2)
        Asrp(3,PD_Param_ID) = -sclfa*DCOS(del_u)*eb(3)
! S term
        PD_Param_ID = PD_Param_ID + 1
        Asrp(1,PD_Param_ID) = -sclfa*DSIN(del_u)*eb(1)
        Asrp(2,PD_Param_ID) = -sclfa*DSIN(del_u)*eb(2)
        Asrp(3,PD_Param_ID) = -sclfa*DSIN(del_u)*eb(3)
Else
        PD_Param_ID = PD_Param_ID
End If

     END IF 

! end of ECOM2 model
! ==================================================================

! use the shadow coefficient for scaling the SRP effect
!-------------------------------------------------------
IF (lambda .lt. 1) THEN
!IF (abs(beta*180.0d0/Pi) .lt. 13.87 .and. del_u*180.0d0/Pi .gt. 167 .and. del_u*180.0d0/Pi .lt. 193) then
!print*,'lambda=, del_u=', lambda, del_u*180/Pi
Asrp(1:3,1) = lambda*Asrp(1:3,1)

END IF


!----------------------------------------------------------------------------------------------

END SUBROUTINE

END MODULE 
