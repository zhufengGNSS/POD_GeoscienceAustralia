MODULE m_orbinfo


! ----------------------------------------------------------------------
! MODULE: m_orbinfo.f90
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling statorbit subroutine
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng
!
! Copyright: GEOSCIENCE AUSTRALIA, AUSTRALIA
!
! Created:      07-05-2019
! ----------------------------------------------------------------------
      IMPLICIT NONE


      Contains

SUBROUTINE orbinfo (mjd, prnnum, rsat, vsat, beta, del_u, yaw, lambda, angX, angY, angZ, &
                    fr, ft, fn )
! ----------------------------------------------------------------------
! SUBROUTINE: orbinfo.f90
! ----------------------------------------------------------------------
! Purpose:
! Compute the required information for orbital analysis in terms of
! satellite argument of latitude, beta angle, attitude yaw angle, angles
! between the solar radiation direction and the surface normal vectors
! ----------------------------------------------------------------------
! Input parameters:
! - mjd:               Time in MJD
! - prnnum:            Satellite PRN number
! - rsat:              Satellite position vector in ICRF
! - vsat:              Satellite velocity vector in ICRF
!
! Output parameters:
! - beta (rad) :       BETA angle 
! - del_u (rad):       Satellite argument of latitude wrt the SUN
! - yaw (rad)  :       Nominal yaw angle  
! - lambda:            Shadow coefficient
! - angX (rad) :       Incident angle between the solar radiation and the X surface normal vector
! - angY (rad) :       Incident angle between the solar radiation and the Y surface normal vector
! - angZ (rad) :       Incident angle between the solar radiation and the Z surface normal vector
! - fr:                SRP-induced acceleration in radial
! - ft:                SRP-induced acceleration in along-track
! - fn:                SRP-induced acceleration in cross-track
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng
!
! Copyright: GEOSCIENCE AUSTRALIA,AUSTRALIA
!
! Created:      07-05-2019
! ----------------------------------------------------------------------
      USE mdl_precision
      USE mdl_num
      USE mdl_planets
      USE m_shadow
      USE mdl_param
      USE mdl_config
      IMPLICIT NONE


!-------------------------------------------------------------------
      INTEGER (KIND = prec_int4) :: prnnum
      REAL (KIND = prec_d), INTENT(IN) :: mjd
      REAL (KIND = prec_d), DIMENSION(3), INTENT(IN) :: rsat, vsat
      REAL (KIND = prec_q), INTENT(OUT) :: angX, angY, angZ
      REAL (KIND = prec_q), INTENT(OUT) :: beta, del_u, yaw
      REAL (KIND = prec_q), INTENT(OUT) :: lambda
      CHARACTER (KIND = 1) :: ECLTYP

!--------------------------------------------------------------------
      DOUBLE PRECISION  JD, Zbody(6)
      INTEGER (KIND = 4) :: i, j
      REAL (KIND = prec_q), DIMENSION(3) :: rbody
      REAL (KIND = prec_q), DIMENSION(3) :: rSun, rMoon
      REAL (KIND = prec_q), DIMENSION(3) :: r_sun1, r_sun2
      REAL (KIND = prec_q) :: u_sun, yaw2
      INTEGER  NTARG_SUN, NTARG_MOON, NCTR
      REAL (KIND = prec_q), DIMENSION(3) :: ed, ey, eb, ex, en, ev, ez, et
      REAL (KIND = prec_q), DIMENSION(3) :: nd, nz, nx
      REAL (KIND = prec_q), DIMENSION(3) :: yy, ed_on
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: Pi, Ps, AU, sclfa
      REAL (KIND = prec_d) :: GM, Ds
      REAL (KIND = prec_d) :: u_sat, i_sat, omega_sat
      REAL (KIND = prec_q) :: fxo,fyo,fzo
      REAL (KIND = prec_q) :: fx,fy,fz
      REAL (KIND = prec_q) :: fr,ft,fn
      REAL (KIND = prec_q), DIMENSION(3) :: fsrp
      REAL (KIND = prec_q), DIMENSION(3) :: FXX(3),FZZ(3),FSP(3),FXD(3),FZD(3)
      REAL (KIND = prec_q), DIMENSION(:),ALLOCATABLE :: srpcoef
      INTEGER (KIND = prec_int2) :: AllocateStatus
      INTEGER              :: PD_Param_ID
      REAL (KIND = 8)      :: II, KN, U
      REAL (KIND = prec_q) :: F0,alpha
! ----------------------------------------------------------------------
! Satellite physical informaiton
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: X_SIDE,Z_SIDE
      REAL (KIND = prec_q) :: AREA
      REAL (KIND = prec_q) :: A_SOLAR
! ----------------------------------------------------------------------
      INTEGER*4 BLKNUM,SVN,REFF,ERM,ANT,GRD,MONTH
      REAL*8  ACCEL(3),SUN(3)
      REAL*8  YSAT(6)
      CHARACTER (LEN=1) :: GNSSid
      CHARACTER (LEN=100) :: fmt_line
! ----------------------------------------------------------------------


! Numerical Constants

      GM = 0.39860044150D+15
      Ps = 4.5567D-6 ! (Nm^-2)
      AU = 1.496d11 ! (m)
      Pi = 4*atan(1.0d0)
! initialize the SRP force
! -------------------------
     DO i=1,3
     fsrp(i)=0.0d0
     END DO
! Identify the GNSS constellation
! -------------------------------
fmt_line = '(A1)'
READ (PRN, fmt_line ) GNSSid


IF (GNSSid == 'G') THEN
! GPS constellation
! -----------------
! I
         if (BLKTYP=='GPS-I')then
         Z_SIDE = 3.020D0
         X_SIDE = 1.728D0
         A_SOLAR= 6.053D0
         F0 = 4.54d-5

! II and IIA
         else if (BLKTYP=='GPS-II' .or. BLKTYP=='GPS-IIA') then
         Z_SIDE = 2.881D0
         X_SIDE = 2.893D0
         A_SOLAR= 11.871D0
         F0 = 8.695d-5

! IIF
         else if(BLKTYP=='GPS-IIF') then
         Z_SIDE = 5.05D0
         X_SIDE = 4.55D0
         A_SOLAR= 22.25D0
         F0 = 16.7d-5

! IIR
         else if (BLKTYP=='GPS-IIR' .or. BLKTYP=='GPS-IIR-A' .or. &
                  BLKTYP=='GPS-IIR-B'.or.BLKTYP=='GPS-IIR-M') then
         Z_SIDE = 4.25D0
         X_SIDE = 4.11D0
         A_SOLAR= 13.92D0
         F0 = 11.15d-5

! III
         else if (BLKTYP=='GPS-IIIA') then
         Z_SIDE = 4.38D0
         X_SIDE = 6.05D0
         A_SOLAR= 22.25D0
         F0 = 11.0d-5

         end if
ELSE IF (GNSSid == 'R') THEN
! GLONASS constellation
! ---------------------
         if(BLKTYP=='GLO' .or. BLKTYP=='GLO-M' .or. BLKTYP == 'GLO-M+' .or.&
            BLKTYP=='GLO-K1A'.or.BLKTYP == 'GLO-K1B')then
         Z_SIDE = 1.6620D0
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0
! GLONASS-K
         if(BLKTYP=='GLO-K1A'.or.BLKTYP == 'GLO-K1B') F0 = 10.0d-5
! GLONASS-M
         if (BLKTYP=='GLO' .or. BLKTYP=='GLO-M' .or. BLKTYP == 'GLO-M+') F0 = 20.9d-5

         end if
ELSE IF (GNSSid == 'E') THEN
! GALILEO constellation
! ---------------------
         if (BLKTYP=='GAL-1' .or. BLKTYP=='GAL-2') then
         Z_SIDE = 3.002D0
         X_SIDE = 1.323D0
         A_SOLAR= 11.0D0
         F0 = 8.35d-5
         end if
ELSE IF (GNSSid == 'C') THEN
! BDS constellation
! -----------------
         if (BLKTYP=='BDS-2M'.or. BLKTYP=='BDS-2G'.or.BLKTYP=='BDS-2I') then
         Z_SIDE = 3.96D0
         X_SIDE = 4.5D0
         A_SOLAR= 22.44D0
! BDS MEO
         if(BLKTYP=='BDS-2M') F0 = 8.35d-5
! BDS IGSO
         if(BLKTYP=='BDS-2I'.or. BLKTYP=='BDS-2G') F0 = 50.1d-5
         end if

ELSE IF (GNSSid == 'J') THEN
! QZSS constellation
! ------------------
         if (BLKTYP=='QZS-1'.or.BLKTYP=='QZS-2I'.or.BLKTYP=='QZS-2G') then
         Z_SIDE = 6.00D0
         X_SIDE = 12.2D0
         A_SOLAR= 40.0D0
         F0 = 50.1d-5 ! Assumed to be the same with BDS/IGSO
         end if
END IF

! ---------------------------------------------------------------------

! Julian Day Number of the input epoch
      JD = mjd + 2400000.5D0

! Center celestial body: Earth (NCTR = 3)
      NCTR = 3
      NTARG_MOON= 10
      NTARG_SUN = 11
! MOON
! Celestial body's (NTARG) Cartesian coordinates w.r.t. Center body
! (NCTR)
      CALL  PLEPH ( JD, NTARG_MOON, NCTR, Zbody )
! Cartesian coordinates of the celestial body in meters: KM to M
      rbody(1) = Zbody(1) * 1000.D0
      rbody(2) = Zbody(2) * 1000.D0
      rbody(3) = Zbody(3) * 1000.D0
      rMoon = rbody
! SUN
      CALL  PLEPH ( JD, NTARG_SUN, NCTR, Zbody )
      rbody(1) = Zbody(1) * 1000.D0
      rbody(2) = Zbody(2) * 1000.D0
      rbody(3) = Zbody(3) * 1000.D0
      rSun = rbody

      CALL shadow (rsat, rSun, rMoon, lambda, ECLTYP )

!! The unit vector ez SAT->EARTH
      ez(1)=-rsat(1)/sqrt(rsat(1)**2+rsat(2)**2+rsat(3)**2)
      ez(2)=-rsat(2)/sqrt(rsat(1)**2+rsat(2)**2+rsat(3)**2)
      ez(3)=-rsat(3)/sqrt(rsat(1)**2+rsat(2)**2+rsat(3)**2)

! The unit vector ed SAT->SUN
! vector)
      Ds=sqrt((rSun(1)-rsat(1))**2+(rSun(2)-rsat(2))**2+(rSun(3)-rsat(3))**2)
      ed(1)=((rSun(1)-rsat(1))/Ds)
      ed(2)=((rSun(2)-rsat(2))/Ds)
      ed(3)=((rSun(3)-rsat(3))/Ds)

! The unit vector ey = ez x ed/|ez x ed|
      CALL productcross (ez,ed,yy)
      ey(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

! The unit vector of x surface
      CALL productcross (ey,ez,ex)

! The unit vector eb = ed x ey
      CALL productcross (ed,ey,eb)

! the orbit normal vector
!------------------------
      ev(1)=vsat(1)/sqrt(vsat(1)**2+vsat(2)**2+vsat(3)**2)
      ev(2)=vsat(2)/sqrt(vsat(1)**2+vsat(2)**2+vsat(3)**2)
      ev(3)=vsat(3)/sqrt(vsat(1)**2+vsat(2)**2+vsat(3)**2)

      CALL productcross (-ez,ev,en)

! the along track direction
  CALL productcross (en,-ez,et)

! computation of the satellite argument of latitude and orbit
! inclination
      CALL kepler_z2k (rsat, vsat, GM, kepler)
      u_sat = kepler(9)*Pi/180.d0
      i_sat = kepler(3)*Pi/180.d0
      omega_sat = kepler(4)*Pi/180.d0
!print*,'ACS POD=', kepler(3), kepler(9), kepler(4)
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

! rotate big omega to make the X-axis of the celestial frame consistent
! with the
! direction of right ascension of ascending node
      CALL R3(omega_sat,R33)
! rotate inclination to make the XY plane of the celestial frame
! consistent with
! the orbit plane
      CALL R1(i_sat,R11)
! convert the sun position in the celestial frame to the orbit plane
      CALL matrix_Rr(R33,rSun,r_sun1)
      CALL matrix_Rr(R11,r_sun1,r_sun2)

      beta  = atan2(r_sun2(3),sqrt(r_sun2(1)**2+r_sun2(2)**2)) ! in rad
!write (*,*) beta*180.0d0/Pi

! yaw angle
      yaw= acos(ex(1)*et(1)+ex(2)*et(2)+ex(3)*et(3))

      IF (beta*180/Pi .gt. 0.d0) yaw= -yaw ! In accordance with the IGS convention (yaw-steering only)


      u_sun = atan2(r_sun2(2),r_sun2(1)) ! in rad
      del_u = u_sat - u_sun
      if (del_u*180/Pi .gt.360.0d0) then
      del_u=del_u-2*Pi
      else if(del_u*180/Pi .lt.0.0d0) then
      del_u=del_u+2*Pi
      end if

! BDS yaw angles for the orbit normal attitude
!     if(prnnum > 300) then
!        if (abs(beta*180.0d0/Pi) < 4.d0) ex = et

!     end if
      if(BLKID == 301)then
      CALL productcross (ez,ev,yy)
      ey(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

      CALL productcross (ed,ey,yy)
      eb(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      eb(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      eb(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

      CALL productcross (ey,eb,yy)
      ed_on(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ed_on(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ed_on(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

      yaw= acos(ed(1)*ed_on(1)+ed(2)*ed_on(2)+ed(3)*ed_on(3))
      IF (beta*180/Pi .gt. 0.d0) yaw= -yaw

      elseif(BLKID == 302 .or. BLKID == 303)then
      if (abs(beta*180.0d0/Pi) < 4.5d0) then
      CALL productcross (ez,ev,yy)
      ey(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ey(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

      CALL productcross (ed,ey,yy)
      eb(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      eb(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      eb(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

      CALL productcross (ey,eb,yy)
      ed_on(1)=yy(1)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ed_on(2)=yy(2)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)
      ed_on(3)=yy(3)/sqrt(yy(1)**2+yy(2)**2+yy(3)**2)

      yaw= acos(ed(1)*ed_on(1)+ed(2)*ed_on(2)+ed(3)*ed_on(3))
      IF (beta*180/Pi .gt. 0.d0) yaw= -yaw

      end if
      end if


!========================================
! angles between ed and each surface(k)
! k=1: +X
!   2: +Y
!   3: +Z
!   4: Solar panel
!=======================================
! SRP force direction

      cosang(1)=ed(1)*ex(1)+ed(2)*ex(2)+ed(3)*ex(3)
      cosang(2)=ed(1)*ey(1)+ed(2)*ey(2)+ed(3)*ey(3)
      cosang(3)=ed(1)*ez(1)+ed(2)*ez(2)+ed(3)*ez(3)
      cosang(4)=ed(1)*ed(1)+ed(2)*ed(2)+ed(3)*ed(3)

!do j =1,4
!if(cosang(j)<0.d0) print*,'cosang < 0.0, j =', j
!end do
nx = ex
nz = ez
nd = ed

if(cosang(3)<0.d0) nz = ez*(-1.d0)

do j = 1,3
FXX(j) = Ps/MASS * X_SIDE * cosang(1) * nx(j)
FZZ(j) = Ps/MASS * Z_SIDE * cosang(3) * nz(j)
FSP(j) = Ps/MASS * A_SOLAR* nd(j)
end do



      angX = acos(cosang(1))
      angY = acos(cosang(2))
      angZ = acos(cosang(3))

sclfa=(AU/Ds)**2

! A PRIORI SRP ACCELERQTION
!---------------------------
! SIMPLE BOX-WING
      if (Flag_BW_cfg == 1 ) then
         fxo=Ps/MASS*(0.3*X_SIDE*cosang(1)*nx(1)+Z_SIDE*cosang(3)*nz(1)+0.7*A_SOLAR*cosang(4)*nd(1))
         fyo=Ps/MASS*(0.3*X_SIDE*cosang(1)*nx(2)+Z_SIDE*cosang(3)*nz(2)+0.7*A_SOLAR*cosang(4)*nd(2))
         fzo=Ps/MASS*(0.3*X_SIDE*cosang(1)*nx(3)+Z_SIDE*cosang(3)*nz(3)+0.7*A_SOLAR*cosang(4)*nd(3))
               
         alpha = sqrt(fxo**2+fyo**2+fzo**2)

! BOX-WING model from the repro3 routine
! --------------------------------------
      else if (Flag_BW_cfg == 2 ) then
         REFF = 0    ! 0: inertial frame,  1: satellite body-fixed frame,
                     ! 2: sun-fixed frame, 3: orbital frame

         YSAT(1:3) = rsat
         YSAT(4:6) = vsat
         CALL SRPFBOXW(REFF,YSAT,rSun,BLKID,SVNID,ACCEL)
         alpha = sqrt(ACCEL(1)**2+ACCEL(2)**2+ACCEL(3)**2)
        
       fxo=ACCEL(1)
       fyo=ACCEL(2)
       fzo=ACCEL(3)
 
      else if (Flag_BW_cfg == 0) then
         alpha = F0/MASS

      else
         alpha = 1.d0
      end if

IF (ECOM_param_glb <= 2) THEN
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


ELSE IF(ECOM_param_glb == 3) THEN
PD_Param_ID = NPARAM_glb
END IF


ALLOCATE (srpcoef(PD_Param_ID), STAT = AllocateStatus)
if (AllocateStatus .ne. 0) then
        print *, "failed to allocate srpcoef"
        goto 100
end if

srpcoef = 0.d0

! ECOM1 model
! ***********************************************************************
      IF (ECOM_param_glb == 1 ) then
      PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
        IF (lambda .lt. 1) srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID)
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
        IF (lambda .lt. 1) srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID)
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
! SIMPLE BOX WING
! *******************************************************************
      ELSE IF (ECOM_param_glb == 3 ) THEN
      PD_Param_ID = 9
      DO PD_Param_ID = 1, 9
         IF (PD_Param_ID == 1) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         IF (lambda .lt. 1) srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*FXX(i)*ed(i)
         END DO
         ELSE IF (PD_Param_ID == 2) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         IF (lambda .lt. 1) srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*FZZ(i)*ed(i)
         END DO
         ELSE IF (PD_Param_ID == 3) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         IF (lambda .lt. 1) srpcoef(PD_Param_ID) = lambda*srpcoef(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*FSP(i)
         END DO
         ELSE IF (PD_Param_ID == 4) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*ey(i)
         END DO
         ELSE IF (PD_Param_ID == 5) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*FXX(i)*eb(i)
         END DO
         ELSE IF (PD_Param_ID == 6) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*FZZ(i)*eb(i)
         END DO
         ELSE IF (PD_Param_ID == 7) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*eb(i)
         END DO
         ELSE IF (PD_Param_ID == 8) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DCOS(del_u)*eb(i)
         END DO
         ELSE IF (PD_Param_ID == 9) THEN
         srpcoef (PD_Param_ID) = ECOM_accel_glb(PD_Param_ID)
         DO i=1,3
         fsrp(i) = fsrp(i) + srpcoef(PD_Param_ID)*sclfa*DSIN(del_u)*eb(i)
         END DO
         END IF
      END DO
      END IF

     fx=-fsrp(1)
     fy=-fsrp(2)
     fz=-fsrp(3)

     fr = fx*(-ez(1))+fy*(-ez(2))+fz*(-ez(3))
     ft = fx*et(1)+fy*et(2)+fz*et(3)
     fn = fx*en(1)+fy*en(2)+fz*en(3)

 100 END SUBROUTINE
END MODULE



