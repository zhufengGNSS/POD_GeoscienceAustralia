
SUBROUTINE force_srp (mjd,GM,prnnum,satsvn,eclpf,srpid,r,v,r_sun,fx,fy,fz )


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
! - mjd          : mjd
! - GM           : the earth gravitational constant  
! - prnnum       : satellite PRN number
! - satsvn       : satellite SVN
! - eclpf        : =0: in the sun illumination; =1 in the satellite eclipse
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
! Created:	01-Nov-2017
! 
! Changes:      10-Oct-2018   Tzupang Tseng: modify the Cannonball model to be 
!                                            similar to the box-wing model in
!                                            preparation for BDS with the ON
!                                            attitude mode.
!
! Copyright:  GEOSCIENCE AUSTRALIA, AUSTRALIA
! ----------------------------------------------------------------------

      USE mdl_precision
      USE mdl_num
      USE m_satinfo
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      INTEGER                           :: eclpf,srpid,ECOM
      INTEGER                           :: prnnum,BLKNUM
      REAL (KIND = prec_q), DIMENSION(3) :: r,v,r_sun
      REAL (KIND = prec_q)               :: fx,fy,fz
      REAL (KIND = prec_q)               :: mjd
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
      REAL (KIND = prec_q), DIMENSION(3) :: fsrp
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(9) :: srpcoef 
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q), DIMENSION(4) :: AREA1,REFL1,DIFU1,ABSP1

      INTEGER              :: zta
      INTEGER              :: i,j,k,m
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
!  BLKNUM = 3 ! for GPS Block IIA testing
! ---------------------------------------------------------------------

! SRP model for different GNSS constellations
      if(prnnum.le.100)then  
      ECOM = 1 
      else
      ECOM = 2
      end if

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

! Initialize the coefficients of the ECOM-based SRP model
      do m=1,9
      srpcoef(m)= 0.0d0
      end do

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
! The unit vector ey = ez x ed, parallel to the rotation axis of solar panel
      CALL cross_product (ez,ed,ey)

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


      if (srpid .eq. 1) then
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
     else if (srpid .eq. 2) then
! initialize the SRP force
! -------------------------
    do i=1,3
     fsrp(i)=0.0d0
     end do

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


     else if (srpid .eq. 3) then
! ECOM2 model
! ***************************************
! a=D(u)*ed+Y(u)*ey+B(u)*eb
!
! D(u)=D0+Dc*cos2*(del_u)+Ds*sin2*(del_u)
!        +Dc*cos4*(del_u)+Ds*sin4*(del_u)
!
! Y(u)=Y0
!
! B(u)=B0+Bc*cos1*(del_u)+Bs*sin1*(del_u)
! ***************************************
! srpcoef(1) = D0
! srpcoef(2) = Y0
! srpcoef(3) = B0
! srpcoef(4) = D2_c
! srpcoef(5) = D4_c
! srpcoef(6) = D2_s
! srpcoef(7) = D4_s
! srpcoef(8) = B1_c 
! srpcoef(9) = B1_s
!=======================================


      if (ECOM .eq. 2) then

      srpcoef(1) = -0.91647d-7
      srpcoef(2) =  0.00811d-7
      srpcoef(3) = -0.00652d-7
      srpcoef(4) =  0.00064d-7
      srpcoef(5) =  0.00073d-7
      srpcoef(6) = -0.00065d-7
      srpcoef(7) = -0.00142d-7
      srpcoef(8) =  0.01236d-7
      srpcoef(9) = -0.00968d-7




! ECOM1 model
! ***************************************
! a=D(u)*ed+Y(u)*ey+B(u)*eb
!
! D(u)=D0+Dc*cos1*(del_u)+Ds*sin1*(del_u)
!        
!
! Y(u)=Y0+Yc*cos1*(del_u)+Ys*sin1*(del_u)
!
! B(u)=B0+Bc*cos1*(del_u)+Bs*sin1*(del_u)
! ***************************************
! srpcoef(1) = D0
! srpcoef(2) = Y0
! srpcoef(3) = B0
! srpcoef(4) = Dc
! srpcoef(5) = Ds
! srpcoef(6) = Yc
! srpcoef(7) = Ys
! srpcoef(8) = Bc
! srpcoef(9) = Bs
!=======================================


      else if (ECOM .eq. 1) then
      srpcoef(1) = -0.91664d-7
      srpcoef(2) =  0.00831d-7
      srpcoef(3) = -0.00617d-7
      srpcoef(4) =  0.00067d-7
      srpcoef(5) =  0.00035d-7
      srpcoef(6) = -0.00051d-7
      srpcoef(7) = -0.00100d-7
      srpcoef(8) =  0.01228d-7
      srpcoef(9) = -0.00894d-7

      end if


! initialize the SRP force
! -------------------------
     do i=1,3
     fsrp(i)=0.0d0
     end do


! A scaling factor is applied to ECOM model as suggested by Bernese
!******************************************************************
     sclfa=(AU/Ds)**2

     do i=1,3

     if(ECOM .eq. 2) then
     fsrp(i)=fsrp(i) + srpcoef(1)*sclfa*ed(i) &
                     + srpcoef(2)*sclfa*ey(i) &
                     + srpcoef(3)*sclfa*eb(i) &
                     + srpcoef(4)*sclfa*ed(i)*DCOS(2.0d0*del_u) &
                     + srpcoef(5)*sclfa*ed(i)*DCOS(4.0d0*del_u) &
                     + srpcoef(6)*sclfa*ed(i)*DSIN(2.0d0*del_u) &
                     + srpcoef(7)*sclfa*ed(i)*DSIN(4.0d0*del_u) &
                     + srpcoef(8)*sclfa*eb(i)*DCOS(1.0d0*del_u) &
                     + srpcoef(9)*sclfa*eb(i)*DSIN(1.0d0*del_u)

     else

     fsrp(i)=fsrp(i) + srpcoef(1)*sclfa*ed(i) &
                     + srpcoef(2)*sclfa*ey(i) &
                     + srpcoef(3)*sclfa*eb(i) &
                     + srpcoef(4)*sclfa*ed(i)*DCOS(del_u) &
                     + srpcoef(5)*sclfa*ed(i)*DSIN(del_u) &
                     + srpcoef(6)*sclfa*ey(i)*DCOS(del_u) &
                     + srpcoef(7)*sclfa*ey(i)*DSIN(del_u) &
                     + srpcoef(8)*sclfa*eb(i)*DCOS(del_u) &
                     + srpcoef(9)*sclfa*eb(i)*DSIN(del_u)

     end if
 
     end do

     fx=fsrp(1)
     fy=fsrp(2)
     fz=fsrp(3)


     end if

! end of ECOM model
! ---------------------------------------------------------

END
