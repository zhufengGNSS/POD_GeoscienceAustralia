
SUBROUTINE force_srp (GM,prnnum,eclpf,srpid,r,v,r_sun,fx,fy,fz )


! ----------------------------------------------------------------------
! SUBROUTINE: force_srp.f90
! ----------------------------------------------------------------------
! Purpose:
! Acceleration due to the solar radiation pressure 
! Computation of SRP acceleration using various SRP models, 
! such as a simply cannonball model (srpid=1),a box-wing model (srpid=2) 
! and a ECOM (D2B1) model (srpid=3) 
! ----------------------------------------------------------------------
! Input arguments:
! - prnnum       : satellite PRN number 
! - eclpf        : =0: in the sun illumination; =1 in the satellite eclipse
! - srpid        : =1: a simply cannonball model; 
!                  =2: box-wing model;
!                  =3: ECOM (D2B1)model
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
! Copyright:  GEOSCIENCE AUSTRALIA
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      INTEGER                           :: eclpf,srpid
      INTEGER                           :: prnnum
!      REAL (KIND = prec_q), INTENT(IN), DIMENSION(3) :: r,v,r_sun
!      REAL (KIND = prec_q), INTENT(OUT) :: fx,fy,fz
      REAL (KIND = prec_q), DIMENSION(3) :: r,v,r_sun
      REAL (KIND = prec_q)               :: fx,fy,fz

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: Ps,AU,Pi
      REAL (KIND = prec_q) :: Cr
      REAL (KIND = prec_q) :: ANG,GM
      REAL (KIND = prec_q) :: Ds,sclfa

      REAL (KIND = prec_q) :: R1(3,3),R3(3,3)
      REAL (KIND = prec_q), DIMENSION(3) :: er,ed,ey,eb
      REAL (KIND = prec_q), DIMENSION(3) :: fsrp
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(9) :: srpcoef 
      INTEGER              :: zta
      INTEGER              :: i,j
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

! computation of the factor zta related to the Sun illumination
      if (eclpf.eq.0) then
         zta=1 ! in the sun light area
      else
         zta=0 ! in the eclipse
         fx=0.0d0
         fy=0.0d0
         fz=0.0d0
!      return
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
          R3(i,j)=1.0d0
         end do
      end do
      R3(1,3)=0.0d0
      R3(2,3)=0.0d0
      R3(3,1)=0.0d0
      R3(3,2)=0.0d0
      R3(3,3)=1.0d0

      do i=1,3
         do j=1,3
          R1(i,j)=1.0d0
         end do
      end do
      R1(1,1)=1.0d0
      R1(1,2)=0.0d0
      R1(1,3)=0.0d0
      R1(2,1)=0.0d0
      R1(2,3)=0.0d0


      CALL iau_RZ(omega_sat,R3)
      CALL iau_RX(i_sat,R1)
      CALL matrix_Rr(R3,r_sun,r_sun1)
      CALL matrix_Rr(R1,r_sun1,r_sun2)

      u_sun = atan2(r_sun2(2),r_sun2(1)) ! in rad

!      if (u_sun*180/Pi .lt. 0.0d0) then
!      u_sun=u_sun+2*Pi
!      end if  

      beta  = atan2(r_sun2(3),sqrt(r_sun2(1)**2+r_sun2(2)**2)) ! in rad
      del_u = u_sat - u_sun


      if (srpid .eq. 1) then
! A simply cannonball model
! *********************************
!  a=-zta*Cr*(A/m)P(1AU/r)^2*Dr(i)
! ********************************

! GPS satellites 
      if (prnnum.ge.1 .and. prnnum.le.32) then

! GPS Block IIR types
      if (prnnum.ge.11 .and. prnnum.le.23 .or. prnnum.eq.2 .or. &
       prnnum.eq.5 .or. prnnum.eq.7 .or. prnnum.eq.28 .or. &
       prnnum.eq.29 .or. prnnum.eq.31) then 
      Z_SIDE = 3.75D0 ! surface-to-mass ratio 
      X_SIDE = 3.05D0
      A_SOLAR= 13.60D0
      MASS   = 1100.0D0
! Block IIF types
      else 
      Z_SIDE = 5.40D0 ! surface-to-mass ratio
      X_SIDE = 5.72D0
      A_SOLAR= 22.25D0
      MASS   = 1555.3D0
      end if
      end if

! GLONASS satellites
      if (prnnum.gt.101 .and. prnnum.le.126)then
      Z_SIDE = 0.877D0 ! surface-to-mass ratio
      X_SIDE = 1.258D0
      A_SOLAR= 23.616D0
      MASS   = 1415.0D0
      end if

! BDS satellites
      if (prnnum.gt.401 .and. prnnum.le.426)then
      Z_SIDE = 3.69D0 ! surface-to-mass ratio
      X_SIDE = 4.5D0
      A_SOLAR= 22.44D0
      MASS   = 1900.0D0
      end if

! ----------------------------------------------------------------------

! The main surface area face toward to the Sun using the SAT->SUN and SAT->EARTH
! vectors
     ANG=acos(ed(1)*er(1)+ed(2)*er(2)+ed(3)*er(3))*180.0d0/Pi
     
     if (abs(ANG) .le. 30.0d0) then
     AREA=Z_SIDE+A_SOLAR
     else
     AREA=X_SIDE+A_SOLAR
     end if 

! Cartesian counterparts (fx,fy,fz) of acceleration fr
      fx = -zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(1)
      fy = -zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(2)
      fz = -zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(3) 

!write(*,*) fx,fy,fz
! end of the cannonball model
!--------------------------------------------------------------------------


! Box-wing model
! **************
     else if (srpid .eq. 2) then
! initialize the SRP force
! -------------------------

     else if (srpid .eq. 3) then
! ECOM (D2B1) model
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
! srpcoef(8) = D1_c 
! srpcoef(9) = D1_c


! GPS satellites
      if (prnnum.ge.1 .and. prnnum.le.32) then

! GPS Block IIR types
      if (prnnum.ge.11 .and. prnnum.le.23 .or. prnnum.eq.2 .or. &
       prnnum.eq.5 .or. prnnum.eq.7 .or. prnnum.eq.28 .or. &
       prnnum.eq.29 .or. prnnum.eq.31) then

       srpcoef(1) = -0.0070d-7
       srpcoef(2) = -0.0020d-7
       srpcoef(3) =  0.0050d-7
       srpcoef(4) = -0.0050d-7
       srpcoef(5) = -0.0084d-7
       srpcoef(6) =  0.0050d-7
       srpcoef(7) =  0.0050d-7
       srpcoef(8) = -0.0008d-7 
       srpcoef(9) =  0.0050d-7

! Block IIF types
      else
      srpcoef(1) = -1.0070d-7
      srpcoef(2) =  0.0015d-7
      srpcoef(3) = -0.0050d-7
      srpcoef(4) = -0.0050d-7
      srpcoef(5) =  0.0030d-7
      srpcoef(6) = -0.0030d-7
      srpcoef(7) = -0.0050d-7
      srpcoef(8) = -0.0090d-7 
      srpcoef(9) = -0.0050d-7

      end if
      end if

! BDS satellites
      if (prnnum.ge.401 .and. prnnum.le.420) then

! BDS IGSO types
      if (prnnum.ge.406 .and. prnnum.le.410) then

      srpcoef(1) =  1.2330d-7
      srpcoef(2) = -0.0002d-7
      srpcoef(3) = -0.0050d-7
      srpcoef(4) = -0.0050d-7
      srpcoef(5) =  0.0005d-7 
      srpcoef(6) = -0.0050d-7
      srpcoef(7) = -0.0050d-7
      srpcoef(8) = -0.0150d-7 
      srpcoef(9) = -0.0050d-7

! BDS MEO types 
      else if (prnnum.ge.411 .and. prnnum.le.415) then
      srpcoef(1) = -1.3000d-7
      srpcoef(2) =  0.0010d-7
      srpcoef(3) = -0.0020d-7
      srpcoef(4) = -0.0250d-7
      srpcoef(5) =  0.0051d-7 
      srpcoef(6) = -0.0100d-7
      srpcoef(7) = -0.0050d-7
      srpcoef(8) =  0.0080d-7 
      srpcoef(9) = -0.0050d-7

      end if
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
     fsrp(i)=fsrp(i) + srpcoef(1)*sclfa*ed(i) &
                     + srpcoef(2)*sclfa*ey(i) &
                     + srpcoef(3)*sclfa*eb(i) &
                     + srpcoef(4)*sclfa*ed(i)*DCOS(2.0d0*del_u) &
                     + srpcoef(5)*sclfa*ed(i)*DCOS(4.0d0*del_u) &
                     + srpcoef(6)*sclfa*ed(i)*DSIN(2.0d0*del_u) &
                     + srpcoef(7)*sclfa*ed(i)*DSIN(4.0d0*del_u) &
                     + srpcoef(8)*sclfa*eb(i)*DCOS(1.0d0*del_u) &
                     + srpcoef(9)*sclfa*eb(i)*DSIN(1.0d0*del_u)
 
     end do

     fx=fsrp(1)
     fy=fsrp(2)
     fz=fsrp(3)


!write(*,*) fx,fy,fz

     end if

! end of ECOM (D2B1) model
! ---------------------------------------------------------

END
