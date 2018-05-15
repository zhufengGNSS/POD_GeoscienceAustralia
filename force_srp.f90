
SUBROUTINE force_srp (GM,prnnum,BLKNUM,eclpf,srpid,r,v,r_sun,beta,del_u,fx,fy,fz )


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
! - beta         : Sun elevation wrt satellite orbit plane
! 
! Output arguments:
! - del_u        : satellite-Sun argument of lattitude       
! - fx,fy,fz:	 : Acceleration's cartesian components (m)
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
      INTEGER                           :: eclpf,srpid,X,ECOM
      INTEGER                           :: prnnum,BLKNUM
      REAL (KIND = prec_q), DIMENSION(3) :: r,v,r_sun
      REAL (KIND = prec_q)               :: fx,fy,fz

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: Ps,AU,Pi
      REAL (KIND = prec_q) :: Cr
      REAL (KIND = prec_q) :: ANG,GM
      REAL (KIND = prec_q) :: Ds,sclfa

      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: surforce(4,3)
      REAL (KIND = prec_q), DIMENSION(3) :: er,ed,ey,eb,ex,en
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
      Cr = 1.4 ! SRP coefficient ranges from 1.3 to 1.5
      AU = 1.4959787066d11 ! (m)
      Pi = 4*atan(1.0d0)
! ---------------------------------------------------------------------

         ! ECOM =1 :ECOM-1 model is applied.
ECOM = 2 ! ECOM =2 :D2B1   model is applied. 

   X = 2 ! X = 1 : forces in the inertial frame.
         ! X = 2 : forces in the orbital frame.
eclpf=0


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


! The unit vector eb = ed x ey
      CALL cross_product (ed,ey,eb)

! The unit vector of x side, which is always illuminated by the
! sun, parallel to the satellite flight/velocity direction
! The unit vector ex
      ex(1)=v(1)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ex(2)=v(2)/sqrt(v(1)**2+v(2)**2+v(3)**2)
      ex(3)=v(3)/sqrt(v(1)**2+v(2)**2+v(3)**2)

! the orbit normal vector
!------------------------

     Call productcross(-er,ex,en)

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
!   sun                 sun
!  r    = R(big_w,i) x r
!   kep                 ICF
!

!      do i=1,3
!         do j=1,3
!          R3(i,j)=1.0d0
!         end do
!      end do
      R33(1,1)=1.0d0
      R33(1,2)=1.0d0
      R33(1,3)=0.0d0
      R33(2,1)=1.0d0
      R33(2,2)=1.0d0
      R33(2,3)=0.0d0
      R33(3,1)=0.0d0
      R33(3,2)=0.0d0
      R33(3,3)=1.0d0


!      do i=1,3
!         do j=1,3
!          R1(i,j)=1.0d0
!         end do
!      end do
      R11(1,1)=1.0d0
      R11(1,2)=0.0d0
      R11(1,3)=0.0d0
      R11(2,1)=0.0d0
      R11(2,2)=1.0d0
      R11(2,3)=1.0d0
      R11(3,1)=0.0d0
      R11(3,2)=1.0d0
      R11(3,3)=1.0d0

! rotate big omega to make the X-axis of the celestial frame consistent with the
! direction of right ascension of ascending node
      CALL R3(omega_sat,R33)

! rotate inclination to make the XY plane of the celestial frame consistent with
! the orbit plane
      CALL R1(i_sat,R11)

! convert the sun position in the celestial frame to the orbit plane
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



      if (srpid .eq. 1) then
! A simply cannonball model
! =================================================================
!  a=-zta*Cr*(A/m)P(1AU/r)^2*Dr(i)
! ==================================================================

! GPS satellites 
      if (prnnum.ge.1 .and. prnnum.le.32) then

! GPS Block IIR types
      if (prnnum.ge.11 .and. prnnum.le.23 .or. prnnum.eq.2 .or. &
       prnnum.eq.5 .or. prnnum.eq.7 .or. prnnum.eq.28 .or. &
       prnnum.eq.29 .or. prnnum.eq.31) then 
      Z_SIDE = 3.75D0 
      X_SIDE = 3.05D0
      A_SOLAR= 13.60D0
      MASS   = 1100.0D0
! Block IIF types
      else 
      Z_SIDE = 5.40D0 
      X_SIDE = 5.72D0
      A_SOLAR= 22.25D0
      MASS   = 1555.3D0
      end if
      end if

! GLONASS satellites
      if (prnnum.gt.101 .and. prnnum.le.126)then
      Z_SIDE = 0.877D0
      X_SIDE = 1.258D0
      A_SOLAR= 23.616D0
      MASS   = 1415.0D0
      end if

! BDS satellites
      if (prnnum.gt.401 .and. prnnum.le.426)then
      Z_SIDE = 3.69D0
      X_SIDE = 4.5D0
      A_SOLAR= 22.44D0
      MASS   = 1900.0D0
      end if


! The main surface area face toward to the Sun using the SAT->SUN and SAT->EARTH
! vectors
     ANG=acos(ed(1)*er(1)+ed(2)*er(2)+ed(3)*er(3))*180.0d0/Pi
     
     if (abs(ANG) .le. 30.0d0) then
     AREA=Z_SIDE+A_SOLAR
     else
     AREA=X_SIDE+A_SOLAR
     end if 

! Cartesian counterparts (fx,fy,fz) of acceleration fsrP
      fsrp(1) = zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(1)
      fsrp(2) = zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(2)
      fsrp(3) = zta*Cr*AREA/MASS*Ps*(AU/Ds)**2*ed(3) 



     if (X .eq. 1)then
! In the inertial frame
!----------------------
     fx=-fsrp(1)
     fy=-fsrp(2)
     fz=-fsrp(3)

     else
! In the orbit frame
!--------------------

     fx=-(fsrp(1)*(-er(1))+fsrp(2)*(-er(2))+fsrp(3)*(-er(3)))
     fy=-(fsrp(1)*ex(1)   +fsrp(2)*ex(2)   +fsrp(3)*ex(3))
     fz=-(fsrp(1)*en(1)   +fsrp(2)*en(2)   +fsrp(3)*en(3))


!     fx=-(fsrp(1)*ed(1)   +fsrp(2)*ed(2)   +fsrp(3)*ed(3))
!     fy=-(fsrp(1)*ey(1)   +fsrp(2)*ey(2)   +fsrp(3)*ey(3))
!     fz=-(fsrp(1)*eb(1)   +fsrp(2)*eb(2)   +fsrp(3)*eb(3))


     end if


!write(*,*) fx,fy,fz
! end of the cannonball model
! =================================================================================



else if (srpid .eq. 2) then

! Box-wing model
! ================================
! initialize the SRP force
! -------------------------
     do i=1,3
     fsrp(i)=0.0d0
     end do


! GPS Block IIR types
      if (prnnum.ge.11 .and. prnnum.le.23 .or. prnnum.eq.2 .or. &
       prnnum.eq.5 .or. prnnum.eq.7 .or. prnnum.eq.28 .or. &
       prnnum.eq.29 .or. prnnum.eq.31) then
      MASS   = 1100.0D0
      
      elseif (prnnum .eq. 4) then
      MASS   = 975.0d0
! Block IIF types
      else
      MASS   = 1555.3D0
      end if


! angles between ed and each surface(k)
! k=1: +X
!   2: +Y
!   3: +Z
!   4: solar panels

     cosang(1)=ed(1)*ex(1)+ed(2)*ex(2)+ed(3)*ex(3)
     cosang(2)=ed(1)*ey(1)+ed(2)*ey(2)+ed(3)*ey(3)
     cosang(3)=ed(1)*er(1)+ed(2)*er(2)+ed(3)*er(3)
     cosang(4)=ed(1)*ed(1)+ed(2)*ed(2)+ed(3)*ed(3)

!write(*,*) cosang(1),cosang(2),cosang(3),cosang(4)

      CALL surfprop(BLKNUM,AREA1,REFL1,DIFU1,ABSP1)
!write(*,*) BLKNUM,AREA1,REFL1,DIFU1,ABSP1 

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
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*er(i))+2*cosang(k)**2*REFL1(k)*er(i)
      end do
        ELSEIF(cosang(k).LT.0D0)THEN
      do i=1,3
     surforce(k,i)=abs(cosang(k))*(ABSP1(k)+DIFU1(k))*(ed(i)+2/3*(-er(i)))+2*cosang(k)**2*REFL1(k)*(-er(i))
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


     if (X .eq. 1)then
! In the inertial frame
!----------------------
     fx=-fsrp(1)
     fy=-fsrp(2)
     fz=-fsrp(3)

     else
! In the orbit frame
!--------------------

     fx=-(fsrp(1)*(-er(1))+fsrp(2)*(-er(2))+fsrp(3)*(-er(3)))
     fy=-(fsrp(1)*ex(1)   +fsrp(2)*ex(2)   +fsrp(3)*ex(3))
     fz=-(fsrp(1)*en(1)   +fsrp(2)*en(2)   +fsrp(3)*en(3))

!     fx=-(fsrp(1)*ed(1)   +fsrp(2)*ed(2)   +fsrp(3)*ed(3))
!     fy=-(fsrp(1)*ey(1)   +fsrp(2)*ey(2)   +fsrp(3)*ey(3))
!     fz=-(fsrp(1)*eb(1)   +fsrp(2)*eb(2)   +fsrp(3)*eb(3))


     end if



!write(*,*) fx,fy,fz
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
! srpcoef(1) = D0
! srpcoef(2) = Y0
! srpcoef(3) = B0
! srpcoef(4) = D2_c
! srpcoef(5) = D4_c
! srpcoef(6) = D2_s
! srpcoef(7) = D4_s
! srpcoef(8) = B1_c 
! srpcoef(9) = B1_s
!=========================================

! GPS satellites
      if (prnnum.ge.1 .and. prnnum.le.32) then
      srpcoef(1) = -0.91698d-7
      srpcoef(2) =  0.00836d-7
      srpcoef(3) = -0.00730d-7
      srpcoef(4) =  0.00053d-7
      srpcoef(5) = -0.00029d-7
      srpcoef(6) = -0.00040d-7
      srpcoef(7) =  0.00010d-7
      srpcoef(8) =  0.01062d-7
      srpcoef(9) = -0.00801d-7

     else

      srpcoef(1) = -1.13184d-7
      srpcoef(2) =  0.00091d-7
      srpcoef(3) =  0.00145d-7
      srpcoef(4) = -0.05658d-7
      srpcoef(5) = -0.00115d-7
      srpcoef(6) = -0.00052d-7
      srpcoef(7) =  0.00805d-7
      srpcoef(8) = -0.02095d-7
      srpcoef(9) = -0.00435d-7



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

     if (X .eq. 1)then
! In the inertial frame
!----------------------
     fx=fsrp(1)
     fy=fsrp(2)
     fz=fsrp(3)

     else
! In the orbit frame
!--------------------

     fx=(fsrp(1)*(-er(1))+fsrp(2)*(-er(2))+fsrp(3)*(-er(3)))
     fy=(fsrp(1)*ex(1)   +fsrp(2)*ex(2)   +fsrp(3)*ex(3))
     fz=(fsrp(1)*en(1)   +fsrp(2)*en(2)   +fsrp(3)*en(3))

!     fx=-(fsrp(1)*ed(1)   +fsrp(2)*ed(2)   +fsrp(3)*ed(3))
!     fy=-(fsrp(1)*ey(1)   +fsrp(2)*ey(2)   +fsrp(3)*ey(3))
!     fz=-(fsrp(1)*eb(1)   +fsrp(2)*eb(2)   +fsrp(3)*eb(3))


     end if

     end if

! end of ECOM (D2B1) model
! ==================================================================

END
