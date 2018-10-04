
SUBROUTINE force_erp (mjd,prnnum,r,vsat,r_sun,fx_erp,fy_erp,fz_erp)


! ----------------------------------------------------------------------
! SUBROUTINE: force_erp.f90
! ----------------------------------------------------------------------
! Purpose:
! Acceleration due to the earth radiation pressure plus antenna thrust effect (simply model) 
! This model is only for GNSS satellites, not for LEO satellites which needs
! more precise model, such as a numerical one.
! ----------------------------------------------------------------------
! Input arguments:
! - prnnum       : satellite PRN number 
! - mjd          : time variable 
! - r            : satellite position vector (m)
! - vsat         : satellite velocity vector
! - r_sun        : Sun position vector
! 
! Output arguments:
! - psi          : angle sat-earth-sun  
! - aa           : Earth irradiance     
! - fx,fy,fz:	 : Accelerations in the inertial frame (m/s^2)
! - fr,ft,fn
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	15-03-2018
!
! Changes:      03-10-2018 Tzupang Tseng: Add the conversion between PRN and SVN
!                                         ADD the antenna thrust info for GNSS 
!
! Copyright:  GEOSCIENCE AUSTRALIA
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE m_satinfo
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
      INTEGER                            :: prnnum
      REAL (KIND = prec_q)               :: mjd  
      REAL (KIND = prec_q), DIMENSION(3) :: r,vsat,r_sun
      REAL (KIND = prec_q)               :: fx_erp,fy_erp,fz_erp
      REAL (KIND = prec_q)               :: fr_erp,ft_erp,fn_erp
!
! Satellite information
!-----------------------------------------------------------------------
      INTEGER (KIND = 4)                 :: satsvn
! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: albedo,Ae,Pi,aa, AU
      REAL (KIND = prec_q) :: Re,c,Esun
      REAL (KIND = prec_q) :: cospsi,psi
      REAL (KIND = prec_q) :: Ds,sclfa,rsat,SS,DDs
      REAL (KIND = prec_q) :: u,v

      REAL (KIND = prec_q), DIMENSION(3) :: er,e_sun,ev,en,ey,ed,ex
      REAL (KIND = prec_q), DIMENSION(4) :: AREA1,REFL1,DIFU1,ABSP1
! ----------------------------------------------------------------------
! Satellite physical informaiton
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: X_SIDE,Z_SIDE
      REAL (KIND = prec_q) :: MASS,AREA
      REAL (KIND = prec_q) :: A_SOLAR

      REAL (KIND = prec_q) :: w, f_ant
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: E1,Evis,Eir
      REAL (KIND = prec_q) :: f_Z_radial,f_solar_radial,f_solar_non_radial
! ----------------------------------------------------------------------
! Sun-related variables
! ----------------------------------------------------------------------
       REAL (KIND = prec_q) :: u_sun,beta,del_u
       REAL (KIND = prec_q), DIMENSION(3) :: r_sun1,r_sun2

! ----------------------------------------------------------------------
! Numerical Constants
      albedo = 0.3 ! albedo coefficient 
      Re = 6.371d6 ! mean earth radian, unit : m 
      Pi = 4*atan(1.0d0)
      Ae = Pi*Re**2
      Esun=1367.d0
      c  = 299792458.d0
      AU = 1.4959787066d11 ! (m)

! ----------------------------------------------------------------------
      CALL satinfo (mjd,prnnum,satsvn)
!print*, mjd,prnnum,satsvn

! GPS constellation
! -----------------
      if(prnnum.le.100)then
! IIR-A/B
         if(satsvn.eq.41 .or.satsvn.eq.51 .or.satsvn.eq.54 .or.satsvn.eq.56 &
           .or.satsvn.ge.43.and.satsvn.le.47 .or. satsvn.ge.59.and.satsvn.le.61)then
            w   = 60.0d0
         MASS   = 1080.0d0
         Z_SIDE = 4.25D0 ! surface-to-mass ratio
         X_SIDE = 4.11D0
         A_SOLAR= 13.92D0

! IIR-M
         elseif(satsvn.eq.52 .or.satsvn.eq.53 .or.satsvn.eq.55 .or.satsvn.eq.57 &
           .or.satsvn.eq.58 .or.satsvn.ge.48.and.satsvn.le.50)then
            w   = 145.0d0
         MASS   = 1080.0d0
         Z_SIDE = 4.25D0 ! surface-to-mass ratio
         X_SIDE = 4.11D0
         A_SOLAR= 13.92D0

! IIF
         elseif(satsvn.ge.62.and.satsvn.le.73)then
            w   = 240.0d0
         MASS   = 1633.0d0
         Z_SIDE = 5.05D0 ! surface-to-mass ratio
         X_SIDE = 4.55D0
         A_SOLAR= 22.25D0

         end if
! GLONASS constellation
! ---------------------
      else if (prnnum .gt. 100 .and. prnnum .le. 200) then
! GLONASS-M L1L/L2L
         if(satsvn.eq.735)then
            w   = 20.0d0
         MASS   = 1415.0d0
         Z_SIDE = 1.6620D0 ! surface-to-mass ratio
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-M L1L/L2M
         elseif(satsvn.eq.715 .or.satsvn.eq.721 .or.satsvn.eq.733.or.satsvn.eq.734.or.satsvn.eq.736)then
            w = 25.0d0
         MASS = 1415.0d0
         Z_SIDE = 1.6620D0 ! surface-to-mass ratio
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-M L1L/L2H
         elseif(satsvn.eq.719)then
            w = 40.0d0
         MASS = 1415.0d0
         Z_SIDE = 1.6620D0 ! surface-to-mass ratio
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-M L1M/L2H
         elseif(satsvn.eq.716 .or.satsvn.eq.720)then
            w = 60.0d0
         MASS = 1415.0d0
         Z_SIDE = 1.6620D0 ! surface-to-mass ratio
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-M L1H/L2M
         elseif(satsvn.eq.717 .or.satsvn.eq.730 .or.satsvn.eq.732)then
            w = 65.0d0
         MASS = 1415.0d0
         Z_SIDE = 1.6620D0 ! surface-to-mass ratio
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-M L1H/L2H
         elseif(satsvn.eq.720 .or.satsvn.eq.731.or.satsvn.eq.747.or.satsvn.eq.851 &
            .or.satsvn.eq.853.or.satsvn.eq.854.or.satsvn.ge.742.and.satsvn.le.745)then
            w = 85.0d0
         MASS = 1415.0d0
         Z_SIDE = 1.6620D0 ! surface-to-mass ratio
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-K
         elseif(satsvn.eq.801)then
            w = 135.0d0
         MASS = 995.0d0
         elseif(satsvn.eq.802)then
            w = 105.0d0
         MASS = 995.0d0
         elseif(satsvn.eq.855)then
            w = 100.0d0
         MASS = 995.0d0
! GLONASS-K end
         else
            w = 50.0d0
         MASS = 1415.0d0
         end if
! GALILEO constellation
! ---------------------
      else if (prnnum .gt. 200 .and. prnnum .le. 300) then
! GALILEO IOV
         if(satsvn.ge.101.and.satsvn.le.104)then
            w   = 130.0d0
         MASS   = 695.0d0
         Z_SIDE = 3.002D0
         X_SIDE = 1.323D0
         A_SOLAR= 11.0D0

! GALILEO FOC
         elseif(satsvn.ge.201.and.satsvn.le.214)then
            w = 265.0d0
         MASS = 707.0d0
         Z_SIDE = 3.002D0
         X_SIDE = 1.323D0
         A_SOLAR= 11.0D0

         end if
! BDS constellation
! -----------------
      else if (prnnum .gt. 300 .and. prnnum .le. 400) then
! BDS MEO
         if(satsvn.ge.12.and.satsvn.le.15)then
            w   = 130.0d0
         MASS   = 800.0d0
         Z_SIDE = 3.96D0 ! surface-to-mass ratio
         X_SIDE = 4.5D0
         A_SOLAR= 22.44D0

! BDS IGSO
         elseif(satsvn.ge.7.and.satsvn.le.10.or.satsvn.eq.5.or.satsvn.eq.17)then
            w = 185.0d0
         MASS = 1400.0d0
         Z_SIDE = 3.96D0 ! surface-to-mass ratio
         X_SIDE = 4.5D0
         A_SOLAR= 22.44D0

         end if
! QZSS constellation
! ------------------
      else if (prnnum .gt. 400 .and. prnnum .le. 500) then
         if(satsvn.eq.1)then
            w = 244.0d0
         MASS = 2000.0d0
         Z_SIDE = 6.00D0 ! surface-to-mass ratio
         X_SIDE = 12.2D0
         A_SOLAR= 40.0D0
         end if
      end if

! --------------------------------------------------------------------

! The first step is to estimate the combination energy of the reflected
! radiation from Sun and the earth thermal radiation (emissivity)

! Radial vector of satellite
! --------------------------
      rsat = sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(1)=r(1)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(2)=r(2)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(3)=r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2)
!print *, 'r=', r(1),r(2),r(3)


! The velocity direction (along-track) of satellite
! -------------------------------------------------------

      ev(1)=vsat(1)/sqrt(vsat(1)**2+vsat(2)**2+vsat(3)**2)
      ev(2)=vsat(2)/sqrt(vsat(1)**2+vsat(2)**2+vsat(3)**2)
      ev(3)=vsat(3)/sqrt(vsat(1)**2+vsat(2)**2+vsat(3)**2)

! Normal to the orbit (cross-track)
! --------------------------------
      Call productcross(er,ev,en)


! The unit vector ed SAT->SUN
! ---------------------------
      DDs=sqrt((r_sun(1)-r(1))**2+(r_sun(2)-r(2))**2+(r_sun(3)-r(3))**2)
      ed(1)=((r_sun(1)-r(1))/DDs)
      ed(2)=((r_sun(2)-r(2))/DDs)
      ed(3)=((r_sun(3)-r(3))/DDs)

! The non-radial force caused by the earth radiation is along the +X direction
! in satellite body frame. 
! Contruct the satellite body-fixed frame
! --------------------------------------
      CALL productcross(-er,ed,ey)
      CALL productcross(er,ey,ex)


! Sun vector from the Earth
! -------------------------

      Ds=sqrt(r_sun(1)**2+r_sun(2)**2+r_sun(3)**2)
      e_sun(1)=r_sun(1)/Ds
      e_sun(2)=r_sun(2)/Ds
      e_sun(3)=r_sun(3)/Ds

!print *, 'r_sun=', r_sun


      CALL productdot(er,e_sun,cospsi)

      psi = acos(cospsi)
! print *, 'psi=', psi 

! Analytical earth radiation model
! --------------------------------

      E1 = (Ae*Esun/rsat**2)

      Evis = E1*((2*albedo)/(3*Pi**2)*((Pi-psi)*cos(psi)+sin(psi)))

      Eir = E1*(1-albedo)/(4*Pi)

      aa = Evis+Eir

!print *, 'aa=',aa

! ----------------------------------------------------------------------


! The second step is to produce the accelerations or forces caused by the earth radiation model.
! In this step, some info of box wing SRP model is introduced, e.g., the surface areas of +Z and the solar panel 
! as well as the optical properties.


!      CALL surfprop(BLKNUM,AREA1,REFL1,DIFU1,ABSP1)

! The radial component
! satellite bus +Z
! ----------------
  u = 0.d0    ! specularity
  v = 0.13d0  ! reflectivity
  sclfa=(AU/Ds)**2
  SS = sclfa/(MASS*c)
   
  f_Z_radial = Z_SIDE*SS*aa*((1+u*v)+2/3*(v-u*v))
!print *, '+Z=',AREA1(3)


! solar panel
! -----------

      if ( psi .lt. Pi/2) then ! solar panel back

      u = 0.5d0
      v = 0.2d0

      f_solar_radial = A_SOLAR*SS*abs(cos(psi))*aa*(1+2/3*(v-u*v)*abs(cos(psi))+(u*v)*cos(2*psi))

      f_solar_non_radial = A_SOLAR*SS*cos(psi)*aa*(2/3*(v-u*v)*sin(psi)+(u*v)*abs(sin(2*psi)))

      else if (psi .ge. Pi/2 ) then ! solar panel front
      u = 0.67d0
      v = 0.24d0

      f_solar_radial = A_SOLAR*SS*abs(cos(psi))*aa*(1+2/3*(v-u*v)*abs(cos(psi))+(u*v)*cos(2*psi))


      f_solar_non_radial = A_SOLAR*SS*cos(psi)*aa*(2/3*(v-u*v)*sin(psi)+(u*v)*abs(sin(2*psi)))
      end if


! The acceleration caused by the satellite antenna thrust
! --------------------------------------------------------
   
      f_ant = w/(MASS*c)

! --------------------------------------------------------

! forces in the inertial frame
!-------------------------------
   fx_erp = -((f_Z_radial+f_solar_radial+f_ant)*er(1)+f_solar_non_radial*ex(1))
   fy_erp = -((f_Z_radial+f_solar_radial+f_ant)*er(2)+f_solar_non_radial*ex(2))
   fz_erp = -((f_Z_radial+f_solar_radial+f_ant)*er(3)+f_solar_non_radial*ex(3))

   
! forces in the orbital frame
! ---------------------------
!   fr_erp = fx_erp*er(1)+fy_erp*er(2)+fz_erp*er(3)
!   ft_erp = fx_erp*ev(1)+fy_erp*ev(2)+fz_erp*ev(3)
!   fn_erp = fx_erp*en(1)+fy_erp*en(2)+fz_erp*en(3)

END SUBROUTINE
