
SUBROUTINE force_ant (mjd,prnnum,satsvn,r,fx,fy,fz)


! ----------------------------------------------------------------------
! SUBROUTINE: force_ant.f90
! ----------------------------------------------------------------------
! Purpose:
! Acceleration due to the antenna thrust effect  
! ----------------------------------------------------------------------
! Input arguments:
! - mjd          : time variable
! - prnnum       : satellite PRN number 
! - satsvn       : satellite SVN number  
! - r            : satellite position vector (m)
! 
! Output arguments:
!
! - fx,fy,fz:	 : Accelerations in the inertial frame (m/s^2)
! 
! ----------------------------------------------------------------------
! Author :	Dr. Tzupang Tseng
!
! Created:	03-12-2018
!
! Changes:       
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
      INTEGER                            :: prnnum
      REAL (KIND = prec_q)               :: mjd  
      REAL (KIND = prec_q), DIMENSION(3) :: r
      REAL (KIND = prec_q)               :: fx,fy,fz
!
! Satellite information
!-----------------------------------------------------------------------
      INTEGER (KIND = 4)                 :: satsvn
! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: rsat
      REAL (KIND = prec_q) :: c
      REAL (KIND = prec_q), DIMENSION(3) :: er
! ----------------------------------------------------------------------
! Satellite physical informaiton
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: X_SIDE,Z_SIDE
      REAL (KIND = prec_q) :: MASS,AREA
      REAL (KIND = prec_q) :: A_SOLAR

      REAL (KIND = prec_q) :: w, f_ant
! ----------------------------------------------------------------------
c  = 299792458.d0
! GPS constellation
! -----------------
      if(prnnum.le.100)then
         MASS   = 1080.0d0
         Z_SIDE = 4.25D0
         X_SIDE = 4.11D0
         A_SOLAR= 13.92D0

! IIR-A/B
         if(satsvn.eq.41 .or.satsvn.eq.51 .or.satsvn.eq.54 .or.satsvn.eq.56 &
           .or.satsvn.ge.43.and.satsvn.le.47 .or. satsvn.ge.59.and.satsvn.le.61)then
            w   = 60.0d0
! IIR-M
         elseif(satsvn.eq.52 .or.satsvn.eq.53 .or.satsvn.eq.55 .or.satsvn.eq.57 &
           .or.satsvn.eq.58 .or.satsvn.ge.48.and.satsvn.le.50)then
            w   = 145.0d0
! IIF
         elseif(satsvn.ge.62.and.satsvn.le.73)then
            w   = 240.0d0
         MASS   = 1633.0d0
         Z_SIDE = 5.05D0
         X_SIDE = 4.55D0
         A_SOLAR= 22.25D0

         end if
! GLONASS constellation
! ---------------------
      else if (prnnum .gt. 100 .and. prnnum .le. 200) then
         MASS   = 1415.0d0
         Z_SIDE = 1.6620D0
         X_SIDE = 4.200D0
         A_SOLAR= 23.616D0

! GLONASS-M L1L/L2L
         if(satsvn.eq.735)then
            w   = 20.0d0
! GLONASS-M L1L/L2M
         elseif(satsvn.eq.715 .or.satsvn.eq.721 .or.satsvn.eq.733.or.satsvn.eq.734.or.satsvn.eq.736)then
            w = 25.0d0
! GLONASS-M L1L/L2H
         elseif(satsvn.eq.719)then
            w = 40.0d0
! GLONASS-M L1M/L2H
         elseif(satsvn.eq.716 .or.satsvn.eq.720)then
            w = 60.0d0
! GLONASS-M L1H/L2M
         elseif(satsvn.eq.717 .or.satsvn.eq.730 .or.satsvn.eq.732)then
            w = 65.0d0
! GLONASS-M L1H/L2H
         elseif(satsvn.eq.720 .or.satsvn.eq.731.or.satsvn.eq.747.or.satsvn.eq.851 &
            .or.satsvn.eq.853.or.satsvn.eq.854.or.satsvn.ge.742.and.satsvn.le.745)then
            w = 85.0d0
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
         Z_SIDE = 3.002D0
         X_SIDE = 1.323D0
         A_SOLAR= 11.0D0

! GALILEO IOV
         if(satsvn.ge.101.and.satsvn.le.104)then
            w   = 130.0d0
         MASS   = 695.0d0
! GALILEO FOC
         elseif(satsvn.ge.201.and.satsvn.le.214)then
            w = 265.0d0
         MASS = 707.0d0
         end if
! BDS constellation
! -----------------
      else if (prnnum .gt. 300 .and. prnnum .le. 400) then
         Z_SIDE = 3.96D0
         X_SIDE = 4.5D0
         A_SOLAR= 22.44D0

! BDS MEO
         if(satsvn.ge.12.and.satsvn.le.15)then
            w   = 130.0d0
         MASS   = 800.0d0

! BDS IGSO
         elseif(satsvn.ge.7.and.satsvn.le.10.or.satsvn.eq.5.or.satsvn.eq.17)then
            w = 185.0d0
         MASS = 1400.0d0

         end if
! QZSS constellation
! ------------------
      else if (prnnum .gt. 400 .and. prnnum .le. 500) then
         if(satsvn.eq.1)then
            w = 244.0d0
         MASS = 2000.0d0
         Z_SIDE = 6.00D0 
         X_SIDE = 12.2D0
         A_SOLAR= 40.0D0
         end if
      end if

! --------------------------------------------------------------------

! Radial vector of satellite
! --------------------------
      rsat = sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(1)=r(1)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(2)=r(2)/sqrt(r(1)**2+r(2)**2+r(3)**2)
      er(3)=r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2)
!print *, 'r=', r(1),r(2),r(3)



! The acceleration caused by the satellite antenna thrust
! --------------------------------------------------------
   
      f_ant = w/(MASS*c)

! --------------------------------------------------------

! forces in the inertial frame
!-------------------------------
   fx = (f_ant)*er(1)
   fy = (f_ant)*er(2)
   fz = (f_ant)*er(3)
!print*,fx, fy, fz
   

END SUBROUTINE
