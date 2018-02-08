MODULE m_pd_gfm


! ----------------------------------------------------------------------
! MODULE: m_pd_gfm
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the following subroutines 
! 
! Subroutines contained within the module:
! - pd_gfm
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	11 December 2017
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE pd_gfm (GM, ae, r, nmax, mmax, Cnm, Snm, dCnm, dSnm , fx,fy,fz, Umatrix)


! ----------------------------------------------------------------------
! SUBROUTINE: pd_gfm
! ----------------------------------------------------------------------
! Purpose:
! Partial derivatives of the acceleration due to the Earth Gravity field 
! (based on a gravity model expressed by a set of spherical harmonic coefficients)
! and the Tides effects (expressed as corrections to the spherical harmonic coefficients)
! ----------------------------------------------------------------------
! Input arguments:
! - GM:				Earth gravity constant  (m^3/sec^2)
! - ae:				Earth radius (meters)
! - r:				Position vector (m) in Terrestrial Reference System (ITRS)
!   				r = [x y z]
! - n_max:          maximum degree expansion
! - m_max:          maximum order expansion (m_max <= n_max)
! - Cnm, Snm:		Spherical Harmonic Coefficients (degree n, order m); dynamic allocatable arrays
! - dCnm, dSnm:		Tidal corrections to the spherical harmonic coefficients (degree n, order m); dynamic allocatable arrays
! 
! Output arguments:
! - fx,fy,fz:		Acceleration's cartesian components in ITRS
! - U: 				Matrix of the partial derivatives of the acceleration in ITRS
!   				given in the form of second partial derivatives of 
!					geopotential Uxx, Uyy, Uzz, Uxy, Uxz, Uyz 
! ----------------------------------------------------------------------
! Remark 1:
!  Cnm, Snm are dynamic allocatable arrays
!  Cnm and Snm arrays are formed into lower triangular matrices.
!  Coefficient Cnm corresponds to the matrix element Cnm(n+1,m+1)
!
! Remark 2:
!  Matrix of second partial derivatives:
!   U = [ Uxx   Uxy   Uxz 
!         Uxy   Uyy   Uyz
!         Uxz   Uyz   Uzz ]  
! 
!  Partial derivatives are noted as:
!   dfx / dy = Uxy,     dfx / dx = Uxx
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	11 December 2017
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_legendre
      USE m_legendre1
      USE m_legendre2
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_q), INTENT(IN) :: GM, ae
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(3) :: r
      INTEGER (KIND = prec_int8), INTENT(IN) :: n_max, m_max
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:), ALLOCATABLE :: Cnm, Snm
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_q), INTENT(OUT) :: fx,fy,fz
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q), DIMENSION(:,:), ALLOCATABLE :: Pnm_norm, dPnm_norm
      REAL (KIND = prec_q) :: phi, lamda, radius, fgrav_3d
      REAL (KIND = prec_q) :: fr, ftheta, flamda
      REAL (KIND = prec_q) :: dV_r, dV_phi, dV_lamda
      INTEGER (KIND = prec_int8) :: n, m, m_limit
      INTEGER (KIND = prec_int8) :: sz_tides, Nmax_tide
      INTEGER (KIND = prec_int2) :: comp_option
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
! ----------------------------------------------------------------------






! ----------------------------------------------------------------------
! Tides corrections Nmax
sz_tides = SIZE (dCnm, DIM=1)
Nmax_tide = sz_tides - 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! computation of spherical coordinates in radians
![lamda,phi,l] = lamda_phi(r);
!rdist = l;
CALL coord_r2sph (r , phi,lamda,radius)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Normalized associated Legendre functions
CALL legendre (phi, n_max, Pnm_norm)

! First-order derivatives of normalized associated Legendre functions
CALL legendre_drv1 (phi, n_max, dPnm_norm)

! Second-order derivatives of the Normalized Associated Legendre functions
![d2Pnm_norm] = Legendre2ord(phi,nmax) ;
CALL legendre2 (phi, nmax, d2Pnm_norm)
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
!% 2nd approach:
! Partial derivatives of potential with respect to the spherical coordinates:
! - dV_r     : partial derivative of geopotential to radius
! - dV_phi   : partial derivative of geopotential to latitude
! - dV_lamda : partial derivative of geopotential to longtitude
dV_r = 0.D0
dV_phi = 0.D0
dV_lamda = 0.D0
Do n = 2 , nmax
    If (n > mmax) Then
        m_limit = mmax
    else
        m_limit = n
    End IF
    if (n > Nmax_tide) then
        Do m = 0 , m_limit    
            dV_r = dV_r         + (n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda));
            dV_phi = dV_phi     + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda)+Snm(n+1,m+1)*sin(m*lamda));
            dV_lamda = dV_lamda + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda)-Cnm(n+1,m+1)*sin(m*lamda));
        end Do
    else
        Do m = 0 , m_limit    
            dV_r = dV_r         + (n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) * (Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda)) ...
                                + (n+1)*((ae/l)^n) * Pnm_norm(n+1,m+1) *(dCnm(n+1,m+1) * cos(m*lamda) +dSnm(n+1,m+1) * sin(m*lamda)) ;
            dV_phi = dV_phi     + ((ae/l)^n) * dPnm_norm(n+1,m+1) * (Cnm(n+1,m+1)*cos(m*lamda) + Snm(n+1,m+1)*sin(m*lamda)) ...
                                + ((ae/l)^n) * dPnm_norm(n+1,m+1) *(dCnm(n+1,m+1)*cos(m*lamda) +dSnm(n+1,m+1)*sin(m*lamda)) ;
            dV_lamda = dV_lamda + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) * (Snm(n+1,m+1)*cos(m*lamda) - Cnm(n+1,m+1)*sin(m*lamda)) ...
                                + m * ((ae/l)^n) * Pnm_norm(n+1,m+1) *(dSnm(n+1,m+1)*cos(m*lamda) -dCnm(n+1,m+1)*sin(m*lamda)) ;
        end Do
    end If
end Do
dV_r = - GM / (l**2) - (GM/l**2) * dV_r
dV_phi = (GM / l) * dV_phi
dV_lamda = (GM / l) * dV_lamda
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Partial derivatives of (r,phi,lamda) with respect to (x,y,z)
! ----------------------------------------------------------------------
PDVrx (1,1:3) = (/ cos(phi)*cos(lamda) ,              cos(phi)*sin(lamda) ,             sin(phi) /) 
PDVrx (2,1:3) = (/ (1.0D0/l)*sin(phi)*cos(lamda) ,    (1.0D0/l)*sin(phi)*sin(lamda),    (-1.0D0/l)*cos(phi) /)
PDVrx (3,1:3) = (/ ( -1.D0/(l*cos(phi)) )*sin(lamda), ( 1.D0/(l*cos(phi)) )*cos(lamda),         0.0D0 /)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Computation of Cartesian counterparts of the acceleration
!fxyz = PDVrx' * [dV_r; dV_phi; dV_lamda];

PDVrx_T = TRANSPOSE (PDVrx)
dV_sph = (/ dV_r, dV_phi, dV_lamda /)
CALL matrix_Rr (PDVrx_T, dV_sph , fxyz)

fx = fxyz(1)
fy = fxyz(2)
fz = fxyz(3)
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
! Geopotential 2nd-order partial derivatives in r,theta,lamda components
! ----------------------------------------------------------------------
! Vrr, VrTheta (Vrt), VrLamda(VrL), V_ThetaTheta (Vtt), V_ThetaLamda (VtL)
! V_LamdaLamda (VLL)
! ----------------------------------------------------------------------
Vrr = 0.D0
Vrt = 0.D0
VrL = 0.D0
Vtt = 0.D0 
VtL = 0.D0 
VLL = 0.D0

Do n = 0 , nmax
    if (n > Nmax_tide) then        
        Do m = 0 , n        
            Vrr_f = (GM / ae**3) * (n+1)*(n+2) * (ae / l)**(n+3)
            Vrr = Vrr + Vrr_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1)

            Vrt_f = - (GM / ae**2) * (n+1) * (ae / l)**(n+2) 					!% Eq.6.24 Revised 2012
            !% Abramowitz and Stegun 1972 (first-order derivative of Pnm)
            !%dPnm = n*tan(phi) * Pnm_norm(n+1,m+1) - (n+m) * (1/cos(phi)) * Pnm_norm(n-1+1,m+1)  ;
            Vrt = Vrt + Vrt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * dPnm_norm(n+1,m+1)
            !%Pnm_dv1 = dPnm_norm(n+1,m+1)

            VrL_f = (GM / ae**2) * (n+1) * (ae / l)**(n+2)
            VrL = VrL + VrL_f * m * ( Cnm(n+1,m+1) * sin(m*lamda) - Snm(n+1,m+1) * cos(m*lamda) ) * Pnm_norm(n+1,m+1)

            Vtt_f = (GM / ae) * (ae / l)**(n+1)
            Vtt = Vtt + Vtt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * d2Pnm_norm(n+1,m+1)
            !%Pnm_dv2 = d2Pnm_norm(n+1,m+1)
    !%         % Abramowitz and Stegun 1972 (first-order derivative of Pnm)
    !%         %dPnm = n*tan(phi) * Pnm_norm(n+1,m+1) - (n+m) * (1/cos(phi)) * Pnm_norm(n-1+1,m+1)  ;
    !%         % Hobson 1931 (second-order derivative of Pnm)
    !%         dPnm = dPnm_norm(n+1,m+1);
    !%         dPnm2 = - tan(phi) * dPnm - ( n*(n+1) - m^2 * (1/cos(phi))^2 ) * Pnm_norm(n+1,m+1) ;
    !%         Vtt = Vtt + Vtt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * dPnm2;  %d2Pnm_norm(n+1,m+1) ;


            VtL_f = (GM / ae) * (ae / l)**(n+1)
            !% Corrected : May 2012
            VtL = VtL + VtL_f * m * ( - Cnm(n+1,m+1) * sin(m*lamda) + Snm(n+1,m+1) * cos(m*lamda) ) * dPnm_norm(n+1,m+1)

            VLL_f = (GM / ae) * (ae / l)**(n+1) 
            VLL = VLL - VLL_f * m**2 * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1) 
        end Do
    else
        Do m = 0 , n        
            Vrr_f = (GM / ae**3) * (n+1)*(n+2) * (ae / l)**(n+3) 
            Vrr = Vrr + Vrr_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1) &
                      + Vrr_f * (dCnm(n+1,m+1) * cos(m*lamda) +dSnm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1)  

            
            Vrt_f = - (GM / ae**2) * (n+1) * (ae / l)**(n+2)  				!% Eq.6.24 Revised 2012
            !% Abramowitz and Stegun 1972 (first-order derivative of Pnm)
            !%dPnm = n*tan(phi) * Pnm_norm(n+1,m+1) - (n+m) * (1/cos(phi)) * Pnm_norm(n-1+1,m+1)  ;
            Vrt = Vrt + Vrt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * dPnm_norm(n+1,m+1) &
                      + Vrt_f * (dCnm(n+1,m+1) * cos(m*lamda) +dSnm(n+1,m+1) * sin(m*lamda) ) * dPnm_norm(n+1,m+1) 

            
            VrL_f = (GM / ae**2) * (n+1) * (ae / l)**(n+2) 
            VrL = VrL + VrL_f * m * ( Cnm(n+1,m+1) * sin(m*lamda) - Snm(n+1,m+1) * cos(m*lamda) ) * Pnm_norm(n+1,m+1) &
                      + VrL_f * m * (dCnm(n+1,m+1) * sin(m*lamda) -dSnm(n+1,m+1) * cos(m*lamda) ) * Pnm_norm(n+1,m+1) 

            
            Vtt_f = (GM / ae) * (ae / l)*(n+1) 
            Vtt = Vtt + Vtt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * d2Pnm_norm(n+1,m+1) &
                      + Vtt_f * (dCnm(n+1,m+1) * cos(m*lamda) +dSnm(n+1,m+1) * sin(m*lamda) ) * d2Pnm_norm(n+1,m+1) 
            !%Pnm_dv2 = d2Pnm_norm(n+1,m+1)
    !%         % Abramowitz and Stegun 1972 (first-order derivative of Pnm)
    !%         %dPnm = n*tan(phi) * Pnm_norm(n+1,m+1) - (n+m) * (1/cos(phi)) * Pnm_norm(n-1+1,m+1)  ;
    !%         % Hobson 1931 (second-order derivative of Pnm)
    !%         dPnm = dPnm_norm(n+1,m+1);
    !%         dPnm2 = - tan(phi) * dPnm - ( n*(n+1) - m^2 * (1/cos(phi))^2 ) * Pnm_norm(n+1,m+1) ;
    !%         Vtt = Vtt + Vtt_f * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * dPnm2;  %d2Pnm_norm(n+1,m+1) ;


            VtL_f = (GM / ae) * (ae / l)**(n+1) 
            !% Corrected : May 2012
            VtL = VtL + VtL_f * m * ( - Cnm(n+1,m+1) * sin(m*lamda) + Snm(n+1,m+1) * cos(m*lamda) ) * dPnm_norm(n+1,m+1) &
                      + VtL_f * m * ( -dCnm(n+1,m+1) * sin(m*lamda) +dSnm(n+1,m+1) * cos(m*lamda) ) * dPnm_norm(n+1,m+1) 

            
            VLL_f = (GM / ae) * (ae / l)**(n+1) 
            VLL = VLL - VLL_f * m**2 * ( Cnm(n+1,m+1) * cos(m*lamda) + Snm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1) &
                      - VLL_f * m**2 * (dCnm(n+1,m+1) * cos(m*lamda) +dSnm(n+1,m+1) * sin(m*lamda) ) * Pnm_norm(n+1,m+1) 
        end Do
    end If
end Do

!% Geopotential 2nd-order partial derivatives in local tangent frame 
!% er,etheta,elamda
Vrtl(1,:) = (/ Vrr , Vrt , VrL /)
Vrtl(2,:) = (/ Vrt , Vtt , VtL /)   
Vrtl(3,:) = (/ VrL , VtL , VLL /) 	 
! ----------------------------------------------------------------------


!% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% % Partial derivatives of (r,theta,lamda) with respect to (x,y,z)
!% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% % PDVrx = [
!% %       sin(theta)*cos(lamda)            sin(theta)*sin(lamda)           cos(theta)
!% %   ( 1/l)*cos(theta)*cos(lamda)     ( 1/l)*cos(theta)*sin(lamda)    (-1/l)*sin(theta)
!% % ( -1/(l*sin(theta)) )*sin(lamda)  ( 1/(l*sin(theta)) )*cos(lamda)         0
!% % ];
!% % Replacement of "theta" with "phi"
!% PDVrx = [
!%       cos(phi)*cos(lamda)            cos(phi)*sin(lamda)          sin(phi)
!%   ( 1/l)*sin(phi)*cos(lamda)     ( 1/l)*sin(phi)*sin(lamda)    (-1/l)*cos(phi)
!% ( -1/(l*cos(phi)) )*sin(lamda)  ( 1/(l*cos(phi)) )*cos(lamda)        0
!% ];
!% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = pi/2 - phi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Partial derivatives of PDVrx with respect to (r,theta,lamda)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdv2_r (1,1:3) = (/ 0.0D0 , 0.0D0 , 0.0D0 /) 
pdv2_r (2,1:3) = (/ (-1.0D0/rdist^2)*(cos(theta)*cos(lamda)), (-1.0D0/rdist^2)*(cos(theta)*sin(lamda)), (1/rdist^2)*sin(theta) /)
pdv2_r (3,1:3) = (/  (1.0D0/rdist^2)*(sin(lamda)/sin(theta)), (-1.0D0/rdist^2)*(cos(lamda)/sin(theta)),             0.0D0      /)
   
pdv2_theta (1,1:3) = (/ cos(theta)*cos(lamda), cos(theta)*sin(lamda), -sin(theta) /)
pdv2_theta (2,1:3) = (/ (-1.0D0/rdist)*(sin(theta)*cos(lamda)), (-1.0D0/rdist)*(sin(theta)*sin(lamda)), (-1.0D0/rdist)*cos(theta) /)
pdv2_theta (3,1:3) = (/  (1.0D0/rdist)*(sin(lamda)*cos(theta)/sin(theta)^2), (-1.0D0/rdist)*(cos(lamda)*cos(theta)/sin(theta)^2), 0.0D0 /)

pdv2_lamda (1,1:3) = (/ -sin(theta)*sin(lamda), sin(theta)*cos(lamda), 0.0D0 /)
pdv2_lamda (2,1:3) = (/ (-1.0D0/rdist)*(cos(theta)*sin(lamda)),  (1.0D0/rdist)*(cos(theta)*cos(lamda)), 0.0D0 /)
pdv2_lamda (3,1:3) = (/ (-1.0D0/rdist)*(cos(lamda)/sin(theta)), (-1.0D0/rdist)*(sin(lamda)/sin(theta)), 0.0D0 /)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Partial derivatives of (Fx,Fy,Fz) with respect to (r,theta,lamda)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% pdv(fxyz / r)
pdvfx_r = Vrr * PDVrx(1,1) + Vrt * PDVrx(2,1) + VrL * PDVrx(3,1) +           &
          dV_r * pdv2_r(1,1) + dV_phi * pdv2_r(2,1) + dV_lamda * pdv2_r(3,1)
      
pdvfy_r = Vrr * PDVrx(1,2) + Vrt * PDVrx(2,2) + VrL * PDVrx(3,2) +           &
          dV_r * pdv2_r(1,2) + dV_phi * pdv2_r(2,2) + dV_lamda * pdv2_r(3,2)
      
pdvfz_r = Vrr * PDVrx(1,3) + Vrt * PDVrx(2,3) + VrL * PDVrx(3,3) +           &
          dV_r * pdv2_r(1,3) + dV_phi * pdv2_r(2,3) + dV_lamda * pdv2_r(3,3)

!% pdv(fxyz / theta)
pdvfx_theta = Vrt * PDVrx(1,1) + Vtt * PDVrx(2,1) + VtL * PDVrx(3,1) +                      &
              dV_r * pdv2_theta(1,1) + dV_phi * pdv2_theta(2,1) + dV_lamda * pdv2_theta(3,1)
      
pdvfy_theta = Vrt * PDVrx(1,2) + Vtt * PDVrx(2,2) + VtL * PDVrx(3,2) +                      &   
              dV_r * pdv2_theta(1,2) + dV_phi * pdv2_theta(2,2) + dV_lamda * pdv2_theta(3,2)

pdvfz_theta = Vrt * PDVrx(1,3) + Vtt * PDVrx(2,3) + VtL * PDVrx(3,3) +                      &
              dV_r * pdv2_theta(1,3) + dV_phi * pdv2_theta(2,3) + dV_lamda * pdv2_theta(3,3)

!% pdv(fxyz / lamda)
pdvfx_lamda = VrL * PDVrx(1,1) + VtL * PDVrx(2,1) + VLL * PDVrx(3,1) +                      &  
              dV_r * pdv2_lamda(1,1) + dV_phi * pdv2_lamda(2,1) + dV_lamda * pdv2_lamda(3,1)

pdvfy_lamda = VrL * PDVrx(1,2) + VtL * PDVrx(2,2) + VLL * PDVrx(3,2) +                      &
              dV_r * pdv2_lamda(1,2) + dV_phi * pdv2_lamda(2,2) + dV_lamda * pdv2_lamda(3,2)
          
pdvfz_lamda = VrL * PDVrx(1,3) + VtL * PDVrx(2,3) + VLL * PDVrx(3,3) +                      & 
              dV_r * pdv2_lamda(1,3) + dV_phi * pdv2_lamda(2,3) + dV_lamda * pdv2_lamda(3,3)
          
!% matrix : pdvFxyz_rtl
pdvFxyz_rtl = [
pdvFxyz_rtl (1,1:3) = (/ pdvfx_r,   pdvfx_theta,   pdvfx_lamda /)
pdvFxyz_rtl (2,1:3) = (/ pdvfy_r,   pdvfy_theta,   pdvfy_lamda /)
pdvFxyz_rtl (3,1:3) = (/ pdvfz_r,   pdvfz_theta,   pdvfz_lamda /)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Geopotential 2nd-order partial derivatives in ITRS X,Y,Z components
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Umatrix = PDVrx' * pdvFxyz_rtl';

! Module mdl_arr.f90
! Allocate arrays
R1 = transp(PDVrx)
R2 = transp(pdvFxyz_rtl)
Call matrix_RxR
Umatrix = R3 

! Rewrite matrix_RxR subroutine in Fortran 2003 for adavnces in dynamic memory allocation (allocatable arrays as dummy arguments)
! or check the Fortran library subroutines/functions re arrays multiplication 



END SUBROUTINE



END Module
